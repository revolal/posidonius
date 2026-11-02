// Tests of the coplanar form of the Kaula tide ((8.19) of docs/tides/kaula_tides.pdf,
// src/effects/tides/kaula/kaula_tide.rs). Expected values (tests/data/kaula_planet/kaula_coplanar_values.json)
// computed independently by a Python implementation of (8.19), with the Hansen coefficients by quadrature.
use posidonius::kaula_tide::{
    self, COPLANAR_INCLINATION, MeanRates, NodalGeometry, love_number_at, nodal_geometry,
    tide, tide_coplanar, tide_single, prepare_spectrum, q_max,
};
use posidonius::{Axes, LoveNumber};
use serde_json::Value;
use std::f64::consts::PI;

const DAY: f64 = 86400.;
const DATA: &str = "tests/data/kaula_planet/";

fn load(name: &str) -> Value {
    let text = std::fs::read_to_string(format!("{DATA}{name}")).unwrap();
    serde_json::from_str(&text).unwrap()
}

fn axes(value: &Value) -> Axes {
    Axes::from(
        value[0].as_f64().unwrap(),
        value[1].as_f64().unwrap(),
        value[2].as_f64().unwrap(),
    )
}

fn f(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

fn vec_f64(value: &Value) -> Vec<f64> {
    value.as_array().unwrap().iter().map(f).collect()
}

fn distance(a: Axes, b: Axes) -> f64 {
    let mut d = a;
    d.sub(&b);
    d.norm()
}

fn angle_difference(a: f64, b: f64) -> f64 {
    ((a - b + PI).rem_euclid(2. * PI) - PI).abs()
}

fn unit(v: Axes) -> Axes {
    let mut u = v;
    u.div(v.norm());
    u
}

// Tabulated spectrum padded to 1024 rows by repeating the last row (tests/common)
fn love_number_from_table(chi: &[f64], rek2: &[f64], imk2: &[f64]) -> LoveNumber {
    let mut x = [0.; 1024];
    let mut re = [0.; 1024];
    let mut im = [0.; 1024];
    for j in 0..1024 {
        let k = j.min(chi.len() - 1);
        x[j] = chi[k];
        re[j] = rek2[k];
        im[j] = imk2[k];
    }
    LoveNumber::new_from(x, re, im, None)
}

// Analytic spectra of the reference file, chi in rad/s: bar k_2(sigma), sigma in rad/day
fn analytic_love(spectrum: &Value) -> impl Fn(f64) -> (f64, f64) + '_ {
    move |sigma: f64| {
        let chi = sigma.abs() / DAY;
        let k2 = f(&spectrum["k2"]);
        let (rek2, imk2) = match spectrum["kind"].as_str().unwrap() {
            "constant_time_lag" => {
                let x = chi * f(&spectrum["dt_s"]);
                (k2 * x.cos(), k2 * x.sin())
            }
            "maxwell" => {
                let x = chi * f(&spectrum["tau_s"]);
                (k2 / (1. + x * x), k2 * x / (1. + x * x))
            }
            kind => panic!("unknown spectrum {kind}"),
        };
        if sigma == 0. {
            (rek2, 0.)
        } else {
            (rek2, -sigma.signum() * imk2)
        }
    }
}

struct Inputs {
    r: Axes,
    v: Axes,
    spin: Axes,
    star_mass: f64,
    planet_mass: f64,
    planet_radius: f64,
    g: f64,
    rates: MeanRates,
}

fn inputs(value: &Value) -> Inputs {
    let rates = &value["mean_rates"];
    Inputs {
        r: axes(&value["r"]),
        v: axes(&value["v"]),
        spin: axes(&value["spin"]),
        star_mass: f(&value["star_mass"]),
        planet_mass: f(&value["planet_mass"]),
        planet_radius: f(&value["planet_radius"]),
        g: f(&value["G"]),
        rates: MeanRates {
            n_bar: f(&rates["nbar"]),
            lambda_dot: f(&rates["lamdot"]),
            node_dot: f(&rates["Omdot"]),
        },
    }
}

fn geometry_of(inp: &Inputs, spin: Axes) -> NodalGeometry {
    nodal_geometry(inp.r, inp.v, spin, inp.g * (inp.star_mass + inp.planet_mass))
}

#[derive(Clone, Copy, PartialEq)]
enum Form {
    Dispatch,
    Coplanar,
    Single,
}

fn run<L: FnMut(usize, usize, i32, f64) -> (f64, f64)>(
    inp: &Inputs,
    geo: &NodalGeometry,
    q: usize,
    form: Form,
    love: L,
) -> kaula_tide::Tide {
    let (r, rates, g, ms, rp) = (inp.r, &inp.rates, inp.g, inp.star_mass, inp.planet_radius);
    match form {
        Form::Dispatch => tide(r, geo, rates, love, true, g, ms, rp, q),
        Form::Coplanar => tide_coplanar(r, geo, rates, love, true, g, ms, rp, q),
        Form::Single => tide_single(r, geo, rates, love, true, g, ms, rp, q),
    }
}

// Tide of the given form with the spectrum of the group (analytic, or tabulated through love_number_at)
fn tide_of_form(inp: &Inputs, geo: &NodalGeometry, spectrum: &Value, q: usize, form: Form) -> kaula_tide::Tide {
    if spectrum["kind"] == "tabulated" {
        let love = love_number_from_table(
            &vec_f64(&spectrum["chi"]),
            &vec_f64(&spectrum["rek2"]),
            &vec_f64(&spectrum["imk2"]),
        );
        let spec = prepare_spectrum(&love);
        let mut hints = [[0usize; 3]; 3];
        run(inp, geo, q, form, |m, p, _q, sigma| {
            love_number_at(&love, &spec, sigma, &mut hints[m][p])
        })
    } else {
        let analytic = analytic_love(spectrum);
        run(inp, geo, q, form, |_, _, _, sigma| analytic(sigma))
    }
}

// Spin of the same rate at obliquity i, tilted about a fixed direction perpendicular to h^
fn tilted_spin(inp: &Inputs, i: f64) -> Axes {
    let h = unit(inp.r.cross(&inp.v));
    let t = unit(h.cross(&Axes::from(0.3, -0.7, 0.2)));
    let w = inp.spin.norm();
    Axes::from(
        w * (i.cos() * h.x + i.sin() * t.x),
        w * (i.cos() * h.y + i.sin() * t.y),
        w * (i.cos() * h.z + i.sin() * t.z),
    )
}

//------------------------------------------------------------------------------ against Python

#[test]
fn coplanar_form_against_python() {
    // (8.19) (tide_coplanar, states with i in {0, 1e-9, 1e-7, 5e-6}) == Python, 72 states:
    // force, torque and their dissipative parts (8.21); true anomaly of the geometry. The choice of the
    // form by tide is tested by dispatch_between_the_two_forms
    let doc = load("kaula_coplanar_values.json");
    let mut worst = [0f64; 5];
    for group in doc["groups"].as_array().unwrap() {
        let spectrum = &group["love_number_spectrum"];
        for state in group["states"].as_array().unwrap() {
            let inp = inputs(&state["inputs"]);
            let geo = geometry_of(&inp, inp.spin);
            let expected = &state["expected"];
            let q = expected["q_max"].as_u64().unwrap() as usize;
            assert_eq!(q, q_max(geo.e));
            let t = tide_of_form(&inp, &geo, spectrum, q, Form::Coplanar);
            let force = axes(&expected["force"]);
            let scale = force.norm();
            let r = inp.r.norm();
            worst[0] = worst[0].max(distance(t.force, force) / scale);
            worst[1] = worst[1].max(distance(t.torque, axes(&expected["dL_spin_dt"])) / (r * scale));
            worst[2] = worst[2].max(distance(t.dissipative_force, axes(&expected["dissipative_force"])) / scale);
            worst[3] = worst[3].max(distance(t.dissipative_torque, axes(&expected["dissipative_torque"])) / (r * scale));
            // f is ill-conditioned for e -> 0 (error ~ 1e-16 / e): compared for e > 1e-3; the force is not
            if geo.e > 1e-3 {
                worst[4] = worst[4].max(angle_difference(geo.true_anomaly, f(&state["derived"]["true_anomaly"])));
            }
        }
    }
    println!(
        "coplanar vs Python: force {:.1e}, torque {:.1e}, dissipative force {:.1e}, dissipative torque {:.1e}, true anomaly {:.1e}",
        worst[0], worst[1], worst[2], worst[3], worst[4]
    );
    assert!(worst[0] < 1e-13 && worst[1] < 1e-13, "{worst:?}");
    assert!(worst[2] < 1e-13 && worst[3] < 1e-13, "{worst:?}");
    assert!(worst[4] < 1e-12, "{worst:?}");
}

//------------------------------------------------------------------------------ against the single form

#[test]
fn coplanar_form_equals_single_form_at_zero_obliquity() {
    // i = 0 (spin along r x v): tide_coplanar == tide_single (all four vectors)
    let doc = load("kaula_coplanar_values.json");
    let mut worst: f64 = 0.;
    let mut count = 0;
    for group in doc["groups"].as_array().unwrap() {
        let spectrum = &group["love_number_spectrum"];
        for state in group["states"].as_array().unwrap() {
            let inp = inputs(&state["inputs"]);
            let spin = tilted_spin(&inp, 0.);
            let geo = geometry_of(&inp, spin);
            let q = q_max(geo.e);
            let a = tide_of_form(&inp, &geo, spectrum, q, Form::Coplanar);
            let b = tide_of_form(&inp, &geo, spectrum, q, Form::Single);
            let scale = b.force.norm();
            let r = inp.r.norm();
            worst = worst
                .max(distance(a.force, b.force) / scale)
                .max(distance(a.torque, b.torque) / (r * scale))
                .max(distance(a.dissipative_force, b.dissipative_force) / scale)
                .max(distance(a.dissipative_torque, b.dissipative_torque) / (r * scale));
            count += 1;
        }
    }
    println!("coplanar == single form at i = 0, {count} states: {worst:.1e}");
    assert!(worst < 1e-14, "{worst:e}");
}

//------------------------------------------------------------------------------ switch

#[test]
fn switch_continuity() {
    // |F_2D - F_3D| <= C(e) i F0, F0 = 3 kappa G M*^2 R^5 / r*^7 (scale of the truncation guarantee), with
    // C(e) <= 55 ((8.20), any |bar k| <= kappa); jump at the switch i = i_s (1 -+ 1e-9) below
    // the truncation tolerance 1e-5 F0. States of kaula_coplanar_values.json, analytic spectra (kappa = k2).
    let doc = load("kaula_coplanar_values.json");
    let mut worst_coefficient: f64 = 0.;
    let mut worst_switch: f64 = 0.;
    for group in doc["groups"].as_array().unwrap() {
        let spectrum = &group["love_number_spectrum"];
        if spectrum["kind"] == "tabulated" {
            continue;
        }
        let kappa = f(&spectrum["k2"]);
        for state in group["states"].as_array().unwrap() {
            let inp = inputs(&state["inputs"]);
            let mut incs = vec![COPLANAR_INCLINATION * (1. - 1e-9), COPLANAR_INCLINATION * (1. + 1e-9)];
            incs.extend([1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]);
            for i in incs {
                let spin = tilted_spin(&inp, i);
                let geo = geometry_of(&inp, spin);
                let q = q_max(geo.e);
                let a = tide_of_form(&inp, &geo, spectrum, q, Form::Coplanar);
                let b = tide_of_form(&inp, &geo, spectrum, q, Form::Single);
                let f0 = 3. * kappa * inp.g * inp.star_mass.powi(2) * inp.planet_radius.powi(5) / geo.r_star.powi(7);
                let coefficient = distance(a.force, b.force) / (f0 * geo.inclination);
                worst_coefficient = worst_coefficient.max(coefficient);
                if (geo.inclination / COPLANAR_INCLINATION - 1.).abs() < 1e-6 {
                    worst_switch = worst_switch.max(distance(a.force, b.force) / f0);
                }
            }
        }
    }
    println!("switch: |F_2D - F_3D| / (F0 i) <= {worst_coefficient:.3}; jump at the switch / F0 <= {worst_switch:.1e}");
    assert!(worst_coefficient <= 55., "{worst_coefficient}");
    assert!(worst_switch <= 1e-5, "{worst_switch:e}");
}

//------------------------------------------------------------------------------ choice of the form

#[test]
fn dispatch_between_the_two_forms() {
    // tide: coplanar for i <= i_s, single form above and near i = pi; a non-rotating planet
    // (s^ = h^, i = 0) takes the coplanar form. Bit-identical to the selected form.
    let doc = load("kaula_coplanar_values.json");
    let group = &doc["groups"][1];
    let spectrum = &group["love_number_spectrum"];
    let inp = inputs(&group["states"][10]["inputs"]);
    let same = |a: &kaula_tide::Tide, b: &kaula_tide::Tide| a == b;
    for (i, expected) in [
        (0., Form::Coplanar),
        (0.5 * COPLANAR_INCLINATION, Form::Coplanar),
        // margin 1e-6: the inclination rebuilt from the vectors is exact to ~1e-16 rad only (1e-9 relative at i_s)
        (COPLANAR_INCLINATION * (1. - 1e-6), Form::Coplanar),
        (COPLANAR_INCLINATION * (1. + 1e-6), Form::Single),
        (0.3, Form::Single),
        (PI - 1e-9, Form::Single),
    ] {
        let spin = tilted_spin(&inp, i);
        let geo = geometry_of(&inp, spin);
        let q = q_max(geo.e);
        let d = tide_of_form(&inp, &geo, spectrum, q, Form::Dispatch);
        let e = tide_of_form(&inp, &geo, spectrum, q, expected);
        assert!(same(&d, &e), "i = {i}: dispatch differs from the expected form");
    }
    let geo = geometry_of(&inp, Axes::new());
    assert_eq!(geo.inclination, 0.);
    let q = q_max(geo.e);
    assert!(same(
        &tide_of_form(&inp, &geo, spectrum, q, Form::Dispatch),
        &tide_of_form(&inp, &geo, spectrum, q, Form::Coplanar)
    ));
}
