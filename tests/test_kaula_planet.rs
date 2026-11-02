// Tests of the planetary Kaula tide in the single form (src/effects/tides/kaula/kaula_tide.rs).
// Expected values (tests/data/kaula_planet) computed independently by a Python implementation of the same
// equations (docs/tides/kaula_tides.pdf, Chapter 8), with the Hansen coefficients by quadrature.
use posidonius::kaula_tide::{
    self, MeanRates, NodalGeometry, RateHistory, averaging_window, eccentricity_functions_g2,
    inclination_functions_f2, love_number_at, mean_motion_about_host, mean_rates,
    nodal_geometry, tide, tide_tabulated, prepare_spectrum, q_max,
    update_rate_history,
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

#[derive(Clone, Copy)]
enum Padding {
    Zeros,
    RepeatLast,
}

// Tabulated spectrum padded to 1024 rows, as by the Python cases (zeros) or tests/common (repetition)
fn love_number_from_table(chi: &[f64], rek2: &[f64], imk2: &[f64], padding: Padding) -> LoveNumber {
    let mut x = [0.; 1024];
    let mut re = [0.; 1024];
    let mut im = [0.; 1024];
    for j in 0..1024 {
        let k = match padding {
            Padding::Zeros if j >= chi.len() => continue,
            Padding::RepeatLast => j.min(chi.len() - 1),
            Padding::Zeros => j,
        };
        x[j] = chi[k];
        re[j] = rek2[k];
        im[j] = imk2[k];
    }
    LoveNumber::new_from(x, re, im, None)
}

// Analytic spectra of the reference files, chi in rad/s: returns bar k_2(sigma), sigma in rad/day
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
            "elastic" => (k2, 0.),
            kind => panic!("unknown spectrum {kind}"),
        };
        if sigma == 0. {
            (rek2, 0.)
        } else {
            (rek2, -sigma.signum() * imk2)
        }
    }
}

fn rates_from(value: &Value) -> MeanRates {
    MeanRates {
        n_bar: f(&value["nbar"]),
        lambda_dot: f(&value["lamdot"]),
        node_dot: f(&value["Omdot"]),
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
    Inputs {
        r: axes(&value["r"]),
        v: axes(&value["v"]),
        spin: axes(&value["spin"]),
        star_mass: f(&value["star_mass"]),
        planet_mass: f(&value["planet_mass"]),
        planet_radius: f(&value["planet_radius"]),
        g: f(&value["G"]),
        rates: rates_from(&value["mean_rates"]),
    }
}

fn geometry_of(inp: &Inputs) -> NodalGeometry {
    nodal_geometry(inp.r, inp.v, inp.spin, inp.g * (inp.star_mass + inp.planet_mass))
}

// Force and torque with the spectrum of the case (analytic or tabulated with the given padding)
fn tide_of(
    inp: &Inputs,
    geo: &NodalGeometry,
    spectrum: &Value,
    q: usize,
    padding: Padding,
) -> kaula_tide::Tide {
    if spectrum["kind"] == "tabulated" {
        let love = love_number_from_table(
            &vec_f64(&spectrum["chi"]),
            &vec_f64(&spectrum["rek2"]),
            &vec_f64(&spectrum["imk2"]),
            padding,
        );
        let spec = prepare_spectrum(&love);
        tide_tabulated(
            inp.r,
            geo,
            &inp.rates,
            &love,
            &spec,
            inp.g,
            inp.star_mass,
            inp.planet_radius,
            q,
        )
    } else {
        let love = analytic_love(spectrum);
        tide(
            inp.r,
            geo,
            &inp.rates,
            |_, _, _, sigma| love(sigma),
            true,
            inp.g,
            inp.star_mass,
            inp.planet_radius,
            q,
        )
    }
}

// Relative errors of force and torque: |dF| / |F| and |dGamma| / (r |F|) (ch. 9, tolerances)
fn errors(tide: &kaula_tide::Tide, expected: &Value, r: Axes) -> (f64, f64) {
    let force = axes(&expected["force"]);
    let torque = axes(&expected["dL_spin_dt"]);
    let scale = force.norm();
    (
        distance(tide.force, force) / scale,
        distance(tide.torque, torque) / (r.norm() * scale),
    )
}

fn check_derived(geo: &NodalGeometry, derived: &Value, label: &str) -> f64 {
    let mut worst: f64 = 0.;
    let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
    worst = worst.max(rel(geo.a, f(&derived["a"])));
    worst = worst.max((geo.e - f(&derived["e"])).abs());
    worst = worst.max(rel(geo.n_osc, f(&derived["n"])));
    worst = worst.max(rel(geo.spin_rate, f(&derived["spin_rate"])));
    worst = worst.max((geo.inclination - f(&derived["obliquity"])).abs());
    worst = worst.max(angle_difference(geo.u_m, f(&derived["u_M"])));
    // omega and M are undefined at e = 0 (only u_M is): compared for e > 1e-8 only
    if geo.e > 1e-8 {
        worst = worst.max(angle_difference(geo.omega, f(&derived["omega_equatorial"])));
        worst = worst.max(angle_difference(geo.mean_anomaly, f(&derived["mean_anomaly"])));
    }
    assert!(worst < 1e-12, "{label}: derived quantities differ by {worst:e}");
    worst
}

//------------------------------------------------------------------------------ reference cases

#[test]
fn reference_cases_force_and_torque() {
    let doc = load("kaula_l2_reference_cases.json");
    for case in doc["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let inp = inputs(&case["inputs"]);
        let geo = geometry_of(&inp);
        check_derived(&geo, &case["derived"], name);
        let table = &case["expected"]["table_tol_1e-5"];
        // q_max of the file: some cases lie on the edges of the table, where e rebuilt from (r, v) may fall on either side
        let q = table["q_max"].as_u64().unwrap() as usize;
        for padding in [Padding::Zeros, Padding::RepeatLast] {
            let tide = tide_of(&inp, &geo, &case["inputs"]["love_number_spectrum"], q, padding);
            let (ef, et) = errors(&tide, table, inp.r);
            println!("{name:40} q_max {q:2}  force {ef:.1e}  torque {et:.1e}");
            assert!(ef < 1e-12 && et < 1e-12, "{name}: force {ef:e}, torque {et:e}");
            // Truncation error against the converged values (|q| <= 70): (8.15), below 1.1e-6
            let (cf, ct) = errors(&tide, &case["expected"]["converged_q70"], inp.r);
            assert!(cf < 1.1e-6 && ct < 1.1e-6, "{name}: truncation {cf:e}, {ct:e}");
            // Elastic case: no dissipative part at all
            if case["inputs"]["love_number_spectrum"]["kind"] == "elastic" {
                assert_eq!(tide.dissipative_force, Axes::from(0., 0., 0.));
                assert_eq!(tide.dissipative_torque, Axes::from(0., 0., 0.));
            }
        }
    }
}

#[test]
fn reference_case_elastic_closed_form() {
    // Elastic limit of (8.17): F = -3 k2 G M*^2 R^5 r^-7 r^ (test V4), torque zero
    let doc = load("kaula_l2_reference_cases.json");
    let case = &doc["cases"][5];
    assert_eq!(case["name"], "elastic_oblique");
    let inp = inputs(&case["inputs"]);
    let geo = geometry_of(&inp);
    let tide = tide_of(&inp, &geo, &case["inputs"]["love_number_spectrum"], 42, Padding::Zeros);
    let r = inp.r.norm();
    let mut expected = inp.r;
    expected.mul(-3. * 0.3 * inp.g * inp.star_mass.powi(2) * inp.planet_radius.powi(5) / r.powi(8));
    let error = distance(tide.force, expected) / expected.norm();
    assert!(error < 2e-6, "elastic closed form: {error:e}");
    assert!(tide.torque.norm() / (r * expected.norm()) < 2e-6);
}

//------------------------------------------------------------------------------ random states

#[test]
fn random_states_force_and_torque() {
    let doc = load("kaula_l2_random_states.json");
    let mut worst = (0f64, 0f64, 0f64);
    for group in doc["groups"].as_array().unwrap() {
        let spectrum = &group["love_number_spectrum"];
        for (k, state) in group["states"].as_array().unwrap().iter().enumerate() {
            let label = format!("{} #{k}", group["spectrum"]);
            let inp = inputs(&state["inputs"]);
            let geo = geometry_of(&inp);
            let derived = &state["derived"];
            worst.2 = worst.2.max(check_derived(&geo, derived, &label));
            for (axis, key) in [(geo.n_hat, "node"), (geo.y_hat, "y_axis"), (geo.s_hat, "spin_axis")] {
                assert!(distance(axis, axes(&derived[key])) < 1e-13, "{label}: {key}");
            }
            let q = state["expected"]["q_max"].as_u64().unwrap() as usize;
            assert_eq!(q, q_max(geo.e), "{label}: q_max");
            for padding in [Padding::Zeros, Padding::RepeatLast] {
                let tide = tide_of(&inp, &geo, spectrum, q, padding);
                let (ef, et) = errors(&tide, &state["expected"], inp.r);
                assert!(ef < 1e-12 && et < 1e-12, "{label}: force {ef:e}, torque {et:e}");
                worst.0 = worst.0.max(ef);
                worst.1 = worst.1.max(et);
            }
        }
    }
    println!(
        "60 random states: force {:.1e}, torque {:.1e}, derived {:.1e}",
        worst.0, worst.1, worst.2
    );
}

//------------------------------------------------------------------------------ functions

#[test]
fn inclination_functions() {
    let doc = load("kaula_planet_unit_values.json");
    for entry in doc["inclination_functions"].as_array().unwrap() {
        let f2 = inclination_functions_f2(f(&entry["i"]));
        for m in 0..3 {
            for p in 0..3 {
                for key in ["F", "F_general"] {
                    let expected = f(&entry[key][m][p]);
                    assert!((f2[m][p] - expected).abs() < 1e-15, "F_2{m}{p}");
                }
            }
        }
    }
}

#[test]
fn eccentricity_functions() {
    let doc = load("kaula_planet_unit_values.json");
    let mut worst: f64 = 0.;
    for entry in doc["eccentricity_functions"].as_array().unwrap() {
        let e = f(&entry["e"]);
        let qm = entry["q_max"].as_u64().unwrap() as usize;
        assert_eq!(qm, q_max(e));
        let g2 = eccentricity_functions_g2(e, qm);
        for p in 0..2 {
            for (k, q) in (-(qm as i32)..=(qm as i32)).enumerate() {
                let value = g2[p][(q + 42) as usize];
                for key in ["tisserand_K40", "quadrature"] {
                    let d = (value - f(&entry[key][p][k])).abs();
                    worst = worst.max(d);
                    assert!(d < 1e-13, "G_2{p}{q}({e}) vs {key}: {d:e}");
                }
            }
        }
    }
    println!("G_2pq: max absolute difference {worst:.1e}");
}

#[test]
fn eccentricity_table_against_tisserand() {
    // Chebyshev tables (8.16) against Tisserand's formula on 4001 eccentricities
    // in [0, 0.5], all p and |q| <= q_max(e); above 0.5 the direct formula is used (identical values)
    let mut worst: f64 = 0.;
    for k in 0..=4000 {
        let e = 0.5 * f64::from(k) / 4000.;
        let qm = q_max(e);
        let table = eccentricity_functions_g2(e, qm);
        let direct = kaula_tide::eccentricity_functions_g2_tisserand(e, qm);
        for p in 0..2 {
            for q in -(qm as i32)..=(qm as i32) {
                let i = (q + 42) as usize;
                worst = worst.max((table[p][i] - direct[p][i]).abs());
            }
        }
    }
    println!("G_2pq table vs Tisserand: max absolute difference {worst:.1e}");
    assert!(worst < 1e-13, "{worst:e}");
    for e in [0.500001, 0.6, 0.8] {
        assert_eq!(eccentricity_functions_g2(e, 42), kaula_tide::eccentricity_functions_g2_tisserand(e, 42));
    }
}

#[test]
fn spectrum_test_table() {
    let doc = load("kaula_l2_reference_cases.json");
    let spectrum = &doc["cases"][7]["inputs"]["love_number_spectrum"];
    assert_eq!(spectrum["kind"], "tabulated");
    let (chi, re, im) = (vec_f64(&spectrum["chi"]), vec_f64(&spectrum["rek2"]), vec_f64(&spectrum["imk2"]));
    let unit = load("kaula_planet_unit_values.json");
    for padding in [Padding::Zeros, Padding::RepeatLast] {
        let love = love_number_from_table(&chi, &re, &im, padding);
        let spec = prepare_spectrum(&love);
        assert_eq!(spec.n_distinct, chi.len());
        assert!(!spec.has_low);
        // chi_1 = 0: no lower bound for the warning
        assert!(!kaula_tide::far_outside_spectrum(&love, &spec, 1e-30));
        assert!(kaula_tide::far_outside_spectrum(&love, &spec, 11. * chi[chi.len() - 1]));
        let mut hint = 0;
        for row in unit["spectrum_test_table"]["evaluations"].as_array().unwrap() {
            let x = f(&row[0]);
            let (r, i) = kaula_tide::evaluate_spectrum(&love, &spec, x, &mut hint);
            assert!((r - f(&row[1])).abs() <= 1e-14 * f(&row[1]).abs(), "rek2({x:e})");
            assert!((i - f(&row[2])).abs() <= 1e-13 * f(&row[2]).abs().max(1e-300), "imk2({x:e})");
            // bar k_2 at sigma of both signs (8.13)
            let sigma = x * DAY;
            let (r_s, i_s) = kaula_tide::evaluate_spectrum(&love, &spec, sigma.abs() / DAY, &mut hint);
            let (kr, ki) = love_number_at(&love, &spec, sigma, &mut hint);
            assert_eq!((kr, ki), (r_s, if sigma == 0. { 0. } else { -i_s }));
            let (_, ki_neg) = love_number_at(&love, &spec, -sigma, &mut hint);
            assert_eq!(ki_neg, if sigma == 0. { 0. } else { i_s });
        }
    }
}

#[test]
fn spectrum_maxwell_cut_and_continuity() {
    let unit = load("kaula_planet_unit_values.json");
    let table = &unit["spectrum_maxwell_cut"];
    let (chi, re, im) = (vec_f64(&table["chi"]), vec_f64(&table["rek2"]), vec_f64(&table["imk2"]));
    let love = love_number_from_table(&chi, &re, &im, Padding::Zeros);
    let spec = prepare_spectrum(&love);
    assert!(spec.has_low);
    for (end, key) in [(spec.high, "high_end"), (spec.low, "low_end")] {
        assert!((end.beta - f(&table[key]["beta"])).abs() < 1e-14, "{key} beta");
        assert!((end.beta_r - f(&table[key]["betaR"])).abs() < 1e-12, "{key} b");
    }
    let mut hint = 0;
    for row in table["evaluations"].as_array().unwrap() {
        let x = f(&row[0]);
        let (r, i) = kaula_tide::evaluate_spectrum(&love, &spec, x, &mut hint);
        // Far above the table rek2 = rek2_N + (rek2_N - rek2_N-1)(...) is a difference of close numbers:
        // absolute tolerance at the rounding level of the table values (k2 = 0.5)
        assert!((r - f(&row[1])).abs() <= 1e-12 * f(&row[1]).abs() + 1e-15 * 0.5, "rek2({x:e}): {r} {}", row[1]);
        assert!((i - f(&row[2])).abs() <= 1e-12 * f(&row[2]).abs().max(1e-300), "imk2({x:e})");
    }
    // Continuity at both ends of the table: relative jump at rounding level
    let n = chi.len();
    for (edge, inner) in [(chi[n - 1], chi[n - 1] * (1. - 1e-12)), (chi[0], chi[0] * (1. + 1e-12))] {
        let outer = if edge == chi[0] { edge * (1. - 1e-12) } else { edge * (1. + 1e-12) };
        let a = kaula_tide::evaluate_spectrum(&love, &spec, inner, &mut hint);
        let b = kaula_tide::evaluate_spectrum(&love, &spec, outer, &mut hint);
        assert!((a.0 - b.0).abs() < 1e-10 * a.0.abs(), "rek2 jump at {edge:e}");
        assert!((a.1 - b.1).abs() < 1e-10 * a.1.abs(), "imk2 jump at {edge:e}");
    }
    // Warning criterion: more than a decade outside the table
    for (x, expected) in [(chi[n - 1], false), (9.9 * chi[n - 1], false), (10.1 * chi[n - 1], true),
                          (chi[0], false), (0.11 * chi[0], false), (0.09 * chi[0], true), (0., false), (chi[100], false)] {
        assert_eq!(kaula_tide::far_outside_spectrum(&love, &spec, x), expected, "chi = {x:e}");
    }
    // imk2(0) = 0 below a spectrum with chi_1 > 0
    let (_, i0) = kaula_tide::evaluate_spectrum(&love, &spec, 0., &mut hint);
    assert_eq!(i0, 0.);
}

#[test]
fn truncation_table() {
    let unit = load("kaula_planet_unit_values.json");
    let edges = vec_f64(&unit["q_max_table"]["edges"]);
    let values: Vec<usize> = unit["q_max_table"]["q_max"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    assert_eq!(edges.as_slice(), &kaula_tide::Q_MAX_ECCENTRICITY_EDGES);
    assert_eq!(values.as_slice(), &kaula_tide::Q_MAX_TABLE);
    assert_eq!(q_max(0.), 2);
    assert_eq!(q_max(0.001), 2);
    assert_eq!(q_max(0.0011), 3);
    assert_eq!(q_max(0.5), 42);
    // Above the table: capped, with a warning
    assert_eq!(q_max(0.7), 42);
}

#[test]
#[should_panic(expected = "unbound orbit")]
fn truncation_unbound_orbit() {
    q_max(1.0);
}

//------------------------------------------------------------------------------ invalid spectra

fn invalid_spectrum(chi: &[f64], re: &[f64], im: &[f64]) {
    let love = love_number_from_table(chi, re, im, Padding::Zeros);
    prepare_spectrum(&love);
}

#[test]
#[should_panic(expected = "strictly increasing")]
fn invalid_spectrum_first_row_of_parameters() {
    // A first row of physical parameters left in the file
    invalid_spectrum(&[90., 1e-8, 1e-7, 1e-6], &[2., 0.9, 0.8, 0.7], &[170., 0.1, 0.2, 0.3]);
}

#[test]
#[should_panic(expected = "strictly increasing")]
fn invalid_spectrum_too_short() {
    invalid_spectrum(&[1e-8, 1e-7], &[0.9, 0.8], &[0.1, 0.2]);
}

#[test]
#[should_panic(expected = "not sorted")]
fn invalid_spectrum_not_sorted() {
    invalid_spectrum(&[1e-8, 1e-7, 1e-6, 1e-9, 1e-5], &[0.9, 0.8, 0.7, 0.6, 0.5], &[0.1, 0.2, 0.3, 0.2, 0.1]);
}

#[test]
#[should_panic(expected = "imk2 < 0")]
fn invalid_spectrum_negative_imk2() {
    invalid_spectrum(&[1e-8, 1e-7, 1e-6], &[0.9, 0.8, 0.7], &[0.1, -0.2, 0.3]);
}

#[test]
#[should_panic(expected = "non-finite")]
fn invalid_spectrum_nan() {
    invalid_spectrum(&[1e-8, 1e-7, 1e-6], &[0.9, f64::NAN, 0.7], &[0.1, 0.2, 0.3]);
}

#[test]
fn valid_spectrum_with_both_paddings() {
    for padding in [Padding::Zeros, Padding::RepeatLast] {
        let love = love_number_from_table(&[0., 1e-8, 1e-7, 1e-6], &[0.9, 0.9, 0.8, 0.7], &[0., 0.1, 0.2, 0.3], padding);
        assert_eq!(prepare_spectrum(&love).n_distinct, 4);
    }
}

//------------------------------------------------------------------------------ mean rates

// Synthetic accepted steps: node precessing uniformly about z, lambda = Omega + u_M linear
fn synthetic_geometry(t: f64, lambda_dot: f64, node_dot: f64, n_osc: f64) -> NodalGeometry {
    let node = 0.3 + node_dot * t;
    let u_m = 1.1 + (lambda_dot - node_dot) * t;
    NodalGeometry {
        n_hat: Axes::from(node.cos(), node.sin(), 0.),
        s_hat: Axes::from(0., 0., 1.),
        n_osc,
        u_m: (u_m + PI).rem_euclid(2. * PI) - PI,
        ..Default::default()
    }
}

#[test]
fn mean_rates_linear_angles() {
    let (lambda_dot, node_dot, n) = (6.2, -2e-3, 6.1);
    let window = 3.0;
    let mut history = RateHistory::default();
    // Start-up: lambda_dot = n_bar = n_osc, node_dot = 0
    let start = mean_rates(&history, 5.0);
    assert_eq!((start.n_bar, start.lambda_dot, start.node_dot), (5.0, 5.0, 0.0));
    let dt = 0.013;
    let mut t = 0.;
    let mut seen_full = false;
    for step in 0..2000 {
        let geo = synthetic_geometry(t, lambda_dot, node_dot, n);
        update_rate_history(&mut history, t, &geo, window);
        // Repeated call at the same time (same accepted step): ignored
        update_rate_history(&mut history, t, &geo, window);
        let rates = mean_rates(&history, 5.0);
        if t < window - dt {
            assert_eq!(rates.lambda_dot, 5.0, "start-up at step {step}");
        } else if t > window + window / 16. {
            seen_full = true;
            assert!((rates.lambda_dot - lambda_dot).abs() < 1e-10 * lambda_dot, "{} at {t}", rates.lambda_dot);
            assert!((rates.node_dot - node_dot).abs() < 1e-10 * lambda_dot, "{}", rates.node_dot);
            assert!((rates.n_bar - n).abs() < 1e-14);
        }
        assert!(history.len <= kaula_tide::HISTORY_CAPACITY);
        t += dt;
    }
    assert!(seen_full);
    // Window suddenly 3 times longer than the buffer: shorter window from the oldest sample, still exact
    let geo = synthetic_geometry(t, lambda_dot, node_dot, n);
    update_rate_history(&mut history, t, &geo, 3. * window);
    let rates = mean_rates(&history, 5.0);
    assert!((rates.lambda_dot - lambda_dot).abs() < 1e-10 * lambda_dot);
}

#[test]
fn averaging_window_limits() {
    let n = 2. * PI;
    // single planet: one period; 3:2: one synodic period; close periods: 10 periods
    assert!((averaging_window(n, &[]) - 1.).abs() < 1e-14);
    let n2 = n / (1.5 * 1.002f64.powf(1.5));
    assert!((averaging_window(n, &[n2]) - 2. * PI / (n - n2)).abs() < 1e-12);
    assert!((averaging_window(n, &[n * 0.999]) - 10.).abs() < 1e-14);
    // unbound bodies (n = 0) are ignored
    assert!((averaging_window(n, &[0., n2]) - 2. * PI / (n - n2)).abs() < 1e-12);
    // circular orbit about the host
    let mu: f64 = 2.9e-4;
    let r: f64 = 0.02;
    let v = (mu / r).sqrt();
    let n_c = mean_motion_about_host(Axes::from(r, 0., 0.), Axes::from(0., v, 0.), mu);
    assert!((n_c - (mu / r.powi(3)).sqrt()).abs() < 1e-14 * n_c);
    assert_eq!(mean_motion_about_host(Axes::from(r, 0., 0.), Axes::from(0., 2. * v, 0.), mu), 0.);
}
