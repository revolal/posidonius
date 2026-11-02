// Tests of the stellar Kaula tide (src/effects/tides/kaula/kaula_tide.rs with D = star, P = planet):
// signed spectra (docs/tides/kaula_tides.pdf, Chapter 9).
use posidonius::constants::{DAY, G, M_EARTH};
use posidonius::kaula_tide::{
    self, MeanRates, evaluate_signed_spectrum, love_number_signed, nodal_geometry,
    prepare_signed_spectrum, q_max, tide, tide_tabulated,
};
use posidonius::{Axes, Integrator, LoveNumber, Particle};
use serde_json::Value;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

const STELLAR_FILE: &str = "input/love_numbers/Aurelie_Stellar/alpha0.51600_P1p2_Ek1p5em6.txt";

// Table padded to 1024 rows by repeating the last row (tests/common)
fn table(x: &[f64], re: &[f64], im: &[f64], spin_ref: Option<f64>) -> LoveNumber {
    let (mut fx, mut fre, mut fim) = ([0.; 1024], [0.; 1024], [0.; 1024]);
    for j in 0..1024 {
        let k = j.min(x.len() - 1);
        (fx[j], fre[j], fim[j]) = (x[k], re[k], im[k]);
    }
    LoveNumber::new_from(fx, fre, fim, spin_ref)
}

// Columns of the stellar files: frequency [rad/s], Im k2, Re k2
fn stellar_columns() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let text = std::fs::read_to_string(STELLAR_FILE).unwrap();
    let (mut x, mut re, mut im) = (vec![], vec![], vec![]);
    for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
        let v: Vec<f64> = line.split_whitespace().map(|s| s.parse().unwrap()).collect();
        x.push(v[0]);
        im.push(v[1]);
        re.push(v[2]);
    }
    (x, re, im)
}

#[test]
fn signed_spectrum_inside_the_table_is_the_linear_interpolation() {
    // Rows reproduced exactly, midpoints linear, and across 0 the interpolation between the two
    // neighbours of 0 is kept (Im(0) != 0)
    let (x, re, im) = stellar_columns();
    let love = table(&x, &re, &im, None);
    let spec = prepare_signed_spectrum(&love);
    assert!(spec.signed && spec.n_distinct == x.len());
    let mut hint = 0;
    for j in 0..x.len() - 1 {
        let (r, i) = evaluate_signed_spectrum(&love, &spec, x[j], &mut hint);
        assert!(r == re[j] && i == im[j], "row {j}");
        let mid = 0.5 * (x[j] + x[j + 1]);
        let (r, i) = evaluate_signed_spectrum(&love, &spec, mid, &mut hint);
        let t = (mid - x[j]) / (x[j + 1] - x[j]);
        let expected = (1. - t) * im[j] + t * im[j + 1];
        assert!((i - expected).abs() <= 1e-15 * expected.abs().max(im[j].abs()), "mid {j}");
        assert!((r - ((1. - t) * re[j] + t * re[j + 1])).abs() <= 1e-15, "mid {j}");
    }
    let k = x.partition_point(|&f| f < 0.);
    let (_, i0) = evaluate_signed_spectrum(&love, &spec, 0., &mut hint);
    let t = -x[k - 1] / (x[k] - x[k - 1]);
    assert!((i0 - ((1. - t) * im[k - 1] + t * im[k])).abs() <= 1e-15 * im[k].abs());
    assert!(i0 != 0.);
}

#[test]
fn signed_spectrum_ends_are_power_laws_of_the_absolute_values() {
    // (8.14) in |x| beyond both ends, imk2 of the sign of x, continuous at the ends
    let x = [-8e-4, -4e-4, -2e-4, -1e-4, 1e-4, 2e-4, 4e-4, 8e-4];
    let im = [-1e-3, -2e-3, -4e-3, -8e-3, 6e-3, 3e-3, 1.5e-3, 1e-3];
    let re = [0.4; 8];
    let love = table(&x, &re, &im, None);
    let spec = prepare_signed_spectrum(&love);
    let mut hint = 0;
    // Negative end: |Im| = 1e-3 (|x| / 8e-4)^beta with beta = ln(1/2) / ln(2) = -1
    for factor in [1.0000001, 2., 10.] {
        let (r, i) = evaluate_signed_spectrum(&love, &spec, x[0] * factor, &mut hint);
        let expected = -1e-3 / factor;
        assert!((i - expected).abs() <= 1e-14 * expected.abs(), "{i:e} vs {expected:e}");
        assert!(r == 0.4);
    }
    // Positive end: |Im| = 1e-3 (x / 8e-4)^beta, beta = ln(1/1.5) / ln(2), within [-1, 0)
    let beta = (1. / 1.5f64).ln() / 2f64.ln();
    for factor in [1.0000001, 3.] {
        let (_, i) = evaluate_signed_spectrum(&love, &spec, x[7] * factor, &mut hint);
        let expected = 1e-3 * factor.powf(beta);
        assert!((i - expected).abs() <= 1e-14 * expected, "{i:e} vs {expected:e}");
    }
    let (_, at_low) = evaluate_signed_spectrum(&love, &spec, x[0], &mut hint);
    let (_, at_high) = evaluate_signed_spectrum(&love, &spec, x[7], &mut hint);
    assert!(at_low == im[0] && at_high == im[7]);
}

#[test]
fn reading_a_with_spin_rescaling() {
    // bar k(sigma) = rek2(x) - i imk2(x) (Omega / Omega_ref)^2 at x = sigma Omega_ref / Omega, (9.1);
    // dissipative: sigma Im bar k < 0 outside the
    // interval of the table that contains 0
    let (x, re, im) = stellar_columns();
    let spin_ref = 2. * PI / (1.2 * DAY);
    let love = table(&x, &re, &im, Some(spin_ref));
    let plain = table(&x, &re, &im, None);
    let spec = prepare_signed_spectrum(&love);
    let k = x.partition_point(|&f| f < 0.);
    let mut checked = 0;
    for period_d in [1.2, 1.9, 1.01] {
        let spin = 2. * PI / period_d; // rad/day
        let ratio = spin_ref * DAY / spin;
        for j in 0..=4000 {
            let sigma = (-3. + 6. * f64::from(j) / 4000.) * spin; // rad/day
            let (mut h1, mut h2) = (0, 0);
            let (kr, ki) = love_number_signed(&love, &spec, sigma, spin, &mut h1);
            let (fr, fi) = evaluate_signed_spectrum(&plain, &spec, sigma / DAY * ratio, &mut h2);
            // Same reading up to the rounding of the rescaled frequency
            assert!((kr - fr).abs() <= 1e-14 * fr.abs(), "sigma {sigma:e}: {kr:e} vs {fr:e}");
            assert!((ki + fi / (ratio * ratio)).abs() <= 1e-12 * fi.abs() / (ratio * ratio), "sigma {sigma:e}: {ki:e} vs {:e}", -fi / (ratio * ratio));
            let xf = sigma / DAY * ratio;
            if !(x[k - 1] < xf && xf < x[k]) {
                assert!(sigma * ki <= 0., "sigma {sigma:e}: Im bar k = {ki:e}");
                checked += 1;
            }
        }
    }
    assert!(checked > 11_000);
}

#[test]
fn hermitian_flag_does_not_change_a_hermitian_spectrum() {
    // With bar k(-sigma) = conj bar k(sigma), evaluating every m = 0 term (signed path) gives the same tide
    // as the evaluation by pairs of Section 8.9, in both forms (coplanar and single)
    let love = |sigma: f64| {
        let w = sigma / 0.37;
        (0.3 + 0.05 / (1. + w * w), -0.02 * w / (1. + w * w))
    };
    let (g, ms, rp) = (0.01720209895f64.powi(2), 0.08, 4.6e-5);
    for (inc, e) in [(0f64, 0.05), (0.3, 0.2), (1.2, 0.4)] {
        let a: f64 = 0.03;
        let n = (g * (ms + 3e-6) / a.powi(3)).sqrt();
        let rp_ = a * (1. - e);
        let v = (g * (ms + 3e-6) * (1. + e) / rp_).sqrt();
        let r = Axes::from(rp_ * 0.6, rp_ * 0.8, 0.);
        let vel = Axes::from(-v * 0.8, v * 0.6, 0.);
        let spin = Axes::from(0., inc.sin() * 1.7 * n, inc.cos() * 1.7 * n);
        let geo = nodal_geometry(r, vel, spin, g * (ms + 3e-6));
        let rates = MeanRates { n_bar: n, lambda_dot: n, node_dot: 0. };
        let q = q_max(e);
        let t1 = tide(r, &geo, &rates, |_, _, _, s| love(s), true, g, ms, rp, q);
        let t2 = tide(r, &geo, &rates, |_, _, _, s| love(s), false, g, ms, rp, q);
        let scale = t1.force.norm();
        let mut d = t1.force;
        d.sub(&t2.force);
        assert!(d.norm() <= 1e-14 * scale, "i = {inc}: {:e}", d.norm() / scale);
        let mut dt = t1.torque;
        dt.sub(&t2.torque);
        assert!(dt.norm() <= 1e-14 * scale * r.norm(), "i = {inc}");
    }
}

#[test]
#[should_panic(expected = "has the sign opposite to the frequency")]
fn incoherent_signed_spectrum_is_rejected() {
    let x = [-3e-4, -2e-4, -1e-4, 1e-4, 2e-4, 3e-4];
    let im = [-6e-3, -4e-3, 2e-3, 2e-3, 4e-3, 6e-3]; // row 2: Im > 0 at a negative frequency
    prepare_signed_spectrum(&table(&x, &[0.4; 6], &im, None));
}

#[test]
#[should_panic(expected = "frequencies of both signs needed")]
fn unsigned_table_is_rejected_as_signed() {
    let x = [0., 1e-5, 2e-5, 4e-5];
    let im = [0., 1e-3, 2e-3, 3e-3];
    prepare_signed_spectrum(&table(&x, &[0.3; 4], &im, None));
}

#[test]
fn prepared_spectrum_kind() {
    let (x, re, im) = stellar_columns();
    assert!(prepare_signed_spectrum(&table(&x, &re, &im, None)).signed);
    let unsigned = kaula_tide::prepare_spectrum(&table(&[0., 1e-5, 2e-5, 4e-5], &[0.3; 4], &[0., 1e-3, 2e-3, 3e-3], None));
    assert!(!unsigned.signed);
}

//------------------------------------------------------------------------------ against Python and the exact response

fn axes(value: &Value) -> Axes {
    let v: Vec<f64> = value.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
    Axes::from(v[0], v[1], v[2])
}

fn distance(a: Axes, b: Axes) -> f64 {
    let mut d = a;
    d.sub(&b);
    d.norm()
}

fn kwok_star_love() -> LoveNumber {
    let (x, re, im) = stellar_columns();
    table(&x, &re, &im, Some(2. * PI / (1.2 * DAY)))
}

#[test]
fn stellar_tide_against_python_and_truth() {
    // 16 random states (D = star, P = planet, Kwok et al. spectrum with rescaling, small precession):
    // the Rust tide at q_max(e) equals the Python single form with the same reading of the spectrum, and at
    // q = 42 the exact linear response computed without Kaula's expansion (tests/data/kaula_stellar)
    let text = fs::read_to_string("tests/data/kaula_stellar/kaula_stellar_values.json").unwrap();
    let data: Value = serde_json::from_str(&text).unwrap();
    let love = kwok_star_love();
    assert_eq!(
        data["spectrum"]["stellar_spectrum_spin_rate"].as_f64().unwrap(),
        2. * PI / (1.2 * DAY)
    );
    let spec = prepare_signed_spectrum(&love);
    let (mut worst_form, mut worst_truth) = (0f64, 0f64);
    for case in data["cases"].as_array().unwrap() {
        let inp = &case["inputs"];
        let (r, v, s) = (axes(&inp["r"]), axes(&inp["v"]), axes(&inp["spin"]));
        let (ms, rs, mp, g) = (
            inp["star_mass"].as_f64().unwrap(),
            inp["star_radius"].as_f64().unwrap(),
            inp["planet_mass"].as_f64().unwrap(),
            inp["G"].as_f64().unwrap(),
        );
        let mr = &inp["mean_rates"];
        let rates = MeanRates {
            n_bar: mr["n_bar"].as_f64().unwrap(),
            lambda_dot: mr["lambda_dot"].as_f64().unwrap(),
            node_dot: mr["node_dot"].as_f64().unwrap(),
        };
        let geo = nodal_geometry(r, v, s, g * (ms + mp));
        let q = case["single_form"]["q_max"].as_u64().unwrap() as usize;
        assert_eq!(q, q_max(geo.e), "{}", case["name"]);
        for (expected, q, worst) in [
            (&case["single_form"], q, &mut worst_form),
            (&case["truth"], kaula_tide::Q_LIMIT, &mut worst_truth),
        ] {
            let t = tide_tabulated(r, &geo, &rates, &love, &spec, g, mp, rs, q);
            let force = axes(&expected["force"]);
            let ef = distance(t.force, force) / force.norm();
            let et = distance(t.torque, axes(&expected["torque"])) / (r.norm() * force.norm());
            *worst = worst.max(ef).max(et);
        }
    }
    println!("stellar tide: Rust vs Python single form {worst_form:.1e}, Rust (q = 42) vs truth {worst_truth:.1e}");
    assert!(worst_form < 1e-12, "{worst_form:e}");
    assert!(worst_truth < 1e-12, "{worst_truth:e}");
}

//------------------------------------------------------------------------------ in Posidonius

const STAR_MASS: f64 = 1.0;
const STAR_RADIUS: f64 = 1.2 * 0.00465047;

fn work_dir() -> PathBuf {
    let dir = PathBuf::from("target/test_kaula_stellar");
    fs::create_dir_all(&dir).unwrap();
    dir
}

// Circular-ish prograde orbit in the plane z = 0 rotated by `inclination` about x, periapsis argument
// `phase`, planet relative to the star
fn planet_state(a: f64, e: f64, inclination: f64, phase: f64, mass: f64) -> (Axes, Axes) {
    let mu = G * (STAR_MASS + mass);
    let p = a * (1. - e * e);
    let (r, vk) = (p / (1. + e), (mu / p).sqrt());
    let (x, y, vx, vy) = (r * phase.cos(), r * phase.sin(), -vk * phase.sin() * (1. + e) / (1. + e), vk * (e + 1.) * phase.cos());
    let (ci, si) = (inclination.cos(), inclination.sin());
    (Axes::from(x, y * ci, y * si), Axes::from(vx, vy * ci, vy * si))
}

// Synthetic signed spectrum, strongly dissipative and asymmetric (Im k2(x) != -Im k2(-x)), Im of the sign of
// x, peak near the main frequency 2 (n - Omega_star) of the integration tests; no spin rescaling
fn dissipative_star_love() -> LoveNumber {
    let tau = 1. / 6.4e-5;
    let (mut x, mut re, mut im) = (vec![], vec![], vec![]);
    for j in 0..1024 {
        let xj = -1e-3 + 2e-3 * (j as f64 + 0.5) / 1024.;
        let w = xj * tau;
        x.push(xj);
        re.push(0.3);
        im.push(0.3 * w / (1. + w * w) * (1. + 0.5 * w.tanh()));
    }
    table(&x, &re, &im, None)
}

fn kaula_star(spin: Axes, love: LoveNumber) -> Particle {
    let zero = Axes::new();
    let mut star = Particle::new(STAR_MASS, STAR_RADIUS, 0.2645, zero, zero, spin);
    star.tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(posidonius::TidalModel::Kaula(
        posidonius::KaulaParameters {
            love_numbers: love,
            ..Default::default()
        },
    )));
    star
}

fn planet(r: Axes, v: Axes, mass: f64, effect: posidonius::TidesEffect) -> Particle {
    let mut planet = Particle::new(mass, 4e-4, 0.5, r, v, Axes::from(0., 0., 1.));
    planet.tides = posidonius::Tides::new(effect);
    planet
}

fn disabled_model() -> posidonius::TidesEffect {
    posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::DisabledModel)
}

fn universe(particles: Vec<Particle>) -> posidonius::Universe {
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    posidonius::Universe::new(0., 1e9, particles, consider_effects)
}

fn host_index(u: &posidonius::Universe) -> usize {
    (0..u.n_particles).find(|&j| u.particles[j].mass == STAR_MASS).unwrap()
}

fn relative(star: &Particle, body: &Particle) -> (Axes, Axes) {
    // r = r_star - r_body, v = v_star - v_body
    let mut r = star.inertial_position;
    r.sub(&body.inertial_position);
    let mut v = star.inertial_velocity;
    v.sub(&body.inertial_velocity);
    (r, v)
}

// Stellar tide of one body at the current state, start-up rates (empty history: n_osc)
fn stellar_tide_alone(star: &Particle, body: &Particle) -> kaula_tide::Tide {
    let (r, v) = relative(star, body);
    let love = kwok_star_love();
    let geo = nodal_geometry(r, v, star.spin, G * (star.mass + body.mass));
    let rates = MeanRates { n_bar: geo.n_osc, lambda_dot: geo.n_osc, node_dot: 0. };
    tide_tabulated(r, &geo, &rates, &love, &prepare_signed_spectrum(&love), G, body.mass, star.radius, q_max(geo.e))
}

#[test]
fn every_body_has_its_own_stellar_tide() {
    // Two planets (OrbitingBody(DisabledModel): stellar tide only), plus a body with
    // TidesEffect::Disabled (no tidal coordinates, no stellar tide). Torque on the star = sum of the
    // torques of each planet alone; force on each planet = -(its own force); power of the star = sum of
    // eq. (8.21)
    let spin = Axes::from(0.05, 0., 2. * PI / 1.9);
    let (r1, v1) = planet_state(0.02, 0.05, 0.1, 0.3, 300. * M_EARTH);
    let (r2, v2) = planet_state(0.045, 0.1, 0.4, 2.0, 20. * M_EARTH);
    let (r3, v3) = planet_state(0.08, 0.0, 0.2, 4.0, 5. * M_EARTH);
    let mut u = universe(vec![
        kaula_star(spin, kwok_star_love()),
        planet(r1, v1, 300. * M_EARTH, disabled_model()),
        planet(r2, v2, 20. * M_EARTH, disabled_model()),
        planet(r3, v3, 5. * M_EARTH, posidonius::TidesEffect::Disabled),
    ]);
    u.calculate_additional_effects(0., false, true, true, posidonius::IgnoreGravityTerms::None);
    u.calculate_denergy_dt();
    let h = host_index(&u);
    let star = u.particles[h];
    let mut torque = Axes::new();
    let mut power = 0.;
    let mut checked = 0;
    for j in (0..u.n_particles).filter(|&j| j != h) {
        let body = u.particles[j];
        let acceleration = body.tides.parameters.output.acceleration;
        if matches!(body.tides.effect, posidonius::TidesEffect::Disabled) {
            assert!(acceleration.norm() == 0., "body without tidal coordinates");
            continue;
        }
        let tide = stellar_tide_alone(&star, &body);
        let mut expected = tide.force;
        expected.mul(-1. / body.mass);
        let err = distance(acceleration, expected) / expected.norm();
        assert!(err < 1e-12, "body {j}: acceleration {err:e}");
        torque.add(&tide.torque);
        let (_, v) = relative(&star, &body);
        power -= tide.dissipative_force.dot(&v) + tide.dissipative_torque.dot(&star.spin);
        checked += 1;
    }
    assert_eq!(checked, 2);
    let measured = star.tides.parameters.output.dangular_momentum_dt;
    let err = distance(measured, torque) / torque.norm();
    println!("torque on the star vs sum of the planets alone {err:.1e}; power {:.3e} vs {power:.3e}", star.tides.parameters.internal.denergy_dt);
    assert!(err < 1e-12, "torque on the star {err:e}");
    let e_power = (star.tides.parameters.internal.denergy_dt - power).abs() / power.abs();
    assert!(e_power < 1e-12, "power dissipated in the star {e_power:e}");
}

fn run(integrator: &mut dyn Integrator, n_steps: usize, name: &str) {
    let path = work_dir().join(format!("{name}_history.bin"));
    let mut writer = posidonius::output::get_universe_history_writer(&path, 0);
    for _ in 0..n_steps {
        let _ = integrator.iterate(&mut writer, true);
    }
}

fn total_angular_momentum(u: &posidonius::Universe) -> (Axes, Axes) {
    let mut orbit = Axes::new();
    let mut spins = Axes::new();
    for particle in &u.particles[..u.n_particles] {
        let mut l = particle.inertial_position.cross(&particle.inertial_velocity);
        l.mul(particle.mass);
        orbit.add(&l);
        spins.add(&particle.angular_momentum);
    }
    (orbit, spins)
}

fn orbital_energy(u: &posidonius::Universe) -> f64 {
    let h = host_index(u);
    let (star, body) = (&u.particles[h], &u.particles[1 - h]);
    let (r, v) = relative(star, body);
    let (m1, m2) = (star.mass, body.mass);
    0.5 * m1 * m2 / (m1 + m2) * v.dot(&v) - G * m1 * m2 / r.norm()
}

// Orbit average of the power F.v and of the torque on the star over the osculating ellipse of the current
// state (mean anomaly uniform, elements and spin fixed), with the mean rates of the pair history
fn orbit_average(u: &posidonius::Universe) -> (f64, Axes) {
    let h = host_index(u);
    let (star, body) = (&u.particles[h], &u.particles[1 - h]);
    let (r0, v0) = relative(star, body);
    let mut spin = star.angular_momentum;
    spin.div(star.moment_of_inertia);
    let mu = G * (star.mass + body.mass);
    let geo0 = nodal_geometry(r0, v0, spin, mu);
    let hv = r0.cross(&v0);
    let mut e_vec = v0.cross(&hv);
    e_vec.div(mu);
    let mut r_hat = r0;
    r_hat.div(r0.norm());
    e_vec.sub(&r_hat);
    let mut p_hat = e_vec;
    p_hat.div(e_vec.norm());
    let mut h_hat = hv;
    h_hat.div(hv.norm());
    let q_hat = h_hat.cross(&p_hat);
    let (a, e) = (geo0.a, geo0.e);
    let love = star.tides.get_kaula().love_numbers;
    let spec = prepare_signed_spectrum(&love);
    let rates = kaula_tide::mean_rates(&body.tides.host_kaula_tide.rate_history, geo0.n_osc);
    let n_grid = 512;
    let (mut power, mut torque) = (0., Axes::new());
    for k in 0..n_grid {
        let m = 2. * PI * k as f64 / n_grid as f64;
        let mut ecc = m;
        for _ in 0..60 {
            ecc -= (ecc - e * ecc.sin() - m) / (1. - e * ecc.cos());
        }
        let (x, y) = (a * (ecc.cos() - e), a * (1. - e * e).sqrt() * ecc.sin());
        let rdot = (mu * a).sqrt() / (a * (1. - e * ecc.cos()));
        let (vx, vy) = (-rdot * ecc.sin(), rdot * (1. - e * e).sqrt() * ecc.cos());
        let comb = |u: f64, w: f64| {
            Axes::from(u * p_hat.x + w * q_hat.x, u * p_hat.y + w * q_hat.y, u * p_hat.z + w * q_hat.z)
        };
        let (r, v) = (comb(x, y), comb(vx, vy));
        let geo = nodal_geometry(r, v, spin, mu);
        let tide = tide_tabulated(r, &geo, &rates, &love, &spec, G, body.mass, star.radius, q_max(geo.e));
        power += tide.force.dot(&v) / n_grid as f64;
        let mut t = tide.torque;
        t.mul(1. / n_grid as f64);
        torque.add(&t);
    }
    (power, torque)
}

#[test]
fn stellar_conservation_and_secular_rates_ias15() {
    check_stellar_conservation_and_secular_rates(0.3, dissipative_star_love(), true);
}

#[test]
fn stellar_conservation_and_secular_rates_ias15_coplanar() {
    check_stellar_conservation_and_secular_rates(0., dissipative_star_love(), true);
}

#[test]
fn stellar_conservation_and_torque_kwok_spectrum() {
    // Kwok et al. spectrum: Re k2 ~ 0.5 and small Im k2, so the periodic energy exchange of the elastic part
    // hides the dissipated energy over the run (no energy check); angular momentum and torque checked
    check_stellar_conservation_and_secular_rates(0.3, kwok_star_love(), false);
    check_stellar_conservation_and_secular_rates(0., kwok_star_love(), false);
}

fn check_stellar_conservation_and_secular_rates(obliquity: f64, love: LoveNumber, check_energy: bool) {
    // Hot Jupiter around a star of spin period 1.9 d (Kwok et al. 2024 "inward" spin), stellar tide only
    // (planet OrbitingBody(DisabledModel)). V9: total angular momentum (orbit + star spin) conserved.
    // Secular rates: measured dE_orb/dt and d(I s_star)/dt = orbit averages of the same model at both ends.
    let mass = 300. * M_EARTH;
    let (a, e) = (0.02, 0.1);
    let (r, v) = planet_state(a, e, 0., 0.4, mass);
    let spin_rate = 2. * PI / 1.9;
    let spin = Axes::from(0., -spin_rate * obliquity.sin(), spin_rate * obliquity.cos());
    let u = universe(vec![kaula_star(spin, love), planet(r, v, mass, disabled_model())]);
    let p_orb = 2. * PI * (a.powi(3) / (G * (STAR_MASS + mass))).sqrt();
    let mut integrator = posidonius::Ias15::new(p_orb / 50., 1e9, 1e9, u);
    integrator.initialize_physical_values();
    let u0 = integrator.universe.clone();
    let (orbit0, spin0) = total_angular_momentum(&u0);
    let t0 = integrator.current_time;
    let n_orbits = 100.;
    // first averages once the history covers a window
    while integrator.current_time < t0 + 10. * p_orb {
        run(&mut integrator, 1, "stellar_secular");
    }
    let (power0, torque0) = orbit_average(&integrator.universe);
    let (energy_a, spin_a, t_a) = (
        orbital_energy(&integrator.universe),
        total_angular_momentum(&integrator.universe).1,
        integrator.current_time,
    );
    while integrator.current_time < t_a + n_orbits * p_orb {
        run(&mut integrator, 1, "stellar_secular");
    }
    let u1 = &integrator.universe;
    let (orbit1, spin1) = total_angular_momentum(u1);
    let mut d_spin = spin1;
    d_spin.sub(&spin0);
    let mut d_total = orbit1;
    d_total.add(&spin1);
    d_total.sub(&orbit0);
    d_total.sub(&spin0);
    let conservation = d_total.norm() / d_spin.norm();
    println!(
        "V9 stellar IAS15 obliquity {obliquity}: |dL_star|/|L_star| = {:.2e}, |dL_tot|/|dL_star| = {conservation:.2e}",
        d_spin.norm() / spin0.norm()
    );
    assert!(conservation < 1e-6, "V9: total angular momentum not conserved: {conservation:e}");

    let (power1, torque1) = orbit_average(u1);
    let elapsed = integrator.current_time - t_a;
    let power_mean = 0.5 * (power0 + power1);
    let mut torque_mean = torque0;
    torque_mean.add(&torque1);
    torque_mean.mul(0.5);
    let energy_rate = (orbital_energy(u1) - energy_a) / elapsed;
    let mut spin_rate_measured = spin1;
    spin_rate_measured.sub(&spin_a);
    spin_rate_measured.mul(1. / elapsed);
    let e_power = (energy_rate - power_mean).abs() / power_mean.abs();
    let e_torque = distance(spin_rate_measured, torque_mean) / torque_mean.norm();
    println!(
        "secular rates, stellar: dE/dt = {energy_rate:.6e} vs <F.v> = {power_mean:.6e} ({e_power:.1e}); d(I s*)/dt vs <Gamma>: {e_torque:.1e}"
    );
    assert!(!check_energy || e_power < 1e-3, "secular energy rate: {e_power:e}");
    assert!(e_torque < 1e-3, "secular spin rate: {e_torque:e}");
    if obliquity == 0. {
        let h = host_index(u1);
        let (r, v) = relative(&u1.particles[h], &u1.particles[1 - h]);
        let geo = nodal_geometry(r, v, u1.particles[h].spin, G * (STAR_MASS + mass));
        assert!(geo.inclination <= kaula_tide::COPLANAR_INCLINATION, "i = {:e}", geo.inclination);
    }
}
