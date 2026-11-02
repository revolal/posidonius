// Integration tests of the planetary Kaula tide (docs/tides/kaula_tides.pdf, Sections 8.10 and 10.9):
//   V9               conservation of the total angular momentum (orbit + spins)
//   secular rates    rates of the orbital energy and of the spin = orbit averages of the same model
//   V3               constant-time-lag limit: same instantaneous force and torque as the ConstantTimeLag model
//   restart          restart from a recovery snapshot (.json and .bin) = uninterrupted run
//   dissipated power denergy_dt of the Kaula branch = -(F_Im . v + Gamma_Im . s), (8.21)
// V9 (IAS15, WHFast), the secular rates, V3 and the restart are run at obliquity 0.3 (single form) and 0
// (coplanar form (8.19)).
use posidonius::constants::{AU, G, M_EARTH, TWO_PI};
use posidonius::kaula_tide::{self, MeanRates};
use posidonius::{Axes, Integrator, LoveNumber, Particle};
use std::fs;
use std::path::{Path, PathBuf};

const DAY: f64 = 86400.;
const R_EARTH: f64 = 6.3781e6 / AU;
const STAR_MASS: f64 = 0.09;

//------------------------------------------------------------------------------ helpers

fn work_dir() -> PathBuf {
    let dir = PathBuf::from("target/test_kaula_planet_integration");
    fs::create_dir_all(&dir).unwrap();
    dir
}

// Keplerian elements -> heliocentric state (planet relative to the star)
fn kepler_state(a: f64, e: f64, i: f64, node: f64, omega: f64, mean_anomaly: f64, mu: f64) -> (Axes, Axes) {
    let mut ecc = mean_anomaly;
    for _ in 0..100 {
        ecc -= (ecc - e * ecc.sin() - mean_anomaly) / (1. - e * ecc.cos());
    }
    let f = 2. * ((1. + e).sqrt() * (ecc / 2.).sin()).atan2((1. - e).sqrt() * (ecc / 2.).cos());
    let p = a * (1. - e * e);
    let r = p / (1. + e * f.cos());
    let (x, y) = (r * f.cos(), r * f.sin());
    let vk = (mu / p).sqrt();
    let (vx, vy) = (-vk * f.sin(), vk * (e + f.cos()));
    let (co, so, ci, si, cw, sw) = (node.cos(), node.sin(), i.cos(), i.sin(), omega.cos(), omega.sin());
    let rot = [
        [co * cw - so * sw * ci, -co * sw - so * cw * ci],
        [so * cw + co * sw * ci, -so * sw + co * cw * ci],
        [sw * si, cw * si],
    ];
    let apply = |u: f64, w: f64| {
        Axes::from(
            rot[0][0] * u + rot[0][1] * w,
            rot[1][0] * u + rot[1][1] * w,
            rot[2][0] * u + rot[2][1] * w,
        )
    };
    (apply(x, y), apply(vx, vy))
}

// Spin of rate `rate` with obliquity `obliquity` with respect to the orbit normal of (r, v)
fn spin_vector(r: Axes, v: Axes, rate: f64, obliquity: f64) -> Axes {
    let h = r.cross(&v);
    let mut h_hat = h;
    h_hat.div(h.norm());
    let mut t = h_hat.cross(&Axes::from(0.3, 1., 0.2));
    t.div(t.norm());
    Axes::from(
        rate * (obliquity.cos() * h_hat.x + obliquity.sin() * t.x),
        rate * (obliquity.cos() * h_hat.y + obliquity.sin() * t.y),
        rate * (obliquity.cos() * h_hat.z + obliquity.sin() * t.z),
    )
}

// Tabulated spectrum, chi in rad/s, 1000 points from 0 to chi_max (log-spaced above 1e-12), zero padding
fn tabulated(love: impl Fn(f64) -> (f64, f64), chi_max: f64) -> LoveNumber {
    let mut x = [0.; 1024];
    let mut re = [0.; 1024];
    let mut im = [0.; 1024];
    let (re0, _) = love(0.);
    re[0] = re0;
    for j in 1..1000 {
        let chi = 10f64.powf(-12. + (chi_max.log10() + 12.) * (j - 1) as f64 / 998.);
        let (r, i) = love(chi);
        x[j] = chi;
        re[j] = r;
        im[j] = i;
    }
    LoveNumber::new_from(x, re, im, None)
}

fn maxwell_spectrum(k2: f64, tau: f64) -> LoveNumber {
    tabulated(|chi| {
        let x = chi * tau;
        (k2 / (1. + x * x), k2 * x / (1. + x * x))
    }, 1.)
}

fn constant_time_lag_spectrum(k2: f64, dt_s: f64) -> LoveNumber {
    // chi dt < pi over the table, so that imk2 >= 0
    tabulated(|chi| (k2 * (chi * dt_s).cos(), k2 * (chi * dt_s).sin()), 1. / dt_s)
}

struct Setup {
    a: f64,
    e: f64,
    obliquity: f64,
    spin_over_n: f64,
    radius: f64,
    mass: f64,
}

impl Default for Setup {
    fn default() -> Self {
        Setup {
            a: 0.0115,
            e: 0.1,
            obliquity: 0.3,
            spin_over_n: 1.7,
            radius: 3. * R_EARTH,
            mass: M_EARTH,
        }
    }
}

fn particles(setup: &Setup, planet_tides: posidonius::Tides, star_tides: posidonius::Tides) -> Vec<Particle> {
    let mu = G * (STAR_MASS + setup.mass);
    let (r, v) = kepler_state(setup.a, setup.e, 0.2, 0.7, 1.1, 0.4, mu);
    let n = (mu / setup.a.powi(3)).sqrt();
    let zero = Axes::new();
    let mut star = Particle::new(STAR_MASS, 0.117 * 0.00465047, 0.447, zero, zero, Axes::from(0., 0., TWO_PI));
    star.tides = star_tides;
    let mut planet = Particle::new(
        setup.mass,
        setup.radius,
        0.339f64.sqrt(),
        r,
        v,
        spin_vector(r, v, setup.spin_over_n * n, setup.obliquity),
    );
    planet.tides = planet_tides;
    vec![star, planet]
}

fn kaula_tides(love: LoveNumber) -> posidonius::Tides {
    posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::Kaula(
        posidonius::KaulaParameters {
            love_numbers: love,
            ..Default::default()
        },
    )))
}

fn disabled_star_tides() -> posidonius::Tides {
    posidonius::Tides::new(posidonius::TidesEffect::CentralBody(posidonius::TidalModel::DisabledModel))
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

fn run(integrator: &mut dyn Integrator, n_steps: usize, name: &str) {
    let path = work_dir().join(format!("{name}_history.bin"));
    let mut writer = posidonius::output::get_universe_history_writer(&path, 0);
    for _ in 0..n_steps {
        let _ = integrator.iterate(&mut writer, true);
    }
}

fn universe_of(integrator: &dyn Integrator) -> &posidonius::Universe {
    let any = integrator.as_any();
    if let Some(i) = any.downcast_ref::<posidonius::Ias15>() {
        &i.universe
    } else if let Some(i) = any.downcast_ref::<posidonius::WHFast>() {
        &i.universe
    } else {
        &any.downcast_ref::<posidonius::LeapFrog>().unwrap().universe
    }
}

fn relative_state(u: &posidonius::Universe) -> (Axes, Axes) {
    let (star, planet) = (&u.particles[0], &u.particles[1]);
    let mut r = planet.inertial_position;
    r.sub(&star.inertial_position);
    let mut v = planet.inertial_velocity;
    v.sub(&star.inertial_velocity);
    (r, v)
}

// The planet stays in the domain of the coplanar form (i <= COPLANAR_INCLINATION) during a coplanar run
fn assert_coplanar(u: &posidonius::Universe) {
    let (r, v) = relative_state(u);
    let geo = kaula_tide::nodal_geometry(r, v, u.particles[1].spin, G * (u.particles[0].mass + u.particles[1].mass));
    assert!(geo.inclination <= kaula_tide::COPLANAR_INCLINATION, "i = {:e}", geo.inclination);
}

// Orbital angular momentum about the barycentre plus the spins
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
    let (r, v) = relative_state(u);
    let (m1, m2) = (u.particles[0].mass, u.particles[1].mass);
    0.5 * m1 * m2 / (m1 + m2) * v.dot(&v) - G * m1 * m2 / r.norm()
}

// Orbit average of the tidal power F.v and of the torque over the osculating ellipse of the current
// state (mean anomaly uniform, elements and spin fixed, Kepler rates), as (8.22)
fn orbit_average(u: &posidonius::Universe) -> (f64, Axes) {
    let (r0, v0) = relative_state(u);
    let planet = &u.particles[1];
    let mut spin = planet.angular_momentum;
    spin.div(planet.moment_of_inertia);
    let mu = G * (STAR_MASS + planet.mass);
    let geo0 = kaula_tide::nodal_geometry(r0, v0, spin, mu);
    // elements of the planet relative to the star in inertial axes
    let h = r0.cross(&v0);
    let mut e_vec = v0.cross(&h);
    e_vec.div(mu);
    let mut r_hat = r0;
    r_hat.div(r0.norm());
    e_vec.sub(&r_hat);
    let mut p_hat = e_vec;
    p_hat.div(e_vec.norm());
    let mut h_hat = h;
    h_hat.div(h.norm());
    let q_hat = h_hat.cross(&p_hat);
    let (a, e) = (geo0.a, geo0.e);
    let kaula = planet.tides.get_kaula();
    let love = kaula.love_numbers;
    let spec = kaula_tide::prepare_spectrum(&love);
    // mean rates of the run (history of the accepted steps), as in the force evaluations
    let rates: MeanRates = kaula_tide::mean_rates(&kaula.rate_history, geo0.n_osc);
    let n_grid = 512;
    let (mut power, mut torque) = (0., Axes::new());
    for k in 0..n_grid {
        let m = TWO_PI * k as f64 / n_grid as f64;
        let mut ecc = m;
        for _ in 0..60 {
            ecc -= (ecc - e * ecc.sin() - m) / (1. - e * ecc.cos());
        }
        let (x, y) = (a * (ecc.cos() - e), a * (1. - e * e).sqrt() * ecc.sin());
        let rdot = (mu * a).sqrt() / (a * (1. - e * ecc.cos()));
        let (vx, vy) = (-rdot * ecc.sin(), rdot * (1. - e * e).sqrt() * ecc.cos());
        let comb = |u: f64, w: f64| {
            Axes::from(
                u * p_hat.x + w * q_hat.x,
                u * p_hat.y + w * q_hat.y,
                u * p_hat.z + w * q_hat.z,
            )
        };
        let (r, v) = (comb(x, y), comb(vx, vy));
        let geo = kaula_tide::nodal_geometry(r, v, spin, mu);
        let tide = kaula_tide::tide_tabulated(
            r,
            &geo,
            &rates,
            &love,
            &spec,
            G,
            STAR_MASS,
            planet.radius,
            kaula_tide::q_max(geo.e),
        );
        power += tide.force.dot(&v) / n_grid as f64;
        let mut t = tide.torque;
        t.mul(1. / n_grid as f64);
        torque.add(&t);
    }
    (power, torque)
}

fn ias15(setup: &Setup, love: LoveNumber, time_step: f64) -> posidonius::Ias15 {
    let u = universe(particles(setup, kaula_tides(love), disabled_star_tides()));
    let mut integrator = posidonius::Ias15::new(time_step, 1e9, 1e9, u);
    integrator.initialize_physical_values();
    integrator
}

fn period(setup: &Setup) -> f64 {
    TWO_PI * (setup.a.powi(3) / (G * (STAR_MASS + setup.mass))).sqrt()
}

//------------------------------------------------------------------------------ V9 and secular rates

#[test]
fn conservation_and_secular_rates_ias15() {
    check_conservation_and_secular_rates_ias15(0.3);
}

#[test]
fn conservation_and_secular_rates_ias15_coplanar() {
    check_conservation_and_secular_rates_ias15(0.);
}

fn check_conservation_and_secular_rates_ias15(obliquity: f64) {
    // Weak tide (a = 0.03 AU, R = 1 R_earth): the spin changes by about 1 % over the run, so that the
    // orbit averages at both ends bracket the measured rates
    let setup = Setup {
        a: 0.03,
        radius: R_EARTH,
        obliquity,
        ..Default::default()
    };
    let p_orb = period(&setup);
    let love = maxwell_spectrum(0.5, 1. / (2. * TWO_PI / (p_orb * DAY)));
    let mut integrator = ias15(&setup, love, p_orb / 50.);
    let u0 = integrator.universe.clone();
    let (orbit0, spin0) = total_angular_momentum(&u0);
    let energy0 = orbital_energy(&u0);
    let (power0, torque0) = orbit_average(&u0);
    let t0 = integrator.current_time;
    let n_orbits = 100.;
    while integrator.current_time < t0 + n_orbits * p_orb {
        run(&mut integrator, 1, "secular");
    }
    let u1 = &integrator.universe;
    let elapsed = integrator.current_time - t0;
    let (orbit1, spin1) = total_angular_momentum(u1);
    let mut total0 = orbit0;
    total0.add(&spin0);
    let mut total1 = orbit1;
    total1.add(&spin1);
    let mut d_spin = spin1;
    d_spin.sub(&spin0);
    let d_total = {
        let mut d = total1;
        d.sub(&total0);
        d
    };
    let conservation = d_total.norm() / d_spin.norm();
    println!(
        "V9 IAS15 obliquity {obliquity}, {n_orbits} orbits: |dL_spin|/|L_tot| = {:.2e}, |dL_tot|/|dL_spin| = {conservation:.2e}",
        d_spin.norm() / total0.norm()
    );
    assert!(conservation < 1e-6, "V9: total angular momentum not conserved: {conservation:e}");

    // Secular rates: rates measured over the run against the trapezoidal mean of the orbit averages at both ends
    let (power1, torque1) = orbit_average(u1);
    let power_mean = 0.5 * (power0 + power1);
    let mut torque_mean = torque0;
    torque_mean.add(&torque1);
    torque_mean.mul(0.5);
    let energy_rate = (orbital_energy(u1) - energy0) / elapsed;
    let mut spin_rate = d_spin;
    spin_rate.mul(1. / elapsed);
    let e_power = (energy_rate - power_mean).abs() / power_mean.abs();
    let e_torque = distance(spin_rate, torque_mean) / torque_mean.norm();
    println!(
        "secular rates: dE/dt = {energy_rate:.6e} vs <F.v> = {power_mean:.6e} ({e_power:.1e}); d(I s)/dt vs <Gamma>: {e_torque:.1e}; change of <F.v> over the run {:.1e}",
        (power1 - power0).abs() / power0.abs()
    );
    assert!(e_power < 1e-3, "secular energy rate: {e_power:e}");
    assert!(e_torque < 1e-3, "secular spin rate: {e_torque:e}");
    if obliquity == 0. {
        assert_coplanar(u1);
    }
}

fn distance(a: Axes, b: Axes) -> f64 {
    let mut d = a;
    d.sub(&b);
    d.norm()
}

#[test]
fn conservation_whfast_and_leapfrog() {
    let setup = Setup::default();
    let p_orb = period(&setup);
    let love = maxwell_spectrum(0.5, 1. / (2. * TWO_PI / (p_orb * DAY)));
    let mut leapfrog_errors = vec![];
    for (name, steps_per_orbit) in [("whfast", 100usize), ("leapfrog", 100), ("leapfrog", 400)] {
        let u = universe(particles(&setup, kaula_tides(love), disabled_star_tides()));
        let (orbit0, spin0) = total_angular_momentum(&u);
        let time_step = p_orb / steps_per_orbit as f64;
        let mut integrator: Box<dyn Integrator> = if name == "whfast" {
            Box::new(posidonius::WHFast::new(
                time_step,
                1e9,
                1e9,
                u,
                posidonius::whfast::CoordinatesType::DemocraticHeliocentric,
            ))
        } else {
            Box::new(posidonius::LeapFrog::new(time_step, 1e9, 1e9, u))
        };
        integrator.initialize_physical_values();
        // 5 orbits: the strong tide of this setup brings the spin near equilibrium in about 10 orbits
        run(integrator.as_mut(), 5 * steps_per_orbit, name);
        let (orbit1, spin1) = total_angular_momentum(universe_of(integrator.as_ref()));
        let mut d_spin = spin1;
        d_spin.sub(&spin0);
        let mut d_total = orbit1;
        d_total.add(&spin1);
        d_total.sub(&orbit0);
        d_total.sub(&spin0);
        let ratio = d_total.norm() / d_spin.norm();
        println!("V9 {name} dt = P/{steps_per_orbit}, 5 orbits: |dL_tot|/|dL_spin| = {ratio:.2e}");
        if name == "whfast" {
            assert!(ratio < 1e-6, "V9 {name}: {ratio:e}");
        } else {
            leapfrog_errors.push(ratio);
        }
    }
    // LeapFrog kicks the spin with the torques at the start and at the middle of the step, the orbit
    // with the torque at the middle only, evaluated at (r_{n+1/2}, v_n): first-order error in dt for
    // a force that depends on the velocity through the elements
    let order = (leapfrog_errors[0] / leapfrog_errors[1]).log(4.);
    println!("V9 leapfrog: order of the angular-momentum error in dt = {order:.2}");
    assert!(order > 0.8, "LeapFrog error does not decrease with dt: order {order}");
}

#[test]
fn conservation_whfast_coplanar() {
    let setup = Setup {
        obliquity: 0.,
        ..Default::default()
    };
    let p_orb = period(&setup);
    let love = maxwell_spectrum(0.5, 1. / (2. * TWO_PI / (p_orb * DAY)));
    let u = universe(particles(&setup, kaula_tides(love), disabled_star_tides()));
    let (orbit0, spin0) = total_angular_momentum(&u);
    let mut integrator = posidonius::WHFast::new(
        p_orb / 100.,
        1e9,
        1e9,
        u,
        posidonius::whfast::CoordinatesType::DemocraticHeliocentric,
    );
    integrator.initialize_physical_values();
    run(&mut integrator, 500, "whfast_coplanar");
    let (orbit1, spin1) = total_angular_momentum(universe_of(&integrator));
    let mut d_spin = spin1;
    d_spin.sub(&spin0);
    let mut d_total = orbit1;
    d_total.add(&spin1);
    d_total.sub(&orbit0);
    d_total.sub(&spin0);
    let ratio = d_total.norm() / d_spin.norm();
    println!("V9 whfast coplanar dt = P/100, 5 orbits: |dL_tot|/|dL_spin| = {ratio:.2e}");
    assert!(ratio < 1e-6, "V9 whfast coplanar: {ratio:e}");
    assert_coplanar(universe_of(&integrator));
}

//------------------------------------------------------------------------------ V3

#[test]
fn constant_time_lag_limit() {
    check_constant_time_lag_limit(0.3);
}

#[test]
fn constant_time_lag_limit_coplanar() {
    check_constant_time_lag_limit(0.);
}

fn check_constant_time_lag_limit(obliquity: f64) {
    // bar k_2 = k2 e^{-i sigma dt} tabulated, against the ConstantTimeLag model of Posidonius with
    // k2 dt = 1.5 sigma_p R^5 / G. e = 0.05 (q_max = 5): the truncation error of the elastic part (bounded
    // by (8.15)) is then negligible against the torque, which is O(sigma dt); at e = 0.01 (q_max = 3) it
    // gives a torque of 1.6e-7 r |F| independent of dt
    let k2 = 0.3;
    let mut errors = vec![];
    let mut powers = vec![];
    for dt_s in [100., 10.] {
        let setup = Setup {
            e: 0.05,
            obliquity,
            ..Default::default()
        };
        let mut kaula = universe(particles(&setup, kaula_tides(constant_time_lag_spectrum(k2, dt_s)), disabled_star_tides()));
        let dissipation_factor = 2. * G * k2 * (dt_s / DAY) / (3. * setup.radius.powi(5));
        let ctl = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(
            posidonius::TidalModel::ConstantTimeLag(posidonius::ConstantTimeLagParameters {
                dissipation_factor,
                dissipation_factor_scale: 1.,
                love_number: k2,
            }),
        ));
        let ctl_star = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
            posidonius::TidalModel::ConstantTimeLag(posidonius::ConstantTimeLagParameters {
                dissipation_factor: 0.,
                dissipation_factor_scale: 1.,
                love_number: 0.,
            }),
        ));
        let mut reference = universe(particles(&setup, ctl, ctl_star));
        let mut outputs = vec![];
        for u in [&mut kaula, &mut reference] {
            u.inertial_to_heliocentric();
            u.calculate_additional_effects(0., false, true, true, posidonius::IgnoreGravityTerms::None);
            u.calculate_denergy_dt();
            let p = &u.particles[1];
            let mut force = p.tides.parameters.output.acceleration;
            force.mul(p.mass);
            outputs.push((force, p.tides.parameters.output.dangular_momentum_dt, p.tides.parameters.internal.denergy_dt));
        }
        let (r, _) = relative_state(&kaula);
        let (f_k, g_k, p_k) = outputs[0];
        let (f_c, g_c, p_c) = outputs[1];
        // (8.21): same dissipated power as denergy_dt of the CTL model (Bolmont et al. 2015 eq. 32) to
        // O((n dt)^2) (0.93 (n dt)^2 in the Python reference)
        let e_power = (p_k - p_c).abs() / p_c.abs();
        if obliquity == 0. {
            assert_coplanar(&kaula);
        }
        println!("V3 obliquity {obliquity}, dt = {dt_s} s: denergy_dt Kaula {p_k:.6e}, CTL {p_c:.6e}, relative difference {e_power:.2e}");
        powers.push(e_power);
        let e_force = distance(f_k, f_c) / f_c.norm();
        let e_torque = distance(g_k, g_c) / g_c.norm();
        let n = (G * (STAR_MASS + setup.mass) / setup.a.powi(3)).sqrt() / DAY;
        println!(
            "V3 dt = {dt_s} s (n dt = {:.1e}): force {e_force:.2e}, torque {e_torque:.2e} (|Gamma|/(r|F|) = {:.1e})",
            n * dt_s,
            g_c.norm() / (r.norm() * f_c.norm())
        );
        // Re k2 = k2 cos(sigma dt) = k2 (1 - (sigma dt)^2 / 2): O((sigma dt)^2) on the force, with
        // sigma dt ~ 7e-3 at dt = 100 s; O(sigma dt) relative on the torque, whose leading term is the
        // same in both models
        assert!(e_force < 3e-5 * (dt_s / 100.).powi(2) + 1e-9, "V3 force: {e_force:e}");
        assert!(e_torque < 2e-3 * (dt_s / 100.), "V3 torque: {e_torque:e}");
        errors.push((e_force, e_torque));
    }
    let order_force = (errors[0].0 / errors[1].0).log10();
    let order_torque = (errors[0].1 / errors[1].1).log10();
    println!("V3: order in dt of the difference: force {order_force:.2}, torque {order_torque:.2}");
    assert!(order_force > 1.8 && order_torque > 0.8);
    let order_power = (powers[0] / powers[1]).log10();
    println!("V3: order in dt of the difference of denergy_dt: {order_power:.2}");
    assert!(powers[0] < 1e-3 && order_power > 0.8);
}

//------------------------------------------------------------------------------ restart

#[test]
fn restart_from_recovery_snapshot() {
    check_restart_from_recovery_snapshot(0.3);
}

#[test]
fn restart_from_recovery_snapshot_coplanar() {
    check_restart_from_recovery_snapshot(0.);
}

fn check_restart_from_recovery_snapshot(obliquity: f64) {
    // .bin and .json snapshots: bit-identical continuation, history included (IAS15 and WHFast). The .json
    // case needs the feature float_roundtrip of serde_json (Cargo.toml; without it about 20 % of the floats
    // are parsed 1 ulp off).
    let setup = Setup {
        obliquity,
        ..Default::default()
    };
    let p_orb = period(&setup);
    let love = maxwell_spectrum(0.5, 1. / (2. * TWO_PI / (p_orb * DAY)));
    for (name, extension, tolerance) in [("ias15", "bin", 0.), ("whfast", "bin", 0.), ("ias15", "json", 0.), ("whfast", "json", 0.)] {
        let make = || -> Box<dyn Integrator> {
            let u = universe(particles(&setup, kaula_tides(love), disabled_star_tides()));
            let mut integrator: Box<dyn Integrator> = if name == "ias15" {
                Box::new(posidonius::Ias15::new(p_orb / 50., 1e9, 1e9, u))
            } else {
                Box::new(posidonius::WHFast::new(
                    p_orb / 50.,
                    1e9,
                    1e9,
                    u,
                    posidonius::whfast::CoordinatesType::DemocraticHeliocentric,
                ))
            };
            integrator.initialize_physical_values();
            integrator
        };
        let mut straight = make();
        run(straight.as_mut(), 300, "restart_straight");
        let mut first = make();
        run(first.as_mut(), 150, "restart_first");
        let path: PathBuf = work_dir().join(format!("restart_{name}_{obliquity}.{extension}"));
        let history = work_dir().join(format!("restart_history_{obliquity}.bin"));
        let mut writer = posidonius::output::get_universe_history_writer(&history, 0);
        first.write_recovery_snapshot(Path::new(&path), &mut writer);
        let mut resumed = posidonius::output::restore_snapshot(&path).unwrap();
        run(resumed.as_mut(), 150, "restart_resumed");
        let (a, b) = (universe_of(straight.as_ref()), universe_of(resumed.as_ref()));
        let (ra, va) = relative_state(a);
        let (rb, vb) = relative_state(b);
        let dr = distance(ra, rb) / ra.norm();
        let dv = distance(va, vb) / va.norm();
        let ds = distance(a.particles[1].angular_momentum, b.particles[1].angular_momentum)
            / a.particles[1].angular_momentum.norm();
        println!("restart {name} .{extension}, obliquity {obliquity}: position {dr:.1e}, velocity {dv:.1e}, spin {ds:.1e}");
        assert!(dr <= tolerance && dv <= tolerance && ds <= tolerance, "restart {name} {extension}: {dr:e} {dv:e} {ds:e}");
        if tolerance == 0. {
            assert_eq!(
                a.particles[1].tides.get_kaula().rate_history,
                b.particles[1].tides.get_kaula().rate_history
            );
        }
    }
}

//------------------------------------------------------------------------------ dissipated power

#[test]
fn denergy_dt_kaula_branch() {
    let setup = Setup::default();
    let p_orb = period(&setup);
    let love = maxwell_spectrum(0.5, 1. / (2. * TWO_PI / (p_orb * DAY)));
    let mut integrator = ias15(&setup, love, p_orb / 50.);
    run(&mut integrator, 10, "denergy");
    let u = &mut integrator.universe;
    u.calculate_denergy_dt();
    let p = &u.particles[1];
    let state = p.tides.get_kaula().tide_state;
    assert!(state.valid);
    let expected = -(state.tide.dissipative_force.dot(&state.velocity) + state.tide.dissipative_torque.dot(&state.spin));
    assert_eq!(p.tides.parameters.internal.denergy_dt, expected);
    assert!(expected > 0.);
    // elastic planet: no dissipated power
    let elastic = tabulated(|_| (0.3, 0.), 1.);
    let mut integrator = ias15(&setup, elastic, p_orb / 50.);
    run(&mut integrator, 3, "denergy_elastic");
    integrator.universe.calculate_denergy_dt();
    assert_eq!(integrator.universe.particles[1].tides.parameters.internal.denergy_dt, 0.);
}
