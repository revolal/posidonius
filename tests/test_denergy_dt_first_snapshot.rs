// At the first historic snapshot the tidal coordinates are not computed yet (distance = 0) and denergy_dt of
// the non-Kaula orbiting bodies was 0/0 = NaN; a DisabledModel body went through the same formula. serde writes NaN as null and the .json
// recovery snapshot could not be read back ("Unknown integrator!"). Now denergy_dt = 0 in both cases.
use posidonius::constants::{AU, G, M_EARTH, TWO_PI};
use posidonius::{Axes, Integrator, Particle, TidalModel, TidesEffect};
use std::fs;
use std::path::PathBuf;

const R_EARTH: f64 = 6.3781e6 / AU;
const STAR_MASS: f64 = 0.09;

fn ctl(dissipation_factor: f64) -> TidalModel {
    TidalModel::ConstantTimeLag(posidonius::ConstantTimeLagParameters {
        dissipation_factor,
        dissipation_factor_scale: 1.,
        love_number: 0.3,
    })
}

fn integrator(star: TidesEffect, planets: Vec<TidesEffect>) -> Box<posidonius::Ias15> {
    let zero = Axes::new();
    let mut star_particle = Particle::new(STAR_MASS, 0.117 * 0.00465047, 0.447, zero, zero, Axes::from(0., 0., TWO_PI));
    star_particle.tides = posidonius::Tides::new(star);
    let mut particles = vec![star_particle];
    for (j, tides) in planets.into_iter().enumerate() {
        let a = 0.0115 * (1. + j as f64);
        let v = (G * (STAR_MASS + M_EARTH) / a).sqrt();
        let mut planet = Particle::new(M_EARTH, R_EARTH, 0.5, Axes::from(a, 0., 0.), Axes::from(0., v, 0.), Axes::from(0., 0., 2.));
        planet.tides = posidonius::Tides::new(tides);
        particles.push(planet);
    }
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    // boxed: the integrator does not fit on the default test thread stack
    Box::new(posidonius::Ias15::new(0.01, 1e9, 1e9, posidonius::Universe::new(0., 1e9, particles, consider_effects)))
}

#[test]
fn denergy_dt_finite_at_first_snapshot() {
    let dir = PathBuf::from("target/test_denergy_dt_first_snapshot");
    fs::create_dir_all(&dir).unwrap();
    let cases = [
        ("disabled_planet", TidesEffect::CentralBody(ctl(1e-3)), vec![TidesEffect::OrbitingBody(ctl(1e-3)), TidesEffect::OrbitingBody(TidalModel::DisabledModel)]),
        ("star_without_dissipation", TidesEffect::CentralBody(ctl(0.)), vec![TidesEffect::OrbitingBody(ctl(1e-3))]),
        ("disabled_star", TidesEffect::CentralBody(TidalModel::DisabledModel), vec![TidesEffect::OrbitingBody(ctl(1e-3))]),
    ];
    for (name, star, planets) in cases {
        let mut integrator = integrator(star, planets);
        integrator.initialize_physical_values();
        let mut writer = posidonius::output::get_universe_history_writer(&dir.join(format!("{name}_history.bin")), 0);
        // first iteration: first historic snapshot (denergy_dt before any force evaluation), then one step
        integrator.iterate(&mut writer, true).unwrap();
        let u = &integrator.universe;
        for particle in &u.particles[..u.n_particles] {
            let p = particle.tides.parameters.internal.denergy_dt;
            println!("first snapshot, {name}: denergy_dt {p:e}");
            assert!(p.is_finite(), "{name}: denergy_dt = {p}");
            if let TidesEffect::OrbitingBody(TidalModel::DisabledModel) = particle.tides.effect {
                assert_eq!(p, 0.);
            }
        }
        // the .json recovery snapshot written at this point is readable
        let path = dir.join(format!("{name}.json"));
        integrator.write_recovery_snapshot(&path, &mut writer);
        // restore in a thread with a larger stack (restoring overflows the default 2 MiB test stack, as
        // planet_tides_rust_vs_python of the original code)
        let readable = std::thread::Builder::new()
            .stack_size(256 << 20)
            .spawn(move || posidonius::output::restore_snapshot(&path).is_ok())
            .unwrap()
            .join()
            .unwrap();
        assert!(readable, "{name}: .json snapshot not readable");
        // later snapshots: the CTL value is computed again (positive dissipated power of the CTL planet)
        integrator.universe.calculate_denergy_dt();
        let p = integrator.universe.particles[1].tides.parameters.internal.denergy_dt;
        assert!(p.is_finite() && p > 0., "{name}: CTL denergy_dt after one step = {p}");
    }
}
