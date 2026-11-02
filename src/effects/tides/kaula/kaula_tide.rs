// Tide raised in a deformed body D by a perturber P (Kaula, degree 2), single form valid for any
// obliquity, and coplanar form. Equation numbers refer to docs/tides/kaula_tides.pdf (Chapter 8),
// where the equations are written for D = planet and P = star: the names below (r_star, c_star, ...)
// follow that notation. The same functions give the stellar tide with D = star and P = planet
// (Chapter 9). All quantities are computed from the state vectors r = r_D - r_P, v = v_D - v_P and the
// spin s of D, in Posidonius units (AU, day, M_sun, rad/day), inertial axes.
use super::LoveNumber;
use crate::Axes;
use crate::constants::{DAY, PI, TWO_PI};
use serde::{Deserialize, Serialize};
use serde_big_array::BigArray;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicBool, Ordering};
use time::OffsetDateTime;
use time::format_description;

// Truncation |q| <= q_max(e): table of Section 8.5, guarantee (8.15) with tol 1e-5
pub const Q_MAX_ECCENTRICITY_EDGES: [f64; 12] = [
    0.001, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50,
];
pub const Q_MAX_TABLE: [usize; 12] = [2, 3, 3, 5, 7, 10, 12, 15, 18, 23, 28, 42];
// Number of terms of Tisserand's series for an absolute error below 1e-13 on G_2pq (Appendix E)
const TISSERAND_TERMS: [usize; 12] = [3, 4, 4, 5, 6, 7, 8, 9, 11, 12, 14, 22];
// Largest q_max of the table: arrays over q are indexed by q + Q_LIMIT
pub const Q_LIMIT: usize = 42;
const N_Q: usize = 2 * Q_LIMIT + 1;
// |s^ x h^| below which the node is undefined and any equatorial direction is taken (8.2)
pub const EPSILON_NODE: f64 = 1e-10;
// Obliquity up to which the coplanar form (8.19) is used (Section 8.7): its error is below C(e) i F0
// with C <= 55 for e <= 0.5 and any |bar k| <= kappa, eq. (8.20), i.e. below the truncation tolerance
// 1e-5 F0 at the switch
pub const COPLANAR_INCLINATION: f64 = 1e-7;
// Mean-rate history: samples stored every window / HISTORY_SAMPLES_PER_WINDOW, ring buffer of
// HISTORY_CAPACITY samples (Sections 8.4.1 and 10.4)
pub const HISTORY_SAMPLES_PER_WINDOW: f64 = 32.;
pub const HISTORY_CAPACITY: usize = 40;

// Chebyshev tables of G_2pq for e <= G_TABLE_ECCENTRICITY, eq. (8.16): H_pq(u) = G_2pq / beta^|q|,
// u = beta^2, interpolated on G_TABLE_INTERVALS equal intervals of u by Chebyshev polynomials of degree
// G_TABLE_DEGREE; absolute error 3e-15 on G
const G_TABLE_ECCENTRICITY: f64 = 0.5;
const G_TABLE_INTERVALS: usize = 8;
const G_TABLE_DEGREE: usize = 8;
// Terms of Tisserand's series at the nodes: the terms grow with |q| before beta^(2k) wins, so the count of
// TISSERAND_TERMS (valid for |q| <= q_max(e)) is not enough for all |q| <= 42; 45 terms
const G_TABLE_TISSERAND_TERMS: usize = 45;

// c_m = (2-m)!/(2+m)! (2 - delta_0m) of eq. (8.6)
const C_M: [f64; 3] = [1.0, 1.0 / 3.0, 1.0 / 12.0];

static WARNED_ECCENTRICITY_OUTSIDE_TABLE: AtomicBool = AtomicBool::new(false);
static WARNED_FREQUENCY_OUTSIDE_SPECTRUM: AtomicBool = AtomicBool::new(false);
static WARNED_FREQUENCY_OUTSIDE_SIGNED_SPECTRUM: AtomicBool = AtomicBool::new(false);

//------------------------------------------------------------------------------
// Data structures

/// Geometry of the perturber seen from the deformed body (eqs. (8.1) to (8.5)).
#[derive(Debug, Copy, Clone, Default, PartialEq)]
pub struct NodalGeometry {
    pub n_hat: Axes,       // N^, node of the star's orbit on the planet's equator (8.2)
    pub y_hat: Axes,       // Y^ = s^ x N^ (8.2)
    pub s_hat: Axes,       // spin axis s^
    pub spin_rate: f64,    // Omega_p = |s|, rad/day
    pub inclination: f64,  // i in [0, pi] (8.3)
    pub a: f64,            // semi-major axis, AU (8.1)
    pub e: f64,            // eccentricity (8.1)
    pub n_osc: f64,        // osculating mean motion, rad/day (8.1)
    pub omega: f64,        // argument of pericentre of the star, counted from N^ (8.4)
    pub mean_anomaly: f64, // M (8.4)
    pub true_anomaly: f64, // f = u - omega (8.4)
    pub u_m: f64,          // u_M = omega + M (8.4)
    pub r_star: f64,       // star-planet distance, AU
    pub r_hat: Axes,       // direction of the star seen from the planet (8.19)
    pub rho: (f64, f64),   // rho = sin(theta*) e^{i lambda_N}, (Re, Im) (8.5)
    pub c_star: f64,       // cos(theta*) (8.5)
}

/// Mean rates of the frequencies (8.11), rad/day.
#[derive(Debug, Copy, Clone, Default, PartialEq)]
pub struct MeanRates {
    pub n_bar: f64,      // mean motion averaged over the window
    pub lambda_dot: f64, // mean rate of lambda = Omega + omega + M
    pub node_dot: f64,   // mean rate of the node Omega
}

/// Fixed-size history of the accepted steps (Section 8.4.1).
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RateHistory {
    pub initialized: bool,
    pub window: f64,      // Delta T of eq. (8.12) at the last accepted step, days
    pub last_time: f64,   // time of the last accepted step, days
    pub lambda: f64,      // lambda accumulated over the accepted steps (full turns kept), rad
    pub node: f64,        // Omega accumulated over the accepted steps, rad
    pub last_u_m: f64,    // u_M at the last accepted step
    pub last_n_hat: Axes, // N^ at the last accepted step
    // Ring buffer of (t, lambda, Omega, n_osc), oldest at index `head`
    #[serde(with = "BigArray")]
    pub samples: [[f64; 4]; HISTORY_CAPACITY],
    pub head: usize,
    pub len: usize,
}

impl Default for RateHistory {
    fn default() -> Self {
        Self {
            initialized: false,
            window: 0.,
            last_time: 0.,
            lambda: 0.,
            node: 0.,
            last_u_m: 0.,
            last_n_hat: Axes::new(),
            samples: [[0.; 4]; HISTORY_CAPACITY],
            head: 0,
            len: 0,
        }
    }
}

/// Power-law extrapolation of one end of the spectrum (8.14).
#[derive(Debug, Copy, Clone, Default, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct SpectrumEnd {
    pub beta: f64,   // exponent of imk2; 0: imk2 and rek2 held constant
    pub beta_r: f64, // exponent b of rek2 - rek2_inf; 0: rek2 held constant
}

/// Spectrum prepared once from `LoveNumber` (Section 8.4.2).
/// Unsigned (planetary tide): chi >= 0, read through eq. (8.13). Signed (stellar tide): frequencies of both
/// signs, read through `love_number_signed`; `low` is then the extrapolation beyond the negative end.
#[derive(Debug, Copy, Clone, Default, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct TideSpectrum {
    pub prepared: bool,
    pub signed: bool,
    pub n_distinct: usize, // N, length of the strictly increasing prefix of the frequencies
    pub high: SpectrumEnd, // extrapolation above chi_N
    pub has_low: bool,     // chi_1 > 0 (unsigned) or signed: extrapolation below chi_1
    pub low: SpectrumEnd,
}

/// Force and torque of the tide, and their dissipative parts (bar k -> i Im bar k).
#[derive(Debug, Copy, Clone, Default, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct Tide {
    pub force: Axes,               // F_{P->D}, eq. (8.17), M_sun AU day^-2
    pub torque: Axes,              // Gamma = d(I s)/dt = -r x F, eq. (8.18), M_sun AU^2 day^-2
    pub dissipative_force: Axes,   // F_Im, eq. (8.21)
    pub dissipative_torque: Axes,  // Gamma_Im = -r x F_Im, eq. (8.21)
}

//------------------------------------------------------------------------------
// Small complex arithmetic for the sums of eq. (8.10)

#[derive(Debug, Copy, Clone, Default, PartialEq)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    #[inline(always)]
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    #[inline(always)]
    fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    #[inline(always)]
    fn scale(self, factor: f64) -> Self {
        Self::new(self.re * factor, self.im * factor)
    }
    #[inline(always)]
    fn add(self, other: Self) -> Self {
        Self::new(self.re + other.re, self.im + other.im)
    }
    #[inline(always)]
    fn mul(self, other: Self) -> Self {
        Self::new(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )
    }
}

//------------------------------------------------------------------------------
// Geometry

fn unit(v: Axes) -> Axes {
    let mut u = v;
    u.div(v.norm());
    u
}

/// Geometry from the state vectors, eqs. (8.1) to (8.5).
/// r = r_D - r_P (AU), v = v_D - v_P (AU/day), spin s of D (rad/day), mu_orb = G (M_D + M_P).
pub fn nodal_geometry(r: Axes, v: Axes, spin: Axes, mu_orb: f64) -> NodalGeometry {
    // Perturber seen from the deformed body
    let mut r_s = r;
    r_s.negate();
    let mut v_s = v;
    v_s.negate();
    let h = r_s.cross(&v_s); // = r x v
    let h_hat = unit(h);
    let spin_rate = spin.norm();
    // A non-rotating body has no equator: any axis gives the same force (frequencies m Omega_p = 0)
    let s_hat = if spin_rate > 0. { unit(spin) } else { h_hat };

    // Node and inclination, eq. (8.2) and eq. (8.3)
    let mut n_raw = s_hat.cross(&h_hat);
    let sin_i = n_raw.norm();
    let inclination = sin_i.atan2(h_hat.dot(&s_hat));
    if sin_i < EPSILON_NODE {
        // Node undefined: inertial axis with the smallest |s_k|
        let (sx, sy, sz) = (s_hat.x.abs(), s_hat.y.abs(), s_hat.z.abs());
        let axis = if sx <= sy && sx <= sz {
            Axes::from(1., 0., 0.)
        } else if sy <= sz {
            Axes::from(0., 1., 0.)
        } else {
            Axes::from(0., 0., 1.)
        };
        n_raw = s_hat.cross(&axis);
    }
    let mut n_perp = s_hat;
    n_perp.mul(n_raw.dot(&s_hat));
    let mut n_hat = n_raw;
    n_hat.sub(&n_perp);
    let n_hat = unit(n_hat);
    let y_hat = s_hat.cross(&n_hat);
    let w_hat = unit(h_hat.cross(&n_hat));

    // Elements, eq. (8.1)
    let r_star = r_s.norm();
    let a = 1.0 / (2.0 / r_star - v_s.dot(&v_s) / mu_orb);
    let mut e_vec = v_s.cross(&h);
    e_vec.div(mu_orb);
    let mut r_hat = r_s;
    r_hat.div(r_star);
    e_vec.sub(&r_hat);
    let e = e_vec.norm();

    // Angles in the orbital plane, counted from the node, eq. (8.4)
    let n_in = unit(w_hat.cross(&h_hat));
    let u = r_hat.dot(&w_hat).atan2(r_hat.dot(&n_in));
    let omega = if e == 0. {
        0.
    } else {
        e_vec.dot(&w_hat).atan2(e_vec.dot(&n_in))
    };
    let f = u - omega;
    let eccentric_anomaly = ((1. - e * e).max(0.).sqrt() * f.sin()).atan2(e + f.cos());
    let mean_anomaly = eccentric_anomaly - e * eccentric_anomaly.sin();

    NodalGeometry {
        n_hat,
        y_hat,
        s_hat,
        spin_rate,
        inclination,
        a,
        e,
        n_osc: (mu_orb / a.powi(3)).sqrt(),
        omega,
        mean_anomaly,
        true_anomaly: f,
        u_m: omega + mean_anomaly,
        r_star,
        r_hat,
        rho: (r_hat.dot(&n_hat), r_hat.dot(&y_hat)),
        c_star: r_hat.dot(&s_hat),
    }
}

//------------------------------------------------------------------------------
// Truncation and Kaula functions

fn eccentricity_bin(e: f64) -> Option<usize> {
    Q_MAX_ECCENTRICITY_EDGES.iter().position(|&edge| e <= edge)
}

/// Truncation of the q-series, table of Section 8.5 (tol 1e-5).
/// Above e = 0.5 the largest value is kept and a warning is printed once: the guarantee
/// (8.15) no longer holds. The expansion does not exist for e >= 1.
pub fn q_max(e: f64) -> usize {
    if e.is_nan() || e >= 1. {
        panic!(
            "[PANIC {} UTC] Kaula tide: unbound orbit (e = {}), the Kaula expansion does not exist!",
            utc_now(),
            e
        );
    }
    match eccentricity_bin(e) {
        Some(bin) => Q_MAX_TABLE[bin],
        None => {
            if !WARNED_ECCENTRICITY_OUTSIDE_TABLE.swap(true, Ordering::Relaxed) {
                println!(
                    "[WARNING {} UTC] Kaula tide: eccentricity {:.3} above 0.5, q_max capped at {}: the truncation error is no longer guaranteed below 1e-5 (printed once)",
                    utc_now(),
                    e,
                    Q_LIMIT
                );
            }
            Q_LIMIT
        }
    }
}

fn utc_now() -> String {
    OffsetDateTime::now_utc()
        .format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap())
        .unwrap()
}

/// Inclination functions F_2mp(i), m, p = 0..2, closed forms of Table C.1.
pub fn inclination_functions_f2(i: f64) -> [[f64; 3]; 3] {
    let (s, c) = i.sin_cos();
    let s2 = s * s;
    [
        [-3. / 8. * s2, 3. / 4. * s2 - 1. / 2., -3. / 8. * s2],
        [
            3. / 4. * s * (1. + c),
            -3. / 4. * (2. * i).sin(),
            -3. / 4. * s * (1. - c),
        ],
        [
            3. / 4. * (1. + c).powi(2),
            3. / 2. * s2,
            3. / 4. * (1. - c).powi(2),
        ],
    ]
}

// Generalised binomial n(n-1)...(n-j+1)/j! for the three values of n needed at l = 2
#[inline(always)]
fn generalised_binomial(n: i32, j: usize) -> f64 {
    let jf = j as f64;
    let sign = if j % 2 == 0 { 1. } else { -1. };
    match n {
        0 => {
            if j == 0 {
                1.
            } else {
                0.
            }
        }
        -2 => sign * (jf + 1.),
        -4 => sign * (jf + 1.) * (jf + 2.) * (jf + 3.) / 6.,
        _ => unreachable!(),
    }
}

fn tisserand_terms(e: f64) -> usize {
    match eccentricity_bin(e) {
        Some(bin) => TISSERAND_TERMS[bin],
        None => {
            // Outside the table (no guarantee): beta^(2K) below 1e-18, plus a margin
            let beta = e / (1. + (1. - e * e).sqrt());
            let k = (18. * 10f64.ln() / (-2. * beta.ln())).ceil() as usize + 10;
            k.clamp(TISSERAND_TERMS[11], 300)
        }
    }
}

// Hansen coefficient G_2pq(e) by Tisserand's formula, eq. (C.3), for p = 0, 1, with k_terms terms
fn hansen_tisserand(p: i32, q: i32, e: f64, k_terms: usize, powers: &mut Vec<f64>) -> f64 {
    let beta = e / (1. + (1. - e * e).sqrt());
    beta.powi(q.abs()) * hansen_tisserand_reduced(p, q, e, k_terms, powers)
}

// H_pq = G_2pq / beta^|q| (Tisserand's formula without its factor beta^|q|), smooth at e = 0
fn hansen_tisserand_reduced(p: i32, q: i32, e: f64, k_terms: usize, powers: &mut Vec<f64>) -> f64 {
    let l = 2;
    let root = (1. - e * e).sqrt();
    let beta = e / (1. + root);
    // e / (2 beta), finite at e = 0
    let ratio = (1. + root) / 2.;
    let x = f64::from(l - 2 * p + q) * ratio;
    // x^r / r! for r up to the largest degree of the polynomials P and Q
    let degree_max = k_terms - 1 + q.unsigned_abs() as usize;
    powers.clear();
    powers.push(1.);
    for r in 1..=degree_max {
        let previous = powers[r - 1];
        powers.push(previous * x / r as f64);
    }
    let beta_2 = beta * beta;
    let mut beta_2k = 1.;
    let mut total = 0.;
    for k in 0..k_terms {
        let (h_p, h_q) = if q > 0 {
            (k + q as usize, k)
        } else {
            (k, k + q.unsigned_abs() as usize)
        };
        let mut p_sum = 0.;
        for r in 0..=h_p {
            let sign = if r % 2 == 0 { 1. } else { -1. };
            p_sum += generalised_binomial(2 * p - 2 * l, h_p - r) * sign * powers[r];
        }
        let mut q_sum = 0.;
        for r in 0..=h_q {
            q_sum += generalised_binomial(-2 * p, h_q - r) * powers[r];
        }
        total += p_sum * q_sum * beta_2k;
        beta_2k *= beta_2;
    }
    let sign = if q % 2 == 0 { 1. } else { -1. };
    sign * (1. + beta_2).powi(l) * total
}

/// Eccentricity functions G_2pq(e), p = 0, 1, |q| <= q_max, by Tisserand's formula (C.3).
/// Index [p][q + Q_LIMIT]; G_22q = G_20,-q. Entries with |q| > q_max are zero.
pub fn eccentricity_functions_g2_tisserand(e: f64, q_max: usize) -> [[f64; N_Q]; 2] {
    let q_max = q_max.min(Q_LIMIT) as i32;
    let k_terms = tisserand_terms(e);
    let mut powers = Vec::with_capacity(k_terms + Q_LIMIT + 1);
    let mut g = [[0.; N_Q]; 2];
    for p in 0..2 {
        for q in -q_max..=q_max {
            g[p][(q + Q_LIMIT as i32) as usize] =
                hansen_tisserand(p as i32, q, e, k_terms, &mut powers);
        }
    }
    g
}

// Chebyshev coefficients of H_pq on each interval of u: [p][q + Q_LIMIT][interval][degree]
struct EccentricityTable {
    u_max: f64,
    coefficients: Vec<[[f64; G_TABLE_DEGREE + 1]; G_TABLE_INTERVALS]>,
}

static ECCENTRICITY_TABLE: OnceLock<EccentricityTable> = OnceLock::new();

// Built once per run from Tisserand's formula at the Chebyshev nodes (first kind) of each interval, with
// G_TABLE_TISSERAND_TERMS terms for all |q| <= 42 (about 50 ms); the same at every run
fn eccentricity_table() -> &'static EccentricityTable {
    ECCENTRICITY_TABLE.get_or_init(|| {
        let e_max = G_TABLE_ECCENTRICITY;
        let beta_max = e_max / (1. + (1. - e_max * e_max).sqrt());
        let u_max = beta_max * beta_max;
        let width = u_max / G_TABLE_INTERVALS as f64;
        let n_nodes = G_TABLE_DEGREE + 1;
        let mut powers = Vec::new();
        let mut coefficients = vec![[[0.; G_TABLE_DEGREE + 1]; G_TABLE_INTERVALS]; 2 * N_Q];
        for interval in 0..G_TABLE_INTERVALS {
            let middle = (interval as f64 + 0.5) * width;
            let nodes: Vec<f64> = (0..n_nodes)
                .map(|j| (PI * (2 * j + 1) as f64 / (2 * n_nodes) as f64).cos())
                .collect();
            for p in 0..2 {
                for q in -(Q_LIMIT as i32)..=(Q_LIMIT as i32) {
                    let values: Vec<f64> = nodes
                        .iter()
                        .map(|&t| {
                            let u = middle + 0.5 * width * t;
                            let beta = u.sqrt();
                            let e = 2. * beta / (1. + u);
                            hansen_tisserand_reduced(p as i32, q, e, G_TABLE_TISSERAND_TERMS, &mut powers)
                        })
                        .collect();
                    // c_k = (2 - delta_k0) / n sum_j H(t_j) T_k(t_j)
                    let row = &mut coefficients[p * N_Q + (q + Q_LIMIT as i32) as usize][interval];
                    for (k, c) in row.iter_mut().enumerate() {
                        let sum: f64 = nodes
                            .iter()
                            .zip(&values)
                            .map(|(&t, &h)| h * (k as f64 * t.acos()).cos())
                            .sum();
                        *c = sum * if k == 0 { 1. } else { 2. } / n_nodes as f64;
                    }
                }
            }
        }
        EccentricityTable { u_max, coefficients }
    })
}

/// Eccentricity functions G_2pq(e), p = 0, 1, |q| <= q_max (C.3): Chebyshev tables for e <= 0.5
/// (absolute error ~3e-15, eq. (8.16)), Tisserand's formula above. Index [p][q + Q_LIMIT];
/// G_22q = G_20,-q. Entries with |q| > q_max are zero.
pub fn eccentricity_functions_g2(e: f64, q_max: usize) -> [[f64; N_Q]; 2] {
    if e > G_TABLE_ECCENTRICITY {
        return eccentricity_functions_g2_tisserand(e, q_max);
    }
    let table = eccentricity_table();
    let q_max = q_max.min(Q_LIMIT);
    let beta = e / (1. + (1. - e * e).sqrt());
    let u = beta * beta;
    let width = table.u_max / G_TABLE_INTERVALS as f64;
    let interval = ((u / width) as usize).min(G_TABLE_INTERVALS - 1);
    let t = (u - (interval as f64 + 0.5) * width) / (0.5 * width);
    let mut beta_powers = [1.; Q_LIMIT + 1];
    for k in 1..=q_max {
        beta_powers[k] = beta_powers[k - 1] * beta;
    }
    let mut g = [[0.; N_Q]; 2];
    for p in 0..2 {
        for q in -(q_max as i32)..=(q_max as i32) {
            let index = (q + Q_LIMIT as i32) as usize;
            // Clenshaw's recurrence for sum_k c_k T_k(t)
            let c = &table.coefficients[p * N_Q + index][interval];
            let (mut b1, mut b2) = (0., 0.);
            for &ck in c[1..].iter().rev() {
                let b0 = ck + 2. * t * b1 - b2;
                b2 = b1;
                b1 = b0;
            }
            g[p][index] = (c[0] + t * b1 - b2) * beta_powers[q.unsigned_abs() as usize];
        }
    }
    g
}

//------------------------------------------------------------------------------
// Love numbers from the tabulated spectra

fn spectrum_end(x: [f64; 3], re: [f64; 3], im: [f64; 3], high: bool) -> SpectrumEnd {
    // x[0] is the end point (chi_N or chi_1), x[1] and x[2] the next points inwards
    let hold = SpectrumEnd::default();
    if !(im[0] > 0. && im[1] > 0.) {
        return hold;
    }
    let mut beta = (im[0] / im[1]).ln() / (x[0] / x[1]).ln();
    if high {
        if beta >= 0. {
            return hold;
        }
        beta = beta.max(-1.);
    } else {
        if beta <= 0. {
            return hold;
        }
        beta = beta.min(1.);
    }
    let (b_low, b_high) = if beta < 2. * beta {
        (beta, 2. * beta)
    } else {
        (2. * beta, beta)
    };
    // Exponent b of rek2 = rek2_inf + C chi^b through the three end points, clamped to [b_low, b_high]
    let d1 = re[0] - re[1];
    let d2 = re[1] - re[2];
    let beta_r = if d1 * d2 <= 0. {
        beta
    } else {
        let g = |b: f64| (x[0].powf(b) - x[1].powf(b)) / (x[1].powf(b) - x[2].powf(b)) - d1 / d2;
        let (g_low, g_high) = (g(b_low), g(b_high));
        if g_low * g_high > 0. {
            if g_low.abs() < g_high.abs() {
                b_low
            } else {
                b_high
            }
        } else {
            let (mut lo, mut hi) = (b_low, b_high);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if g(lo) * g(mid) <= 0. {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            0.5 * (lo + hi)
        }
    };
    SpectrumEnd { beta, beta_r }
}

/// Prepare the spectrum once (Section 8.4.2): strictly increasing prefix of
/// the frequencies, padding check, extrapolation exponents at both ends.
/// Panics on an invalid spectrum.
pub fn prepare_spectrum(love: &LoveNumber) -> TideSpectrum {
    let (chi, rek2, imk2) = love.spectrum();
    let len = chi.len();
    let mut n = 1;
    while n < len && chi[n] > chi[n - 1] {
        n += 1;
    }
    let invalid = |reason: String| -> ! {
        panic!(
            "[PANIC {} UTC] Kaula tide: invalid Love-number spectrum: {}",
            utc_now(),
            reason
        )
    };
    if n < 3 {
        invalid(format!(
            "only {n} strictly increasing frequencies at the start of the table (at least 3 needed; a first row of physical parameters must be removed)"
        ));
    }
    for j in 0..n {
        if !(chi[j].is_finite() && rek2[j].is_finite() && imk2[j].is_finite()) {
            invalid(format!("non-finite value at row {j}"));
        }
        if imk2[j] < 0. {
            invalid(format!("imk2 < 0 at row {j} (imk2 >= 0 is a lag)"));
        }
    }
    if chi[0] < 0. {
        invalid(format!("negative frequency chi_1 = {}", chi[0]));
    }
    // After the prefix, only padding: zeros, or repetitions of the last row
    for j in n..len {
        if chi[j] != 0. && chi[j] != chi[n - 1] {
            invalid(format!(
                "frequencies not sorted: row {j} (chi = {:e}) after the increasing prefix of {n} rows",
                chi[j]
            ));
        }
    }
    let high = spectrum_end(
        [chi[n - 1], chi[n - 2], chi[n - 3]],
        [rek2[n - 1], rek2[n - 2], rek2[n - 3]],
        [imk2[n - 1], imk2[n - 2], imk2[n - 3]],
        true,
    );
    let has_low = chi[0] > 0.;
    let low = if has_low {
        spectrum_end(
            [chi[0], chi[1], chi[2]],
            [rek2[0], rek2[1], rek2[2]],
            [imk2[0], imk2[1], imk2[2]],
            false,
        )
    } else {
        SpectrumEnd::default()
    };
    TideSpectrum {
        prepared: true,
        signed: false,
        n_distinct: n,
        high,
        has_low,
        low,
    }
}

fn extrapolate(chi: f64, x: [f64; 2], re: [f64; 2], im0: f64, end: &SpectrumEnd) -> (f64, f64) {
    // eq. (8.14), from the end point x[0] and its neighbour x[1]
    let imk2 = if end.beta != 0. {
        im0 * (chi / x[0]).powf(end.beta)
    } else {
        im0
    };
    let rek2 = if end.beta_r == 0. {
        re[0]
    } else {
        let (x0, x1) = (x[0].powf(end.beta_r), x[1].powf(end.beta_r));
        re[0] + (re[0] - re[1]) * (chi.powf(end.beta_r) - x0) / (x0 - x1)
    };
    (rek2, imk2)
}

/// True if chi lies more than a decade outside the tabulated spectrum (above 10 chi_N, or below chi_1 / 10 when
/// chi_1 > 0; chi = 0 excluded): (8.14) is then used far from the table, with no guaranteed accuracy.
pub fn far_outside_spectrum(love: &LoveNumber, spec: &TideSpectrum, chi: f64) -> bool {
    let (x, _, _) = love.spectrum();
    chi > 10. * x[spec.n_distinct - 1] || (spec.has_low && chi > 0. && chi < 0.1 * x[0])
}

// Warning printed once per run when a forcing frequency is more than a decade outside the spectrum
fn warn_far_outside_spectrum(love: &LoveNumber, spec: &TideSpectrum, chi: f64) {
    if far_outside_spectrum(love, spec, chi) && !WARNED_FREQUENCY_OUTSIDE_SPECTRUM.swap(true, Ordering::Relaxed) {
        let (x, _, _) = love.spectrum();
        println!(
            "[WARNING {} UTC] Kaula tide: forcing frequency {:e} rad/s more than a decade outside the tabulated spectrum [{:e}, {:e}] rad/s: Love number extrapolated, accuracy not guaranteed (printed once)",
            utc_now(),
            chi,
            x[0],
            x[spec.n_distinct - 1]
        );
    }
}

/// (rek2, imk2) at chi >= 0 in rad/s: linear interpolation inside the table, eq. (8.14) outside.
/// `hint` is the interval found by the previous call (any value is valid; it only speeds up the search).
pub fn evaluate_spectrum(
    love: &LoveNumber,
    spec: &TideSpectrum,
    chi: f64,
    hint: &mut usize,
) -> (f64, f64) {
    let (x_all, re_all, im_all) = love.spectrum();
    let n = spec.n_distinct;
    let (x, re, im) = (&x_all[..n], &re_all[..n], &im_all[..n]);
    if chi >= x[n - 1] {
        warn_far_outside_spectrum(love, spec, chi);
        return extrapolate(
            chi,
            [x[n - 1], x[n - 2]],
            [re[n - 1], re[n - 2]],
            im[n - 1],
            &spec.high,
        );
    }
    if chi <= x[0] {
        warn_far_outside_spectrum(love, spec, chi);
        if !spec.has_low || chi == x[0] {
            return (re[0], im[0]);
        }
        if chi == 0. {
            if spec.low.beta == 0. {
                return (re[0], im[0]);
            }
            let rek2 = if spec.low.beta_r == 0. {
                re[0]
            } else {
                let (x0, x1) = (x[0].powf(spec.low.beta_r), x[1].powf(spec.low.beta_r));
                re[0] - (re[0] - re[1]) * x0 / (x0 - x1)
            };
            return (rek2, 0.);
        }
        return extrapolate(chi, [x[0], x[1]], [re[0], re[1]], im[0], &spec.low);
    }
    // Interval x[j] < chi <= x[j + 1]
    let j = if *hint + 1 < n && x[*hint] < chi && chi <= x[*hint + 1] {
        *hint
    } else {
        x.partition_point(|&value| value < chi) - 1
    };
    *hint = j;
    let t = (chi - x[j]) / (x[j + 1] - x[j]);
    (
        (1. - t) * re[j] + t * re[j + 1],
        (1. - t) * im[j] + t * im[j + 1],
    )
}

/// bar k_2(sigma) = rek2(|sigma|) - i sgn(sigma) imk2(|sigma|), eq. (8.13), with chi = |sigma| / DAY (Section 10.3).
/// sigma in rad/day. Returns (Re bar k_2, Im bar k_2).
pub fn love_number_at(
    love: &LoveNumber,
    spec: &TideSpectrum,
    sigma: f64,
    hint: &mut usize,
) -> (f64, f64) {
    let (rek2, imk2) = evaluate_spectrum(love, spec, sigma.abs() / DAY, hint);
    if sigma == 0. {
        (rek2, 0.)
    } else {
        (rek2, -sigma.signum() * imk2)
    }
}

/// Prepare a signed spectrum once (stellar tide, Section 9.2): strictly increasing prefix
/// x_1 < ... < x_N with x_1 < 0 < x_N, padding check, Im k2 of the sign of the frequency at every row,
/// extrapolation exponents of (8.14) in |x| at both ends (held constant when an end has fewer than
/// 3 points of its sign). Panics on an invalid spectrum.
pub fn prepare_signed_spectrum(love: &LoveNumber) -> TideSpectrum {
    let (x, rek2, imk2) = love.spectrum();
    let len = x.len();
    let mut n = 1;
    while n < len && x[n] > x[n - 1] {
        n += 1;
    }
    let invalid = |reason: String| -> ! {
        panic!(
            "[PANIC {} UTC] Kaula tide: invalid signed Love-number spectrum: {}",
            utc_now(),
            reason
        )
    };
    if n < 3 {
        invalid(format!(
            "only {n} strictly increasing frequencies at the start of the table (at least 3 needed)"
        ));
    }
    for j in 0..n {
        if !(x[j].is_finite() && rek2[j].is_finite() && imk2[j].is_finite()) {
            invalid(format!("non-finite value at row {j}"));
        }
        if imk2[j] * x[j] < 0. {
            invalid(format!(
                "row {j}: imk2 = {:e} has the sign opposite to the frequency {:e}",
                imk2[j], x[j]
            ));
        }
    }
    if !(x[0] < 0. && x[n - 1] > 0.) {
        invalid(format!(
            "frequencies of both signs needed (first {:e}, last {:e} rad/s)",
            x[0],
            x[n - 1]
        ));
    }
    for j in n..len {
        if x[j] != 0. && x[j] != x[n - 1] {
            invalid(format!(
                "frequencies not sorted: row {j} (x = {:e}) after the increasing prefix of {n} rows",
                x[j]
            ));
        }
    }
    let high = if x[n - 3] > 0. {
        spectrum_end(
            [x[n - 1], x[n - 2], x[n - 3]],
            [rek2[n - 1], rek2[n - 2], rek2[n - 3]],
            [imk2[n - 1], imk2[n - 2], imk2[n - 3]],
            true,
        )
    } else {
        SpectrumEnd::default()
    };
    // Negative end in |x|: |x| grows outwards, as above the positive end
    let low = if x[2] < 0. {
        spectrum_end(
            [-x[0], -x[1], -x[2]],
            [rek2[0], rek2[1], rek2[2]],
            [-imk2[0], -imk2[1], -imk2[2]],
            true,
        )
    } else {
        SpectrumEnd::default()
    };
    TideSpectrum {
        prepared: true,
        signed: true,
        n_distinct: n,
        high,
        has_low: true,
        low,
    }
}

// Warning printed once per run when a forcing frequency falls outside a signed spectrum
fn warn_outside_signed_spectrum(x_first: f64, x_last: f64, x: f64) {
    if !WARNED_FREQUENCY_OUTSIDE_SIGNED_SPECTRUM.swap(true, Ordering::Relaxed) {
        println!(
            "[WARNING {} UTC] Kaula tide: forcing frequency {:e} rad/s outside the signed spectrum [{:e}, {:e}] rad/s: Love number extrapolated (8.14), accuracy not guaranteed (printed once)",
            utc_now(),
            x,
            x_first,
            x_last
        );
    }
}

/// (rek2, imk2) of a signed spectrum at x in rad/s, any sign: linear interpolation inside the table,
/// also across 0 (imk2(0) != 0 kept); (8.14) in |x| outside, with the sign of imk2 restored.
/// `hint` as in `evaluate_spectrum`.
pub fn evaluate_signed_spectrum(
    love: &LoveNumber,
    spec: &TideSpectrum,
    x_eval: f64,
    hint: &mut usize,
) -> (f64, f64) {
    let (x_all, re_all, im_all) = love.spectrum();
    let n = spec.n_distinct;
    let (x, re, im) = (&x_all[..n], &re_all[..n], &im_all[..n]);
    if x_eval >= x[n - 1] {
        if x_eval > x[n - 1] {
            warn_outside_signed_spectrum(x[0], x[n - 1], x_eval);
        }
        return extrapolate(
            x_eval,
            [x[n - 1], x[n - 2]],
            [re[n - 1], re[n - 2]],
            im[n - 1],
            &spec.high,
        );
    }
    if x_eval <= x[0] {
        if x_eval < x[0] {
            warn_outside_signed_spectrum(x[0], x[n - 1], x_eval);
        }
        let (rek2, imk2) =
            extrapolate(-x_eval, [-x[0], -x[1]], [re[0], re[1]], -im[0], &spec.low);
        return (rek2, -imk2);
    }
    // Interval x[j] < x_eval <= x[j + 1]
    let j = if *hint + 1 < n && x[*hint] < x_eval && x_eval <= x[*hint + 1] {
        *hint
    } else {
        x.partition_point(|&value| value < x_eval) - 1
    };
    *hint = j;
    let t = (x_eval - x[j]) / (x[j + 1] - x[j]);
    (
        (1. - t) * re[j] + t * re[j + 1],
        (1. - t) * im[j] + t * im[j + 1],
    )
}

/// bar k_2(sigma) from a signed spectrum (stellar tide), eq. (9.1): the
/// frequency of the file is sigma, and bar k_2 = conj k_file = rek2 - i imk2 (imk2 has the sign of the
/// frequency: Im bar k_2 < 0 for a lag at sigma > 0, as eq. (8.13)). With a spectrum spin rate Omega_ref
/// (rescaling of commit ad904f5), the file is read at sigma Omega_ref / Omega and imk2 is multiplied by
/// (Omega / Omega_ref)^2. sigma and spin_rate (Omega) in rad/day. Returns (Re bar k_2, Im bar k_2).
pub fn love_number_signed(
    love: &LoveNumber,
    spec: &TideSpectrum,
    sigma: f64,
    spin_rate: f64,
    hint: &mut usize,
) -> (f64, f64) {
    let (x, scale) = match love.spectrum_spin_rate() {
        Some(reference) => {
            if !(spin_rate > 0.) {
                panic!(
                    "[PANIC {} UTC] Kaula tide: spin rate {:e} rad/day, the rescaling of the spectrum by the spin needs a rotating body",
                    utc_now(),
                    spin_rate
                );
            }
            let spin = spin_rate / DAY;
            (sigma / DAY * reference / spin, (spin / reference).powi(2))
        }
        None => (sigma / DAY, 1.),
    };
    let (rek2, imk2) = evaluate_signed_spectrum(love, spec, x, hint);
    (rek2, -imk2 * scale)
}

//------------------------------------------------------------------------------
// Mean rates

fn wrap_angle(x: f64) -> f64 {
    // (x + pi) mod 2 pi - pi, with the sign convention of Python's %
    (x + PI).rem_euclid(TWO_PI) - PI
}

/// Osculating mean motion about the host (rad/day), 0 for an unbound body (8.12).
/// r, v: body relative to the host.
pub fn mean_motion_about_host(r: Axes, v: Axes, mu: f64) -> f64 {
    let inverse_a = 2. / r.norm() - v.dot(&v) / mu;
    if inverse_a > 0. {
        (mu * inverse_a.powi(3)).sqrt()
    } else {
        0.
    }
}

/// Averaging window of a planet, eq. (8.12): longest synodic period with the other bodies orbiting
/// the host, clamped to [P, 10 P]. Zero mean motions (unbound bodies) are ignored.
pub fn averaging_window(n_planet: f64, n_others: &[f64]) -> f64 {
    let period = TWO_PI / n_planet;
    let mut synodic = 0.;
    for &n_other in n_others.iter().filter(|&&n_other| n_other > 0.) {
        synodic = f64::max(synodic, TWO_PI / (n_planet - n_other).abs());
    }
    synodic.max(period).min(10. * period)
}

impl RateHistory {
    fn sample(&self, k: usize) -> [f64; 4] {
        // k = 0 is the oldest sample
        self.samples[(self.head + k) % HISTORY_CAPACITY]
    }
}

/// Append an accepted step to the history (Section 8.4.1): accumulate lambda
/// and Omega, store a sample every window / 32. Called once per accepted step; a second call at
/// the same time is ignored.
pub fn update_rate_history(history: &mut RateHistory, t: f64, geo: &NodalGeometry, window: f64) {
    if history.initialized && t == history.last_time {
        return;
    }
    if !history.initialized {
        history.lambda = geo.u_m;
        history.node = 0.;
        history.initialized = true;
    } else {
        let previous = history.last_n_hat;
        let d_node = geo
            .s_hat
            .dot(&previous.cross(&geo.n_hat))
            .atan2(previous.dot(&geo.n_hat));
        history.node += d_node;
        history.lambda += wrap_angle(d_node + geo.u_m - history.last_u_m);
    }
    history.last_n_hat = geo.n_hat;
    history.last_u_m = geo.u_m;
    history.last_time = t;
    history.window = window;

    let store = history.len == 0
        || t - history.sample(history.len - 1)[0] >= window / HISTORY_SAMPLES_PER_WINDOW;
    if store {
        let entry = [t, history.lambda, history.node, geo.n_osc];
        if history.len < HISTORY_CAPACITY {
            let index = (history.head + history.len) % HISTORY_CAPACITY;
            history.samples[index] = entry;
            history.len += 1;
        } else {
            history.samples[history.head] = entry;
            history.head = (history.head + 1) % HISTORY_CAPACITY;
        }
    }
}

/// Mean rates from the history (read only), Section 8.4.1.
/// lambda and Omega at t - window are interpolated linearly between the stored samples; n_bar is the
/// mean of the stored n from there on. Start-up (history shorter than the window):
/// lambda_dot = n_bar = n_osc, node_dot = 0. Buffer full but shorter than the window: oldest sample.
pub fn mean_rates(history: &RateHistory, n_osc: f64) -> MeanRates {
    let start_up = MeanRates {
        n_bar: n_osc,
        lambda_dot: n_osc,
        node_dot: 0.,
    };
    if !history.initialized || history.len == 0 || history.window <= 0. {
        return start_up;
    }
    let (t1, lambda1, node1) = (history.last_time, history.lambda, history.node);
    let t_window = t1 - history.window;
    // Last stored sample at or before t - window
    let mut older = None;
    for k in 0..history.len {
        if history.sample(k)[0] <= t_window {
            older = Some(k);
        } else {
            break;
        }
    }
    let (t0, lambda0, node0, first) = match older {
        Some(k) => {
            let b0 = history.sample(k);
            let b1 = if k + 1 < history.len {
                history.sample(k + 1)
            } else {
                [t1, lambda1, node1, 0.]
            };
            if b1[0] > b0[0] {
                let w = (t_window - b0[0]) / (b1[0] - b0[0]);
                (
                    t_window,
                    b0[1] + (b1[1] - b0[1]) * w,
                    b0[2] + (b1[2] - b0[2]) * w,
                    k,
                )
            } else {
                (b0[0], b0[1], b0[2], k)
            }
        }
        None => {
            if history.len < HISTORY_CAPACITY {
                return start_up;
            }
            // Window longer than the buffer: shorter window from the oldest sample
            let b0 = history.sample(0);
            (b0[0], b0[1], b0[2], 0)
        }
    };
    if t1 <= t0 {
        return start_up;
    }
    let mut n_sum = 0.;
    for k in first..history.len {
        n_sum += history.sample(k)[3];
    }
    MeanRates {
        n_bar: n_sum / (history.len - first) as f64,
        lambda_dot: (lambda1 - lambda0) / (t1 - t0),
        node_dot: (node1 - node0) / (t1 - t0),
    }
}

/// Frequency of the term (m, p, q), eq. (8.11), rad/day.
#[inline(always)]
pub fn frequency(m: usize, p: usize, q: i32, rates: &MeanRates, spin_rate: f64) -> f64 {
    let k = m as i32 - 2 + 2 * p as i32;
    f64::from(2 - 2 * p as i32) * rates.lambda_dot + f64::from(q) * rates.n_bar
        + f64::from(k) * rates.node_dot
        - m as f64 * spin_rate
}

//------------------------------------------------------------------------------
// Gradient, force and torque

// e^{i q M} for |q| <= q_max, index q + Q_LIMIT, by successive products of e^{i M} (Section 8.9)
fn phases_q(mean_anomaly: f64, q_max: usize) -> [Complex; N_Q] {
    let (s, c) = mean_anomaly.sin_cos();
    let step = Complex::new(c, s);
    let mut phase = [Complex::new(1., 0.); N_Q];
    for q in 1..=q_max.min(Q_LIMIT) {
        phase[Q_LIMIT + q] = phase[Q_LIMIT + q - 1].mul(step);
        phase[Q_LIMIT - q] = phase[Q_LIMIT + q].conj();
    }
    phase
}

// Love numbers of the m = 0 terms, [p][q + Q_LIMIT], for the p flagged in `needed`. Hermitian spectrum
// (planetary tide, Section 8.9): one evaluation per pair, exact since bar k(-sigma) = conj bar k(sigma) by
// construction (8.13) and sigma_{0,1,-q} = -sigma_{0,1,q}, sigma_{0,2,q} = -sigma_{0,0,-q} hold exactly in
// floating point. Signed spectrum (stellar tide), without this symmetry: every term is evaluated.
fn love_numbers_m0<L: FnMut(usize, usize, i32, f64) -> (f64, f64)>(
    love: &mut L,
    rates: &MeanRates,
    spin_rate: f64,
    q_max: usize,
    needed: [bool; 3],
    hermitian: bool,
) -> [[(f64, f64); N_Q]; 3] {
    let q_max = q_max.min(Q_LIMIT) as i32;
    let index = |q: i32| (q + Q_LIMIT as i32) as usize;
    let mut k = [[(0., 0.); N_Q]; 3];
    if !hermitian {
        for p in 0..3 {
            if needed[p] {
                for q in -q_max..=q_max {
                    k[p][index(q)] = love(0, p, q, frequency(0, p, q, rates, spin_rate));
                }
            }
        }
        return k;
    }
    if needed[0] || needed[2] {
        for q in -q_max..=q_max {
            k[0][index(q)] = love(0, 0, q, frequency(0, 0, q, rates, spin_rate));
        }
        if needed[2] {
            for q in -q_max..=q_max {
                let (re, im) = k[0][index(-q)];
                k[2][index(q)] = (re, -im);
            }
        }
    }
    if needed[1] {
        for q in 0..=q_max {
            let (re, im) = love(0, 1, q, frequency(0, 1, q, rates, spin_rate));
            k[1][index(q)] = (re, im);
            k[1][index(-q)] = (re, -im);
        }
    }
    k
}

/// Gradient of T at the perturber in the nodal frame, eqs. (8.6) to (8.10), split into the parts due to
/// Re bar k and to i Im bar k (8.21). `love(sigma)` returns (Re bar k_2, Im bar k_2) at sigma
/// in rad/day; `hermitian`: bar k(-sigma) = conj bar k(sigma) (planetary spectra, see love_numbers_m0). Returns [(V_X + i V_Y, V_Z) of the Re part, same for the Im part], in AU day^-2
/// per unit mass of the perturber.
pub fn tidal_gradient_nodal<L: FnMut(usize, usize, i32, f64) -> (f64, f64)>(
    geo: &NodalGeometry,
    rates: &MeanRates,
    mut love: L,
    hermitian: bool,
    gm_perturber: f64,
    radius: f64,
    q_max: usize,
) -> [((f64, f64), f64); 2] {
    let q_max = q_max.min(Q_LIMIT);
    let f2 = inclination_functions_f2(geo.inclination);
    let g2 = eccentricity_functions_g2(geo.e, q_max);
    let c = geo.c_star;
    let p3 = [
        0.5 * (5. * c.powi(3) - 3. * c),
        1.5 * (5. * c * c - 1.),
        15. * c,
        15.,
    ];
    let rho = Complex::new(geo.rho.0, geo.rho.1);
    let rho_bar = rho.conj();
    let one = Complex::new(1., 0.);
    let rho_pow = [one, rho, rho.mul(rho), rho.mul(rho).mul(rho)];
    let rho_bar_pow = [one, rho_bar, rho_bar.mul(rho_bar)];
    // e^{i f_2m}: f_20 = f_22 = 0, f_21 = -pi/2
    let parity = [one, Complex::new(0., -1.), one];
    // e^{i (2-2p) u_M}
    let mut phase_p = [one; 3];
    for (p, phase) in phase_p.iter_mut().enumerate() {
        let (s, co) = ((2. - 2. * p as f64) * geo.u_m).sin_cos();
        *phase = Complex::new(co, s);
    }
    // e^{i q M} by successive products of e^{i M} (rounding error ~ |q| 1e-16)
    let phase_q = phases_q(geo.mean_anomaly, q_max);
    let scale = gm_perturber * radius.powi(5) / geo.a.powi(3);

    // m = 0: one evaluation per pair for a hermitian spectrum (love_numbers_m0)
    let love_m0 = love_numbers_m0(&mut love, rates, geo.spin_rate, q_max, [f2[0][0] != 0., f2[0][1] != 0., f2[0][2] != 0.], hermitian);

    let mut v_xy = [Complex::default(); 2];
    let mut v_z = [0.; 2];
    for m in 0..3 {
        for p in 0..3 {
            let f_mp = f2[m][p];
            if f_mp == 0. {
                continue;
            }
            for q in -(q_max as i32)..=(q_max as i32) {
                // G_22q = G_20,-q
                let g = if p < 2 {
                    g2[p][(q + Q_LIMIT as i32) as usize]
                } else {
                    g2[0][(Q_LIMIT as i32 - q) as usize]
                };
                let k_mpq = scale * C_M[m] * f_mp * g;
                if k_mpq == 0. {
                    continue;
                }
                let (re_k, im_k) = if m == 0 {
                    love_m0[p][(q + Q_LIMIT as i32) as usize]
                } else {
                    love(m, p, q, frequency(m, p, q, rates, geo.spin_rate))
                };
                let phase = phase_p[p]
                    .mul(phase_q[(q + Q_LIMIT as i32) as usize])
                    .mul(parity[m]);
                // E_mpq = bar k e^{i Psi}, eq. (8.8), split into Re bar k and i Im bar k
                let e_parts = [phase.scale(re_k), Complex::new(0., im_k).mul(phase)];
                for part in 0..2 {
                    let e_mpq = e_parts[part];
                    v_z[part] +=
                        (3 - m) as f64 * k_mpq * p3[m] * rho_bar_pow[m].mul(e_mpq).re;
                    if m == 0 {
                        v_xy[part] = v_xy[part].add(rho.scale(k_mpq * p3[1] * e_mpq.re));
                    } else {
                        let first = rho_pow[m + 1].mul(e_mpq.conj()).scale(p3[m + 1]);
                        let second = rho_bar_pow[m - 1]
                            .mul(e_mpq)
                            .scale(((3 - m) * (4 - m)) as f64 * p3[m - 1]);
                        v_xy[part] = v_xy[part].add(
                            first
                                .add(second.scale(-1.))
                                .scale(0.5 * k_mpq),
                        );
                    }
                }
            }
        }
    }
    let r4 = geo.r_star.powi(4);
    [
        ((v_xy[0].re / r4, v_xy[0].im / r4), v_z[0] / r4),
        ((v_xy[1].re / r4, v_xy[1].im / r4), v_z[1] / r4),
    ]
}

/// Gradient of T at the perturber for a coplanar orbit, eq. (8.19) (the single form eq. (8.10) at i = 0,
/// only (m, p) = (0, 1) and (2, 0)), split into the parts due to Re bar k and to i Im bar k (8.21).
/// `love(m, p, q, sigma)` and `hermitian` as in `tidal_gradient_nodal`. Returns [(radial, azimuthal) of the Re part,
/// same for the Im part], along r^* and s^ x r^*, in AU day^-2 per unit mass of the perturber.
pub fn tidal_gradient_coplanar<L: FnMut(usize, usize, i32, f64) -> (f64, f64)>(
    geo: &NodalGeometry,
    rates: &MeanRates,
    mut love: L,
    hermitian: bool,
    gm_perturber: f64,
    radius: f64,
    q_max: usize,
) -> [(f64, f64); 2] {
    let q_max = q_max.min(Q_LIMIT);
    let g2 = eccentricity_functions_g2(geo.e, q_max);
    // Love numbers of (0, 1, q) (one per pair for a hermitian spectrum) and phases by successive products
    let love_m0 = love_numbers_m0(&mut love, rates, geo.spin_rate, q_max, [false, true, false], hermitian);
    let phase_q = phases_q(geo.mean_anomaly, q_max);
    let (s2f, c2f) = (2. * geo.true_anomaly).sin_cos();
    let phase_2m_2f = phases_q(geo.mean_anomaly, 2)[Q_LIMIT + 2].mul(Complex::new(c2f, -s2f));
    // w_r and w_phi of eq. (8.19), [Re part, Im part]
    let mut w_r = [0.; 2];
    let mut w_phi = [0.; 2];
    for q in -(q_max as i32)..=(q_max as i32) {
        let index = (q + Q_LIMIT as i32) as usize;
        let (g_20q, g_21q) = (g2[0][index], g2[1][index]);
        // m = 0, p = 1: 3/4 G_21q Re(bar k_01q e^{i q M}), sigma = q n
        let phase = phase_q[index];
        if g_21q != 0. {
            let (re_k, im_k) = love_m0[1][index];
            let (s, c) = (phase.im, phase.re);
            w_r[0] += 0.75 * g_21q * re_k * c;
            w_r[1] -= 0.75 * g_21q * im_k * s;
        }
        // m = 2, p = 0: Z_q = bar k_20q e^{i((2+q) M - 2 f)}, sigma = 2 lambda_dot + q n - 2 Omega_p
        if g_20q != 0. {
            let (re_k, im_k) = love(2, 0, q, frequency(2, 0, q, rates, geo.spin_rate));
            // e^{i((2+q) M - 2 f)} = e^{i q M} e^{i(2M - 2f)}
            let z = phase.mul(phase_2m_2f);
            let (s, c) = (z.im, z.re);
            // 9/4 G_20q Re Z_q and -3/2 G_20q Im Z_q
            w_r[0] += 2.25 * g_20q * re_k * c;
            w_phi[0] -= 1.5 * g_20q * re_k * s;
            w_r[1] -= 2.25 * g_20q * im_k * s;
            w_phi[1] -= 1.5 * g_20q * im_k * c;
        }
    }
    let scale = gm_perturber * radius.powi(5) / (geo.a.powi(3) * geo.r_star.powi(4));
    [
        (w_r[0] * scale, w_phi[0] * scale),
        (w_r[1] * scale, w_phi[1] * scale),
    ]
}

fn coplanar_to_inertial(geo: &NodalGeometry, gradient: (f64, f64), perturber_mass: f64) -> Axes {
    // eq. (8.19): F = M_P [radial r^* + azimuthal (s^ x r^*)]
    let (radial, azimuthal) = gradient;
    let phi_hat = geo.s_hat.cross(&geo.r_hat);
    Axes::from(
        perturber_mass * (radial * geo.r_hat.x + azimuthal * phi_hat.x),
        perturber_mass * (radial * geo.r_hat.y + azimuthal * phi_hat.y),
        perturber_mass * (radial * geo.r_hat.z + azimuthal * phi_hat.z),
    )
}

fn nodal_to_inertial(geo: &NodalGeometry, gradient: ((f64, f64), f64), perturber_mass: f64) -> Axes {
    // eq. (8.17): F = M_P [Re(V_X + i V_Y) N^ + Im(V_X + i V_Y) Y^ + V_Z s^]
    let ((v_x, v_y), v_z) = gradient;
    Axes::from(
        perturber_mass * (v_x * geo.n_hat.x + v_y * geo.y_hat.x + v_z * geo.s_hat.x),
        perturber_mass * (v_x * geo.n_hat.y + v_y * geo.y_hat.y + v_z * geo.s_hat.y),
        perturber_mass * (v_x * geo.n_hat.z + v_y * geo.y_hat.z + v_z * geo.s_hat.z),
    )
}

fn torque_from_force(r: Axes, force: Axes) -> Axes {
    // eq. (8.18): Gamma = -r x F
    let mut torque = r.cross(&force);
    torque.negate();
    torque
}

fn assemble_tide(r: Axes, force_re: Axes, dissipative_force: Axes) -> Tide {
    let mut force = force_re;
    force.add(&dissipative_force);
    Tide {
        force,
        torque: torque_from_force(r, force),
        dissipative_force,
        dissipative_torque: torque_from_force(r, dissipative_force),
    }
}

/// Force of P on D (8.17), torque on the spin of D (8.18) and their dissipative parts
/// (8.21), inertial axes, with the single form eq. (8.10) (any obliquity). r = r_D - r_P.
/// `love(m, p, q, sigma)` returns (Re bar k_2, Im bar k_2); `hermitian` as in `tidal_gradient_nodal`.
#[allow(clippy::too_many_arguments)]
pub fn tide_single<L: FnMut(usize, usize, i32, f64) -> (f64, f64)>(
    r: Axes,
    geo: &NodalGeometry,
    rates: &MeanRates,
    love: L,
    hermitian: bool,
    g: f64,
    perturber_mass: f64,
    radius: f64,
    q_max: usize,
) -> Tide {
    let [real_part, imaginary_part] =
        tidal_gradient_nodal(geo, rates, love, hermitian, g * perturber_mass, radius, q_max);
    assemble_tide(
        r,
        nodal_to_inertial(geo, real_part, perturber_mass),
        nodal_to_inertial(geo, imaginary_part, perturber_mass),
    )
}

/// Same with the coplanar form eq. (8.19) (valid for i <= COPLANAR_INCLINATION).
#[allow(clippy::too_many_arguments)]
pub fn tide_coplanar<L: FnMut(usize, usize, i32, f64) -> (f64, f64)>(
    r: Axes,
    geo: &NodalGeometry,
    rates: &MeanRates,
    love: L,
    hermitian: bool,
    g: f64,
    perturber_mass: f64,
    radius: f64,
    q_max: usize,
) -> Tide {
    let [real_part, imaginary_part] =
        tidal_gradient_coplanar(geo, rates, love, hermitian, g * perturber_mass, radius, q_max);
    assemble_tide(
        r,
        coplanar_to_inertial(geo, real_part, perturber_mass),
        coplanar_to_inertial(geo, imaginary_part, perturber_mass),
    )
}

/// Tide: coplanar form eq. (8.19) if i <= COPLANAR_INCLINATION (prograde orbit in the
/// equator), single form eq. (8.10) otherwise (Section 8.7).
#[allow(clippy::too_many_arguments)]
pub fn tide<L: FnMut(usize, usize, i32, f64) -> (f64, f64)>(
    r: Axes,
    geo: &NodalGeometry,
    rates: &MeanRates,
    love: L,
    hermitian: bool,
    g: f64,
    perturber_mass: f64,
    radius: f64,
    q_max: usize,
) -> Tide {
    if geo.inclination <= COPLANAR_INCLINATION {
        tide_coplanar(r, geo, rates, love, hermitian, g, perturber_mass, radius, q_max)
    } else {
        tide_single(r, geo, rates, love, hermitian, g, perturber_mass, radius, q_max)
    }
}

/// Tide with the tabulated spectrum, unsigned (love_number_at) or signed (love_number_signed, spin rate
/// of D), with one search hint per (m, p), since the frequencies of consecutive q are close.
#[allow(clippy::too_many_arguments)]
pub fn tide_tabulated(
    r: Axes,
    geo: &NodalGeometry,
    rates: &MeanRates,
    love: &LoveNumber,
    spec: &TideSpectrum,
    g: f64,
    perturber_mass: f64,
    radius: f64,
    q_max: usize,
) -> Tide {
    let mut hints = [[0usize; 3]; 3];
    if spec.signed {
        let spin_rate = geo.spin_rate;
        tide(
            r,
            geo,
            rates,
            |m, p, _q, sigma| love_number_signed(love, spec, sigma, spin_rate, &mut hints[m][p]),
            false,
            g,
            perturber_mass,
            radius,
            q_max,
        )
    } else {
        tide(
            r,
            geo,
            rates,
            |m, p, _q, sigma| love_number_at(love, spec, sigma, &mut hints[m][p]),
            true,
            g,
            perturber_mass,
            radius,
            q_max,
        )
    }
}
