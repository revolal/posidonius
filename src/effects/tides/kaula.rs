// Implemented by Alexandre Revol alexandre.revol@unige.ch
use super::super::super::constants::{G, TWO_PI, DAY};
use super::super::super::tools;
use super::super::super::{Particle};
use super::super::super::{Axes};

// Strucure declare the intput parameters needed for the Kaula calculation
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct KaulaParameters {
    // The tidal Love numbers function of the excitation frequency
    // the data are stored in 2d array 32 by 32.
    pub love_number_excitation_frequency: [[f64;32];32], // The excitation frequency 
    pub real_part_love_number: [[f64;32];32], // The real part of the complex love number
    pub imaginary_part_love_number: [[f64;32];32], // The imaginary part of the complex love number
    pub num_datapoints: f64, // The number of point stored in the 32 by 32 array
}

// Compute the tidal force applied by the host particle to the particle
// Called from common.rs
pub fn calculate_tidal_force(tidal_host_particle: &Particle, particle: &mut Particle, number_data_point:f64, excitation_frequency_tab:[[f64;32];32], real_love_number_tab:[[f64;32];32], imag_love_number_tab:[[f64;32];32]) -> Axes {
    // --- The spherical coordinate --- //
    // The following elements correspond to the coordinate in the spherical coordinate
    // The coplanar distance is the radial distance projected in the x-y plane
    // The theta angle is the angle btw the radial distance vector with respect to the z axis
    // The phi angle is the angle btw the coplanar distance with respect to the x axis
    // ---
    let radial_distance:f64 = particle.tides.parameters.internal.distance;
    let coplanar_distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) ).sqrt();
    let cos_theta:f64 = particle.tides.coordinates.position.z / radial_distance;
    let sin_theta:f64 = coplanar_distance / radial_distance;
    let cos_phi:f64 = particle.tides.coordinates.position.x / coplanar_distance;
    let sin_phi:f64 = particle.tides.coordinates.position.y / coplanar_distance;

    // --- The Keplerian orbital elements --- //
    let gm = G*(tidal_host_particle.mass+particle.mass);
    let keplerian_elements = tools::calculate_keplerian_orbital_elements( gm, particle.heliocentric_position, particle.heliocentric_velocity);
    // let (a, _q, eccentricity, inclination,longitude_perihelion, longitude_of_ascending_node, mean_anomaly, _orbital_period) = keplerian_elements;
    // println!("\n|\t sma {} \n|\t ecc {} \n|\t incl {} \n|\t long_peri {} {} \n|\t long_asc_node {} {}  \n|\t mean_ano {}", a, eccentricity, inclination, longitude_perihelion, TWO_PI+longitude_perihelion, longitude_of_ascending_node, TWO_PI+longitude_of_ascending_node, mean_anomaly);
    // panic!("At the disco");
    // --- The spherical component of the tidal force --- //
    // The radial component is the force applicated through the radial axis
    // The normal component act on the co longitude axis
    // The orthogonal component act on the co latitude axis
    // ---
    let (radial_component_of_the_tidal_force, radial_component_of_the_tidal_force_secular) = calculate_radial_component_of_the_tidal_force(tidal_host_particle, particle, keplerian_elements, number_data_point, excitation_frequency_tab, real_love_number_tab, imag_love_number_tab );
    let (normal_component_of_the_tidal_force, normal_component_of_the_tidal_force_secular) = calculate_normal_component_of_the_tidal_force(tidal_host_particle, particle, keplerian_elements, number_data_point, excitation_frequency_tab, real_love_number_tab, imag_love_number_tab );
    let (orthogonal_component_of_the_tidal_force, orthogonal_component_of_the_tidal_force_secular) = calculate_orthogonal_component_of_the_tidal_force(tidal_host_particle, particle, keplerian_elements, number_data_point, excitation_frequency_tab, real_love_number_tab, imag_love_number_tab );

    // --- The cartesian tidal force --- // computed by projection of the spherical coordinates
    let tidal_force_x = radial_component_of_the_tidal_force
                                *sin_theta *cos_phi 
                            + normal_component_of_the_tidal_force
                                *cos_theta *cos_phi 
                            - orthogonal_component_of_the_tidal_force 
                                *sin_phi ;
    let tidal_force_y = radial_component_of_the_tidal_force
                                *sin_theta *sin_phi
                            + normal_component_of_the_tidal_force
                                *cos_theta *sin_phi
                            + orthogonal_component_of_the_tidal_force 
                                *cos_phi ;
    let tidal_force_z = radial_component_of_the_tidal_force
                                *cos_theta
                            - normal_component_of_the_tidal_force
                                *sin_theta ;
    // --- The tidal torque --- // general 
    // particle.tides.parameters.internal.kaula_tidal_force.x = radial_component_of_the_tidal_force
    //                             *sin_theta *cos_phi 
    //                         + normal_component_of_the_tidal_force
    //                             *cos_theta *cos_phi 
    //                         - orthogonal_component_of_the_tidal_force
    //                             *sin_phi ;
    //                         particle.tides.parameters.internal.kaula_tidal_force.y = radial_component_of_the_tidal_force
    //                             *sin_theta *sin_phi
    //                         + normal_component_of_the_tidal_force
    //                             *cos_theta *sin_phi
    //                         + orthogonal_component_of_the_tidal_force
    //                             *cos_phi ;
    //                         particle.tides.parameters.internal.kaula_tidal_force.z = radial_component_of_the_tidal_force
    //                             *cos_theta
    //                         - normal_component_of_the_tidal_force
    //                             *sin_theta ;

    // // --- The tidal force --- // secular
    // let tidal_force_x = radial_component_of_the_tidal_force_secular
    //                             *sin_theta *cos_phi 
    //                         + normal_component_of_the_tidal_force_secular
    //                             *cos_theta *cos_phi 
    //                         - orthogonal_component_of_the_tidal_force_secular 
    //                             *sin_phi ;
    // let tidal_force_y = radial_component_of_the_tidal_force_secular
    //                             *sin_theta *sin_phi
    //                         + normal_component_of_the_tidal_force_secular
    //                             *cos_theta *sin_phi
    //                         + orthogonal_component_of_the_tidal_force_secular
    //                             *cos_phi ;
    // let tidal_force_z = radial_component_of_the_tidal_force_secular
    //                             *cos_theta
    //                         - normal_component_of_the_tidal_force_secular
    //                             *sin_theta ;
    // --- The tidal torque --- // return the secular part of the tidal torque (simplified from the rapid varying phases)
    particle.tides.parameters.internal.kaula_tidal_force.x = radial_component_of_the_tidal_force_secular
                                *sin_theta *cos_phi 
                            + normal_component_of_the_tidal_force_secular
                                *cos_theta *cos_phi 
                            - orthogonal_component_of_the_tidal_force_secular
                                *sin_phi ;
    particle.tides.parameters.internal.kaula_tidal_force.y = radial_component_of_the_tidal_force_secular
                                *sin_theta *sin_phi
                            + normal_component_of_the_tidal_force_secular
                                *cos_theta *sin_phi
                            + orthogonal_component_of_the_tidal_force_secular 
                                *cos_phi ;
    particle.tides.parameters.internal.kaula_tidal_force.z = radial_component_of_the_tidal_force_secular
                                *cos_theta
                            - normal_component_of_the_tidal_force_secular
                                *sin_theta ;


    return Axes{x: tidal_force_x, y: tidal_force_y, z: tidal_force_z}
}

// ------------------------------------- //
// --- Calculate tidal force modules --- //
// Function which compute the components of the tidal torque

// --- The radial (e_{r}) component of tidal force
fn calculate_radial_component_of_the_tidal_force(
    tidal_host_particle: & Particle, particle: &Particle, keplerian_elements:(f64, f64, f64, f64, f64, f64, f64, f64), 
    number_data_point:f64, excitation_frequency_tab:[[f64;32];32], real_love_number_tab:[[f64;32];32], imag_love_number_tab:[[f64;32];32] ) -> (f64,f64) {
        // --- The keplerian elements
        // --- 
        let ( semi_major_axis, _perihelion_distance, eccentricity, _inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period 
        ) = keplerian_elements;
        let obliquity = tools::calculate_inclination_orbital_equatorial_plane(particle.heliocentric_position, particle.heliocentric_velocity, particle.spin);
        let argument_perihelion = longitude_perihelion -longitude_of_ascending_node ;
        let orbital_frequency = TWO_PI / (orbital_period *DAY); // The orbital mean motion in [rad.s^-1]
        let spin = particle.norm_spin_vector_2.sqrt() / DAY; // The planetary spin in [rad.s^-1]

        // ---
        // --- local quantities needed --- //
        let cste:f64 = -(G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6)*particle.heliocentric_distance) ;

        // - Kaula inclination and eccentricity functions
        // let eccentricity_function_g_2pq = calculate_eccentricty_function_g_2pq(eccentricity);

        // --- summation --- //
        let radial_force:f64;
        let radial_force_secular:f64;
        // let case_2d:bool = true;

        if obliquity <= 1.0e-8 { // --- if the inclination btw the equatorial plane and the orbital plane is negligible
            if eccentricity == 0.0 { // --- If circular coplanar orbit
                let frequ_2010 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0., 1., 0., spin, orbital_frequency);
                let frequ_2200 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., 0., 0., spin, orbital_frequency);
                let (rek2_2010, _imk2_2010) = calculate_kaula_numbers(frequ_2010, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                let (_rek2_2200, imk2_2200) = calculate_kaula_numbers(frequ_2200, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                radial_force = cste *( (3_f64/4_f64) *rek2_2010 +(9_f64/4_f64) *imk2_2200 ) ;
                radial_force_secular = radial_force;

            } else { //  --- If eccentric orbit 
                let eccentricity_function_g_20q = calculate_eccentricty_function_g_20q(eccentricity);
                let eccentricity_function_g_21q = calculate_eccentricty_function_g_21q(eccentricity);
                let mut sum_over_q:f64 = 0.;
                let mut sum_over_q_secular:f64 = 0.;

                let order = match () { // Select the order of the summation q over the eccentricity function G_lpq
                    _ if (eccentricity > 0.3 ) => 7, // if ecc > 0.20 take q btw -7 and 7
                    _ if (eccentricity > 0.25) => 6, // if ecc > 0.20 take q btw -6 and 6
                    _ if (eccentricity > 0.2 ) => 5, // if ecc > 0.20 take q btw -5 and 5
                    _ if (eccentricity > 0.15) => 4, // if ecc > 0.15 take q btw -4 and 4
                    _ if (eccentricity > 0.1 ) => 3, // if ecc > 0.10 take q btw -3 and 3 
                    _ if (eccentricity <= 0.1 ) => 2, // if ecc > 0.10 take q btw -2 and 2
                    _ => 0,
                };
                // let q_max = [7-order, 2*order+1];
                let q_min =7-order;
                let q_max =2*order+1;

                // for q in -7..8{
                // for (q, g_201q) in eccentricity_function_g_20q.iter().zip(eccentricity_function_g_21q.iter()).skip(q_max[0]).take(q_max[1]).enumerate(){
                for (q, g_201q) in eccentricity_function_g_20q.iter().zip(eccentricity_function_g_21q.iter()).skip(q_min).take(q_max).enumerate(){
                    // let frequ_201q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0., 1., q as f64, spin, orbital_frequency);
                    // let frequ_220q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                    let frequ_201q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0., 1., q as f64 -order as f64, spin, orbital_frequency);
                    let frequ_220q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., 0., q as f64 -order as f64, spin, orbital_frequency);
                    let (rek2_201q, imk2_201q) = calculate_kaula_numbers(frequ_201q, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                    let (rek2_220q, imk2_220q) = calculate_kaula_numbers(frequ_220q, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );

                    let mut sum_over_j_1:f64=0.;
                    let mut sum_over_j_3:f64=0.;

                    let mut sum_over_j_1_secular:f64=0.;
                    let mut sum_over_j_3_secular:f64=0.;
                    // for j in -7..8{
                    for (j, g_201j) in eccentricity_function_g_20q.iter().zip(eccentricity_function_g_21q.iter()).skip(q_min).take(q_max).enumerate(){
                        // let alpha_qj = alpha_pqkj(0., q as f64, 0., j as f64, mean_anomaly, argument_perihelion );
                        let alpha_qj = alpha_pqkj(0., q as f64 -order as f64, 0., j as f64 -order as f64, mean_anomaly, argument_perihelion );
                        // let g_21j:f64 = eccentricity_function_g_2pq[1][(j+7) as usize];
                        // let g_20j:f64 = eccentricity_function_g_2pq[0][(j+7) as usize];
                        let (g_20j, g_21j) = g_201j;

                        sum_over_j_1 = sum_over_j_1 +g_21j*(
                            alpha_qj.cos()*rek2_201q 
                            -alpha_qj.sin()*imk2_201q
                        );
                        sum_over_j_3 = sum_over_j_3 +g_20j*(
                            alpha_qj.cos()*rek2_220q 
                            -alpha_qj.sin()*imk2_220q
                        );
                        if q==j {
                            sum_over_j_1_secular = sum_over_j_1_secular +g_21j *rek2_201q ;
                            sum_over_j_3_secular = sum_over_j_3_secular +g_20j *rek2_220q ;
                        }
                    } // end loop over j
                    // let g_21q:f64 = eccentricity_function_g_2pq[1][(q+7) as usize];
                    // let g_20q:f64 = eccentricity_function_g_2pq[0][(q+7) as usize];
                    let (g_20q, g_21q) = g_201q;
                    sum_over_q = sum_over_q + (
                        (3_f64/4_f64) *g_21q *sum_over_j_1 
                        +(9_f64/4_f64)*g_20q *sum_over_j_3
                    );
                    sum_over_q_secular = sum_over_q_secular + (
                        (3_f64/4_f64) *g_21q *sum_over_j_1_secular
                        +(9_f64/4_f64)*g_20q *sum_over_j_3_secular
                    );
                }
                radial_force = cste *sum_over_q;
                radial_force_secular = cste *sum_over_q_secular;
            }
        } else { // --- if there s an non-negligible inclination (3D case)
            let eccentricity_function_g_2pq = calculate_eccentricty_function_g_2pq(eccentricity);
            let inclination_function_f_2mp = calculate_inclination_function_f_2mp(obliquity);

            let heliocentric_r = particle.heliocentric_distance;
            let heliocentric_varphi: f64=0.;
            let spin_angle:f64=0.;
            let host_particle_star_2:f64 = tidal_host_particle.mass.powi(2);
            let particle_radius_5:f64 = particle.radius.powi(5);
            let semi_major_axis_6:f64 = semi_major_axis.powi(6);
            let _cot_theta:f64=0.;

            let inclination_function_f_20p = calculate_inclination_function_f_20p(obliquity);
            let inclination_function_f_21p = calculate_inclination_function_f_21p(obliquity);
            let inclination_function_f_22p = calculate_inclination_function_f_22p(obliquity);
            let eccentricity_function_g_2pq = calculate_eccentricty_function_g_2pq(eccentricity);
            // let eccentricity_function_g_20q = calculate_eccentricty_function_g_20q(eccentricity);
            // let eccentricity_function_g_21q = calculate_eccentricty_function_g_21q(eccentricity);
            // let eccentricity_function_g_22q = calculate_eccentricty_function_g_22q(eccentricity);

            let order = match () { // Select the order of the summation q over the eccentricity function G_lpq
                _ if (eccentricity > 0.3 ) => 7, // if ecc > 0.20 take q btw -7 and 7
                _ if (eccentricity > 0.25) => 6, // if ecc > 0.20 take q btw -6 and 6
                _ if (eccentricity > 0.2 ) => 5, // if ecc > 0.20 take q btw -5 and 5
                _ if (eccentricity > 0.15) => 4, // if ecc > 0.15 take q btw -4 and 4
                _ if (eccentricity > 0.1 ) => 3, // if ecc > 0.10 take q btw -3 and 3 
                _ if (eccentricity <= 0.1 ) => 2, // if ecc > 0.10 take q btw -2 and 2
                _ => 0,
            };
            let q_min = 7-order;
            let q_max = 2*order+1;

            let mut sum_over_p_m0:f64=0.;
            let mut sum_over_p_m1:f64=0.;
            let mut sum_over_p_m2:f64=0.;
            let mut sum_over_p_m0_s:f64=0.;
            let mut sum_over_p_m1_s:f64=0.;
            let mut sum_over_p_m2_s:f64=0.;
            for (p, f_2mp) in inclination_function_f_20p.iter().zip(inclination_function_f_21p.iter().zip(inclination_function_f_22p.iter())).enumerate() {
                // println!("|\t p = {:?}", p);

                let (f_20p, (f_21p,f_22p) )= f_2mp;

                let tmp_p:f64 = p as f64;
                let f_20p_2 = f_20p.powi(2);
                let f_21p_2 = f_21p.powi(2);
                let f_22p_2 = f_22p.powi(2);

                let mut sum_over_q_m0:f64=0.;
                let mut sum_over_q_m1:f64=0.;
                let mut sum_over_q_m2:f64=0.;
                let mut sum_over_q_m0_s:f64=0.;
                let mut sum_over_q_m1_s:f64=0.;
                let mut sum_over_q_m2_s:f64=0.;

                let sum_g_2pq = eccentricity_function_g_2pq[p];

                for (q, g_2pq) in sum_g_2pq.iter().skip(q_min).take(q_max).enumerate() {
                    
                    let tmp_q:f64 = q as f64 -order as f64;
                    let g_2pq_2 = g_2pq.powi(2);

                    let frequ_20pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0., tmp_p, tmp_q, spin, orbital_frequency);
                    let frequ_21pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(1., tmp_p, tmp_q, spin, orbital_frequency);
                    let frequ_22pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., tmp_p, tmp_q, spin, orbital_frequency);
                    let (rek2_20pq, imk2_20pq) = calculate_kaula_numbers(frequ_20pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                    let (rek2_21pq, imk2_21pq) = calculate_kaula_numbers(frequ_21pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                    let (rek2_22pq, imk2_22pq) = calculate_kaula_numbers(frequ_22pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );

                    let mut sum_over_k_m0:f64=0.;
                    let mut sum_over_k_m1:f64=0.;
                    let mut sum_over_k_m2:f64=0.;
                    let mut sum_over_k_m0_s:f64=0.;
                    let mut sum_over_k_m1_s:f64=0.;
                    let mut sum_over_k_m2_s:f64=0.;

                    for (k, f_2mk) in inclination_function_f_20p.iter().zip(inclination_function_f_21p.iter().zip(inclination_function_f_22p.iter())).enumerate() {

                        let (f_20k, (f_21k,f_22k) )= f_2mk;

                        let tmp_k:f64 = k as f64;
                        let mut sum_over_j_m0:f64=0.;
                        let mut sum_over_j_m1:f64=0.;
                        let mut sum_over_j_m2:f64=0.;
                        let mut sum_over_j_m0_s:f64=0.;
                        let mut sum_over_j_m1_s:f64=0.;
                        let mut sum_over_j_m2_s:f64=0.;
                        let sum_g_2kj = eccentricity_function_g_2pq[k];

                        for (j, g_2kj) in sum_g_2kj.iter().skip(q_min).take(q_max).enumerate() {
                            
                            let tmp_j:f64 = j as f64 -order as f64;
                            let phase_alpha = alpha_pqkj(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, argument_perihelion );
                            
                            let cos_alpha:f64 = phase_alpha.cos();
                            let sin_alpha:f64 = phase_alpha.sin();

                            sum_over_j_m0 += g_2kj * (cos_alpha *rek2_20pq -sin_alpha*imk2_20pq) ;
                            sum_over_j_m1 += g_2kj * (cos_alpha *rek2_21pq -sin_alpha*imk2_21pq) ;
                            sum_over_j_m2 += g_2kj * (cos_alpha *rek2_22pq -sin_alpha*imk2_22pq) ;
                            
                        }
                        sum_over_k_m0 += f_20k *sum_over_j_m0 ;
                        sum_over_k_m1 += f_21k *sum_over_j_m1 ;
                        sum_over_k_m2 += f_22k *sum_over_j_m2 ;
                        sum_over_k_m0_s += f_20k *sum_over_j_m0_s ;
                        sum_over_k_m1_s += f_21k *sum_over_j_m1_s ;
                        sum_over_k_m2_s += f_22k *sum_over_j_m2_s ;
                    }
                    sum_over_q_m0 += g_2pq *sum_over_k_m0;
                    sum_over_q_m1 += g_2pq *sum_over_k_m1;
                    sum_over_q_m2 += g_2pq *sum_over_k_m2;
                    sum_over_q_m0_s += g_2pq_2 *rek2_20pq ;
                    sum_over_q_m1_s += g_2pq_2 *rek2_21pq ;
                    sum_over_q_m2_s += g_2pq_2 *rek2_22pq ;
                }
                sum_over_p_m0 += f_20p *sum_over_q_m0;
                sum_over_p_m1 += f_21p *sum_over_q_m1;
                sum_over_p_m2 += f_22p *sum_over_q_m2;
                sum_over_p_m0_s += f_20p_2 *sum_over_q_m0_s;
                sum_over_p_m1_s += f_21p_2 *sum_over_q_m1_s;
                sum_over_p_m2_s += f_22p_2 *sum_over_q_m2_s;
                
            }

            let cste = (G *host_particle_star_2 *particle_radius_5 *3.0_f64) / (semi_major_axis_6 *heliocentric_r);
            radial_force =  (sum_over_p_m0 +(1_f64/3_f64)*sum_over_p_m1 + (1_f64/12_f64)*sum_over_p_m2) *cste  ;
            radial_force_secular = (sum_over_p_m0_s +(1_f64/3_f64)*sum_over_p_m1_s + (1_f64/12_f64)*sum_over_p_m2_s) *cste  ;

        } // End if 3D
    return (radial_force, radial_force_secular) ;
}

// --- The Normal (e_{\theta}) component of tidal force
fn calculate_normal_component_of_the_tidal_force(
    tidal_host_particle: &Particle, particle: &Particle, keplerian_elements:(f64, f64, f64, f64, f64, f64, f64, f64), 
    number_data_point:f64, excitation_frequency_tab:[[f64;32];32], real_love_number_tab:[[f64;32];32], imag_love_number_tab:[[f64;32];32] ) -> (f64,f64) {
        // --- The keplerian elements
        // --- 
        let ( semi_major_axis, _perihelion_distance, eccentricity, _inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period 
        ) = keplerian_elements;
        let obliquity = tools::calculate_inclination_orbital_equatorial_plane(particle.heliocentric_position, particle.heliocentric_velocity, particle.spin);
        let argument_perihelion = longitude_perihelion -longitude_of_ascending_node ;
        let orbital_frequency = TWO_PI / (orbital_period *DAY); // The orbital mean motion in [rad.s^-1]
        let spin = particle.norm_spin_vector_2.sqrt() / DAY; // The planetary spin in [rad.s^-1]

        // println!("|\t --- Check orbital elements --- ");

        // println!("| \t sma {:?}\n| \t ecc {:?}\n| \t incl {:?}\n| \t long_asc {:?}\n| \t long_peri {:?}\n| \t mean_ano {:?}\n| \t mean_long {:?} ",
        // semi_major_axis, eccentricity, tools::modulus(_inclination,TWO_PI), tools::modulus(longitude_of_ascending_node,TWO_PI),  tools::modulus(longitude_perihelion,TWO_PI),  tools::modulus(mean_anomaly,TWO_PI),  tools::modulus(longitude_perihelion+mean_anomaly,TWO_PI) ); 

        // println!("| \t obliquity {:?} rad {:?} deg  ",obliquity, obliquity*360./TWO_PI );
        // println!("| \t inclination {:?} rad {:?} deg  ",_inclination, _inclination*360./TWO_PI );
        // println!("| \t spin {:?} {:?}", spin, TWO_PI/spin/DAY);
        // println!("| \t orbital frequency {:?} | \t period {:?}", orbital_frequency, orbital_period);

        // panic!("At the disco");

        let normal_force:f64;
        let normal_force_secular:f64;

        if obliquity <= 1.0e-8 {
            normal_force = 0.;
            normal_force_secular = 0.;
        } else { // 3D case
            let heliocentric_r = particle.heliocentric_distance;
            let heliocentric_varphi: f64=0.;
            let spin_angle:f64=0.;
            let host_particle_star_2:f64 = tidal_host_particle.mass.powi(2);
            let particle_radius_5:f64 = particle.radius.powi(5);
            let semi_major_axis_6:f64 = semi_major_axis.powi(6);
            let _cot_theta:f64=0.;

            let inclination_function_f_20p = calculate_inclination_function_f_20p(obliquity);
            let inclination_function_f_21p = calculate_inclination_function_f_21p(obliquity);
            let inclination_function_f_22p = calculate_inclination_function_f_22p(obliquity);
            let eccentricity_function_g_2pq = calculate_eccentricty_function_g_2pq(eccentricity);
            // let eccentricity_function_g_20q = calculate_eccentricty_function_g_20q(eccentricity);
            // let eccentricity_function_g_21q = calculate_eccentricty_function_g_21q(eccentricity);
            // let eccentricity_function_g_22q = calculate_eccentricty_function_g_22q(eccentricity);

            let order = match () { // Select the order of the summation q over the eccentricity function G_lpq
                _ if (eccentricity > 0.3 ) => 7, // if ecc > 0.20 take q btw -7 and 7
                _ if (eccentricity > 0.25) => 6, // if ecc > 0.20 take q btw -6 and 6
                _ if (eccentricity > 0.2 ) => 5, // if ecc > 0.20 take q btw -5 and 5
                _ if (eccentricity > 0.15) => 4, // if ecc > 0.15 take q btw -4 and 4
                _ if (eccentricity > 0.1 ) => 3, // if ecc > 0.10 take q btw -3 and 3 
                _ if (eccentricity <= 0.1 ) => 2, // if ecc > 0.10 take q btw -2 and 2
                _ => 0,
            };
            let q_min = 7-order;
            let q_max = 2*order+1;

            // term m = 0
            let mut sum_over_p:f64=0.;
            let mut sum_over_p_s:f64=0.;
            for (p, f_20p) in inclination_function_f_20p.iter().enumerate() {
                // println!("|\t p = {:?}", p);
                let tmp_p:f64 = p as f64;
                let mut sum_over_q:f64=0.;
                let mut sum_over_q_s:f64=0.;
                let sum_g_2pq = eccentricity_function_g_2pq[p];

                for (q, g_2pq) in sum_g_2pq.iter().skip(q_min).take(q_max).enumerate() {
                    
                    let tmp_q:f64 = q as f64 -order as f64;
                    // println!("|\t \t q = {:?} {:?}", q, tmp_q);
                    let frequ_20pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0., tmp_p, tmp_q, spin, orbital_frequency);
                    let mut sum_over_k:f64=0.;
                    let mut sum_over_k_s:f64=0.;

                    for (k, f_21k) in inclination_function_f_21p.iter().enumerate() {
                        // println!("|\t \t \t k = {:?}", k);
                        let tmp_k:f64 = k as f64;
                        let mut sum_over_j:f64=0.;
                        let mut sum_over_j_s:f64=0.;
                        let sum_g_2kj = eccentricity_function_g_2pq[k];

                        for (j, g_2kj) in sum_g_2kj.iter().skip(q_min).take(q_max).enumerate() {
                            
                            let tmp_j:f64 = j as f64 -order as f64;
                            // println!("|\t \t \t \t j = {:?} {:?}", j, tmp_j);
                            let phase_beta:f64 =compute_phase_beta(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );
                            let phase_alpha_1:f64 =compute_phase_alpha_1(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );

                            let (rek2_20pq, imk2_20pq) = calculate_kaula_numbers(frequ_20pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                            
                            let cos_alpha_1:f64 = phase_alpha_1.cos();
                            let sin_alpha_1:f64 = phase_alpha_1.sin();
                            let cos_beta:f64 = phase_beta.cos();
                            let sin_beta:f64 = phase_beta.sin();

                            sum_over_j = sum_over_j +g_2kj *( 0.5_f64*(cos_alpha_1 *rek2_20pq -sin_alpha_1 *imk2_20pq) +1.5_f64*(cos_beta *rek2_20pq -sin_beta *imk2_20pq) );
                            sum_over_j_s = sum_over_j_s +g_2kj *( 0.5_f64*rek2_20pq +1.5_f64*rek2_20pq  );
                            
                        }
                        sum_over_k = sum_over_k + f_21k *sum_over_j ;
                        sum_over_k_s = sum_over_k_s + f_21k *sum_over_j_s ;
                    }
                    sum_over_q = sum_over_q + g_2pq *sum_over_k;
                    sum_over_q_s = sum_over_q_s + g_2pq *sum_over_k_s;
                }
                sum_over_p = sum_over_p + f_20p *sum_over_q;
                sum_over_p_s = sum_over_p_s + f_20p *sum_over_q_s;
            }
            let _term_m0:f64 = sum_over_p ;
            let _term_m0_s:f64 = sum_over_p_s ;

            // term m = 1
            sum_over_p=0.;
            sum_over_p_s=0.;
            for (p, f_21p) in inclination_function_f_21p.iter().enumerate() {
                
                let tmp_p:f64 = p as f64;
                let mut sum_over_q:f64=0.;
                let mut sum_over_q_s:f64=0.;
                let sum_g_2pq = eccentricity_function_g_2pq[p];

                for (q, g_2pq) in sum_g_2pq.iter().skip(q_min).take(q_max).enumerate() {

                    let tmp_q:f64 = q as f64 -order as f64;
                    let frequ_21pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(1., tmp_p, tmp_q, spin, orbital_frequency);

                    let mut sum_over_k_m2:f64=0.;
                    let mut sum_over_k_m0:f64=0.;
                    let mut sum_over_k_m2_s:f64=0.;
                    let mut sum_over_k_m0_s:f64=0.;

                    for (k, f_22_20_k) in inclination_function_f_22p.iter().zip(inclination_function_f_20p.iter()).enumerate() {

                        let (f_22k, f_20k) = f_22_20_k; // is sum over p for each m

                        let tmp_k:f64 = k as f64;
                        let mut sum_over_j_p2:f64=0.;
                        let mut sum_over_j_p0:f64=0.;
                        let mut sum_over_j_p2_s:f64=0.;
                        let mut sum_over_j_p0_s:f64=0.;
                        let sum_g_2kj = eccentricity_function_g_2pq[k];

                        for (j, g_2kj) in sum_g_2kj.iter().skip(q_min).take(q_max).enumerate() {

                            let tmp_j:f64 = j as f64 -order as f64;

                            // let phase_beta:f64 =compute_phase_beta(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_ascending_node );
                            let phase_alpha_1:f64 =compute_phase_alpha_1(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );
                            let phase_alpha_2:f64 =compute_phase_alpha_2(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );

                            let (rek2_21pq, imk2_21pq) = calculate_kaula_numbers(frequ_21pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                            
                            let cos_alpha_1:f64 = phase_alpha_1.cos();
                            let sin_alpha_1:f64 = phase_alpha_1.sin();
                            let cos_alpha_2:f64 = phase_alpha_2.cos();
                            let sin_alpha_2:f64 = phase_alpha_2.sin();

                            sum_over_j_p2 = sum_over_j_p2 +g_2kj *(cos_alpha_1 *rek2_21pq -sin_alpha_1 *imk2_21pq );
                            sum_over_j_p0 = sum_over_j_p0 +g_2kj *(cos_alpha_2 *rek2_21pq -sin_alpha_2 *imk2_21pq );
                            sum_over_j_p2_s = sum_over_j_p2_s +g_2kj *rek2_21pq ;
                            sum_over_j_p0_s = sum_over_j_p0_s +g_2kj *rek2_21pq ;
                        }
                        sum_over_k_m2 = sum_over_k_m2 + f_22k *sum_over_j_p2 ;
                        sum_over_k_m0 = sum_over_k_m0 + f_20k *sum_over_j_p0 ;
                        sum_over_k_m2_s = sum_over_k_m2_s + f_22k *sum_over_j_p2_s ;
                        sum_over_k_m0_s = sum_over_k_m0_s + f_20k *sum_over_j_p0_s ;
                    }
                    sum_over_q = sum_over_q + g_2pq * (0.5_f64 *sum_over_k_m2 -3_f64*sum_over_k_m0);
                    sum_over_q_s = sum_over_q_s + g_2pq * (0.5_f64 *sum_over_k_m2_s -3_f64*sum_over_k_m0_s);
                }
                sum_over_p = sum_over_p + f_21p *sum_over_q;
                sum_over_p_s = sum_over_p_s + f_21p *sum_over_q_s;
            }
            let _term_m1:f64 = (1_f64/3_f64) *sum_over_p ;
            let _term_m1_s:f64 = (1_f64/3_f64) *sum_over_p_s ;

            // term m = 2
            sum_over_p=0.;
            sum_over_p_s=0.;
            for (p, f_22p) in inclination_function_f_22p.iter().enumerate() {
                
                let tmp_p:f64 = p as f64;
                let mut sum_over_q:f64=0.;
                let mut sum_over_q_s:f64=0.;
                let sum_g_2pq = eccentricity_function_g_2pq[p];

                for (q, g_2pq) in sum_g_2pq.iter().skip(q_min).take(q_max).enumerate() {

                    let tmp_q:f64 = q as f64 -order as f64;
                    let frequ_22pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., tmp_p, tmp_q, spin, orbital_frequency);

                    let mut sum_over_k:f64=0.;
                    let mut sum_over_k_s:f64=0.;

                    for (k, f_21k) in inclination_function_f_21p.iter().enumerate() {

                        let tmp_k:f64 = k as f64;
                        let mut sum_over_j:f64=0.;
                        let mut sum_over_j_s:f64=0.;

                        let sum_g_2kj = eccentricity_function_g_2pq[k];

                        for (j, g_2kj) in sum_g_2kj.iter().skip(q_min).take(q_max).enumerate() {

                            let tmp_j:f64 = j as f64 -order as f64;

                            // let phase_beta:f64 =compute_phase_beta(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_ascending_node );
                            let phase_alpha_2:f64 =compute_phase_alpha_2(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );

                            let (rek2_22pq, imk2_22pq) = calculate_kaula_numbers(frequ_22pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );

                            let cos_alpha_2:f64 = phase_alpha_2.cos();
                            let sin_alpha_2:f64 = phase_alpha_2.sin();

                            sum_over_j = sum_over_j +g_2kj *(cos_alpha_2 *rek2_22pq -sin_alpha_2 *imk2_22pq );
                            sum_over_j_s = sum_over_j_s +g_2kj *rek2_22pq ;
                        }
                        sum_over_k = sum_over_k + f_21k *sum_over_j ;
                        sum_over_k_s = sum_over_k_s + f_21k *sum_over_j_s ;
                    }
                    sum_over_q = sum_over_q + g_2pq * sum_over_k ;
                    sum_over_q_s = sum_over_q_s + g_2pq * sum_over_k_s ;
                }
                sum_over_p = sum_over_p  + f_22p *sum_over_q;
                sum_over_p_s = sum_over_p_s  + f_22p *sum_over_q_s;
            }
            let _term_m2:f64 = - (1_f64/6_f64) *sum_over_p ;
            let _term_m2_s:f64 = - (1_f64/6_f64) *sum_over_p_s ;

            let cste:f64 = (G *host_particle_star_2 *particle_radius_5) / (semi_major_axis_6 *heliocentric_r) ;
            
            normal_force = cste * (_term_m0 + _term_m1 +_term_m2);
            normal_force_secular = cste * (_term_m0_s + _term_m1_s +_term_m2_s);

        }// End if case 3D
    return (normal_force, normal_force_secular) ;
}

// --- The Ortho-radial (the e_{\varphi}) component of tidal force
fn calculate_orthogonal_component_of_the_tidal_force(
    tidal_host_particle: &Particle, particle: &Particle, keplerian_elements:(f64, f64, f64, f64, f64, f64, f64, f64),
    number_data_point:f64, excitation_frequency_tab:[[f64;32];32], real_love_number_tab:[[f64;32];32],  imag_love_number_tab:[[f64;32];32] ) ->(f64,f64) {
        // --- The keplerian elements
        // --- 
        let ( semi_major_axis, _perihelion_distance, eccentricity, _inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period 
        ) = keplerian_elements;
        let argument_perihelion = longitude_perihelion -longitude_of_ascending_node ;
        let orbital_frequency = TWO_PI / (orbital_period *DAY); // The orbital mean motion in [rad.s^-1]
        let spin = particle.norm_spin_vector_2.sqrt() / DAY; // The planetary spin in [rad.s^-1]
        let obliquity = tools::calculate_inclination_orbital_equatorial_plane(particle.heliocentric_position, particle.heliocentric_velocity, particle.spin);

        // --- local quantities needed --- //

        // --- The cartesian angles --- //
        // --- 
        let distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) +particle.tides.coordinates.position.z.powi(2) ).sqrt();
        let coplanar_distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) ).sqrt();
        let sin_theta:f64 = coplanar_distance / distance;

        

        // Kaula eccentricity function
        // let eccentricity_function_g_2pq = calculate_eccentricty_function_g_2pq(eccentricity);

        // --- summation --- //
        let orthogonal_force:f64;
        let orthogonal_force_secular:f64;

        if obliquity <= 1.0e-8 {
            let cste:f64 = -(G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6)*particle.heliocentric_distance*sin_theta) ;
            if eccentricity == 0. {
                let frequ_2200 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., 0., 0., spin, orbital_frequency);
                let (_rek2_2200, imk2_2200) = calculate_kaula_numbers(frequ_2200, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                orthogonal_force = cste *(3_f64/2_f64) *imk2_2200;
                orthogonal_force_secular = orthogonal_force;
            } else {
                let eccentricity_function_g_20q = calculate_eccentricty_function_g_20q(eccentricity);
                let mut sum_over_q:f64=0.;
                let mut sum_over_q_secular:f64=0.;

                let order = match () { // Select the order of the summation q over the eccentricity function G_lpq
                    _ if (eccentricity > 0.3 ) => 7, // if ecc > 0.20 take q btw -7 and 7
                    _ if (eccentricity > 0.25) => 6, // if ecc > 0.20 take q btw -6 and 6
                    _ if (eccentricity > 0.2 ) => 5, // if ecc > 0.20 take q btw -5 and 5
                    _ if (eccentricity > 0.15) => 4, // if ecc > 0.15 take q btw -4 and 4
                    _ if (eccentricity > 0.1 ) => 3, // if ecc > 0.10 take q btw -3 and 3 
                    _ if (eccentricity <= 0.1 ) => 2, // if ecc > 0.10 take q btw -2 and 2
                    _ => 0,
                };
                // let q_max = [7-order, 2*order+1];
                let q_min = 7-order;
                let q_max = 2*order+1;
                
                // for q in -7..8{
                for (q, g_20q) in eccentricity_function_g_20q.iter().skip(q_min).take(q_max).enumerate(){
                    // let frequ_220q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                    let frequ_220q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., 0., q as f64 -order as f64, spin, orbital_frequency);
                    let (rek2_220q, imk2_220q) = calculate_kaula_numbers(frequ_220q, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );

                    let mut sum_over_j_2:f64=0.;
                    let mut sum_over_j_2_secular:f64=0.;
                    // for j in -7..8{
                    for (j, g_20j) in eccentricity_function_g_20q.iter().skip(q_min).take(q_max).enumerate(){
                        // let alpha_qj = alpha_pqkj(0., q as f64, 0., j as f64, mean_anomaly, argument_perihelion );
                        let alpha_qj = alpha_pqkj(0., q as f64 -order as f64, 0., j as f64 -order as f64, mean_anomaly, argument_perihelion );
                        // let g_20j = eccentricity_function_g_2pq[0][(j+7) as usize];

                        sum_over_j_2 = sum_over_j_2 +g_20j*(
                                alpha_qj.sin()*rek2_220q 
                                +alpha_qj.cos()*imk2_220q
                                );
                        if q==j { // --- the secular part of the force 
                            sum_over_j_2_secular = sum_over_j_2_secular +g_20j *imk2_220q;
                        }
                    } 
                    // let g_20q = eccentricity_function_g_2pq[0][(q+7) as usize];

                    sum_over_q = sum_over_q +g_20q *(3_f64/2_f64) *sum_over_j_2 ;
                    sum_over_q_secular = sum_over_q_secular +g_20q *(3_f64/2_f64) *sum_over_j_2_secular ;
                }
                orthogonal_force = cste *sum_over_q;
                orthogonal_force_secular = cste *sum_over_q_secular;
            }
        } else {
            let heliocentric_varphi: f64=0.;
            let spin_angle:f64=0.;
            let _cot_theta:f64=0.;

            let inclination_function_f_21p = calculate_inclination_function_f_21p(obliquity);
            let inclination_function_f_22p = calculate_inclination_function_f_22p(obliquity);
            let inclination_function_f_30p = calculate_inclination_function_f_30p(obliquity);
            let inclination_function_f_31p = calculate_inclination_function_f_31p(obliquity);
            let inclination_function_f_32p = calculate_inclination_function_f_32p(obliquity);
            let inclination_function_f_33p = calculate_inclination_function_f_33p(obliquity);
            let eccentricity_function_g_2pq = calculate_eccentricty_function_g_2pq(eccentricity);
            let eccentricity_function_g_3pq = calculate_eccentricty_function_g_3pq(eccentricity);

            let order = match () { // Select the order of the summation q over the eccentricity function G_lpq
                _ if (eccentricity > 0.3 ) => 7, // if ecc > 0.20 take q btw -7 and 7
                _ if (eccentricity > 0.25) => 6, // if ecc > 0.20 take q btw -6 and 6
                _ if (eccentricity > 0.2 ) => 5, // if ecc > 0.20 take q btw -5 and 5
                _ if (eccentricity > 0.15) => 4, // if ecc > 0.15 take q btw -4 and 4
                _ if (eccentricity > 0.1 ) => 3, // if ecc > 0.10 take q btw -3 and 3 
                _ if (eccentricity <= 0.1 ) => 2, // if ecc > 0.10 take q btw -2 and 2
                _ => 0,
            };
            let q_min = 7-order;
            let q_max = 2*order+1;

            let mut sum_over_p:f64=0.;
            let mut sum_over_p_s:f64=0.;

            // term m = 1
            sum_over_p=0.;
            sum_over_p_s=0.;
            let c1 = 1_f64/6_f64;
            // let c2 = -(4_f64) / (15_f64).sqrt();
            for (p, f_21p) in inclination_function_f_21p.iter().enumerate() {
                
                let tmp_p:f64 = p as f64;
                let mut sum_over_q:f64=0.;
                let mut sum_over_q_s:f64=0.;
                let sum_g_2pq = eccentricity_function_g_2pq[p];

                for (q, g_2pq) in sum_g_2pq.iter().skip(q_min).take(q_max).enumerate() {

                    let tmp_q:f64 = q as f64 -order as f64;
                    let frequ_21pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(1., tmp_p, tmp_q, spin, orbital_frequency);

                    let mut sum_over_k_term1:f64=0.;
                    let mut sum_over_k_term2:f64=0.;
                    let mut sum_over_k_term1_s:f64=0.;
                    let mut sum_over_k_term2_s:f64=0.;

                    for (k, f_30_32_k) in inclination_function_f_30p.iter().zip(inclination_function_f_32p.iter()).enumerate() {
                        // println!("|\t integer k {:?}", k);
                        let (f_30k, f_32k) = f_30_32_k; // is sum over p for each m

                        let tmp_k:f64 = k as f64;
                        let mut sum_over_j_term1:f64=0.;
                        let mut sum_over_j_term2:f64=0.;
                        let mut sum_over_j_term1_s:f64=0.;
                        let mut sum_over_j_term2_s:f64=0.;
                        let sum_g_3kj = eccentricity_function_g_3pq[k];

                        for (j, g_3kj) in sum_g_3kj.iter().skip(q_min).take(q_max).enumerate() {

                            let tmp_j:f64 = j as f64 -order as f64;

                            // let phase_beta:f64 =compute_phase_beta(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_ascending_node );
                            let phase_alpha_3:f64 =compute_phase_alpha_3(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );
                            let phase_alpha_4:f64 =compute_phase_alpha_4(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );
                            
                            let (rek2_21pq, imk2_21pq) = calculate_kaula_numbers(frequ_21pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                            
                            let cos_alpha_1:f64 = (phase_alpha_3).cos();
                            let sin_alpha_1:f64 = (phase_alpha_3).sin();
                            let cos_alpha_2:f64 = (phase_alpha_4).cos();
                            let sin_alpha_2:f64 = (phase_alpha_4).sin();

                            sum_over_j_term1 = sum_over_j_term1 + g_3kj*(sin_alpha_1 *rek2_21pq - cos_alpha_1 *imk2_21pq);
                            sum_over_j_term2 = sum_over_j_term2 + g_3kj*(sin_alpha_2 *rek2_21pq - cos_alpha_2 *imk2_21pq);
                            sum_over_j_term1_s = sum_over_j_term1_s - g_3kj *imk2_21pq;
                            sum_over_j_term2_s = sum_over_j_term2_s - g_3kj *imk2_21pq;

                        }
                        sum_over_k_term1 = sum_over_k_term1 + f_32k *sum_over_j_term1 ;
                        sum_over_k_term2 = sum_over_k_term2 + f_30k *sum_over_j_term2 ;
                        sum_over_k_term1_s = sum_over_k_term1_s + f_32k *sum_over_j_term1_s ;
                        sum_over_k_term2_s = sum_over_k_term2_s + f_30k *sum_over_j_term2_s ;
                    }
                    sum_over_q = sum_over_q + g_2pq*( c1*sum_over_k_term1 + sum_over_k_term2 ) ;
                    sum_over_q_s = sum_over_q_s + g_2pq*( c1*sum_over_k_term1_s + sum_over_k_term2_s ) ;
                }
                sum_over_p = sum_over_p + f_21p *sum_over_q ;
                sum_over_p_s = sum_over_p_s + f_21p *sum_over_q_s ;
            }
            let _term_m1:f64 = -(4_f64) / (15_f64).sqrt() * sum_over_p ;
            let _term_m1_s:f64 = -(4_f64) / (15_f64).sqrt() * sum_over_p_s ;

            // term m = 2
            sum_over_p=0.;
            sum_over_p_s=0.;
            // let c1 = -(5_f64) / (48_f64*6_f64.sqrt());
            // let c3 = 2_f64;
            for (p, f_21p) in inclination_function_f_22p.iter().enumerate() {
                
                let tmp_p:f64 = p as f64;
                let mut sum_over_q:f64=0.;
                let mut sum_over_q_s:f64=0.;
                let sum_g_2pq = eccentricity_function_g_2pq[p];

                for (q, g_2pq) in sum_g_2pq.iter().skip(q_min).take(q_max).enumerate() {

                    let tmp_q:f64 = q as f64 -order as f64;
                    let frequ_22pq:f64 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2., tmp_p, tmp_q, spin, orbital_frequency);

                    let mut sum_over_k_term1:f64=0.;
                    let mut sum_over_k_term2:f64=0.;
                    let mut sum_over_k_term1_s:f64=0.;
                    let mut sum_over_k_term2_s:f64=0.;

                    for (k, f_31_33_k) in inclination_function_f_31p.iter().zip(inclination_function_f_33p.iter()).enumerate() {

                        let (f_31k, f_33k) = f_31_33_k; // is sum over p for each m

                        let tmp_k:f64 = k as f64;
                        let mut sum_over_j_term1:f64=0.;
                        let mut sum_over_j_term2:f64=0.;
                        let mut sum_over_j_term1_s:f64=0.;
                        let mut sum_over_j_term2_s:f64=0.;
                        let sum_g_3kj = eccentricity_function_g_3pq[k];

                        for (j, g_3kj) in sum_g_3kj.iter().skip(q_min).take(q_max).enumerate() {

                            let tmp_j:f64 = j as f64 -order as f64;

                            // let phase_beta:f64 =compute_phase_beta(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_ascending_node );
                            let phase_alpha_3:f64 =compute_phase_alpha_3(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );
                            let phase_alpha_4:f64 =compute_phase_alpha_4(tmp_p, tmp_q, tmp_k, tmp_j, mean_anomaly, spin_angle, argument_perihelion, longitude_of_ascending_node, heliocentric_varphi );

                            let (rek2_22pq, imk2_22pq) = calculate_kaula_numbers(frequ_22pq, number_data_point, real_love_number_tab, imag_love_number_tab, excitation_frequency_tab );
                            
                            let cos_alpha_1:f64 = (phase_alpha_3).cos();
                            let sin_alpha_1:f64 = (phase_alpha_3).sin();
                            let cos_alpha_2:f64 = (phase_alpha_4).cos();
                            let sin_alpha_2:f64 = (phase_alpha_4).sin();

                            sum_over_j_term1 = sum_over_j_term1 + g_3kj*(sin_alpha_1 *rek2_22pq - cos_alpha_1 *imk2_22pq);
                            sum_over_j_term2 = sum_over_j_term2 + g_3kj*(sin_alpha_2 *rek2_22pq - cos_alpha_2 *imk2_22pq);
                            sum_over_j_term1_s = sum_over_j_term1_s - g_3kj *imk2_22pq;
                            sum_over_j_term2_s = sum_over_j_term2_s - g_3kj *imk2_22pq;

                        }
                        sum_over_k_term1 = sum_over_k_term1 + f_33k *sum_over_j_term1 ;
                        sum_over_k_term2 = sum_over_k_term2 + f_31k *sum_over_j_term2 ;
                        sum_over_k_term1_s = sum_over_k_term1_s + f_33k *sum_over_j_term1_s ;
                        sum_over_k_term2_s = sum_over_k_term2_s + f_31k *sum_over_j_term2_s ;
                    }
                    sum_over_q = sum_over_q + g_2pq*( sum_over_k_term1 + 2_f64*sum_over_k_term2 ) ;
                    sum_over_q_s = sum_over_q_s + g_2pq*( sum_over_k_term1_s + 2_f64*sum_over_k_term2_s ) ;
                }
                sum_over_p = sum_over_p + f_21p *sum_over_q ;
                sum_over_p_s = sum_over_p_s + f_21p *sum_over_q_s ;
            }
            let _term_m2:f64 = -(5_f64) / (48_f64*6_f64.sqrt()) *sum_over_p ;
            let _term_m2_s:f64 = -(5_f64) / (48_f64*6_f64.sqrt()) *sum_over_p_s ;

            let cste:f64 = -(G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(7)) ;
            orthogonal_force = cste *(_term_m1+_term_m2);
            orthogonal_force_secular = cste *(_term_m1_s+_term_m2_s); // cste *sum_over_p_secular;
        } // End 3D case
    return (orthogonal_force, orthogonal_force_secular);
}

// -------------------------------------------------- //
// --- Calculate tidal torque due to tidal forces --- //
pub fn calculate_torque_due_to_tides(particle: &Particle) -> Axes {

    let position_x : f64 = particle.tides.coordinates.position.x ;
    let position_y : f64 = particle.tides.coordinates.position.y ;
    let position_z : f64 = particle.tides.coordinates.position.z ;
    let tidal_force_x : f64 = particle.tides.parameters.internal.kaula_tidal_force.x ;
    let tidal_force_y : f64 = particle.tides.parameters.internal.kaula_tidal_force.y ;
    let tidal_force_z : f64 = particle.tides.parameters.internal.kaula_tidal_force.z ; 

    // Let the torque be the cross product of the radial distance vector and the tidal force vector
    let torque_due_to_tides_x = position_y *tidal_force_z
                                    -position_z *tidal_force_y;
    let torque_due_to_tides_y = position_z *tidal_force_x
                                    -position_x *tidal_force_z;
    let torque_due_to_tides_z = position_x *tidal_force_y
                                    -position_y *tidal_force_x;
    return Axes{x: torque_due_to_tides_x, y: torque_due_to_tides_y, z: torque_due_to_tides_z} ;
}

// ------------- //
// --- Tools --- //

// --- Inpuit Love numbers

// pub fn find_input_love_numbers(){
//     particle.tides.parameters.internal.num_datapoints = 1.0;
// }


// --- Calculate the associated Legendre polynomials P_{l=2}^{m}
fn calculate_associated_legendre_polynomials_p_2m(x: f64) -> [f64;3]{
    let mut p_2m = [0.0;3];
    let x_2 = x*x ;

    p_2m[0] = 0.5*(3.*x_2 -1.);
    p_2m[1] = -3.*x*((1.-x_2).sqrt()); //p_2m[1] = -3.*((1.-x_2).sqrt());
    p_2m[2] = 3.*(1.-x_2);

    return p_2m;
}

// --- Calculate the tidal excitation frequency mode \sigma_{2mpq}
fn calculate_tidal_excitation_frequency_mode_sigma_2mpq(m:f64, p:f64, q:f64, spin:f64, orbital_frequency: f64) -> f64{
    return (2.0_f64 -2.0_f64*p +q)*orbital_frequency -m*spin;
}

// --- Calculate the phases of the 2-Kaula transformed tidal forces
fn alpha_pqkj(p:f64, q:f64, k:f64, j:f64, mean_anomaly:f64, argument_perihelion:f64 ) ->f64{
    return (2_f64*p -2_f64*k +j -q) *mean_anomaly +2_f64 *(p -k) *argument_perihelion ;
}

// --- Calculate the phases of the 2-Kaula transformed tidal forces
fn compute_phase_beta(p:f64, q:f64, k:f64, j:f64, mean_anomaly:f64, spin_angle:f64, argument_perihelion:f64, longitude_ascending_node:f64, heliocentric_varphi:f64 ) ->f64{
    return (4_f64 - 2_f64*p -2_f64*k +q -j)*mean_anomaly +(4_f64 -2_f64 *(k+p)) *argument_perihelion +spin_angle -longitude_ascending_node +heliocentric_varphi;
}
fn compute_phase_alpha_1(p:f64, q:f64, k:f64, j:f64, mean_anomaly:f64, spin_angle:f64, argument_perihelion:f64, longitude_ascending_node:f64, heliocentric_varphi:f64 ) ->f64{
    return (2_f64*k -2_f64*p +q -j)*mean_anomaly +2_f64 *(k -p) *argument_perihelion +spin_angle -longitude_ascending_node -heliocentric_varphi;
}
fn compute_phase_alpha_2(p:f64, q:f64, k:f64, j:f64, mean_anomaly:f64, spin_angle:f64, argument_perihelion:f64, longitude_ascending_node:f64, heliocentric_varphi:f64 ) ->f64{
    return (2_f64*k -2_f64*p +q -j)*mean_anomaly +2_f64 *(k -p) *argument_perihelion -spin_angle +longitude_ascending_node +heliocentric_varphi;
}
fn compute_phase_alpha_3(p:f64, q:f64, k:f64, j:f64, mean_anomaly:f64, spin_angle:f64, argument_perihelion:f64, longitude_ascending_node:f64, heliocentric_varphi:f64 ) ->f64{
    return (2_f64*k -2_f64*p +q -j -1_f64)*mean_anomaly +2_f64 *(k -p -1_f64) *argument_perihelion +spin_angle -longitude_ascending_node -heliocentric_varphi;
}
fn compute_phase_alpha_4(p:f64, q:f64, k:f64, j:f64, mean_anomaly:f64, spin_angle:f64, argument_perihelion:f64, longitude_ascending_node:f64, heliocentric_varphi:f64 ) ->f64{
    return (2_f64*k -2_f64*p +q -j -1_f64)*mean_anomaly +2_f64 *(k -p -1_f64) *argument_perihelion -spin_angle +longitude_ascending_node -heliocentric_varphi;
}

// --- Find the real part and the imaginary part of the Love number associated to the excitation frequenccy wk2
fn calculate_kaula_numbers(wk2:f64, nm_data:f64, real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> (f64,f64){
    let nm_data = nm_data as usize;
    let mut line_up:usize=0;
    let mut column_up:usize=0;
    let mut line_down:usize=0;
    let mut column_down:usize=0;
    let mut in_the_file:bool=true;

    let (last_line, last_column):(usize, usize) = match nm_data%32{
        0=> (31, nm_data/32-1 ),
        _=> (nm_data%32-1, nm_data/32 ), 
    };
    let case_oneline_last_col:bool = match last_line{
        0=>true,
        1=>true,
        _=>false,
    };
    let (parity, w_k2):(bool, f64) = match (){
        _ if wk2 < 0.0 => (true, wk2.abs()),
        _ => (false, wk2.abs()),
    };

    let mut ctrl = true;
    let mut im_k2 = 0.;
    let mut re_k2 = 0.;
    if w_k2 > love_number_excitation_frequency[last_column][last_line-1] { // --- Check beyond the last value of the 2D table --- //
        if case_oneline_last_col {
            im_k2 = imaginary_part_love_number[last_column][last_line]  
                + (w_k2 - love_number_excitation_frequency[last_column][last_line]) 
                * (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column-1][31])
                /(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]);
            re_k2 = real_part_love_number[last_column][last_line]       
                + (w_k2 - love_number_excitation_frequency[last_column][last_line]) 
                * (real_part_love_number[last_column][last_line] - real_part_love_number[last_column-1][31])
                /(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]) ;
        } else {
            let delta = (w_k2 - love_number_excitation_frequency[last_column][last_line-1]) 
                / (love_number_excitation_frequency[last_column][last_line-1] - love_number_excitation_frequency[last_column][last_line-2]);
            im_k2 = delta
                *( imaginary_part_love_number[last_column][last_line-1] -imaginary_part_love_number[last_column][last_line-2] ) 
                +imaginary_part_love_number[last_column][last_line-1];
            re_k2 = delta
                *( real_part_love_number[last_column][last_line-1] -real_part_love_number[last_column][last_line-2] )
                +real_part_love_number[last_column][last_line-1];
            if im_k2<0. {
                im_k2 = imaginary_part_love_number[last_column][last_line]
                    + (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column][last_line-1])/2.;
            }
            if re_k2<0. {
                re_k2 = real_part_love_number[last_column][last_line] 
                    + (real_part_love_number[last_column][last_line] - real_part_love_number[last_column][last_line-1])/2.;
            }
        }
        column_up = last_column;
        line_up = last_line;
        column_down = 0;
        line_down = 0;
        in_the_file = false;
    } else if w_k2 > love_number_excitation_frequency[last_column-1][31] { // --- If is in the last column --- //
        if  w_k2 < love_number_excitation_frequency[last_column][0] { // --- Check is the value is btw the two last column --- //
            let delta = (w_k2 - love_number_excitation_frequency[last_column-1][31]) 
                / (love_number_excitation_frequency[last_column][0] - love_number_excitation_frequency[last_column-1][31]);
            im_k2 = (1.-delta)* imaginary_part_love_number[last_column-1][31]
                +delta*(imaginary_part_love_number[last_column][0]);
            re_k2 = (1.-delta)* real_part_love_number[last_column-1][31]
                +delta*(real_part_love_number[last_column][0]);
            column_up  = last_column;
            line_up    = 0;
            column_down  = last_column-1;
            line_down    = 31;
        } else { // --- Check into the last column --- //
            for frequency2 in 0..last_line+1{ // --- loop for from 00 to last line-1 --- //
                if  love_number_excitation_frequency[last_column][frequency2] >= w_k2 {
                    if  love_number_excitation_frequency[last_column][frequency2] == w_k2 {
                        im_k2 = imaginary_part_love_number[last_column][frequency2];
                        re_k2 = real_part_love_number[last_column][frequency2];
                        column_up  = last_column;
                        line_up    = frequency2;
                        column_down  = last_column;
                        line_down    = frequency2;
                        break;
                    } else {
                        let delta = (w_k2 - love_number_excitation_frequency[last_column][frequency2-1]) 
                            / (love_number_excitation_frequency[last_column][frequency2] -love_number_excitation_frequency[last_column][frequency2-1]);
                        im_k2 = (1.-delta)* imaginary_part_love_number[last_column][frequency2-1]
                            +delta*(imaginary_part_love_number[last_column][frequency2]);
                        re_k2 = (1.-delta)* real_part_love_number[last_column][frequency2-1] 
                            +delta*(real_part_love_number[last_column][frequency2]);
                        column_up  = last_column;
                        line_up    = frequency2;
                        column_down  = last_column;
                        line_down    = frequency2-1;
                        break;
                    }
                }
            }
        } // --- End check last column --- //
    } else { // --- Check All other column --- //
        for frequency1 in 0..last_column{ // --- Loop into the column without the last one (loop from 0 to last_column-1) --- //
            if ctrl && love_number_excitation_frequency[frequency1][31] >= w_k2 { // --- If the last value of the col frequency1 is GT w2, w2 is in this col --- //
                if ctrl && love_number_excitation_frequency[frequency1][0] > w_k2 { // --- Check is the fisrt value is GT w2, then the value is btw the two lat colunm --- //
                    let delta = (w_k2 - love_number_excitation_frequency[frequency1-1][31]) 
                        / (love_number_excitation_frequency[frequency1][0] - love_number_excitation_frequency[frequency1-1][31]);
                    im_k2 = (1.-delta)* imaginary_part_love_number[frequency1-1][31]
                        +delta*(imaginary_part_love_number[frequency1][0]);
                    re_k2 = (1.-delta)* real_part_love_number[frequency1-1][31]
                        +delta*(real_part_love_number[frequency1][0]);
                    column_up  = frequency1;
                    line_up    = 0;
                    column_down  = frequency1-1;
                    line_down    = 31;
                    ctrl = false;
                } else if ctrl {
                    for frequency2 in 0..love_number_excitation_frequency.len(){ // --- loop into the 32 elements of each line --- // 
                        if ctrl && love_number_excitation_frequency[frequency1][frequency2] >= w_k2 { 
                            if  love_number_excitation_frequency[frequency1][frequency2] == w_k2 {
                                im_k2 = imaginary_part_love_number[frequency1][frequency2];
                                re_k2 = real_part_love_number[frequency1][frequency2];
                                column_up  = frequency1;
                                line_up    = frequency2;
                                column_down  = frequency1;
                                line_down    = frequency2;
                                ctrl = false;
                                break;
                            }
                            else if ctrl { 
                                let delta = (w_k2 - love_number_excitation_frequency[frequency1][frequency2-1]) 
                                    / (love_number_excitation_frequency[frequency1][frequency2] -love_number_excitation_frequency[frequency1][frequency2-1]);
                                im_k2 = (1.-delta)* imaginary_part_love_number[frequency1][frequency2-1] 
                                    +delta*(imaginary_part_love_number[frequency1][frequency2]);
                                re_k2 = (1.-delta)* real_part_love_number[frequency1][frequency2-1]
                                    +delta*(real_part_love_number[frequency1][frequency2]);
                                column_up  = frequency1;
                                line_up    = frequency2;
                                column_down  = frequency1;
                                line_down    = frequency2-1;
                                ctrl = false;
                                break;
                            }
                        }
                        if !ctrl {break;}
                    }
                }
            }
            if !ctrl {break;}
        }
    }
    if parity {
        im_k2 = -im_k2;
    }
    if ((w_k2!=0. && im_k2==0.) || w_k2>love_number_excitation_frequency[column_up][line_up] || w_k2<love_number_excitation_frequency[column_down][line_down]) && in_the_file {
        panic!("=== Error  in tide.rs :: fail interpolation Imk2 === found btw column&line {} {} < {} {} in frequency: {} < {} < {} \t the ImK2 {}", 
        column_down, line_down, column_up, line_up, love_number_excitation_frequency[column_down][line_down], w_k2, love_number_excitation_frequency[column_up][line_up], im_k2);
    }
    return (-re_k2, im_k2)
}

// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_2mp(inclination: f64) -> [[f64; 3]; 3]{

    let mut inclination_function = [[0.;3];3];
    // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3
    let sin_incl = inclination.sin();
    let sin_incl_2 = sin_incl*sin_incl;
    let cos_incl = inclination.cos();

    inclination_function[0][0] = -(3./8.)*sin_incl_2;
    inclination_function[0][1] = (3./4.)*sin_incl_2 -(1./2.);
    inclination_function[0][2] = -(3./8.)*sin_incl_2;

    inclination_function[1][0] = (3./4.)*sin_incl*(1. +cos_incl);
    inclination_function[1][1] = -(3./2.)*sin_incl*cos_incl;
    inclination_function[1][2] = (3./4.)*sin_incl*(cos_incl -1.);

    inclination_function[2][0] = (3./4.)*(1. +cos_incl).powi(2);
    inclination_function[2][1] = (3./2.)*sin_incl.powi(2);
    inclination_function[2][2] = (3./4.)*(1. -cos_incl).powi(2);
    
    return inclination_function;
}

// F2mp
// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_20p(inclination: f64) -> [f64; 3]{
    let mut inclination_function = [0.;3];
    // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3
    let sin_incl = inclination.sin();
    let sin_incl_2 = sin_incl*sin_incl;

    inclination_function[0] = -(3./8.)*sin_incl_2;
    inclination_function[1] = (3./4.)*sin_incl_2 -(1./2.);
    inclination_function[2] = -(3./8.)*sin_incl_2;

    return inclination_function;
}
// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_21p(inclination: f64) -> [f64; 3]{
    let mut inclination_function = [0.;3];
    let sin_incl = inclination.sin();
    let cos_incl = inclination.cos();

    inclination_function[0] = (3./4.)*sin_incl*(1. +cos_incl);
    inclination_function[1] = -(3./2.)*sin_incl*cos_incl;
    inclination_function[2] = (3./4.)*sin_incl*(cos_incl -1.);

    return inclination_function;
}
// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_22p(inclination: f64) -> [f64; 3]{
    let mut inclination_function = [0.;3];
    let sin_incl = inclination.sin();
    let cos_incl = inclination.cos();

    inclination_function[0] = (3./4.)*(1. +cos_incl).powi(2);
    inclination_function[1] = (3./2.)*sin_incl.powi(2);
    inclination_function[2] = (3./4.)*(1. -cos_incl).powi(2);

    return inclination_function;
}
// F3pmp
// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_30p(inclination: f64) -> [f64; 4]{
    let mut inclination_function = [0.;4];
    let sin_incl = inclination.sin();
    let sin_3_factor_1 = (5./16.)*sin_incl.powi(3);
    let sin_factor_2 = (3./4.)*sin_incl;

    inclination_function[0] = - sin_3_factor_1;
    inclination_function[1] = sin_3_factor_1 -sin_factor_2;
    inclination_function[2] = - sin_3_factor_1 +sin_factor_2;
    inclination_function[3] = - inclination_function[0];

    return inclination_function;
}
// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_31p(inclination: f64) -> [f64; 4]{
    let mut inclination_function = [0.;4];
    let sin_incl = inclination.sin();
    let sin_2_factor_1 = (15./16.)*sin_incl.powi(2);
    let sin_2_factor_43 = sin_2_factor_1 - (3./4.);
    let cos_incl = inclination.cos();
    let cos_incl_factor_1 = 1.0+3.0*cos_incl;
    let cos_incl_factor_2 = 1.0-3.0*cos_incl;

    inclination_function[0] = - sin_2_factor_1 *cos_incl_factor_1;
    inclination_function[1] = sin_2_factor_43 *cos_incl_factor_1;
    inclination_function[2] = sin_2_factor_43 *cos_incl_factor_2;
    inclination_function[3] = - sin_2_factor_1 *cos_incl_factor_2;

    return inclination_function;
}
// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_32p(inclination: f64) -> [f64; 4]{
    let mut inclination_function = [0.;4];
    let sin_incl = inclination.sin();
    let sin_2_factor_1 = (15./8.)*sin_incl.powi(2);
    let cos_incl = inclination.cos();
    let cos_incl_2 = cos_incl.powi(2);

    inclination_function[0] = sin_2_factor_1 *(1.0+cos_incl).powi(2);
    inclination_function[1] = sin_2_factor_1 *(1.0-2.0*cos_incl-3.0*cos_incl_2);
    inclination_function[2] = - sin_2_factor_1 *(1.0+2.0*cos_incl-3.0*cos_incl_2);
    inclination_function[3] = - sin_2_factor_1 *(1.0-cos_incl).powi(2);

    return inclination_function;
}
// --- The inclination function of the Kaula 1961 developpement Flmp(i)
fn calculate_inclination_function_f_33p(inclination: f64) -> [f64; 4]{
    let mut inclination_function = [0.;4];
    let sin_incl = inclination.sin();
    let sin_2_factor_1 = (48./8.)*sin_incl.powi(2);
    let cos_incl = inclination.cos();
    let cos_incl_factor_1 = 1.0+cos_incl;
    let cos_incl_factor_2 = 1.0-cos_incl;

    inclination_function[0] = (15./8.)*cos_incl_factor_1.powi(3);
    inclination_function[1] = sin_2_factor_1 *cos_incl_factor_1;
    inclination_function[2] = sin_2_factor_1 *cos_incl_factor_2;
    inclination_function[3] = (15./8.)*cos_incl_factor_2.powi(3);

    return inclination_function;
}


// --- The eccentricity function of the Kaula 1961 developpement Glpq(e) from the Hensen coeff [?] (see Cayley table [?] )
fn calculate_eccentricty_function_g_2pq(eccentricity: f64) -> [[f64; 15]; 3]{

    let mut eccentricity_function = [[0.0;15];3];
    // in this notation:
    // [0] is -7
    // [1] is -6
    // [2] is -5
    // [3] is -4
    // [4] is -3
    // [5] is -2
    // [6] is -1
    // [7] is 0
    // [8] is 1
    // [9] is 2
    // [10] is 3   
    // [11] is 4
    // [12] is 5
    // [13] is 6
    // [14] is 7

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc*ecc;
    let ecc_3:f64=ecc_2*ecc;
    let ecc_4:f64=ecc_3*ecc;
    let ecc_5:f64=ecc_4*ecc;
    let ecc_6:f64=ecc_5*ecc;
    let ecc_7:f64=ecc_6*ecc;

    eccentricity_function[0][0] = (15625./129024.)*ecc_7 ;
    eccentricity_function[0][1] = (4./45.)*ecc_6 ;
    eccentricity_function[0][2] = (81./1280.)*ecc_5     +(81./2048.)*ecc_7 ;
    eccentricity_function[0][3] = (1./24.)*ecc_4        +(7./240.)*ecc_6 ;
    eccentricity_function[0][4] = (1./48.)*ecc_3        +(11./768.)*ecc_5       +(313./30720.)*ecc_7 ;
    eccentricity_function[0][5] = 0.;
    eccentricity_function[0][6] = -0.5*ecc              +(1./16.)*ecc_3         -(5./384.)*ecc_5        -(143./18432.)*ecc_7 ;
    eccentricity_function[0][7] = 1.                    -(5./2.)*ecc_2          +(13./16.)*ecc_4        -(35./288.)*ecc_6 ;
    eccentricity_function[0][8] = (7./2.)*ecc           -(123./16.)*ecc_3       +(489./128.)*ecc_5      -(1763./2048.)*ecc_7 ;
    eccentricity_function[0][9] = (17./2.)*ecc_2        -(115./16.)*ecc_4       +(601./48.)*ecc_6 ;
    eccentricity_function[0][10] = (845./48.)*ecc_3     -(32525./768.)*ecc_5    +(208225./6144.)*ecc_7 ;
    eccentricity_function[0][11] = (533./16.)*ecc_4     -(13827./160.)*ecc_6 ;
    eccentricity_function[0][12] = (228347./3840.)*ecc_5 -(3071075./18432.)*ecc_7 ;
    eccentricity_function[0][13] = (73369./720.)*ecc_6 ;
    eccentricity_function[0][14] = (12144273./71680.)*ecc_7 ;

    eccentricity_function[1][0] = (432091./30720.)*ecc_7 ;
    eccentricity_function[1][1] = (3167./320.)*ecc_6 ;
    eccentricity_function[1][2] = (1773./256.)*ecc_5    -(4987./6144.)*ecc_7 ;
    eccentricity_function[1][3] = (77./16.)*ecc_4       +(129./160.)*ecc_6 ;
    eccentricity_function[1][4] = (53./16.)*ecc_3       +(393./256.)*ecc_5      +(24753./10240.)*ecc_7;
    eccentricity_function[1][5] = (9./4.)*ecc_2         +(7./4.)*ecc_4          +(141./64.)*ecc_6 ;
    eccentricity_function[1][6] = (3./2.)*ecc           +(27./16.)*ecc_3        +(261./128.)*ecc_5      +(14309./6144.)*ecc_7 ;
    eccentricity_function[1][7] = (1.-ecc_2).powf(-3./2.) ;
    eccentricity_function[1][8] = eccentricity_function[1][6] ;
    eccentricity_function[1][9] = eccentricity_function[1][5] ;
    eccentricity_function[1][10] = eccentricity_function[1][4] ;
    eccentricity_function[1][11] = eccentricity_function[1][3] ;
    eccentricity_function[1][12] = eccentricity_function[1][2] ;
    eccentricity_function[1][13] = eccentricity_function[1][1] ;
    eccentricity_function[1][14] = eccentricity_function[1][0] ;

    eccentricity_function[2][0] = eccentricity_function[0][14];
    eccentricity_function[2][1] = eccentricity_function[0][13];
    eccentricity_function[2][2] = eccentricity_function[0][12];
    eccentricity_function[2][3] = eccentricity_function[0][11];
    eccentricity_function[2][4] = eccentricity_function[0][10];
    eccentricity_function[2][5] = eccentricity_function[0][9] ;
    eccentricity_function[2][6] = eccentricity_function[0][8] ;
    eccentricity_function[2][7] = eccentricity_function[0][7] ;
    eccentricity_function[2][8] = eccentricity_function[0][6] ;
    eccentricity_function[2][9] = eccentricity_function[0][5] ;
    eccentricity_function[2][10] =eccentricity_function[0][4] ;
    eccentricity_function[2][11] =eccentricity_function[0][3] ;
    eccentricity_function[2][12] =eccentricity_function[0][2] ;
    eccentricity_function[2][13] =eccentricity_function[0][1] ;
    eccentricity_function[2][14] =eccentricity_function[0][0] ;

    return eccentricity_function;
}

// --- The eccentricity function of the Kaula 1961 developpement G20q(e) from the Hensen coeff [?] (see Cayley table [?] )
fn calculate_eccentricty_function_g_20q(eccentricity: f64) -> [f64; 15]{

    let mut eccentricity_function = [0.0;15];
    // in this notation:
    // [0] is -7
    // [1] is -6
    // [2] is -5
    // [3] is -4
    // [4] is -3
    // [5] is -2
    // [6] is -1
    // [7] is 0
    // [8] is 1
    // [9] is 2
    // [10] is 3   
    // [11] is 4
    // [12] is 5
    // [13] is 6
    // [14] is 7

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc*ecc;
    let ecc_3:f64=ecc_2*ecc;
    let ecc_4:f64=ecc_3*ecc;
    let ecc_5:f64=ecc_4*ecc;
    let ecc_6:f64=ecc_5*ecc;
    let ecc_7:f64=ecc_6*ecc;

    eccentricity_function[0] = (15625./129024.)*ecc_7 ;
    eccentricity_function[1] = (4./45.)*ecc_6 ;
    eccentricity_function[2] = (81./1280.)*ecc_5     +(81./2048.)*ecc_7 ;
    eccentricity_function[3] = (1./24.)*ecc_4        +(7./240.)*ecc_6 ;
    eccentricity_function[4] = (1./48.)*ecc_3        +(11./768.)*ecc_5       +(313./30720.)*ecc_7 ;
    eccentricity_function[5] = 0.;
    eccentricity_function[6] = -0.5*ecc              +(1./16.)*ecc_3         -(5./384.)*ecc_5        -(143./18432.)*ecc_7 ;
    eccentricity_function[7] = 1.                    -(5./2.)*ecc_2          +(13./16.)*ecc_4        -(35./288.)*ecc_6 ;
    eccentricity_function[8] = (7./2.)*ecc           -(123./16.)*ecc_3       +(489./128.)*ecc_5      -(1763./2048.)*ecc_7 ;
    eccentricity_function[9] = (17./2.)*ecc_2        -(115./16.)*ecc_4       +(601./48.)*ecc_6 ;
    eccentricity_function[10] = (845./48.)*ecc_3     -(32525./768.)*ecc_5    +(208225./6144.)*ecc_7 ;
    eccentricity_function[11] = (533./16.)*ecc_4     -(13827./160.)*ecc_6 ;
    eccentricity_function[12] = (228347./3840.)*ecc_5 -(3071075./18432.)*ecc_7 ;
    eccentricity_function[13] = (73369./720.)*ecc_6 ;
    eccentricity_function[14] = (12144273./71680.)*ecc_7 ;

    return eccentricity_function;
}

// --- The eccentricity function of the Kaula 1961 developpement G21q(e) from the Hensen coeff [?] (see Cayley table [?] )
fn calculate_eccentricty_function_g_21q(eccentricity: f64) -> [f64; 15]{

    let mut eccentricity_function = [0.0;15];
    // in this notation:
    // [0] is -7
    // [1] is -6
    // [2] is -5
    // [3] is -4
    // [4] is -3
    // [5] is -2
    // [6] is -1
    // [7] is 0
    // [8] is 1
    // [9] is 2
    // [10] is 3   
    // [11] is 4
    // [12] is 5
    // [13] is 6
    // [14] is 7

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc*ecc;
    let ecc_3:f64=ecc_2*ecc;
    let ecc_4:f64=ecc_3*ecc;
    let ecc_5:f64=ecc_4*ecc;
    let ecc_6:f64=ecc_5*ecc;
    let ecc_7:f64=ecc_6*ecc;

    eccentricity_function[0] = (432091./30720.)*ecc_7 ;
    eccentricity_function[1] = (3167./320.)*ecc_6 ;
    eccentricity_function[2] = (1773./256.)*ecc_5    -(4987./6144.)*ecc_7 ;
    eccentricity_function[3] = (77./16.)*ecc_4       +(129./160.)*ecc_6 ;
    eccentricity_function[4] = (53./16.)*ecc_3       +(393./256.)*ecc_5      +(24753./10240.)*ecc_7;
    eccentricity_function[5] = (9./4.)*ecc_2         +(7./4.)*ecc_4          +(141./64.)*ecc_6 ;
    eccentricity_function[6] = (3./2.)*ecc           +(27./16.)*ecc_3        +(261./128.)*ecc_5      +(14309./6144.)*ecc_7 ;
    eccentricity_function[7] = (1.-ecc_2).powf(-3./2.) ;
    eccentricity_function[8] = eccentricity_function[6] ;
    eccentricity_function[9] = eccentricity_function[5] ;
    eccentricity_function[10] = eccentricity_function[4] ;
    eccentricity_function[11] = eccentricity_function[3] ;
    eccentricity_function[12] = eccentricity_function[2] ;
    eccentricity_function[13] = eccentricity_function[1] ;
    eccentricity_function[14] = eccentricity_function[0] ;

    return eccentricity_function;
}

// --- The eccentricity function of the Kaula 1961 developpement Glpq(e) from the Hensen coeff [?] (see Cayley table [?] )
fn calculate_eccentricty_function_g_3pq(eccentricity: f64) -> [[f64; 15]; 4]{

    let mut eccentricity_function = [[0.0;15];4];
    // in this notation:
    // [0] is -7
    // [1] is -6
    // [2] is -5
    // [3] is -4
    // [4] is -3
    // [5] is -2
    // [6] is -1
    // [7] is 0
    // [8] is 1
    // [9] is 2
    // [10] is 3   
    // [11] is 4
    // [12] is 5
    // [13] is 6
    // [14] is 7

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc.powi(2);
    let ecc_3:f64=ecc.powi(3);
    let ecc_4:f64=ecc.powi(4);
    let ecc_5:f64=ecc.powi(5);
    let ecc_6:f64=ecc.powi(6);
    let ecc_7:f64=ecc.powi(7);

    eccentricity_function[0][0] = (8./315.)*ecc_7 ;
    eccentricity_function[0][1] = (81./5120.)*ecc_6 ;
    eccentricity_function[0][2] = (1./120.)*ecc_5 +(13./1440.)*ecc_7 ;
    eccentricity_function[0][3] = (1./384.)*ecc_4 +(1./384.)*ecc_6 ;
    eccentricity_function[0][4] = 0.;
    eccentricity_function[0][5] = (1./8.)*ecc_2 +(1./48.)*ecc_4 +(55./3072.)*ecc_6 ;
    eccentricity_function[0][6] = -ecc +(5./4.)*ecc_3 -(7./48.)*ecc_5 +(23./288.)*ecc_7 ;
    eccentricity_function[0][7] = 1. -6.*ecc_2 +(423./64.)*ecc_4 -(125./64.)*ecc_6 ;
    eccentricity_function[0][8] = 5.*ecc -22.*ecc_3 +(607./24.)*ecc_5 -(98./9.)*ecc_7 ;
    eccentricity_function[0][9] = (127./8.)*ecc_2 -(3065./48.)*ecc_4 +(243805./3072.)*ecc_6 ;
    eccentricity_function[0][10] = (163./4.)*ecc_3 -(2577./16.)*ecc_5 +(1089./5.)*ecc_7 ;
    eccentricity_function[0][11] = (35413./384.)*ecc_4 -(709471./1920.)*ecc_6 ;
    eccentricity_function[0][12] = (23029./120.)*ecc_5 -(35614./45.)*ecc_7 ;
    eccentricity_function[0][13] = (385095./1024.)*ecc_6 ;
    eccentricity_function[0][14] = (44377./63.)*ecc_7 ;

    eccentricity_function[1][0] = (16337./2240.)*ecc_7 ;
    eccentricity_function[1][1] = (48203./9240.)*ecc_6 ;
    eccentricity_function[1][2] = (899./240.)*ecc_5 +(2441./480.)*ecc_7 ;
    eccentricity_function[1][3] = (343./128.)*ecc_4 +(2819./640.)*ecc_6 ;
    eccentricity_function[1][4] = (23./12.)*ecc_3  +(89./24.)*ecc_5 +(5663./960.)*ecc_7 ;
    eccentricity_function[1][5] = (11./8.)*ecc_2 +(49./16.)*ecc_4 +(15665./3072.)*ecc_6 ;
    eccentricity_function[1][6] = ecc*(1.-ecc_2).powf(-5./2.) ;
    eccentricity_function[1][7] = 1. +(1./2.)*ecc_2 +(239./64.)*ecc_4 -(3323./576.)*ecc_6 ;
    eccentricity_function[1][8] = 3.*ecc +(11./4.)*ecc_3 +(245./48.)*ecc_5 +(463./64.)*ecc_7 ;
    eccentricity_function[1][9] = (53./8.)*ecc_2 +(39./16.)*ecc_4 +(7041./1024.)*ecc_6 ;
    eccentricity_function[1][10] = (163./4.)*ecc_3 -(2577./16.)*ecc_5 +(1089./5.)*ecc_7 ;
    eccentricity_function[1][11] = (35413./384.)*ecc_4 -(709471./1920.)*ecc_6 ;
    eccentricity_function[1][12] = (23029./120.)*ecc_5 -(35614./45.)*ecc_7 ;
    eccentricity_function[1][13] = (385095./1024.)*ecc_6 ;
    eccentricity_function[1][14] = (44377./63.)*ecc_7 ;

    eccentricity_function[2][0] = eccentricity_function[1][14];
    eccentricity_function[2][1] = eccentricity_function[1][13];
    eccentricity_function[2][2] = eccentricity_function[1][12];
    eccentricity_function[2][3] = eccentricity_function[1][11];
    eccentricity_function[2][4] = eccentricity_function[1][10];
    eccentricity_function[2][5] = eccentricity_function[1][9] ;
    eccentricity_function[2][6] = eccentricity_function[1][8] ;
    eccentricity_function[2][7] = eccentricity_function[1][7] ;
    eccentricity_function[2][8] = eccentricity_function[1][6] ;
    eccentricity_function[2][9] = eccentricity_function[1][5] ;
    eccentricity_function[2][10] =eccentricity_function[1][4] ;
    eccentricity_function[2][11] =eccentricity_function[1][3] ;
    eccentricity_function[2][12] =eccentricity_function[1][2] ;
    eccentricity_function[2][13] =eccentricity_function[1][1] ;
    eccentricity_function[2][14] =eccentricity_function[1][0] ;

    eccentricity_function[3][0] = eccentricity_function[0][14];
    eccentricity_function[3][1] = eccentricity_function[0][13];
    eccentricity_function[3][2] = eccentricity_function[0][12];
    eccentricity_function[3][3] = eccentricity_function[0][11];
    eccentricity_function[3][4] = eccentricity_function[0][10];
    eccentricity_function[3][5] = eccentricity_function[0][9] ;
    eccentricity_function[3][6] = eccentricity_function[0][8] ;
    eccentricity_function[3][7] = eccentricity_function[0][7] ;
    eccentricity_function[3][8] = eccentricity_function[0][6] ;
    eccentricity_function[3][9] = eccentricity_function[0][5] ;
    eccentricity_function[3][10] =eccentricity_function[0][4] ;
    eccentricity_function[3][11] =eccentricity_function[0][3] ;
    eccentricity_function[3][12] =eccentricity_function[0][2] ;
    eccentricity_function[3][13] =eccentricity_function[0][1] ;
    eccentricity_function[3][14] =eccentricity_function[0][0] ;

    return eccentricity_function;
}

// --- The eccentricity function of the Kaula 1961 developpement G30q(e) from the Hensen coeff [?] (see Cayley table [?] )
fn calculate_eccentricty_function_g_30q(eccentricity: f64) -> [f64; 15]{

    let mut eccentricity_function = [0.0;15];
    // in this notation:
    // [0] is -7
    // [1] is -6
    // [2] is -5
    // [3] is -4
    // [4] is -3
    // [5] is -2
    // [6] is -1
    // [7] is 0
    // [8] is 1
    // [9] is 2
    // [10] is 3   
    // [11] is 4
    // [12] is 5
    // [13] is 6
    // [14] is 7

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc.powi(2);
    let ecc_3:f64=ecc.powi(3);
    let ecc_4:f64=ecc.powi(4);
    let ecc_5:f64=ecc.powi(5);
    let ecc_6:f64=ecc.powi(6);
    let ecc_7:f64=ecc.powi(7);

    eccentricity_function[0] = (8./315.)*ecc_7 ;
    eccentricity_function[1] = (81./5120.)*ecc_6 ;
    eccentricity_function[2] = (1./120.)*ecc_5 +(13./1440.)*ecc_7 ;
    eccentricity_function[3] = (1./384.)*ecc_4 +(1./384.)*ecc_6 ;
    eccentricity_function[4] = 0.;
    eccentricity_function[5] = (1./8.)*ecc_2 +(1./48.)*ecc_4 +(55./3072.)*ecc_6 ;
    eccentricity_function[6] = -ecc +(5./4.)*ecc_3 -(7./48.)*ecc_5 +(23./288.)*ecc_7 ;
    eccentricity_function[7] = 1. -6.*ecc_2 +(423./64.)*ecc_4 -(125./64.)*ecc_6 ;
    eccentricity_function[8] = 5.*ecc -22.*ecc_3 +(607./24.)*ecc_5 -(98./9.)*ecc_7 ;
    eccentricity_function[9] = (127./8.)*ecc_2 -(3065./48.)*ecc_4 +(243805./3072.)*ecc_6 ;
    eccentricity_function[10] = (163./4.)*ecc_3 -(2577./16.)*ecc_5 +(1089./5.)*ecc_7 ;
    eccentricity_function[11] = (35413./384.)*ecc_4 -(709471./1920.)*ecc_6 ;
    eccentricity_function[12] = (23029./120.)*ecc_5 -(35614./45.)*ecc_7 ;
    eccentricity_function[13] = (385095./1024.)*ecc_6 ;
    eccentricity_function[14] = (44377./63.)*ecc_7 ;

    return eccentricity_function;
}
// --- The eccentricity function of the Kaula 1961 developpement G31q(e) from the Hensen coeff [?] (see Cayley table [?] )
fn calculate_eccentricty_function_g_31q(eccentricity: f64) -> [f64; 15]{

    let mut eccentricity_function = [0.0;15];
    // in this notation:
    // [0] is -7
    // [1] is -6
    // [2] is -5
    // [3] is -4
    // [4] is -3
    // [5] is -2
    // [6] is -1
    // [7] is 0
    // [8] is 1
    // [9] is 2
    // [10] is 3   
    // [11] is 4
    // [12] is 5
    // [13] is 6
    // [14] is 7

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc.powi(2);
    let ecc_3:f64=ecc.powi(3);
    let ecc_4:f64=ecc.powi(4);
    let ecc_5:f64=ecc.powi(5);
    let ecc_6:f64=ecc.powi(6);
    let ecc_7:f64=ecc.powi(7);

    eccentricity_function[0] = (16337./2240.)*ecc_7 ;
    eccentricity_function[1] = (48203./9240.)*ecc_6 ;
    eccentricity_function[2] = (899./240.)*ecc_5 +(2441./480.)*ecc_7 ;
    eccentricity_function[3] = (343./128.)*ecc_4 +(2819./640.)*ecc_6 ;
    eccentricity_function[4] = (23./12.)*ecc_3  +(89./24.)*ecc_5 +(5663./960.)*ecc_7 ;
    eccentricity_function[5] = (11./8.)*ecc_2 +(49./16.)*ecc_4 +(15665./3072.)*ecc_6 ;
    eccentricity_function[6] = ecc*(1.-ecc_2).powf(-5./2.) ;
    eccentricity_function[7] = 1. +(1./2.)*ecc_2 +(239./64.)*ecc_4 -(3323./576.)*ecc_6 ;
    eccentricity_function[8] = 3.*ecc +(11./4.)*ecc_3 +(245./48.)*ecc_5 +(463./64.)*ecc_7 ;
    eccentricity_function[9] = (53./8.)*ecc_2 +(39./16.)*ecc_4 +(7041./1024.)*ecc_6 ;
    eccentricity_function[10] = (163./4.)*ecc_3 -(2577./16.)*ecc_5 +(1089./5.)*ecc_7 ;
    eccentricity_function[11] = (35413./384.)*ecc_4 -(709471./1920.)*ecc_6 ;
    eccentricity_function[12] = (23029./120.)*ecc_5 -(35614./45.)*ecc_7 ;
    eccentricity_function[13] = (385095./1024.)*ecc_6 ;
    eccentricity_function[14] = (44377./63.)*ecc_7 ;

    return eccentricity_function;
}
