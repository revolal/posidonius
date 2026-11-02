use serde::{Deserialize, Serialize};
use serde_big_array::BigArray;

// Tabulated Love-number spectrum of a body. The spectrum is prepared and read by kaula_tide:
// unsigned for the planetary tide (prepare_spectrum, love_number_at), signed for the stellar tide
// (prepare_signed_spectrum, love_number_signed).
#[derive(Serialize, Copy, Deserialize, PartialEq, Debug, Clone)]
#[serde(default, deny_unknown_fields)]
pub struct LoveNumber {
    #[serde(with = "BigArray")]
    spectrum_excitation_frequency: [f64; 1024],
    #[serde(with = "BigArray")]
    spectrum_real_part: [f64; 1024],
    #[serde(with = "BigArray")]
    spectrum_imaginary_part: [f64; 1024],

    // If we are simulating stellar tide, then we need the value of the spectrum spin rate,
    // which was an input value of stellar spin rate from the computation of the spectrum.
    stellar_spectrum_spin_rate: Option<f64>, // Specified the initial spin rate of the star, or `None` if planetary tide.
}

impl Default for LoveNumber {
    fn default() -> Self {
        Self {
            spectrum_excitation_frequency: [0.; 1024],
            spectrum_real_part: [0.; 1024],
            spectrum_imaginary_part: [0.; 1024],
            stellar_spectrum_spin_rate: None,
        }
    }
}

impl LoveNumber {
    pub fn new_from(
        spectrum_excitation_frequency: [f64; 1024],
        spectrum_real_part: [f64; 1024],
        spectrum_imaginary_part: [f64; 1024],
        stellar_spectrum_spin_rate: Option<f64>,
    ) -> Self {
        Self {
            spectrum_excitation_frequency,
            spectrum_real_part,
            spectrum_imaginary_part,
            stellar_spectrum_spin_rate,
        }
    }

    /// Read access to the tabulated spectrum: (chi in rad/s, rek2, imk2), used by kaula_tide.
    pub(crate) fn spectrum(&self) -> (&[f64; 1024], &[f64; 1024], &[f64; 1024]) {
        (
            &self.spectrum_excitation_frequency,
            &self.spectrum_real_part,
            &self.spectrum_imaginary_part,
        )
    }

    /// Spin rate of the body for which the spectrum was computed (rad/s), `None` for no rescaling.
    pub(crate) fn spectrum_spin_rate(&self) -> Option<f64> {
        self.stellar_spectrum_spin_rate
    }
}
