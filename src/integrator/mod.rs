mod leapfrog;
mod ias15;
pub mod whfast;
pub mod output;

pub use self::leapfrog::*;
pub use self::ias15::*;
pub use self::whfast::WHFast;

use std::io::{BufWriter};
use std::fs::File;
use std::path::Path;
use std::any::Any;


pub trait Integrator {
    fn as_any(&self) -> &dyn Any;
    fn get_n_historic_snapshots(&self) -> usize;
    // fn set_n_historic_snapshots(&self,new_n_historic_snapshots:usize);
    fn get_new_n_historic_snapshots(&self, historic_snapshot_period: f64) -> usize;
    fn get_n_particles(&self) -> usize;
    fn get_current_time(&self) -> f64;
    fn set_time_limit(&mut self, time_limit: f64);
    fn set_snapshot_periods(&mut self, historic_snapshot_period: f64, recovery_snapshot_period: f64);
    fn initialize_physical_values(&mut self);
    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool, just_resumed:bool, new_n_historic_snapshots:usize) -> Result<bool, String>;
    fn write_recovery_snapshot(&mut self, snapshot_path: &Path, universe_history_writer: &mut BufWriter<File>);
}


