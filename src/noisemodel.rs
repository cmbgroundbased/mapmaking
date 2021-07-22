//! Given the total optical load `T_{sys}` and the `fs`, `NoiseModel` returns a TOD of white noise

use rand_distr::{Distribution, Normal};
pub struct NoiseModel {
    _sigma_k: f32,
    _f_knee: f32,
    _slope: f32,
    _seed: i32,
    _samples: usize,
    noise_tod: Vec<f32>,
}


impl NoiseModel {
    pub fn new(t_sys_k:f32, bandwidth_hz: f32, tau_s: f32, f_knee:f32, slope: f32, seed:i32, samples:usize) -> Self {
        let den = bandwidth_hz * tau_s;
        let sigma_k = t_sys_k / den.sqrt();
      
        // Starting with only white noise
        let mut _rng = rand::thread_rng();
        let normal = Normal::new(0.0, sigma_k).unwrap();
        let white_noise_tod: Vec<f32> = normal.sample_iter(_rng).take(samples).collect();

        NoiseModel{
            _sigma_k: sigma_k,
            _f_knee: f_knee,
            _slope: slope,
            _seed: seed,
            _samples: samples,
            noise_tod: white_noise_tod,
        }

    }
}

impl NoiseModel {
    pub fn get_noise_tod(&self) -> Vec<f32> {
        return self.noise_tod.clone();
    }

}