/*!
# Strip MapMaking library
 Lorem ipsum dolor sit amet, consectetur adipiscing elit. Curabitur facilisis consectetur arcu. Etiam semper, sem sit amet lacinia dignissim, mauris eros rutrum massa, a imperdiet orci urna vel elit. Nulla at sagittis lacus. Curabitur eu gravida turpis. Mauris blandit porta orci. Aliquam fringilla felis a sem aliquet rhoncus. Suspendisse porta, mi vel euismod porta, mi ex cursus diam, quis iaculis sapien massa eget massa. Fusce sit amet neque vel turpis interdum tempus et non nisl. Nunc aliquam nunc vitae justo accumsan pretium. Morbi eget urna quis ex pellentesque molestie. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Integer vehicula vehicula tortor sit amet dignissim. Duis finibus, felis ut fringilla tincidunt, mi lectus fermentum eros, ut laoreet justo lacus id urna.

Duis iaculis faucibus mollis. Maecenas dignissim efficitur ex. Sed pulvinar justo a arcu lobortis imperdiet. Suspendisse placerat venenatis volutpat. Aenean eu nulla vitae libero porta dignissim ut sit amet ante. Vestibulum porttitor sodales nibh, nec imperdiet tortor accumsan quis. Ut sagittis arcu eu efficitur varius. Etiam at ex condimentum, volutpat ipsum sed, posuere nibh. Sed posuere fringilla mi in commodo. Ut sodales, elit volutpat finibus dapibus, dui lacus porttitor enim, ac placerat erat ligula quis ipsum. Morbi sagittis et nisl mollis fringilla. Praesent commodo faucibus erat, nec congue lectus finibus vitae. Sed eu ipsum in lorem congue vehicula. 

# Using the Strip MapMaking

Duis iaculis faucibus mollis. Maecenas dignissim efficitur ex. Sed pulvinar justo a arcu lobortis imperdiet. Suspendisse placerat venenatis volutpat. Aenean eu nulla vitae libero porta dignissim ut sit amet ante. Vestibulum porttitor sodales nibh, nec imperdiet tortor accumsan quis. Ut sagittis arcu eu efficitur varius. Etiam at ex condimentum, volutpat ipsum sed, posuere nibh. Sed posuere fringilla mi in commodo. Ut sodales, elit volutpat finibus dapibus, dui lacus porttitor enim, ac placerat erat ligula quis ipsum. Morbi sagittis et nisl mollis fringilla. Praesent commodo faucibus erat, nec congue lectus finibus vitae. Sed eu ipsum in lorem congue vehicula. 
*/
extern crate rustfft;
pub mod directory;
pub mod iteratorscustom;
pub mod sky;
pub mod noisemodel;
pub mod conjugategradient;
pub mod conjugategradient2;
pub mod plot_suite;

// use std::f64::consts::PI;
// use rustfft::num_traits::Zero;
use std::fs::File;
use std::io::Write;
// use std::ops::Index;
use colored::Colorize;

use iteratorscustom::FloatIterator;
// use ndarray::iter::Windows;
use num::ToPrimitive;

use num::complex::Complex32;
// use num::traits::sign;
//use rand_distr::{Distribution, Normal};
use sky::Sky;
use noisemodel::NoiseModel;
// use fftw::array::AlignedVec;
// use fftw::plan::*;
// use fftw::types::*;

use rustfft::FftPlanner;
use rustfft::num_complex::Complex;
// use rustfft::num_traits::Zero;

//use plot_suite::plot_vector;
//use gnuplot::*;

use crate::conjugategradient2::conjgrad2;
//use crate::conjugategradient::conjgrad;


#[derive(Debug)]
pub struct Obs {
    start: String,
    stop: String,
    detector: Vec<String>,
    mc_id: u8,
    alpha: f32,
    f_knee: f32,
    pix: Vec<Vec<i32>>,
    tod: Vec<Vec<f32>>,
    sky_t: Vec<f32>,
}

// ```
// Documentation
// Creation function
// ```
impl Obs {
    pub fn new(
        start: String, 
        stop: String, 
        detector: Vec<String>, 
        mc_id: u8, 
        alpha: f32, 
        f_knee: f32, 
        pix: Vec<Vec<i32>>, 
        tod: & mut Vec<Vec<f32>>,
        sky: Vec<f32>) -> Self 
        {

            for (i, j) in tod.into_iter().zip(pix.iter()){
                //let noise = NoiseModel::new(50.0, 7e9, 1.0/20.0, 0.1, 1.0, 123, i.len());
                //let tod_noise = noise.get_noise_tod();
                for (n, (k, l)) in i.into_iter().zip(j.iter()).enumerate(){
                    let t_sky = sky[match l.to_usize() {Some(p) => p, None=>0}];
                    //let r = tod_noise[n];
                    // *k = 0.56 * (*k) + r + t_sky; 
                    *k = 0.56 * (*k) + t_sky; 
                }
            }
            
            return Obs {
                start,
                stop,
                detector,
                mc_id,
                alpha,
                f_knee,
                pix,
                tod: tod.clone(),
                sky_t: sky,
            }
    }
}

// The `get` methods
impl Obs {
    pub fn get_start(&self) -> &String {
        &self.start
    }

    pub fn get_stop(&self) -> &String {
        &self.stop
    }

    pub fn get_detector(&self) -> &Vec<String> {
        &self.detector
    }

    pub fn get_mcid(&self) -> &u8 {
        &self.mc_id
    }

    pub fn get_pix(&self) -> &Vec<Vec<i32>> {
        &self.pix
    }
    pub fn get_tod(&self) -> &Vec<Vec<f32>> {
        &self.tod
    }
}

// Mitigation of the systematic effects
// Starting from the binning, to the
// implementation of a de_noise model
impl Obs {
    pub fn binning(&self) {

        println!("");
        println!("Start {}", "binning".bright_blue().bold());

        const NSIDE: usize = 128;
        const NUM_PIX: usize = NSIDE*NSIDE*12;


        // MEMORY!!!!!
        let mut signal_map: Vec<f32> = Vec::new();
        let mut hit_map: Vec<f64>   = Vec::new();
        for _i in 0..NUM_PIX {
            signal_map.push(0.0);
            hit_map.push(0.0)
        }

        let vect_pix = &self.pix;
        let tod = &self.tod;
        

        for (i, j) in vect_pix.iter().zip(tod.iter()) {            
            let mut iterator: usize = 0;
            for pix in i.iter() {

                let pixel = match pix.to_usize(){Some(p) => p, None=>0};
                hit_map[pixel] += 1.0;
                signal_map[pixel] += match j.get(iterator) {Some(s) => s.clone(), None=>0.0};
                iterator += 1;
            }
        }

        let vec_signal = signal_map;
        let vec_hit    = hit_map;

        println!("{}", "COMPLETED".bright_green());

        /***PRINT ON FILE */
        println!("");
        let id_number = self.get_mcid();
        let file_name = format!("binned_{}.dat", id_number);

        println!("Print maps on file: {}", file_name.bright_green().bold());

        let mut f = File::create(file_name).unwrap();

        let hit: Vec<String> = vec_hit.iter().map(|a| a.to_string()).collect();
        let sig: Vec<String> = vec_signal.iter().map(|a| a.to_string()).collect();

        for (i,j) in hit.iter().zip(sig.iter()) {
            writeln!(f, "{}\t{}",i, j).unwrap();
        }

        drop(vec_signal);
        drop(vec_hit);
        println!("{}", "COMPLETED".bright_green());
    }
}


impl Obs {

    pub fn dummy_denoise(&self) {
        let tod  = &self.tod;
        let pix  = &self.pix;
        
        let baseline_len: usize = 2 * 20; 
        let mut signal_map: Vec<f32> = Vec::new();
        let mut hit_map: Vec<f64>   = Vec::new();

        const NSIDE: usize = 128;
        const NUM_PIX: usize = NSIDE*NSIDE*12;

        for _i in 0..NUM_PIX {
            signal_map.push(0.0);
            hit_map.push(0.0)
        }
        
        for (tod_t, pix_t) in tod.iter().zip(pix.iter()) {
            let idx_max =  tod_t.len() / baseline_len;

            let mut tod_tmp: Vec<f32> = Vec::new();
            for idx in 0..(idx_max) {
                let start = baseline_len * idx;
                let stop = baseline_len * (idx+1);
                let mut sum: f32 = 0.0;

                for t in start..stop {
                    sum += tod_t[t];
                }
                for tod_idx in start..stop {
                    tod_tmp.push(tod_t[tod_idx] - sum/200.0);
                }
            }
            
            let mut iterator: usize = 0;
            for pix in pix_t.iter() {
                let pixel = match pix.to_usize(){Some(p) => p, None=>0};
                hit_map[pixel] += 1.0;
                signal_map[pixel] += match tod_tmp.get(iterator) {Some(s) => s.clone(), None=>0.0};
                iterator += 1;
            }
            
               
            
        }

        let vec_signal = signal_map;
        let vec_hit    = hit_map;
        println!("{}", "AVG REDUCTION COMPLETED".bright_green());

        /***PRINT ON FILE */
        println!("");
        let id_number = self.get_mcid();
        let file_name = format!("mappe_{}.dat", id_number);

        println!("Print maps on file: {}", file_name.bright_green().bold());

        let mut f = File::create(file_name).unwrap();

        let hit: Vec<String> = vec_hit.iter().map(|a| a.to_string()).collect();
        let sig: Vec<String> = vec_signal.iter().map(|a| a.to_string()).collect();

        for (i,j) in hit.iter().zip(sig.iter()) {
            writeln!(f, "{}\t{}",i, j).unwrap();
        }

        println!("{}", "WRITE MAP COMPLETED".bright_green());

    }
}

pub fn fn_noise_prior(f: f32, alpha: f32, f_k: f32, sigma: f32) -> f32 {
    let mut _np: f32 = 0.0;
    if f > 0.0 {

        _np = sigma*sigma * f32::powf(  1.0 + f_k/(10.0-f), alpha.clone());
    } else {
        _np = sigma*sigma * f32::powf(  1.0 + f_k/(10.0+f), alpha.clone());
    }

    _np
} 

/// Design Kaiser window from parameters.
///
/// The length depends on the parameters given, and it's always odd.
// fn kaiser(atten: f32, freq: Vec<f32>) -> Vec<f32> {
//     use crate::misc::bessel_i0 as bessel;

//     let beta: f32;
//     if atten > 50. {
//         beta = 0.1102 * (atten - 8.7);
//     } else if atten < 21. {
//         beta = 0.;
//     } else {
//         beta = 0.5842 * (atten - 21.).powf(0.4) + 0.07886 * (atten - 21.);
//     }

//     // Filter length, we want an odd length
//     let mut length: i32 = ((atten - 8.) / (2.285 * freq.ceil() as i32 + 1;
//     if length % 2 == 0 {
//         length += 1;
//     }

//     let mut window: Signal = Vec::with_capacity(length as usize);

//     for n in -(length - 1) / 2..=(length - 1) / 2 {
//         let n = n as f32;
//         let m = length as f32;
//         window.push(bessel(beta * (1. - (n / (m / 2.)).powi(2)).sqrt()) / bessel(beta))
//     }

//     debug!(
//         "Kaiser window design finished, beta: {}, length: {}",
//         beta, length
//     );

//     window
// }




pub fn denoise(tod: Vec<f32>, _alpha: f32, _f_k: f32, _sigma: f32, _fs: f32) -> Vec<f32> {
    
    let k0 = 3.141592 / ( tod.len() as f32);

    let sin_window: Vec<f32> = (0..tod.len()).map(|i| f32::sin(k0 * (i as f32))).collect();
    let cos_window: Vec<f32> = (0..tod.len()).map(|i| f32::cos(k0 * (i as f32))).collect();

    let mut input: Vec<Complex<f32>> = tod.iter().zip(sin_window.iter()).map(|x| Complex32::new(
        *x.0 *x.1, 
        0.0)).
        collect();
   
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(tod.len());
    fft.process(&mut input);

    let mut noise_p: Vec<f32> = Vec::new();
    let freq_iter: FloatIterator = FloatIterator::new(-10.0, 10.0, match input.len().to_u32() {Some(p) => p, None => 0});
    let mut freq: Vec<f32> = Vec::new();

    for f in freq_iter {
        noise_p.push(fn_noise_prior(f, _alpha, _f_k, _sigma));
        freq.push(f);
    }
 
    // Denoise
    let mut tod_corrected: Vec<Complex32> = input.iter().zip(noise_p.iter()).zip(cos_window.iter()).map(|((a, b), c)| Complex32::new( (a.re / b) * c, 0.0)).collect();

    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(tod.len());

    ifft.process(&mut tod_corrected);

    let tod_real: Vec<f32> = tod_corrected.iter().map(|f| f.re).collect();

    tod_real
}

impl Obs {

    pub fn gls_denoise(&self, tol: f32, maxiter: usize, nside: usize){
        
        println!("{}", "gls_denoise process in execution...".blue());

        const NUM_PIX: usize = 12*128*128;

        let tods = self.get_tod();
        let pixs = self.get_pix();
        let mut b: Vec<f32> = Vec::new();
        for _i in 0..NUM_PIX {
            b.push(0.0);
        }

        let _x = b.clone();

        for (i, j) in tods.iter().zip(pixs.iter()) {

            let tod = denoise(i.clone(), 8.0/3.0, 7.0, 1.0, 20.0);
            let (map, _) = bin_map(tod.clone(), &j.clone(), nside);
            for i in 0..NUM_PIX {
                b[i] += map[i];
            }
        }

        fn a() -> Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>> {
            Box::new(|_x: Vec<f32>, puntamenti:Vec<Vec<i32>>|  {
                let mut res: Vec<f32> = Vec::new();
                for _i in 0..NUM_PIX {
                    res.push(0.0);
                }
                for i_det in puntamenti.iter(){
                    let mut _tmp: Vec<f32> = Vec::new();
                    
                    for pix in i_det.iter() {
                        let pix_id = match pix.to_usize(){Some(p) => p, None => 0};
                        _tmp.push(_x[pix_id]);
                        
                    }

                    let tmp_denoise = denoise(_tmp, 8.0/3.0, 7.0, 1.0, 20.0);
                    let (map, _) = bin_map(tmp_denoise, i_det, 128);

                    for i in 0..NUM_PIX{
                        res[i] += map[i];
                    }

                }
                res 
            })
        }

        fn p() -> Box<dyn Fn(Vec<f32>) -> Vec<f32>> {
            Box::new(|x| x)
        }


        let _map = conjgrad2(a(), b, tol, maxiter, p(), pixs);

        
        /*PRINT ON FILE
        ************************************************************************/
        println!("");
        let id_number = self.get_mcid();

        let file_name = format!("gls_denoise_{}.dat", id_number);

        println!("Print maps on file: {}", file_name.bright_green().bold());

        let mut f = File::create(file_name).unwrap();

        let sig: Vec<String> = _map.iter().map(|a| a.to_string()).collect();

        for i in sig.iter() {
            writeln!(f, "{}", i).unwrap();
        }
        println!("{}", "WRITE MAP COMPLETED".bright_green());
        /*
        ************************************************************************/
    }
} // End of GLS_DENOISE

pub fn bin_map(tod: Vec<f32>, pix: &Vec<i32>, nside: usize) -> (Vec<f32>, Vec<i32>) {

    let num_pixs: usize = 12*nside*nside;

    let mut signal_map: Vec<f32> = Vec::new();
    let mut hit_map: Vec<i32> = Vec::new();
    for _i in 0..num_pixs {
        signal_map.push(0.0);
        hit_map.push(0);
    }

    let mut index: usize = 0;
    for i in 0..tod.len() {
        let pix_id = match pix[i].to_usize() {Some(p) => p, None => 0};
        signal_map[pix_id] = tod[index];
        hit_map[pix_id] += 1;
        index += 1
               
    }
    
    (signal_map, hit_map)

}


impl Obs {
    pub fn atm_mitigation(&self, baselines_length: usize, maxiter: usize, tol: f32, nside: usize){ // It does not work...
        
        println!("{}", "atm_mitigation process in execution...".blue());

        const NUM_PIX: usize = 12*128*128;

        let tods = self.get_tod();
        let pixs = self.get_pix();
        
        let mut b: Vec<f32> = Vec::new();
        let mut x: Vec<f32> = Vec::new();

        let b_l = baselines_length*20; // In seconds
        let b_l_f32 = match (baselines_length*20).to_f32(){Some(p) => p, None => 0.0}; 

        for _i in 0..NUM_PIX {
            b.push(0.0);
            x.push(0.0);
        }

        for (i, j) in tods.iter().zip(pixs.iter()) {
            let mut new_tod: Vec<f32> = Vec::new();
            

            for c in i.chunks(b_l) {
                let mut avg: f32 = c.iter().map(|i| i).sum();
                avg /= b_l_f32;
                let mut new_c: Vec<f32> = c.iter().map(|i| i-avg).collect();
                new_tod.append(&mut new_c);
            } 
            let (map, _hit) = bin_map(new_tod.clone(), &j.clone(), nside);
            //let map: Vec<f32> = map.iter().zip(hit.iter()).map(|p| *p.0 / match (*p.1).to_f32(){Some(p) => p, None=>0.0}).collect();
            for i in 0..NUM_PIX {
                b[i] += map[i];
            }
        }
        
        fn a() -> Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>> {
            Box::new(|_x: Vec<f32>, puntamenti:Vec<Vec<i32>>|  { 
                let mut res: Vec<f32> = Vec::new();
                for _i in 0..NUM_PIX {
                    res.push(0.0);
                }
                for i_det in puntamenti.iter(){
                    let mut _tmp: Vec<f32> = Vec::new();
                    for pix in i_det.iter() {
                        let pix_id = match pix.to_usize(){Some(p) => p, None => 0};
                        _tmp.push(_x[pix_id]);                      
                    }
                    let mut tmp_denoised: Vec<f32> = Vec::new();
                    
                    /*HARDCODED!!!!!!**************** */
                    let baselines_length = 1*20;
                    /******************************** */
                    
                    let b_l = baselines_length*20;
                    let b_l_f32 = match (baselines_length*20).to_f32(){Some(p) => p, None => 0.0};

                    for c in _tmp.chunks(b_l) {
                        let mut avg: f32 = c.iter().map(|i| i).sum();
                        avg /= b_l_f32;
                        let mut new_c: Vec<f32> = c.iter().map(|i| i-avg).collect();
                        tmp_denoised.append(&mut new_c);
                    } 

                    let (map, _) = bin_map(tmp_denoised, i_det, 128);

                    for i in 0..NUM_PIX{
                        res[i] += map[i];
                    }
                }
               _x 
            })
        }

        fn p() -> Box<dyn Fn(Vec<f32>) -> Vec<f32>> {
            Box::new(|x| x)
        }

        let _map = conjgrad2(a(), b, tol, maxiter, p(), pixs);
        
        /*PRINT ON FILE
        ************************************************************************/
        println!("");
        let id_number = self.get_mcid();
        let file_name = format!("atm_destripe_{}.dat", id_number);

        println!("Print maps on file: {}", file_name.bright_green().bold());

        let mut f = File::create(file_name).unwrap();

        let sig: Vec<String> = _map.iter().map(|a| a.to_string()).collect();

        for i in sig.iter() {
            writeln!(f, "{}", i).unwrap();
        }

        println!("{}", "WRITE MAP COMPLETED".bright_green());
        /*
        ************************************************************************/


        
    }
}
