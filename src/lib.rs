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
pub mod threadpool2;
pub mod misc;
pub mod noisemodel;
// pub mod conjugategradient;
pub mod conjugategradient2;
pub mod plot_suite;
use threadpool2::ThreadPool;
use std::fs::File;
// use std::f64::consts::PI;
// use rustfft::num_traits::Zero;
use std::io::Write;
use std::sync::mpsc;
use std::sync::{Arc, Barrier};
use std::vec;
use std::thread::sleep;
// use std::ops::Index;
use colored::Colorize;




use iteratorscustom::FloatIterator;
// use ndarray::iter::Windows;
use num::ToPrimitive;

use num::complex::Complex32;
// use num::traits::sign;
//use rand_distr::{Distribution, Normal};
// use fftw::array::AlignedVec;
// use fftw::plan::*;
// use fftw::types::*;

use rustfft::FftPlanner;
use rustfft::num_complex::Complex;
use noisemodel::NoiseModel;
// use rustfft::num_traits::Zero;

//use plot_suite::plot_vector;
//use gnuplot::{AxesCommon, Caption, Color, Figure};

use crate::conjugategradient2::conjgrad2;
//use crate::conjugategradient::conjgrad;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct Obs <'a> {
    start: String,
    stop: String,
    detector: Vec<String>,
    mc_id: u8,
    alpha: f32,
    f_knee: f32,
    pix: Vec<Vec<i32>>,
    tod: Vec<Vec<f32>>,
    sky_t: Vec<f32>,
    phantom: PhantomData<&'a f32>,
}

// ```
// Documentation
// Creation function
// ```
impl <'a> Obs <'a> {
    pub fn new(
        start: String, 
        stop: String, 
        detector: Vec<String>, 
        mc_id: u8, 
        alpha: f32, 
        f_knee: f32, 
        tod: Vec<Vec<f32>>, 
        sky: Vec<f32>,
        pix: Vec<Vec<i32>> ) -> Self 
        {
            let mut tod_final: Vec<Vec<f32>> = Vec::new();

            for (i, j) in tod.iter().zip(pix.iter()){
                let noise = NoiseModel::new(50.0, 7e9, 1.0/20.0, 0.1, 1.0, 123, i.len());
                let tod_noise = noise.get_noise_tod();
                let mut tmp: Vec<f32> = Vec::new();
                for (_n, (k, l)) in i.into_iter().zip(j.iter()).enumerate(){
                    let t_sky = sky[match l.to_usize() {Some(p) => p, None=>0}];
                    let r = tod_noise[_n];
                    tmp.push(k + t_sky + r)
                }
                tod_final.push(tmp);
            }
            
            return Obs {
                start,
                stop,
                detector,
                mc_id,
                alpha,
                f_knee,
                pix,
                tod: tod_final,
                sky_t: sky,
                phantom: PhantomData,
            }
    }
}

// The `get` methods
impl <'a> Obs <'a>{
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
impl <'a> Obs <'a>{
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

pub fn fn_noise_prior(f: f32, alpha: f32, f_k: f32, sigma: f32, _n: f32) -> f32 {
    let mut _np: f32 = 0.0;
    if f > 0.0 {
        let _np_g = f32::exp( -((10.0 - f) * (10.0-f))   /   (2.0 * 0.0002));
        _np = sigma * f32::powf(  1.0 + f_k/(10.0-f), alpha.clone()) + 8E6 * _np_g ;
    } else {
        //sigma * (1+f_knee/(-1.0*f[f<0]))**alpha + 8e8 *np.exp(-f[f<0]**2/(2*0.0002))  

        let _np_g = f32::exp( -((10.0 + f) * (10.0-f))   /   (2.0 * 0.0002));
        _np = sigma*sigma * f32::powf(  1.0 + f_k/(10.0+f), alpha.clone()) + 8E6 * _np_g;
    }

    _np
} 

// Design Kaiser window from parameters.
// The length depends on the parameters given, and it's always odd.
pub fn kaiser(beta: f32, length: i32) -> Vec<f32> {
    use crate::misc::bessel_i0 as bessel;

    let mut window: Vec<f32> = Vec::with_capacity(length as usize);

    let start: f32 = (-(length - 1) / 2) as f32;
    let end: f32 = ((length - 1) / 2) as f32;

    let n_idx = iteratorscustom::FloatIterator::new(start, end, match length.to_u32(){Some(p) => p, None => 0});

    for n in n_idx {
        let m = length as f32;
        window.push(bessel(beta * (1. - (n / (m / 2.)).powi(2)).sqrt()) / bessel(beta))
    }

    window
}


pub fn denoise(tod: Vec<f32>, _alpha: f32, _f_k: f32, _sigma: f32, _fs: f32) -> Vec<f32> {
    
    let kaiser_win: Vec<f32> = kaiser(5.0, match tod.len().to_i32(){Some(p)=>p, None=> 0});
    let mut input: Vec<Complex<f32>> = tod.iter().zip(kaiser_win.iter()).map(|x| Complex32::new(
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
        noise_p.push(fn_noise_prior(f, _alpha, _f_k, _sigma, match tod.len().to_f32(){Some(p)=>p, None=>0.0} ));
        freq.push(f);
    }
 
    // Denoise
    let mut tod_denoised: Vec<Complex32> = input.iter().zip(noise_p.iter()).map(|(a, b)| {
        let (_, angle) = a.to_polar();
        let module = a.norm()/b;
        Complex32::from_polar(module, angle)
    }).collect();

    // let mut fg = Figure::new();
    // fg.axes2d().points(0..input.len(), input.iter().map(|t| t.norm()).collect::<Vec<f32>>(), &[Caption("FFT - raw")]).set_x_log(Some(10.0)).set_y_log(Some(10.0)).
    //             points(0..tod_denoised.len(), tod_denoised.iter().map(|f| f.norm()).collect::<Vec<f32>>(), &[Caption("FFT - denoised")]).set_x_log(Some(10.0)).set_y_log(Some(10.0)).
    //             lines(0..noise_p.len(), noise_p, &[Caption("Noiser prior"), Color("red")]);
    // fg.show().unwrap();

    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(tod.len());

    ifft.process(&mut tod_denoised);
    let tod_denoised: Vec<Complex32> = tod_denoised.iter().map(|c|  {
        let (_, angle) = c.to_polar();
        let module = c.norm() / (tod.len() as f32);
        Complex32::from_polar(module, angle)
    }).collect();

    
    let tod_real: Vec<f32> = tod_denoised.iter().zip(kaiser_win.iter()).map(|t| t.0.re  / t.1 ).collect();

    // let mut fg = Figure::new();
    // fg.axes2d().lines(0..tod_real.len(), tod_real.clone(), &[Caption("TOD - denoised")]);
    // fg.show().unwrap();

    tod_real

}


pub fn get_b(tod: &Vec<f32>, pix: &Vec<i32>, nside: usize) -> Vec<f32> {
    let mut b: Vec<f32> = vec![0.0; 12*128*128];
    let tod_n = denoise(tod.clone(), 8.0/3.0, 7.0, 1.0, 20.0);
    let (map, _) = bin_map(tod_n.clone(), &pix, nside);
    for i in 0..12*nside*nside {
        b[i] += map[i];
    }
    // let mut a = 1.0;
    // loop     {

    //     a +=1.0;
    //     a -=1.0;
    //     a = f32::sqrt(a);

    //     if a == 0.0 {
    //         break;
    //     }
    // }

    b

}

pub fn a() -> Box<dyn Fn(&Vec<f32>, &Vec<i32>) -> Vec<f32>> {
    Box::new(|_x: &Vec<f32>, pointings: &Vec<i32>|  {
        
        let mut res: Vec<f32> = vec![0.0; 12*128*128];        
        let mut _tmp: Vec<f32> = vec![0.0; pointings.len()];                 
        let mut index: usize = 0;

        for pix in pointings.iter() {
            let pix_id = match pix.to_usize(){Some(p) => p, None => 0};
            _tmp[index] += _x[pix_id];
            index += 1;
            
        }
        
        let tmp_denoised = denoise(_tmp, 8.0/3.0, 7.0, 1.0, 20.0);
        let (map, _) = bin_map(tmp_denoised, pointings, 128);
        for i in 0..12*128*128{
            res[i] += map[i];
        }
                          
        res

    })
}


pub fn p() -> Box<dyn Fn(Vec<f32>) -> Vec<f32>> {
    Box::new(|x| x)
}



impl <'a> Obs <'a>{

    pub fn gls_denoise(&self, _tol: f32, _maxiter: usize, _nside: usize){
        
        println!("{}", "Execution of the gls_denoise".bright_blue().bold());

        const NUM_PIX: usize = 12*128*128;
        let _x: Vec<f32> = vec![0.0; NUM_PIX];

        let (tx, rx) = mpsc::channel();


        let tods = &self.tod;
        let pixs = &self.pix;

        println!("Len TODs: {}", tods.len()); // 55 * n_hour

        let mut partial_maps: Vec<Vec<f32>> = Vec::new();
       
        let my_pool_b = ThreadPool::new(16);
        for idx_th in 0..tods.len() {
            
            let t = tods[idx_th].clone();
            let p = pixs[idx_th].clone();
            let tx = tx.clone();
            my_pool_b.execute(move || {
                let b = get_b(&t, &p, 128);
              
                tx.send(b).expect("channel will be there waiting for the pool");
             
            });            
        }

        for _ in 0..tods.len(){
            partial_maps.push(rx.recv().unwrap());
        }

        let mut b: Vec<f32> = vec![0.0; 12*128*128];

        for i in partial_maps.iter(){
            for (n,j) in i.iter().enumerate(){
                b[n] += j;
            }
        }

        // let mut partial_maps: Vec<Vec<f32>> = Vec::new();
        // let (tx_a, rx_a) = mpsc::channel();
        // let my_pool_a = ThreadPool::new(16);
        // for idx_th in 0..tods.len() {
        //     let p = pixs[idx_th].clone();
        //     let b = b.clone();
        //     let tx_a = tx_a.clone();
        //     my_pool_a.execute(move || {
        //         let p = p.clone();
        //         //let b = b.clone();
        //         let map = conjgrad2(a(), b, 1e-4, 2, &p);
        //         tx_a.send(map).expect("channel will be there waiting for the pool");
        //     });


        // }
        // for _ in 0..tods.len(){
        //     partial_maps.push(rx_a.recv().unwrap());
        // }

        // let mut map: Vec<f32> = vec![0.0; 12*128*128];

        // for i in partial_maps.iter(){
        //     for (n,j) in i.iter().enumerate(){
        //         map[n] += j;
        //     }
        // }
        
        /***PRINT ON FILE */
        println!("");
        let id_number = self.get_mcid();
        let file_name = format!("gls_{}.dat", id_number);
        println!("Print maps on file: {}", file_name.bright_green().bold());
        let mut f = File::create(file_name).unwrap();
        let sig: Vec<String> = b.iter().map(|a| a.to_string()).collect();
        for i in sig.iter() {
            writeln!(f, "{}",i).unwrap();
        }
        println!("{}", "COMPLETED".bright_green());          
    }


} // End of GLS_DENOISE

pub fn bin_map(tod: Vec<f32>, pix: &Vec<i32>, nside: usize) -> (Vec<f32>, Vec<i32>) {

    let num_pixs: usize = 12*nside*nside;

    let mut signal_map: Vec<f32> = vec![0.0; num_pixs];
    let mut hit_map: Vec<i32> = vec![0; num_pixs];

       
    let mut iterator: usize = 0;
    for i in pix.iter() {

        let pixel = match i.to_usize(){Some(p)=> p, None=> 0};
        hit_map[pixel] += 1;
        signal_map[pixel] += tod[iterator];
        iterator += 1;
    }

    (signal_map, hit_map)

}