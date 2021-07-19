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
pub mod threadpool;
pub mod misc;
pub mod noisemodel;
pub mod plot_suite;
pub mod conjugategradient;
use threadpool::ThreadPool;
use std::{fs::File, io::Write, sync::mpsc, usize, vec};
use colored::Colorize;
use conjugategradient::conjgrad;
use iteratorscustom::FloatIterator;
use num::{ToPrimitive, complex::Complex32};
use rustfft::{FftPlanner, num_complex::Complex};
// use rustfft::algorithm::Radix4;
// use noisemodel::NoiseModel;
// use std::time::Instant;

// use gnuplot::*;

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
                //let noise = NoiseModel::new(50.0, 7e9, 1.0/20.0, 0.1, 1.0, 123, i.len());
                //let tod_noise = noise.get_noise_tod();
                let mut tmp: Vec<f32> = Vec::new();
                for (_n, (k, l)) in i.into_iter().zip(j.iter()).enumerate(){
                    let t_sky = sky[match l.to_usize() {Some(p) => p, None=>0}];
                    //let r = tod_noise[_n];
                    tmp.push(0.54*k + t_sky);
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

        
        let mut signal_maps: Vec<Vec<f32>> = Vec::new();
        let mut hit_maps: Vec<Vec<i32>>    = Vec::new();

        let pix = &self.pix;
        let tod = &self.tod;
        let num_threads = num_cpus::get();
        let bin_pool = ThreadPool::new(num_threads);
        let (tx, rx) = mpsc::channel();

        for i in 0..tod.len() {
            
            let t = tod[i].clone();
            let p = pix[i].clone();
            
            let tx = tx.clone();
            bin_pool.execute(move ||{
                  
                let (sig_par, hit_par) = bin_map(t.clone(), p.clone(), 128);
                
                tx.send((sig_par, hit_par)).unwrap();
            });
        }

        for _i in 0..tod.len() {
            let rec = rx.recv().unwrap();
            signal_maps.push(rec.0.clone());
            hit_maps.push(rec.1.clone());
        }

        let mut final_sig: Vec<f32> = vec![0.0; NUM_PIX];
        let mut final_hit: Vec<i32> = vec![0; NUM_PIX];

        
        for idx in 0..signal_maps.len() {
            let s = signal_maps[idx].clone();
            let h = hit_maps[idx].clone();

            for pidx in 0..NUM_PIX{
                let signal = s[pidx];
                let hits = h[pidx];
                final_sig[pidx] += signal;
                final_hit[pidx] += hits;
            }
        }
        
        println!("{}", "COMPLETED".bright_green());

        /***PRINT ON FILE */
        println!("");
        let id_number = self.get_mcid();
        let file_name = format!("binned_{}.dat", id_number);

        println!("Print maps on file: {}", file_name.bright_green().bold());

        let mut f = File::create(file_name).unwrap();

        let hit: Vec<String> = final_hit.iter().map(|a| a.to_string()).collect();
        let sig: Vec<String> = final_sig.iter().map(|a| a.to_string()).collect();

        for (i,j) in hit.iter().zip(sig.iter()) {
            writeln!(f, "{}\t{}",i, j).unwrap();
        }

        drop(final_hit);
        drop(final_sig);
        drop(signal_maps);
        drop(hit_maps);
        drop(bin_pool);
        drop(tod);
        drop(pix);
        println!("{}", "COMPLETED".bright_green());
    }
}

// GLS DENOISE
impl <'a> Obs <'a>{

    pub fn gls_denoise(&self, _tol: f32, _maxiter: usize, _nside: usize){
        
        println!("{}", "Execution of the gls_denoise".bright_blue().bold());

        const NUM_PIX: usize = 12*128*128;
        let _x: Vec<f32> = vec![0.0; NUM_PIX];

        let (tx, rx) = mpsc::channel();


        let tods = &self.tod;
        let pixs = &self.pix;

        //println!("Len TODs: {}", tods.len()); // 55 * n_hour

        let mut partial_maps: Vec<Vec<f32>> = Vec::new();
       
        let num_threads = num_cpus::get();
        let my_pool_b = ThreadPool::new(num_threads);
        for idx_th in 0..tods.len() {
            
            let t = tods[idx_th].clone();
            let p = pixs[idx_th].clone();
            let tx = tx.clone();
            my_pool_b.execute(move || {
                let b = get_b(t, p, 128);
              
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

        let map = conjgrad(a(), b, _tol, _maxiter, p(), pixs.clone());

        /***PRINT ON FILE */
        println!("");
        let id_number = self.get_mcid();
        let file_name = format!("gls_{}.dat", id_number);
        println!("Print maps on file: {}", file_name.bright_green().bold());
        let mut f = File::create(file_name).unwrap();
        let sig: Vec<String> = map.iter().map(|a| a.to_string()).collect();
        for i in sig.iter() {
            writeln!(f, "{}",i).unwrap();
        }
        println!("{}", "COMPLETED".bright_green());          
    
    }// End of GLS_DENOISE

} 


// UTILS functions
pub fn fn_noise_prior(f: f32, alpha: f32, f_k: f32, sigma: f32, _n: f32) -> f32 {
    let mut _np: f32 = 0.0;
    if f > 0.0 {
      
        let _np_g = f32::exp( -((10.0 - f) * (10.0-f))   /   (2.0 * 0.0002));
        _np = sigma * f32::powf(  1.0 + f_k/(10.0-f), alpha.clone()) + 8E8 * _np_g ;

    } else {

        let _np_g = f32::exp( -((10.0 + f) * (10.0-f))   /   (2.0 * 0.0002));
        _np = sigma*sigma * f32::powf(  1.0 + f_k/(10.0+f), alpha.clone()) + 8E8 * _np_g;
    }
    _np
} 

pub fn kaiser(beta: f32, length: i32) -> Vec<f32> {
    use crate::misc::bessel_i0 as bessel;

    let mut window: Vec<f32> = Vec::new();

    let start: f32 = (-(length - 1) / 2) as f32;
    let end: f32 = ((length - 1) / 2) as f32;

    let n_idx = iteratorscustom::FloatIterator::new(start, end, match length.to_u32(){Some(p) => p, None => 0});
     
    for n in n_idx {
        let m = length as f32;
        window.push(bessel(beta * (1. - (n / (m / 2.)).powi(2)).sqrt()) / bessel(beta))
    }// 70-80
 
    window
}

pub fn hann(length: i32) -> Vec<f32> {
    let mut window: Vec<f32> = Vec::new();
    let n_idx = iteratorscustom::FloatIterator::new(0.0, length as f32, match length.to_u32(){Some(p) => p, None => 0});
    
    for n in n_idx {
        window.push(f32::powi( f32::sin(3.141592 * n / (length as f32)), 2));
    }

    window
}

pub fn denoise(tod: Vec<f32>, _alpha: f32, _f_k: f32, _sigma: f32, _fs: f32) -> Vec<f32> {
    
    let win: Vec<f32> = kaiser(5.0, match tod.len().to_i32(){Some(p)=>p, None=> 0});//70!!!!!
    //let win: Vec<f32> = hann( match tod.len().to_i32(){Some(p)=>p, None=> 0});


    // let mut fg = Figure::new();
    // fg.axes2d().lines(0..win.len(), win.clone(), &[Caption("Kaiser")]).
    //             lines(0..win.len(), win_h, &[Caption("Hanning")]);
    // fg.show().unwrap();

    //let now = Instant::now();
    let mut input: Vec<Complex<f32>> = tod.iter().zip(win.iter()).map(|x| Complex32::new(
        *x.0 *x.1, 
        0.0)).
        collect(); // ~0
    //println!("Denoise process: {}", now.elapsed().as_millis());

    
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(tod.len()); // 10
    

    //let fft = Radix4::new()
    
    fft.process(&mut input); // 50

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
    // fg.axes2d().points(0..input.len()/2, input.iter().map(|t| t.re).collect::<Vec<f32>>(), &[Caption("FFT - raw")]).set_x_log(Some(10.0)).set_y_log(Some(10.0)).
    //             points(0..tod_denoised.len()/2, tod_denoised.iter().map(|f| f.re).collect::<Vec<f32>>(), &[Caption("FFT - denoised")]).set_x_log(Some(10.0)).set_y_log(Some(10.0)).
    //             lines(0..noise_p.len()/2, noise_p, &[Caption("Noise model")]);
    // fg.show().unwrap();
    
    let mut planner = FftPlanner::new(); // 60
    let ifft = planner.plan_fft_inverse(tod.len());

    ifft.process(&mut tod_denoised);
    
    let tod_denoised: Vec<Complex32> = tod_denoised.iter().map(|c|  {
        let (module, angle) = c.to_polar();
        let module2 = module / (tod.len() as f32);
        Complex32::from_polar(module2, angle)
    }).collect();
    
    let tod_real: Vec<f32> = tod_denoised.iter().zip(   win.iter()).map(|t| t.0.re  / t.1 ).collect();
    
    // let mut fg = Figure::new();
    // fg.axes2d().lines(0..tod_real.len(), tod_real.clone(), &[Caption("TOD denoised")]).lines(0..tod.len(), tod, &[Caption("TOD RAW")]);
    // fg.show().unwrap();

    tod_real
    
}

pub fn get_b(tod: Vec<f32>, pix: Vec<i32>, nside: usize) -> Vec<f32> {
    let mut b: Vec<f32> = vec![0.0; 12*128*128];
    let tod_n = denoise(tod.clone(), 4.0/3.0, 7.0, 30.0, 20.0);
    let (map, _) = bin_map(tod_n.clone(), pix, nside);
    for i in 0..12*nside*nside {
        b[i] += map[i];
    }
    b
}

fn a() -> Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>> {
    Box::new(|_x: Vec<f32>, pointings: Vec<Vec<i32>>|  {
        let mut temp_maps: Vec<Vec<f32>> = Vec::new();
        let mut res: Vec<f32> = vec![0.0; 12*128*128];

        let num_threads = num_cpus::get();
        let pool_denoise = ThreadPool::new(num_threads);
        let (tx, rx) = mpsc::channel();

        for i_det in pointings.iter() {
            let mut tmp: Vec<f32> = Vec::new();
            let tx = tx.clone();
            let point_i = i_det.clone();
            let x = _x.clone();
            pool_denoise.execute(move ||{
                for i in point_i.iter() {
                    tmp.push(x[*i as usize]);
                }
                let tmp_denoised = denoise(tmp.clone(), 4.0/3.0, 7.0, 30.0, 20.0);
                let (map, _) = bin_map(tmp_denoised.clone(), point_i.clone(), 128);
                // let mut final_map = Vec::new();
                // for (m, h) in map.iter().zip(hit.iter()) {
                //     final_map.push(m.clone()/(h.clone() as f32));
                // }
                tx.send(map).unwrap(); 
            });
        }

        for _i in 0..pointings.len() {
            temp_maps.push(rx.recv().unwrap());
        }
        
        for map in temp_maps.iter() {
            res = res.iter().zip(map.iter()).map(|(r, m)| r+m).collect::<Vec<f32>>();
        }
                                   
        res

    })
}

pub fn p() -> Box<dyn Fn(Vec<f32>) -> Vec<f32>> {

    Box::new(|m| m.iter().map(|m| { 1.0*m } ).collect::<Vec<f32>>())
}

pub fn bin_map(tod: Vec<f32>, pix: Vec<i32>, nside: usize) -> (Vec<f32>, Vec<i32>) {

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