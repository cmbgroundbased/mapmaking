/*!
# Strip MapMaking library
 Lorem ipsum dolor sit amet, consectetur adipiscing elit. Curabitur facilisis consectetur arcu. Etiam semper, sem sit amet lacinia dignissim, mauris eros rutrum massa, a imperdiet orci urna vel elit. Nulla at sagittis lacus. Curabitur eu gravida turpis. Mauris blandit porta orci. Aliquam fringilla felis a sem aliquet rhoncus. Suspendisse porta, mi vel euismod porta, mi ex cursus diam, quis iaculis sapien massa eget massa. Fusce sit amet neque vel turpis interdum tempus et non nisl. Nunc aliquam nunc vitae justo accumsan pretium. Morbi eget urna quis ex pellentesque molestie. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Integer vehicula vehicula tortor sit amet dignissim. Duis finibus, felis ut fringilla tincidunt, mi lectus fermentum eros, ut laoreet justo lacus id urna.

Duis iaculis faucibus mollis. Maecenas dignissim efficitur ex. Sed pulvinar justo a arcu lobortis imperdiet. Suspendisse placerat venenatis volutpat. Aenean eu nulla vitae libero porta dignissim ut sit amet ante. Vestibulum porttitor sodales nibh, nec imperdiet tortor accumsan quis. Ut sagittis arcu eu efficitur varius. Etiam at ex condimentum, volutpat ipsum sed, posuere nibh. Sed posuere fringilla mi in commodo. Ut sodales, elit volutpat finibus dapibus, dui lacus porttitor enim, ac placerat erat ligula quis ipsum. Morbi sagittis et nisl mollis fringilla. Praesent commodo faucibus erat, nec congue lectus finibus vitae. Sed eu ipsum in lorem congue vehicula. 

# Using the Strip MapMaking

Duis iaculis faucibus mollis. Maecenas dignissim efficitur ex. Sed pulvinar justo a arcu lobortis imperdiet. Suspendisse placerat venenatis volutpat. Aenean eu nulla vitae libero porta dignissim ut sit amet ante. Vestibulum porttitor sodales nibh, nec imperdiet tortor accumsan quis. Ut sagittis arcu eu efficitur varius. Etiam at ex condimentum, volutpat ipsum sed, posuere nibh. Sed posuere fringilla mi in commodo. Ut sodales, elit volutpat finibus dapibus, dui lacus porttitor enim, ac placerat erat ligula quis ipsum. Morbi sagittis et nisl mollis fringilla. Praesent commodo faucibus erat, nec congue lectus finibus vitae. Sed eu ipsum in lorem congue vehicula. 
*/

pub mod directory;
pub mod iteratorscustom;
pub mod sky;
pub mod noisemodel;
pub mod conjugategradient;

use std::fs::File;
use std::io::Write;
use colored::Colorize;
use num::ToPrimitive;
use num::complex::Complex32;
//use rand_distr::{Distribution, Normal};
use sky::Sky;
use noisemodel::NoiseModel;
use fftw::array::AlignedVec;
use fftw::plan::*;
use fftw::types::*;
use gnuplot::*;

use crate::conjugategradient::conjgrad;



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
        tod: & mut Vec<Vec<f32>>) -> Self 
        {
            let my_sky = Sky::new();
            let t_map = my_sky.get_t_sky();


            for (i, j) in tod.into_iter().zip(pix.iter()){
                let noise = NoiseModel::new(50.0, 7e9, 1.0/20.0, 0.1, 1.0, 123, i.len());
                let tod_noise = noise.get_noise_tod();
                for (n, (k, l)) in i.into_iter().zip(j.iter()).enumerate(){
                    let t_sky = t_map[match l.to_usize() {Some(p) => p, None=>0}];
                    let r = tod_noise[n];
                    *k = 0.56*(*k) + r + t_sky;
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
// Starting from the dummy binning, to the
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

        // This could be ok...
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
        let file_name = format!("mappe_{}.dat", id_number);

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

impl Obs {

    pub fn gls_denoise(&self, tol: f32, maxiter: usize, nside: usize){
              
        const NUM_PIX: usize = 12*128*128;

        let tods = self.get_tod();
        let pixs = self.get_pix();
        let mut b: Vec<f32> = Vec::new();
        for i in 0..NUM_PIX {
            b.push(0.0);
        }

        for (i, j) in tods.iter().zip(pixs.iter()) {
            let tod = denoise(i.clone(), 5.0/11.0, 0.1, 0.003);
            let (map, hit) = bin_map(tod, j.clone(), nside); // ERROREEEEE!!!!
            for i in 0..NUM_PIX {
                b[i] += map[i];
            }
            
            
        }

        // Check the iteration on the detectors. Is it correct here or do I have to do it 
        // into the CG function? In my opinion the CG fun has to manage only one TOD

        let mut a = |x: Vec<f32>| {
            
            let mut res: Vec<f32> = Vec::new();
            for _i in 0..NUM_PIX {
                res.push(0.0);
            }
            let mut tod: Vec<Vec<f32>> = Vec::new();

            for i_det in pixs.iter() {
                let mut tmp: Vec<f32> = Vec::new();
                for pix in i_det {
                    let idx = match pix.to_usize() {Some(p) => p, None=> 0,};
                    tmp.push(x[idx]);
                }
                let mut _tmp = denoise(tmp, 5.0/11.0, 0.1, 0.003);
                let mut _tmp_b = _tmp.clone();
                tod.push(_tmp);
                let (map, _) = bin_map(_tmp_b, i_det.clone(), nside);
                for i in 0..NUM_PIX {
                    res[i] += map[i];
                }
            }
            return res;
            
        };

        let mut precon = |x: Vec<f32>, y: Vec<f32> | {
            let mut out: Vec<f32> = Vec::new();
            for i in 0..x.len() {
                out.push(y[i]);
            }
            out 
         };

        let mut map = conjgrad(a, b, 1E-10, 100, precon);

    }
}

pub fn bin_map(tod: Vec<f32>, pix: Vec<i32>, nside: usize) -> (Vec<f32>, Vec<i32>) {

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


pub fn denoise(tod: Vec<f32>, alpha: f32, f_k: f32, sigma: f32) -> Vec<f32> {
    
    // Do the complex to complex C2C

    let mut c2r: C2RPlan64 = C2RPlan::aligned(&[tod.len()], Flag::MEASURE).unwrap();
    let mut a = AlignedVec::new(tod.len()/2 + 1);
    let mut b: AlignedVec<f64> = AlignedVec::new(tod.len());

    for i in 0..(tod.len() / 2 + 1) {
        a[i] = c64::new(match tod[i].to_f64(){Some(p)=> p, None => 0.0}, 0.0);
    }
    c2r.c2r(&mut a, &mut b).unwrap();

    let mut index: f64 = 0.5;
    let mut y: Vec<f64> = Vec::new();
    let mut x: Vec<f64> = Vec::new();

    for i in 1..(tod.len()/2+1) {
        x.push(index.log10());  
        let b = b.to_vec();
        y.push(b[i].log10());
        index += 1.0;
    }
    
    // Filter it out by the nois prior

    let tod_corrected: Vec<f32> = Vec::new();

    tod_corrected
    // let mut fg = gnuplot::Figure::new();


    // fg.axes2d().lines_points(&x, &y, &[Caption("A line"), Color("black")]);
    
    // fg.show().unwrap();

}