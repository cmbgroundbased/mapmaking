pub mod iteratorscustom;

use rustfft::{FftPlanner};
use rustfft::num_complex::Complex;
use rustfft::num_traits::{ToPrimitive};
use iteratorscustom::FloatIterator;
use ndarray_npy::NpzReader;
use ndarray::Array1;
use std::fs::File;
use gnuplot::*;

fn main() {

    /* ONLY ATMOSPHERE TOD */
    let mut a = NpzReader::new(File::open("/home/algebrato/Progetti/mapmaking/one_tod/2022010721/43000008.0/043/I0.npz").unwrap()).unwrap();
    let c: Array1<f32> = a.by_name("arr_0.npy").unwrap();
    let tod = c.to_vec();

    let mut input: Vec<Complex<f32>> = Vec::new();
    
    for i in tod.iter() {
        input.push(Complex::<f32>::new(*i, 0.0));
    }

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft(tod.len(), rustfft::FftDirection::Forward);
    fft.process(&mut input);

    /* The input vector has to be ordered!! Something like fftshift */

    let mut new_input: Vec<Complex<f32>> = Vec::new(); 
    for _i in 0..input.len() {
        new_input.push(Complex::new(0.0, 0.0));
    }

    for i in 0..input.len()/2 {
        let idx_n = input.len()/2 - i;

        new_input[i] = input[idx_n];
    }

    for i in 1..input.len()/2 {
        let idx = i + input.len()/2;
        let idx_n = input.len() - i;
        new_input[idx] = input[idx_n];
    }

    let input = new_input;

    let mut noise_prior: Vec<f32> = Vec::new();
    let iterator_f = FloatIterator::new(-10.0, 10.0, match input.len().to_u32() {Some(p) => p, None => 0});

    let mut freq: Vec<f32> = Vec::new();
    for f in iterator_f{
        noise_prior.push(noise_prior_f(f, 5.0/3.0, 2.0, 30.0));
        freq.push(f);

    }

    let mut fg = Figure::new();
    // let avg: f32 =input.iter().map(|f| f.norm()).sum::<f32>() / (match input.len().to_f32() {Some(p) => p, None=> 0.0});

    fg.axes2d().
        points(freq.clone(), input.iter().map(|f| f.norm()).collect::<Vec<f32>>(), &[Caption("TOD FFT"), Color("red"), PointSymbol('O'), PointSize(0.5)]).
        lines(freq.clone(), noise_prior, &[Caption("Noise Prior"), Color("black")]).
        set_y_log(Some(10.0));//.set_x_range(Fix(0.01), Fix(10.0)).set_x_log(Some(10.0));
    fg.show().unwrap();

}

pub fn noise_prior_f(f: f32, alpha: f32, f_k: f32, sigma: f32) -> f32 {
    let mut _np: f32 = 0.0;
    if f > 0.0 {

        _np = sigma*sigma * f32::powf(  1.0 + f_k/(10.0-f), alpha.clone());
    } else {
        _np = sigma*sigma * f32::powf(  1.0 + f_k/(10.0+f), alpha.clone());
    }

    _np
} 