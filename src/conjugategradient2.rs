// extern crate rayon;
// use rayon::prelude::*;

pub fn conjgrad2(a: Box<dyn Fn(&Vec<f32>, &Vec<i32>) -> Vec<f32>>, b: Vec<f32>, _tol: f32, maxiter:usize, pixes: &Vec<i32>) -> Vec<f32> {
    // println!("{}", "Starting CG...");
    let n = b.len();
 
    let mut x: Vec<f32> = vec![0.0; b.len()];
    let mut r: Vec<f32> = b.clone();
    let mut p: Vec<f32> = r.clone(); 

    let mut ap: Vec<f32> = a(&p, pixes);
    let mut pap = p.iter().zip(ap.iter()).map(|i|  i.0 * i.1).sum::<f32>();
    
    let rr = r.iter().map(|r| r*r).sum::<f32>();
    let alpha: f32 = rr / pap;

    for i in 0..r.len(){
        r[i] = &r[i] - alpha*&ap[i]; 
    }

    let r_0: f32 = r.iter().map(|r| r*r).sum::<f32>();

    for _i in 1..maxiter {

        let rr_i: f32 = r.iter().map(|r| r*r).sum::<f32>();
        let beta = rr_i / rr;

        for i in 0..p.len(){
            p[i] = &r[i] + beta*&p[i];
        }

        let p = r.iter().zip(p.iter()).map(|f| f.0 + beta*f.1).collect::<Vec<f32>>();

        ap = a(&p, pixes);
        pap = p.iter().zip(ap.iter()).map(|i| i.0 * i.1).sum::<f32>(); 
        
        let alpha = rr_i / pap;

        for i in 0..n {
            x[i] += alpha * p[i];
        }

        for i in 0..n {
            r[i] -=  ap[i];
        }

        let rr = rr_i;

        println!("Residuals: {}", rr/r_0);
    }

    x

}
