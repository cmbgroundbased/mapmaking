use num::Float;
//use approx;
use colored::*;

pub fn conjgrad(A: fn(f32) -> f32, b: Vec<f32>, tol: f32, maxiter: usize, precon: fn(Vec<f32>, Vec<f32>)->Vec<f32>) -> Vec<f32> {
    
    // let mut data = CGData::new(b.len());
    let _n_iter: usize = 0;

    let mut x: Vec<f32> = Vec::with_capacity(b.len());
    let mut norm_b: f32 = 0.0;

    for i in b.iter(){
        norm_b += i*i;
    }

    let norm_b = Float::sqrt(norm_b);

    // Use approx!
    if norm_b == 0.0 {
        for i in 0..b.len() {
            x[i] = 0.0;
        }
        println!("{}", "Zero Vector".yellow());
        x
    } else {

        // x = A(r)
        let mut r: Vec<f32> = Vec::with_capacity(b.len());
        for i in 0..b.len() {
            r[i] = 0.0;
        }

        for i in 0..b.len() {
            x[i] = A(r[i]);
        }

        // -r = -1.0 * r
        for i in 0..b.len() {        
            r[i] = -1.0 * r[i];
        }

        // axpy (y += a * x) axpy(1.0, b, r)
        for i in 0..b.len(){
            r[i] += 1.0 * b[i];
        }

        let r_0: f32 = Float::sqrt(r.iter().zip(r.iter()).map(|(x, y)| x * y).sum());

        if r_0 < tol {
            println!("{}", "Error 2".red());
            x
        } else {

            let mut z: Vec<f32> = Vec::with_capacity(b.len());
            let mut p: Vec<f32> = Vec::with_capacity(b.len());

            for i in 0..b.len() {
                z[i] = 0.0;
                p[i] = 0.0;
            }

            let mut z = precon(z, r.clone());

            for i in 0..b.len() {
                p[i] = z[i];
            }

            let mut r = r.clone();
            for iter in 0..maxiter {
                let mut ap: Vec<f32> = p.iter().map(|p| A(*p)).collect();
                let mut gamma:f32 = r.clone().iter().zip(z.iter()).map(|(x, y)| x * y).sum();
                let mut pAp: f32 = p.iter().zip(ap.iter()).map(|(x, y)| x * y).sum();
                let mut alpha = gamma / pAp;
                for i in 0..b.len(){
                    x[i] += alpha * p[i];
                }
                for i in 0..b.len(){
                    r[i] += -1.0 * alpha * ap[i];
                }
                
                let res: f32 = r.iter().zip(r.iter()).map(|(x, y)| x * y).sum();
                let res = res / r_0;
                let res = Float::sqrt(res);

                println!("CG ratio residuals: {}", res);
                
                z = precon(z, r.clone());

                let beta: f32 = z.iter().zip(r.iter()).map(|(x, y)| x * y).sum();
                let beta = beta / gamma;

                let p: Vec<f32> = p.iter().map(|p| beta * p).collect();

                let mut p = p.clone();

                for i in 0..b.len() {
                    p[i] += 1.0 * z[i];
                }

                // implement the auto-rerutn when the condition: res < tol is 
                // satisfied.
            }
            x
        }
    }
}