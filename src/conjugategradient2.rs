extern crate rayon;

use rayon::prelude::*;

pub fn conjgrad2(a: Box<dyn Fn(&Vec<f32>, &Vec<Vec<i32>>) -> Vec<f32>>, b: Vec<f32>, _tol: f32, maxiter:usize, _pr: Box<dyn Fn(Vec<f32>) -> Vec<f32>> , pixes: &Vec<Vec<i32>>) -> Vec<f32> {
    println!("{}", "Starting CG...");
    //let n = b.len();
 
    let mut x: Vec<f32> = vec![0.0; b.len()];
    let r: Vec<f32> = b.clone();
    let mut p: Vec<f32> = r.clone();


    
    let mut ap: Vec<f32> = a(&p, &pixes);

    //let mut pap: f32 = p.iter().zip(ap.iter()).map(|i| i.0 * i.1).sum::<f32>();
    let mut pap = p.par_iter().zip(ap.par_iter()).map(|i|  i.0 * i.1).sum::<f32>();
    
    //let rr: f32 = r.iter().map(|r| r*r).sum();
    let rr = r.par_iter().map(|r| r*r).sum::<f32>();

    println!("rr: {}", rr);

    let alpha: f32 = rr / pap;

    println!("{}", "Starting iteration...");

    // for i in 0..r.len(){
    //     r[i] = &r[i] - alpha*&ap[i]; 
    // }

    let r = r.par_iter().zip(ap.par_iter()).map(|r| r.0 - alpha*r.1).collect::<Vec<f32>>();


    let r_0: f32 = r.par_iter().map(|r| r*r).sum::<f32>();
    println!("r_0: {}", rr);

    for _i in 1..maxiter {

        let rr_i: f32 = r.par_iter().map(|r| r*r).sum::<f32>();
        println!("rr_i: {}", rr_i);
        
        let beta = rr_i / rr;
        println!("beta: {}", beta);

        for i in 0..p.len(){
            p[i] = &r[i] + beta*&p[i];
        }
        let p = r.par_iter().zip(p.par_iter()).map(|f| f.0 + beta*f.1).collect::<Vec<f32>>();

        ap = a(&p, &pixes);
        //pap = p.clone().iter().zip(ap.iter()).map(|i| i.0 * i.1).sum::<f32>(); 

        pap = p.par_iter().zip(ap.par_iter()).map(|i| i.0 * i.1).sum::<f32>(); 

        println!("pap: {}", pap);
        
        let alpha = rr_i / pap;

        // for i in 0..n {
        //     x[i] += alpha * p[i];
        // }

        x = x.par_iter().zip(p.par_iter()).map(|f| f.0 + alpha*f.1).collect::<Vec<f32>>();

        // for i in 0..n {
        //     r[i] -=  ap[i];
        // }

        let r = r.par_iter().zip(ap.par_iter()).map(|f| f.0 - f.1).collect::<Vec<f32>>();

        let rr = rr_i;

        println!("Iteration: {} -- Residuals: {} -- RR: {}", _i, r.iter().map(|r| r*r).sum::<f32>()/r_0, rr);
    }

    x

}
