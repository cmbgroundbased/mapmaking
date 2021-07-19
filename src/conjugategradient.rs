pub mod linalgebra;
use linalgebra::dot_prod;
// use std::time::{Duration, Instant, SystemTime};

pub fn conjgrad(a: Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>>, 
                b: Vec<f32>, _tol: f32, maxiter:usize, 
                pr: Box<dyn Fn(Vec<f32>)->Vec<f32>>, pixs: Vec<Vec<i32>>) -> Vec<f32> {
    
   
    let mut x = vec![0.0; b.len()];
    let mut r = b.clone();
    let z = pr(r.clone()); 

    let mut p = z.clone();

    let r_0 = dot_prod(r.clone(), z.clone());
   
    for _i in 0..maxiter {

        let z = pr(r.clone());
        let gamma = dot_prod(r.clone(), z.clone());

        // 19 sec !!!
        //let now = Instant::now(); 
        let ap = a(p.clone(), pixs.clone());  // E' UN SUICIDIO!!!!!!!!!!!!!
        //println!("A(x) = {}", now.elapsed().as_secs());

        let alpha = gamma / dot_prod(p.clone(), ap.clone());
        
        for i in 0..x.len() {
            x[i] = x[i] + alpha * p[i]
        }

        for i in 0..r.len() {
            r[i] = r[i] - alpha*ap[i];
        }

        let znew = pr(r.clone());

        let beta = dot_prod(r.clone(), znew.clone()) / gamma;
        
        
        for i in 0..p.len(){
            p[i] = znew[i] + beta * p[i];
        }

        println!("Res: {}", (beta*gamma)/r_0); 
        
        if (beta*gamma/r_0) <= _tol {
            break;
        }

        //rold = rnew;

    }
    x
    
}