pub fn conjgrad(a: Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>>, b: Vec<f32>, tol: f32, maxiter:usize, pr: Box<dyn Fn(Vec<f32>) -> Vec<f32>> , pixes: &Vec<Vec<i32>>) -> Vec<f32> {
  
    let mut x: Vec<f32> = b.iter().map(|f| *f * 0.0).collect();
    let mut r = b.clone();
    let n: usize = b.len();


    let z = pr(r.clone()); //*************** */
    let mut _rz: f32 = 0.0;

    for i in 0..n {
        _rz += r[i]*z[i];
    }

    //let rz: f32 = r.clone().iter().zip(z.iter()).map(|c| c.0 * c.1).sum();

    let rz0 = _rz;
    let p: Vec<f32> = z;
    let _n_iter: usize = 0;

 
    for _i in 0..maxiter {
        let ap: Vec<f32> = a(p.clone(), pixes.clone());
        let pap: f32 = p.clone().iter().zip(ap.iter()).map(|c| c.0 * c.1).sum();
        let alpha = _rz / pap;


        println!("Alpha: {}", alpha);

        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }

        let z = pr(r.clone());

        let next_rz: f32 = r.clone().iter().zip(z.iter()).map(|c| c.0 * c.1).sum();
        let err: f32 = next_rz / rz0;
        
        /* DEBUG
        *******************************/
        println!("CG error: {}", err);
        /******************************
        */
        
        let beta: f32 = next_rz / _rz; 
        let _rz: f32 = next_rz;
        let mut p = p.clone();
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    x
}
