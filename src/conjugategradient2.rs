pub fn conjgrad2(a: Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>>, b: Vec<f32>, tol: f32, maxiter:usize, pr: Box<dyn Fn(Vec<f32>) -> Vec<f32>> , pixes: &Vec<Vec<i32>>) -> Vec<f32> {
    let n = b.len();
    let _n_iter: usize = 0;

    let mut r: Vec<f32> = Vec::new();
    for _i in 0..n {
        r.push(0.0);
    }
    
     let mut x = a(r.clone(), pixes.clone());

    for i in 0..n {
        r[i] += 1.0 * b[i];
    }

    let r_0: f32 = r.iter().map(|i| i * i).sum();
    let z: Vec<f32> = pr(r.clone()); // **************

    let p: Vec<f32>  = z.clone();

    for _i in 1..maxiter{
        let ap = a(p.clone(), pixes.clone());

        let mut gamma: f32 = 0.0;
        for i in 0..n{
            gamma += r[i] * z[i];
        }

        let mut pap: f32 = 0.0;
        for i in 0..n {
            pap += p[i] * ap[i];
        }

        let alpha = gamma / pap;

        println!("Alpha: {}", alpha);
        println!("Gamma: {}", gamma);
        println!("PAP  : {}", pap); // Rimane sempre uguale
        


        if alpha >= f32::MAX {
            panic!();
        }

        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }

        let mut residual: f32 = 0.0;
        for i in 0..n{
            residual += r[i] * r[i];
        }

        println!("Resid: {}", residual/r_0);
        println!("*********");

        let z = pr(r.clone());

        let mut beta: f32 = 0.0;
        for i in 0..n {
            beta += z[i]*r[i];
        }
        let beta = beta / gamma;

        

        let mut new_p: Vec<f32> = Vec::new();
        for i in 0..n {
            new_p.push(beta*p[i])
        }

        let mut p: Vec<f32> = Vec::new();
        for i in 0..n {
            p.push(new_p[i] + z[i]);
        }
    }
    x

}
