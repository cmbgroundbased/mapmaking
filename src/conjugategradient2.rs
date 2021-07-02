pub fn conjgrad2(a: Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>>, b: Vec<f32>, tol: f32, maxiter:usize, pr: Box<dyn Fn(Vec<f32>) -> Vec<f32>> , pixes: &Vec<Vec<i32>>) -> Vec<f32> {
    let n = b.len();
 
    let mut x: Vec<f32> = vec![0.0; n];
    let mut r: Vec<f32> = b.clone();
    let mut p: Vec<f32> = r.clone();
    
    let mut ap: Vec<f32> = a(p.clone(), pixes.clone());
    let mut pap = p.clone().iter().zip(ap.iter()).map(|i| i.0 * i.1).sum::<f32>(); 
    let mut rr: f32 = r.iter().map(|r| r*r).sum();

    let mut alpha = rr / pap;

    for i in 0..r.len(){
        r[i] -= alpha*ap[i]; 
    }

    let r_0: f32 = r.iter().map(|r| r*r).sum::<f32>();

    for _i in 1..maxiter {

         
        let rr_i: f32 = r.iter().map(|r| r*r).sum();
        
        let beta = rr_i / rr;

        for i in 0..p.len(){
            p[i] = r[i] + beta*p[i];
        }

        ap = a(p.clone(), pixes.clone());
        pap = p.clone().iter().zip(ap.iter()).map(|i| i.0 * i.1).sum::<f32>(); 
        
        let alpha = rr_i / pap;

        for i in 0..n {
            x[i] += alpha * p[i];
        }

        for i in 0..n {
            r[i] -=  ap[i];
        }

        rr = rr_i.clone();

        println!("Iteration: {} -- Residuals: {}", _i, r.iter().map(|r| r*r).sum::<f32>()/r_0);
    }

    x

}
