use crate::threadpool;
use threadpool::ThreadPool;
use std::sync::mpsc;

pub fn dot_prod(a: Vec<f32>, b: Vec<f32>) -> f32{
    let my_pool = ThreadPool::new(16);
    let (tx, rx) = mpsc::channel();
    
    assert!(a.len() == b.len());
    
    for i in 0..a.len(){
        let tx = tx.clone();
        let a = a[i].clone();
        let b = b[i].clone();
        my_pool.execute(move || {
            let p = a * b;
            tx.send(p).expect("channel will be there waiting for the pool");
        });
    }

    let mut dot_res = 0.0;
    let mut cop = Vec::new();
    for _ in 0..a.len(){
        cop.push(rx.recv().unwrap());
    }
    
    for i in 0..cop.len(){
        dot_res += cop[i];
    }

    dot_res
    
}