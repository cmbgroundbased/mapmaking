extern crate mapmaking;
use approx;
use mapmaking::iteratorscustom::FloatIterator;
extern crate npy_derive;
extern crate npy;


#[test]
fn test_float_iterator() {
    let float_iter = FloatIterator::new(0.0, 0.9, 9);
    for i in float_iter.enumerate() {
        let n = i.0 as f32;
        approx::assert_relative_eq!(i.1, 0.0f32 + 0.1f32 * n);
    }
}

use crate::mapmaking::conjugategradient::conjgrad;

#[test]
fn test_conjugate_gradient() {
    let b = vec![1.0, 2.0];

    fn mul_m() -> Box< dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>> {
        Box::new(|_x: Vec<f32>, _p: Vec<Vec<i32>>| {
            let a = vec![4.0, 1.0, 1.0, 3.0];
            let mut x = vec![0.0; 2];
            for i in 0..2{
                for j in 0..2{
                    x[i] += a[i*2+j] * _x[j];
                }
            }
            x
        })
    }

    let pixel = vec![vec![0, 1, 2]];


    fn precc() -> Box<dyn Fn(Vec<f32>) -> Vec<f32>> {
        Box::new(|x| x.iter().map(|x| 1.0*x).collect::<Vec<f32>>())
    } 

    let res = conjgrad(mul_m(), b, 1e-5, 100, precc(), pixel);

    println!("{:?}", res);

    assert_eq!(vec![0.09090909, 0.6363636], res);
}