extern crate mapmaking;
use approx;
use mapmaking::iteratorscustom::FloatIterator;
use mapmaking::conjugategradient::conjgrad;
use mapmaking::conjugategradient2::conjgrad2;
use num::ToPrimitive;
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

#[test]
fn test_conjugate_gradient() {

    let b: Vec<f32> = vec![2.0, 2.0, 2.0, 2.0];
    let pixes: Vec<Vec<i32>> = Vec::new();

    fn a() -> Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>> {
        Box::new(|_x: Vec<f32>, _y: Vec<Vec<i32>>| {
            
            let mut res: Vec<f32> = Vec::new();
            
            for _i in 0..4 {
                res.push(0.0);
            }

            for i in 0..4 {
                res[i] = 1.0; //* match i.to_f32() {Some(p) => p, None => 0.0} + 1.0; 
            }

            res  
        })   
    }
 
    fn p() -> Box<dyn Fn(Vec<f32>) -> Vec<f32>> {
        Box::new(|x| x)
    }

    let res = conjgrad(a(), b, 1E-10, 5, p(), &pixes);

    assert_eq!(res, vec![1.0, 0.5, 0.33, 0.25] );

}

#[test]
fn test_conjugate_gradient2() {

    let b: Vec<f32> = vec![1.0, 1.0, 1.0, 1.0];
    let pixes: Vec<Vec<i32>> = Vec::new();

    fn a() -> Box<dyn Fn(Vec<f32>, Vec<Vec<i32>>) -> Vec<f32>> {
        Box::new(|_x: Vec<f32>, _y: Vec<Vec<i32>>| {
            
            let mut res: Vec<f32> = Vec::new();
            
            for _i in 0..4 {
                res.push(0.0);
            }

            for i in 0..4 {
                res[i] = 1.0 * match i.to_f32() {Some(p) => p, None => 0.0} + 1.0; 
            }

            res  
        })   
    }
 
    fn p() -> Box<dyn Fn(Vec<f32>) -> Vec<f32>> {
        Box::new(|x| x)
    }

    let res = conjgrad2(a(), b, 1E-10, 5, p(), &pixes);

    assert_eq!(res, vec![1.0, 0.5, 0.33, 0.25] );

}