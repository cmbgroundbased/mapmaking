extern crate mapmaking;
use approx;
extern crate npy_derive;
extern crate npy;
//use mapmaking::directory::DirStruct;
use mapmaking::iteratorscustom::FloatIterator;

#[test]
fn test_float_iterator() {
    let float_iter = FloatIterator::new(0.0, 0.9, 9);
    for i in float_iter.enumerate() {
        let n = i.0 as f64;
        approx::assert_relative_eq!(i.1, 0.0f64 + 0.1f64 * n);
    }
}