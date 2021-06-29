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