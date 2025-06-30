
use dgeom::euclidean::*;

#[test]
fn test_euclid3() {
    let e = EuclideanSpace::<3>::reference_frame();

    let v = e.from_local(&[1.0, 2.0, 3.0]);

    // let v1 = e.clone().from_local([1.0, 2.0, 3.0]);

    // let t = Translation::<3>::zero();
    // let t2 = t.multiply(&t);
    // assert_eq!(t2.raw, [0.0, 0.0, 0.0]);
}