use std::f64::consts::PI;

use geometrica::{euclidean::*, manifold::ChartTransform};

#[test]
fn test_euclid3() {
    {
        let e = EuclideanSpace::<3>::reference_frame();
        let v = e.from_local(&[1.0, 2.0, 3.0]);
        let r = Rotation3D::from_axis_angle(&e.basis[0], PI / 2.0);
        let m = Motion3D::from_rotation(&r, &e.origin);
        let b = m.transform(&e);

        println!("v: {:?}", e.to_local(&v));
        println!("b: {:?}", b.to_local(&v));
    }

    // let v1 = e.clone().from_local([1.0, 2.0, 3.0]);

    // let t = Translation::<3>::zero();
    // let t2 = t.multiply(&t);
    // assert_eq!(t2.raw, [0.0, 0.0, 0.0]);
}
