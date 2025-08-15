use geometrica::euclidean::*;
use geometrica::kinematics::*;

type BodyFrame3D<K> = OrthonormalAffineFrame<3, EuclideanSpace<3, K>>;
type RigidBodyKinematics3D<K> = Kinematics<6, BodyFrame3D<K>, Motion3D<K>>;

struct Value {
    raw: f64,
}
struct HasRef<'a> {
    value: &'a Value,
}

impl HasRef<'_> {
    fn freeze(&self) -> f64 {
        self.value.raw
    }
}

#[test]
fn test_rigid_body_kinematics() {
    let reference = EuclideanSpace::<3>::reference_frame();
    let body = RigidBodyKinematics3D::<f64>::stationary(reference);

    // let vel: RigidBodyKinematics3D<f64>::Velocity = body.velocity;

    let mut value: Value = Value { raw: 0.0 };

    let has_ref = HasRef { value: &value };

    // value.raw = 1.0;

    print!("Value: {}", has_ref.value.raw);

    let raw = has_ref.freeze();

    value.raw = 2.0;

    println!("Raw value: {}", raw);
    println!("Value after freeze: {}", value.raw);

    // let body = kinematics.transform(&reference);
}
