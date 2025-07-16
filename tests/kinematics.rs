use geometrica::euclidean::{AffineSpace, EuclideanSpace, Motion3D};
use geometrica::kinematics::LieGroupKinematics;
use geometrica::*;

type RigidBodyKinematics3D<K: Scalar> = LieGroupKinematics<6, Motion3D<K>>;

#[test]
fn test_rigid_body_kinematics() {
    let reference = EuclideanSpace::<3>::reference_frame();
    let kinematics = RigidBodyKinematics3D::<f64>::zero();

    let body = kinematics.transform(&reference);
}
