use geometrica::euclidean::{AffineSpace, EuclideanSpace, Motion3D};
use geometrica::kinematics::LieGroupKinematics;

type RigidBodyKinematics3D<K = f64> = LieGroupKinematics<3, Motion3D<K>>;

#[test]
fn test_rigid_body_kinematics() {
    // Create a RigidBodyKinematics instance
    // let reference = EuclideanSpace::<3>::reference_frame();
    // let kinematics =
    //     RigidBodyKinematics3D::new(Motion3D::identity(), Motion3D::LieAlgebra::identity());

    // // Test position method
    // let position = kinematics.position();
    // assert_eq!(position.raw, [0.0, 0.0, 0.0]);

    // // Test with_velocity method
    // let velocity_result = kinematics.with_velocity(|v| v.raw);
    // assert_eq!(velocity_result, [0.0, 0.0, 0.0]);
}
