// use geometrica::affine_space::*;
use geometrica::euclidean::*;
use geometrica::kinematics::*;

#[test]
fn test_particle_motion() {
    struct Origin;
    struct Particle;
    impl EuclideanSpace<3> for Origin {}
    impl EuclideanSpace<3> for Particle {}

    struct Reference;
    impl OrthonormalFrame<3> for Reference {}

    type ParticleKinematics = Kinematics1<3, Vector<3>, Origin, Particle>;

    let mut particle = ParticleKinematics::zero();

    let force: ForceOf<ParticleKinematics> = MomentumOf::<ParticleKinematics>::zero();

    let frame = Kinematics0::<6, Motion3D<f64>, Reference>::zero();

    // let particle = kinematics.transform(&reference);
}

// #[test]
// fn test_rigid_body_rotation() {
//     struct Body;

//     type BodyFrame3D<K = f64> = OrthonormalBasis<3, Vector<3, K>>;
//     type RigidBodyKinematics3D<K = f64> = Kinematics<3, BodyFrame3D<K>, Rotation3D<K>>;

//     let reference = EuclideanSpace::<3>::reference_basis();
//     let body = RigidBodyKinematics3D::stationary(reference);

//     let force: ForceOf<RigidBodyKinematics3D> = MomentumOf::<RigidBodyKinematics3D>::zero();

//     // let body = kinematics.transform(&reference);
// }

// #[test]
// fn test_rigid_body_motion() {
//     type BodyFrame3D<K = f64> = OrthonormalAffineFrame<3, EuclideanSpace<3, K>>;
//     type RigidBodyKinematics3D<K = f64> = Kinematics<6, BodyFrame3D<K>, Motion3D<K>>;

//     let reference = EuclideanSpace::<3>::reference_frame();
//     let body = RigidBodyKinematics3D::stationary(reference);

//     // let v: &VelocityOf<RigidBodyKinematics3D> = &body.velocity;

//     // let motion = body.action;
//     // let position = body.point;
//     // let velocity = body.velocity;
// }
