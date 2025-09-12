use crate::euclidean::*;
use crate::kinematics::*;

// pub type ECI<K = f64> = OrthonormalAffineFrame<3, EuclideanSpace<3, K>>;

pub struct EarthCenter;
impl EuclideanSpace<3> for EarthCenter {}

pub struct ECI;
impl OrthonormalFrame<3> for ECI {}

pub struct ECEF;
impl OrthonormalFrame<3> for ECEF {}

pub fn ecef<K: Scalar>(time: K) -> Kinematics1<6, Motion3D<K>, ECI, ECEF> {
    // !todo!("Implement ECEF at time {}", time);
    unimplemented!()
}

// impl<K: Scalar> EuclideanSpace<3, K> for ECI {}

// impl<K: Scalar> ECI<K> {
//     pub fn as_reference() -> Self {
//         EuclideanSpace::<3, K>::reference_frame()
//     }
// }

// pub type ECEF<K = f64> =
//     Kinematics<6, OrthonormalAffineFrame<3, EuclideanSpace<3, K>>, Motion3D<K>>;

// impl<K: Scalar> ECEF<K> {
//     pub fn at_time(t: K) -> Self {
//         !todo!("Implement ECEF at time {}", t);
//     }
// }
