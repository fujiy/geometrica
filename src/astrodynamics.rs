use crate::euclidean::*;
use crate::kinematics::*;

pub type ECI<K = f64> = OrthonormalAffineFrame<3, EuclideanSpace<3, K>>;

impl<K: Scalar> ECI<K> {
    pub fn as_reference() -> Self {
        EuclideanSpace::<3, K>::reference_frame()
    }
}

pub type ECEF<K = f64> =
    Kinematics<6, OrthonormalAffineFrame<3, EuclideanSpace<3, K>>, Motion3D<K>>;

impl<K: Scalar> ECEF<K> {
    pub fn at_time(t: K) -> Self {
        !todo!("Implement ECEF at time {}", t);
    }
}
