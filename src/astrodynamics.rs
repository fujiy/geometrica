
use crate::manifold::Scalar;
use crate::euclidean::{AffineSpace, AffineFrame, EuclideanSpace};


pub struct ECI<K: Scalar = f64>(AffineFrame<3, EuclideanSpace<3, K>>);

impl<K: Scalar> ECI<K> {
    pub fn as_reference() -> Self {
        ECI(EuclideanSpace::<3, K>::reference_frame())
    }
}


pub struct ECEF<K: Scalar = f64>(AffineFrame<3, EuclideanSpace<3, K>>);

impl<K: Scalar> ECEF<K> {
    
    pub fn at_time(t: K) -> Self {
        !todo!("Implement ECEF at time {}", t);
    }
}