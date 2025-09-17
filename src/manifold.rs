use generic_array::{ArrayLength, GenericArray};
use nalgebra::{DimName, RealField};
pub use num_traits::{Num, One, Zero};
pub use std::ops::{Add, Neg};

pub trait Scalar: RealField + Copy + Num + Zero + One + Neg {}
impl<T: RealField + Copy + Num + Zero + One + Neg> Scalar for T {}
pub trait Manifold: Sized {
    type Field: Scalar;
    type Dim: DimName;
}

pub trait Chart<M: Manifold> {
    // type InducedVectorField: VectorField<N, M = M>;

    fn to_local(&self, point: &M) -> GenericArray<M::Field, M::Dim>
    where
        M::Dim: ArrayLength;
    fn from_local(&self, components: &GenericArray<M::Field, M::Dim>) -> M
    where
        M::Dim: ArrayLength;
}
