use nalgebra as na;
pub use num_traits::{Num, One, Zero};

use crate::linear_space::LinearSpace;

pub trait Scalar: na::RealField + Copy + Num + Zero + One {}
impl<T: na::RealField + Copy + Num> Scalar for T {}

pub trait Manifold {
    type Field: Scalar;
}

pub trait LieGroup<const N: usize>: Manifold {
    fn identity() -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn inverse(&self) -> Self;

    fn coords_in(&self, chart: &impl Chart<N, Self>) -> [Self::Field; N]
    where
        Self: Sized,
    {
        chart.to_coords(&self)
    }
    fn new(&self, chart: &impl Chart<N, Self>, coords: &[Self::Field; N]) -> Self
    where
        Self: Sized,
    {
        chart.element_from(coords)
    }
}

pub trait Torsor<const N: usize, G: LieGroup<N>> {}

pub trait Chart<const N: usize, G: LieGroup<N>> {
    fn to_coords(&self, of: &G) -> [G::Field; N];
    fn element_from(&self, coords: &[G::Field; N]) -> G;
}

pub struct LieAlgebra<const N: usize, G: LieGroup<N>> {
    raw: na::SVector<G::Field, N>,
}

impl<const N: usize, G: LieGroup<N>> Manifold for LieAlgebra<N, G> {
    type Field = G::Field;
}

impl<const N: usize, G: LieGroup<N>> LinearSpace<N> for LieAlgebra<N, G> {
    type DualSpace = LieAlgebraDual<N, G>;

    fn _get_raw(&self) -> &na::SVector<G::Field, N> {
        &self.raw
    }

    fn _from_raw(raw: na::SVector<G::Field, N>) -> Self {
        Self { raw }
    }
}

pub struct LieAlgebraDual<const N: usize, G: LieGroup<N>> {
    raw: na::SVector<G::Field, N>,
}

impl<const N: usize, G: LieGroup<N>> Manifold for LieAlgebraDual<N, G> {
    type Field = G::Field;
}

impl<const N: usize, G: LieGroup<N>> LinearSpace<N> for LieAlgebraDual<N, G> {
    type DualSpace = LieAlgebra<N, G>;

    fn _get_raw(&self) -> &na::SVector<G::Field, N> {
        &self.raw
    }

    fn _from_raw(raw: na::SVector<G::Field, N>) -> Self {
        Self { raw }
    }
}
