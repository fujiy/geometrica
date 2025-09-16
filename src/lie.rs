use nalgebra as na;

use crate::manifold::Manifold;

pub trait LieGroup<const N: usize>: Manifold<N> {
    fn identity() -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn inverse(&self) -> Self;
}

pub struct LieAlgebra<const N: usize, G: LieGroup<N>> {
    raw: na::SVector<G::Field, N>,
}

impl<const N: usize, G: LieGroup<N>> LieAlgebra<N, G> {
    pub fn zero() -> Self {
        LieAlgebra {
            raw: na::SVector::zeros(),
        }
    }
}

pub struct LieAlgebraDual<const N: usize, G: LieGroup<N>> {
    raw: na::SVector<G::Field, N>,
}

pub struct Torsor<const N: usize, G: LieGroup<N>> {
    from_origin: G,
}

// pub trait GroupAction<const N: usize, const D: usize, M: Manifold<N>>: LieGroup<D> {
//     fn act_on(&self, point: &M) -> M;
// }

// pub trait ChartTransformation<const N: usize, const D: usize, C: Chart<N>> {
//     type Transformed: Chart<N, M = C::M>;

//     fn transform(&self, chart: &C) -> Self::Transformed;
//     // fn between(chart: &C, transformed: &Self::Transformed) -> Self;
// }

// pub trait ChartDifference<const N: usize, const D: usize, C: Chart<N>> {
//     type Other: Chart<N, M = C::M>;

//     fn difference(chart: &C, other: &Self::Other) -> Self;
// }
