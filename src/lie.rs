use crate::manifold::Manifold;
use nalgebra::OVector;
pub use nalgebra::{DefaultAllocator, allocator::Allocator};

pub trait LieGroup: Manifold {
    fn identity() -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn inverse(&self) -> Self;
}

pub struct LieAlgebra<G: LieGroup>
where
    DefaultAllocator: Allocator<G::Dim>,
{
    raw: OVector<G::Field, G::Dim>,
}

impl<G: LieGroup> LieAlgebra<G>
where
    DefaultAllocator: Allocator<G::Dim>,
{
    pub fn zero() -> Self {
        LieAlgebra {
            raw: OVector::zeros(),
        }
    }
}

pub struct LieAlgebraDual<G: LieGroup>
where
    DefaultAllocator: Allocator<G::Dim>,
{
    raw: OVector<G::Field, G::Dim>,
}

pub struct Torsor<G: LieGroup> {
    pub _from_origin: G,
}

impl<G: LieGroup> Manifold for Torsor<G> {
    type Field = G::Field;
    type Dim = G::Dim;
}

pub trait GroupAction<M: Manifold>: LieGroup {
    fn act_on(&self, point: &M) -> M;
}

impl<G: LieGroup> GroupAction<Torsor<G>> for G {
    fn act_on(&self, point: &Torsor<G>) -> Torsor<G> {
        Torsor {
            _from_origin: self.multiply(&point._from_origin),
        }
    }
}

// pub trait ChartTransformation<const N: usize, const D: usize, C: Chart<N>> {
//     type Transformed: Chart<N, M = C::M>;

//     fn transform(&self, chart: &C) -> Self::Transformed;
//     // fn between(chart: &C, transformed: &Self::Transformed) -> Self;
// }

// pub trait ChartDifference<const N: usize, const D: usize, C: Chart<N>> {
//     type Other: Chart<N, M = C::M>;

//     fn difference(chart: &C, other: &Self::Other) -> Self;
// }
