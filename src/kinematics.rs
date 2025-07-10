// use crate::manifold::{LieGroup, Manifold, TangentBundle, TangentVector};
use crate::manifold::*;

pub type Kinematics<const N: usize, M> = TangentBundle<N, M>;

impl<const N: usize, M: Manifold<N>> Kinematics<N, M> {
    pub fn with_velocity<F, R>(&self, f: F) -> R
    where
        F: for<'point> FnOnce(TangentVector<'point, N, M>) -> R,
    {
        self.with_vector(f)
    }
}

pub struct LieGroupKinematics<const N: usize, G: LieGroup<N>> {
    pub point: G,
    pub velocity: G::LieAlgebra,
}

impl<const N: usize, G: LieGroup<N>> LieGroupKinematics<N, G> {
    pub fn new(point: G, velocity: G::LieAlgebra) -> Self {
        Self { point, velocity }
    }

    pub fn zero() -> Self {
        Self {
            point: G::identity(),
            velocity: G::LieAlgebra::zero(),
        }
    }
}

pub struct MovingChart<const N: usize, C: Chart<N>> {
    pub chart: C,
    pub velocity: C::InducedVectorField,
}
