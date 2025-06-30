

use crate::manifold::{Manifold, TangentBundle, TangentVector, LieGroup};



pub type Kinematics<const N: usize, M> = TangentBundle<N, M>;

impl<const N: usize, M: Manifold<N>> Kinematics<N, M> {
    pub fn position(&self) -> &M {
        &self.point
    }
    pub fn with_velocity<F, R>(&self, f: F) -> R where
        F: for<'point> FnOnce(TangentVector<'point, N, M>) -> R {
        self.with_vector(f)
    }
}

pub struct LieGroupKinematics<const N: usize, G: LieGroup<N>> {
    pub point: G,
    pub velocity: G::LieAlgebra,
}

impl<const N: usize, G: LieGroup<N>> LieGroupKinematics<N, G> {
    
}
