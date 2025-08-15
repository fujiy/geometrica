// use crate::manifold::{LieGroup, Manifold, TangentBundle, TangentVector};
use crate::manifold::*;

// pub struct Velocity<'a, const D: usize, G: LieGroup<D>> {
//     position: &'a G,
//     velocity: &'a LieAlgebra<D, G>,
// }

pub trait HasVelocity {
    type Velocity;
}
pub trait HasAcceleration {
    type Acceleration;
}

pub struct Kinematics<const D: usize, M: Torsor<D, G>, G: LieGroup<D>> {
    pub point: M,
    action: G,
    pub velocity: LieAlgebra<D, G>,
}

impl<const D: usize, M: Torsor<D, G>, G: LieGroup<D>> HasVelocity for Kinematics<D, M, G> {
    type Velocity = LieAlgebra<D, G>;
}

impl<const D: usize, G: LieGroup<D>, M: Torsor<D, G>> Kinematics<D, M, G> {
    pub fn new(point: M, velocity: LieAlgebra<D, G>) -> Self {
        Self {
            point: point,
            action: G::identity(),
            velocity,
        }
    }

    pub fn stationary(point: M) -> Self {
        Self {
            point: point,
            action: G::identity(),
            velocity: <Self as HasVelocity>::Velocity::zero(),
        }
    }

    // pub fn velocity_from_local<const N: usize, Q: Torsor<N, H>, H: LieGroup<N>>(
    //     &self,
    //     components: &[Q::Field; N],
    // ) -> Kinematics<N, Q, H>
    // where
    //     M: Chart<N, Q>,
    // {
    //     Kinematics {
    //         point: self.point.from_local(components),
    //         action: self.action.clone(),
    //         velocity: LieAlgebra::<N, H>::zero(),
    //     }
    // }
}
