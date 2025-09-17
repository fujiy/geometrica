// use crate::manifold::{LieGroup, Manifold, TangentBundle, TangentVector};
use crate::lie::{LieAlgebra, LieGroup, Torsor};
use nalgebra::{DefaultAllocator, allocator::Allocator};

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

pub struct Kinematics<G: LieGroup>
where
    DefaultAllocator: Allocator<G::Dim>,
{
    pub point: Torsor<G>,
    pub velocity: LieAlgebra<G>,
}

impl<G: LieGroup> HasVelocity for Kinematics<G>
where
    DefaultAllocator: Allocator<G::Dim>,
{
    type Velocity = LieAlgebra<G>;
}

impl<G: LieGroup> Kinematics<G>
where
    DefaultAllocator: Allocator<G::Dim>,
{
    pub fn new(point: Torsor<G>, velocity: LieAlgebra<G>) -> Self {
        Self { point, velocity }
    }

    pub fn stationary(point: Torsor<G>) -> Self {
        Self {
            point: point,
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
