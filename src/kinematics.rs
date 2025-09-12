use std::marker::PhantomData;

pub use crate::linear_space::LinearSpace;

// use crate::manifold::{LieGroup, Manifold, TangentBundle, TangentVector};
use crate::lie_group::{LieAlgebra, LieAlgebraDual, LieGroup, Manifold, Torsor};

// pub struct Velocity<'a, const D: usize, G: LieGroup<D>> {
//     position: &'a G,
//     velocity: &'a LieAlgebra<D, G>,
// }

pub trait HasVelocity {
    type Velocity;
    type Momentum;
}
pub type VelocityOf<T> = <T as HasVelocity>::Velocity;
pub type MomentumOf<T> = <T as HasVelocity>::Momentum;
pub type ForceOf<T> = MomentumOf<T>;

pub trait HasAcceleration {
    type Acceleration;
}
pub type AccelerationOf<T> = <T as HasAcceleration>::Acceleration;

pub struct Kinematics0<const N: usize, G: LieGroup<N>, From: Torsor<N, G>, To: Torsor<N, G> = From>
{
    transformation: G,
    _from: PhantomData<From>,
    _to: PhantomData<To>,
}

impl<const N: usize, G: LieGroup<N>, From: Torsor<N, G>, To: Torsor<N, G>>
    Kinematics0<N, G, From, To>
{
    pub fn zero() -> Self {
        Self {
            transformation: G::identity(),
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

pub struct Kinematics1<const N: usize, G: LieGroup<N>, From: Torsor<N, G>, To: Torsor<N, G> = From>
{
    transformation: G,
    velocity: VelocityOf<Self>,
    _from: PhantomData<From>,
    _to: PhantomData<To>,
}

impl<const N: usize, G: LieGroup<N>, From: Torsor<N, G>, To: Torsor<N, G>>
    Kinematics1<N, G, From, To>
{
    pub fn zero() -> Self {
        Self {
            transformation: G::identity(),
            velocity: VelocityOf::<Self>::zero(),
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

impl<const N: usize, G: LieGroup<N>, From: Torsor<N, G>, To: Torsor<N, G>> HasVelocity
    for Kinematics1<N, G, From, To>
{
    type Velocity = LieAlgebra<N, G>;
    type Momentum = LieAlgebraDual<N, G>;
}

// impl<const D: usize, G: LieGroup<D>, M: Torsor<D, G>> Kinematics<D, M, G> {
//     pub fn absolute(origin: M, displacement: G, velocity: LieAlgebra<D, G>) -> Self {
//         Self {
//             anchor: Anchor::Fixed { point: origin },
//             displacement,
//             velocity,
//         }
//     }

//     pub fn stationary(origin: M) -> Self {
//         Self::absolute(origin, G::identity(), LieAlgebra::zero())
//     }

//     pub fn relative_to(
//         parent: Arc<Kinematics<D, M, G>>,
//         displacement: G,
//         velocity: LieAlgebra<D, G>,
//     ) -> Self {
//         Self {
//             anchor: Anchor::Follow { parent },
//             displacement,
//             velocity,
//         }
//     }
// }
