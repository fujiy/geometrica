use nalgebra as na;

use crate::linear::group::{SpecialEuclideanGroup2D, SpecialEuclideanGroup3D};
use crate::linear::group::{SpecialOrthogonalGroup2D, SpecialOrthogonalGroup3D};
// use crate::linear::space::impl_vector_ops;

pub use crate::lie::{LieGroup, Torsor};
pub use crate::linear::{AffineSpace, InnerProductSpace, LinearSpace};
pub use crate::manifold::{Chart, Manifold, Scalar};

// Vector space associated with Euclidean space

pub struct Vector<const N: usize, K = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K: Scalar> Manifold<N> for Vector<N, K> {
    type Field = K;
}
impl<const N: usize, K: Scalar> LinearSpace<N> for Vector<N, K> {
    type DualSpace = Covector<N, K>;

    fn _get_raw(&self) -> &na::SVector<K, N> {
        &self.raw
    }

    fn _from_raw(raw: na::SVector<K, N>) -> Self {
        Self { raw }
    }
}
impl<const N: usize, K: Scalar> InnerProductSpace<N> for Vector<N, K> {
    fn dot(&self, other: &Self) -> <Self as Manifold<N>>::Field {
        self.raw.dot(&other.raw)
    }
}

crate::impl_vector_ops!(Vector);
pub struct Covector<const N: usize, K = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K: Scalar> Manifold<N> for Covector<N, K> {
    type Field = K;
}
impl<const N: usize, K: Scalar> LinearSpace<N> for Covector<N, K> {
    type DualSpace = Vector<N, K>;

    fn _get_raw(&self) -> &na::SVector<K, N> {
        &self.raw
    }

    fn _from_raw(raw: na::SVector<K, N>) -> Self {
        Self { raw }
    }
}

// Euclidean space

pub type EuclideanSpace<const N: usize, K = f64> = AffineSpace<N, Vector<N, K>>;

// impl<const N: usize, K: Scalar> EuclideanSpace<N, K> {
//     pub fn reference_frame() -> OrthonormalAffineFrame<N, Self> {
//         OrthonormalAffineFrame::new(
//             EuclideanSpace {
//                 raw: na::SVector::from_element(K::zero()),
//             },
//             OrthonormalBasis::new(std::array::from_fn(|i| {
//                 let mut repr = na::SVector::from_element(K::zero());
//                 repr[i] = K::one();
//                 Vector { raw: repr }
//             })),
//         )
//     }
// }

// pub struct PolarChart<const N: usize, K: Scalar = f64> {
//     base: AffineFrame<N, EuclideanSpace<N, K>>,
// }

// Translation group T(N)

pub type Translation<const N: usize, K = f64> = Vector<N, K>;

// Rotation group or special orthogonal group SO(N)

// pub trait Rotation<const N: usize, K: Scalar = f64>: Manifold<N> + LieGroup<N> {}

pub type Rotation2D<K = f64> = SpecialOrthogonalGroup2D<Vector<2, K>>;
pub type Rotation3D<K = f64> = SpecialOrthogonalGroup3D<Vector<3, K>>;

// Motion group or special Euclidean group SE(N)

pub type Motion2D<K = f64> = SpecialEuclideanGroup2D<Vector<2, K>>;
pub type Motion3D<K = f64> = SpecialEuclideanGroup3D<Vector<3, K>>;
