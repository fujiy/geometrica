use nalgebra::{DefaultAllocator, DimName, OVector, U1, U2, U3, allocator::Allocator};

use crate::linear::{SpecialEuclideanGroup, SpecialOrthogonalGroup};
// use crate::linear::space::impl_vector_ops;

pub use crate::lie::{LieGroup, Torsor};
pub use crate::linear::{AffineSpace, InnerProductSpace, LinearSpace};
pub use crate::manifold::{Chart, Manifold, Scalar};

// Vector space associated with Euclidean space

pub struct Vector<N: DimName, K: Scalar = f64>
where
    DefaultAllocator: Allocator<N, U1>,
{
    raw: OVector<K, N>,
}

impl<N: DimName, K: Scalar> Manifold for Vector<N, K>
where
    DefaultAllocator: Allocator<N, U1>,
{
    type Field = K;
    type Dim = N;
}
impl<N: DimName, K: Scalar> LinearSpace for Vector<N, K>
where
    DefaultAllocator: Allocator<N, U1>,
{
    type DualSpace = Covector<N, K>;

    fn _get_raw(&self) -> &OVector<K, N> {
        &self.raw
    }

    fn _from_raw(raw: OVector<K, N>) -> Self {
        Self { raw }
    }
}
impl<N: DimName, K: Scalar> InnerProductSpace for Vector<N, K>
where
    DefaultAllocator: Allocator<N, U1>,
{
    fn dot(&self, other: &Self) -> <Self as Manifold>::Field {
        self.raw.dot(&other.raw)
    }
}

pub struct Covector<N: DimName, K: Scalar = f64>
where
    DefaultAllocator: Allocator<N, U1>,
{
    raw: OVector<K, N>,
}

impl<N: DimName, K: Scalar> Manifold for Covector<N, K>
where
    DefaultAllocator: Allocator<N, U1>,
{
    type Field = K;
    type Dim = N;
}

impl<N: DimName, K: Scalar> LinearSpace for Covector<N, K>
where
    DefaultAllocator: Allocator<N, U1>,
{
    type DualSpace = Vector<N, K>;

    fn _get_raw(&self) -> &OVector<K, N> {
        &self.raw
    }

    fn _from_raw(raw: OVector<K, N>) -> Self {
        Self { raw }
    }
}

// Euclidean space

pub type EuclideanSpace<N: DimName, K = f64> = AffineSpace<Vector<N, K>>;

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

pub type Translation<N: DimName, K = f64> = Vector<N, K>;

// Rotation group or special orthogonal group SO(N)

// pub trait Rotation<const N: usize, K: Scalar = f64>: Manifold<N> + LieGroup<N> {}

pub type Rotation<N: DimName, K = f64> = SpecialOrthogonalGroup<Vector<N, K>>;
pub type Rotation2D<K = f64> = Rotation<U2, K>;
pub type Rotation3D<K = f64> = Rotation<U3, K>;

// Motion group or special Euclidean group SE(N)

pub type Motion<N: DimName, K = f64> = SpecialEuclideanGroup<Vector<N, K>>;
pub type Motion2D<K = f64> = Motion<U2, K>;
pub type Motion3D<K = f64> = Motion<U3, K>;
