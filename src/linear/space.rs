use crate::lie::Torsor;
use crate::manifold::Manifold;
use generic_array::{ArrayLength, GenericArray, sequence::GenericSequence};
use nalgebra::{OMatrix, OVector, U1};

pub use nalgebra::{DefaultAllocator, allocator::Allocator};

pub trait LinearSpace: Manifold
where
    DefaultAllocator: Allocator<Self::Dim, U1>,
{
    type DualSpace: LinearSpace<DualSpace = Self, Field = Self::Field, Dim = Self::Dim>;

    fn zero() -> Self {
        Self::_from_raw(OVector::zeros())
    }
    fn scale(&self, scalar: <Self as Manifold>::Field) -> Self {
        Self::_from_raw(self._get_raw() * scalar)
    }

    fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold>::Field {
        self._get_raw().dot(&dual._get_raw())
    }
    fn dual_basis(basis: &GenericArray<Self, Self::Dim>) -> GenericArray<Self::DualSpace, Self::Dim>
    where
        Self::Dim: ArrayLength,
        DefaultAllocator: Allocator<Self::Dim, Self::Dim>,
    {
        let raws: GenericArray<OVector<Self::Field, Self::Dim>, Self::Dim> =
            GenericArray::generate(|i| basis[i]._get_raw().clone());
        let matrix = OMatrix::<Self::Field, Self::Dim, Self::Dim>::from_columns(&raws);
        let inverse = matrix.try_inverse().expect("Basis is singular");

        let inverse_t = inverse.transpose();
        GenericArray::<Self::DualSpace, Self::Dim>::generate(|i| {
            Self::DualSpace::_from_raw(inverse_t.column(i).into_owned())
        })
    }

    // fn _induced_basis(basis: &Basis<N, Self>) -> [na::SVector<<Self as Manifold<N>>::Field, N>; N];

    fn _get_raw(&self) -> &OVector<Self::Field, Self::Dim>;
    fn _from_raw(raw: OVector<Self::Field, Self::Dim>) -> Self;
}

pub trait InnerProductSpace: LinearSpace
where
    DefaultAllocator: Allocator<Self::Dim>,
{
    fn dot(&self, other: &Self) -> Self::Field;
}

pub type AffineSpace<V: LinearSpace> = Torsor<V>;
