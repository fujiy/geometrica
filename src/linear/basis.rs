use crate::lie::{GroupAction, LieGroup, Torsor};
use crate::linear::InnerProductSpace;
use crate::linear::group::orthogonal::SORepr;
use crate::linear::group::{DimOfGL, DimOfSO, GeneralLinearGroup, SpecialOrthogonalGroup};
use crate::linear::space::LinearSpace;
use crate::manifold::{Chart, Manifold};
use generic_array::{ArrayLength, GenericArray};
use nalgebra::{DefaultAllocator, allocator::Allocator};
use nalgebra::{Dim, OVector};

pub type LinearBasis<V = f64> = Torsor<GeneralLinearGroup<V>>;

impl<V: LinearSpace> Chart<V> for LinearBasis<V>
where
    DefaultAllocator: Allocator<V::Dim> + Allocator<V::Dim, V::Dim>,
    V::Dim: DimOfGL,
{
    fn from_local(&self, components: &GenericArray<V::Field, V::Dim>) -> V
    where
        V::Dim: ArrayLength,
    {
        self._from_origin
            .inverse()
            .act_on(&V::_from_raw(OVector::from_column_slice(
                components.as_slice(),
            )))
    }

    fn to_local(&self, point: &V) -> GenericArray<V::Field, V::Dim>
    where
        V::Dim: ArrayLength,
    {
        GenericArray::from_slice(self._from_origin.act_on(point)._get_raw().as_slice()).clone()
    }
}

pub type OrthonormalLinearBasis<V = f64> = Torsor<SpecialOrthogonalGroup<V>>;

impl<V: InnerProductSpace> Chart<V> for OrthonormalLinearBasis<V>
where
    DefaultAllocator: Allocator<V::Dim>,
    V::Dim: DimOfSO + SORepr<V>,
{
    fn from_local(&self, components: &GenericArray<V::Field, V::Dim>) -> V
    where
        V::Dim: ArrayLength,
    {
        self._from_origin
            .inverse()
            .act_on(&V::_from_raw(OVector::from_column_slice(
                components.as_slice(),
            )))
    }

    fn to_local(&self, point: &V) -> GenericArray<V::Field, V::Dim>
    where
        V::Dim: ArrayLength,
    {
        GenericArray::from_slice(self._from_origin.act_on(point)._get_raw().as_slice()).clone()
    }
}
