use nalgebra as na;
use num_traits::{Num, One, Zero};
use std::marker::PhantomData;

use crate::linear_space::{AffineSpace, Basis};

pub use crate::linear_space::LinearSpace;

pub trait Scalar: na::RealField + Copy + Num + Zero + One {}
impl<T: na::RealField + Copy + Num> Scalar for T {}

#[derive(Clone, Copy, Default)]
struct InvariantLifetime<'id>(PhantomData<*mut &'id ()>);

pub trait Manifold<const N: usize>: Sized {
    type Field: Scalar;
}

// Tangent Bundle

pub struct TangentBundle<const N: usize, M: Manifold<N>> {
    pub point: M,
    vector_raw: na::SVector<M::Field, N>,
}

impl<const N: usize, M: Manifold<N>> TangentBundle<N, M> {
    pub fn with_vector<F, R>(&self, f: F) -> R
    where
        F: for<'point> FnOnce(TangentVector<'point, N, M>) -> R,
    {
        let local = TangentVector {
            raw: self.vector_raw,
            _marker: InvariantLifetime(PhantomData),
        };
        f(local)
    }
}

// Cotangent Bundle

pub struct CotangentBundle<const N: usize, M: Manifold<N>> {
    pub point: M,
    vector_raw: na::SVector<M::Field, N>,
}

impl<const N: usize, M: Manifold<N>> CotangentBundle<N, M> {
    pub fn with_vector<F, R>(&self, f: F) -> R
    where
        F: for<'point> FnOnce(CotangentVector<'point, N, M>) -> R,
    {
        let local = CotangentVector {
            raw: self.vector_raw,
            _marker: InvariantLifetime(PhantomData),
        };
        f(local)
    }
}

// Tangent Vector

pub struct TangentVector<'point, const N: usize, M: Manifold<N>> {
    raw: na::SVector<M::Field, N>,
    _marker: InvariantLifetime<'point>,
}

impl<'point, const N: usize, M: Manifold<N>> TangentVector<'point, N, M> {}

impl<'point, const N: usize, M: Manifold<N>> Manifold<N> for TangentVector<'point, N, M> {
    type Field = M::Field;
}

impl<'point, const N: usize, M: Manifold<N>> LieGroup<N> for TangentVector<'point, N, M> {
    fn identity() -> Self {
        TangentVector {
            raw: na::SVector::zeros(),
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        TangentVector {
            raw: self.raw + other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn inverse(&self) -> Self {
        TangentVector {
            raw: -self.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }
}

impl<'point, const N: usize, M: Manifold<N>> Torsor<N, TangentVector<'point, N, M>>
    for TangentVector<'point, N, M>
{
    fn act(&self, group_element: &Self) -> Self {
        TangentVector {
            raw: self.raw + group_element.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn difference(&self, other: &Self) -> Self {
        TangentVector {
            raw: self.raw - other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }
}

impl<'point, const N: usize, M: Manifold<N>> AffineSpace<N> for TangentVector<'point, N, M> {
    type VectorSpace = Self;
}

impl<'point, const N: usize, M: Manifold<N>> LinearSpace<N> for TangentVector<'point, N, M> {
    type DualSpace = CotangentVector<'point, N, M>;

    fn zero() -> Self {
        TangentVector {
            raw: na::SVector::zero(),
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn scale(&self, scalar: M::Field) -> Self {
        TangentVector {
            raw: self.raw * scalar,
            _marker: InvariantLifetime(PhantomData),
        }
    }
    fn pair_with(&self, dual: &Self::DualSpace) -> M::Field {
        self.raw.dot(&dual.raw)
    }
    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N] {
        let raws: [na::SVector<M::Field, N>; N] = std::array::from_fn(|i| basis[i].raw.clone());
        let matrix = na::SMatrix::<M::Field, N, N>::from_columns(&raws);
        let inverse = matrix.try_inverse().expect("Basis is singular");

        let inverse_t = inverse.transpose();
        std::array::from_fn(|i| CotangentVector {
            raw: inverse_t.column(i).into_owned(),
            _marker: InvariantLifetime(PhantomData),
        })
    }

    fn add(&self, other: &Self) -> Self {
        TangentVector {
            raw: self.raw + other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn _induced_basis(
        basis: &Basis<N, Self>,
    ) -> [nalgebra::SVector<<Self as Manifold<N>>::Field, N>; N] {
        std::array::from_fn(|i| basis.basis[i].raw)
    }
}

// Cotangent Vector

pub struct CotangentVector<'point, const N: usize, M: Manifold<N>> {
    raw: na::SVector<M::Field, N>,
    _marker: InvariantLifetime<'point>,
}

impl<'point, const N: usize, M: Manifold<N>> CotangentVector<'point, N, M> {}

impl<'point, const N: usize, M: Manifold<N>> Manifold<N> for CotangentVector<'point, N, M> {
    type Field = M::Field;
}

impl<'point, const N: usize, M: Manifold<N>> LieGroup<N> for CotangentVector<'point, N, M> {
    fn identity() -> Self {
        CotangentVector {
            raw: na::SVector::zeros(),
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        CotangentVector {
            raw: self.raw + other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn inverse(&self) -> Self {
        CotangentVector {
            raw: -self.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }
}

impl<'point, const N: usize, M: Manifold<N>> Torsor<N, CotangentVector<'point, N, M>>
    for CotangentVector<'point, N, M>
{
    fn act(&self, group_element: &Self) -> Self {
        CotangentVector {
            raw: self.raw + group_element.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn difference(&self, other: &Self) -> Self {
        CotangentVector {
            raw: self.raw - other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }
}

impl<'point, const N: usize, M: Manifold<N>> AffineSpace<N> for CotangentVector<'point, N, M> {
    type VectorSpace = Self;
}

impl<'point, const N: usize, M: Manifold<N>> LinearSpace<N> for CotangentVector<'point, N, M> {
    type DualSpace = TangentVector<'point, N, M>;

    fn zero() -> Self {
        CotangentVector {
            raw: na::SVector::zero(),
            _marker: InvariantLifetime(PhantomData),
        }
    }
    fn add(&self, other: &Self) -> Self {
        CotangentVector {
            raw: self.raw + other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }
    fn scale(&self, scalar: M::Field) -> Self {
        CotangentVector {
            raw: self.raw * scalar,
            _marker: InvariantLifetime(PhantomData),
        }
    }
    fn pair_with(&self, dual: &Self::DualSpace) -> M::Field {
        self.raw.dot(&dual.raw)
    }
    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N] {
        let raws: [na::SVector<M::Field, N>; N] = std::array::from_fn(|i| basis[i].raw.clone());
        let matrix = na::SMatrix::<M::Field, N, N>::from_columns(&raws);
        let inverse = matrix.try_inverse().expect("Basis is singular");
        let inverse_t = inverse.transpose();
        std::array::from_fn(|i| TangentVector {
            raw: inverse_t.column(i).into_owned(),
            _marker: InvariantLifetime(PhantomData),
        })
    }

    fn _induced_basis(
        basis: &Basis<N, Self>,
    ) -> [nalgebra::SVector<<Self as Manifold<N>>::Field, N>; N] {
        std::array::from_fn(|i| basis.basis[i].raw)
    }
}

pub trait VectorField<const N: usize> {
    type M: Manifold<N>;

    fn at(&self, point: Self::M) -> TangentBundle<N, Self::M>;
}
pub trait CovectorField<const N: usize> {
    type M: Manifold<N>;

    fn at(&self, point: Self::M) -> CotangentBundle<N, Self::M>;
}

pub struct ZeroVectorField<const N: usize, M: Manifold<N>> {
    _marker: PhantomData<M>,
}

impl<const N: usize, M: Manifold<N>> ZeroVectorField<N, M> {
    pub fn new() -> Self {
        ZeroVectorField {
            _marker: PhantomData,
        }
    }
}

impl<const N: usize, M: Manifold<N>> VectorField<N> for ZeroVectorField<N, M> {
    type M = M;

    fn at(&self, point: Self::M) -> TangentBundle<N, Self::M> {
        TangentBundle {
            point: point,
            vector_raw: na::SVector::zeros(),
        }
    }
}
pub struct ZeroCovectorField<const N: usize, M: Manifold<N>> {
    _marker: PhantomData<M>,
}

impl<const N: usize, M: Manifold<N>> ZeroCovectorField<N, M> {
    pub fn new() -> Self {
        ZeroCovectorField {
            _marker: PhantomData,
        }
    }
}

impl<const N: usize, M: Manifold<N>> CovectorField<N> for ZeroCovectorField<N, M> {
    type M = M;
    fn at(&self, point: Self::M) -> CotangentBundle<N, Self::M> {
        CotangentBundle {
            point: point,
            vector_raw: na::SVector::zeros(),
        }
    }
}

pub trait Chart<const N: usize, M: Manifold<N>> {
    type InducedVectorField: VectorField<N, M = M>;

    fn to_local(&self, point: &M) -> [M::Field; N];
    fn from_local(&self, components: &[M::Field; N]) -> M;

    fn _induced_basis(&self, point: &M) -> [na::SVector<M::Field, N>; N];

    fn from_local_vector(
        &self,
        point: M,
        vector_components: &[M::Field; N],
    ) -> TangentBundle<N, M> {
        let basis_raw = self._induced_basis(&point);
        let basis: Basis<N, TangentVector<'_, N, M>> =
            Basis::new(std::array::from_fn(|i| TangentVector {
                raw: basis_raw[i],
                _marker: InvariantLifetime(PhantomData),
            }));
        let vector = basis.from_local(vector_components);
        TangentBundle {
            point: point,
            vector_raw: vector.raw,
        }
    }

    fn with_local_vector<F, R>(&self, point_and_vector: &TangentBundle<N, M>, f: F) -> R
    where
        F: for<'point> FnOnce(&'point [M::Field; N]) -> R,
    {
        point_and_vector.with_vector(|local| {
            let basis_raw = self._induced_basis(&point_and_vector.point);
            let basis: Basis<N, TangentVector<'_, N, M>> =
                Basis::new(std::array::from_fn(|i| TangentVector {
                    raw: basis_raw[i],
                    _marker: InvariantLifetime(PhantomData),
                }));
            let vector_components = basis.to_local(&local);
            f(&vector_components)
        })
    }

    // fn transform<'a, const D: usize, G: GroupAction<N, D, Self::M>>(
    //     &'a self,
    //     transformation: G,
    // ) -> TransformedChart<'a, N, D, Self, G> {
    //     TransformedChart {
    //         original_chart: self,
    //         transformation,
    //     }
    // }
}

pub trait LieGroup<const N: usize>: Manifold<N> {
    fn identity() -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn inverse(&self) -> Self;
}

pub struct LieAlgebra<const N: usize, G: LieGroup<N>> {
    raw: na::SVector<G::Field, N>,
}

impl<const N: usize, G: LieGroup<N>> LieAlgebra<N, G> {
    pub fn zero() -> Self {
        LieAlgebra {
            raw: na::SVector::zeros(),
        }
    }
}

pub struct LieAlgebraDual<const N: usize, G: LieGroup<N>> {
    raw: na::SVector<G::Field, N>,
}

pub trait Torsor<const N: usize, G: LieGroup<N>>: Manifold<N> {
    fn act(&self, group_element: &G) -> Self;
    fn difference(&self, other: &Self) -> G;
}

// pub trait GroupAction<const N: usize, const D: usize, M: Manifold<N>>: LieGroup<D> {
//     fn act_on(&self, point: &M) -> M;
// }

// pub trait ChartTransformation<const N: usize, const D: usize, C: Chart<N>> {
//     type Transformed: Chart<N, M = C::M>;

//     fn transform(&self, chart: &C) -> Self::Transformed;
//     // fn between(chart: &C, transformed: &Self::Transformed) -> Self;
// }

// pub trait ChartDifference<const N: usize, const D: usize, C: Chart<N>> {
//     type Other: Chart<N, M = C::M>;

//     fn difference(chart: &C, other: &Self::Other) -> Self;
// }
