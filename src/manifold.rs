use std::marker::PhantomData;
use num_traits::{Num, Zero, One};
use nalgebra as na;

use crate::linear_space::{AffineSpace, AffineFrame, LinearSpace, Basis,};

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
    pub fn with_vector<F, R>(&self, f: F) -> R where
        F: for<'point> FnOnce(TangentVector<'point, N, M>) -> R {
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
    pub fn with_vector<F, R>(&self, f: F) -> R where
        F: for<'point> FnOnce(CotangentVector<'point, N, M>) -> R {
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

impl<'point, const N: usize, M: Manifold<N>> TangentVector<'point, N, M> {
    
}

impl<'point, const N: usize, M: Manifold<N>> Manifold<N> for TangentVector<'point, N, M> {
    type Field = M::Field;
}

impl<'point, const N: usize, M: Manifold<N>> AffineSpace<N> for TangentVector<'point, N, M> {
    type VectorSpace = Self; 

    fn translate(&self, other: &Self::VectorSpace) -> Self {
        TangentVector {
            raw: self.raw + other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn difference(&self, other: &Self) -> Self::VectorSpace {
        TangentVector {
            raw: self.raw - other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn reference_frame() -> AffineFrame<N, Self> {
        AffineFrame {
            origin: TangentVector { 
                raw: na::SVector::from_element(M::Field::zero()),
                _marker: InvariantLifetime(PhantomData),
            },
            basis: Self::VectorSpace::reference_basis(),
        }
    }
}

impl<'point, const N: usize, M: Manifold<N>> LinearSpace<N> for TangentVector<'point, N, M> {
    type DualSpace = CotangentVector<'point, N, M>;

    fn zero() -> Self {
        TangentVector {
            raw: na::SVector::zero(),
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn add(&self, other: &Self) -> Self {
        TangentVector {
            raw: self.raw + other.raw,
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
    fn reference_basis() -> Basis<N, Self> {
        let vectors = std::array::from_fn(|i| {
            let mut repr = na::SVector::from_element(M::Field::zero());
            repr[i] = M::Field::one();
            TangentVector { raw: repr, _marker: InvariantLifetime(PhantomData) }
        });
        Basis {
            dual_basis: Self::dual_basis(&vectors),
            basis: vectors, 
        }
    }
}


// Cotangent Vector

pub struct CotangentVector<'point, const N: usize, M: Manifold<N>> {
    raw: na::SVector<M::Field, N>,
    _marker: InvariantLifetime<'point>,
}

impl<'point, const N: usize, M: Manifold<N>> CotangentVector<'point, N, M> {
    
}

impl<'point, const N: usize, M: Manifold<N>> Manifold<N> for CotangentVector<'point, N, M> {
    type Field = M::Field;
}

impl<'point, const N: usize, M: Manifold<N>> AffineSpace<N> for CotangentVector<'point, N, M> {
    type VectorSpace = Self;

    fn translate(&self, other: &Self::VectorSpace) -> Self {
        CotangentVector {
            raw: self.raw + other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn difference(&self, other: &Self) -> Self::VectorSpace {
        CotangentVector {
            raw: self.raw - other.raw,
            _marker: InvariantLifetime(PhantomData),
        }
    }

    fn reference_frame() -> AffineFrame<N, Self> {
        AffineFrame {
            origin: CotangentVector { 
                raw: na::SVector::from_element(M::Field::zero()),
                _marker: InvariantLifetime(PhantomData),
            },
            basis: Self::VectorSpace::reference_basis(),
        }
    }
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
    fn reference_basis() -> Basis<N, Self> {
        let vectors = std::array::from_fn(|i| {
            let mut repr = na::SVector::from_element(M::Field::zero());
            repr[i] = M::Field::one();
            CotangentVector { raw: repr, _marker: InvariantLifetime(PhantomData) }
        });
        Basis {
            dual_basis: Self::dual_basis(&vectors),
            basis: vectors, 
        }
    }

}






pub trait Chart<const N: usize> {
    type M: Manifold<N>;

    fn to_local(&self, point: &Self::M) -> [<Self::M as Manifold<N>>::Field; N];
    fn from_local(&self, components: &[<Self::M as Manifold<N>>::Field; N]) -> Self::M;
    fn indiced_basis(&self, point: &Self::M) -> [[<Self::M as Manifold<N>>::Field; N]; N];

    fn with_induced_basis<F, R>(&self, point_and_vector: &TangentBundle<N, Self::M>, f: F) -> R where 
        F: for<'point> FnOnce(&Basis<N, TangentVector<'point, N, Self::M>>, TangentVector<'point, N, Self::M>) -> R {
        point_and_vector.with_vector(|local| {
            let basis_raw = self.indiced_basis(&point_and_vector.point);
            let marker = InvariantLifetime(PhantomData);
            let base_vectors = std::array::from_fn(|i| {
                let mut raw = na::SVector::from_element(<Self::M as Manifold<N>>::Field::zero());
                for j in 0..N {
                    raw[j] = basis_raw[i][j];
                }
                TangentVector { raw, _marker: marker }
            });
            let basis = Basis::new(base_vectors);
            f(&basis, local)
        })
    }
}




pub trait LieGroup<const N: usize>: Manifold<N> {
    type LieAlgebra;

    fn identity() -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn inverse(&self) -> Self;
}


pub trait GroupAction<const N: usize, M: Manifold<N>>: LieGroup<N> {
    fn act_on(&self, point: &M) -> M;
}

