use nalgebra as na;

use crate::manifold::{Chart, Manifold, VectorField, ZeroVectorField};

// Affine Space

pub trait AffineSpace<const N: usize>: Manifold<N> {
    type VectorSpace: LinearSpace<N, Field = <Self as Manifold<N>>::Field>;

    fn translate(&self, other: &Self::VectorSpace) -> Self;
    fn difference(&self, other: &Self) -> Self::VectorSpace;
}

// Affine Frame for affine space

pub struct AffineFrame<const N: usize, A: AffineSpace<N>> {
    pub origin: A,
    pub basis: Basis<N, A::VectorSpace>,
}

impl<const N: usize, A: AffineSpace<N>> AffineFrame<N, A> {
    pub fn new(origin: A, basis: Basis<N, A::VectorSpace>) -> Self {
        AffineFrame { origin, basis }
    }
}

impl<const N: usize, A: AffineSpace<N>> Chart<N> for AffineFrame<N, A> {
    type M = A;
    type InducedVectorField = AffineVectorField<N, A>;

    fn to_local(&self, point: &Self::M) -> [<Self::M as Manifold<N>>::Field; N] {
        self.basis.to_local(&point.difference(&self.origin))
    }

    fn from_local(&self, point: &[<Self::M as Manifold<N>>::Field; N]) -> Self::M {
        self.origin.translate(&self.basis.from_local(point))
    }

    fn _induced_basis(
        &self,
        _point: &Self::M,
    ) -> [na::SVector<<Self::M as Manifold<N>>::Field, N>; N] {
        A::VectorSpace::_induced_basis(&self.basis)
    }
}

pub struct AffineVectorField<const N: usize, A: AffineSpace<N>> {
    pub affine_group: AffineGroup<N, A>,
}

impl<const N: usize, A: AffineSpace<N>> VectorField<N> for AffineVectorField<N, A> {
    type M = A;

    fn at(&self, point: Self::M) -> crate::manifold::TangentBundle<N, Self::M> {
        !todo!("Implement tangent bundle for affine vector field");
    }
}

struct AffineGroup<const N: usize, A: AffineSpace<N>> {
    pub frame: AffineFrame<N, A>,
}

// Linear Space

pub trait LinearSpace<const N: usize>: AffineSpace<N, VectorSpace = Self> {
    type DualSpace: LinearSpace<N, DualSpace = Self, Field = <Self as Manifold<N>>::Field>;

    fn zero() -> Self;
    fn scale(&self, scalar: <Self as Manifold<N>>::Field) -> Self;

    fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold<N>>::Field;
    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N];

    fn add(&self, other: &Self) -> Self {
        self.translate(other)
    }

    fn _induced_basis(basis: &Basis<N, Self>) -> [na::SVector<<Self as Manifold<N>>::Field, N>; N];
}

pub trait Tensor<const N: usize, const R: usize, const S: usize> {
    type V: LinearSpace<N>;
}

// struct Tensor2<const N: usize, K: Scalar> {
//     raw: na::SMatrix<K, N, N>,
// }

// impl<const N: usize, K: Scalar> Tensor2<N, K> {
//     fn from_columns(columns: [&Vector<N, K>; N]) -> Self {
//         let raw = na::SMatrix::<K, N, N>::from_columns(&columns.map(|v| v.raw));
//         Tensor2 {
//             raw,
//         }
//     }

//     fn inverse(&self) -> Self {
//         Tensor2 {
//             raw: self.raw.try_inverse().expect("Matrix is not invertible"),
//         }
//     }
// }

// impl<const N: usize, V: LinearSpace<N>> Tensor<N, 2, 0> for Tensor2<N, V> {
//     type V = V;
// }

// Basis for Linear Space

pub struct Basis<const N: usize, V: LinearSpace<N>> {
    pub basis: [V; N],
    pub dual_basis: [V::DualSpace; N], // inverse_matrix: OnceCell<Tensor2<N, V>>, // Cache for the inverse matrix
}

impl<const N: usize, V: LinearSpace<N>> Basis<N, V> {
    pub fn new(vectors: [V; N]) -> Self {
        Basis {
            dual_basis: V::dual_basis(&vectors),
            basis: vectors,
        }
    }
    pub fn new_with_dual(vectors: [V; N], duals: [V::DualSpace; N]) -> Self {
        Basis {
            basis: vectors,
            dual_basis: duals,
        }
    }
}

impl<const N: usize, V: LinearSpace<N>> Chart<N> for Basis<N, V> {
    type M = V;
    type InducedVectorField = ZeroVectorField<N, V>;

    fn to_local(&self, point: &Self::M) -> [<Self::M as Manifold<N>>::Field; N] {
        std::array::from_fn(|i| self.dual_basis[i].pair_with(point))
    }

    fn from_local(&self, components: &[<Self::M as Manifold<N>>::Field; N]) -> Self::M {
        components
            .iter()
            .enumerate()
            .fold(V::zero(), |v, (i, c)| v.add(&self.basis[i].scale(*c)))
    }
    fn _induced_basis(
        &self,
        _point: &Self::M,
    ) -> [nalgebra::SVector<<Self::M as Manifold<N>>::Field, N>; N] {
        V::_induced_basis(self)
    }
}

impl<const N: usize, V: LinearSpace<N>> std::ops::Index<usize> for Basis<N, V> {
    type Output = V;

    fn index(&self, index: usize) -> &Self::Output {
        &self.basis[index]
    }
}

macro_rules! define_vector_and_covector {
    ($vector_name:ident, $covector_name:ident) => {
        define_vector!($vector_name, $covector_name);
        define_vector!($covector_name, $vector_name);
    };
}

macro_rules! impl_vector_and_covector {
    ($vector_name:ident, $covector_name:ident) => {
        impl_vector!($vector_name, $covector_name);
        impl_vector!($covector_name, $vector_name);
    };
}

macro_rules! define_vector {
    ($vector_name:ident, $covector_name:ident) => {
        #[derive(Clone)]
        pub struct $vector_name<const N: usize, K: Scalar = f64> {
            pub raw: na::SVector<K, N>,
        }

        impl_vector!($vector_name, $covector_name);
    };
}

macro_rules! impl_vector {
    ($vector_name:ident, $covector_name:ident) => {
        impl<const N: usize, K: Scalar> Manifold<N> for $vector_name<N, K> {
            type Field = K;
        }

        impl<const N: usize, K: Scalar> AffineSpace<N> for $vector_name<N, K> {
            type VectorSpace = Self;

            fn translate(&self, other: &Self) -> Self {
                Self {
                    raw: self.raw + other.raw,
                }
            }

            fn difference(&self, other: &Self) -> Self::VectorSpace {
                Self {
                    raw: self.raw - other.raw,
                }
            }
        }

        impl<const N: usize, K: Scalar> LinearSpace<N> for $vector_name<N, K> {
            type DualSpace = $covector_name<N, K>;

            fn zero() -> Self {
                Self {
                    raw: na::SVector::zeros(),
                }
            }

            fn add(&self, other: &Self) -> Self {
                Self {
                    raw: &self.raw + &other.raw,
                }
            }

            fn scale(&self, scalar: K) -> Self {
                Self {
                    raw: &self.raw * scalar,
                }
            }

            fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold<N>>::Field {
                self.raw.dot(&dual.raw)
            }

            fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N] {
                let raws: [na::SVector<K, N>; N] = std::array::from_fn(|i| basis[i].raw.clone());
                let matrix = na::SMatrix::<K, N, N>::from_columns(&raws);
                let inverse = matrix.try_inverse().expect("Basis is singular");

                let inverse_t = inverse.transpose();
                std::array::from_fn(|i| $covector_name {
                    raw: inverse_t.column(i).into_owned(),
                })
            }

            fn _induced_basis(basis: &Basis<N, Self>) -> [na::SVector<K, N>; N] {
                std::array::from_fn(|i| basis.basis[i].raw)
            }
        }

        impl<'a, const N: usize, K: Scalar> std::ops::Add<&'a $vector_name<N, K>>
            for &'a $vector_name<N, K>
        {
            type Output = $vector_name<N, K>;

            fn add(self, rhs: Self) -> Self::Output {
                $vector_name {
                    raw: self.raw + rhs.raw,
                }
            }
        }
        impl<'a, const N: usize, K: Scalar> std::ops::Sub<&'a $vector_name<N, K>>
            for &'a $vector_name<N, K>
        {
            type Output = $vector_name<N, K>;

            fn sub(self, rhs: Self) -> Self::Output {
                $vector_name {
                    raw: self.raw - rhs.raw,
                }
            }
        }
        impl<'a, const N: usize, K: Scalar> std::ops::Neg for &'a $vector_name<N, K> {
            type Output = $vector_name<N, K>;

            fn neg(self) -> Self::Output {
                $vector_name { raw: -self.raw }
            }
        }
        impl<'a, const N: usize, K: Scalar> std::ops::Mul<K> for &'a $vector_name<N, K> {
            type Output = $vector_name<N, K>;

            fn mul(self, scalar: K) -> Self::Output {
                $vector_name {
                    raw: self.raw * scalar,
                }
            }
        }
    };
}
