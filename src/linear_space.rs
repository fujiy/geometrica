
use crate::manifold::{Chart, Manifold};


// Affine Space

pub trait AffineSpace<const N: usize>: Manifold<N> {
    type VectorSpace: LinearSpace<N, Field = <Self as Manifold<N>>::Field>;

    fn translate(&self, other: &Self::VectorSpace) -> Self;
    fn difference(&self, other: &Self) -> Self::VectorSpace;

    fn reference_frame() -> AffineFrame<N, Self>;
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

    fn to_local(&self, point: &Self::M) -> [<Self::M as Manifold<N>>::Field; N] {
        self.basis.to_local(&point.difference(&self.origin))
    }

    fn from_local(&self, point: &[<Self::M as Manifold<N>>::Field; N]) -> Self::M {
        self.origin.translate(&self.basis.from_local(point))
    }

    fn indiced_basis(&self, _point: &Self::M) -> [[<Self::M as Manifold<N>>::Field; N]; N] {
        self.basis.indiced_basis()
    }
}


// Linear Space

pub trait LinearSpace<const N: usize> : AffineSpace<N, VectorSpace = Self>
    {
    type DualSpace: LinearSpace<N, DualSpace=Self, Field = <Self as Manifold<N>>::Field>;

    fn zero() -> Self;
    fn scale(&self, scalar: <Self as Manifold<N>>::Field) -> Self;

    fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold<N>>::Field;
    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N];

    fn reference_basis() -> Basis<N, Self>;

    fn add(&self, other: &Self) -> Self {
        self.translate(other)
    }
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
    pub dual_basis: [V::DualSpace; N]
    // inverse_matrix: OnceCell<Tensor2<N, V>>, // Cache for the inverse matrix
}

impl<const N: usize, V: LinearSpace<N>> Basis<N, V> {
    pub fn new(vectors: [V; N]) -> Self {
        Basis {
            dual_basis: V::dual_basis(&vectors),
            basis: vectors,
        }
    }

    fn indiced_basis(&self) -> [[<V as Manifold<N>>::Field; N]; N] {
        let reference = V::reference_basis();
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                reference.dual_basis[j].pair_with(&self.basis[i])
            })
        })
    }
}

impl<const N: usize, V: LinearSpace<N>> Chart<N> for Basis<N, V> {
    type M = V;
    
    fn to_local(&self, point: &Self::M) -> [<Self::M as Manifold<N>>::Field; N] {
        std::array::from_fn(|i| {
            self.dual_basis[i].pair_with(point)
        })
    }

    fn from_local(&self, components: &[<Self::M as Manifold<N>>::Field; N]) -> Self::M {
        components.iter().enumerate().fold(
            V::zero(),
            |v, (i, c)| {
                v.add(&self.basis[i].scale(*c))
            }
        )
    }

    fn indiced_basis(&self, _point: &Self::M) -> [[<Self::M as Manifold<N>>::Field; N]; N] {
        self.indiced_basis()
    }

}
