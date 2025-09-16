use nalgebra as na;

use super::space::LinearSpace;

pub trait Tensor<const N: usize, const R: usize, const S: usize> {
    type V: LinearSpace<N>;
}

pub struct Tensor1_1<const N: usize, V: LinearSpace<N>> {
    pub raw: na::SMatrix<V::Field, N, N>,
}

impl<const N: usize, V: LinearSpace<N>> Tensor1_1<N, V> {
    // fn apply(&self, vector: &V) -> V::DualSpace {
    //     V::DualSpace {
    //         raw: self.raw * vector.raw,
    //     }
    // }
}

impl<const N: usize, V: LinearSpace<N>> Tensor<N, 1, 1> for Tensor1_1<N, V> {
    type V = V;
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
