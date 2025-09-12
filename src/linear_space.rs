use nalgebra as na;
use typenum::{U0, U1, U2, Unsigned};

use crate::lie_group::{
    Chart, LieAlgebra, LieAlgebraDual, LieGroup, Manifold, Scalar, Torsor, Zero,
};

// Linear Space

pub trait LinearSpace<const N: usize>: Manifold + Sized {
    type DualSpace: LinearSpace<N, DualSpace = Self, Field = Self::Field>;

    fn scale(&self, scalar: Self::Field) -> Self {
        Self::_from_raw(self._get_raw() * scalar)
    }

    fn pair_with(&self, dual: &Self::DualSpace) -> Self::Field {
        self._get_raw().dot(&dual._get_raw())
    }
    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N] {
        let raws: [na::SVector<Self::Field, N>; N] =
            std::array::from_fn(|i| basis[i]._get_raw().clone());
        let matrix = na::SMatrix::<Self::Field, N, N>::from_columns(&raws);
        let inverse = matrix.try_inverse().expect("Basis is singular");

        let inverse_t = inverse.transpose();
        std::array::from_fn(|i| Self::DualSpace::_from_raw(inverse_t.column(i).into_owned()))
    }

    fn zero() -> Self {
        Self::_from_raw(na::SVector::zeros())
    }
    fn add(&self, other: &Self) -> Self {
        Self::_from_raw(self._get_raw() + other._get_raw())
    }
    fn difference(&self, other: &Self) -> Self {
        Self::_from_raw(self._get_raw() - other._get_raw())
    }

    fn _get_raw(&self) -> &na::SVector<Self::Field, N>;
    fn _from_raw(raw: na::SVector<Self::Field, N>) -> Self;
}

// impl<'a, const N: usize, K: Scalar> std::ops::Sub<&'a $vector_name<N, K>>
//     for &'a $vector_name<N, K>
// {
//     type Output = $vector_name<N, K>;

//     fn sub(self, rhs: Self) -> Self::Output {
//         $vector_name {
//             raw: self.raw - rhs.raw,
//         }
//     }
// }
// impl<'a, const N: usize, K: Scalar> std::ops::Neg for &'a $vector_name<N, K> {
//     type Output = $vector_name<N, K>;

//     fn neg(self) -> Self::Output {
//         $vector_name { raw: -self.raw }
//     }
// }
// impl<'a, const N: usize, K: Scalar> std::ops::Mul<K> for &'a $vector_name<N, K> {
//     type Output = $vector_name<N, K>;

//     fn mul(self, scalar: K) -> Self::Output {
//         $vector_name {
//             raw: self.raw * scalar,
//         }
//     }
// }

pub trait InnerProductSpace<const N: usize>: LinearSpace<N> {
    fn flat(&self) -> Self::DualSpace;
    fn sharp(covector: &Self::DualSpace) -> Self;
    fn dot(&self, other: &Self) -> Self::Field {
        self.pair_with(&other.flat())
    }

    fn random_unit() -> Self {
        // Self::_from_raw(na::SVector::new_random())
        todo!("Implement random unit vector generation")
    }
}

trait OrientedSpace<const N: usize>: LinearSpace<N> {}

// impl<const N: usize, V: LinearSpace<N>> InnerProductSpace<N> for V
// where
//     V::DualSpace: InnerProductSpace<N>,
// {
//     fn flat(&self) -> Self::DualSpace {
//         Self::DualSpace::sharp(self)
//     }

//     fn sharp(covector: &Self::DualSpace) -> Self {
//         covector.flat()
//     }
// }

pub trait Tensor<const N: usize, V: LinearSpace<N>, R: Unsigned, S: Unsigned> {
    // fn product<R_: Unsigned, S_: Unsigned, T: Tensor<N, V, R_, S_>>(
    //     &self,
    //     other: &T,
    // ) -> dyn Tensor<N, V, R + R_, S + S_>;
}

impl<const N: usize, V: LinearSpace<N>> Tensor<N, V, U0, U0> for V::Field {}

pub struct Tensor02<const N: usize, V: LinearSpace<N>> {
    raw: na::SMatrix<V::Field, N, N>,
}
pub struct Tensor11<const N: usize, V: LinearSpace<N>> {
    raw: na::SMatrix<V::Field, N, N>,
}

pub struct Tensor20<const N: usize, V: LinearSpace<N>> {
    raw: na::SMatrix<V::Field, N, N>,
}

impl<const N: usize, V: LinearSpace<N>> Tensor<N, V, U1, U1> for Tensor11<N, V> {}
impl<const N: usize, V: LinearSpace<N>> Tensor<N, V, U0, U2> for Tensor02<N, V> {}
impl<const N: usize, V: LinearSpace<N>> Tensor<N, V, U2, U0> for Tensor20<N, V> {}

impl<const N: usize, V: LinearSpace<N>> Tensor11<N, V> {
    // pub fn from_local(
    //     input_basis: Basis<N, V>,
    //     output_basis: Basis<N, V::DualSpace>,
    //     components: &[V::Field; N],
    // ) -> Self {
    //     let input_matrix = na::SMatrix::<V::Field, N, N>::from_columns(
    //         &input_basis.basis.map(|v| v._get_raw().clone()),
    //     );
    //     let output_matrix = na::SMatrix::<V::Field, N, N>::from_columns(
    //         &output_basis.dual_basis.map(|v| v._get_raw().clone()),
    //     );
    //     let component_matrix = na::SMatrix::<V::Field, N, N>::from_column_slice(components);
    //     Tensor11 {
    //         raw: output_matrix * component_matrix * input_matrix,
    //     }
    // }

    pub fn apply(&self, vector: &V) -> V::DualSpace {
        V::DualSpace::_from_raw(self.raw * vector._get_raw())
    }
}

impl<const N: usize, V: LinearSpace<N>> Tensor<N, V, U1, U0> for V {}

impl<const N: usize, V: LinearSpace<N>> Tensor<N, V::DualSpace, U0, U1> for V {}

type Tensor2<const N: usize, V: InnerProductSpace<N>> = Tensor11<N, V>;

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

// General Linear Group GL(N)

type GL<const N: usize, V> = GeneralLinearGroup<N, V>;

pub struct GeneralLinearGroup<const N: usize, V: LinearSpace<N>> {
    raw: na::SMatrix<V::Field, N, N>,
}

impl<const N: usize, V: LinearSpace<N>> GeneralLinearGroup<N, V> {
    pub fn apply(&self, vector: &V) -> V {
        V::_from_raw(self.raw * vector._get_raw())
    }
}

macro_rules! impl_general_linear_group {
    ($N:literal, $D:literal) => {
        impl<V: LinearSpace<$N>> Manifold for GeneralLinearGroup<$N, V> {
            type Field = V::Field;
        }
        impl<V: LinearSpace<$N>> LieGroup<$D> for GeneralLinearGroup<$N, V> {
            fn identity() -> Self {
                Self {
                    raw: na::SMatrix::<V::Field, $N, $N>::identity(),
                }
            }

            fn multiply(&self, other: &Self) -> Self {
                Self {
                    raw: self.raw * other.raw,
                }
            }

            fn inverse(&self) -> Self {
                Self {
                    raw: self.raw.try_inverse().expect("Matrix is not invertible"),
                }
            }
        }

        impl<V: LinearSpace<$N>, T: Basis<$N>> Torsor<$D, GeneralLinearGroup<$N, V>> for T {}

        impl<V: LinearSpace<$N> + LieGroup<$N>> Chart<$N, V> for GeneralLinearGroup<$N, V> {
            fn to_coords(&self, of: &V) -> [V::Field; $N] {
                (self.inverse().raw * of._get_raw()).into()
            }

            fn element_from(&self, coords: &[V::Field; $N]) -> V {
                V::_from_raw(self.raw * na::SVector::<V::Field, $N>::from_column_slice(coords))
            }
        }
    };
}

impl_general_linear_group!(1, 1);
impl_general_linear_group!(2, 4);
impl_general_linear_group!(3, 9);

pub trait Basis<const N: usize> {
    type Space: LinearSpace<N>;
}

// Special Orthogonal Group SO(N)

pub type SO2<V> = SpecialOrthogonalGroup2D<V>;
pub type SO3<V> = SpecialOrthogonalGroup3D<V>;

pub struct SpecialOrthogonalGroup2D<V: LinearSpace<2>> {
    angle: V::Field,
}

impl<V: LinearSpace<2>> SpecialOrthogonalGroup2D<V> {
    fn from_angle(angle: V::Field) -> Self {
        Self { angle: angle }
    }

    fn to_angle(&self) -> V::Field {
        self.angle
    }
}

impl<V: LinearSpace<2>> Manifold for SpecialOrthogonalGroup2D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<2>> LieGroup<1> for SpecialOrthogonalGroup2D<V> {
    fn identity() -> Self {
        Self {
            angle: V::Field::zero(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Self {
            angle: self.angle + other.angle,
        }
    }

    fn inverse(&self) -> Self {
        Self { angle: -self.angle }
    }
}

impl<V: LinearSpace<2>, T: OrthonormalBasis<2>> Torsor<1, SpecialOrthogonalGroup2D<V>> for T {}

impl<V: LinearSpace<2> + LieGroup<2>> Chart<2, V> for SpecialOrthogonalGroup2D<V> {
    fn element_from(&self, coords: &[V::Field; 2]) -> V {
        let rotation_matrix = na::Rotation2::new(self.angle);
        V::_from_raw(rotation_matrix * na::SVector::<V::Field, 2>::from_column_slice(coords))
    }
    fn to_coords(&self, of: &V) -> [V::Field; 2] {
        let rotation_matrix = na::Rotation2::new(-self.angle);
        (rotation_matrix * of._get_raw()).into()
    }
}

pub struct SpecialOrthogonalGroup3D<V: LinearSpace<3>> {
    quaternion: na::UnitQuaternion<V::Field>,
}

impl<V: LinearSpace<3>> Manifold for SpecialOrthogonalGroup3D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<3>> LieGroup<3> for SpecialOrthogonalGroup3D<V> {
    fn identity() -> Self {
        Self {
            quaternion: na::UnitQuaternion::identity(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Self {
            quaternion: self.quaternion * other.quaternion,
        }
    }

    fn inverse(&self) -> Self {
        Self {
            quaternion: self.quaternion.inverse(),
        }
    }
}

impl<V: LinearSpace<3>, T: OrthonormalBasis<3>> Torsor<3, SpecialOrthogonalGroup3D<V>> for T {}

pub trait OrthonormalBasis<const N: usize> {
    type Space: InnerProductSpace<N>;
}
