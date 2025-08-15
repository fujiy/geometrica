use nalgebra as na;

use crate::manifold::{Chart, LieGroup, Manifold, Torsor, VectorField, ZeroVectorField};

// Affine Space

pub trait AffineSpace<const N: usize>: Torsor<N, Self::VectorSpace> {
    type VectorSpace: LinearSpace<N, Field = <Self as Manifold<N>>::Field> + LieGroup<N>;
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

// generic_const_exprs pls...
// impl<const N: usize, A: AffineSpace<N>> Manifold<{ N + N * N }> for AffineFrame<N, A> {
//    type Field = A::Field;
// }
impl<A: AffineSpace<1>> Manifold<2> for AffineFrame<1, A> {
    type Field = A::Field;
}
impl<A: AffineSpace<2>> Manifold<6> for AffineFrame<2, A> {
    type Field = A::Field;
}
impl<A: AffineSpace<3>> Manifold<12> for AffineFrame<3, A> {
    type Field = A::Field;
}

impl<const N: usize, A: AffineSpace<N>> Chart<N, A> for AffineFrame<N, A> {
    type InducedVectorField = AffineVectorField<N, A>;

    fn to_local(&self, point: &A) -> [A::Field; N] {
        self.basis.to_local(&point.difference(&self.origin))
    }

    fn from_local(&self, point: &[A::Field; N]) -> A {
        self.origin.act(&self.basis.from_local(point))
    }

    fn _induced_basis(&self, _point: &A) -> [na::SVector<A::Field, N>; N] {
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

pub struct OrthogonalAffineFrame<
    const N: usize,
    A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>,
> {
    pub origin: A,
    pub basis: OrthogonalBasis<N, A::VectorSpace>,
}

impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>>
    OrthogonalAffineFrame<N, A>
{
    pub fn new(origin: A, basis: OrthogonalBasis<N, A::VectorSpace>) -> Self {
        OrthogonalAffineFrame { origin, basis }
    }
}

impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>> Into<AffineFrame<N, A>>
    for OrthogonalAffineFrame<N, A>
{
    fn into(self) -> AffineFrame<N, A> {
        AffineFrame {
            origin: self.origin,
            basis: self.basis.into(),
        }
    }
}

// generic_const_exprs pls...
// impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>> Manifold<{ N + N * (N + 1) / 2 }> for OrthogonalAffineFrame<N, A> {
//    type Field = A::Field;
// }
impl<A: AffineSpace<1, VectorSpace: InnerProductSpace<1>>> Manifold<3>
    for OrthogonalAffineFrame<1, A>
{
    type Field = A::Field;
}
impl<A: AffineSpace<2, VectorSpace: InnerProductSpace<2>>> Manifold<6>
    for OrthogonalAffineFrame<2, A>
{
    type Field = A::Field;
}
impl<A: AffineSpace<3, VectorSpace: InnerProductSpace<3>>> Manifold<12>
    for OrthogonalAffineFrame<3, A>
{
    type Field = A::Field;
}

impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>> Chart<N, A>
    for OrthogonalAffineFrame<N, A>
{
    type InducedVectorField = AffineVectorField<N, A>;

    fn to_local(&self, point: &A) -> [A::Field; N] {
        self.basis.to_local(&point.difference(&self.origin))
    }

    fn from_local(&self, point: &[A::Field; N]) -> A {
        self.origin.act(&self.basis.from_local(point))
    }

    fn _induced_basis(&self, _point: &A) -> [na::SVector<A::Field, N>; N] {
        A::VectorSpace::_induced_basis(&self.basis.0)
    }
}

pub struct OrthonormalAffineFrame<
    const N: usize,
    A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>,
> {
    pub origin: A,
    pub basis: OrthonormalBasis<N, A::VectorSpace>,
}

impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>>
    OrthonormalAffineFrame<N, A>
{
    pub fn new(origin: A, basis: OrthonormalBasis<N, A::VectorSpace>) -> Self {
        OrthonormalAffineFrame { origin, basis }
    }
}

impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>> Into<AffineFrame<N, A>>
    for OrthonormalAffineFrame<N, A>
{
    fn into(self) -> AffineFrame<N, A> {
        AffineFrame {
            origin: self.origin,
            basis: self.basis.into(),
        }
    }
}

impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>>
    Into<OrthogonalAffineFrame<N, A>> for OrthonormalAffineFrame<N, A>
{
    fn into(self) -> OrthogonalAffineFrame<N, A> {
        OrthogonalAffineFrame {
            origin: self.origin,
            basis: self.basis.into(),
        }
    }
}

// generic_const_exprs pls...
// impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>> Manifold<{ N + N * (N - 1) / 2 }> for OrthonormalAffineFrame<N, A> {
//    type Field = A::Field;
// }
impl<A: AffineSpace<1, VectorSpace: InnerProductSpace<1>>> Manifold<1>
    for OrthonormalAffineFrame<1, A>
{
    type Field = A::Field;
}
impl<A: AffineSpace<2, VectorSpace: InnerProductSpace<2>>> Manifold<3>
    for OrthonormalAffineFrame<2, A>
{
    type Field = A::Field;
}
impl<A: AffineSpace<3, VectorSpace: InnerProductSpace<3>>> Manifold<6>
    for OrthonormalAffineFrame<3, A>
{
    type Field = A::Field;
}

impl<const N: usize, A: AffineSpace<N, VectorSpace: InnerProductSpace<N>>> Chart<N, A>
    for OrthonormalAffineFrame<N, A>
{
    type InducedVectorField = AffineVectorField<N, A>;

    fn to_local(&self, point: &A) -> [A::Field; N] {
        self.basis.to_local(&point.difference(&self.origin))
    }

    fn from_local(&self, components: &[A::Field; N]) -> A {
        self.origin.act(&self.basis.from_local(components))
    }

    fn _induced_basis(&self, _point: &A) -> [na::SVector<A::Field, N>; N] {
        A::VectorSpace::_induced_basis(&self.basis.0.0)
    }
}

// Linear Topological Space

pub trait LinearSpace<const N: usize>: AffineSpace<N, VectorSpace = Self> {
    type DualSpace: LinearSpace<N, DualSpace = Self, Field = <Self as Manifold<N>>::Field>;

    fn zero() -> Self;
    fn scale(&self, scalar: <Self as Manifold<N>>::Field) -> Self;

    fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold<N>>::Field;
    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N];

    fn add(&self, other: &Self) -> Self;

    fn _induced_basis(basis: &Basis<N, Self>) -> [na::SVector<<Self as Manifold<N>>::Field, N>; N];
}

pub trait InnerProductSpace<const N: usize>: LinearSpace<N> {
    fn dot(&self, other: &Self) -> <Self as Manifold<N>>::Field;
}

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

// generic_const_exprs pls...
// impl<const N: usize, V: LinearSpace<N>> Manifold<{ N * N }> for Basis<N, V> {
//     type Field = V::Field;
// }
impl<V: LinearSpace<1>> Manifold<1> for Basis<1, V> {
    type Field = V::Field;
}
impl<V: LinearSpace<2>> Manifold<4> for Basis<2, V> {
    type Field = V::Field;
}
impl<V: LinearSpace<3>> Manifold<9> for Basis<3, V> {
    type Field = V::Field;
}

impl<const N: usize, V: LinearSpace<N>> Chart<N, V> for Basis<N, V> {
    type InducedVectorField = ZeroVectorField<N, V>;

    fn to_local(&self, point: &V) -> [V::Field; N] {
        std::array::from_fn(|i| self.dual_basis[i].pair_with(point))
    }

    fn from_local(&self, components: &[V::Field; N]) -> V {
        components
            .iter()
            .enumerate()
            .fold(V::zero(), |v, (i, c)| v.add(&self.basis[i].scale(*c)))
    }
    fn _induced_basis(&self, _point: &V) -> [nalgebra::SVector<V::Field, N>; N] {
        V::_induced_basis(self)
    }
}

impl<const N: usize, V: LinearSpace<N>> std::ops::Index<usize> for Basis<N, V> {
    type Output = V;

    fn index(&self, index: usize) -> &Self::Output {
        &self.basis[index]
    }
}

pub struct OrthogonalBasis<const N: usize, V: InnerProductSpace<N>>(Basis<N, V>);

impl<const N: usize, V: InnerProductSpace<N>> OrthogonalBasis<N, V> {
    pub fn new(vectors: [V; N]) -> Self {
        let duals = V::dual_basis(&vectors);
        OrthogonalBasis(Basis::new_with_dual(vectors, duals))
    }
}

// generic_const_exprs pls...
// impl<const N: usize, V: InnerProductSpace<N>> Manifold<{ N * (N + 1) / 2 }> for OrthogonalBasis<N, V> {
//    type Field = V::Field;
// }
impl<V: InnerProductSpace<1>> Manifold<1> for OrthogonalBasis<1, V> {
    type Field = V::Field;
}
impl<V: InnerProductSpace<2>> Manifold<3> for OrthogonalBasis<2, V> {
    type Field = V::Field;
}
impl<V: InnerProductSpace<3>> Manifold<6> for OrthogonalBasis<3, V> {
    type Field = V::Field;
}

impl<const N: usize, V: InnerProductSpace<N>> Chart<N, V> for OrthogonalBasis<N, V> {
    type InducedVectorField = ZeroVectorField<N, V>;

    fn to_local(&self, point: &V) -> [V::Field; N] {
        self.0.to_local(point)
    }

    fn from_local(&self, components: &[V::Field; N]) -> V {
        self.0.from_local(components)
    }

    fn _induced_basis(&self, _point: &V) -> [na::SVector<V::Field, N>; N] {
        V::_induced_basis(&self.0)
    }
}

impl<const N: usize, V: InnerProductSpace<N>> std::ops::Index<usize> for OrthogonalBasis<N, V> {
    type Output = V;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0.basis[index]
    }
}

impl<const N: usize, V: InnerProductSpace<N>> Into<Basis<N, V>> for OrthogonalBasis<N, V> {
    fn into(self) -> Basis<N, V> {
        self.0
    }
}

impl<const N: usize, V: InnerProductSpace<N>> AsRef<Basis<N, V>> for OrthogonalBasis<N, V> {
    fn as_ref(&self) -> &Basis<N, V> {
        &self.0
    }
}

pub struct OrthonormalBasis<const N: usize, V: InnerProductSpace<N>>(OrthogonalBasis<N, V>);

impl<const N: usize, V: InnerProductSpace<N>> OrthonormalBasis<N, V> {
    pub fn new(vectors: [V; N]) -> Self {
        let duals = V::dual_basis(&vectors);
        OrthonormalBasis(OrthogonalBasis(Basis::new_with_dual(vectors, duals)))
    }

    pub fn basis(&self) -> &[V; N] {
        &self.0.0.basis
    }
}

// generic_const_exprs pls...
// impl<const N: usize, V: InnerProductSpace<N>> Manifold<{ N * (N - 1) / 2 }> for OrthonormalBasis<N, V> {
//     type Field = V::Field;
// }
impl<V: InnerProductSpace<1>> Manifold<1> for OrthonormalBasis<1, V> {
    type Field = V::Field;
}
impl<V: InnerProductSpace<2>> Manifold<2> for OrthonormalBasis<2, V> {
    type Field = V::Field;
}
impl<V: InnerProductSpace<3>> Manifold<3> for OrthonormalBasis<3, V> {
    type Field = V::Field;
}

impl<const N: usize, V: InnerProductSpace<N>> Chart<N, V> for OrthonormalBasis<N, V> {
    type InducedVectorField = ZeroVectorField<N, V>;

    fn to_local(&self, point: &V) -> [V::Field; N] {
        self.0.to_local(point)
    }

    fn from_local(&self, components: &[V::Field; N]) -> V {
        self.0.from_local(components)
    }

    fn _induced_basis(&self, _point: &V) -> [na::SVector<V::Field, N>; N] {
        V::_induced_basis(&self.0.0)
    }
}

impl<const N: usize, V: InnerProductSpace<N>> std::ops::Index<usize> for OrthonormalBasis<N, V> {
    type Output = V;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0.0.basis[index]
    }
}

impl<const N: usize, V: InnerProductSpace<N>> Into<Basis<N, V>> for OrthonormalBasis<N, V> {
    fn into(self) -> Basis<N, V> {
        self.0.0
    }
}

impl<const N: usize, V: InnerProductSpace<N>> Into<OrthogonalBasis<N, V>>
    for OrthonormalBasis<N, V>
{
    fn into(self) -> OrthogonalBasis<N, V> {
        self.0
    }
}

impl<const N: usize, V: InnerProductSpace<N>> AsRef<Basis<N, V>> for OrthonormalBasis<N, V> {
    fn as_ref(&self) -> &Basis<N, V> {
        &self.0.0
    }
}

impl<const N: usize, V: InnerProductSpace<N>> AsRef<OrthogonalBasis<N, V>>
    for OrthonormalBasis<N, V>
{
    fn as_ref(&self) -> &OrthogonalBasis<N, V> {
        &self.0
    }
}

macro_rules! define_vector_and_covector {
    ($vector_name:ident, $covector_name:ident) => {
        define_vector!($vector_name, $covector_name);
        define_vector!($covector_name, $vector_name);
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

        impl<const N: usize, K: Scalar> LieGroup<N> for $vector_name<N, K> {
            fn identity() -> Self {
                Self {
                    raw: na::SVector::zeros(),
                }
            }

            fn multiply(&self, other: &Self) -> Self {
                Self {
                    raw: self.raw + other.raw,
                }
            }
            fn inverse(&self) -> Self {
                Self { raw: -self.raw }
            }
        }

        impl<const N: usize, K: Scalar> Torsor<N, Self> for $vector_name<N, K> {
            fn act(&self, other: &Self) -> Self {
                Self {
                    raw: self.raw + other.raw,
                }
            }

            fn difference(&self, other: &Self) -> Self {
                Self {
                    raw: self.raw - other.raw,
                }
            }
        }

        impl<const N: usize, K: Scalar> AffineSpace<N> for $vector_name<N, K> {
            type VectorSpace = Self;
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
