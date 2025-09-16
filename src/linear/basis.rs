use crate::lie::Torsor;

pub use super::group::GeneralLinearGroup;

pub type Basis<const N: usize, V = f64> = Torsor<N, GeneralLinearGroup<N, V>>;
