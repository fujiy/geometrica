pub mod affine;
pub mod euclidean;
pub mod general;
pub mod orthogonal;

pub use euclidean::{DimOfSE, SE, SpecialEuclideanGroup};
pub use general::{DimOfGL, GL, GeneralLinearGroup};
pub use orthogonal::{DimOfSO, SO, SpecialOrthogonalGroup};
