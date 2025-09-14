use num_complex::{Complex32, Complex64};

use crate::{
    base::{BasePoly, UnivariateMarker},
    scalar_traits::{BasicScalar, NonIntegerScalar},
    storage_traits::BasicStorage,
};

pub(crate) type C<S> = <<S as BasicStorage>::T as NonIntegerScalar>::ToComplexScalar;

pub(crate) type R<S> = <C<S> as BasicScalar>::RealPartScalar;

pub type Poly32 = BasePoly<Vec<Complex32>, UnivariateMarker>;

pub type Poly64 = BasePoly<Vec<Complex64>, UnivariateMarker>;

pub type RealPoly32 = BasePoly<Vec<f32>, UnivariateMarker>;

pub type RealPoly64 = BasePoly<Vec<f64>, UnivariateMarker>;
