use crate::vf64::Aligned4xf64;
use crate::r3::{dot_product, norm_sq};
use num_traits::Zero;

///
/// The cartesian harmonic potential
///   u(v1, v2) = (1/2) k |v1 - v2|^2
pub fn harmonic_cartesian(
    v1: &[Aligned4xf64; 3], v2: &[Aligned4xf64; 3], k: f64
) -> Aligned4xf64{
    let mut u = Aligned4xf64::zero();
    for (&x1i, &x2i) in v1.iter().zip(v2.iter()){
        let dxi = x1i - x2i;
        u += dxi*dxi;
    };
    u *= 0.5 * k;

    u
}

/// Evaluate the angular harmonic potential
///       u(v1, v2) = (1/2) k (theta - theta0)^2
/// where cos(theta) = v1.v2 / (|v1| |v2|), and theta0 is within [0, pi]
pub fn harmonic_angular(
    v1: &[Aligned4xf64; 3], v2: &[Aligned4xf64; 3], theta0: f64, k: f64
){
    let mut u = Aligned4xf64::zero();
    let dot = dot_product(v1, v2);
    let (n1, n2) = (norm_sq(v1).map(f64::sqrt), norm_sq(v2).map(f64::sqrt));
    let cos_theta = dot * (n1 * n2).map(f64::recip );
    let theta = cos_theta.map(f64::acos);


}

