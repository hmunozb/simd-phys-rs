pub mod fields;
use alga::general::{ClosedAdd};
use alga::general::{ClosedMul, ClosedSub};
use core::ops::{Sub, Mul};
use crate::vf64::Aligned4xf64;
use num_traits::Zero;
use num_traits::float::FloatCore;

pub fn dot_product<T>(
    a: &[T; 3], b: &[T; 3])  -> T
where T: Copy + ClosedMul<T> + ClosedAdd<T>
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

pub fn norm_sq<T>( a: &[T; 3] ) -> T
where T: Copy + ClosedMul<T> + ClosedAdd<T>
{
    dot_product(a, a)
}

pub fn cross_product<T>(
    a: &[T; 3],
    b: &[T; 3],
    c: &mut [T; 3]
)
where T: Copy + ClosedMul<T> + ClosedSub<T>
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/// Analytically evaluate the matrix exponential of the antisymmetric matrix represented
/// in the SO(3) Lie algebra basis by a.
/// Equivalently: the magnitude of a is the angle of rotation about the unit vector in the direction of a.
/// and the corresponding rotation matrix is written to exp_a in row-major order
///
pub fn cross_exponential_3d(
    a: &[Aligned4xf64; 3],
    exp_a: &mut [Aligned4xf64; 9]
){
    //let mut omega_mat = *e;


    let (ax_sq, ay_sq, az_sq) = (a[0]*a[0], a[1]*a[1], a[2]*a[2]);
    let a_sq = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    let mut a_norm = Aligned4xf64::default();
    for (v_norm, &v_sq) in a_norm.dat.iter_mut().zip(a_sq.dat.iter()){
        *v_norm = v_sq.sqrt();
    }
    //todo: use parallelizable versions of cos, sin
    let rc_a_norm = a_norm.map(|v| v.recip());

    let mut c_arr = Aligned4xf64::zero();
    let mut s_arr = Aligned4xf64::zero();
    for( (c, s), &a ) in c_arr.dat.iter_mut().zip(s_arr.dat.iter_mut()).zip(a_norm.dat.iter()){
        let sc = f64::sin_cos(a);
        *s = sc.0; *c = sc.1;
    }
//    let c_arr = a_norm.map(f64::cos);
//    let s_arr = a_norm.map(f64::sin);
    //protect against small-norm division
    let a_norm = a_norm.map(|v| if v.abs() < f64::epsilon(){ 1.0} else { v});

    let mut omega_mat : [Aligned4xf64; 9] = [Aligned4xf64::zero(); 9];
    omega_mat[1] = -a[2] ; omega_mat[2] =  a[1];
    omega_mat[3] =  a[2] ; omega_mat[5] = -a[0];
    omega_mat[6] = -a[1] ; omega_mat[7] =  a[0];

    let mut omega_sq_mat : [Aligned4xf64; 9] = [Aligned4xf64::zero(); 9];
    omega_sq_mat[0] = -(ay_sq + az_sq);  omega_sq_mat[1] = a[0] * a[1];  omega_sq_mat[2] = a[0] * a[2];
    omega_sq_mat[3] = omega_sq_mat[1]; omega_sq_mat[4] = -(ax_sq + az_sq); omega_sq_mat[5] = a[1] * a[2];
    omega_sq_mat[6] = omega_sq_mat[2]; omega_sq_mat[7] = omega_sq_mat[5]; omega_sq_mat[8] = -(ax_sq + ay_sq);


    let sinc_a_norm = s_arr * rc_a_norm;
    for v in omega_mat.iter_mut(){
        *v *= sinc_a_norm;
    }
    let cosc_a_norm =(-c_arr + 1.0) * rc_a_norm * rc_a_norm;
    for v in omega_sq_mat.iter_mut(){
        *v *= cosc_a_norm;
    }

    //let mut exp_a = omega_mat;
    exp_a[0] += 1.0; exp_a[4] += 1.0; exp_a[8] += 1.0;
    for (v, w) in exp_a.iter_mut().zip(omega_sq_mat.iter()){
        *v += *w;
    }

}

#[cfg(test)]
mod tests{
    use std;

    use crate::vf64::Aligned4xf64;
    use num_traits::{Zero, FloatConst};
    use crate::r3::cross_exponential_3d;

    #[test]
    fn test_cross_exponential_3d(){
        println!("Testing...");
        let a = [Aligned4xf64::zero(), Aligned4xf64::zero(),
            Aligned4xf64::from(1.0*f64::PI())];
        let mut e = [Aligned4xf64::zero(); 9];

        cross_exponential_3d(&a, &mut e);
        for b in e.chunks_exact(3){
            println!("{}\t{}\t{}", b[0].dat[0], b[1].dat[0], b[2].dat[0]);
        }
    }
}