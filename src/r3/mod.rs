use alga::general::Ring;
use alga::general::{ClosedMul, ClosedSub};
use core::ops::{Sub, Mul};

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