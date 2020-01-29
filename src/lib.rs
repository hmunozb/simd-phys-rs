//#![no_std]
/// Support for SIMD-optimizable computation with stable Rust
///

pub mod aligned;
pub mod vf64;
pub mod vc64;
pub mod r3;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
