/// Builds a 'profile' of `b` in `64`-bit blocks, and compressed `a` into a `[0,1,2,3]` alphabet.
///
/// Returns a bitpacked `B` indicating which chars of `b` equal a given char of `a`.
///
/// TODO: 32bit variant.
pub trait Profile: Clone + Copy + std::fmt::Debug {
    /// Encoding for a single character in `a`.
    type A;
    /// Encoding for 64 characters in `b`.
    type B;

    fn encode_a(a: &[u8]) -> Vec<Self::A>;
    fn encode_b(b: &[u8; 64]) -> Self::B;
    /// Given the encoding of an `a` and the encoding for 64 `b`s,
    /// return a bitmask of which characters of `b` equal the corresponding character of `a`.
    fn eq(ca: &Self::A, cb: &Self::B) -> u64;
}
