use sassy::*;

fn main() {
    // Normal > IUPAC
    let mut seq = [b'g'; 32];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    let mut result = [0u32; 32];
    packed_nibbles(&seq, b"", &mut result);

    // IUPAC > normal
    let mut seq = [b'g'; 32];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    seq[2] = b'C';
    let mut result = [0u32; 32];
    packed_nibbles(&seq, b"Y", &mut result);
}
