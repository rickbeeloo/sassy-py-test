use sassy::*;

fn main() {
    // Normal > IUPAC
    let mut seq = [b'g'; 32];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    let mut result = [0u32; 32];
    match_bases_packed_nibbles_defaults(&seq, b"", &mut result);
    print_match_positions(&result, b"ATGC");

    // IUPAC > normal
    let mut seq = [b'g'; 32];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    seq[2] = b'C';
    let mut result = [0u32; 32];
    match_bases_packed_nibbles_defaults(&seq, b"Y", &mut result);
    print_match_positions(&result, b"ATGCY");
}
