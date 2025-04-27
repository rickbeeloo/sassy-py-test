use sassy::*;

fn main() {
    // Normal > IUPAC
    let mut seq = [b'g'; 32];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    let bases = b"ATGC";
    let mut result = vec![];
    unsafe { match_bases(&seq, bases, &mut result) };
    print_match_positions(&result, bases);

    // IUPAC > normal
    let mut seq = [b'g'; 32];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    seq[2] = b'C';
    let bases = b"Y";
    unsafe { match_bases(&seq, bases, &mut result) };
    print_match_positions(&result, bases);
}
