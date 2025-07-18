#![allow(non_snake_case)]

use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::profiles::{Ascii, Dna, Iupac};
use crate::search::{self, Strand};

pub enum SearcherType {
    Ascii(search::Searcher<Ascii>),
    Dna(search::Searcher<Dna>),
    Iupac(search::Searcher<Iupac>),
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Match {
    pub text_start: usize,
    pub text_end: usize,
    pub pattern_start: usize,
    pub pattern_end: usize,
    pub cost: i32,
    /// 0 = Fwd, 1 = Rc
    pub strand: u8,
}

impl From<search::Match> for Match {
    fn from(m: search::Match) -> Self {
        Match {
            text_start: m.text_start,
            text_end: m.text_end,
            pattern_start: m.pattern_start,
            pattern_end: m.pattern_end,
            cost: m.cost,
            strand: match m.strand {
                Strand::Fwd => 0,
                Strand::Rc => 1,
            },
        }
    }
}

/// Create a new `Searcher` instance.
///
/// `alphabet`: one of "ascii", "dna", "iupac" (case-insensitive).
/// `rc`: whether to also search the reverse-complement strand.
/// `alpha`: overhang parameter. Pass `NAN` to disable.
///
/// Returns a pointer to an opaque `Searcher` object, or panics on error.
#[unsafe(no_mangle)]
pub extern "C" fn sassy_searcher(
    alphabet: *const c_char,
    rc: bool,
    alpha: f32,
) -> *mut SearcherType {
    assert!(!alphabet.is_null(), "Alphabet pointer must not be null");
    let c_str = unsafe { CStr::from_ptr(alphabet) };
    let alphabet_str = c_str.to_str().unwrap();
    let alpha_opt = if alpha.is_nan() { None } else { Some(alpha) };

    let searcher = match alphabet_str.to_ascii_lowercase().as_str() {
        "ascii" => SearcherType::Ascii(search::Searcher::<Ascii>::new(rc, alpha_opt)),
        "dna" => SearcherType::Dna(search::Searcher::<Dna>::new(rc, alpha_opt)),
        "iupac" => SearcherType::Iupac(search::Searcher::<Iupac>::new(rc, alpha_opt)),
        _ => panic!("Unsupported alphabet: {}", alphabet_str),
    };

    Box::into_raw(Box::new(searcher))
}

/// Free a `Searcher` previously created with `sassy_searcher`.
#[unsafe(no_mangle)]
pub extern "C" fn sassy_searcher_free(ptr: *mut SearcherType) {
    if ptr.is_null() {
        panic!("Pointer to SearcherType must not be null");
    }
    unsafe {
        drop(Box::from_raw(ptr));
    }
}

/// Search for `pattern` in `text` allowing up to `k` edits.
///
/// `out_matches` will point to a newly allocated rust `Vec` of `Match` results. The
/// function returns the number of matches found.
/// Matches should be freed using `sassy_matches_free`.
#[unsafe(no_mangle)]
pub extern "C" fn search(
    searcher: *mut SearcherType,
    pattern: *const u8,
    pattern_len: usize,
    text: *const u8,
    text_len: usize,
    k: usize,
    out_matches: *mut *mut Match,
) -> usize {
    if searcher.is_null() || pattern.is_null() || text.is_null() || out_matches.is_null() {
        panic!("Pointers in search() must not be null");
    }

    let searcher = unsafe { &mut *searcher };
    let pattern = unsafe { slice::from_raw_parts(pattern, pattern_len) };
    let text = unsafe { slice::from_raw_parts(text, text_len) };

    let matches_vec: Vec<search::Match> = match searcher {
        SearcherType::Ascii(s) => s.search(pattern, &text, k),
        SearcherType::Dna(s) => s.search(pattern, &text, k),
        SearcherType::Iupac(s) => s.search(pattern, &text, k),
    };

    let mut c_matches: Vec<Match> = matches_vec.into_iter().map(Match::from).collect();
    // Ensure len == capacity so our free routine is safe even if it assumes this.
    c_matches.shrink_to_fit();
    let len = c_matches.len();
    let ptr_matches = c_matches.as_mut_ptr();
    std::mem::forget(c_matches);
    unsafe {
        *out_matches = ptr_matches;
    }
    len
}

/// Free a match array previously obtained from `sassy_searcher_search`.
#[unsafe(no_mangle)]
pub extern "C" fn sassy_matches_free(ptr: *mut Match, len: usize) {
    assert!(!ptr.is_null(), "Pointer to matches must not be null");
    unsafe {
        Vec::from_raw_parts(ptr, len, len);
    }
}
