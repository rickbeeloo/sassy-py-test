#![allow(non_snake_case)]

use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::profiles::{Ascii, Dna, Iupac};
use crate::search::{self, Match, Strand};

#[repr(C)]
#[cbindgen::opaque]
pub enum SearcherType {
    Ascii(search::Searcher<Ascii>),
    Dna(search::Searcher<Dna>),
    Iupac(search::Searcher<Iupac>),
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CMatch {
    pub pattern_start: i32,
    pub text_start: i32,
    pub pattern_end: i32,
    pub text_end: i32,
    pub cost: i32,
    pub strand: u8, // 0 = Fwd, 1 = Rc
}

impl From<Match> for CMatch {
    fn from(m: Match) -> Self {
        CMatch {
            pattern_start: m.start.0,
            text_start: m.start.1,
            pattern_end: m.end.0,
            text_end: m.end.1,
            cost: m.cost,
            strand: match m.strand {
                Strand::Fwd => 0,
                Strand::Rc => 1,
            },
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn sassy_searcher(
    alphabet: *const c_char,
    rc: bool,
    alpha: f32,
) -> *mut SearcherType {
    if alphabet.is_null() {
        return std::ptr::null_mut();
    }
    let c_str = unsafe { CStr::from_ptr(alphabet) };
    let Ok(alphabet_str) = c_str.to_str() else {
        return std::ptr::null_mut();
    };
    let alpha_opt = if alpha.is_nan() { None } else { Some(alpha) };

    let searcher = match alphabet_str.to_ascii_lowercase().as_str() {
        "ascii" => SearcherType::Ascii(search::Searcher::<Ascii>::new(rc, alpha_opt)),
        "dna" => SearcherType::Dna(search::Searcher::<Dna>::new(rc, alpha_opt)),
        "iupac" => SearcherType::Iupac(search::Searcher::<Iupac>::new(rc, alpha_opt)),
        _ => return std::ptr::null_mut(),
    };

    Box::into_raw(Box::new(searcher))
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn sassy_searcher_free(ptr: *mut SearcherType) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(ptr));
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn search(
    ptr: *mut SearcherType,
    pattern: *const u8,
    pattern_len: usize,
    text: *const u8,
    text_len: usize,
    k: usize,
    out_matches: *mut *mut CMatch,
) -> usize {
    if ptr.is_null() || pattern.is_null() || text.is_null() || out_matches.is_null() {
        // Return MAX for error, or should we do 0?
        return usize::MAX;
    }

    let searcher = unsafe { &mut *ptr };
    let pattern_slice = unsafe { slice::from_raw_parts(pattern, pattern_len) };
    let text_slice = unsafe { slice::from_raw_parts(text, text_len) };

    let matches_vec: Vec<Match> = match searcher {
        SearcherType::Ascii(s) => s.search(pattern_slice, &text_slice, k),
        SearcherType::Dna(s) => s.search(pattern_slice, &text_slice, k),
        SearcherType::Iupac(s) => s.search(pattern_slice, &text_slice, k),
    };

    let mut c_matches: Vec<CMatch> = matches_vec.into_iter().map(CMatch::from).collect();
    // Ensure len == capacity so our free routine is safe even if it assumes this.
    c_matches.shrink_to_fit();
    let len = c_matches.len();
    let ptr_matches = c_matches.as_mut_ptr();
    // Transfer buffer ownership
    std::mem::forget(c_matches);

    unsafe {
        *out_matches = ptr_matches;
    }
    len
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn sassy_matches_free(ptr: *mut CMatch, len: usize) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        Vec::from_raw_parts(ptr, len, len);
    }
}
