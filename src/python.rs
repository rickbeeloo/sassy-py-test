use crate::profiles::{Ascii, Dna, Iupac};
use crate::search::{Match, Searcher, Strand};
use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "sassy")]
fn sassy(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySearcher>()?;
    Ok(())
}

enum SearcherType {
    Ascii(Searcher<Ascii>),
    Dna(Searcher<Dna>),
    Iupac(Searcher<Iupac>),
}

#[pyclass]
#[doc = "A reusable searcher object for fast sequence search."]
pub struct PySearcher {
    searcher: SearcherType,
}

#[pymethods]
impl PySearcher {
    #[new]
    #[pyo3(signature = (alphabet, rc=true, alpha=None))]
    fn new(alphabet: &str, rc: bool, alpha: Option<f64>) -> PyResult<Self> {
        let searcher = match alphabet.to_lowercase().as_str() {
            "ascii" => {
                let s = Searcher::<Ascii>::new(false, alpha.map(|a| a as f32));
                SearcherType::Ascii(s)
            }
            "dna" => {
                let s = Searcher::<Dna>::new(rc, alpha.map(|a| a as f32));
                SearcherType::Dna(s)
            }
            "iupac" => {
                let s = Searcher::<Iupac>::new(rc, alpha.map(|a| a as f32));
                SearcherType::Iupac(s)
            }
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Unsupported alphabet: {}",
                    alphabet
                )));
            }
        };

        Ok(PySearcher { searcher })
    }

    #[pyo3(signature = (query, text, k))]
    #[doc = "Search for a query in a text. Returns a list of PyMatch."]
    fn search(&mut self, query: Vec<u8>, text: Vec<u8>, k: usize) -> Vec<Match> {
        match &mut self.searcher {
            SearcherType::Ascii(searcher) => searcher.search(&query, &text, k),
            SearcherType::Dna(searcher) => searcher.search(&query, &text, k),
            SearcherType::Iupac(searcher) => searcher.search(&query, &text, k),
        }
    }
}

#[pymethods]
impl Match {
    #[getter]
    fn pattern_start(&self) -> i32 {
        self.start.0
    }

    #[getter]
    fn text_start(&self) -> i32 {
        self.start.1
    }

    #[getter]
    fn pattern_end(&self) -> i32 {
        self.end.0
    }

    #[getter]
    fn text_end(&self) -> i32 {
        self.end.1
    }

    #[getter]
    fn cost(&self) -> i32 {
        self.cost
    }

    #[getter]
    fn strand(&self) -> &'static str {
        match self.strand {
            Strand::Fwd => "+",
            Strand::Rc => "-",
        }
    }

    #[getter]
    fn cigar(&self) -> String {
        self.cigar.to_string()
    }
}
