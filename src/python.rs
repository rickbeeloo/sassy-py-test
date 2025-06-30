use crate::profiles::{Ascii, Dna, Iupac};
use crate::search::{Match, Searcher, Strand};
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

/// A Python module implemented in Rust.
#[pymodule]
fn sassy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PySearcher>()?;
    m.add_class::<PyMatch>()?;
    m.add_function(wrap_pyfunction!(search_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(search_sequences, m)?)?;
    Ok(())
}

#[pyclass]
#[doc = "A reusable searcher object for fast sequence search."]
pub struct PySearcher {
    searcher: Option<Box<dyn std::any::Any + Send + Sync>>,
    alphabet: String,
    no_rc: bool,
}

#[pymethods]
impl PySearcher {
    #[new]
    #[pyo3(signature = (alphabet, no_rc, alpha=None))]
    fn new(alphabet: &str, no_rc: bool, alpha: Option<f64>) -> PyResult<Self> {
        let searcher = match alphabet.to_lowercase().as_str() {
            "ascii" => {
                let s = Searcher::<Ascii>::new(!no_rc, alpha.map(|a| a as f32));
                Some(Box::new(s) as Box<dyn std::any::Any + Send + Sync>)
            }
            "dna" => {
                let s = Searcher::<Dna>::new(!no_rc, alpha.map(|a| a as f32));
                Some(Box::new(s) as Box<dyn std::any::Any + Send + Sync>)
            }
            "iupac" => {
                let s = Searcher::<Iupac>::new(!no_rc, alpha.map(|a| a as f32));
                Some(Box::new(s) as Box<dyn std::any::Any + Send + Sync>)
            }
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Unsupported alphabet: {}",
                    alphabet
                )));
            }
        };

        Ok(PySearcher {
            searcher,
            alphabet: alphabet.to_string(),
            no_rc,
        })
    }

    #[pyo3(text_signature = "(self, query, text, k)")]
    #[doc = "Search for a query in a text. Returns a list of PyMatch."]
    fn search(&mut self, query: Vec<u8>, text: Vec<u8>, k: usize) -> PyResult<Vec<PyMatch>> {
        match self.alphabet.to_lowercase().as_str() {
            "ascii" => {
                let searcher = self
                    .searcher
                    .as_mut()
                    .unwrap()
                    .downcast_mut::<Searcher<Ascii>>()
                    .unwrap();
                let matches = searcher.search(&query, &text, k);
                Ok(matches.into_iter().map(PyMatch::from).collect())
            }
            "dna" => {
                let searcher = self
                    .searcher
                    .as_mut()
                    .unwrap()
                    .downcast_mut::<Searcher<Dna>>()
                    .unwrap();
                let matches = searcher.search(&query, &text, k);
                Ok(matches.into_iter().map(PyMatch::from).collect())
            }
            "iupac" => {
                let searcher = self
                    .searcher
                    .as_mut()
                    .unwrap()
                    .downcast_mut::<Searcher<Iupac>>()
                    .unwrap();
                let matches = searcher.search(&query, &text, k);
                Ok(matches.into_iter().map(PyMatch::from).collect())
            }
            _ => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Unsupported alphabet: {}",
                self.alphabet
            ))),
        }
    }
}

#[pyclass]
#[derive(Clone)]
#[doc = "A match result object."]
pub struct PyMatch {
    #[pyo3(get)]
    pub start: (i32, i32),
    #[pyo3(get)]
    pub end: (i32, i32),
    #[pyo3(get)]
    pub cost: i32,
    #[pyo3(get)]
    pub strand: String,
    #[pyo3(get)]
    pub cigar: String,
}

impl From<Match> for PyMatch {
    fn from(m: Match) -> Self {
        PyMatch {
            start: (m.start.0, m.start.1),
            end: (m.end.0, m.end.1),
            cost: m.cost,
            strand: match m.strand {
                Strand::Fwd => "forward".to_string(),
                Strand::Rc => "reverse".to_string(),
            },
            cigar: m.cigar.to_string(),
        }
    }
}

#[pyfunction]
#[pyo3(text_signature = "(query, text, k, alphabet, no_rc, alpha=None)")]
#[doc = "Search for a single query in a single text. Returns a list of PyMatch."]
fn search_sequence(
    query: Vec<u8>,
    text: Vec<u8>,
    k: usize,
    alphabet: &str,
    no_rc: bool,
    alpha: Option<f64>,
) -> PyResult<Vec<PyMatch>> {
    let mut searcher = PySearcher::new(alphabet, no_rc, alpha)?;
    searcher.search(query, text, k)
}

#[pyfunction]
#[pyo3(text_signature = "(queries, texts, k, alphabet, no_rc, alpha=None)")]
#[doc = "Search multiple queries against multiple texts. Returns a list of lists of PyMatch."]
fn search_sequences(
    queries: Vec<Vec<u8>>,
    texts: Vec<Vec<u8>>,
    k: usize,
    alphabet: &str,
    no_rc: bool,
    alpha: Option<f64>,
) -> PyResult<Vec<Vec<PyMatch>>> {
    let mut searcher = PySearcher::new(alphabet, no_rc, alpha)?;
    let mut results = Vec::new();

    for query in &queries {
        let mut query_results = Vec::new();
        for text in &texts {
            let matches = searcher.search(query.clone(), text.clone(), k)?;
            query_results.extend(matches);
        }
        results.push(query_results);
    }

    Ok(results)
}
