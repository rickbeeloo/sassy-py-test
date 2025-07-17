use crate::search::OwnedStaticText;
use needletail::{FastxReader, parse_fastx_file};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Minimum number of reference bases to hold in batch
const DEFAULT_BATCH_BYTES: usize = 256 * 1024; // 256 KB

pub type SharedRecord = Arc<(String, OwnedStaticText)>;

#[derive(Clone)]
pub struct Query {
    pub id: String,
    pub seq: Vec<u8>,
}

pub type QueryArc = Arc<Query>;

#[derive(Clone)]
pub struct BatchItem {
    pub query: QueryArc,
    pub record: SharedRecord,
}

pub type Batch = Vec<BatchItem>;

struct RecordState {
    reader: Box<dyn FastxReader + Send>,
    /// Index of the next pattern for the *current* record.
    next_pattern_idx: usize,
    /// Current reference record, reused for all patterns until we go to the next
    current_record: Option<SharedRecord>,
}

/// Thread-safe iterator giving *batches* of ((query_id, query_seq), record) pairs.
/// each batch tries to "fill" the batch_byte_limit
pub struct RecordIterator {
    queries: Vec<QueryArc>,
    state: Mutex<RecordState>,
    batch_byte_limit: usize,
}

impl RecordIterator {
    /// Create a new iterator over `fasta_path`, going through `queries`.
    /// `max_batch_bytes` controls how many texts are bundled together.
    pub fn new<P: AsRef<Path>>(
        fasta_path: P,
        queries: Vec<Query>,
        max_batch_bytes: Option<usize>,
    ) -> Self {
        let reader = parse_fastx_file(fasta_path).expect("valid fasta");
        // Just empty state when we create the iterator
        let state = RecordState {
            reader,
            next_pattern_idx: 0,
            current_record: None,
        };
        Self {
            queries: queries.into_iter().map(Arc::new).collect(),
            state: Mutex::new(state),
            batch_byte_limit: max_batch_bytes.unwrap_or(DEFAULT_BATCH_BYTES),
        }
    }

    /// Pull the next batch, or returns None if empty
    pub fn next_batch(&self) -> Option<Batch> {
        let mut guard = self.state.lock().unwrap();
        let mut batch: Batch = Vec::new();
        let mut bytes_in_batch = 0usize;

        loop {
            // Make sure we have a current record, just so we can unwrap
            if guard.current_record.is_none() {
                match guard.reader.next() {
                    Some(Ok(rec)) => {
                        let id = String::from_utf8(rec.id().to_owned()).expect("valid UTF-8 id");
                        let seq = rec.seq().into_owned();
                        let static_text = OwnedStaticText::new(seq);
                        guard.current_record = Some(Arc::new((id, static_text)));
                        guard.next_pattern_idx = 0; // fresh record -> pattern index reset
                    }
                    Some(Err(e)) => panic!("Error reading FASTA record: {e}"),
                    None => {
                        // Done, reached end
                        return if batch.is_empty() { None } else { Some(batch) };
                    }
                }
            }

            let rec_arc = guard.current_record.as_ref().unwrap().clone();
            let rec_len = rec_arc.1.text.len();

            // If no sapce left for next record, we return current batch
            if !batch.is_empty() && bytes_in_batch + rec_len > self.batch_byte_limit {
                break; // return current batch – keep state for next call
            }

            // Add next pattern
            let pat_arc = self.queries[guard.next_pattern_idx].clone();
            let item = BatchItem {
                query: pat_arc,
                record: rec_arc.clone(),
            };
            batch.push(item);

            // Only count the reference length once per batch (we ignore pattern lengths)
            if guard.next_pattern_idx == 0 {
                bytes_in_batch += rec_len;
            }

            // Advance pattern index, but if we dont have more queries we reset to empty
            // record so we pull new sequence in next batch
            guard.next_pattern_idx += 1;
            if guard.next_pattern_idx == self.queries.len() {
                guard.current_record = None;
            }

            // If byte limit reached exactly, finish batch.
            if bytes_in_batch >= self.batch_byte_limit {
                break;
            }
        }

        Some(batch)
    }
}

// The iterator is `Send + Sync` because all its interior mutability is behind a mutex.
unsafe impl Send for RecordIterator {}
unsafe impl Sync for RecordIterator {}
