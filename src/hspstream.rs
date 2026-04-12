//! Rust equivalent of blast_hspstream.c — HSP collection and streaming.
//! Replaces the C vtable + mutex approach with Rust's type system.

use std::sync::Mutex;

/// A single High-Scoring Pair (alignment hit).
#[derive(Debug, Clone)]
pub struct Hsp {
    pub score: i32,
    pub num_ident: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_offset: i32,
    pub query_end: i32,
    pub subject_offset: i32,
    pub subject_end: i32,
    pub context: i32,
    pub num_gaps: i32,
}

/// A list of HSPs for one query-subject pair.
#[derive(Debug, Clone)]
pub struct HspList {
    pub oid: i32,
    pub hsps: Vec<Hsp>,
    pub best_evalue: f64,
}

impl HspList {
    pub fn new(oid: i32) -> Self {
        HspList {
            oid,
            hsps: Vec::new(),
            best_evalue: f64::MAX,
        }
    }

    pub fn add_hsp(&mut self, hsp: Hsp) {
        if hsp.evalue < self.best_evalue {
            self.best_evalue = hsp.evalue;
        }
        self.hsps.push(hsp);
    }

    pub fn sort_by_score(&mut self) {
        self.hsps.sort_by(|a, b| b.score.cmp(&a.score));
    }
}

/// Results for one query: a collection of HspLists from different subjects.
#[derive(Debug, Clone)]
pub struct HitList {
    pub hsp_lists: Vec<HspList>,
}

impl HitList {
    pub fn new() -> Self {
        HitList {
            hsp_lists: Vec::new(),
        }
    }

    pub fn add_hsp_list(&mut self, list: HspList) {
        self.hsp_lists.push(list);
    }

    pub fn sort_by_evalue(&mut self) {
        self.hsp_lists.sort_by(|a, b| {
            a.best_evalue.partial_cmp(&b.best_evalue).unwrap_or(std::cmp::Ordering::Equal)
        });
    }
}

/// Complete search results for all queries.
#[derive(Debug, Clone)]
pub struct HspResults {
    pub hitlists: Vec<Option<HitList>>,
}

impl HspResults {
    pub fn new(num_queries: i32) -> Self {
        HspResults {
            hitlists: (0..num_queries).map(|_| None).collect(),
        }
    }

    pub fn sort_by_evalue(&mut self) {
        for hl in self.hitlists.iter_mut().flatten() {
            hl.sort_by_evalue();
        }
    }
}

/// Thread-safe HSP stream for collecting results during parallel search.
pub struct HspStream {
    results: Mutex<HspResults>,
}

impl HspStream {
    pub fn new(num_queries: i32) -> Self {
        HspStream {
            results: Mutex::new(HspResults::new(num_queries)),
        }
    }

    /// Write an HspList to the stream (thread-safe).
    pub fn write(&self, query_index: i32, hsp_list: HspList) {
        let mut results = self.results.lock().unwrap();
        let idx = query_index as usize;
        if idx < results.hitlists.len() {
            let hitlist = results.hitlists[idx].get_or_insert_with(HitList::new);
            hitlist.add_hsp_list(hsp_list);
        }
    }

    /// Consume the stream and return the collected results.
    pub fn into_results(self) -> HspResults {
        self.results.into_inner().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hsp_stream() {
        let stream = HspStream::new(2);

        let mut list1 = HspList::new(0);
        list1.add_hsp(Hsp {
            score: 100, num_ident: 50, bit_score: 91.5, evalue: 1e-20,
            query_offset: 0, query_end: 50, subject_offset: 76, subject_end: 126,
            context: 0, num_gaps: 0,
        });
        stream.write(0, list1);

        let mut list2 = HspList::new(5);
        list2.add_hsp(Hsp {
            score: 30, num_ident: 12, bit_score: 24.3, evalue: 1.7,
            query_offset: 10, query_end: 22, subject_offset: 100, subject_end: 112,
            context: 0, num_gaps: 0,
        });
        stream.write(0, list2);

        let results = stream.into_results();
        assert!(results.hitlists[0].is_some());
        assert!(results.hitlists[1].is_none());
        let hl = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hl.hsp_lists.len(), 2);
        assert_eq!(hl.hsp_lists[0].oid, 0);
        assert_eq!(hl.hsp_lists[1].oid, 5);
    }

    /// Helper to create an HSP with given score and evalue.
    fn make_hsp(score: i32, evalue: f64) -> Hsp {
        Hsp {
            score,
            num_ident: score / 2,
            bit_score: score as f64 * 0.9,
            evalue,
            query_offset: 0,
            query_end: score,
            subject_offset: 0,
            subject_end: score,
            context: 0,
            num_gaps: 0,
        }
    }

    #[test]
    fn test_hsp_stream_ordering() {
        // Add HSP lists for multiple subjects in arbitrary order,
        // then sort by evalue and verify OID ordering matches evalue ordering.
        let stream = HspStream::new(1);

        // Add subjects out of OID order with varying evalues
        let oids_and_evalues = vec![
            (7, 1e-5),
            (2, 1e-20),
            (5, 1e-10),
            (0, 1e-2),
            (9, 1e-30),
        ];

        for &(oid, evalue) in &oids_and_evalues {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(50, evalue));
            stream.write(0, list);
        }

        let mut results = stream.into_results();
        let hitlist = results.hitlists[0].as_mut().unwrap();

        // Sort by OID to verify we can retrieve them in OID order
        hitlist.hsp_lists.sort_by_key(|l| l.oid);
        let oids: Vec<i32> = hitlist.hsp_lists.iter().map(|l| l.oid).collect();
        assert_eq!(oids, vec![0, 2, 5, 7, 9]);

        // Also verify sort_by_evalue produces the correct evalue ordering
        hitlist.sort_by_evalue();
        for i in 1..hitlist.hsp_lists.len() {
            assert!(
                hitlist.hsp_lists[i].best_evalue >= hitlist.hsp_lists[i - 1].best_evalue,
                "HspLists should be sorted by evalue: {} >= {}",
                hitlist.hsp_lists[i].best_evalue,
                hitlist.hsp_lists[i - 1].best_evalue
            );
        }
    }

    #[test]
    fn test_hsp_stream_empty() {
        let stream = HspStream::new(3);
        let results = stream.into_results();

        // All query slots should be None
        assert_eq!(results.hitlists.len(), 3);
        for hitlist in &results.hitlists {
            assert!(hitlist.is_none(), "Empty stream should have no results");
        }
    }

    #[test]
    fn test_hsp_stream_single_subject() {
        let stream = HspStream::new(1);

        let mut list = HspList::new(42);
        list.add_hsp(make_hsp(75, 1e-15));
        list.add_hsp(make_hsp(30, 1e-3));
        stream.write(0, list);

        let results = stream.into_results();

        assert!(results.hitlists[0].is_some());
        let hitlist = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 42);
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 2);
        // best_evalue should track the minimum
        assert_eq!(hitlist.hsp_lists[0].best_evalue, 1e-15);
    }
}
