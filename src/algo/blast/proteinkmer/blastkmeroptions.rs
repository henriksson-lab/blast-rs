#[derive(Debug, Clone, PartialEq)]
pub struct BlastKmerOptions {
    thresh: f64,
    min_hits: i32,
    num_target_seqs: i32,
}

impl Default for BlastKmerOptions {
    fn default() -> Self {
        Self::new()
    }
}

impl BlastKmerOptions {
    pub fn new() -> Self {
        Self {
            thresh: 0.10,
            min_hits: 1,
            num_target_seqs: 500,
        }
    }

    pub fn validate(&self) -> bool {
        if self.thresh <= 0.0 || self.thresh > 1.0 {
            return false;
        }

        if self.min_hits < 0 {
            return false;
        }

        if self.num_target_seqs < 0 {
            return false;
        }

        true
    }

    pub fn thresh(&self) -> f64 {
        self.thresh
    }

    pub fn set_thresh(&mut self, thresh: f64) {
        self.thresh = thresh;
    }

    pub fn min_hits(&self) -> i32 {
        self.min_hits
    }

    pub fn set_min_hits(&mut self, min_hits: i32) {
        self.min_hits = min_hits;
    }

    pub fn num_target_seqs(&self) -> i32 {
        self.num_target_seqs
    }

    pub fn set_num_target_seqs(&mut self, num_target_seqs: i32) {
        self.num_target_seqs = num_target_seqs;
    }
}
