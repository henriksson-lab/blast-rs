use crate::algo::blast::api::query_data::{CBioseq, CSeqAlignSet, CSeqId};
pub struct PsiCdMsa;
pub struct PsiCdMsaCell;

pub struct PsiCdMsaCellData {
    pub wfreqs: *mut f64,
    pub iobsr: f64,
}

pub struct PsiMsaDimensions;
pub struct PsiBlastOptions;
pub struct PsiDiagnosticsRequest;

pub const K_ALPHABET_SIZE: usize = 28;

pub struct IPssmInputCdd {
    pub base: IPssmInputBase,
}

pub struct IPssmInputBase {
    pub gap_existence: i32,
    pub gap_extension: i32,
}

pub struct CCddInputData {
    pub query_data: Vec<u8>,
    pub query_title: String,
    pub db_name: String,
    pub seqalign_set: Option<Box<CSeqAlignSet>>,
    pub hits: Vec<Box<CddHit>>,
    pub cdd_data: PsiCdMsa,
    pub msa_dimensions: PsiMsaDimensions,
    pub msa_data: Vec<PsiCdMsaCell>,
    pub msa: *mut *mut PsiCdMsaCell,
    pub opts: PsiBlastOptions,
    pub matrix_name: String,
    pub diagnostics_request: *mut PsiDiagnosticsRequest,
    pub min_evalue: f64,
    pub query_bioseq: Option<Box<CBioseq>>,
    pub gap_existence: i32,
    pub gap_extension: i32,
}

impl CCddInputData {
    pub fn x_validate_hits(&self) -> bool {
        for hit in &self.hits {
            hit.validate();
        }
        true
    }
}

pub struct CddRange {
    pub from: i32,
    pub to: i32,
}

pub struct CddHitSegment {
    pub query_range: CddRange,
    pub subject_range: CddRange,
    pub msa_data: Vec<PsiCdMsaCellData>,
    pub w_freqs_data: Vec<f64>,
}

impl CddHitSegment {
    pub fn adjust_ranges(&mut self, d_from: i32, d_to: i32) {
        self.query_range.from += d_from;
        self.query_range.to += d_to;
        self.subject_range.from += d_from;
        self.subject_range.to += d_to;
    }

    pub fn is_empty(&self) -> bool {
        self.query_range.from > self.query_range.to
            || self.subject_range.from > self.subject_range.to
    }

    pub fn get_length(&self) -> i32 {
        self.query_range.to - self.query_range.from
    }

    pub fn validate(&self) -> bool {
        debug_assert!(self.query_range.from >= 0 && self.query_range.to >= 0);
        debug_assert!(self.subject_range.from >= 0 && self.subject_range.to >= 0);

        let query_length = self.query_range.to - self.query_range.from;
        let subject_length = self.subject_range.to - self.subject_range.from;
        if query_length != subject_length {
            return false;
        }

        debug_assert_eq!(
            self.w_freqs_data.len(),
            subject_length.max(0) as usize * K_ALPHABET_SIZE
        );
        debug_assert_eq!(self.msa_data.len(), subject_length.max(0) as usize);

        for data in &self.msa_data {
            debug_assert!(!data.wfreqs.is_null());
        }

        true
    }
}

pub enum CddHitApplyTo {
    Query,
    Subject,
}

pub struct CddHit {
    pub subject_id: Option<Box<CSeqId>>,
    pub evalue: f64,
    pub msa_idx: i32,
    pub segment_list: Vec<Box<CddHitSegment>>,
}

impl CddHit {
    pub fn get_length(&self) -> i32 {
        if self.is_empty() {
            return 0;
        }

        let mut result = 0;
        for segment in &self.segment_list {
            result += segment.get_length();
        }

        result
    }

    pub fn is_empty(&self) -> bool {
        if self.segment_list.is_empty() {
            return true;
        }

        for segment in &self.segment_list {
            if !segment.is_empty() {
                return false;
            }
        }

        true
    }

    pub fn validate(&self) -> bool {
        debug_assert!(self.subject_id.is_some());

        for segment in &self.segment_list {
            segment.validate();
        }

        true
    }

    pub fn intersect_with_ranges(&mut self, ranges: &[CddRange], app: CddHitApplyTo) {
        let mut range_index = 0usize;
        let mut new_segments: Vec<Box<CddHitSegment>> = Vec::new();

        for segment in std::mem::take(&mut self.segment_list) {
            if range_index >= ranges.len() {
                break;
            }

            let (seg_from, seg_to) = if matches!(&app, CddHitApplyTo::Subject) {
                (segment.subject_range.from, segment.subject_range.to)
            } else {
                (segment.query_range.from, segment.query_range.to)
            };

            while range_index < ranges.len() && ranges[range_index].to < seg_from {
                range_index += 1;
            }

            if range_index == ranges.len() {
                break;
            }

            let intersection_from = seg_from.max(ranges[range_index].from);
            let intersection_to = seg_to.min(ranges[range_index].to);

            if intersection_from == seg_from && intersection_to == seg_to {
                new_segments.push(segment);
                continue;
            }

            if intersection_to < intersection_from {
                continue;
            }

            while range_index < ranges.len() && ranges[range_index].from < seg_to {
                let d_from = seg_from.max(ranges[range_index].from) - seg_from;
                let d_to = seg_to.min(ranges[range_index].to) - seg_to;

                let mut new_seg = Box::new(CddHitSegment {
                    query_range: CddRange {
                        from: segment.query_range.from,
                        to: segment.query_range.to,
                    },
                    subject_range: CddRange {
                        from: segment.subject_range.from,
                        to: segment.subject_range.to,
                    },
                    msa_data: Vec::new(),
                    w_freqs_data: Vec::new(),
                });
                new_seg.adjust_ranges(d_from, d_to);
                debug_assert!(!new_seg.is_empty());
                new_segments.push(new_seg);

                range_index += 1;
            }
        }

        new_segments.sort_by(|a, b| a.subject_range.from.cmp(&b.subject_range.from));
        self.segment_list = new_segments;
    }

    pub fn intersect_with_hit(&mut self, hit: &CddHit, app: CddHitApplyTo) {
        let mut ranges = Vec::with_capacity(hit.segment_list.len());
        for segment in &hit.segment_list {
            let range = if matches!(&app, CddHitApplyTo::Query) {
                CddRange {
                    from: segment.query_range.from,
                    to: segment.query_range.to,
                }
            } else {
                CddRange {
                    from: segment.subject_range.from,
                    to: segment.subject_range.to,
                }
            };
            ranges.push(range);
        }

        ranges.sort_by(|a, b| a.from.cmp(&b.from).then_with(|| a.to.cmp(&b.to)));
        self.intersect_with_ranges(&ranges, app);
    }

    pub fn subtract(&mut self, hit: &CddHit) {
        if self.is_empty() || hit.is_empty() {
            return;
        }

        let from = hit.segment_list.first().unwrap().query_range.from;
        let to = hit.segment_list.last().unwrap().query_range.to;

        if self.segment_list.first().unwrap().query_range.from >= to
            || self.segment_list.last().unwrap().query_range.to <= from
        {
            return;
        }

        let mut old_segments = std::mem::take(&mut self.segment_list).into_iter();
        let mut new_segments: Vec<Box<CddHitSegment>> = Vec::new();
        let mut current = None;

        for segment in old_segments.by_ref() {
            if segment.query_range.to <= from {
                new_segments.push(segment);
            } else {
                current = Some(segment);
                break;
            }
        }

        let Some(mut segment) = current else {
            self.segment_list = new_segments;
            return;
        };

        if segment.query_range.from > to {
            new_segments.push(segment);
            new_segments.extend(old_segments);
            self.segment_list = new_segments;
            return;
        }

        if segment.query_range.to > to {
            let mut new_seg;

            if segment.query_range.from < from {
                new_seg = Box::new(CddHitSegment {
                    query_range: CddRange {
                        from: segment.query_range.from,
                        to: segment.query_range.to,
                    },
                    subject_range: CddRange {
                        from: segment.subject_range.from,
                        to: segment.subject_range.to,
                    },
                    msa_data: Vec::new(),
                    w_freqs_data: Vec::new(),
                });

                let d_to = from - segment.query_range.to;
                debug_assert!(d_to < 0);
                segment.adjust_ranges(0, d_to);
                debug_assert!(segment.query_range.from < segment.query_range.to);
                new_segments.push(segment);
            } else {
                new_seg = segment;
            }

            let d_from = to - new_seg.query_range.from;
            debug_assert!(d_from >= 0);
            new_seg.adjust_ranges(d_from, 0);
            debug_assert!(new_seg.query_range.from < new_seg.query_range.to);
            new_segments.push(new_seg);
            new_segments.extend(old_segments);
        } else {
            if segment.query_range.from < from {
                let d_to = from - segment.query_range.to;
                debug_assert!(d_to < 0);
                segment.adjust_ranges(0, d_to);
                debug_assert!(segment.query_range.from < segment.query_range.to);
                new_segments.push(segment);
            }

            let mut next_segment = None;
            for segment in old_segments.by_ref() {
                if segment.query_range.to <= to {
                    continue;
                }
                next_segment = Some(segment);
                break;
            }

            if let Some(mut segment) = next_segment {
                if segment.query_range.from < to {
                    let d_from = to - segment.query_range.from;
                    debug_assert!(d_from > 0);
                    segment.adjust_ranges(d_from, 0);
                    debug_assert!(segment.query_range.from < segment.query_range.to);
                    new_segments.push(segment);
                }

                new_segments.extend(old_segments);
            }
        }

        self.segment_list = new_segments;
    }
}
