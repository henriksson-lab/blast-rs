//! Rust equivalent of blast_filter.c — query sequence masking/filtering.
//! Implements DUST (nucleotide) and SEG (protein) low-complexity filtering.

use std::collections::VecDeque;

/// A masked region in a sequence.
#[derive(Debug, Clone, Copy)]
pub struct MaskedRegion {
    pub start: i32,
    pub end: i32,
}

/// Collection of masked regions for a query.
#[derive(Debug, Clone, Default)]
pub struct MaskLoc {
    pub regions: Vec<MaskedRegion>,
}

impl MaskLoc {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, start: i32, end: i32) {
        self.regions.push(MaskedRegion { start, end });
    }

    /// Check if a position is masked.
    pub fn is_masked(&self, pos: i32) -> bool {
        self.regions.iter().any(|r| pos >= r.start && pos < r.end)
    }

    /// Apply masking to a sequence by replacing masked positions with sentinel value.
    pub fn apply(&self, sequence: &mut [u8], sentinel: u8) {
        for region in &self.regions {
            for pos in region.start..region.end {
                if (pos as usize) < sequence.len() {
                    sequence[pos as usize] = sentinel;
                }
            }
        }
    }
}

const DUST_DEFAULT_LEVEL: u32 = 20;
const DUST_DEFAULT_WINDOW: usize = 64;
const DUST_DEFAULT_LINKER: usize = 1;
const TRIPLET_MASK: u8 = 0x3f;

#[derive(Clone, Copy)]
struct PerfectInterval {
    start: usize,
    end: usize,
    score: u32,
    len: usize,
}

struct DustTriplets {
    triplet_list: VecDeque<u8>,
    perfects: Vec<PerfectInterval>,
    start: usize,
    stop: usize,
    max_size: usize,
    low_k: u8,
    l: usize,
    thresholds: Vec<u32>,
    c_w: [u8; 64],
    c_v: [u8; 64],
    r_w: u32,
    r_v: u32,
    num_diff: u32,
}

impl DustTriplets {
    fn new(window: usize, low_k: u8, thresholds: Vec<u32>) -> Self {
        Self {
            triplet_list: VecDeque::new(),
            perfects: Vec::new(),
            start: 0,
            stop: 0,
            max_size: window.saturating_sub(2),
            low_k,
            l: 0,
            thresholds,
            c_w: [0; 64],
            c_v: [0; 64],
            r_w: 0,
            r_v: 0,
            num_diff: 0,
        }
    }

    fn add_triplet_info(r: &mut u32, counts: &mut [u8; 64], t: u8) {
        *r += counts[t as usize] as u32;
        counts[t as usize] += 1;
    }

    fn rem_triplet_info(r: &mut u32, counts: &mut [u8; 64], t: u8) {
        counts[t as usize] -= 1;
        *r -= counts[t as usize] as u32;
    }

    fn needs_processing(&self) -> bool {
        let count = self.stop - self.l;
        count < self.triplet_list.len() && 10 * self.r_w > self.thresholds[count]
    }

    fn shift_high(&mut self, t: u8) -> bool {
        let s = self.triplet_list.pop_back().unwrap();
        Self::rem_triplet_info(&mut self.r_w, &mut self.c_w, s);
        if self.c_w[s as usize] == 0 {
            self.num_diff -= 1;
        }
        self.start += 1;

        self.triplet_list.push_front(t);
        if self.c_w[t as usize] == 0 {
            self.num_diff += 1;
        }
        Self::add_triplet_info(&mut self.r_w, &mut self.c_w, t);
        self.stop += 1;

        if self.num_diff <= 1 {
            self.perfects.insert(
                0,
                PerfectInterval {
                    start: self.start,
                    end: self.stop + 1,
                    score: 0,
                    len: 0,
                },
            );
            false
        } else {
            true
        }
    }

    fn shift_window(&mut self, t: u8) -> bool {
        if self.triplet_list.len() >= self.max_size {
            if self.num_diff <= 1 {
                return self.shift_high(t);
            }

            let s = self.triplet_list.pop_back().unwrap();
            Self::rem_triplet_info(&mut self.r_w, &mut self.c_w, s);
            if self.c_w[s as usize] == 0 {
                self.num_diff -= 1;
            }

            if self.l == self.start {
                self.l += 1;
                Self::rem_triplet_info(&mut self.r_v, &mut self.c_v, s);
            }

            self.start += 1;
        }

        self.triplet_list.push_front(t);
        if self.c_w[t as usize] == 0 {
            self.num_diff += 1;
        }
        Self::add_triplet_info(&mut self.r_w, &mut self.c_w, t);
        Self::add_triplet_info(&mut self.r_v, &mut self.c_v, t);

        if self.c_v[t as usize] > self.low_k {
            let mut off = self.triplet_list.len() - (self.l - self.start) - 1;
            loop {
                let triplet = self.triplet_list[off];
                Self::rem_triplet_info(&mut self.r_v, &mut self.c_v, triplet);
                self.l += 1;
                if triplet == t {
                    break;
                }
                off -= 1;
            }
        }

        self.stop += 1;

        if self.triplet_list.len() >= self.max_size && self.num_diff <= 1 {
            self.perfects.clear();
            self.perfects.push(PerfectInterval {
                start: self.start,
                end: self.stop + 1,
                score: 0,
                len: 0,
            });
            false
        } else {
            true
        }
    }

    fn find_perfect(&mut self) {
        let mut counts = self.c_v;
        let mut count = self.stop - self.l;
        let mut score = self.r_v;
        let mut perfect_idx = 0usize;
        let mut max_perfect_score = 0u32;
        let mut max_len = 0usize;
        let mut pos = self.l as isize - 1;

        for triplet in self.triplet_list.iter().skip(count) {
            let cnt = counts[*triplet as usize];
            Self::add_triplet_info(&mut score, &mut counts, *triplet);

            if cnt > 0 && score * 10 > self.thresholds[count] {
                while perfect_idx < self.perfects.len()
                    && pos >= 0
                    && (pos as usize) <= self.perfects[perfect_idx].start
                {
                    let perfect = self.perfects[perfect_idx];
                    if max_perfect_score == 0
                        || max_len * perfect.score as usize
                            > max_perfect_score as usize * perfect.len
                    {
                        max_perfect_score = perfect.score;
                        max_len = perfect.len;
                    }
                    perfect_idx += 1;
                }

                if max_perfect_score == 0
                    || score as usize * max_len >= max_perfect_score as usize * count
                {
                    max_perfect_score = score;
                    max_len = count;
                    if pos >= 0 {
                        self.perfects.insert(
                            perfect_idx,
                            PerfectInterval {
                                start: pos as usize,
                                end: self.stop + 1,
                                score: max_perfect_score,
                                len: count,
                            },
                        );
                    }
                }
            }

            count += 1;
            pos -= 1;
        }
    }
}

fn blastna_to_ncbi2na(base: u8) -> u8 {
    match base {
        b'C' | b'c' | 1 => 1,
        b'G' | b'g' | 2 => 2,
        b'T' | b't' | b'U' | b'u' | 3 => 3,
        b'N' | b'n' => 0,
        _ if base <= 15 => base & 0x3,
        _ => 0,
    }
}

fn save_masked_regions(
    mask: &mut MaskLoc,
    perfects: &mut Vec<PerfectInterval>,
    wstart: usize,
    start: usize,
    linker: usize,
) {
    if let Some(bounds) = perfects.last().copied() {
        if bounds.start < wstart {
            let start_pos = (bounds.start + start) as i32;
            let end_pos = (bounds.end + start) as i32;
            if let Some(last) = mask.regions.last_mut() {
                if last.end as usize + linker >= start_pos as usize {
                    last.end = last.end.max(end_pos);
                } else {
                    mask.add(start_pos, end_pos);
                }
            } else {
                mask.add(start_pos, end_pos);
            }

            while perfects.last().is_some_and(|p| p.start < wstart) {
                perfects.pop();
            }
        }
    }
}

/// Faithful port of NCBI's symmetric DUST interval finder (`symdust.cpp`).
///
/// `level`, `window`, and `linker` correspond to BLAST's real DUST settings.
/// The returned intervals are half-open `[start, end)` regions over the input.
pub fn dust_filter(sequence: &[u8], level: u32, window: usize, linker: usize) -> MaskLoc {
    let level = if (2..=64).contains(&level) {
        level
    } else {
        DUST_DEFAULT_LEVEL
    };
    let window = if (8..=64).contains(&window) {
        window
    } else {
        DUST_DEFAULT_WINDOW
    };
    let linker = if (1..=32).contains(&linker) {
        linker
    } else {
        DUST_DEFAULT_LINKER
    };
    let low_k = (level / 5) as u8;
    let mut thresholds = Vec::with_capacity(window.saturating_sub(2));
    thresholds.push(1);
    for i in 1..window.saturating_sub(2) {
        thresholds.push(i as u32 * level);
    }

    let mut mask = MaskLoc::new();
    if sequence.len() < 3 {
        return mask;
    }

    let mut start = 0usize;
    let stop = sequence.len() - 1;
    while stop > start + 2 {
        let mut triplets = DustTriplets::new(window, low_k, thresholds.clone());
        let mut t =
            (blastna_to_ncbi2na(sequence[start]) << 2) + blastna_to_ncbi2na(sequence[start + 1]);
        let mut next_pos = start + triplets.stop + 2;
        let mut done = false;

        while !done && next_pos <= stop {
            save_masked_regions(
                &mut mask,
                &mut triplets.perfects,
                triplets.start,
                start,
                linker,
            );

            t = ((t << 2) & TRIPLET_MASK) + (blastna_to_ncbi2na(sequence[next_pos]) & 0x3);
            next_pos += 1;

            if triplets.shift_window(t) {
                if triplets.needs_processing() {
                    triplets.find_perfect();
                }
            } else {
                while next_pos <= stop {
                    save_masked_regions(
                        &mut mask,
                        &mut triplets.perfects,
                        triplets.start,
                        start,
                        linker,
                    );
                    t = ((t << 2) & TRIPLET_MASK) + (blastna_to_ncbi2na(sequence[next_pos]) & 0x3);
                    if triplets.shift_window(t) {
                        done = true;
                        break;
                    }
                    next_pos += 1;
                }
            }
        }

        let mut wstart = triplets.start;
        while !triplets.perfects.is_empty() {
            save_masked_regions(&mut mask, &mut triplets.perfects, wstart, start, linker);
            wstart += 1;
        }

        if triplets.start > 0 {
            start += triplets.start;
        } else {
            break;
        }
    }

    mask
}

#[cfg_attr(not(test), allow(dead_code))]
/// Merge overlapping masked regions.
fn merge_regions(mask: &mut MaskLoc) {
    if mask.regions.len() <= 1 {
        return;
    }
    mask.regions.sort_by_key(|r| r.start);
    let mut merged = Vec::new();
    let mut current = mask.regions[0];
    for r in &mask.regions[1..] {
        if r.start <= current.end {
            current.end = current.end.max(r.end);
        } else {
            merged.push(current);
            current = *r;
        }
    }
    merged.push(current);
    mask.regions = merged;
}

const SEG_DEFAULT_WINDOW: usize = 12;
const SEG_DEFAULT_HICUT: f64 = 2.5;
const SEG_DEFAULT_MAXTRIM: usize = 50;
const SEG_DEFAULT_MAXBOGUS: usize = 2;
const SEG_LN_20: f64 = 2.995_732_273_553_991;
const SEG_VALID_AA_CODES: [u8; 20] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];
const LOG_WIN10: [f64; 11] = [
    0.0,
    -2.302_585_09,
    -1.609_437_91,
    -1.203_982_804,
    -0.916_290_73,
    -0.693_147_8,
    -0.510_825_623,
    -0.356_674_944,
    -0.223_143_55,
    -0.105_360_515,
    0.0,
];

#[derive(Clone, Copy)]
struct SegParameters {
    window: usize,
    locut: f64,
    hicut: f64,
    maxtrim: usize,
    maxbogus: usize,
    overlaps: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SegSegment {
    begin: usize,
    end: usize,
}

fn normalize_seg_params(window: usize, locut: f64, hicut: f64) -> SegParameters {
    let window = if window == 0 {
        SEG_DEFAULT_WINDOW
    } else {
        window
    };
    let locut = locut.max(0.0);
    let hicut = hicut.max(locut);
    let maxbogus = SEG_DEFAULT_MAXBOGUS.min(window);
    SegParameters {
        window,
        locut,
        hicut,
        maxtrim: SEG_DEFAULT_MAXTRIM,
        maxbogus,
        overlaps: false,
    }
}

fn seg_ascii_to_ncbistdaa(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| crate::input::aminoacid_to_ncbistdaa(b))
        .collect()
}

fn seg_alpha_index(code: u8) -> Option<usize> {
    SEG_VALID_AA_CODES.iter().position(|&valid| valid == code)
}

fn seg_state_and_bogus(seq: &[u8]) -> (Vec<i32>, usize) {
    let mut composition = [0i32; 20];
    let mut bogus = 0usize;
    for &aa in seq {
        if let Some(idx) = seg_alpha_index(aa) {
            composition[idx] += 1;
        } else {
            bogus += 1;
        }
    }
    let mut state: Vec<i32> = composition.into_iter().filter(|&count| count > 0).collect();
    state.sort_unstable_by(|a, b| b.cmp(a));
    (state, bogus)
}

fn seg_entropy(state: &[i32]) -> f64 {
    let total: i32 = state.iter().copied().sum();
    if total == 0 {
        return 0.0;
    }

    let mut ent = 0.0;
    if total == 10 {
        for &count in state {
            ent += count as f64 * LOG_WIN10[count as usize] / crate::math::NCBIMATH_LN2;
        }
    } else {
        let total_f = total as f64;
        for &count in state {
            let count_f = count as f64;
            ent += count_f * (count_f / total_f).ln() / crate::math::NCBIMATH_LN2;
        }
    }
    (ent / total as f64).abs()
}

fn seg_window_entropy(seq: &[u8]) -> (f64, Vec<i32>, usize) {
    let (state, bogus) = seg_state_and_bogus(seq);
    (seg_entropy(&state), state, bogus)
}

fn seg_seq_entropy(seq: &[u8], window: usize, maxbogus: usize) -> Option<Vec<f64>> {
    if window > seq.len() {
        return None;
    }

    let downset = (window + 1) / 2 - 1;
    let upset = window - downset;
    let first = downset;
    let last = seq.len() - upset;
    let mut values = vec![-1.0; seq.len()];

    for center in first..=last {
        let start = center - downset;
        let end = start + window;
        let (entropy, _, bogus) = seg_window_entropy(&seq[start..end]);
        if bogus <= maxbogus {
            values[center] = entropy;
        }
    }

    Some(values)
}

fn seg_find_low(i: usize, limit: usize, hicut: f64, entropies: &[f64]) -> usize {
    let mut j = i;
    loop {
        if entropies[j] == -1.0 || entropies[j] > hicut {
            return j + 1;
        }
        if j == limit {
            return limit;
        }
        j -= 1;
    }
}

fn seg_find_high(i: usize, limit: usize, hicut: f64, entropies: &[f64]) -> usize {
    let mut j = i;
    while j <= limit {
        if entropies[j] == -1.0 || entropies[j] > hicut {
            return j.saturating_sub(1);
        }
        j += 1;
    }
    limit
}

fn seg_ln_perm(state: &[i32], window_length: usize) -> f64 {
    let mut ans = crate::math::ln_factorial(window_length as i32);
    for &count in state {
        ans -= crate::math::ln_factorial(count);
    }
    ans
}

fn seg_ln_ass(state: &[i32], alphasize: usize) -> f64 {
    let mut ans = crate::math::ln_factorial(alphasize as i32);
    if state.is_empty() {
        return ans;
    }

    let mut total = alphasize as i32;
    let mut class = 1i32;
    let mut prev = state[0];
    let mut i = 1usize;

    loop {
        if i == alphasize {
            ans -= crate::math::ln_factorial(class);
            break;
        } else if i < state.len() && state[i] == prev {
            class += 1;
            i += 1;
            continue;
        } else {
            total -= class;
            ans -= crate::math::ln_factorial(class);
            if i >= state.len() {
                ans -= crate::math::ln_factorial(total);
                break;
            }
            class = 1;
            prev = state[i];
            i += 1;
        }
    }

    ans
}

fn seg_get_prob(state: &[i32], total: usize) -> f64 {
    seg_ln_ass(state, 20) + seg_ln_perm(state, total) - total as f64 * SEG_LN_20
}

fn seg_trim(seq: &[u8], leftend: &mut usize, rightend: &mut usize, params: SegParameters) {
    let mut lend = 0usize;
    let mut rend = seq.len().saturating_sub(1);
    let minlen = (seq.len().saturating_sub(params.maxtrim)).max(1);
    let mut minprob = 1.0f64;

    for len in (minlen + 1..=seq.len()).rev() {
        for start in 0..=seq.len() - len {
            let end = start + len;
            let (state, _) = seg_state_and_bogus(&seq[start..end]);
            let prob = seg_get_prob(&state, len);
            if prob < minprob {
                minprob = prob;
                lend = start;
                rend = end - 1;
            }
        }
    }

    *leftend += lend;
    *rightend -= seq.len() - rend - 1;
}

fn seg_seq(seq: &[u8], params: SegParameters, segs: &mut Vec<SegSegment>, offset: usize) {
    if params.window == 0 || seq.is_empty() || params.window > seq.len() {
        return;
    }

    let Some(entropies) = seg_seq_entropy(seq, params.window, params.maxbogus) else {
        return;
    };

    let downset = (params.window + 1) / 2 - 1;
    let upset = params.window - downset;
    let first = downset;
    let last = seq.len() - upset;
    let mut lowlim = first;
    let mut i = first;

    while i <= last {
        if entropies[i] != -1.0 && entropies[i] <= params.locut {
            let loi = seg_find_low(i, lowlim, params.hicut, &entropies);
            let hii = seg_find_high(i, last, params.hicut, &entropies);

            let mut leftend = loi - downset;
            let mut rightend = hii + upset - 1;
            seg_trim(
                &seq[leftend..=rightend],
                &mut leftend,
                &mut rightend,
                params,
            );

            if i + upset - 1 < leftend {
                let lend = loi - downset;
                let rend = leftend - 1;
                seg_seq(&seq[lend..=rend], params, segs, offset + lend);
            }

            segs.push(SegSegment {
                begin: leftend + offset,
                end: rightend + offset,
            });

            i = hii.max(rightend + downset);
            lowlim = i + 1;
        }
        i += 1;
    }
}

fn seg_merge_segments(seq_len: usize, segs: &mut Vec<SegSegment>) {
    if segs.is_empty() {
        return;
    }

    segs.sort_unstable_by_key(|seg| seg.begin);
    let hilenmin = 0usize;

    if seq_len - 1 - segs[0].end < hilenmin {
        segs[0].end = seq_len - 1;
    }

    let mut merged = Vec::with_capacity(segs.len());
    let mut current = segs[0];

    for &next in segs.iter().skip(1) {
        if current.begin <= next.end + hilenmin + 1 {
            current.begin = current.begin.min(next.begin);
            current.end = current.end.max(next.end);
        } else {
            merged.push(current);
            current = next;
        }
    }

    if current.begin < hilenmin {
        current.begin = 0;
    }
    merged.push(current);
    *segs = merged;
}

/// Faithful port of NCBI's SEG low-complexity masking over NCBIstdaa input.
///
/// The returned intervals are half-open `[start, end)` regions.
pub fn seg_filter_ncbistdaa(sequence: &[u8], window: usize, locut: f64, hicut: f64) -> MaskLoc {
    let params = normalize_seg_params(window, locut, hicut);
    let mut mask = MaskLoc::new();
    if sequence.len() < params.window {
        return mask;
    }

    let mut segs = Vec::new();
    seg_seq(sequence, params, &mut segs, 0);
    if params.overlaps {
        seg_merge_segments(sequence.len(), &mut segs);
    }

    for seg in segs {
        if seg.begin <= seg.end {
            mask.add(seg.begin as i32, seg.end as i32 + 1);
        }
    }

    mask
}

/// SEG filtering over ASCII amino-acid input.
///
/// The wrapper preserves the crate's historical ASCII-oriented API while
/// delegating to the faithful NCBIstdaa implementation with the standard
/// NCBI `hicut` default.
pub fn seg_filter(sequence: &[u8], window: usize, locut: f64) -> MaskLoc {
    seg_filter_ncbistdaa(
        &seg_ascii_to_ncbistdaa(sequence),
        window,
        locut,
        SEG_DEFAULT_HICUT.max(locut),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mask_loc() {
        let mut mask = MaskLoc::new();
        mask.add(5, 10);
        assert!(!mask.is_masked(4));
        assert!(mask.is_masked(5));
        assert!(mask.is_masked(9));
        assert!(!mask.is_masked(10));
    }

    #[test]
    fn test_apply_mask() {
        let mut seq = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let mut mask = MaskLoc::new();
        mask.add(2, 5);
        mask.apply(&mut seq, 14); // 14 = N in BLASTNA
        assert_eq!(seq, vec![0, 1, 14, 14, 14, 1, 2, 3]);
    }

    #[test]
    fn test_merge_regions() {
        let mut mask = MaskLoc::new();
        mask.add(0, 5);
        mask.add(3, 8);
        mask.add(10, 15);
        merge_regions(&mut mask);
        assert_eq!(mask.regions.len(), 2);
        assert_eq!(mask.regions[0].start, 0);
        assert_eq!(mask.regions[0].end, 8);
        assert_eq!(mask.regions[1].start, 10);
    }

    /// Low-complexity poly-A sequence should be flagged by DUST.
    #[test]
    fn test_dust_low_complexity_sequence() {
        // 64-base poly-A — maximally repetitive
        let seq: Vec<u8> = vec![0u8; 64]; // all A (encoded as 0)
        let mask = dust_filter(&seq, 20, 64, 1);
        assert!(
            !mask.regions.is_empty(),
            "Poly-A sequence should produce at least one masked region"
        );
        assert_eq!(mask.regions[0].start, 0);
        assert_eq!(mask.regions[0].end, 63);
    }

    /// A de Bruijn sequence over A/C/G/T contains each 3-mer exactly once, so DUST
    /// should leave it unmasked.
    #[test]
    fn test_dust_normal_sequence() {
        fn debruijn(k: usize, n: usize) -> Vec<usize> {
            fn db(t: usize, p: usize, k: usize, n: usize, a: &mut [usize], out: &mut Vec<usize>) {
                if t > n {
                    if n.is_multiple_of(p) {
                        out.extend_from_slice(&a[1..=p]);
                    }
                } else {
                    a[t] = a[t - p];
                    db(t + 1, p, k, n, a, out);
                    for j in (a[t - p] + 1)..k {
                        a[t] = j;
                        db(t + 1, t, k, n, a, out);
                    }
                }
            }

            let mut a = vec![0; k * n + 1];
            let mut out = Vec::new();
            db(1, 1, k, n, &mut a, &mut out);
            out
        }

        let seq: Vec<u8> = debruijn(4, 3).into_iter().map(|x| x as u8).collect();
        let mask = dust_filter(&seq, 20, 64, 1);
        assert!(
            mask.regions.is_empty(),
            "High-complexity sequence should produce no masked regions"
        );
    }

    /// Empty and very short sequences must not crash.
    #[test]
    fn test_dust_empty_sequence() {
        let empty: Vec<u8> = vec![];
        let mask = dust_filter(&empty, 20, 64, 1);
        assert!(mask.regions.is_empty());

        let short = vec![0u8; 3];
        let mask = dust_filter(&short, 20, 64, 1);
        assert!(mask.regions.is_empty());

        let single = vec![1u8];
        let mask = dust_filter(&single, 20, 64, 1);
        assert!(mask.regions.is_empty());
    }

    /// Verify specific masked positions for a known input.
    #[test]
    fn test_dust_mask_positions() {
        // Build: 30 bases of poly-A, then 60 bases of high-complexity sequence.
        let window = 16;
        let level = 20;
        let linker = 1;
        let mut seq = vec![0u8; 30]; // 30× A — maximally repetitive
        for i in 0..60 {
            seq.push(((i * 7 + 3) % 4) as u8);
        }

        let mask = dust_filter(&seq, level, window, linker);
        // The poly-A region should be masked
        assert!(
            mask.is_masked(0),
            "Position 0 (inside poly-A) should be masked"
        );
        assert!(
            mask.is_masked(15),
            "Position 15 (middle of poly-A) should be masked"
        );
        assert!(
            !mask.is_masked(75),
            "Position 75 (well inside high-complexity region) should not be masked"
        );
    }

    /// Real DUST parameters should affect masking the way BLAST exposes them.
    #[test]
    fn test_dust_with_different_parameters() {
        let mut seq: Vec<u8> = (0..40).map(|i| if i % 2 == 0 { 0 } else { 1 }).collect();
        seq.extend([
            0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 1, 1, 0, 1, 2, 0, 1, 3, 0, 2, 1, 0, 2, 2, 0, 2, 3, 0,
            3, 1, 0, 3, 2, 0, 3, 3, 1, 1, 1, 2, 1, 1, 3, 1, 2, 2, 1, 2, 3, 1, 3, 2, 1, 3, 3, 2, 2,
            2, 3, 2, 2, 3,
        ]);

        let mask_default = dust_filter(&seq, 20, 64, 1);
        let mask_small_window = dust_filter(&seq, 20, 8, 1);
        assert!(
            !mask_default.regions.is_empty(),
            "Default DUST settings should mask the low-complexity prefix"
        );
        assert!(
            mask_small_window.regions.is_empty(),
            "A much smaller DUST window should leave this mixed sequence unmasked"
        );
    }

    #[test]
    fn test_dust_accepts_ascii_nucleotides() {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        let mask = dust_filter(&seq, 20, 32, 1);
        assert_eq!(mask.regions.len(), 1);
        assert_eq!(mask.regions[0].start, 0);
        assert_eq!(mask.regions[0].end, 31);
    }

    #[test]
    fn test_seg_low_complexity_sequence() {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        let mask = seg_filter(&seq, 12, 2.2);
        assert!(
            !mask.regions.is_empty(),
            "poly-A protein sequence should trigger SEG masking"
        );
        assert_eq!(mask.regions[0].start, 0);
        assert!(mask.regions[0].end >= 12);
    }

    #[test]
    fn test_seg_high_complexity_sequence() {
        let seq = b"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY".to_vec();
        let mask = seg_filter(&seq, 12, 2.2);
        assert!(
            mask.regions.is_empty(),
            "diverse protein sequence should not be masked by default SEG"
        );
    }

    #[test]
    fn test_seg_ascii_and_ncbistdaa_match() {
        let ascii = b"AAAAAAAAAAAAQQQQQQQQQQQQ".to_vec();
        let ncbi: Vec<u8> = ascii
            .iter()
            .map(|&b| crate::input::aminoacid_to_ncbistdaa(b))
            .collect();
        let ascii_mask = seg_filter(&ascii, 12, 2.2);
        let ncbi_mask = seg_filter_ncbistdaa(&ncbi, 12, 2.2, 2.5);
        assert_eq!(ascii_mask.regions.len(), ncbi_mask.regions.len());
        assert_eq!(
            ascii_mask
                .regions
                .iter()
                .map(|r| (r.start, r.end))
                .collect::<Vec<_>>(),
            ncbi_mask
                .regions
                .iter()
                .map(|r| (r.start, r.end))
                .collect::<Vec<_>>()
        );
    }
}
