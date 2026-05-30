//! Flat query-anchored alignment view (`-outfmt 4`,
//! `eFlatQueryAnchoredNoIdentities`).
//!
//! NCBI's `CDisplaySeqalign` with the `eMergeAlign` flag merges every HSP for a
//! query into a single multiple-sequence alignment anchored on the query
//! coordinate axis (`CAlnMix::Merge` with `fMinGap | fQuerySeqMergeOnly |
//! fFillUnalignedRegions`). The query is row 0; each HSP becomes one subject
//! row projected onto the shared query axis. Insertions relative to the query
//! (query gap columns) from different HSPs are overlaid into the same columns,
//! their count at each junction being the maximum over all HSPs (`fMinGap`).
//!
//! Each block of `line_len` (default 60) merged columns is printed as:
//! `<id padded to maxIdLen+2><start padded to maxStartLen+2><sequence>  <end>`
//! Rows whose aligned range does not intersect a block are omitted from it; a
//! row that contributes no advancing residue in a block (all blank/gap) has its
//! start and end coordinates suppressed (printed blank). This mirrors
//! `showalign.cpp:x_DisplayRowDataSet` / `x_DisplaySequenceLine`.

use std::io::Write;

/// One HSP's alignment data, in the query's plus orientation.
pub struct FlatHsp<'a> {
    /// Subject row label (accession / Seq-id string, e.g. `BP722512.1`).
    pub subject_id: &'a str,
    /// 1-based query start (<= query_end; query is the plus-strand anchor).
    pub query_start: i32,
    /// 1-based subject start in alignment direction.
    pub subject_start: i32,
    /// 1-based subject end in alignment direction (< subject_start on minus).
    pub subject_end: i32,
    /// Aligned query string (uppercase letters with `-` for gaps).
    pub qseq: &'a str,
    /// Aligned subject string (uppercase letters with `-` for gaps).
    pub sseq: &'a str,
}

// A merged-axis column is either a query match column (carrying a 0-based query
// position) or an insert column sitting between query positions.
#[derive(Clone, Copy)]
enum Column {
    Match(usize), // 0-based query position
    Insert,       // an inserted (query-gap) column between query positions
}

/// Render the flat query-anchored block for a single query.
///
/// `query_index` is the 1-based query ordinal (`Query_1`, `Query_2`, ...).
/// `query_len` is the ungapped query length (number of query residues).
/// `query_seq`, when supplied, holds the full query residue letters (uppercase
/// ASCII, 0-based by query position). NCBI's merge view prints the actual query
/// residues across the entire displayed span — including columns no subject HSP
/// covers — so the anchor row is filled from `query_seq` when available. When
/// `None`, the anchor row is reconstructed from HSP `qseq` fragments only (which
/// is sufficient for a single contiguous HSP, e.g. blastn).
/// Returns without writing anything if there are no HSPs.
pub fn format_flat_query_anchored<W: Write>(
    writer: &mut W,
    query_index: usize,
    query_len: usize,
    query_seq: Option<&[u8]>,
    hsps: &[FlatHsp<'_>],
    line_len: usize,
) -> std::io::Result<()> {
    if hsps.is_empty() || query_len == 0 {
        return Ok(());
    }
    let line_len = line_len.max(1);

    // --- 1. For each HSP, find the maximum run of inserted (query-gap) subject
    // chars at each query junction. junction j sits *before* query position j
    // (0-based); j in 0..=query_len. The merged column layout reserves, at each
    // junction, the maximum insertion length seen over all HSPs (`fMinGap`).
    let mut max_insert = vec![0usize; query_len + 1];

    for hsp in hsps {
        let qbytes: Vec<char> = hsp.qseq.chars().collect();
        let sbytes: Vec<char> = hsp.sseq.chars().collect();
        let n = qbytes.len().min(sbytes.len());
        // 0-based query cursor (advances on non-gap query chars).
        let mut qpos = (hsp.query_start - 1) as usize;
        let mut pending_insert = 0usize;
        for i in 0..n {
            if qbytes[i] == '-' {
                // Insertion relative to query: subject char with no query base.
                // It belongs to the junction before the current qpos.
                pending_insert += 1;
            } else {
                // Match/mismatch column at query position `qpos`.
                if pending_insert != 0 {
                    max_insert[qpos] = max_insert[qpos].max(pending_insert);
                    pending_insert = 0;
                }
                qpos += 1;
            }
        }
        // Trailing insertions (subject extends past query end of this HSP)
        // attach to the junction after the last matched query pos.
        if pending_insert != 0 {
            max_insert[qpos] = max_insert[qpos].max(pending_insert);
        }
    }

    // --- 2. Build the global merged column layout.
    let mut columns: Vec<Column> = Vec::new();
    for q in 0..query_len {
        for _ in 0..max_insert[q] {
            columns.push(Column::Insert);
        }
        columns.push(Column::Match(q));
    }
    for _ in 0..max_insert[query_len] {
        columns.push(Column::Insert);
    }

    // Map: (kind) -> column index, for fast row projection.
    // match_col[q] = column index of query pos q.
    let mut match_col = vec![0usize; query_len];
    // insert_base[j] = first column index of junction j's insert run.
    let mut insert_base = vec![0usize; query_len + 1];
    {
        let mut cursor = 0usize;
        for q in 0..query_len {
            insert_base[q] = cursor;
            cursor += max_insert[q];
            match_col[q] = cursor;
            cursor += 1;
        }
        insert_base[query_len] = cursor;
    }
    let total_cols = columns.len();

    // --- 3. Build the query (anchor) row and each subject row over total_cols.
    // ' ' = blank (no coverage), '-' = gap, letter = residue.
    // The query anchor is contiguous across the displayed span: every match
    // column carries a query letter. When the full query residues are supplied
    // (`query_seq`) they are used directly (matching NCBI even where no subject
    // HSP covers a column); otherwise the query letters are reconstructed from
    // the HSP `qseq` fragments. Insert columns in the query row are gaps.
    let mut query_row = vec![b' '; total_cols];
    if let Some(qseq) = query_seq {
        // Full query residues available: fill every match column from the actual
        // query sequence, matching NCBI's anchor row across the whole span.
        for q in 0..query_len {
            if let Some(&c) = qseq.get(q) {
                query_row[match_col[q]] = c;
            }
        }
    } else {
        // Reconstruct query residues from HSP qseq (first HSP covering each pos).
        for hsp in hsps.iter() {
            let qbytes: Vec<char> = hsp.qseq.chars().collect();
            let mut qpos = (hsp.query_start - 1) as usize;
            for &qc in &qbytes {
                if qc != '-' {
                    if qpos < query_len {
                        let col = match_col[qpos];
                        if query_row[col] == b' ' {
                            query_row[col] = qc as u8;
                        }
                    }
                    qpos += 1;
                }
            }
        }
    }
    // Query row: insert columns are gaps within the query's covered span.
    let q_lo = hsps
        .iter()
        .map(|h| (h.query_start - 1) as usize)
        .min()
        .unwrap_or(0);
    let q_hi = hsps
        .iter()
        .map(|h| query_end_pos(h))
        .max()
        .unwrap_or(0);
    for q in q_lo..=q_hi.min(query_len.saturating_sub(1)) {
        for c in insert_base[q]..(insert_base[q] + max_insert[q]) {
            query_row[c] = b'-';
        }
    }
    // trailing junction inserts inside span
    if q_hi + 1 <= query_len {
        for c in insert_base[q_hi + 1]..(insert_base[q_hi + 1] + max_insert[q_hi + 1]) {
            query_row[c] = b'-';
        }
    }

    // Subject rows.
    struct Row {
        id: String,
        seq: Vec<u8>,       // total_cols, ' ' for blank
        first_col: usize,   // first non-blank column
        last_col: usize,    // last non-blank column
        // For each merged column index, the 1-based subject coordinate of the
        // residue placed there (only meaningful where a residue/gap-with-coord
        // exists). We track coords via a parallel walk instead.
        // Per-column subject coordinate after consuming this column's residue.
        coord_at: Vec<i32>, // subject coordinate of the residue in this column (0 if blank/gap)
        dir: i32,
    }

    let mut rows: Vec<Row> = Vec::with_capacity(hsps.len());
    for hsp in hsps.iter() {
        let dir: i32 = if hsp.subject_end >= hsp.subject_start { 1 } else { -1 };
        let mut seq = vec![b' '; total_cols];
        let mut coord_at = vec![0i32; total_cols];

        // Walk the HSP in alignment order to assign subject coordinates and
        // place residues into merged columns. We must place inserts
        // left-justified into the shared insert run for their junction.
        let qbytes: Vec<char> = hsp.qseq.chars().collect();
        let sbytes: Vec<char> = hsp.sseq.chars().collect();
        let nn = qbytes.len().min(sbytes.len());
        let mut qpos = (hsp.query_start - 1) as usize;
        let mut scoord = hsp.subject_start;
        let mut insert_slot: usize = 0;
        let mut cur_junction: usize = qpos;
        for i in 0..nn {
            let qc = qbytes[i];
            let sc = sbytes[i];
            if qc == '-' {
                // insertion at junction before qpos
                if cur_junction != qpos {
                    cur_junction = qpos;
                    insert_slot = 0;
                }
                let col = insert_base[qpos] + insert_slot;
                insert_slot += 1;
                if sc == '-' {
                    seq[col] = b'-';
                } else {
                    seq[col] = sc as u8;
                    coord_at[col] = scoord;
                    scoord += dir;
                }
            } else {
                let col = match_col[qpos];
                if sc == '-' {
                    seq[col] = b'-';
                } else {
                    seq[col] = sc as u8;
                    coord_at[col] = scoord;
                    scoord += dir;
                }
                qpos += 1;
                insert_slot = 0;
                cur_junction = qpos;
            }
        }

        // Determine covered span (first/last non-blank column). Inserts and
        // matches define the span; blanks remain outside.
        let first_col = seq.iter().position(|&b| b != b' ').unwrap_or(0);
        let last_col = seq.iter().rposition(|&b| b != b' ').unwrap_or(0);
        // Within the row's covered span, any column not carrying a residue is a
        // gap (`-`) — this includes insert columns created by *other* HSPs
        // where this HSP has no insertion. Outside the span, columns stay blank
        // (spaces). Mirrors NCBI's per-row gap fill within `rowRng`.
        if seq.iter().any(|&b| b != b' ') {
            for col in &mut seq[first_col..=last_col] {
                if *col == b' ' {
                    *col = b'-';
                }
            }
        }
        rows.push(Row {
            id: hsp.subject_id.to_string(),
            seq,
            first_col,
            last_col,
            coord_at,
            dir,
        });
    }

    // --- 4. Compute column widths.
    let max_id_len = std::iter::once(format!("Query_{}", query_index).len())
        .chain(rows.iter().map(|r| r.id.len()))
        .max()
        .unwrap_or(0);

    // maxStartLen: width of the widest start coordinate that will be printed.
    // Compute over all blocks for query + subjects. The query start at each
    // block is its 1-based query position; subject start is its coord.
    let span_first = rows
        .iter()
        .filter(|r| r.seq.iter().any(|&b| b != b' '))
        .map(|r| r.first_col)
        .min()
        .unwrap_or(0);
    let span_last = rows
        .iter()
        .filter(|r| r.seq.iter().any(|&b| b != b' '))
        .map(|r| r.last_col)
        .max()
        .unwrap_or(total_cols.saturating_sub(1));

    let mut max_start_len = 1usize;
    // query start coords per block:
    let mut block = 0usize;
    while span_first + block * line_len <= span_last {
        let cstart = span_first + block * line_len;
        let cend = (cstart + line_len).min(span_last + 1);
        // query start in this block = 1-based query position of first match col
        if let Some(qs) = first_query_pos_in(&columns, cstart, cend) {
            max_start_len = max_start_len.max((qs + 1).to_string().len());
        }
        for r in &rows {
            if r.last_col < cstart || r.first_col >= cend {
                continue;
            }
            if let Some(c) = (cstart..cend).find(|&c| r.coord_at[c] != 0) {
                max_start_len = max_start_len.max(r.coord_at[c].unsigned_abs().to_string().len());
            }
        }
        block += 1;
    }

    let id_pad = max_id_len + 2; // k_IdStartMargin
    // start field width = maxStartLen + k_StartSequenceMargin(2)

    // --- 5. Emit blocks.
    // NCBI's merge view only displays the columns actually covered by some
    // subject row; the query (anchor) row is shown over that same span and
    // unaligned query prefixes/suffixes are not printed. The first block begins
    // at the first covered merged column (so block boundaries are anchored on
    // the covered region, not on absolute column 0), and the last block ends at
    // the last covered column. `span_first`/`span_last` were computed above.
    let query_label = format!("Query_{}", query_index);
    // NCBI's merge view prints one leading blank line before the first block;
    // blocks are separated by a single blank line, and there is no trailing
    // blank after the last block (the per-query stats footer supplies the
    // following blank lines).
    writeln!(writer)?;
    let num_blocks = (span_last + 1 - span_first).div_ceil(line_len);
    // query coordinate cursor (1-based) advances by query residues per block.
    let mut block = 0usize;
    while span_first + block * line_len <= span_last {
        let cstart = span_first + block * line_len;
        let cend = (cstart + line_len).min(span_last + 1);

        // Query (anchor) row — always printed.
        let q_start_1b = first_query_pos_in(&columns, cstart, cend).map(|p| p + 1);
        let q_end_1b = last_query_pos_in(&columns, cstart, cend).map(|p| p + 1);
        write_row(
            writer,
            &query_label,
            id_pad,
            max_start_len,
            q_start_1b,
            q_end_1b,
            &query_row[cstart..cend],
        )?;

        for r in &rows {
            if r.last_col < cstart || r.first_col >= cend {
                continue; // row does not intersect this block
            }
            // subject start = coord of first residue-bearing column in block
            let start_coord = (cstart..cend).find(|&c| r.coord_at[c] != 0).map(|c| r.coord_at[c]);
            let end_coord = (cstart..cend).rev().find(|&c| r.coord_at[c] != 0).map(|c| r.coord_at[c]);
            let _ = r.dir;
            write_row(
                writer,
                &r.id,
                id_pad,
                max_start_len,
                start_coord.map(|c| c.unsigned_abs() as usize),
                end_coord.map(|c| c.unsigned_abs() as usize),
                &r.seq[cstart..cend],
            )?;
        }
        if block + 1 < num_blocks {
            writeln!(writer)?;
        }
        block += 1;
    }

    Ok(())
}

fn query_end_pos(h: &FlatHsp<'_>) -> usize {
    // last 0-based query position covered = query_start-1 + (#non-gap qseq)-1
    let nq = h.qseq.chars().filter(|&c| c != '-').count();
    ((h.query_start - 1) as usize) + nq.saturating_sub(1)
}

fn first_query_pos_in(columns: &[Column], cstart: usize, cend: usize) -> Option<usize> {
    columns[cstart..cend].iter().find_map(|c| match c {
        Column::Match(q) => Some(*q),
        _ => None,
    })
}

fn last_query_pos_in(columns: &[Column], cstart: usize, cend: usize) -> Option<usize> {
    columns[cstart..cend].iter().rev().find_map(|c| match c {
        Column::Match(q) => Some(*q),
        _ => None,
    })
}

#[allow(clippy::too_many_arguments)]
fn write_row<W: Write>(
    writer: &mut W,
    id: &str,
    id_pad: usize,
    max_start_len: usize,
    start: Option<usize>,
    end: Option<usize>,
    seq: &[u8],
) -> std::io::Result<()> {
    // id left-justified, padded to id_pad
    write!(writer, "{:<width$}", id, width = id_pad)?;
    // start field: start number then pad to (max_start_len + 2). When start is
    // None (empty/all-gap row in block) the whole field is blanks.
    match start {
        Some(s) => {
            let s_str = s.to_string();
            write!(writer, "{}", s_str)?;
            for _ in 0..(max_start_len - s_str.len() + 2) {
                write!(writer, " ")?;
            }
        }
        None => {
            for _ in 0..(max_start_len + 2) {
                write!(writer, " ")?;
            }
        }
    }
    // sequence (blanks already encoded as spaces)
    writer.write_all(seq)?;
    // 2-space margin then end coordinate (suppressed when start suppressed)
    write!(writer, "  ")?;
    match end {
        Some(e) => writeln!(writer, "{}", e)?,
        None => writeln!(writer)?,
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(hsps: &[FlatHsp<'_>], qlen: usize, line: usize) -> String {
        let mut buf = Vec::new();
        format_flat_query_anchored(&mut buf, 1, qlen, None, hsps, line).unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn self_hit_no_gaps() {
        let hsps = [FlatHsp {
            subject_id: "S1",
            query_start: 1,
            subject_start: 1,
            subject_end: 8,
            qseq: "ACGTACGT",
            sseq: "ACGTACGT",
        }];
        let out = run(&hsps, 8, 60);
        // NCBI's merge view emits a leading blank line before the first block,
        // so the query/subject rows are not at index 0.
        let query_line = out.lines().find(|l| l.starts_with("Query_1")).unwrap();
        assert!(query_line.contains("ACGTACGT"));
        assert!(query_line.trim_end().ends_with("8"));
        let s1 = out.lines().find(|l| l.starts_with("S1")).unwrap();
        assert!(s1.contains("ACGTACGT"));
    }

    #[test]
    fn second_subject_partial_with_query_gap() {
        // S2 starts later and the query has an inserted gap from S1.
        let hsps = [
            FlatHsp {
                subject_id: "S1",
                query_start: 1,
                subject_start: 1,
                subject_end: 7,
                qseq: "ACGT-ACG",
                sseq: "ACGTAACG",
            },
            FlatHsp {
                subject_id: "S2",
                query_start: 5,
                subject_start: 100,
                subject_end: 103,
                qseq: "ACG",
                sseq: "ACG",
            },
        ];
        let out = run(&hsps, 7, 60);
        // query row should contain a '-' from S1's insertion
        let q = out.lines().find(|l| l.starts_with("Query_1")).unwrap();
        assert!(q.contains('-'), "query row should carry merged gap: {q}");
        // S2 row appears, coordinates 100..102
        let s2 = out.lines().find(|l| l.starts_with("S2")).unwrap();
        assert!(s2.contains("100"));
    }

    #[test]
    fn wraps_at_line_len() {
        let q: String = "ACGT".repeat(30); // 120
        let hsps = [FlatHsp {
            subject_id: "S1",
            query_start: 1,
            subject_start: 1,
            subject_end: 120,
            qseq: &q,
            sseq: &q,
        }];
        let out = run(&hsps, 120, 60);
        let query_lines = out.lines().filter(|l| l.starts_with("Query_1")).count();
        assert_eq!(query_lines, 2);
    }
}
