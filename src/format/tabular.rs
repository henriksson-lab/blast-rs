//! BLAST tabular output format (-outfmt 6/7).
//!
//! Default columns: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore

use std::collections::HashMap;
use std::io::Write;

/// Strip a trailing `.N` version suffix from an accession-like string.
/// `AAC44598.1` → `AAC44598`; `gi|123|ref|NP_839091.2|` is left untouched
/// because the trailing `|` after the digits prevents the strip.
pub fn strip_accession_version(s: &str) -> &str {
    if let Some(dot) = s.rfind('.') {
        let suffix = &s[dot + 1..];
        if !suffix.is_empty() && suffix.bytes().all(|b| b.is_ascii_digit()) {
            return &s[..dot];
        }
    }
    s
}

pub const DEFAULT_TABULAR_COLUMNS: &str =
    "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore";

// blast-rs: local parser helper for `-outfmt` token expansion; NCBI handles
// this inside its command-line/output-format parser, not as a standalone helper.
fn expand_column_token<'a>(token: &'a str, columns: &mut Vec<&'a str>) {
    if token == "std" {
        columns.extend(DEFAULT_TABULAR_COLUMNS.split_whitespace());
    } else {
        columns.push(token);
    }
}

fn remove_column_token(token: &str, columns: &mut Vec<&str>) {
    if token == "std" {
        for col in DEFAULT_TABULAR_COLUMNS.split_whitespace() {
            columns.retain(|existing| *existing != col);
        }
    } else {
        columns.retain(|existing| *existing != token);
    }
}

fn apply_column_token<'a>(token: &'a str, columns: &mut Vec<&'a str>) {
    if let Some(remove) = token.strip_prefix('-').filter(|remove| !remove.is_empty()) {
        remove_column_token(remove, columns);
    } else {
        expand_column_token(token, columns);
    }
}

fn is_known_column_token(token: &str) -> bool {
    let token = token.strip_prefix('-').unwrap_or(token);
    token == "std" || field_display_name(token) != "unknown field"
}

fn normalize_column_tokens(cols: &mut Vec<&str>) {
    let mut seen = std::collections::HashSet::new();
    cols.retain(|col| seen.insert(*col));
    if cols.is_empty() {
        cols.extend(DEFAULT_TABULAR_COLUMNS.split_whitespace());
    }
}

pub fn expanded_column_tokens(columns: &str) -> Vec<&str> {
    expanded_column_tokens_checked(columns)
        .unwrap_or_else(|_| DEFAULT_TABULAR_COLUMNS.split_whitespace().collect())
}

fn expanded_column_tokens_checked(columns: &str) -> std::io::Result<Vec<&str>> {
    let mut cols = Vec::new();
    for col in columns.split_whitespace() {
        if !col.starts_with("delim=") && is_known_column_token(col) {
            apply_column_token(col, &mut cols);
        }
    }
    normalize_column_tokens(&mut cols);
    Ok(cols)
}

/// A single alignment hit for tabular output.
#[derive(Clone)]
pub struct TabularHit {
    pub query_id: String,
    pub query_gi: Option<String>,
    pub query_acc: Option<String>,
    pub query_accver: Option<String>,
    pub subject_id: String,
    /// Full BLAST Seq-id chain (e.g. `gi|N|ref|acc.ver|`, `pir||name`) for
    /// `sseqid` / `sallseqid` tabular columns. When `None`, the formatter
    /// falls back to `subject_id`.
    pub subject_seqid: Option<String>,
    pub subject_gi: Option<String>,
    pub subject_acc: Option<String>,
    pub subject_accver: Option<String>,
    pub subject_title: String,
    pub pct_identity: f64,
    pub align_len: i32,
    pub mismatches: i32,
    pub gap_opens: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub subject_start: i32,
    pub subject_end: i32,
    /// True only when the subject coordinates/sequence are a plain nucleotide
    /// sequence, as in blastn-style output. Translated subject programs report
    /// protein-frame alignments and NCBI prints N/A for nucleotide-only fields.
    pub subject_is_plain_nucleotide: bool,
    pub evalue: f64,
    pub bit_score: f64,
    pub query_len: i32,
    pub subject_len: i32,
    pub raw_score: i32,
    pub qseq: Option<String>,
    pub sseq: Option<String>,
    pub qframe: i32,
    pub sframe: i32,
    pub subject_taxids: Vec<i32>,
    pub subject_sci_name: String,
    pub subject_common_name: String,
    pub subject_blast_name: String,
    pub subject_kingdom: String,
    pub num_ident: i32,
    pub num_positives: i32,
    /// Number of HSPs in this HSP's linked-set chain (NCBI `BlastHSP::num`).
    /// 1 for singletons (pairwise prints `Expect = X`); >=2 for sum-stats
    /// chains (pairwise prints `Expect(N) = X`).
    pub num_links: i32,
    /// Composition adjustment method actually applied per HSP, mirroring
    /// `BlastHSP::comp_adjustment_method`. Drives the `, Method: ...` suffix
    /// on the pairwise Score line.
    /// 0 = none, 1 = Composition-based stats., 2 = Compositional matrix adjust.
    pub comp_adjust_method: u8,
}

/// Format an e-value matching BLAST tabular output (-outfmt 6).
/// Uses NStr::DoubleToString(evalue, 2, fDoubleScientific) which is
/// essentially %.2e for small values, transitioning to fixed for larger ones.
/// The fixed-notation threshold uses `< 0.0005` so values like `9.98e-04`
/// (which round to `0.001` at 3 decimal places) render as `0.001` —
/// matching NCBI's tabular output behavior (`%.3f` rounding).
pub fn format_evalue(val: f64) -> String {
    // NCBI's tabular formatter (`align_format/tabular.cpp:1335`) overrides
    // `GetScoreString`'s small-evalue branch with `DoubleToString(evalue, 2,
    // fDoubleScientific)` for `1.0e-180 <= evalue < 0.0009`. We must match
    // that 0.0009 boundary, not 0.0005, so 5.23e-04 renders as "5.23e-04"
    // rather than rounding into the "0.001" fixed-point branch.
    if val < 1.0e-180 {
        "0.0".to_string()
    } else if val < 0.0009 {
        // Scientific notation with 2 decimal places, C-style zero-padded exponent
        let s = format!("{:.2e}", val);
        // Rust: "2.01e-4" → need "2.01e-04"
        // Fix: ensure exponent has at least 2 digits
        if let Some(e_pos) = s.find('e') {
            let (mantissa, exp_part) = s.split_at(e_pos);
            let exp_str = &exp_part[1..]; // skip 'e'
            let (sign, digits) = if let Some(rest) = exp_str.strip_prefix('-') {
                ("-", rest)
            } else if let Some(rest) = exp_str.strip_prefix('+') {
                ("", rest)
            } else {
                ("", exp_str)
            };
            if digits.len() < 2 {
                format!(
                    "{}e{}{:02}",
                    mantissa,
                    sign,
                    digits.parse::<u32>().unwrap_or(0)
                )
            } else {
                s
            }
        } else {
            s
        }
    } else if val < 0.1 {
        // Fixed with 3 decimal places
        format!("{:.3}", val)
    } else if val < 1.0 {
        // Fixed with 2 decimal places
        format!("{:.2}", val)
    } else if val < 10.0 {
        // Fixed with 1 decimal place
        format!("{:.1}", val)
    } else {
        // Integer
        format!("{:.0}", val)
    }
}

/// Format a bit score matching BLAST reference style.
/// `> 99999`: scientific 3 decimal places (NCBI `%5.3le` → `"1.000e+05"`);
/// `> 99.9`: integer (e.g. `"198"`); else: one decimal (e.g. `"56.0"`).
///
/// NCBI `align_format_util.cpp::GetScoreString` (~line 987):
/// ```c
/// if (bit_score > 99999)      snprintf(..., "%5.3le", bit_score);
/// else if (bit_score > 99.9)  snprintf(..., "%3.0ld", (long)bit_score);
/// else                         snprintf(..., "%4.1lf", bit_score);
/// ```
/// Rust's `{:.3e}` produces `"1.000e5"` (no sign, no zero-padded exponent),
/// so we have to manually shape the exponent to match `%5.3le`.
pub fn format_bitscore(val: f64) -> String {
    if val > 99999.0 {
        // NCBI `%5.3le` format: 3 decimal places, explicit `e` with sign and
        // zero-padded 2-digit exponent (e.g. `1.234e+05`).
        let raw = format!("{val:.3e}");
        if let Some(e_idx) = raw.find('e') {
            let (mantissa, exp) = raw.split_at(e_idx);
            let exp_digits = &exp[1..];
            let (sign, digits) = if let Some(rest) = exp_digits.strip_prefix('-') {
                ('-', rest)
            } else {
                ('+', exp_digits.strip_prefix('+').unwrap_or(exp_digits))
            };
            let n: u32 = digits.parse().unwrap_or(0);
            format!("{mantissa}e{sign}{n:02}")
        } else {
            raw
        }
    } else if val > 99.9 {
        format!("{:.0}", val as i64)
    } else {
        format!("{:4.1}", val)
    }
}

/// Get a field value by column name for custom tabular output.
pub fn get_field(hit: &TabularHit, column: &str) -> String {
    get_field_with_qcovs(hit, column, None)
}

pub fn field_display_name(column: &str) -> &'static str {
    match column {
        "qseqid" => "query id",
        "qgi" => "query gi",
        "qacc" => "query acc.",
        "qaccver" => "query acc.ver",
        "sseqid" => "subject id",
        "sgi" => "subject gi",
        "sallgi" => "all subject gis",
        "sacc" => "subject acc.",
        "saccver" => "subject acc.ver",
        "sallseqid" => "all subject ids",
        "sallacc" => "all subject acc.",
        "stitle" => "subject title",
        "salltitles" => "all subject title(s)",
        "pident" => "% identity",
        "length" => "alignment length",
        "mismatch" => "mismatches",
        "gapopen" => "gap opens",
        "qstart" => "q. start",
        "qend" => "q. end",
        "sstart" => "s. start",
        "send" => "s. end",
        "evalue" => "evalue",
        "bitscore" => "bit score",
        "score" => "score",
        "nident" => "identical",
        "positive" => "positives",
        "gaps" => "gaps",
        "ppos" => "% positives",
        "qlen" => "query length",
        "slen" => "subject length",
        "qcovs" => "% query coverage per subject",
        "qcovhsp" => "% query coverage per hsp",
        "qcovus" => "% query coverage per unique subject",
        "staxid" => "subject tax id",
        "staxids" => "unique subject tax ids",
        "ssciname" | "sscinames" => "subject sci name",
        "scomname" | "scomnames" => "subject common name",
        "sblastname" | "sblastnames" => "subject blast name",
        "sskingdom" | "sskingdoms" => "subject super kingdom",
        "qseq" => "aligned part of query sequence",
        "sseq" => "aligned part of subject sequence",
        "btop" => "BTOP",
        "sstrand" => "subject strand",
        "qframe" => "query frame",
        "sframe" => "subject frame",
        "frames" => "query/sbjct frames",
        _ => "unknown field",
    }
}

fn get_field_with_qcovs(hit: &TabularHit, column: &str, qcovs: Option<i32>) -> String {
    match column {
        "qseqid" => hit.query_id.clone(),
        "qacc" => hit.query_acc.as_ref().unwrap_or(&hit.query_id).clone(),
        "qaccver" => hit
            .query_accver
            .as_ref()
            .or(hit.query_acc.as_ref())
            .unwrap_or(&hit.query_id)
            .clone(),
        "qgi" => hit.query_gi.as_deref().unwrap_or("0").to_string(),
        "sseqid" | "sallseqid" => hit
            .subject_seqid
            .as_deref()
            .unwrap_or(hit.subject_id.as_str())
            .to_string(),
        "sacc" | "sallacc" => {
            // NCBI's `-outfmt 6 sacc` returns the bare accession when one
            // is available from the ASN.1 Seq-id with a separate `version`
            // field (Refseq, GenBank, etc.). For DBs whose IDs are parsed
            // from titles (e.g. celegans's `NC_003279.8` — the `.8` is part
            // of the accession itself, not a separable version), NCBI keeps
            // the `.8`. Mirror by stripping only when we have an explicit
            // versioned accession source (`subject_accver`); otherwise pass
            // subject_id through verbatim.
            if let Some(acc) = hit.subject_acc.as_ref() {
                acc.clone()
            } else if let Some(accver) = hit.subject_accver.as_ref() {
                strip_accession_version(accver).to_string()
            } else {
                hit.subject_id.clone()
            }
        }
        "saccver" => hit
            .subject_accver
            .as_ref()
            .unwrap_or(&hit.subject_id)
            .clone(),
        "sgi" | "sallgi" => hit.subject_gi.as_deref().unwrap_or("0").to_string(),
        "stitle" | "salltitles" => {
            if hit.subject_title.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_title.clone()
            }
        }
        "pident" => format!("{:.3}", hit.pct_identity),
        "length" => hit.align_len.to_string(),
        "mismatch" => {
            if let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) {
                // NCBI counts mismatches against the unmasked alignment; for
                // translated programs we may carry SEG soft-mask info as
                // lowercase ASCII in qseq. Case-fold so a lowercase-vs-uppercase
                // pair is recognized as an identity, not a mismatch.
                qseq.bytes()
                    .zip(sseq.bytes())
                    .filter(|&(q, s)| {
                        q != b'-' && s != b'-' && q.to_ascii_uppercase() != s.to_ascii_uppercase()
                    })
                    .count()
                    .to_string()
            } else {
                hit.mismatches.to_string()
            }
        }
        "gapopen" => hit.gap_opens.to_string(),
        "qstart" => hit.query_start.to_string(),
        "qend" => hit.query_end.to_string(),
        "sstart" => hit.subject_start.to_string(),
        "send" => hit.subject_end.to_string(),
        "evalue" => format_evalue(hit.evalue),
        "bitscore" => format_bitscore(hit.bit_score),
        "score" => hit.raw_score.to_string(),
        "nident" => hit.num_ident.to_string(),
        "positive" => hit.num_positives.to_string(),
        "gaps" => {
            if let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) {
                qseq.bytes()
                    .chain(sseq.bytes())
                    .filter(|&base| base == b'-')
                    .count()
                    .to_string()
            } else {
                (hit.align_len - hit.num_ident - hit.mismatches)
                    .max(0)
                    .to_string()
            }
        }
        "ppos" => {
            if hit.align_len > 0 {
                format!(
                    "{:.2}",
                    100.0 * hit.num_positives as f64 / hit.align_len as f64
                )
            } else {
                "0.00".to_string()
            }
        }
        "qlen" => hit.query_len.to_string(),
        "slen" => hit.subject_len.to_string(),
        "qcovs" => {
            if let Some(qcovs) = qcovs {
                qcovs.to_string()
            } else if hit.query_len > 0 {
                let cov = 100.0 * (hit.query_end - hit.query_start + 1).abs() as f64
                    / hit.query_len as f64;
                format_query_coverage(cov)
            } else {
                "0".to_string()
            }
        }
        "qcovus" => {
            // NCBI emits `qcovus` only for blastn-style plain nucleotide
            // subjects. Protein alignments can spell words like ACGT, so this
            // must not be inferred from the rendered sequence alphabet.
            if !hit.subject_is_plain_nucleotide {
                "N/A".to_string()
            } else if let Some(qcovs) = qcovs {
                qcovs.to_string()
            } else if hit.query_len > 0 {
                let cov = 100.0 * (hit.query_end - hit.query_start + 1).abs() as f64
                    / hit.query_len as f64;
                format_query_coverage(cov)
            } else {
                "0".to_string()
            }
        }
        "qcovhsp" => {
            // NCBI `Blast_HSPGetQueryCoverage` (blast_hits.c:1034): the coverage
            // is the QUERY SPAN of the HSP over the query length, NOT align_len:
            //   pct = 100 * (query.end - query.offset) / query_length
            //   if (pct < 99) pct += 0.5;   // round half up, but only below 99
            //   return (int)pct;            // truncate (so 99.x stays 99)
            // align_len includes subject-gap columns and (for translated
            // searches) uses different units than query_length, so it is wrong.
            ncbi_query_coverage(
                (hit.query_end - hit.query_start + 1).unsigned_abs() as i64,
                hit.query_len as i64,
            )
            .to_string()
        }
        // Taxonomy ID fields — from .nto database file
        "staxid" => hit
            .subject_taxids
            .first()
            .map(|t| t.to_string())
            .unwrap_or_else(|| "N/A".to_string()),
        "staxids" => {
            if hit.subject_taxids.is_empty() {
                "N/A".to_string()
            } else {
                let mut unique_taxids = Vec::new();
                for taxid in &hit.subject_taxids {
                    if !unique_taxids.contains(taxid) {
                        unique_taxids.push(*taxid);
                    }
                }
                unique_taxids
                    .iter()
                    .map(|t| t.to_string())
                    .collect::<Vec<_>>()
                    .join(";")
            }
        }
        // Taxonomy name fields — from taxdb.bti/taxdb.btd
        "ssciname" | "sscinames" => {
            if hit.subject_sci_name.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_sci_name.clone()
            }
        }
        "scomname" | "scomnames" => {
            if hit.subject_common_name.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_common_name.clone()
            }
        }
        "sblastname" | "sblastnames" => {
            if hit.subject_blast_name.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_blast_name.clone()
            }
        }
        "sskingdom" | "sskingdoms" => {
            if hit.subject_kingdom.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_kingdom.clone()
            }
        }
        // Sequence hash
        // NCBI's tabular qseq/sseq is always uppercase (the soft-mask
        // lowercase-display is a pairwise-only convention from `showalign.cpp`).
        // Our translated engines may leave soft-mask info in the alignment as
        // lowercase ASCII; force-uppercase here so tabular matches NCBI.
        "qseq" => hit.qseq.as_deref().unwrap_or("N/A").to_ascii_uppercase(),
        "sseq" => hit.sseq.as_deref().unwrap_or("N/A").to_ascii_uppercase(),
        "btop" => format_btop(hit),
        "sstrand" => {
            // Same applicability as qcovus: only plain nucleotide subject
            // coordinates carry plus/minus strand in NCBI tabular output.
            if !hit.subject_is_plain_nucleotide {
                "N/A".to_string()
            } else if hit.subject_start <= hit.subject_end {
                "plus".to_string()
            } else {
                "minus".to_string()
            }
        }
        // Frame fields. NCBI emits `1` for both query and subject in
        // blastp / blastn / psiblast (non-translated). Translated programs
        // (blastx / tblastn / tblastx) emit the actual frame for the
        // translated side and `0` for the un-translated side. We detect
        // "both zero" — only true for non-translated programs — and emit
        // 1/1 in that case; otherwise pass the recorded frame through.
        "qframe" => default_frame_for_nontranslated(hit.qframe, hit.sframe),
        "sframe" => default_frame_for_nontranslated(hit.sframe, hit.qframe),
        "frames" => format!(
            "{}/{}",
            default_frame_for_nontranslated(hit.qframe, hit.sframe),
            default_frame_for_nontranslated(hit.sframe, hit.qframe)
        ),
        _ => "N/A".to_string(),
    }
}

fn default_frame_for_nontranslated(this: i32, other: i32) -> String {
    if this == 0 && other == 0 {
        "1".to_string()
    } else {
        this.to_string()
    }
}

fn format_btop(hit: &TabularHit) -> String {
    let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) else {
        return "N/A".to_string();
    };

    let mut out = String::new();
    let mut matches = 0;
    for (q, s) in qseq.bytes().zip(sseq.bytes()) {
        let q = q.to_ascii_uppercase();
        let s = s.to_ascii_uppercase();
        if q == s && q != b'-' {
            matches += 1;
        } else {
            if matches > 0 {
                out.push_str(&matches.to_string());
                matches = 0;
            }
            out.push(q as char);
            out.push(s as char);
        }
    }
    if matches > 0 {
        out.push_str(&matches.to_string());
    }
    out
}

// blast-rs: native interval merge for qcovs/qcovus aggregation across emitted
// hits; NCBI computes equivalent coverage during formatting state assembly.
fn compute_qcovs_by_query_subject(hits: &[TabularHit]) -> HashMap<(&str, &str), i32> {
    // (query_len, list of (start, end) intervals on this query×subject pair)
    type QCovBuckets<'a> = HashMap<(&'a str, &'a str), (i32, Vec<(i32, i32)>)>;
    let mut intervals: QCovBuckets = HashMap::new();
    for hit in hits {
        let start = hit.query_start.min(hit.query_end);
        let end = hit.query_start.max(hit.query_end);
        intervals
            .entry((hit.query_id.as_str(), hit.subject_id.as_str()))
            .or_insert_with(|| (hit.query_len, Vec::new()))
            .1
            .push((start, end));
    }

    let mut out = HashMap::new();
    for (key, (query_len, mut ranges)) in intervals {
        if query_len <= 0 || ranges.is_empty() {
            out.insert(key, 0);
            continue;
        }
        ranges.sort_unstable_by_key(|&(start, end)| (start, end));
        let mut covered = 0;
        let mut cur = ranges[0];
        for &(start, end) in &ranges[1..] {
            if start <= cur.1 + 1 {
                cur.1 = cur.1.max(end);
            } else {
                covered += cur.1 - cur.0 + 1;
                cur = (start, end);
            }
        }
        covered += cur.1 - cur.0 + 1;
        let cov = rounded_query_coverage(100.0 * covered as f64 / query_len as f64);
        out.insert(key, cov);
    }
    out
}

/// Write tabular output with custom columns.
/// `columns` is a space-separated list of field names.
pub fn format_tabular_custom<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    columns: &str,
) -> std::io::Result<()> {
    format_tabular_custom_with_delimiter(writer, hits, columns, "\t")
}

pub fn format_tabular_custom_with_delimiter<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    columns: &str,
    default_delimiter: &str,
) -> std::io::Result<()> {
    let mut delimiter = default_delimiter;
    let mut delimiter_seen = false;
    let mut cols = Vec::new();
    for col in columns.split_whitespace() {
        if let Some(value) = col.strip_prefix("delim=") {
            if !delimiter_seen {
                delimiter = match value {
                    "" => default_delimiter,
                    value => value,
                };
                delimiter_seen = true;
            }
        } else if is_known_column_token(col) {
            apply_column_token(col, &mut cols);
        }
    }
    normalize_column_tokens(&mut cols);
    let qcovs_by_subject = if cols.iter().any(|&c| c == "qcovs" || c == "qcovus") {
        compute_qcovs_by_query_subject(hits)
    } else {
        HashMap::new()
    };
    for hit in hits {
        let qcovs = qcovs_by_subject
            .get(&(hit.query_id.as_str(), hit.subject_id.as_str()))
            .copied();
        let fields: Vec<String> = cols
            .iter()
            .map(|c| get_field_with_qcovs(hit, c, qcovs))
            .collect();
        writeln!(writer, "{}", fields.join(delimiter))?;
    }
    Ok(())
}

/// Write tabular output (outfmt 6) for a set of hits.
pub fn format_tabular<W: Write>(writer: &mut W, hits: &[TabularHit]) -> std::io::Result<()> {
    format_tabular_custom_with_delimiter(writer, hits, DEFAULT_TABULAR_COLUMNS, "\t")
}

fn format_query_coverage(cov: f64) -> String {
    rounded_query_coverage(cov).to_string()
}

fn rounded_query_coverage(cov: f64) -> i32 {
    // NCBI-style half-up rounding via `BLAST_Nint`, clamped to 100.
    (crate::math::blast_nint(cov.min(100.0)) as i32).min(100)
}

/// Port of NCBI `Blast_HSPGetQueryCoverage` (`blast_hits.c:1034`), used for the
/// `qcovhsp` column. `query_span` is the inclusive query span of the HSP
/// (`qend - qstart + 1`); `query_len` the full query length. The `< 99`
/// guard on the `+0.5` plus the integer truncation means high coverage values
/// (e.g. 99.6%) truncate to 99 rather than rounding to 100.
fn ncbi_query_coverage(query_span: i64, query_len: i64) -> i32 {
    if query_len <= 0 {
        return 0;
    }
    let mut pct = 100.0 * query_span as f64 / query_len as f64;
    if pct < 99.0 {
        pct += 0.5;
    }
    pct as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(qseq: Option<&str>, sseq: Option<&str>) -> TabularHit {
        TabularHit {
            query_id: "q1".to_string(),
            query_gi: None,
            query_acc: None,
            query_accver: None,
            subject_id: "s1".to_string(),
            subject_seqid: None,
            subject_gi: None,
            subject_acc: None,
            subject_accver: None,
            subject_title: "s1 synthetic title".to_string(),
            pct_identity: 95.0,
            align_len: 50,
            mismatches: 2,
            gap_opens: 1,
            query_start: 1,
            query_end: 50,
            subject_start: 10,
            subject_end: 59,
            subject_is_plain_nucleotide: true,
            evalue: 1e-20,
            bit_score: 80.0,
            query_len: 100,
            subject_len: 200,
            raw_score: 120,
            qseq: qseq.map(|s| s.to_string()),
            sseq: sseq.map(|s| s.to_string()),
            qframe: 1,
            sframe: 0,
            subject_taxids: vec![9606],
            subject_sci_name: "Homo sapiens".to_string(),
            subject_common_name: "human".to_string(),
            subject_blast_name: "primates".to_string(),
            subject_kingdom: "E".to_string(),
            num_ident: 47,
            num_positives: 47,
            num_links: 1,
            comp_adjust_method: 0,
        }
    }

    #[test]
    fn test_qseq_sseq_returned_when_present() {
        let hit = make_hit(Some("ACGTACGT"), Some("ACGT-CGT"));
        assert_eq!(get_field(&hit, "qseq"), "ACGTACGT");
        assert_eq!(get_field(&hit, "sseq"), "ACGT-CGT");
    }

    #[test]
    fn test_qseq_sseq_na_when_none() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "qseq"), "N/A");
        assert_eq!(get_field(&hit, "sseq"), "N/A");
    }

    #[test]
    fn test_subject_title_fields() {
        let mut hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "stitle"), "s1 synthetic title");
        assert_eq!(get_field(&hit, "salltitles"), "s1 synthetic title");
        hit.subject_title.clear();
        assert_eq!(get_field(&hit, "stitle"), "N/A");
        assert_eq!(get_field(&hit, "salltitles"), "N/A");
    }

    #[test]
    fn test_custom_format_includes_aligned_seqs() {
        let hit = make_hit(Some("AAACCC"), Some("AAA-CC"));
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "qseqid qseq sseq").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output.trim(), "q1\tAAACCC\tAAA-CC");
    }

    #[test]
    fn test_custom_delimiter_values_are_literal() {
        let hit = make_hit(None, None);
        let hits = vec![hit];

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "delim=tab qseqid sseqid length").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1tabs1tab50");

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, r"delim=\t qseqid sseqid length").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), r"q1\ts1\t50");

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "delim=space qseqid sseqid length").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1spaces1space50");

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "delim= qseqid sseqid length").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1\ts1\t50");

        let mut buf = Vec::new();
        format_tabular_custom_with_delimiter(&mut buf, &hits, "delim= qseqid sseqid length", ",")
            .unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1,s1,50");

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "delim=, delim=tab qseqid sseqid length").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1,s1,50");

        let mut buf = Vec::new();
        format_tabular_custom_with_delimiter(
            &mut buf,
            &hits,
            "delim=tab delim=, qseqid sseqid length",
            ",",
        )
        .unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1tabs1tab50");
    }

    #[test]
    fn test_custom_columns_match_blast_parser_edges() {
        let hit = make_hit(None, None);
        let hits = vec![hit];

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "qseqid bogus sseqid").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1\ts1");

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "std qlen").unwrap();
        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(
            fields.len(),
            DEFAULT_TABULAR_COLUMNS.split_whitespace().count() + 1
        );
        assert_eq!(fields.last().copied(), Some("100"));

        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "std qaccver saccver qlen").unwrap();
        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(
            fields.len(),
            DEFAULT_TABULAR_COLUMNS.split_whitespace().count() + 1
        );
        assert_eq!(fields.last().copied(), Some("100"));
    }

    #[test]
    fn test_custom_columns_deduplicate_known_fields() {
        let hit = make_hit(None, None);
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "qseqid sseqid qseqid length length").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap().trim(), "q1\ts1\t50");
    }

    #[test]
    fn test_raw_score_field() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "score"), "120");
    }

    #[test]
    fn test_frame_fields() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "qframe"), "1");
        assert_eq!(get_field(&hit, "sframe"), "0");
        assert_eq!(get_field(&hit, "frames"), "1/0");
    }

    #[test]
    fn qcovus_and_sstrand_use_subject_molecule_type_not_letters() {
        let mut protein_like_hit = make_hit(Some("ACGT"), Some("ACGT"));
        protein_like_hit.subject_is_plain_nucleotide = false;

        assert_eq!(get_field(&protein_like_hit, "qcovus"), "N/A");
        assert_eq!(get_field(&protein_like_hit, "sstrand"), "N/A");

        let mut nucleotide_hit = protein_like_hit.clone();
        nucleotide_hit.subject_is_plain_nucleotide = true;
        assert_eq!(get_field(&nucleotide_hit, "qcovus"), "50");
        assert_eq!(get_field(&nucleotide_hit, "sstrand"), "plus");
    }

    #[test]
    fn test_taxid_fields() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "staxid"), "9606");
        assert_eq!(get_field(&hit, "staxids"), "9606");
    }

    #[test]
    fn test_taxid_multiple() {
        let mut hit = make_hit(None, None);
        hit.subject_taxids = vec![9606, 9605];
        assert_eq!(get_field(&hit, "staxid"), "9606");
        assert_eq!(get_field(&hit, "staxids"), "9606;9605");
    }

    #[test]
    fn test_taxids_are_unique_preserving_first_seen_order() {
        let mut hit = make_hit(None, None);
        hit.subject_taxids = vec![9606, 9605, 9606, 9604, 9605];
        assert_eq!(get_field(&hit, "staxid"), "9606");
        assert_eq!(get_field(&hit, "staxids"), "9606;9605;9604");
    }

    #[test]
    fn test_tax_name_fields() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "ssciname"), "Homo sapiens");
        assert_eq!(get_field(&hit, "scomname"), "human");
        assert_eq!(get_field(&hit, "sblastname"), "primates");
        assert_eq!(get_field(&hit, "sskingdom"), "E");
    }

    #[test]
    fn test_tax_name_empty_fallback() {
        let mut hit = make_hit(None, None);
        hit.subject_sci_name = String::new();
        assert_eq!(get_field(&hit, "ssciname"), "N/A");
    }

    #[test]
    fn test_taxid_empty() {
        let mut hit = make_hit(None, None);
        hit.subject_taxids = vec![];
        assert_eq!(get_field(&hit, "staxid"), "N/A");
        assert_eq!(get_field(&hit, "staxids"), "N/A");
    }

    #[test]
    fn test_tabular_all_standard_fields() {
        let hit = make_hit(None, None);
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular(&mut buf, &hits).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(
            fields.len(),
            12,
            "standard tabular format must have 12 fields"
        );
        assert_eq!(fields[0], "q1"); // qseqid
        assert_eq!(fields[1], "s1"); // sseqid
        assert_eq!(fields[2], "95.000"); // pident
        assert_eq!(fields[3], "50"); // length
        assert_eq!(fields[4], "2"); // mismatch
        assert_eq!(fields[5], "1"); // gapopen
        assert_eq!(fields[6], "1"); // qstart
        assert_eq!(fields[7], "50"); // qend
        assert_eq!(fields[8], "10"); // sstart
        assert_eq!(fields[9], "59"); // send
                                     // evalue and bitscore are formatted strings
        assert!(!fields[10].is_empty(), "evalue field must be present");
        assert!(!fields[11].is_empty(), "bitscore field must be present");
    }

    #[test]
    fn test_tabular_custom_columns() {
        let hit = make_hit(Some("ACGTACGT"), Some("ACGT-CGT"));
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "qseq sseq qframe sframe score staxid").unwrap();
        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(fields.len(), 6);
        assert_eq!(fields[0], "ACGTACGT"); // qseq
        assert_eq!(fields[1], "ACGT-CGT"); // sseq
        assert_eq!(fields[2], "1"); // qframe
        assert_eq!(fields[3], "0"); // sframe
        assert_eq!(fields[4], "120"); // score
        assert_eq!(fields[5], "9606"); // staxid
    }

    #[test]
    fn test_tabular_custom_columns_remove_fields() {
        let hit = make_hit(Some("ACGTACGT"), Some("ACGT-CGT"));
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "std -evalue -bitscore qseq").unwrap();
        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        let cols = expanded_column_tokens("std -evalue -bitscore qseq");

        assert_eq!(
            cols,
            vec![
                "qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "qseq",
            ]
        );
        assert_eq!(fields.len(), cols.len());
        assert_eq!(fields.last().copied(), Some("ACGTACGT"));
    }

    #[test]
    fn test_btop_normalizes_lowercase_masked_residues() {
        let hit = make_hit(Some("acGt-A"), Some("ACgtTA"));
        assert_eq!(get_field(&hit, "btop"), "4-T1");
    }

    #[test]
    fn test_positive_fields_use_positive_count_not_identity_count() {
        let mut hit = make_hit(Some("ARND"), Some("AKND"));
        hit.align_len = 4;
        hit.num_ident = 3;
        hit.num_positives = 4;

        assert_eq!(get_field(&hit, "nident"), "3");
        assert_eq!(get_field(&hit, "positive"), "4");
        assert_eq!(get_field(&hit, "ppos"), "100.00");
    }

    #[test]
    fn test_query_coverage_rounds_half_up() {
        let mut hit = make_hit(None, None);
        // Query span = qend-qstart+1 = 20 over qlen = 32 gives 62.5%; since
        // 62.5 < 99 NCBI adds 0.5 (-> 63.0) then truncates to 63. align_len is
        // set differently to confirm qcovhsp uses the query SPAN, not align_len.
        hit.query_len = 32;
        hit.query_start = 1;
        hit.query_end = 20;
        hit.align_len = 25;
        assert_eq!(get_field(&hit, "qcovhsp"), "63");
    }

    #[test]
    fn test_tabular_multiple_hits() {
        let hit1 = make_hit(None, None);
        let mut hit2 = make_hit(None, None);
        hit2.query_id = "q2".to_string();
        hit2.subject_id = "s2".to_string();
        let mut hit3 = make_hit(None, None);
        hit3.query_id = "q3".to_string();
        hit3.subject_id = "s3".to_string();
        let hits = vec![hit1, hit2, hit3];
        let mut buf = Vec::new();
        format_tabular(&mut buf, &hits).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.trim().split('\n').collect();
        assert_eq!(lines.len(), 3, "3 hits should produce 3 lines");
        assert!(lines[0].starts_with("q1\ts1\t"));
        assert!(lines[1].starts_with("q2\ts2\t"));
        assert!(lines[2].starts_with("q3\ts3\t"));
    }

    #[test]
    fn test_tabular_evalue_formatting() {
        // Very small evalue: scientific notation
        assert_eq!(format_evalue(1e-20), "1.00e-20");
        // Tiny evalue below threshold: "0.0"
        assert_eq!(format_evalue(0.0), "0.0");
        // Medium small: fixed with 3 decimals
        assert_eq!(format_evalue(0.005), "0.005");
        // Near 1: fixed with 2 decimals
        assert_eq!(format_evalue(0.5), "0.50");
        // Moderate: fixed with 1 decimal
        assert_eq!(format_evalue(5.0), "5.0");
        // Large: integer
        assert_eq!(format_evalue(100.0), "100");
        // Single-digit exponent should be zero-padded
        assert_eq!(format_evalue(1e-5), "1.00e-05");
    }

    #[test]
    fn test_tabular_identity_calculation() {
        let mut hit = make_hit(None, None);
        // 47 identities out of 50 alignment length = 94.0%
        hit.num_ident = 47;
        hit.align_len = 50;
        hit.pct_identity = 100.0 * 47.0 / 50.0; // 94.0
        assert_eq!(get_field(&hit, "pident"), "94.000");

        // Perfect identity
        hit.pct_identity = 100.0;
        assert_eq!(get_field(&hit, "pident"), "100.000");

        // Low identity
        hit.pct_identity = 33.333;
        assert_eq!(get_field(&hit, "pident"), "33.333");
    }
}
