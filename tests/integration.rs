//! Integration tests ported from the previous blast-rs implementation.
//!
//! Tests build in-memory BLAST databases, run searches, and validate results.

use tempfile::TempDir;

use blast_rs::db::DbType;
use blast_rs::{
    BlastDbBuilder, SequenceEntry, SearchParams,
    blastp, blastn, blastx, tblastn, tblastx,
    parse_fasta, reverse_complement, six_frame_translate,
};

// ── Helpers ──────────────────────────────────────────────────────────────────

fn build_protein_db(entries: Vec<SequenceEntry>) -> (TempDir, blast_rs::db::BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("testdb");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "test protein db");
    for e in entries { builder.add(e); }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    (tmp, db)
}

fn build_nucleotide_db(entries: Vec<SequenceEntry>) -> (TempDir, blast_rs::db::BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("testdb");
    let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "test nt db");
    for e in entries { builder.add(e); }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    (tmp, db)
}

fn protein_entry(acc: &str, title: &str, seq: &str) -> SequenceEntry {
    SequenceEntry {
        title: title.to_string(),
        accession: acc.to_string(),
        sequence: seq.as_bytes().to_vec(),
        taxid: None,
    }
}

fn nt_entry(acc: &str, title: &str, seq: &str) -> SequenceEntry {
    SequenceEntry {
        title: title.to_string(),
        accession: acc.to_string(),
        sequence: seq.as_bytes().to_vec(),
        taxid: None,
    }
}

// ── BLASTP tests ─────────────────────────────────────────────────────────────

#[test]
fn blastp_exact_match() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "exact match protein", seq),
        protein_entry("P002", "unrelated protein", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, seq.as_bytes(), &params);

    assert!(!results.is_empty(), "should find at least one hit");
    let best = &results[0];
    let hsp = &best.hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01,
        "exact match should be 100% identity, got {:.1}%", hsp.percent_identity());
    assert_eq!(hsp.alignment_length, seq.len());
}

#[test]
fn blastp_no_hit_for_unrelated() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "all alanine", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    ]);
    let params = SearchParams::blastp()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"WWWWWWWWWWWWWWWWWWWWWWWWWWWW", &params);
    assert!(results.is_empty(), "unrelated sequences should not produce hits at strict evalue");
}

#[test]
fn blastp_finds_similar_sequence() {
    let query   = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let subject = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let mutated = "MKFLILLFNILCLFPVLAAENHGVSMNAS"; // D→E
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "near identical", subject),
        protein_entry("P002", "one mismatch", mutated),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(results.len() >= 2, "should find both similar sequences, got {}", results.len());
}

#[test]
fn blastp_max_target_seqs_limits_results() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let entries: Vec<_> = (0..10)
        .map(|i| protein_entry(&format!("P{:03}", i), "copy", query))
        .collect();
    let (_tmp, db) = build_protein_db(entries);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .max_target_seqs(3)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(results.len() <= 3,
        "max_target_seqs=3 should limit results to at most 3, got {}", results.len());
}

#[test]
fn blastp_empty_query() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"", &params);
    assert!(results.is_empty(), "empty query should produce no results");
}

#[test]
fn blastp_short_query_below_word_size() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"MK", &params);
    let _ = results; // should not panic
}

#[test]
fn blastp_single_residue_query() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"M", &params);
    let _ = results; // should not panic
}

#[test]
fn blastp_all_twenty_amino_acids() {
    let query = "ACDEFGHIKLMNPQRSTVWY";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "all aa", query),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(!results.is_empty());
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01);
    assert_eq!(hsp.alignment_length, 20);
}

// ── BLASTN tests ─────────────────────────────────────────────────────────────

#[test]
fn blastn_exact_match() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "exact nt", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert!(!results.is_empty(), "should find exact nucleotide match");
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01);
}

#[test]
fn blastn_reverse_complement_hit() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let rc = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "forward strand", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("both");

    let results = blastn(&db, rc.as_bytes(), &params);
    assert!(!results.is_empty(), "should find hit on reverse complement strand");
}

#[test]
fn blastn_no_hit_unrelated() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "seq A", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC", &params);
    assert!(results.is_empty(), "unrelated nucleotide sequences should not match");
}

#[test]
fn blastn_mismatch_scoring() {
    let seq     = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let mutated = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACAATCGTAAGGCCTTAGCAGTCAATGC"; // G→A at pos 76
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "original", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, mutated.as_bytes(), &params);
    assert!(!results.is_empty());
    let hsp = &results[0].hsps[0];
    assert!(hsp.percent_identity() > 90.0, "one mismatch in ~100bp should be >90% identity");
    assert!(hsp.percent_identity() < 100.0, "should not be 100% with mismatch");
}

#[test]
fn blastn_empty_query() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"", &params);
    assert!(results.is_empty(), "empty query should produce no results");
}

#[test]
fn blastn_short_query_below_word_size() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"ATGCGT", &params);
    let _ = results; // should not panic
}

#[test]
fn blastn_plus_strand_only() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let rc = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "forward only", seq),
    ]);

    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("plus");

    let results = blastn(&db, rc.as_bytes(), &params);
    assert!(results.is_empty(), "plus-strand-only should not find reverse complement hit");
}

#[test]
fn blastn_alignment_strings_are_ascii() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert!(!results.is_empty());
    let hsp = &results[0].hsps[0];
    for &b in &hsp.query_aln {
        assert!(b == b'-' || b.is_ascii_alphabetic(),
            "query_aln byte {} is not ASCII letter or gap", b);
    }
    for &b in &hsp.subject_aln {
        assert!(b == b'-' || b.is_ascii_alphabetic(),
            "subject_aln byte {} is not ASCII letter or gap", b);
    }
}

// ── BLASTX test ─────────────────────────────────────────────────────────────

#[test]
fn blastx_finds_translated_hit() {
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTT";
    let protein = "MKFLILLF";

    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target protein", protein),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert!(!results.is_empty(), "blastx should find the translated protein");
    let hsp = &results[0].hsps[0];
    assert!(hsp.percent_identity() > 80.0);
    assert!(hsp.query_frame != 0, "blastx HSP should have a non-zero query frame");
}

#[test]
fn blastx_empty_nt_query() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAAD"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, b"", &params);
    assert!(results.is_empty());
}

// ── TBLASTN test ────────────────────────────────────────────────────────────

#[test]
fn tblastn_finds_protein_in_nt_db() {
    let protein_query = "MKFLILLF";
    let nt_subject = "ATGAAATTTCTGATTCTGCTGTTT";

    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "coding region", nt_subject),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastn(&db, protein_query.as_bytes(), &params);
    assert!(!results.is_empty(), "tblastn should find protein in translated nucleotide db");
    let hsp = &results[0].hsps[0];
    assert!(hsp.subject_frame != 0, "tblastn HSP should have non-zero subject frame");
}

#[test]
fn tblastn_empty_protein_query() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "coding", "ATGAAATTTCTGATTCTGCTGTTT"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastn(&db, b"", &params);
    assert!(results.is_empty());
}

// ── TBLASTX test ────────────────────────────────────────────────────────────

#[test]
fn tblastx_translated_vs_translated() {
    let nt_seq = "ATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "coding nt", nt_seq),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastx(&db, nt_seq.as_bytes(), &params);
    assert!(!results.is_empty(), "tblastx should find self-hit in translated mode");
    let hsp = &results[0].hsps[0];
    assert!(hsp.query_frame != 0, "tblastx should set query frame");
    assert!(hsp.subject_frame != 0, "tblastx should set subject frame");
}

// ── Database round-trip tests ────────────────────────────────────────────────

#[test]
fn protein_db_roundtrip() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "test seq", seq),
    ]);
    assert_eq!(db.num_oids, 1);
}

#[test]
fn nucleotide_db_roundtrip() {
    let seq = "ATGGCTAGCGATCGATCGATCGATCG";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "test nt", seq),
    ]);
    assert_eq!(db.num_oids, 1);
    assert_eq!(db.get_seq_len(0), seq.len() as u32);
}

#[test]
fn db_multiple_sequences() {
    let entries = vec![
        protein_entry("P001", "first", "MKFLILLFNILCLFPVLAAD"),
        protein_entry("P002", "second", "NHGVSMNASQRDHFKLAEV"),
        protein_entry("P003", "third", "ACDEFGHIKLMNPQRSTVWY"),
    ];
    let (_tmp, db) = build_protein_db(entries);
    assert_eq!(db.num_oids, 3);
}

// ── Reverse complement test ─────────────────────────────────────────────────

#[test]
fn reverse_complement_roundtrip() {
    let seq = b"ATGGCTAGCGATCG";
    let rc = reverse_complement(seq);
    let rc2 = reverse_complement(&rc);
    assert_eq!(&rc2, seq, "double reverse complement should return original");
}

// ── Six-frame translation test ──────────────────────────────────────────────

#[test]
fn six_frame_translate_produces_six_frames() {
    let seq = b"ATGGCTAGCGATCGATCGATCGATCG";
    let frames = six_frame_translate(seq);
    assert_eq!(frames.len(), 6);
    let frame_nums: Vec<i32> = frames.iter().map(|f| f.frame).collect();
    assert_eq!(frame_nums, vec![1, 2, 3, -1, -2, -3]);
    // Frame +1 starts with M (ATG = Met)
    assert_eq!(frames[0].protein[0], b'M');
}

// ── FASTA parsing edge cases ─────────────────────────────────────────────────

#[test]
fn parse_fasta_with_blank_lines() {
    let input = b">seq1\n\nACGT\n\nTGCA\n\n>seq2\nAAAA\n";
    let seqs = parse_fasta(input);
    assert_eq!(seqs.len(), 2);
    assert_eq!(&seqs[0].1, b"ACGTTGCA".as_slice());
    assert_eq!(&seqs[1].1, b"AAAA".as_slice());
}

#[test]
fn parse_fasta_no_trailing_newline() {
    let input = b">seq1\nACGT";
    let seqs = parse_fasta(input);
    assert_eq!(seqs.len(), 1);
    assert_eq!(&seqs[0].1, b"ACGT".as_slice());
}

// ── Multi-query search ──────────────────────────────────────────────────────

#[test]
fn multi_query_fasta() {
    let fasta = b">q1\nMKFLILLFNILCLFPVLAAD\n>q2\nNHGVSMNASQRDHFKLAEV\n";
    let queries = parse_fasta(fasta);
    assert_eq!(queries.len(), 2);

    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "match1", "MKFLILLFNILCLFPVLAAD"),
        protein_entry("P002", "match2", "NHGVSMNASQRDHFKLAEV"),
    ]);
    let params = SearchParams::blastp().evalue(10.0).num_threads(1)
        .filter_low_complexity(false).comp_adjust(0);

    for (title, seq) in &queries {
        let results = blastp(&db, seq, &params);
        assert!(!results.is_empty(), "query '{}' should find a hit", title);
    }
}

// ── Edge cases ──────────────────────────────────────────────────────────────

#[test]
fn blastp_short_subject_in_db() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "tiny", "MK"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    let _ = results; // should not panic
}

#[test]
fn blastp_stop_codon_in_sequence() {
    let query   = "MKFLILLFNILCLFPVLAAD";
    let subject = "MKFLILLF*ILCLFPVLAAD";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "with stop", subject),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    let _ = results; // should not panic
}

#[test]
fn blastn_query_with_ambiguous_bases() {
    let query   = "ATGCGTACCTGAAAGCTTCAGNACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", subject),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    let _ = results; // should not panic
}

#[test]
fn db_single_base_nucleotide() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "single", "A"),
    ]);
    assert_eq!(db.num_oids, 1);
}

#[test]
fn db_single_residue_protein() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "single", "M"),
    ]);
    assert_eq!(db.num_oids, 1);
}

// ── Large query test ────────────────────────────────────────────────────────

#[test]
fn blastp_large_query_1000aa() {
    let aa = b"ACDEFGHIKLMNPQRSTVWY";
    let query: Vec<u8> = (0..1000).map(|i| aa[i % aa.len()]).collect();
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", std::str::from_utf8(&query).unwrap()),
        protein_entry("P002", "unrelated", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
    ]);
    let params = SearchParams::blastp().evalue(10.0).num_threads(1)
        .filter_low_complexity(false).comp_adjust(0);
    let results = blastp(&db, &query, &params);
    assert!(!results.is_empty(), "should find self-match for 1000aa query");
}

// ── Blastn with multiple mismatches ─────────────────────────────────────────

#[test]
fn blastn_multiple_mismatches() {
    let query   = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let mut subj_bytes = subject.as_bytes().to_vec();
    subj_bytes[10] = b'T';
    subj_bytes[30] = b'A';
    subj_bytes[50] = b'G';
    subj_bytes[70] = b'C';
    subj_bytes[90] = b'T';
    let subj_str = String::from_utf8(subj_bytes).unwrap();

    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "5 mismatches", &subj_str),
    ]);
    let params = SearchParams::blastn().evalue(10.0).num_threads(1)
        .filter_low_complexity(false);
    let results = blastn(&db, query.as_bytes(), &params);
    assert!(!results.is_empty(), "should find hit with 5 mismatches in 100bp");
    let hsp = &results[0].hsps[0];
    assert!(hsp.percent_identity() > 90.0, "95% identity expected, got {:.1}%", hsp.percent_identity());
    assert!(hsp.percent_identity() < 100.0, "should not be 100% with mismatches");
}

// ── Multi-subject index tests ───────────────────────────────────────────────

#[test]
fn blastp_multi_subject_finds_correct_hit() {
    // 5 different proteins; only one matches the query
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "unrelated1", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
        protein_entry("P002", "unrelated2", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        protein_entry("P003", "the target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
        protein_entry("P004", "unrelated3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        protein_entry("P005", "unrelated4", "GGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
    ]);
    let params = SearchParams::blastp().evalue(10.0).num_threads(1)
        .filter_low_complexity(false).comp_adjust(0);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    assert!(!results.is_empty(), "should find the matching subject");
    // The matching subject should be OID 2 (0-indexed)
    assert_eq!(results[0].subject_oid, 2, "should match OID 2 (the target)");
}

#[test]
fn blastp_multi_subject_finds_multiple_hits() {
    // Two subjects match, three don't
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "match exact", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
        protein_entry("P002", "unrelated", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
        protein_entry("P003", "match similar", "MKFLILLFNILCLFPVLAAENHGVSMNAS"), // one mismatch
        protein_entry("P004", "unrelated", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
    ]);
    let params = SearchParams::blastp().evalue(10.0).num_threads(1)
        .filter_low_complexity(false).comp_adjust(0);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    assert!(results.len() >= 2, "should find at least 2 hits, got {}", results.len());

    let oids: Vec<u32> = results.iter().map(|r| r.subject_oid).collect();
    assert!(oids.contains(&0), "should find OID 0 (exact match)");
    assert!(oids.contains(&2), "should find OID 2 (similar match)");
}

#[test]
fn blastn_multi_subject_finds_correct_hit() {
    let target = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "decoy1", "TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC"),
        nt_entry("N002", "decoy2", "AAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGG"),
        nt_entry("N003", "target", target),
        nt_entry("N004", "decoy3", "GGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATT"),
    ]);
    let params = SearchParams::blastn().evalue(10.0).num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, target.as_bytes(), &params);
    assert!(!results.is_empty(), "should find the target in multi-subject DB");
    assert_eq!(results[0].subject_oid, 2, "should match OID 2 (the target)");
    assert!((results[0].hsps[0].percent_identity() - 100.0).abs() < 0.01);
}

#[test]
fn blastn_multi_subject_many_sequences() {
    // 20 sequences, only one matches
    let target = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let decoy = "TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC";
    let mut entries: Vec<SequenceEntry> = (0..19)
        .map(|i| nt_entry(&format!("N{:03}", i), "decoy", decoy))
        .collect();
    entries.push(nt_entry("N019", "target", target));

    let (_tmp, db) = build_nucleotide_db(entries);
    assert_eq!(db.num_oids, 20);

    let params = SearchParams::blastn().evalue(10.0).num_threads(1)
        .filter_low_complexity(false);
    let results = blastn(&db, target.as_bytes(), &params);
    assert!(!results.is_empty(), "should find target in 20-sequence DB");
    assert_eq!(results[0].subject_oid, 19, "target should be OID 19");
}

#[test]
fn blastp_multi_subject_20_sequences() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let entries: Vec<SequenceEntry> = (0..20)
        .map(|i| protein_entry(&format!("P{:03}", i), "copy", query))
        .collect();

    let (_tmp, db) = build_protein_db(entries);
    assert_eq!(db.num_oids, 20);

    let params = SearchParams::blastp().evalue(10.0).num_threads(1)
        .filter_low_complexity(false).comp_adjust(0);
    let results = blastp(&db, query.as_bytes(), &params);
    // Should find hits in many of the 20 identical subjects
    assert!(results.len() >= 10, "should find many hits in 20-sequence DB, got {}", results.len());
}

#[test]
fn blastx_multi_subject_protein_db() {
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTT"; // encodes MKFLILLF
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "decoy", "WWWWWWWWWWWWWWWWWWWW"),
        protein_entry("P002", "target", "MKFLILLF"),
        protein_entry("P003", "decoy", "AAAAAAAAAAAAAAAAAAA"),
    ]);
    let params = SearchParams::blastp().evalue(10.0).num_threads(1)
        .filter_low_complexity(false).comp_adjust(0);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert!(!results.is_empty(), "blastx should find target in multi-subject protein DB");
    assert_eq!(results[0].subject_oid, 1, "should match OID 1 (the target protein)");
}

#[test]
fn tblastn_multi_subject_nt_db() {
    let protein_query = "MKFLILLF";
    let nt_target = "ATGAAATTTCTGATTCTGCTGTTT"; // encodes MKFLILLF
    let nt_decoy = "TTTGGGCCCAAATTTGGGCCCAAA";

    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "decoy1", nt_decoy),
        nt_entry("N002", "decoy2", nt_decoy),
        nt_entry("N003", "target", nt_target),
    ]);
    let params = SearchParams::blastp().evalue(10.0).num_threads(1)
        .filter_low_complexity(false).comp_adjust(0);

    let results = tblastn(&db, protein_query.as_bytes(), &params);
    assert!(!results.is_empty(), "tblastn should find protein in multi-subject nt DB");
    assert_eq!(results[0].subject_oid, 2, "should match OID 2 (the coding sequence)");
}
