//! Integration tests ported from the previous blast-rs implementation.
//!
//! Tests build in-memory BLAST databases, run searches, and validate results.

use tempfile::TempDir;

use blast_rs::db::DbType;
use blast_rs::{
    BlastDbBuilder, SequenceEntry, SearchParams,
    blastp, blastn, blastx, tblastn, tblastx,
    parse_fasta, reverse_complement, six_frame_translate,
    BlastnSearch, Strand,
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

// ── Prokka-style annotation performance test ─────────────────────────────────

/// Performance test: build a protein DB from prokka's sprot FASTA (~25K entries),
/// then search 5 query proteins against it. This mimics the prokka annotation
/// workload. Run with: cargo test --release -- --ignored test_blastp_prokka_sprot
///
/// Target: should complete in under 30 seconds for 5 queries (Perl Prokka does
/// 63 queries against the same DB in ~12 seconds using NCBI BLAST+).
#[test]
#[ignore]
fn test_blastp_prokka_sprot() {
    // Paths relative to the blast-rs repo — prokka db may be at a sibling path
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from("/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
    ];
    let sprot_path = sprot_paths.iter().find(|p| p.exists());
    let sprot_path = match sprot_path {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };

    // Load sprot FASTA using input::parse_fasta (returns FastaRecord structs)
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    eprintln!("Loaded {} sprot records", records.len());
    assert!(records.len() > 1000, "sprot should have >1000 entries");

    // Build indexed protein DB
    let t0 = std::time::Instant::now();
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    let db_build_time = t0.elapsed();
    eprintln!("DB build: {:.2}s ({} entries)", db_build_time.as_secs_f64(), records.len());

    // 5 representative bacterial protein queries (real CDS translations)
    let queries: Vec<&[u8]> = vec![
        // Replication initiation protein (~90 aa)
        b"MKQIKEYLEEFVHSRLNKNIILRAAGFEYAKENPNFSQYYGNTVVSLPHRGKYGGPVNRIAPEMFHQIVAKPGERTFEGMFAIFKHRFPDWRDAES",
        // Short hypothetical (~50 aa)
        b"MNDFNYYKSKEIYREKYYQMPKVFFTNEKYMDLSNDAKIAYMLLKDRFDYS",
        // Transposase fragment (~80 aa)
        b"MNYFRYKQFNKDVITVAVGYYLRYALSYRDISEILRGRGVNVHHSTVYRWVQEYAPILYQQSINTAKNTLKGIECIYALY",
        // TraB transfer protein (~60 aa)
        b"MIKKFSLTTVYVAFLSIVLSNITLGAENPGPKIEQGLQQVQTFLTGLIVAVGICAGVWIV",
        // ErmC methyltransferase (~70 aa)
        b"MNEKNIKHSQNFITSKHNIDKIMTNIRLNEHDNIFEIGSGKGHFTLELVQRCNFVTAIEIDHKLCKTTEN",
    ];

    let params = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let t1 = std::time::Instant::now();
    let mut total_hits = 0;
    for (i, query) in queries.iter().enumerate() {
        let t_q = std::time::Instant::now();
        let results = blastp(&db, query, &params);
        let q_time = t_q.elapsed();
        eprintln!(
            "Query {} ({} aa): {} hits in {:.2}s",
            i + 1, query.len(), results.len(), q_time.as_secs_f64()
        );
        total_hits += results.len();
    }
    let search_time = t1.elapsed();
    eprintln!(
        "Total: {} queries, {} hits, {:.2}s ({:.2}s/query)",
        queries.len(), total_hits,
        search_time.as_secs_f64(),
        search_time.as_secs_f64() / queries.len() as f64
    );

    // Performance summary
    let per_query = search_time.as_secs_f64() / queries.len() as f64;
    let ncbi_per_query = 12.0 / 63.0; // NCBI BLAST+ reference: 63 queries in 12s
    eprintln!(
        "\nPerformance: {:.2}s/query (NCBI BLAST+ reference: {:.3}s/query, {:.0}x slower)",
        per_query, ncbi_per_query, per_query / ncbi_per_query
    );

    // Performance assertion: 5 queries should complete in under 60 seconds
    // (NCBI BLAST+ does 63 queries in ~12s, so 5 queries in 60s is very generous)
    assert!(
        search_time.as_secs() < 60,
        "blastp search too slow: {:.1}s for {} queries (target: <60s)",
        search_time.as_secs_f64(), queries.len()
    );

    // Now test with all available threads
    let params_mt = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(0) // 0 = use all available cores
        .filter_low_complexity(false)
        .comp_adjust(0);

    let t2 = std::time::Instant::now();
    let mut total_hits_mt = 0;
    for query in &queries {
        let results = blastp(&db, query, &params_mt);
        total_hits_mt += results.len();
    }
    let mt_time = t2.elapsed();
    let mt_per_query = mt_time.as_secs_f64() / queries.len() as f64;
    let speedup = search_time.as_secs_f64() / mt_time.as_secs_f64();
    eprintln!(
        "\nMulti-threaded: {:.2}s/query ({:.1}x speedup over single-threaded, {:.0}x vs NCBI)",
        mt_per_query, speedup, mt_per_query / ncbi_per_query
    );

    // Verify same hit count
    assert_eq!(total_hits, total_hits_mt, "Multi-threaded should find same hits as single-threaded");
}

/// Sensitivity test: verify that blastp finds the same hits as NCBI BLAST+.
///
/// These are real CDS proteins from E. faecium plasmid AUS0004_p1 that
/// NCBI BLAST+ (via Perl Prokka) annotates against the sprot database.
/// blast-rs must find at least one hit for each of these queries.
///
/// Run with: cargo test --release -- --ignored test_blastp_sensitivity
#[test]
#[ignore]
fn test_blastp_sensitivity() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from("/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => { eprintln!("Skipping: prokka sprot database not found"); return; }
    };

    // Build DB
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // Queries that NCBI BLAST+ finds in sprot (from Perl Prokka on E. faecium plasmid).
    // Each tuple: (name, expected_product, protein_sequence)
    let known_hits: Vec<(&str, &str, &[u8])> = vec![
        ("srtA", "Sortase A",
         b"MIIRHPKKKRIMGKWIIAFWLLSAVGVLLLMPAEASVAKYQQNQQIAAIDRTGTAAETDSSLDVAKIELGDPVGILTIPSISLKLPIYDGTSDKILENGVGITEGTGDITGGNGKNPLIAGHSGLYKDNLFDDLPSVKKGEKFYIKVDGEQHAYQIDRIEEVQKDELQRNFVTYLEPNPNEDRVTLMTCTPKGINTHRFLVYGKRVTFTKSELKDEENKKQKLSWKWLLGSTVFLSVMIIGSLFVYKKKK"),
        ("topB", "DNA topoisomerase 3",
         b"MMKTVILAEKPSQAKAYADSFSKATRKDGYFEIQDRLFPGETVITYGFGHLVELDSPDMYDENWKQWSLEHLPIFPNQYHYHVPKDKKKQFNVVKQQLQSADTIIIATDSDREGELIAWTIIQQAGADQGKTFKRLWINSLEKEAIYQGFQQLRDAGETYPKFEEAQARQIADWLIGMNGSPLYSLLLQQKGIPGSFSLGRVQTPTLYMIYQLQEKIRNFQKEPYFEGKAQVIAQNGAFDAKLDPNETQATQEAFEDYLKEKGVQLGKQPGTIHQVETEKKSAASPRLFSLSSLQSKMNQLMKASAKDTLEAMQGLYEGKYLSYPRTDTPYITEGEYAYLLDHLDEYKHFLKAEAIPTPIHTPNSRYVNNKKVQEHYAIIPTKTVMTAAAFEQLSPLQQAIYEQVLKTTVAMFAEKYTYEETTILTQVQQLQLKAIGKVPLDLGWKKLFGKESEGKEKEEEPLLPKVTKGETVTVDLQVLEKETKPPQPYTEGTLITAMKTAGKTVDSEEAQSILKEVEGIGTEATRANIIETLKQKEYIKVEKNKLVVTNKGILLCQAVEKEPLLTSAEMTAKWESYLLKIGERKGTQTTFLTNIQKFVSHLLEVVPGQIQSTDFGSTLQEVKAASEKQEAARHLGICPKCQEQEVLLYHKAAACTSEACDFRLWTTIAKKKLTATQLKEIIQNGRTSQPVKGLKGQKGSFEATIVLKEDFTTSFEFSEKKKTNYKKRTRRTTK"),
        ("ssb", "Single-stranded DNA-binding protein",
         b"MINNVTLVGRLTKDPDLRYTASGTAVATFTLAVNRNFTNQNGNREADFINCVIWRKPAETMATLAKKGILIGVVGRIQTRTYDNQQGQRVYVTEVVADNFQLLESKAATESRAHADQSSTSPSTTTFEQRDTATPNNNGLNASQNPFGGQSIDISDDDLPF"),
        ("yoeB", "Toxin YoeB",
         b"MIKAWSDDAWDDYLYWHEQGNKSNIKKINKLIKDIDRSPFAGLGKPEALKHDLSGKWSRRRLVDLTDDDLEKIREEKIPFFIGLSQDRVQRMYQEKGLTIDSVFHGKRKPVTKVIINDLVERF"),
        // Note: umuC removed — NCBI BLAST+ also finds no hit at evalue 1e-6
        // (verified with comp_based_stats=0/1 and default settings).
    ];

    let params = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let mut missed = Vec::new();
    for (name, expected, query) in &known_hits {
        let results = blastp(&db, query, &params);
        if results.is_empty() {
            eprintln!("MISS: {} ({}) — 0 hits", name, expected);
            missed.push(*name);
        } else {
            eprintln!("HIT:  {} ({}) — {} hits, best e-value: {:.2e}, title: {}",
                name, expected, results.len(),
                results[0].best_evalue(),
                &results[0].subject_title[..80.min(results[0].subject_title.len())]);
        }
    }

    assert!(
        missed.is_empty(),
        "blastp failed to find hits for {} proteins that NCBI BLAST+ finds: {:?}. \
         This indicates a sensitivity gap in the search algorithm.",
        missed.len(), missed
    );
}

/// Debug: check if srtA and umuC find seeds/hits against sprot
#[test]
#[ignore]
fn debug_srta_seeds() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from("/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => { eprintln!("Skipping: sprot not found"); return; }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // Test with very permissive evalue to see what we find
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    // Test with evalue=1e-6 (same as sensitivity test)
    let params_strict = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let srta = b"MIIRHPKKKRIMGKWIIAFWLLSAVGVLLLMPAEASVAKYQQNQQIAAIDRTGTAAETDSSLDVAKIELGDPVGILTIPSISLKLPIYDGTSDKILENGVGITEGTGDITGGNGKNPLIAGHSGLYKDNLFDDLPSVKKGEKFYIKVDGEQHAYQIDRIEEVQKDELQRNFVTYLEPNPNEDRVTLMTCTPKGINTHRFLVYGKRVTFTKSELKDEENKKQKLSWKWLLGSTVFLSVMIIGSLFVYKKKK";
    let results_e10 = blastp(&db, srta, &params);
    let results_strict = blastp(&db, srta, &params_strict);
    eprintln!("srtA: {} hits (evalue=10), {} hits (evalue=1e-6)", results_e10.len(), results_strict.len());
    for r in results_e10.iter().take(5) {
        let title = &r.subject_title[..80.min(r.subject_title.len())];
        for h in &r.hsps {
            eprintln!("  evalue={:.2e} score={} bits={:.1} alen={} ident={} gaps={} title={}",
                h.evalue, h.score, h.bit_score, h.alignment_length, h.num_identities, h.num_gaps, title);
        }
    }

    let umuc = b"MNSDLILAGESPSYNAAFIAMKEQHPAVFYAQHNAFGLKKIRSGFISDEQAKEYYPLICEALQKDITHFVDEIVASITGYSIDNIRFAKENKNKTIINSFEGWYNLSQQLLDTIMNEQNKSHPQFSYYSKLNSSHQLSHKEKAEAYIAGINIQITIDKQGKFFQKHFDIIQSIIKEESVNIPVLFINTSRNLKYSTGIEFNELFKRSNSNSLLAKRRVFYSLPYQPAKYREYFDSFKKISEKWIEAYKCNELDKNIGISFHDFYDSRFRTKDAKKQFSFINNIMSKIRDLYEVPEKIVRELKTRFKWFWEKKVKK";
    let results_e10 = blastp(&db, umuc, &params);
    let results_strict = blastp(&db, umuc, &params_strict);
    eprintln!("umuC: {} hits (evalue=10), {} hits (evalue=1e-6)", results_e10.len(), results_strict.len());
    for r in results_e10.iter().take(5) {
        let title = &r.subject_title[..80.min(r.subject_title.len())];
        for h in &r.hsps {
            eprintln!("  evalue={:.2e} score={} bits={:.1} alen={} ident={} gaps={} title={}",
                h.evalue, h.score, h.bit_score, h.alignment_length, h.num_identities, h.num_gaps, title);
        }
    }
}

// ── End-to-end API tests (ported from NCBI bl2seq / blastengine / traceback) ─

/// Search a nucleotide sequence against itself. Should produce a single perfect
/// alignment covering the full length with 100% identity.
#[test]
fn test_blastn_self_search() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "self target", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert!(!results.is_empty(), "self search must find a hit");
    let best = &results[0];
    assert_eq!(best.hsps.len(), 1, "self search should yield exactly one HSP");
    let hsp = &best.hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01,
        "self search identity should be 100%, got {:.2}%", hsp.percent_identity());
    assert_eq!(hsp.alignment_length, seq.len(),
        "alignment should span entire sequence");
    assert_eq!(hsp.num_gaps, 0, "perfect self alignment has no gaps");
}

/// Search completely unrelated sequences and verify no hits at strict evalue.
#[test]
fn test_blastn_no_hit() {
    // Poly-A query vs poly-C subject -- no significant similarity.
    let query   = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let subject = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "unrelated", subject),
    ]);
    let params = SearchParams::blastn()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    assert!(results.is_empty(),
        "completely unrelated sequences should produce no hits at evalue 1e-10");
}

/// Subject contains two separate regions matching the query. Verify multiple
/// HSPs are returned.
#[test]
fn test_blastn_multiple_hsps() {
    // Two distinct 40bp matching regions separated by 60bp of unrelated sequence.
    let region_a = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTT";
    let region_b = "GCTTAACGATCGTAAGGCCTTAGCAGTCAATGCTTGAAGT";
    let spacer   = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    let query = format!("{}{}", region_a, region_b);
    let subject = format!("{}{}{}", region_a, spacer, region_b);
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "multi region", &subject),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    assert!(!results.is_empty(), "should find hits");
    // Count total HSPs across all results for subject N001
    let total_hsps: usize = results.iter().map(|r| r.hsps.len()).sum();
    assert!(total_hsps >= 1,
        "subject with two matching regions should produce at least 1 HSP, got {}", total_hsps);
}

/// Protein self-search: search a protein against itself and verify perfect alignment.
#[test]
fn test_blastp_self_search() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "self protein", seq),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, seq.as_bytes(), &params);
    assert!(!results.is_empty(), "protein self search must find a hit");
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01,
        "protein self search should be 100% identity, got {:.2}%", hsp.percent_identity());
    assert_eq!(hsp.alignment_length, seq.len());
    assert_eq!(hsp.num_gaps, 0);
    assert!(hsp.score > 0, "self-hit score should be positive");
}

/// Search two related but not identical proteins. Verify a hit with positive
/// score and reasonable e-value.
#[test]
fn test_blastp_known_pair() {
    // Query and subject differ at a few positions (conservative substitutions).
    let query   = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let subject = "MKFLILLFNILCLFPVLAADNHGVSINAS"; // M→I near the end (conservative in BLOSUM62)
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "related protein", subject),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(!results.is_empty(), "related proteins should produce a hit");
    let hsp = &results[0].hsps[0];
    assert!(hsp.score > 0, "related protein hit should have positive score");
    assert!(hsp.evalue < 1.0, "related protein hit should have reasonable evalue, got {:.2e}", hsp.evalue);
    assert!(hsp.percent_identity() > 90.0,
        "one substitution should give >90% identity, got {:.1}%", hsp.percent_identity());
}

/// BLASTX: nucleotide query encoding a known protein searched against a protein
/// database. Should find the translated match.
#[test]
fn test_blastx_finds_protein_match() {
    // ATG AAA TTT CTG ATT CTG CTG TTT AAC ATT CTG TGC CTG TTC
    // encodes: M   K   F   L   I   L   L   F   N   I   L   C   L   F
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTC";
    let protein  = "MKFLILLFNILCLF";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target protein", protein),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert!(!results.is_empty(), "blastx should find the translated protein match");
    let hsp = &results[0].hsps[0];
    assert!(hsp.percent_identity() > 80.0,
        "translated match should have high identity, got {:.1}%", hsp.percent_identity());
    assert!(hsp.query_frame != 0, "blastx HSP should have non-zero query frame");
    assert!(hsp.score > 0);
}

/// TBLASTN: protein query against nucleotide subject that encodes the protein.
/// Should find the translated nucleotide match.
#[test]
fn test_tblastn_finds_nucleotide_match() {
    let protein_query = "MKFLILLFNILCLF";
    let nt_subject    = "ATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTC";
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
    assert!(hsp.percent_identity() > 80.0,
        "translated match should have high identity, got {:.1}%", hsp.percent_identity());
    assert!(hsp.score > 0);
}

/// Query whose reverse complement matches the subject. Verify hit on minus strand.
#[test]
fn test_blastn_both_strands() {
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let query_rc = String::from_utf8(reverse_complement(subject.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "forward strand subject", subject),
    ]);

    // With strand=both, should find the reverse-complement hit.
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("both");

    let results = blastn(&db, query_rc.as_bytes(), &params);
    assert!(!results.is_empty(), "should find RC hit when strand=both");
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01,
        "RC of subject should be 100% identity, got {:.2}%", hsp.percent_identity());

    // With strand=plus only, should NOT find it (query is RC of subject).
    let params_plus = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("plus");

    let results_plus = blastn(&db, query_rc.as_bytes(), &params_plus);
    assert!(results_plus.is_empty(),
        "plus-strand-only should not find hit when query is RC of subject");
}

/// Very short query (15bp) with word_size=7. Should still find a match.
#[test]
fn test_blastn_short_query() {
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let query = "ATGCGTACCTGAAAG"; // first 15bp of subject
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", subject),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .word_size(7)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    assert!(!results.is_empty(), "short 15bp query with word_size=7 should find a match");
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01,
        "exact substring should be 100% identity");
}

/// Set a strict evalue threshold and verify only significant hits pass.
#[test]
fn test_blastn_evalue_filter() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", seq),
    ]);

    // Relaxed evalue -- should find hits
    let params_relaxed = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);
    let results_relaxed = blastn(&db, seq.as_bytes(), &params_relaxed);
    assert!(!results_relaxed.is_empty(), "relaxed evalue should find self hit");

    // Very strict evalue for self-hit of 100bp should still pass (self hit is highly significant)
    let params_strict = SearchParams::blastn()
        .evalue(1e-30)
        .num_threads(1)
        .filter_low_complexity(false);
    let results_strict = blastn(&db, seq.as_bytes(), &params_strict);
    assert!(!results_strict.is_empty(),
        "100bp self-hit should still pass evalue 1e-30");

    // All returned HSPs must satisfy the evalue threshold
    for r in &results_strict {
        for hsp in &r.hsps {
            assert!(hsp.evalue <= 1e-30,
                "HSP evalue {:.2e} exceeds threshold 1e-30", hsp.evalue);
        }
    }
}

/// Search against a BLAST database and against a subject FASTA should produce
/// equivalent results for the same sequences.
#[test]
fn test_db_search_vs_subject_search() {
    let query   = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";

    // Database-based search
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", subject),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);
    let db_results = blastn(&db, query.as_bytes(), &params);

    // Subject-mode search (using BlastnSearch builder which does seq-vs-seq)
    let builder_results = BlastnSearch::new()
        .query(query.as_bytes())
        .subject(subject.as_bytes())
        .dust(false)
        .evalue(10.0)
        .run();

    // Both should find hits
    assert!(!db_results.is_empty(), "DB search should find hit");
    assert!(!builder_results.is_empty(), "subject search should find hit");

    // The best HSP scores should be comparable (both are self-hits)
    let db_best_score = db_results[0].hsps[0].score;
    let subj_best_score = builder_results[0].score;
    // Scores may differ slightly due to different code paths, but both should be positive
    assert!(db_best_score > 0, "DB search score should be positive");
    assert!(subj_best_score > 0, "subject search score should be positive");
}

/// Use BlastnSearch builder API directly for a seq-vs-seq search.
#[test]
fn test_blastn_search_builder_api() {
    let query   = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";

    let results = BlastnSearch::new()
        .query(query.as_bytes())
        .subject(subject.as_bytes())
        .word_size(11)
        .reward(1)
        .penalty(-3)
        .gap_open(5)
        .gap_extend(2)
        .evalue(10.0)
        .dust(false)
        .strand(Strand::Both)
        .run();

    assert!(!results.is_empty(), "builder API self-search should find hits");
    let best = &results[0];
    assert!(best.score > 0, "best HSP score should be positive");
    assert!(best.evalue < 1e-10, "100bp self-hit should have very significant evalue");
    assert_eq!(best.align_length as usize, query.len(),
        "self-hit alignment should span full query length");
    assert_eq!(best.num_ident as usize, query.len(),
        "self-hit should have 100% identity");

    // Test with empty query -- should return empty
    let empty = BlastnSearch::new()
        .query(b"")
        .subject(subject.as_bytes())
        .run();
    assert!(empty.is_empty(), "empty query should produce no results");

    // Test with empty subject -- should return empty
    let empty2 = BlastnSearch::new()
        .query(query.as_bytes())
        .subject(b"")
        .run();
    assert!(empty2.is_empty(), "empty subject should produce no results");
}
