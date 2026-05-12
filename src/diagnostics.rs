//! Rust equivalent of blast_diagnostics.c — search diagnostics and statistics.
//! Tracks counts of lookup hits, extensions, and other search phases.

/// Diagnostics from the ungapped search phase.
#[derive(Debug, Clone, Default)]
pub struct UngappedStats {
    /// Number of successful lookup table hits
    pub lookup_hits: i64,
    /// Number of sequences with at least one lookup hit
    pub num_seqs_lookup_hits: i32,
    /// Number of initial words found and extended
    pub init_extends: i32,
    /// Number of successful initial extensions (HSPs saved)
    pub good_init_extends: i32,
    /// Number of sequences with at least one HSP after ungapped stage
    pub num_seqs_passed: i32,
}

/// Diagnostics from the gapped search phase.
#[derive(Debug, Clone, Default)]
pub struct GappedStats {
    /// Number of sequences that passed the ungapped phase
    pub seqs_ungapped_passed: i32,
    /// Number of sequences with gapped extensions
    pub num_seqs_passed: i32,
    /// Number of gapped extensions performed
    pub extensions: i32,
    /// Number of gapped extensions producing significant HSPs
    pub good_extensions: i32,
}

/// Raw cutoff diagnostics (`BlastRawCutoffs` in `blast_diagnostics.h`).
#[derive(Debug, Clone, Default)]
pub struct RawCutoffs {
    pub x_drop_ungapped: i32,
    pub x_drop_gap: i32,
    pub x_drop_gap_final: i32,
    pub ungapped_cutoff: i32,
    pub cutoff_score: i32,
}

/// Complete search diagnostics.
#[derive(Debug, Clone, Default)]
pub struct Diagnostics {
    pub ungapped: UngappedStats,
    pub gapped: GappedStats,
    pub cutoffs: RawCutoffs,
}

/// NCBI `EBlastSeverity` values (`blast_message.h:55`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum BlastSeverity {
    Info = 1,
    Warning = 2,
    Error = 3,
    Fatal = 4,
}

/// Structure to enclose the origin of an error message or warning.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SMessageOrigin {
    pub filename: String,
    pub lineno: u32,
}

/// Structure to hold a message from the BLAST core.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastMessage {
    pub severity: BlastSeverity,
    pub context: i32,
    pub message: String,
    pub origin: Option<SMessageOrigin>,
    pub next: Option<Box<BlastMessage>>,
}

pub const K_BLAST_MESSAGE_NO_CONTEXT: i32 = -1;
pub const BLASTERR_MEMORY: i16 = 50;
pub const BLASTERR_INVALIDPARAM: i16 = 75;
pub const BLASTERR_IDEALSTATPARAMCALC: i16 = 100;
pub const BLASTERR_REDOALIGNMENTCORE_NOTSUPPORTED: i16 = 101;
pub const BLASTERR_INVALIDQUERIES: i16 = 102;
pub const BLASTERR_INTERRUPTED: i16 = 103;
pub const BLASTERR_NOVALIDKARLINALTSCHUL: i16 = 104;
pub const BLASTERR_SUBJECT_LENGTH_INVALID: i16 = 203;
pub const BLASTERR_SEQSRC: i16 = 300;
pub const BLASTERR_DB_MEMORY_MAP: i16 = 400;
pub const BLASTERR_DB_TOO_MANY_OPEN_FILES: i16 = 401;

pub const K_BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS: &str =
    "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or \
     its translation. Please verify the query sequence(s) and/or filtering options";

impl Diagnostics {
    pub fn new() -> Self {
        blast_diagnostics_init()
    }

    /// Record a lookup hit.
    pub fn add_lookup_hit(&mut self) {
        self.ungapped.lookup_hits += 1;
    }

    /// Record a successful ungapped extension.
    pub fn add_ungapped_hit(&mut self) {
        self.ungapped.good_init_extends += 1;
    }

    /// Record a successful gapped extension.
    pub fn add_gapped_hit(&mut self) {
        self.gapped.good_extensions += 1;
    }

    /// Summary string for logging.
    pub fn summary(&self) -> String {
        format!(
            "Lookup hits: {}, Ungapped extensions: {}/{}, Gapped extensions: {}/{}",
            self.ungapped.lookup_hits,
            self.ungapped.good_init_extends,
            self.ungapped.init_extends,
            self.gapped.good_extensions,
            self.gapped.extensions,
        )
    }
}

/// Port of `SMessageOriginNew` (`blast_message.c:48`).
pub fn s_message_origin_new(filename: &str, lineno: u32) -> Option<SMessageOrigin> {
    if filename.is_empty() {
        return None;
    }
    Some(SMessageOrigin {
        filename: filename.to_string(),
        lineno,
    })
}

/// Port of `SMessageOriginFree` (`blast_message.c:70`).
pub fn s_message_origin_free(_msgo: Option<SMessageOrigin>) -> Option<SMessageOrigin> {
    None
}

/// Port of `Blast_MessageFree` (`blast_message.c:80`).
pub fn blast_message_free(blast_msg: &mut Option<Box<BlastMessage>>) -> Option<Box<BlastMessage>> {
    let mut var_msg = blast_msg.take();
    while let Some(mut msg) = var_msg {
        msg.message.clear();
        msg.origin = s_message_origin_free(msg.origin.take());
        var_msg = msg.next.take();
    }
    None
}

/// Port of `Blast_MessageWrite` (`blast_message.c:102`).
pub fn blast_message_write(
    blast_msg: &mut Option<Box<BlastMessage>>,
    severity: BlastSeverity,
    context: i32,
    message: &str,
) -> i16 {
    let new_msg = Box::new(BlastMessage {
        severity,
        context,
        message: message.to_string(),
        origin: None,
        next: None,
    });
    append_blast_message(blast_msg, new_msg);
    0
}

/// Port of `Blast_MessagePost` (`blast_message.c:136`).
pub fn blast_message_post(blast_msg: Option<&BlastMessage>) -> i16 {
    let Some(msg) = blast_msg else {
        return 1;
    };
    eprint!("{}", msg.message);
    0
}

/// Port of `Blast_PerrorEx` (`blast_message.c:154`).
pub fn blast_perror_ex(
    msg: &mut Option<Box<BlastMessage>>,
    error_code: i16,
    file_name: Option<&str>,
    lineno: i32,
    context: i32,
) {
    let Some((message, severity)) = blast_error_code_to_message(error_code) else {
        return;
    };
    let origin = if let Some(file_name) = file_name {
        if lineno > 0 {
            s_message_origin_new(file_name, lineno as u32)
        } else {
            None
        }
    } else {
        None
    };
    let new_msg = Box::new(BlastMessage {
        severity,
        context,
        message,
        origin,
        next: None,
    });
    append_blast_message(msg, new_msg);
}

fn append_blast_message(list: &mut Option<Box<BlastMessage>>, new_msg: Box<BlastMessage>) {
    let mut cursor = list;
    loop {
        match cursor {
            Some(node) => {
                cursor = &mut node.next;
            }
            None => {
                *cursor = Some(new_msg);
                break;
            }
        }
    }
}

fn blast_error_code_to_message(error_code: i16) -> Option<(String, BlastSeverity)> {
    let pair = match error_code {
        BLASTERR_IDEALSTATPARAMCALC => (
            "Failed to calculate ideal Karlin-Altschul parameters",
            BlastSeverity::Error,
        ),
        BLASTERR_REDOALIGNMENTCORE_NOTSUPPORTED => (
            "Composition based statistics or Smith-Waterman not supported for your program type",
            BlastSeverity::Error,
        ),
        BLASTERR_INTERRUPTED => (
            "BLAST search interrupted at user's request",
            BlastSeverity::Info,
        ),
        BLASTERR_NOVALIDKARLINALTSCHUL => (
            K_BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS,
            BlastSeverity::Error,
        ),
        BLASTERR_MEMORY => ("Out of memory", BlastSeverity::Fatal),
        BLASTERR_INVALIDPARAM => ("Invalid argument to function", BlastSeverity::Fatal),
        BLASTERR_INVALIDQUERIES => (
            "search cannot proceed due to errors in all contexts/frames of query sequences",
            BlastSeverity::Fatal,
        ),
        BLASTERR_SUBJECT_LENGTH_INVALID => (
            "The average subject length is too short",
            BlastSeverity::Fatal,
        ),
        BLASTERR_SEQSRC => (
            "search cannot proceed due to errors retrieving sequences from databases",
            BlastSeverity::Fatal,
        ),
        BLASTERR_DB_MEMORY_MAP => ("Database memory map file error", BlastSeverity::Fatal),
        BLASTERR_DB_TOO_MANY_OPEN_FILES => (
            "Too many open files, please raise the open file limit",
            BlastSeverity::Fatal,
        ),
        0 => return None,
        _ => {
            return Some((
                format!("Unknown error code {}", error_code),
                BlastSeverity::Error,
            ));
        }
    };
    Some((pair.0.to_string(), pair.1))
}

/// 1-1 port of `Blast_DiagnosticsInit` (`blast_diagnostics.c:75`).
pub fn blast_diagnostics_init() -> Diagnostics {
    Diagnostics::default()
}

/// 1-1 port of `Blast_DiagnosticsInitMT` (`blast_diagnostics.c:89`).
/// Rust callers synchronize externally, so this is the same initialized value.
pub fn blast_diagnostics_init_mt() -> Diagnostics {
    blast_diagnostics_init()
}

/// 1-1 port of `Blast_DiagnosticsCopy` (`blast_diagnostics.c:48`).
pub fn blast_diagnostics_copy(diagnostics: Option<&Diagnostics>) -> Option<Diagnostics> {
    let Some(src) = diagnostics else {
        return None;
    };
    let mut dst = Diagnostics::default();
    dst.ungapped = src.ungapped.clone();
    dst.gapped = src.gapped.clone();
    dst.cutoffs = src.cutoffs.clone();
    Some(dst)
}

/// 1-1 port of `Blast_DiagnosticsFree` (`blast_diagnostics.c:35`).
/// Dropping the value releases owned Rust storage; clearing the slot mirrors
/// the C function returning `NULL`.
pub fn blast_diagnostics_free(diagnostics: &mut Option<Diagnostics>) -> Option<Diagnostics> {
    *diagnostics = None;
    None
}

/// 1-1 port of `Blast_UngappedStatsUpdate` (`blast_diagnostics.c:102`).
pub fn blast_ungapped_stats_update(
    ungapped_stats: Option<&mut UngappedStats>,
    total_hits: i32,
    extended_hits: i32,
    saved_hits: i32,
) {
    let Some(stats) = ungapped_stats else { return };
    if total_hits == 0 {
        return;
    }
    stats.lookup_hits += total_hits as i64;
    stats.num_seqs_lookup_hits += 1;
    stats.init_extends += extended_hits;
    stats.good_init_extends += saved_hits;
    if saved_hits > 0 {
        stats.num_seqs_passed += 1;
    }
}

/// 1-1 port of `Blast_DiagnosticsUpdate` (`blast_diagnostics.c:119`).
pub fn blast_diagnostics_update(global: &mut Diagnostics, local: Option<&Diagnostics>) {
    let Some(local) = local else { return };

    global.ungapped.lookup_hits += local.ungapped.lookup_hits;
    global.ungapped.num_seqs_lookup_hits += local.ungapped.num_seqs_lookup_hits;
    global.ungapped.init_extends += local.ungapped.init_extends;
    global.ungapped.good_init_extends += local.ungapped.good_init_extends;
    global.ungapped.num_seqs_passed += local.ungapped.num_seqs_passed;

    global.gapped.seqs_ungapped_passed += local.gapped.seqs_ungapped_passed;
    global.gapped.extensions += local.gapped.extensions;
    global.gapped.good_extensions += local.gapped.good_extensions;
    global.gapped.num_seqs_passed += local.gapped.num_seqs_passed;

    global.cutoffs = local.cutoffs.clone();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diagnostics() {
        let mut diag = Diagnostics::new();
        diag.add_lookup_hit();
        diag.add_lookup_hit();
        diag.add_ungapped_hit();
        assert_eq!(diag.ungapped.lookup_hits, 2);
        assert_eq!(diag.ungapped.good_init_extends, 1);
        assert!(!diag.summary().is_empty());
    }

    #[test]
    fn test_diagnostics_default_zeros() {
        let diag = Diagnostics::new();
        assert_eq!(diag.ungapped.lookup_hits, 0);
        assert_eq!(diag.ungapped.num_seqs_lookup_hits, 0);
        assert_eq!(diag.ungapped.init_extends, 0);
        assert_eq!(diag.ungapped.good_init_extends, 0);
        assert_eq!(diag.ungapped.num_seqs_passed, 0);
        assert_eq!(diag.gapped.seqs_ungapped_passed, 0);
        assert_eq!(diag.gapped.num_seqs_passed, 0);
        assert_eq!(diag.gapped.extensions, 0);
        assert_eq!(diag.gapped.good_extensions, 0);
        assert_eq!(diag.cutoffs.cutoff_score, 0);
    }

    #[test]
    fn test_diagnostics_gapped_tracking() {
        let mut diag = Diagnostics::new();
        diag.add_gapped_hit();
        diag.add_gapped_hit();
        diag.add_gapped_hit();
        assert_eq!(diag.gapped.good_extensions, 3);
        // Ungapped should remain zero
        assert_eq!(diag.ungapped.good_init_extends, 0);
    }

    #[test]
    fn test_diagnostics_mixed_tracking() {
        let mut diag = Diagnostics::new();
        for _ in 0..100 {
            diag.add_lookup_hit();
        }
        for _ in 0..10 {
            diag.add_ungapped_hit();
        }
        for _ in 0..3 {
            diag.add_gapped_hit();
        }
        assert_eq!(diag.ungapped.lookup_hits, 100);
        assert_eq!(diag.ungapped.good_init_extends, 10);
        assert_eq!(diag.gapped.good_extensions, 3);
    }

    #[test]
    fn test_diagnostics_summary_format() {
        let mut diag = Diagnostics::new();
        diag.add_lookup_hit();
        diag.ungapped.init_extends = 5;
        diag.add_ungapped_hit();
        diag.gapped.extensions = 2;
        diag.add_gapped_hit();
        let s = diag.summary();
        assert!(s.contains("Lookup hits: 1"));
        assert!(s.contains("Ungapped extensions: 1/5"));
        assert!(s.contains("Gapped extensions: 1/2"));
    }

    #[test]
    fn test_ungapped_stats_clone() {
        let stats = UngappedStats {
            lookup_hits: 42,
            num_seqs_passed: 7,
            ..Default::default()
        };
        let cloned = stats.clone();
        assert_eq!(cloned.lookup_hits, 42);
        assert_eq!(cloned.num_seqs_passed, 7);
    }

    #[test]
    fn test_gapped_stats_clone() {
        let stats = GappedStats {
            extensions: 10,
            good_extensions: 3,
            ..Default::default()
        };
        let cloned = stats.clone();
        assert_eq!(cloned.extensions, 10);
        assert_eq!(cloned.good_extensions, 3);
    }

    #[test]
    fn test_diagnostics_direct_field_mutation() {
        let mut diag = Diagnostics::new();
        diag.ungapped.num_seqs_lookup_hits = 5;
        diag.ungapped.num_seqs_passed = 3;
        diag.gapped.num_seqs_passed = 2;
        diag.gapped.extensions = 8;
        assert_eq!(diag.ungapped.num_seqs_lookup_hits, 5);
        assert_eq!(diag.ungapped.num_seqs_passed, 3);
        assert_eq!(diag.gapped.num_seqs_passed, 2);
        assert_eq!(diag.gapped.extensions, 8);
    }

    #[test]
    fn test_blast_diagnostics_update_accumulates_stats_and_copies_cutoffs() {
        let mut global = blast_diagnostics_init();
        let mut local = blast_diagnostics_init();
        local.ungapped.lookup_hits = 11;
        local.ungapped.num_seqs_lookup_hits = 2;
        local.ungapped.init_extends = 7;
        local.ungapped.good_init_extends = 3;
        local.ungapped.num_seqs_passed = 1;
        local.gapped.seqs_ungapped_passed = 4;
        local.gapped.extensions = 5;
        local.gapped.good_extensions = 2;
        local.gapped.num_seqs_passed = 1;
        local.cutoffs.cutoff_score = 42;

        blast_diagnostics_update(&mut global, Some(&local));

        assert_eq!(global.ungapped.lookup_hits, 11);
        assert_eq!(global.ungapped.num_seqs_lookup_hits, 2);
        assert_eq!(global.ungapped.init_extends, 7);
        assert_eq!(global.ungapped.good_init_extends, 3);
        assert_eq!(global.ungapped.num_seqs_passed, 1);
        assert_eq!(global.gapped.seqs_ungapped_passed, 4);
        assert_eq!(global.gapped.extensions, 5);
        assert_eq!(global.gapped.good_extensions, 2);
        assert_eq!(global.gapped.num_seqs_passed, 1);
        assert_eq!(global.cutoffs.cutoff_score, 42);
    }

    #[test]
    fn test_blast_diagnostics_copy_and_free() {
        let mut diag = blast_diagnostics_init_mt();
        diag.ungapped.lookup_hits = 9;
        let copied = blast_diagnostics_copy(Some(&diag)).expect("copy");
        assert_eq!(copied.ungapped.lookup_hits, 9);
        assert!(blast_diagnostics_copy(None).is_none());

        let mut slot = Some(diag);
        assert!(blast_diagnostics_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn test_blast_message_helpers_match_c_lifecycle() {
        assert!(s_message_origin_new("", 10).is_none());
        let origin = s_message_origin_new("blast_message.c", 154).expect("origin");
        assert_eq!(origin.filename, "blast_message.c");
        assert_eq!(origin.lineno, 154);
        assert!(s_message_origin_free(Some(origin)).is_none());

        let mut messages = None;
        assert_eq!(
            blast_message_write(
                &mut messages,
                BlastSeverity::Warning,
                K_BLAST_MESSAGE_NO_CONTEXT,
                "first"
            ),
            0
        );
        assert_eq!(
            blast_message_write(&mut messages, BlastSeverity::Error, 3, "second"),
            0
        );
        let first = messages.as_ref().expect("first message");
        assert_eq!(first.message, "first");
        assert_eq!(first.severity, BlastSeverity::Warning);
        assert_eq!(first.next.as_ref().expect("second").message, "second");
        assert_eq!(blast_message_post(None), 1);

        blast_perror_ex(
            &mut messages,
            BLASTERR_INTERRUPTED,
            Some("blast_engine.c"),
            99,
            4,
        );
        let third = messages
            .as_ref()
            .unwrap()
            .next
            .as_ref()
            .unwrap()
            .next
            .as_ref()
            .expect("third message");
        assert_eq!(third.severity, BlastSeverity::Info);
        assert_eq!(third.context, 4);
        assert_eq!(third.message, "BLAST search interrupted at user's request");
        assert_eq!(third.origin.as_ref().unwrap().filename, "blast_engine.c");

        assert!(third.next.is_none());
        blast_perror_ex(&mut messages, 0, Some("unused.c"), 1, 0);
        let third_after_zero = messages
            .as_ref()
            .unwrap()
            .next
            .as_ref()
            .unwrap()
            .next
            .as_ref()
            .expect("third message");
        assert!(third_after_zero.next.is_none());

        assert!(blast_message_free(&mut messages).is_none());
        assert!(messages.is_none());
    }
}
