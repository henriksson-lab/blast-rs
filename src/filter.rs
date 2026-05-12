//! Rust equivalent of blast_filter.c — query sequence masking/filtering.
//! Implements DUST (nucleotide) and SEG (protein) low-complexity filtering.

use std::collections::VecDeque;

use crate::program;
use crate::queryinfo::QueryInfo;
use crate::util::{BlastSequenceBlk, SSeqRange};

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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastSeqLoc {
    pub ssl: SSeqRange,
    pub next: Option<Box<BlastSeqLoc>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastMaskLoc {
    pub masks: Vec<Option<Box<BlastSeqLoc>>>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SReadQualityOptions {
    pub frac_ambig: f64,
    pub entropy: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SDustOptions {
    pub level: i32,
    pub window: i32,
    pub linker: i32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SSegOptions {
    pub window: i32,
    pub locut: f64,
    pub hicut: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SRepeatFilterOptions {
    pub database: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SWindowMaskerOptions {
    pub taxid: i32,
    pub database: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EFilterOptions {
    ESeg,
    EDust,
    ERepeats,
    EDustRepeats,
    EEmpty,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SBlastFilterOptions {
    pub mask_at_hash: bool,
    pub dust_options: Option<SDustOptions>,
    pub seg_options: Option<SSegOptions>,
    pub repeat_filter_options: Option<SRepeatFilterOptions>,
    pub window_masker_options: Option<SWindowMaskerOptions>,
    pub read_quality_options: Option<SReadQualityOptions>,
}

pub const K_DUST_LEVEL: i32 = 20;
pub const K_DUST_WINDOW: i32 = 64;
pub const K_DUST_LINKER: i32 = 1;
pub const K_SEG_WINDOW: i32 = 12;
pub const K_SEG_LOCUT: f64 = 2.2;
pub const K_SEG_HICUT: f64 = 2.5;
pub const K_DEFAULT_REPEAT_FILTER_DB: &str = "repeat/repeat_9606";
pub const K_NUCL_MASK: u8 = 14;
pub const K_PROT_MASK: u8 = crate::encoding::NCBISTDAA_X;

pub fn sdust_options_free(dust_options: &mut Option<SDustOptions>) -> Option<SDustOptions> {
    *dust_options = None;
    None
}

pub fn sdust_options_new(dust_options: &mut Option<SDustOptions>) -> i16 {
    *dust_options = Some(SDustOptions {
        level: K_DUST_LEVEL,
        window: K_DUST_WINDOW,
        linker: K_DUST_LINKER,
    });
    0
}

pub fn sseg_options_free(seg_options: &mut Option<SSegOptions>) -> Option<SSegOptions> {
    *seg_options = None;
    None
}

pub fn sseg_options_new(seg_options: &mut Option<SSegOptions>) -> i16 {
    *seg_options = Some(SSegOptions {
        window: K_SEG_WINDOW,
        locut: K_SEG_LOCUT,
        hicut: K_SEG_HICUT,
    });
    0
}

pub fn srepeat_filter_options_free(
    repeat_options: &mut Option<SRepeatFilterOptions>,
) -> Option<SRepeatFilterOptions> {
    if let Some(options) = repeat_options {
        options.database = None;
    }
    *repeat_options = None;
    None
}

pub fn srepeat_filter_options_new(repeat_options: &mut Option<SRepeatFilterOptions>) -> i16 {
    *repeat_options = Some(SRepeatFilterOptions {
        database: Some(K_DEFAULT_REPEAT_FILTER_DB.to_string()),
    });
    0
}

pub fn srepeat_filter_options_reset_db(
    repeat_options: &mut Option<SRepeatFilterOptions>,
    db: &str,
) -> i16 {
    if repeat_options.is_none() && srepeat_filter_options_new(repeat_options) != 0 {
        return 1;
    }

    if let Some(options) = repeat_options {
        options.database = Some(db.to_string());
    }
    0
}

pub fn swindow_masker_options_free(
    winmask_options: &mut Option<SWindowMaskerOptions>,
) -> Option<SWindowMaskerOptions> {
    if let Some(options) = winmask_options {
        options.database = None;
    }
    *winmask_options = None;
    None
}

pub fn swindow_masker_options_new(winmask_options: &mut Option<SWindowMaskerOptions>) -> i16 {
    *winmask_options = Some(SWindowMaskerOptions {
        taxid: 0,
        database: None,
    });
    0
}

pub fn swindow_masker_options_reset_db(
    winmask_options: &mut Option<SWindowMaskerOptions>,
    db: Option<&str>,
) -> i16 {
    if winmask_options.is_none() && swindow_masker_options_new(winmask_options) != 0 {
        return 1;
    }

    if let Some(options) = winmask_options {
        options.database = db.map(str::to_string);
    }
    0
}

pub fn sread_quality_options_free(
    read_quality_options: &mut Option<SReadQualityOptions>,
) -> Option<SReadQualityOptions> {
    *read_quality_options = None;
    None
}

pub fn sread_quality_options_new(read_quality_options: &mut Option<SReadQualityOptions>) -> i16 {
    *read_quality_options = Some(SReadQualityOptions {
        frac_ambig: 0.5,
        entropy: 16,
    });
    0
}

pub fn sblast_filter_options_free(
    filter_options: &mut Option<SBlastFilterOptions>,
) -> Option<SBlastFilterOptions> {
    if let Some(options) = filter_options {
        sdust_options_free(&mut options.dust_options);
        sseg_options_free(&mut options.seg_options);
        srepeat_filter_options_free(&mut options.repeat_filter_options);
        swindow_masker_options_free(&mut options.window_masker_options);
        sread_quality_options_free(&mut options.read_quality_options);
    }
    *filter_options = None;
    None
}

pub fn sblast_filter_options_new(
    filter_options: &mut Option<SBlastFilterOptions>,
    filter_type: EFilterOptions,
) -> i16 {
    let mut options = SBlastFilterOptions {
        mask_at_hash: false,
        dust_options: None,
        seg_options: None,
        repeat_filter_options: None,
        window_masker_options: None,
        read_quality_options: None,
    };

    if filter_type == EFilterOptions::ESeg {
        sseg_options_new(&mut options.seg_options);
    }
    if matches!(
        filter_type,
        EFilterOptions::EDust | EFilterOptions::EDustRepeats
    ) {
        sdust_options_new(&mut options.dust_options);
    }
    if matches!(
        filter_type,
        EFilterOptions::ERepeats | EFilterOptions::EDustRepeats
    ) {
        srepeat_filter_options_new(&mut options.repeat_filter_options);
    }

    *filter_options = Some(options);
    0
}

fn filter_string_is_nucleotide_program(program_number: program::ProgramType) -> bool {
    program_number == program::BLASTN || program_number == program::MAPPING
}

fn filter_string_token_is_keyword(token: &str) -> bool {
    matches!(
        token,
        "d" | "dust"
            | "r"
            | "repeat"
            | "repeats"
            | "s"
            | "seg"
            | "m"
            | "f"
            | "false"
            | "none"
            | "t"
            | "true"
            | "yes"
            | "windowmasker"
            | "window_masker"
            | "window-masker"
            | "wm"
            | "windowmasker_taxid"
            | "read_quality"
            | "read-quality"
            | "rq"
    )
}

fn filter_string_take_i32(tokens: &[String], index: &mut usize) -> Option<i32> {
    let next = tokens.get(*index)?;
    if filter_string_token_is_keyword(next) {
        return None;
    }
    let value = next.parse().ok()?;
    *index += 1;
    Some(value)
}

fn filter_string_take_f64(tokens: &[String], index: &mut usize) -> Option<f64> {
    let next = tokens.get(*index)?;
    if filter_string_token_is_keyword(next) {
        return None;
    }
    let value = next.parse().ok()?;
    *index += 1;
    Some(value)
}

fn filter_string_take_text(tokens: &[String], index: &mut usize) -> Option<String> {
    let next = tokens.get(*index)?;
    if filter_string_token_is_keyword(next) {
        return None;
    }
    *index += 1;
    Some(next.clone())
}

fn filter_string_set_default_filter(
    options: &mut SBlastFilterOptions,
    program_number: program::ProgramType,
) {
    if filter_string_is_nucleotide_program(program_number) {
        if options.dust_options.is_none() {
            let mut dust = None;
            sdust_options_new(&mut dust);
            options.dust_options = dust;
        }
    } else if options.seg_options.is_none() {
        let mut seg = None;
        sseg_options_new(&mut seg);
        options.seg_options = seg;
    }
}

/// Port of NCBI internal `s_SafeStrCat` (`blast_filter.c:301`).
pub fn s_safe_str_cat(dest: &mut String, dest_size: usize, string2append: &str) -> Option<usize> {
    let needed = dest.len().checked_add(string2append.len())?;
    if needed >= dest_size {
        return None;
    }
    dest.push_str(string2append);
    Some(dest.len())
}

/// Port of NCBI internal `s_LoadOptionsToBuffer` (`blast_filter.c:54`).
pub fn s_load_options_to_buffer(
    buffer: &mut String,
    dest_size: usize,
    option_name: &str,
    option_values: &[String],
) -> i16 {
    if option_name.is_empty() {
        return crate::util::BLASTERR_INVALIDPARAM;
    }

    let mut fields = Vec::with_capacity(option_values.len() + 1);
    fields.push(option_name);
    fields.extend(option_values.iter().map(String::as_str));

    for field in fields {
        if !buffer.is_empty() && s_safe_str_cat(buffer, dest_size, " ").is_none() {
            return crate::util::BLASTERR_INVALIDPARAM;
        }
        if s_safe_str_cat(buffer, dest_size, field).is_none() {
            return crate::util::BLASTERR_INVALIDPARAM;
        }
    }
    0
}

/// Port of NCBI internal `s_ParseSegOptions` (`blast_filter.c:242`).
pub fn s_parse_seg_options(tokens: &[String], index: &mut usize) -> Result<SSegOptions, i16> {
    let mut seg = SSegOptions {
        window: K_SEG_WINDOW,
        locut: K_SEG_LOCUT,
        hicut: K_SEG_HICUT,
    };

    if *index >= tokens.len() || filter_string_token_is_keyword(&tokens[*index]) {
        return Ok(seg);
    }

    let saved = *index;
    let Some(window) = filter_string_take_i32(tokens, index) else {
        *index = saved;
        return Err(crate::util::BLASTERR_INVALIDPARAM);
    };
    let Some(locut) = filter_string_take_f64(tokens, index) else {
        *index = saved;
        return Err(crate::util::BLASTERR_INVALIDPARAM);
    };
    let Some(hicut) = filter_string_take_f64(tokens, index) else {
        *index = saved;
        return Err(crate::util::BLASTERR_INVALIDPARAM);
    };

    if window <= 0 || locut < 0.0 || hicut < 0.0 || locut > hicut {
        *index = saved;
        return Err(crate::util::BLASTERR_INVALIDPARAM);
    }

    seg.window = window;
    seg.locut = locut;
    seg.hicut = hicut;
    Ok(seg)
}

/// Port of NCBI `BlastFilteringOptionsToString` (`blast_filter.c`).
pub fn blast_filtering_options_to_string(
    program_number: program::ProgramType,
    filter_options: Option<&SBlastFilterOptions>,
) -> Option<String> {
    let options = filter_options?;
    let mut buffer = String::new();
    let dest_size = usize::MAX;

    if let Some(dust) = options.dust_options {
        if dust
            == (SDustOptions {
                level: K_DUST_LEVEL,
                window: K_DUST_WINDOW,
                linker: K_DUST_LINKER,
            })
        {
            s_load_options_to_buffer(&mut buffer, dest_size, "D", &[]);
        } else {
            s_load_options_to_buffer(
                &mut buffer,
                dest_size,
                "D",
                &[
                    dust.level.to_string(),
                    dust.window.to_string(),
                    dust.linker.to_string(),
                ],
            );
        }
    }

    if let Some(seg) = options.seg_options {
        if seg
            == (SSegOptions {
                window: K_SEG_WINDOW,
                locut: K_SEG_LOCUT,
                hicut: K_SEG_HICUT,
            })
        {
            s_load_options_to_buffer(&mut buffer, dest_size, "S", &[]);
        } else {
            s_load_options_to_buffer(
                &mut buffer,
                dest_size,
                "S",
                &[
                    seg.window.to_string(),
                    seg.locut.to_string(),
                    seg.hicut.to_string(),
                ],
            );
        }
    }

    if let Some(repeat) = &options.repeat_filter_options {
        match repeat.database.as_deref() {
            Some(K_DEFAULT_REPEAT_FILTER_DB) | None => {
                s_load_options_to_buffer(&mut buffer, dest_size, "R", &[]);
            }
            Some(database) => {
                s_load_options_to_buffer(&mut buffer, dest_size, "R", &[database.to_string()]);
            }
        }
    }

    if let Some(window_masker) = &options.window_masker_options {
        if window_masker.taxid != 0 {
            s_load_options_to_buffer(
                &mut buffer,
                dest_size,
                "windowmasker_taxid",
                &[window_masker.taxid.to_string()],
            );
        }
        if let Some(database) = window_masker.database.as_deref() {
            s_load_options_to_buffer(
                &mut buffer,
                dest_size,
                "windowmasker",
                &[database.to_string()],
            );
        }
        if window_masker.taxid == 0 && window_masker.database.is_none() {
            s_load_options_to_buffer(&mut buffer, dest_size, "windowmasker", &[]);
        }
    }

    if let Some(read_quality) = options.read_quality_options {
        s_load_options_to_buffer(
            &mut buffer,
            dest_size,
            "read_quality",
            &[
                read_quality.frac_ambig.to_string(),
                read_quality.entropy.to_string(),
            ],
        );
    }

    if options.mask_at_hash {
        s_load_options_to_buffer(&mut buffer, dest_size, "m", &[]);
    }

    if buffer.is_empty() {
        if filter_string_is_nucleotide_program(program_number) {
            Some("F".to_string())
        } else {
            Some("none".to_string())
        }
    } else {
        Some(buffer)
    }
}

/// Port of NCBI `BlastFilteringOptionsFromString` (`blast_filter.c`).
pub fn blast_filtering_options_from_string(
    program_number: program::ProgramType,
    instructions: Option<&str>,
    filtering_options: &mut Option<SBlastFilterOptions>,
) -> i16 {
    let Some(instructions) = instructions else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    let normalized = instructions.trim().to_ascii_lowercase();
    let tokens: Vec<String> = normalized
        .split(|ch: char| ch.is_ascii_whitespace() || matches!(ch, ',' | ';' | '(' | ')'))
        .map(|token| token.trim_matches('"').trim_matches('\''))
        .filter(|token| !token.is_empty())
        .map(str::to_string)
        .collect();

    if tokens.is_empty()
        || tokens
            .iter()
            .all(|token| matches!(token.as_str(), "f" | "false" | "none"))
    {
        return sblast_filter_options_new(filtering_options, EFilterOptions::EEmpty);
    }

    let mut parsed = None;
    let status = sblast_filter_options_new(&mut parsed, EFilterOptions::EEmpty);
    if status != 0 {
        return status;
    }
    let Some(options) = parsed.as_mut() else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    let mut saw_filter_instruction = false;
    let mut index = 0;
    while index < tokens.len() {
        let token = tokens[index].as_str();
        index += 1;

        match token {
            "f" | "false" | "none" => {}
            "t" | "true" | "yes" => {
                saw_filter_instruction = true;
                filter_string_set_default_filter(options, program_number);
            }
            "m" => {
                saw_filter_instruction = true;
                options.mask_at_hash = true;
            }
            "d" | "dust" => {
                if !filter_string_is_nucleotide_program(program_number) {
                    return crate::options::BLASTERR_OPTION_PROGRAM_INVALID;
                }
                saw_filter_instruction = true;
                let mut dust = SDustOptions {
                    level: K_DUST_LEVEL,
                    window: K_DUST_WINDOW,
                    linker: K_DUST_LINKER,
                };
                let saved = index;
                if let (Some(level), Some(window), Some(linker)) = (
                    filter_string_take_i32(&tokens, &mut index),
                    filter_string_take_i32(&tokens, &mut index),
                    filter_string_take_i32(&tokens, &mut index),
                ) {
                    dust = SDustOptions {
                        level,
                        window,
                        linker,
                    };
                } else {
                    index = saved;
                }
                options.dust_options = Some(dust);
            }
            "r" | "repeat" | "repeats" => {
                if !filter_string_is_nucleotide_program(program_number) {
                    return crate::options::BLASTERR_OPTION_PROGRAM_INVALID;
                }
                saw_filter_instruction = true;
                let database = filter_string_take_text(&tokens, &mut index)
                    .unwrap_or_else(|| K_DEFAULT_REPEAT_FILTER_DB.to_string());
                options.repeat_filter_options = Some(SRepeatFilterOptions {
                    database: Some(database),
                });
            }
            "s" | "seg" => {
                if filter_string_is_nucleotide_program(program_number) {
                    return crate::options::BLASTERR_OPTION_PROGRAM_INVALID;
                }
                saw_filter_instruction = true;
                match s_parse_seg_options(&tokens, &mut index) {
                    Ok(seg) => options.seg_options = Some(seg),
                    Err(status) => return status,
                }
            }
            "windowmasker" | "window_masker" | "window-masker" | "wm" => {
                saw_filter_instruction = true;
                let mut winmask = SWindowMaskerOptions {
                    taxid: 0,
                    database: None,
                };
                if let Some(taxid) = filter_string_take_i32(&tokens, &mut index) {
                    winmask.taxid = taxid;
                } else {
                    winmask.database = filter_string_take_text(&tokens, &mut index);
                }
                options.window_masker_options = Some(winmask);
            }
            "windowmasker_taxid" => {
                saw_filter_instruction = true;
                let Some(taxid) = filter_string_take_i32(&tokens, &mut index) else {
                    return crate::util::BLASTERR_INVALIDPARAM;
                };
                options.window_masker_options = Some(SWindowMaskerOptions {
                    taxid,
                    database: None,
                });
            }
            "read_quality" | "read-quality" | "rq" => {
                saw_filter_instruction = true;
                let mut read_quality = SReadQualityOptions {
                    frac_ambig: 0.5,
                    entropy: 16,
                };
                let saved = index;
                if let (Some(frac_ambig), Some(entropy)) = (
                    filter_string_take_f64(&tokens, &mut index),
                    filter_string_take_i32(&tokens, &mut index),
                ) {
                    read_quality = SReadQualityOptions {
                        frac_ambig,
                        entropy,
                    };
                } else {
                    index = saved;
                }
                options.read_quality_options = Some(read_quality);
            }
            _ => return crate::util::BLASTERR_INVALIDPARAM,
        }
    }

    if !saw_filter_instruction {
        return crate::util::BLASTERR_INVALIDPARAM;
    }

    let status = sblast_filter_options_validate(program_number, Some(options));
    if status != 0 {
        return status;
    }
    *filtering_options = parsed;
    0
}

pub fn s_merge_dust_options(
    opt1: Option<&SDustOptions>,
    opt2: Option<&SDustOptions>,
) -> Option<SDustOptions> {
    match (opt1, opt2) {
        (None, None) => None,
        (Some(options), None) | (None, Some(options)) => Some(*options),
        (Some(first), Some(second)) => Some(SDustOptions {
            level: if first.level != K_DUST_LEVEL {
                first.level
            } else {
                second.level
            },
            window: if first.window != K_DUST_WINDOW {
                first.window
            } else {
                second.window
            },
            linker: if first.linker != K_DUST_LINKER {
                first.linker
            } else {
                second.linker
            },
        }),
    }
}

pub fn s_merge_seg_options(
    opt1: Option<&SSegOptions>,
    opt2: Option<&SSegOptions>,
) -> Option<SSegOptions> {
    match (opt1, opt2) {
        (None, None) => None,
        (Some(options), None) | (None, Some(options)) => Some(*options),
        (Some(first), Some(second)) => Some(SSegOptions {
            window: if first.window != K_SEG_WINDOW {
                first.window
            } else {
                second.window
            },
            locut: if first.locut != K_SEG_LOCUT {
                first.locut
            } else {
                second.locut
            },
            hicut: if first.hicut != K_SEG_HICUT {
                first.hicut
            } else {
                second.hicut
            },
        }),
    }
}

pub fn s_merge_repeat_options(
    opt1: Option<&SRepeatFilterOptions>,
    opt2: Option<&SRepeatFilterOptions>,
) -> Option<SRepeatFilterOptions> {
    match (opt1, opt2) {
        (None, None) => None,
        (Some(options), None) | (None, Some(options)) => Some(options.clone()),
        (Some(_), Some(second)) => Some(second.clone()),
    }
}

pub fn s_merge_window_masker_options(
    opt1: Option<&SWindowMaskerOptions>,
    opt2: Option<&SWindowMaskerOptions>,
) -> Option<SWindowMaskerOptions> {
    let have1 = opt1.is_some_and(|options| options.database.is_some() || options.taxid != 0);
    let have2 = opt2.is_some_and(|options| options.database.is_some() || options.taxid != 0);

    match (have1, have2) {
        (false, false) => None,
        (true, false) => opt1.cloned(),
        (false, true) | (true, true) => opt2.cloned(),
    }
}

pub fn sblast_filter_options_merge(
    combined: &mut Option<SBlastFilterOptions>,
    opt1: Option<&SBlastFilterOptions>,
    opt2: Option<&SBlastFilterOptions>,
) -> i16 {
    *combined = None;
    if opt1.is_none() && opt2.is_none() {
        return 0;
    }

    let mut merged = None;
    let status = sblast_filter_options_new(&mut merged, EFilterOptions::EEmpty);
    if status != 0 {
        return status;
    }

    if let Some(options) = &mut merged {
        options.mask_at_hash = opt1.is_some_and(|options| options.mask_at_hash)
            || opt2.is_some_and(|options| options.mask_at_hash);
        options.dust_options = s_merge_dust_options(
            opt1.and_then(|options| options.dust_options.as_ref()),
            opt2.and_then(|options| options.dust_options.as_ref()),
        );
        options.seg_options = s_merge_seg_options(
            opt1.and_then(|options| options.seg_options.as_ref()),
            opt2.and_then(|options| options.seg_options.as_ref()),
        );
        options.repeat_filter_options = s_merge_repeat_options(
            opt1.and_then(|options| options.repeat_filter_options.as_ref()),
            opt2.and_then(|options| options.repeat_filter_options.as_ref()),
        );
        options.window_masker_options = s_merge_window_masker_options(
            opt1.and_then(|options| options.window_masker_options.as_ref()),
            opt2.and_then(|options| options.window_masker_options.as_ref()),
        );
        options.read_quality_options = opt2
            .and_then(|options| options.read_quality_options)
            .or_else(|| opt1.and_then(|options| options.read_quality_options));
    }

    *combined = merged;
    0
}

pub fn sblast_filter_options_no_filtering(filter_options: Option<&SBlastFilterOptions>) -> bool {
    let Some(options) = filter_options else {
        return true;
    };

    options.dust_options.is_none()
        && options.seg_options.is_none()
        && options.repeat_filter_options.is_none()
        && options.window_masker_options.is_none()
        && options.read_quality_options.is_none()
}

pub fn sblast_filter_options_mask_at_hash(filter_options: Option<&SBlastFilterOptions>) -> bool {
    filter_options.is_some_and(|options| options.mask_at_hash)
}

/// Port of NCBI `SBlastFilterOptionsValidate` (`blast_options.c:441`).
pub fn sblast_filter_options_validate(
    program_number: crate::program::ProgramType,
    filter_options: Option<&SBlastFilterOptions>,
) -> i16 {
    let Some(options) = filter_options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if let Some(repeat_options) = &options.repeat_filter_options {
        if program_number != program::BLASTN && program_number != program::MAPPING {
            return crate::options::BLASTERR_OPTION_PROGRAM_INVALID;
        }
        if repeat_options.database.as_deref().is_none_or(str::is_empty) {
            return crate::util::BLASTERR_INVALIDPARAM;
        }
    }

    if options.dust_options.is_some()
        && program_number != program::BLASTN
        && program_number != program::MAPPING
    {
        return crate::options::BLASTERR_OPTION_PROGRAM_INVALID;
    }

    if options.seg_options.is_some() && filter_string_is_nucleotide_program(program_number) {
        return crate::options::BLASTERR_OPTION_PROGRAM_INVALID;
    }

    0
}

/// Port of NCBI `BlastSeqLocNew` (`blast_filter.c:608`).
pub fn blast_seq_loc_new(head: &mut Option<Box<BlastSeqLoc>>, start: i32, stop: i32) -> i16 {
    let node = Box::new(BlastSeqLoc {
        ssl: SSeqRange {
            left: start,
            right: stop,
        },
        next: None,
    });
    blast_seq_loc_append(head, Some(node))
}

/// Port of NCBI `BlastSeqLocAppend` (`blast_filter.c:621`).
pub fn blast_seq_loc_append(
    head: &mut Option<Box<BlastSeqLoc>>,
    node: Option<Box<BlastSeqLoc>>,
) -> i16 {
    if node.is_none() {
        return 0;
    }

    if head.is_none() {
        *head = node;
        return 0;
    }

    let mut tail = head.as_mut();
    while let Some(current) = tail {
        if current.next.is_none() {
            current.next = node;
            return 0;
        }
        tail = current.next.as_mut();
    }
    0
}

/// Port of NCBI internal `s_BlastSeqLocNodeDup` (`blast_filter.c:651`).
pub fn s_blast_seq_loc_node_dup(node: Option<&BlastSeqLoc>) -> Option<Box<BlastSeqLoc>> {
    let node = node?;
    Some(Box::new(BlastSeqLoc {
        ssl: node.ssl,
        next: None,
    }))
}

/// Port of NCBI internal `s_BlastSeqLocLen` (`blast_filter.c:664`).
pub fn s_blast_seq_loc_len(head: Option<&BlastSeqLoc>) -> usize {
    let mut length = 0;
    let mut cursor = head;
    while let Some(node) = cursor {
        length += 1;
        cursor = node.next.as_deref();
    }
    length
}

/// Port of NCBI internal `s_BlastSeqLocListToArrayOfPointers` (`blast_filter.c:681`).
pub fn s_blast_seq_loc_list_to_array_of_pointers(head: Option<&BlastSeqLoc>) -> Vec<&BlastSeqLoc> {
    let mut array = Vec::with_capacity(s_blast_seq_loc_len(head));
    let mut cursor = head;
    while let Some(node) = cursor {
        array.push(node);
        cursor = node.next.as_deref();
    }
    array
}

/// Port of NCBI `BlastSeqLocListReverse` (`blast_filter.c:704`).
pub fn blast_seq_loc_list_reverse(head: &mut Option<Box<BlastSeqLoc>>) -> i16 {
    let mut previous = None;
    let mut current = head.take();

    while let Some(mut node) = current {
        let next = node.next.take();
        node.next = previous;
        previous = Some(node);
        current = next;
    }

    *head = previous;
    0
}

/// Port of NCBI `BlastSeqLocNodeFree` (`blast_filter.c:727`).
pub fn blast_seq_loc_node_free(node: &mut Option<Box<BlastSeqLoc>>) -> Option<Box<BlastSeqLoc>> {
    if let Some(node) = node {
        node.next = None;
    }
    *node = None;
    None
}

/// Port of NCBI `BlastSeqLocFree` (`blast_filter.c:737`).
pub fn blast_seq_loc_free(head: &mut Option<Box<BlastSeqLoc>>) -> Option<Box<BlastSeqLoc>> {
    while head.is_some() {
        let next = head.as_mut().and_then(|node| node.next.take());
        *head = next;
    }
    None
}

/// Port of NCBI `BlastSeqLoc_RestrictToInterval` (`blast_setup.c:1030`).
pub fn blast_seq_loc_restrict_to_interval(
    mask: Option<&mut Option<Box<BlastSeqLoc>>>,
    from: i32,
    mut to: i32,
) {
    to = to.max(0);
    let Some(mask) = mask else {
        return;
    };
    if mask.is_none() || (from == 0 && to == 0) {
        return;
    }

    let mut old = mask.take();
    let mut kept: Option<Box<BlastSeqLoc>> = None;
    while let Some(mut node) = old {
        old = node.next.take();
        node.ssl.left = 0.max(node.ssl.left - from);
        node.ssl.right = node.ssl.right.min(to) - from;
        if node.ssl.left <= node.ssl.right {
            let _ = blast_seq_loc_append(&mut kept, Some(node));
        }
    }
    *mask = kept;
}

/// Port of NCBI `s_SeqLocListInvert` (`blast_nalookup.c:330`).
/// Returns the complement intervals between a non-empty, ordered mask list.
pub fn s_seq_loc_list_invert(
    locations: Option<&BlastSeqLoc>,
    length: i32,
) -> Option<Box<BlastSeqLoc>> {
    let mut cursor = locations?;
    let mut retval = None;

    let mut start = 0;
    let mut stop = 0.max(cursor.ssl.left - 1);
    if stop - start > 2 {
        let _ = blast_seq_loc_new(&mut retval, start, stop);
    }

    loop {
        start = cursor.ssl.right + 1;
        if let Some(next) = cursor.next.as_deref() {
            stop = next.ssl.left - 1;
            append_seq_loc_tail(&mut retval, start, stop);
            cursor = next;
        } else {
            stop = length - 1;
            append_seq_loc_tail(&mut retval, start, stop);
            break;
        }
    }

    retval
}

fn append_seq_loc_tail(head: &mut Option<Box<BlastSeqLoc>>, start: i32, stop: i32) {
    if start > stop {
        return;
    }
    if head.is_none() {
        let _ = blast_seq_loc_new(head, start, stop);
    } else {
        let _ = blast_seq_loc_append(
            head,
            Some(Box::new(BlastSeqLoc {
                ssl: SSeqRange {
                    left: start,
                    right: stop,
                },
                next: None,
            })),
        );
    }
}

/// Port of NCBI `BlastSeqLocListDup` (`blast_filter.c:747`).
pub fn blast_seq_loc_list_dup(head: Option<&BlastSeqLoc>) -> Option<Box<BlastSeqLoc>> {
    let mut copied = None;
    let mut cursor = head;
    while let Some(node) = cursor {
        let _ = blast_seq_loc_append(&mut copied, s_blast_seq_loc_node_dup(Some(node)));
        cursor = node.next.as_deref();
    }
    copied
}

/// Port of NCBI `BlastMaskLocNew` (`blast_filter.c:760`).
pub fn blast_mask_loc_new(num_contexts: usize) -> Option<BlastMaskLoc> {
    if num_contexts == 0 {
        return None;
    }
    Some(BlastMaskLoc {
        masks: vec![None; num_contexts],
    })
}

/// Port of NCBI `BlastMaskLocDup` (`blast_filter.c:770`).
pub fn blast_mask_loc_dup(mask_loc: Option<&BlastMaskLoc>) -> Option<BlastMaskLoc> {
    let mask_loc = mask_loc?;
    let mut duplicate = BlastMaskLoc {
        masks: Vec::with_capacity(mask_loc.masks.len()),
    };
    for list in &mask_loc.masks {
        duplicate
            .masks
            .push(blast_seq_loc_list_dup(list.as_deref()));
    }
    Some(duplicate)
}

/// Port of NCBI `BlastMaskLocFree` (`blast_filter.c:789`).
pub fn blast_mask_loc_free(mask_loc: &mut Option<BlastMaskLoc>) -> Option<BlastMaskLoc> {
    if let Some(mask_loc) = mask_loc {
        for list in &mut mask_loc.masks {
            blast_seq_loc_free(list);
        }
        mask_loc.masks.clear();
    }
    *mask_loc = None;
    None
}

/// Port of NCBI `BlastMaskLocDNAToProtein` (`blast_filter.c:806`).
///
/// C gets each original nucleotide query length through
/// `BlastQueryInfoGetQueryLength`; Rust callers pass those lengths explicitly
/// because `QueryInfo` stores only per-context translated lengths.
pub fn blast_mask_loc_dna_to_protein(
    mask_loc: Option<&mut BlastMaskLoc>,
    query_info: &QueryInfo,
    dna_lengths: &[i32],
) -> i16 {
    let Some(mask_loc) = mask_loc else {
        return 0;
    };
    assert_eq!(mask_loc.masks.len(), query_info.contexts.len());
    assert!(dna_lengths.len() >= query_info.num_queries.max(0) as usize);

    for seq_index in 0..query_info.num_queries.max(0) as usize {
        let ctx_idx = crate::util::NUM_FRAMES * seq_index;
        let dna_length = dna_lengths[seq_index];
        let mut dna_seqlocs: Vec<Option<Box<BlastSeqLoc>>> = (0..crate::util::NUM_FRAMES)
            .map(|context| mask_loc.masks[ctx_idx + context].take())
            .collect();

        for context in 0..crate::util::NUM_FRAMES {
            let frame = crate::util::blast_context_to_frame_blastx(context as u32);
            let fallback = dna_seqlocs[0].as_deref();
            let mut iter = dna_seqlocs[context].as_deref().or(fallback);

            while let Some(node) = iter {
                let mut from;
                let mut to;
                if frame < 0 {
                    from = (dna_length + frame - node.ssl.right) / crate::util::CODON_LENGTH as i32;
                    to = (dna_length + frame - node.ssl.left) / crate::util::CODON_LENGTH as i32;
                } else {
                    from = (node.ssl.left - frame + 1) / crate::util::CODON_LENGTH as i32;
                    to = (node.ssl.right - frame + 1) / crate::util::CODON_LENGTH as i32;
                }

                if from < 0 {
                    from = 0;
                }
                if to < 0 {
                    to = 0;
                }

                let query_length = query_info.contexts[ctx_idx + context].query_length;
                if from >= query_length {
                    from = query_length - 1;
                }
                if to >= query_length {
                    to = query_length - 1;
                }

                blast_seq_loc_append(
                    &mut mask_loc.masks[ctx_idx + context],
                    Some(Box::new(BlastSeqLoc {
                        ssl: SSeqRange {
                            left: from,
                            right: to,
                        },
                        next: None,
                    })),
                );
                iter = node.next.as_deref();
            }
        }

        for list in &mut dna_seqlocs {
            blast_seq_loc_free(list);
        }
    }

    0
}

/// Port of NCBI `BlastMaskLocProteinToDNA` (`blast_filter.c:892`).
pub fn blast_mask_loc_protein_to_dna(
    mask_loc: Option<&mut BlastMaskLoc>,
    query_info: &QueryInfo,
    dna_lengths: &[i32],
) -> i16 {
    let Some(mask_loc) = mask_loc else {
        return 0;
    };
    assert_eq!(mask_loc.masks.len(), query_info.contexts.len());
    assert!(dna_lengths.len() >= query_info.num_queries.max(0) as usize);

    for index in 0..query_info.num_queries.max(0) as usize {
        let frame_start = index * crate::util::NUM_FRAMES;
        let dna_length = dna_lengths[index];

        for frame_index in frame_start..frame_start + crate::util::NUM_FRAMES {
            let frame = crate::util::blast_context_to_frame_blastx(
                (frame_index % crate::util::NUM_FRAMES) as u32,
            );
            let mut loc = mask_loc.masks[frame_index].as_deref_mut();
            while let Some(node) = loc {
                let (mut from, mut to) = if frame < 0 {
                    (
                        dna_length - crate::util::CODON_LENGTH as i32 * node.ssl.right + frame - 2,
                        dna_length - crate::util::CODON_LENGTH as i32 * node.ssl.left + frame,
                    )
                } else {
                    (
                        crate::util::CODON_LENGTH as i32 * node.ssl.left + frame - 1,
                        crate::util::CODON_LENGTH as i32 * node.ssl.right + frame + 1,
                    )
                };

                if from < 0 {
                    from = 0;
                }
                if to < 0 {
                    to = 0;
                }
                if from >= dna_length {
                    from = dna_length - 1;
                }
                if to >= dna_length {
                    to = dna_length - 1;
                }

                node.ssl.left = from;
                node.ssl.right = to;
                loc = node.next.as_deref_mut();
            }
        }
    }

    0
}

/// Port of NCBI `BlastIsReverseStrand` macro use in filtering paths.
///
/// Nucleotide query contexts alternate plus/minus strands, so odd contexts are
/// reverse. Protein/translated contexts are not treated as nucleotide strands
/// in `blast_filter.c`.
pub fn blast_is_reverse_strand(is_nucl: bool, context: usize) -> bool {
    is_nucl && context % crate::util::NUM_STRANDS == 1
}

/// Port of NCBI `BlastSeqLocCombine` (`blast_filter.c:971`).
pub fn blast_seq_loc_combine(mask_loc: &mut Option<Box<BlastSeqLoc>>, link_value: i32) {
    let mut ranges = Vec::new();
    let mut cursor = mask_loc.as_deref();
    while let Some(node) = cursor {
        ranges.push(node.ssl);
        cursor = node.next.as_deref();
    }

    if ranges.len() < 2 {
        return;
    }

    ranges.sort_by_key(|range| (range.left, range.right));
    let mut merged: Vec<SSeqRange> = Vec::with_capacity(ranges.len());
    for range in ranges {
        if range.left > range.right {
            continue;
        }

        if let Some(last) = merged.last_mut() {
            if range.left <= last.right.saturating_add(link_value).saturating_add(1) {
                last.right = last.right.max(range.right);
                continue;
            }
        }
        merged.push(range);
    }

    blast_seq_loc_free(mask_loc);
    for range in merged {
        blast_seq_loc_append(
            mask_loc,
            Some(Box::new(BlastSeqLoc {
                ssl: range,
                next: None,
            })),
        );
    }
}

/// Port of NCBI `BLAST_ComplementMaskLocations` (`blast_filter.c:1016`).
pub fn blast_complement_mask_locations(
    program_number: crate::program::ProgramType,
    query_info: &QueryInfo,
    mask_loc: Option<&BlastMaskLoc>,
    complement_mask: Option<&mut Option<Box<BlastSeqLoc>>>,
) -> i16 {
    let Some(complement_mask) = complement_mask else {
        return -1;
    };
    *complement_mask = None;

    let is_nucl = program_number == program::BLASTN || program_number == program::MAPPING;
    for context in 0..query_info.contexts.len() {
        let context_info = &query_info.contexts[context];
        if !context_info.is_valid {
            continue;
        }

        let start_offset = context_info.query_offset;
        let end_offset = context_info.query_length + start_offset - 1;
        assert!(start_offset <= end_offset);

        let mask = mask_loc
            .and_then(|mask_loc| mask_loc.masks.get(context))
            .and_then(|mask| mask.as_deref());

        let Some(mask) = mask else {
            append_seq_loc_tail(complement_mask, start_offset, end_offset);
            continue;
        };

        let reverse = blast_is_reverse_strand(is_nucl, context);
        let mut locs = s_blast_seq_loc_list_to_array_of_pointers(Some(mask));
        if reverse {
            locs.reverse();
        }

        let mut first = true;
        let mut last_interval_open = true;
        let mut left = 0;

        for loc in locs {
            let (filter_start, filter_end) = if reverse {
                (end_offset - loc.ssl.right, end_offset - loc.ssl.left)
            } else {
                (start_offset + loc.ssl.left, start_offset + loc.ssl.right)
            };

            if first {
                last_interval_open = true;
                first = false;

                if filter_start > start_offset {
                    left = start_offset;
                } else {
                    left = filter_end + 1;
                    continue;
                }
            }

            let right = filter_start - 1;
            append_seq_loc_tail(complement_mask, left, right);
            if filter_end >= end_offset {
                last_interval_open = false;
                break;
            } else {
                left = filter_end + 1;
            }
        }

        if last_interval_open {
            append_seq_loc_tail(complement_mask, left, end_offset);
        }
    }

    0
}

/// Port of NCBI `BlastSeqLocReverse` (`blast_filter.c:1173`).
pub fn blast_seq_loc_reverse(mut masks: Option<&mut BlastSeqLoc>, query_length: i32) {
    while let Some(node) = masks {
        let old_left = node.ssl.left;
        let old_right = node.ssl.right;
        node.ssl.left = query_length - 1 - old_right;
        node.ssl.right = query_length - 1 - old_left;
        masks = node.next.as_deref_mut();
    }
}

/// Port of NCBI `BlastSetUp_MaskQuery` (`blast_filter.c:1349`).
pub fn blast_setup_mask_query(
    query_blk: &mut BlastSequenceBlk,
    query_info: &QueryInfo,
    filter_maskloc: &BlastMaskLoc,
    program_number: crate::program::ProgramType,
) {
    let is_nucl = program_number == program::BLASTN || program_number == program::MAPPING;
    if !filter_maskloc.masks.iter().any(Option::is_some) {
        return;
    }

    let Some(sequence) = query_blk.sequence.as_mut() else {
        return;
    };
    let total_length = query_info
        .contexts
        .last()
        .map(|context| context.query_offset + context.query_length + 2)
        .unwrap_or(0)
        .max(0) as usize;

    query_blk.sequence_start_nomask = query_blk
        .sequence_start
        .as_ref()
        .map(|sequence_start| sequence_start[..total_length.min(sequence_start.len())].to_vec());
    query_blk.sequence_nomask = query_blk
        .sequence_start_nomask
        .as_ref()
        .map(|sequence_start| {
            if sequence_start.len() > 1 {
                sequence_start[1..].to_vec()
            } else {
                Vec::new()
            }
        });
    query_blk.nomask_allocated = true;

    for context in 0..query_info.contexts.len() {
        let context_info = &query_info.contexts[context];
        if !context_info.is_valid {
            continue;
        }

        let query_length = context_info.query_length.max(0) as usize;
        let context_offset = context_info.query_offset.max(0) as usize;
        if context_offset >= sequence.len() {
            continue;
        }
        let end = context_offset
            .saturating_add(query_length)
            .min(sequence.len());
        blast_mask_the_residues(
            &mut sequence[context_offset..end],
            is_nucl,
            filter_maskloc
                .masks
                .get(context)
                .and_then(|mask| mask.as_deref()),
            blast_is_reverse_strand(is_nucl, context),
            0,
        );
    }
}

/// Port of NCBI `Blast_MaskTheResidues` (`blast_filter.c:1306`).
pub fn blast_mask_the_residues(
    buffer: &mut [u8],
    is_na: bool,
    mut mask_loc: Option<&BlastSeqLoc>,
    reverse: bool,
    offset: i32,
) {
    let masking_letter = if is_na { K_NUCL_MASK } else { K_PROT_MASK };
    let length = buffer.len() as i32;

    while let Some(node) = mask_loc {
        let (mut start, mut stop) = if reverse {
            (length - 1 - node.ssl.right, length - 1 - node.ssl.left)
        } else {
            (node.ssl.left, node.ssl.right)
        };

        start -= offset;
        stop -= offset;

        assert!(start >= 0);
        assert!(start < length);
        assert!(stop < length);

        for index in start..=stop {
            buffer[index as usize] = masking_letter;
        }

        mask_loc = node.next.as_deref();
    }
}

fn append_mask_loc_regions(seq_loc: &mut Option<Box<BlastSeqLoc>>, mask: &MaskLoc) {
    for region in &mask.regions {
        if region.start < region.end {
            blast_seq_loc_append(
                seq_loc,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange {
                        left: region.start,
                        right: region.end - 1,
                    },
                    next: None,
                })),
            );
        }
    }
}

/// Port of NCBI `BlastSetUp_Filter` (`blast_filter.c:1121`).
pub fn blast_setup_filter(
    program_number: crate::program::ProgramType,
    sequence: &[u8],
    offset: i32,
    filter_options: Option<&SBlastFilterOptions>,
    seqloc_retval: &mut Option<Box<BlastSeqLoc>>,
) -> i16 {
    *seqloc_retval = None;

    let status = sblast_filter_options_validate(program_number, filter_options);
    if status != 0 {
        return status;
    }
    let filter_options = filter_options.expect("validated non-null filter options");

    if let Some(seg_options) = filter_options.seg_options {
        let mut seg_params = seg_parameters_new_aa();
        seg_params.overlaps = true;
        if seg_options.window > 0 {
            seg_params.window = seg_options.window as usize;
        }
        if seg_options.locut > 0.0 {
            seg_params.locut = seg_options.locut;
        }
        if seg_options.hicut > 0.0 {
            seg_params.hicut = seg_options.hicut;
        }

        let (status, mask) = seq_buffer_seg(sequence, offset, Some(seg_params));
        if status != 0 {
            return status;
        }
        append_mask_loc_regions(seqloc_retval, &mask);
        seg_parameters_free(Some(seg_params));
    }

    if let Some(read_quality_options) = filter_options.read_quality_options {
        let mut seq_range = None;
        let status =
            filter_queries_for_mapping(sequence, offset, &read_quality_options, &mut seq_range);
        if status != 0 {
            return status;
        }
        if let Some(range) = seq_range {
            blast_seq_loc_append(
                seqloc_retval,
                Some(Box::new(BlastSeqLoc {
                    ssl: range,
                    next: None,
                })),
            );
        }
    }

    0
}

/// Port of NCBI internal `s_GetFilteringLocationsForOneContext`
/// (`blast_filter.c:1191`).
///
/// Rust keeps lower-case masks outside `BlastSequenceBlk`, so callers pass the
/// optional mask table explicitly; this function consumes the context entry
/// just like C detaches `query_blk->lcase_mask->seqloc_array[context]`.
pub fn s_get_filtering_locations_for_one_context(
    query_blk: &mut BlastSequenceBlk,
    query_info: &QueryInfo,
    context: usize,
    program_number: crate::program::ProgramType,
    filter_options: Option<&SBlastFilterOptions>,
    filter_out: &mut Option<Box<BlastSeqLoc>>,
    lcase_mask: Option<&mut BlastMaskLoc>,
) -> i16 {
    let Some(context_info) = query_info.contexts.get(context) else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    let context_offset = context_info.query_offset.max(0) as usize;
    let Some(sequence) = query_blk.sequence.as_ref() else {
        *filter_out = None;
        return 0;
    };
    if context_offset > sequence.len() {
        *filter_out = None;
        return crate::util::BLASTERR_INVALIDPARAM;
    }

    if !context_info.is_valid {
        *filter_out = None;
        return 0;
    }

    let query_length = context_info.query_length.max(0) as usize;
    let end = context_offset
        .saturating_add(query_length)
        .min(sequence.len());
    let status = blast_setup_filter(
        program_number,
        &sequence[context_offset..end],
        0,
        filter_options,
        filter_out,
    );
    if status != 0 {
        return status;
    }

    let is_nucl = program_number == program::BLASTN || program_number == program::MAPPING;
    if blast_is_reverse_strand(is_nucl, context) {
        blast_seq_loc_reverse(filter_out.as_deref_mut(), context_info.query_length);
    }

    if let Some(lcase_mask) = lcase_mask {
        if context < lcase_mask.masks.len() {
            let lcase_mask_slp = lcase_mask.masks[context].take();
            blast_seq_loc_append(filter_out, lcase_mask_slp);
        }
    }

    blast_seq_loc_combine(filter_out, 0);
    0
}

/// Port of NCBI `BlastSetUp_GetFilteringLocations` (`blast_filter.c:1261`).
pub fn blast_setup_get_filtering_locations(
    query_blk: &mut BlastSequenceBlk,
    query_info: &QueryInfo,
    program_number: crate::program::ProgramType,
    filter_options: Option<&SBlastFilterOptions>,
    filter_maskloc: &mut Option<BlastMaskLoc>,
    mut lcase_mask: Option<&mut BlastMaskLoc>,
) -> i16 {
    let num_contexts = query_info.contexts.len();
    *filter_maskloc = blast_mask_loc_new(num_contexts);

    let Some(mask_loc) = filter_maskloc.as_mut() else {
        return if num_contexts == 0 {
            0
        } else {
            crate::util::BLASTERR_INVALIDPARAM
        };
    };

    for context in 0..num_contexts {
        let mut filter_per_context = None;
        let status = s_get_filtering_locations_for_one_context(
            query_blk,
            query_info,
            context,
            program_number,
            filter_options,
            &mut filter_per_context,
            lcase_mask.as_deref_mut(),
        );
        if status != 0 {
            return status;
        }
        mask_loc.masks[context] = filter_per_context;
    }

    0
}

/// Port of NCBI `Blast_MaskUnsupportedAA` (`blast_filter.c:1335`).
pub fn blast_mask_unsupported_aa(seq: &mut [u8], min_invalid: u8) {
    for residue in seq {
        if *residue >= min_invalid {
            *residue = K_PROT_MASK;
        }
    }
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

const NUM_DIMERS: usize = 1 << 4;

/// Port of NCBI internal `s_FindDimerEntropy` (`jumper.c:4482`).
pub fn s_find_dimer_entropy(sequence: &[u8]) -> i32 {
    let mut counts = [0i32; NUM_DIMERS];
    let mut num = 0i32;
    let mut sum = 0.0;

    for pair in sequence.windows(2) {
        let base_1 = pair[0];
        let base_2 = pair[1];
        if (base_1 & 0xfc) == 0 && (base_2 & 0xfc) == 0 {
            let dimer = ((base_1 << 2) | base_2) as usize;
            counts[dimer] += 1;
            num += 1;
        }
    }

    if num > 0 {
        for &count in &counts {
            if count != 0 {
                sum += count as f64 * ((count as f64) / (num as f64)).ln();
            }
        }
    }

    (-sum * (1.0 / 16.0f64.ln()) + 0.5) as i32
}

fn s_mask_sequence(offset: i32, length: i32, seq_loc: &mut Option<SSeqRange>) -> i16 {
    *seq_loc = Some(SSeqRange {
        left: offset,
        right: offset + length - 1,
    });
    0
}

/// Port-shaped equivalent of NCBI `FilterQueriesForMapping` (`jumper.c:4531`).
pub fn filter_queries_for_mapping(
    sequence: &[u8],
    offset: i32,
    options: &SReadQualityOptions,
    seq_loc: &mut Option<SSeqRange>,
) -> i16 {
    let length = sequence.len() as i32;
    if length <= 0 {
        return 0;
    }

    let mut num = 0i32;
    for &base in sequence {
        if base & 0xfc != 0 {
            num += 1;
        }
    }

    if (num as f64) / (length as f64) > options.frac_ambig {
        return s_mask_sequence(offset, length, seq_loc);
    }

    let entropy = s_find_dimer_entropy(sequence);
    if entropy <= options.entropy {
        return s_mask_sequence(offset, length, seq_loc);
    }

    0
}

/// Port of NCBI `GetCutoffScore` (`jumper.c:4563`).
pub fn get_cutoff_score(query_length: i32) -> i32 {
    if query_length <= 20 {
        query_length
    } else if query_length <= 34 {
        20
    } else if query_length < 200 {
        (0.6 * query_length as f64) as i32
    } else {
        120
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
        let mut t = (crate::encoding::blastna_or_iupac_to_ncbi2na_base(sequence[start]) << 2)
            + crate::encoding::blastna_or_iupac_to_ncbi2na_base(sequence[start + 1]);
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

            t = ((t << 2) & TRIPLET_MASK)
                + (crate::encoding::blastna_or_iupac_to_ncbi2na_base(sequence[next_pos]) & 0x3);
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
                    t = ((t << 2) & TRIPLET_MASK)
                        + (crate::encoding::blastna_or_iupac_to_ncbi2na_base(sequence[next_pos])
                            & 0x3);
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
    -std::f64::consts::LN_10,
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

#[derive(Clone, Debug, PartialEq)]
struct SSequence {
    seq: Vec<u8>,
    parent_seq: Vec<u8>,
    length: usize,
    parent_length: usize,
    start: usize,
    left: usize,
    right: usize,
    bogus: usize,
    dash: bool,
    composition: Vec<i32>,
    state: Vec<i32>,
    entropy: f64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SSeg {
    begin: usize,
    end: usize,
    next: Option<Box<SSeg>>,
}

/// Port-shaped constructor for NCBI `s_SSequenceNew`.
#[allow(dead_code)]
fn s_ssequence_new() -> SSequence {
    SSequence {
        seq: Vec::new(),
        parent_seq: Vec::new(),
        length: 0,
        parent_length: 0,
        start: 0,
        left: 0,
        right: 0,
        bogus: 0,
        dash: false,
        composition: Vec::new(),
        state: Vec::new(),
        entropy: 0.0,
    }
}

/// Rust ownership equivalent of NCBI `s_SSequenceFree`.
#[allow(dead_code)]
fn s_ssequence_free(seq: Option<SSequence>) -> Option<SSequence> {
    if let Some(mut seq) = seq {
        seq.seq.clear();
        seq.composition.clear();
        seq.state.clear();
    }
    None
}

/// Rust ownership equivalent of NCBI `s_SegFree`.
#[allow(dead_code)]
fn s_seg_free(seg: Option<Box<SSeg>>) -> Option<Box<SSeg>> {
    let mut seg = seg;
    while let Some(mut current) = seg {
        seg = current.next.take();
    }
    None
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

/// Port-shaped constructor for NCBI `SegParametersNewAa`.
#[allow(dead_code)]
fn seg_parameters_new_aa() -> SegParameters {
    let window = SEG_DEFAULT_WINDOW;
    let locut = 2.2;
    let hicut = SEG_DEFAULT_HICUT;
    normalize_seg_params(window, locut, hicut)
}

/// Port-shaped validation/defaulting helper for NCBI `s_SegParametersCheck`.
#[allow(dead_code)]
fn s_seg_parameters_check(params: &mut SegParameters) {
    let window = params.window;
    let locut = params.locut;
    let hicut = params.hicut;
    let checked = normalize_seg_params(window, locut, hicut);
    *params = checked;
}

/// Rust ownership drops SEG parameters; this exists to keep the audit mapping
/// aligned with NCBI `SegParametersFree`.
#[allow(dead_code)]
fn seg_parameters_free(_params: Option<SegParameters>) {}

fn seg_alpha_index(code: u8) -> Option<usize> {
    SEG_VALID_AA_CODES.iter().position(|&valid| valid == code)
}

#[allow(dead_code)]
fn s_has_dash(seq: &[u8]) -> bool {
    seq.contains(&b'-')
}

#[allow(dead_code)]
fn s_state_cmp(a: i32, b: i32) -> std::cmp::Ordering {
    b.cmp(&a)
}

#[allow(dead_code)]
fn s_lnfact(value: i32) -> f64 {
    crate::math::ln_factorial(value)
}

#[allow(dead_code)]
fn s_aa20alpha_std() -> &'static [u8; 20] {
    let alphabet = &SEG_VALID_AA_CODES;
    debug_assert_eq!(alphabet.len(), 20);
    debug_assert_eq!(alphabet[0], 1);
    debug_assert_eq!(alphabet[19], 24);
    alphabet
}

fn s_comp_on(seq: &[u8]) -> ([i32; 20], usize) {
    let mut composition = [0i32; 20];
    let mut bogus = 0usize;
    for &aa in seq {
        if let Some(idx) = seg_alpha_index(aa) {
            composition[idx] += 1;
        } else {
            bogus += 1;
        }
    }
    (composition, bogus)
}

fn s_state_on(seq: &[u8]) -> (Vec<i32>, usize) {
    let (composition, bogus) = s_comp_on(seq);
    let mut state: Vec<i32> = composition.into_iter().filter(|&count| count > 0).collect();
    state.sort_unstable_by(|a, b| b.cmp(a));
    (state, bogus)
}

/// Port-shaped constructor for NCBI `s_OpenWin`.
#[allow(dead_code)]
fn s_open_win(parent: &SSequence, start: usize, length: usize) -> Option<SSequence> {
    if start.checked_add(length)? > parent.seq.len() {
        return None;
    }

    let seq = parent.seq[start..start + length].to_vec();
    let (composition, bogus) = s_comp_on(&seq);
    let (state, _) = s_state_on(&seq);
    Some(SSequence {
        parent_seq: parent.seq.clone(),
        parent_length: parent.seq.len(),
        start,
        length,
        left: start,
        right: start + length,
        dash: seq.contains(&b'-'),
        composition: composition.to_vec(),
        state,
        entropy: -2.0,
        bogus,
        seq,
    })
}

fn s_entropy(state: &[i32]) -> f64 {
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

/// Port of NCBI `s_DecrementSV`.
#[allow(dead_code)]
fn s_decrement_sv(sv: &mut Vec<i32>, class: i32) {
    for i in 0..sv.len() {
        if sv[i] == class && (i + 1 == sv.len() || sv[i + 1] < class) {
            sv[i] -= 1;
            break;
        }
    }
    sv.retain(|&count| count > 0);
    sv.sort_unstable_by(|a, b| b.cmp(a));
}

/// Port of NCBI `s_IncrementSV`.
#[allow(dead_code)]
fn s_increment_sv(sv: &mut Vec<i32>, class: i32) {
    if class <= 0 {
        sv.push(1);
    } else if let Some(count) = sv.iter_mut().find(|count| **count == class) {
        *count += 1;
    } else {
        sv.push(class + 1);
    }
    sv.sort_unstable_by(|a, b| b.cmp(a));
}

/// Port-shaped window shifter for NCBI `s_ShiftWin1`.
#[allow(dead_code)]
fn s_shift_win1(win: &mut SSequence) -> bool {
    if win.parent_seq.is_empty() || win.start + win.length >= win.parent_length {
        return false;
    }
    if win.seq.is_empty() || win.length == 0 {
        return false;
    }

    let outgoing = win.seq[0];
    let incoming = win.parent_seq[win.start + win.length];
    if let Some(idx) = seg_alpha_index(outgoing) {
        let class = win.composition.get(idx).copied().unwrap_or(0);
        if let Some(count) = win.composition.get_mut(idx) {
            *count -= 1;
        }
        s_decrement_sv(&mut win.state, class);
    } else {
        win.bogus = win.bogus.saturating_sub(1);
    }

    if let Some(idx) = seg_alpha_index(incoming) {
        if win.composition.len() < SEG_VALID_AA_CODES.len() {
            win.composition.resize(SEG_VALID_AA_CODES.len(), 0);
        }
        let class = win.composition[idx];
        win.composition[idx] += 1;
        s_increment_sv(&mut win.state, class);
    } else {
        win.bogus += 1;
    }

    win.start += 1;
    win.left = win.start;
    win.right = win.start + win.length;
    win.seq.remove(0);
    win.seq.push(incoming);
    win.dash = win.seq.contains(&b'-');
    if win.entropy > -2.0 {
        win.entropy = s_entropy(&win.state);
    }
    true
}

/// Rust ownership equivalent of NCBI `s_CloseWin`.
#[allow(dead_code)]
fn s_close_win(win: Option<SSequence>) -> Option<SSequence> {
    if let Some(mut win) = win {
        win.state.clear();
        win.composition.clear();
        win.seq.clear();
        win.parent_seq.clear();
    }
    None
}

/// Port of NCBI `s_EntropyOn`.
#[allow(dead_code)]
fn s_entropy_on(win: &mut SSequence) {
    if win.state.is_empty() {
        let (state, bogus) = s_state_on(&win.seq);
        win.state = state;
        win.bogus = bogus;
    }
    win.entropy = s_entropy(&win.state);
}

fn seg_window_entropy(seq: &[u8]) -> (f64, Vec<i32>, usize) {
    let (state, bogus) = s_state_on(seq);
    (s_entropy(&state), state, bogus)
}

fn s_seq_entropy(seq: &[u8], window: usize, maxbogus: usize) -> Option<Vec<f64>> {
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

fn s_find_low(i: usize, limit: usize, hicut: f64, entropies: &[f64]) -> usize {
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

fn s_find_high(i: usize, limit: usize, hicut: f64, entropies: &[f64]) -> usize {
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

fn s_ln_ass(state: &[i32], alphasize: usize) -> f64 {
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

fn s_get_prob(state: &[i32], total: usize) -> f64 {
    s_ln_ass(state, 20) + seg_ln_perm(state, total) - total as f64 * SEG_LN_20
}

fn s_trim(seq: &[u8], leftend: &mut usize, rightend: &mut usize, params: SegParameters) {
    let mut lend = 0usize;
    let mut rend = seq.len().saturating_sub(1);
    let minlen = (seq.len().saturating_sub(params.maxtrim)).max(1);
    let mut minprob = 1.0f64;

    for len in (minlen + 1..=seq.len()).rev() {
        for start in 0..=seq.len() - len {
            let end = start + len;
            let (state, _) = s_state_on(&seq[start..end]);
            let prob = s_get_prob(&state, len);
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

fn s_seg_seq(seq: &[u8], params: SegParameters, segs: &mut Vec<SegSegment>, offset: usize) {
    if params.window == 0 || seq.is_empty() || params.window > seq.len() {
        return;
    }

    let Some(entropies) = s_seq_entropy(seq, params.window, params.maxbogus) else {
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
            let loi = s_find_low(i, lowlim, params.hicut, &entropies);
            let hii = s_find_high(i, last, params.hicut, &entropies);

            let mut leftend = loi - downset;
            let mut rightend = hii + upset - 1;
            s_trim(
                &seq[leftend..=rightend],
                &mut leftend,
                &mut rightend,
                params,
            );

            if i + upset - 1 < leftend {
                let lend = loi - downset;
                let rend = leftend - 1;
                s_seg_seq(&seq[lend..=rend], params, segs, offset + lend);
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

fn s_merge_segs(seq_len: usize, segs: &mut Vec<SegSegment>) {
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
    s_seg_seq(sequence, params, &mut segs, 0);
    if params.overlaps {
        s_merge_segs(sequence.len(), &mut segs);
    }

    for seg in segs {
        if seg.begin <= seg.end {
            mask.add(seg.begin as i32, seg.end as i32 + 1);
        }
    }

    mask
}

/// Port-shaped Rust entry point for NCBI `SeqBufferSeg`; returns the same
/// half-open mask representation used by this crate.
#[allow(dead_code)]
fn seq_buffer_seg(
    sequence: &[u8],
    offset: i32,
    mut params: Option<SegParameters>,
) -> (i16, MaskLoc) {
    let params = params.get_or_insert_with(seg_parameters_new_aa);
    s_seg_parameters_check(params);
    let mut mask = seg_filter_ncbistdaa(sequence, params.window, params.locut, params.hicut);
    if offset != 0 {
        for region in &mut mask.regions {
            region.start += offset;
            region.end += offset;
        }
    }
    (0, mask)
}

#[allow(dead_code)]
fn s_segs_to_blast_seq_loc(segs: &[SegSegment], offset: i32) -> MaskLoc {
    let mut mask = MaskLoc::new();
    for seg in segs {
        if seg.begin <= seg.end {
            mask.add(seg.begin as i32 + offset, seg.end as i32 + offset + 1);
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
        &crate::encoding::encode_ncbistdaa_sequence(sequence),
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

    #[test]
    fn test_mapping_read_quality_helpers_mask_ambiguous_or_low_entropy_reads() {
        assert_eq!(s_find_dimer_entropy(&[0, 1, 2, 3, 0, 1]), 2);
        assert_eq!(s_find_dimer_entropy(&[0, 0, 0, 0, 0]), 0);

        let options = SReadQualityOptions {
            frac_ambig: 0.25,
            entropy: 0,
        };

        let mut loc = None;
        assert_eq!(
            filter_queries_for_mapping(&[0, 1, 0xfc, 0xfd], 7, &options, &mut loc),
            0
        );
        assert_eq!(loc, Some(SSeqRange { left: 7, right: 10 }));

        let mut low_entropy = None;
        assert_eq!(
            filter_queries_for_mapping(&[0, 0, 0, 0], 3, &options, &mut low_entropy),
            0
        );
        assert_eq!(low_entropy, Some(SSeqRange { left: 3, right: 6 }));

        let mut kept = None;
        assert_eq!(
            filter_queries_for_mapping(&[0, 1, 2, 3, 0, 1], 0, &options, &mut kept),
            0
        );
        assert!(kept.is_none());

        assert_eq!(get_cutoff_score(20), 20);
        assert_eq!(get_cutoff_score(34), 20);
        assert_eq!(get_cutoff_score(100), 60);
        assert_eq!(get_cutoff_score(200), 120);
    }

    #[test]
    fn translated_filter_option_lifecycle_matches_c_defaults() {
        let mut dust = None;
        assert_eq!(sdust_options_new(&mut dust), 0);
        assert_eq!(
            dust,
            Some(SDustOptions {
                level: K_DUST_LEVEL,
                window: K_DUST_WINDOW,
                linker: K_DUST_LINKER,
            })
        );
        assert!(sdust_options_free(&mut dust).is_none());
        assert!(dust.is_none());

        let mut seg = None;
        assert_eq!(sseg_options_new(&mut seg), 0);
        assert_eq!(
            seg,
            Some(SSegOptions {
                window: K_SEG_WINDOW,
                locut: K_SEG_LOCUT,
                hicut: K_SEG_HICUT,
            })
        );
        assert!(sseg_options_free(&mut seg).is_none());
        assert!(seg.is_none());

        let mut read_quality = None;
        assert_eq!(sread_quality_options_new(&mut read_quality), 0);
        assert_eq!(
            read_quality,
            Some(SReadQualityOptions {
                frac_ambig: 0.5,
                entropy: 16,
            })
        );
        assert!(sread_quality_options_free(&mut read_quality).is_none());
        assert!(read_quality.is_none());
    }

    #[test]
    fn translated_repeat_and_window_masker_reset_db() {
        let mut repeat = None;
        assert_eq!(srepeat_filter_options_new(&mut repeat), 0);
        assert_eq!(
            repeat
                .as_ref()
                .and_then(|options| options.database.as_deref()),
            Some(K_DEFAULT_REPEAT_FILTER_DB)
        );
        assert_eq!(
            srepeat_filter_options_reset_db(&mut repeat, "repeat/custom"),
            0
        );
        assert_eq!(
            repeat
                .as_ref()
                .and_then(|options| options.database.as_deref()),
            Some("repeat/custom")
        );
        assert!(srepeat_filter_options_free(&mut repeat).is_none());
        assert!(repeat.is_none());

        let mut winmask = None;
        assert_eq!(swindow_masker_options_new(&mut winmask), 0);
        assert_eq!(
            winmask.as_ref(),
            Some(&SWindowMaskerOptions {
                taxid: 0,
                database: None,
            })
        );
        assert_eq!(
            swindow_masker_options_reset_db(&mut winmask, Some("window/db")),
            0
        );
        assert_eq!(
            winmask
                .as_ref()
                .and_then(|options| options.database.as_deref()),
            Some("window/db")
        );
        assert_eq!(swindow_masker_options_reset_db(&mut winmask, None), 0);
        assert_eq!(
            winmask
                .as_ref()
                .and_then(|options| options.database.as_deref()),
            None
        );
        assert!(swindow_masker_options_free(&mut winmask).is_none());
        assert!(winmask.is_none());
    }

    #[test]
    fn translated_blast_filter_options_new_and_free() {
        let mut seg_filter = None;
        assert_eq!(
            sblast_filter_options_new(&mut seg_filter, EFilterOptions::ESeg),
            0
        );
        let seg_filter_ref = seg_filter.as_ref().expect("seg filter");
        assert!(!seg_filter_ref.mask_at_hash);
        assert!(seg_filter_ref.seg_options.is_some());
        assert!(seg_filter_ref.dust_options.is_none());
        assert!(seg_filter_ref.repeat_filter_options.is_none());

        let mut dust_repeats_filter = None;
        assert_eq!(
            sblast_filter_options_new(&mut dust_repeats_filter, EFilterOptions::EDustRepeats),
            0
        );
        let dust_repeats_ref = dust_repeats_filter.as_ref().expect("dust repeats filter");
        assert!(dust_repeats_ref.dust_options.is_some());
        assert!(dust_repeats_ref.repeat_filter_options.is_some());
        assert!(dust_repeats_ref.seg_options.is_none());

        assert!(sblast_filter_options_free(&mut dust_repeats_filter).is_none());
        assert!(dust_repeats_filter.is_none());
        assert!(sblast_filter_options_free(&mut seg_filter).is_none());
        assert!(seg_filter.is_none());
    }

    #[test]
    fn translated_filter_option_merge_helpers_prefer_non_defaults() {
        let default_dust = SDustOptions {
            level: K_DUST_LEVEL,
            window: K_DUST_WINDOW,
            linker: K_DUST_LINKER,
        };
        let custom_dust = SDustOptions {
            level: 15,
            window: K_DUST_WINDOW,
            linker: 4,
        };
        assert_eq!(
            s_merge_dust_options(Some(&custom_dust), Some(&default_dust)),
            Some(SDustOptions {
                level: 15,
                window: K_DUST_WINDOW,
                linker: 4,
            })
        );

        let default_seg = SSegOptions {
            window: K_SEG_WINDOW,
            locut: K_SEG_LOCUT,
            hicut: K_SEG_HICUT,
        };
        let custom_seg = SSegOptions {
            window: 20,
            locut: K_SEG_LOCUT,
            hicut: 3.0,
        };
        assert_eq!(
            s_merge_seg_options(Some(&default_seg), Some(&custom_seg)),
            Some(custom_seg)
        );

        let first_repeat = SRepeatFilterOptions {
            database: Some("repeat/first".to_string()),
        };
        let second_repeat = SRepeatFilterOptions {
            database: Some("repeat/second".to_string()),
        };
        assert_eq!(
            s_merge_repeat_options(Some(&first_repeat), Some(&second_repeat)),
            Some(second_repeat.clone())
        );

        let empty_winmask = SWindowMaskerOptions {
            taxid: 0,
            database: None,
        };
        let taxid_winmask = SWindowMaskerOptions {
            taxid: 9606,
            database: Some("window/human".to_string()),
        };
        assert!(s_merge_window_masker_options(Some(&empty_winmask), None).is_none());
        assert_eq!(
            s_merge_window_masker_options(Some(&empty_winmask), Some(&taxid_winmask)),
            Some(taxid_winmask)
        );
    }

    #[test]
    fn translated_blast_filter_options_merge_and_query_helpers() {
        let first = SBlastFilterOptions {
            mask_at_hash: true,
            dust_options: Some(SDustOptions {
                level: 15,
                window: K_DUST_WINDOW,
                linker: K_DUST_LINKER,
            }),
            seg_options: None,
            repeat_filter_options: None,
            window_masker_options: None,
            read_quality_options: Some(SReadQualityOptions {
                frac_ambig: 0.5,
                entropy: 16,
            }),
        };
        let second = SBlastFilterOptions {
            mask_at_hash: false,
            dust_options: Some(SDustOptions {
                level: K_DUST_LEVEL,
                window: 32,
                linker: 2,
            }),
            seg_options: Some(SSegOptions {
                window: K_SEG_WINDOW,
                locut: 1.9,
                hicut: K_SEG_HICUT,
            }),
            repeat_filter_options: Some(SRepeatFilterOptions {
                database: Some("repeat/second".to_string()),
            }),
            window_masker_options: Some(SWindowMaskerOptions {
                taxid: 9606,
                database: None,
            }),
            read_quality_options: Some(SReadQualityOptions {
                frac_ambig: 0.25,
                entropy: 12,
            }),
        };

        let mut combined = None;
        assert_eq!(
            sblast_filter_options_merge(&mut combined, Some(&first), Some(&second)),
            0
        );
        let combined_ref = combined.as_ref().expect("combined filter");
        assert!(sblast_filter_options_mask_at_hash(Some(combined_ref)));
        assert!(!sblast_filter_options_no_filtering(Some(combined_ref)));
        assert_eq!(
            combined_ref.dust_options,
            Some(SDustOptions {
                level: 15,
                window: 32,
                linker: 2,
            })
        );
        assert_eq!(combined_ref.seg_options, second.seg_options);
        assert_eq!(
            combined_ref
                .repeat_filter_options
                .as_ref()
                .and_then(|options| options.database.as_deref()),
            Some("repeat/second")
        );
        assert_eq!(
            combined_ref.window_masker_options,
            second.window_masker_options
        );
        assert_eq!(
            combined_ref.read_quality_options,
            second.read_quality_options
        );

        let mut empty_merge = Some(first);
        assert_eq!(sblast_filter_options_merge(&mut empty_merge, None, None), 0);
        assert!(empty_merge.is_none());
        assert!(sblast_filter_options_no_filtering(None));
        assert!(!sblast_filter_options_mask_at_hash(None));

        let read_quality_only = SBlastFilterOptions {
            mask_at_hash: false,
            dust_options: None,
            seg_options: None,
            repeat_filter_options: None,
            window_masker_options: None,
            read_quality_options: Some(SReadQualityOptions {
                frac_ambig: 0.5,
                entropy: 16,
            }),
        };
        assert!(!sblast_filter_options_no_filtering(Some(
            &read_quality_only
        )));
    }

    #[test]
    fn translated_blast_filter_options_validate_matches_c_rules() {
        let mut dust_repeats = None;
        assert_eq!(
            sblast_filter_options_new(&mut dust_repeats, EFilterOptions::EDustRepeats),
            0
        );
        let mut dust_repeats = dust_repeats.expect("dust/repeats filter");
        assert_eq!(
            sblast_filter_options_validate(crate::program::BLASTN, Some(&dust_repeats)),
            0
        );
        assert_eq!(
            sblast_filter_options_validate(crate::program::MAPPING, Some(&dust_repeats)),
            0
        );
        assert_eq!(
            sblast_filter_options_validate(crate::program::BLASTP, Some(&dust_repeats)),
            crate::options::BLASTERR_OPTION_PROGRAM_INVALID
        );

        dust_repeats
            .repeat_filter_options
            .as_mut()
            .expect("repeat options")
            .database = Some(String::new());
        assert_eq!(
            sblast_filter_options_validate(crate::program::BLASTN, Some(&dust_repeats)),
            crate::util::BLASTERR_INVALIDPARAM
        );

        let mut seg = None;
        assert_eq!(sblast_filter_options_new(&mut seg, EFilterOptions::ESeg), 0);
        let seg = seg.expect("seg filter");
        assert_eq!(
            sblast_filter_options_validate(crate::program::BLASTP, Some(&seg)),
            0
        );
        assert_eq!(
            sblast_filter_options_validate(crate::program::BLASTN, Some(&seg)),
            crate::options::BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert_eq!(
            sblast_filter_options_validate(crate::program::MAPPING, Some(&seg)),
            crate::options::BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert_eq!(
            sblast_filter_options_validate(crate::program::BLASTN, None),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_filtering_options_string_helpers_round_trip() {
        let mut options = None;
        assert_eq!(
            blast_filtering_options_from_string(
                crate::program::BLASTN,
                Some("D 22 70 2 R repeat/custom m windowmasker_taxid 9606"),
                &mut options,
            ),
            0
        );
        let options_ref = options.as_ref().expect("nucleotide filters");
        assert_eq!(
            options_ref.dust_options,
            Some(SDustOptions {
                level: 22,
                window: 70,
                linker: 2,
            })
        );
        assert_eq!(
            options_ref
                .repeat_filter_options
                .as_ref()
                .and_then(|repeat| repeat.database.as_deref()),
            Some("repeat/custom")
        );
        assert!(options_ref.mask_at_hash);
        assert_eq!(
            options_ref.window_masker_options,
            Some(SWindowMaskerOptions {
                taxid: 9606,
                database: None,
            })
        );

        let serialized =
            blast_filtering_options_to_string(crate::program::BLASTN, options.as_ref())
                .expect("filter string");
        assert!(serialized.contains("D 22 70 2"));
        assert!(serialized.contains("R repeat/custom"));
        assert!(serialized.contains("windowmasker_taxid 9606"));
        assert!(serialized.contains('m'));

        let mut round_tripped = None;
        assert_eq!(
            blast_filtering_options_from_string(
                crate::program::BLASTN,
                Some(&serialized),
                &mut round_tripped,
            ),
            0
        );
        assert_eq!(round_tripped, options);

        let mut protein = None;
        assert_eq!(
            blast_filtering_options_from_string(
                crate::program::BLASTP,
                Some("seg 14 2.1 2.8 read_quality 0.25 12"),
                &mut protein,
            ),
            0
        );
        let protein_ref = protein.as_ref().expect("protein filters");
        assert_eq!(
            protein_ref.seg_options,
            Some(SSegOptions {
                window: 14,
                locut: 2.1,
                hicut: 2.8,
            })
        );
        assert_eq!(
            protein_ref.read_quality_options,
            Some(SReadQualityOptions {
                frac_ambig: 0.25,
                entropy: 12,
            })
        );

        let mut invalid = None;
        assert_eq!(
            blast_filtering_options_from_string(crate::program::BLASTP, Some("dust"), &mut invalid),
            crate::options::BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert_eq!(
            blast_filtering_options_from_string(
                crate::program::BLASTN,
                Some("unknown"),
                &mut invalid
            ),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_blast_seq_loc_and_mask_loc_lifecycle() {
        let mut list = None;
        assert_eq!(blast_seq_loc_new(&mut list, 4, 9), 0);
        assert_eq!(
            list.as_ref().map(|node| node.ssl),
            Some(SSeqRange { left: 4, right: 9 })
        );

        let mut second = None;
        assert_eq!(blast_seq_loc_new(&mut second, 20, 25), 0);
        assert_eq!(blast_seq_loc_new(&mut second, 30, 35), 0);
        assert_eq!(blast_seq_loc_append(&mut list, second), 0);
        assert_eq!(blast_seq_loc_new(&mut list, 40, 45), 0);
        assert_eq!(s_blast_seq_loc_len(list.as_deref()), 4);

        let array = s_blast_seq_loc_list_to_array_of_pointers(list.as_deref());
        assert_eq!(array.len(), 4);
        assert_eq!(array[0].ssl, SSeqRange { left: 4, right: 9 });
        assert_eq!(
            array[1].ssl,
            SSeqRange {
                left: 20,
                right: 25,
            }
        );
        assert_eq!(
            array[2].ssl,
            SSeqRange {
                left: 30,
                right: 35,
            }
        );
        assert_eq!(
            array[3].ssl,
            SSeqRange {
                left: 40,
                right: 45,
            }
        );

        let duplicated_node = s_blast_seq_loc_node_dup(list.as_deref()).expect("duplicated node");
        assert_eq!(duplicated_node.ssl, SSeqRange { left: 4, right: 9 });
        assert!(duplicated_node.next.is_none());

        let duplicated_list = blast_seq_loc_list_dup(list.as_deref()).expect("duplicated list");
        assert_eq!(s_blast_seq_loc_len(Some(&duplicated_list)), 4);
        assert_eq!(duplicated_list.ssl, SSeqRange { left: 4, right: 9 });
        assert_eq!(
            duplicated_list.next.as_ref().map(|node| node.ssl),
            Some(SSeqRange {
                left: 20,
                right: 25
            })
        );
        assert_ne!(
            duplicated_list
                .next
                .as_ref()
                .map(|node| &**node as *const BlastSeqLoc),
            list.as_ref()
                .and_then(|node| node.next.as_ref())
                .map(|node| &**node as *const BlastSeqLoc)
        );

        assert_eq!(blast_seq_loc_list_reverse(&mut list), 0);
        assert_eq!(
            list.as_ref().map(|node| node.ssl),
            Some(SSeqRange {
                left: 40,
                right: 45,
            })
        );

        let mut mask_loc = blast_mask_loc_new(2);
        let mask = mask_loc.as_mut().expect("mask location");
        mask.masks[0] = list;
        let duplicate = blast_mask_loc_dup(Some(mask)).expect("duplicated mask location");
        assert_eq!(duplicate.masks.len(), 2);
        assert_eq!(s_blast_seq_loc_len(duplicate.masks[0].as_deref()), 4);
        assert!(duplicate.masks[1].is_none());

        assert!(blast_mask_loc_new(0).is_none());
        assert_eq!(blast_mask_loc_free(&mut mask_loc), None);
        assert!(mask_loc.is_none());

        let mut single = Some(duplicated_node);
        assert_eq!(blast_seq_loc_node_free(&mut single), None);
        assert!(single.is_none());
    }

    #[test]
    fn translated_blast_seq_loc_restrict_to_interval_matches_c_rewrite() {
        let mut list = None;
        assert_eq!(blast_seq_loc_new(&mut list, 4, 9), 0);
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange {
                        left: 12,
                        right: 18
                    },
                    next: None
                }))
            ),
            0
        );
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange {
                        left: 25,
                        right: 30
                    },
                    next: None
                }))
            ),
            0
        );

        blast_seq_loc_restrict_to_interval(Some(&mut list), 10, 20);
        let ranges: Vec<SSeqRange> = s_blast_seq_loc_list_to_array_of_pointers(list.as_deref())
            .iter()
            .map(|node| node.ssl)
            .collect();
        assert_eq!(ranges, vec![SSeqRange { left: 2, right: 8 }]);

        let before = list.clone();
        blast_seq_loc_restrict_to_interval(Some(&mut list), 0, 0);
        assert_eq!(list, before);

        blast_seq_loc_restrict_to_interval(Some(&mut list), 100, 110);
        assert!(list.is_none());

        blast_seq_loc_restrict_to_interval(None, 0, 1);
    }

    #[test]
    fn translated_blast_seq_loc_combine_and_mask_unsupported_aa() {
        let mut list = None;
        assert_eq!(blast_seq_loc_new(&mut list, 30, 35), 0);
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange { left: 5, right: 10 },
                    next: None,
                })),
            ),
            0
        );
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange {
                        left: 12,
                        right: 14,
                    },
                    next: None,
                })),
            ),
            0
        );
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange {
                        left: 40,
                        right: 42,
                    },
                    next: None,
                })),
            ),
            0
        );

        blast_seq_loc_combine(&mut list, 1);
        let combined = s_blast_seq_loc_list_to_array_of_pointers(list.as_deref());
        let ranges: Vec<SSeqRange> = combined.iter().map(|node| node.ssl).collect();
        assert_eq!(
            ranges,
            vec![
                SSeqRange { left: 5, right: 14 },
                SSeqRange {
                    left: 30,
                    right: 35,
                },
                SSeqRange {
                    left: 40,
                    right: 42,
                },
            ]
        );

        let mut protein = vec![
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_U,
            crate::encoding::NCBISTDAA_STOP,
            crate::encoding::NCBISTDAA_J,
        ];
        blast_mask_unsupported_aa(&mut protein, crate::encoding::NCBISTDAA_U);
        assert_eq!(
            protein,
            vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_X,
            ]
        );
    }

    #[test]
    fn translated_blast_mask_the_residues_uses_c_offsets_and_letters() {
        let mut list = None;
        assert_eq!(blast_seq_loc_new(&mut list, 3, 5), 0);
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange { left: 8, right: 9 },
                    next: None,
                })),
            ),
            0
        );

        let mut protein = vec![1u8; 12];
        blast_mask_the_residues(&mut protein, false, list.as_deref(), false, 2);
        assert_eq!(&protein[1..=3], &[K_PROT_MASK; 3]);
        assert_eq!(&protein[6..=7], &[K_PROT_MASK; 2]);
        assert_eq!(protein[0], 1);
        assert_eq!(protein[8], 1);

        let mut nucleotide = vec![0u8; 10];
        blast_mask_the_residues(&mut nucleotide, true, list.as_deref(), true, 0);
        assert_eq!(&nucleotide[4..=6], &[K_NUCL_MASK; 3]);
        assert_eq!(&nucleotide[0..=1], &[K_NUCL_MASK; 2]);
        assert_eq!(nucleotide[2], 0);
        assert_eq!(nucleotide[7], 0);
    }

    #[test]
    fn translated_blast_setup_filter_runs_seg_and_read_quality_paths() {
        let seg_options = SBlastFilterOptions {
            mask_at_hash: false,
            dust_options: None,
            seg_options: Some(SSegOptions {
                window: 12,
                locut: 2.2,
                hicut: 2.5,
            }),
            repeat_filter_options: None,
            window_masker_options: None,
            read_quality_options: None,
        };
        let low_complexity = vec![crate::encoding::NCBISTDAA_A; 32];
        let mut seg_locs = None;
        assert_eq!(
            blast_setup_filter(
                crate::program::BLASTP,
                &low_complexity,
                5,
                Some(&seg_options),
                &mut seg_locs,
            ),
            0
        );
        let seg_ranges: Vec<SSeqRange> =
            s_blast_seq_loc_list_to_array_of_pointers(seg_locs.as_deref())
                .iter()
                .map(|node| node.ssl)
                .collect();
        assert!(!seg_ranges.is_empty());
        assert!(seg_ranges[0].left >= 5);
        assert!(seg_ranges[0].right >= seg_ranges[0].left);

        let read_quality_options = SBlastFilterOptions {
            mask_at_hash: false,
            dust_options: None,
            seg_options: None,
            repeat_filter_options: None,
            window_masker_options: None,
            read_quality_options: Some(SReadQualityOptions {
                frac_ambig: 0.25,
                entropy: 16,
            }),
        };
        let mut read_quality_locs = None;
        assert_eq!(
            blast_setup_filter(
                crate::program::MAPPING,
                &[0xfc, 0xfd, 0, 1],
                7,
                Some(&read_quality_options),
                &mut read_quality_locs,
            ),
            0
        );
        assert_eq!(
            read_quality_locs.as_ref().map(|node| node.ssl),
            Some(SSeqRange { left: 7, right: 10 })
        );

        let mut invalid = Some(Box::new(BlastSeqLoc {
            ssl: SSeqRange { left: 1, right: 2 },
            next: None,
        }));
        assert_eq!(
            blast_setup_filter(
                crate::program::BLASTN,
                b"AAAA",
                0,
                Some(&seg_options),
                &mut invalid
            ),
            crate::options::BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert!(invalid.is_none());
    }

    #[test]
    fn translated_get_filtering_locations_one_context_combines_lcase_masks() {
        let query_info = QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 2,
                query_length: 6,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 6,
        };
        let mut query_blk = BlastSequenceBlk {
            sequence: Some(vec![99, 99, 0, 0, 0, 0, 0, 0, 88]),
            ..BlastSequenceBlk::default()
        };
        let options = SBlastFilterOptions {
            mask_at_hash: false,
            dust_options: None,
            seg_options: None,
            repeat_filter_options: None,
            window_masker_options: None,
            read_quality_options: Some(SReadQualityOptions {
                frac_ambig: 1.0,
                entropy: 16,
            }),
        };
        let mut lcase_mask = blast_mask_loc_new(1).expect("lcase mask");
        assert_eq!(blast_seq_loc_new(&mut lcase_mask.masks[0], 1, 2), 0);
        let mut filter_out = None;

        assert_eq!(
            s_get_filtering_locations_for_one_context(
                &mut query_blk,
                &query_info,
                0,
                crate::program::BLASTP,
                Some(&options),
                &mut filter_out,
                Some(&mut lcase_mask),
            ),
            0
        );

        let ranges: Vec<SSeqRange> =
            s_blast_seq_loc_list_to_array_of_pointers(filter_out.as_deref())
                .iter()
                .map(|node| node.ssl)
                .collect();
        assert_eq!(ranges, vec![SSeqRange { left: 0, right: 5 }]);
        assert!(lcase_mask.masks[0].is_none());

        let mut invalid_out = Some(Box::new(BlastSeqLoc {
            ssl: SSeqRange { left: 9, right: 10 },
            next: None,
        }));
        assert_eq!(
            s_get_filtering_locations_for_one_context(
                &mut query_blk,
                &query_info,
                9,
                crate::program::BLASTP,
                Some(&options),
                &mut invalid_out,
                None,
            ),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_blast_setup_get_filtering_locations_fills_context_table() {
        let query_info = QueryInfo {
            num_queries: 1,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: 6,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 7,
                    query_length: 6,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: -1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 14,
                    query_length: 4,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: false,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 6,
        };
        let mut query_blk = BlastSequenceBlk {
            sequence: Some(vec![
                0, 0, 0, 0, 0, 0, 99, 0xfc, 0xfd, 0, 1, 2, 3, 99, 5, 5, 5, 5,
            ]),
            ..BlastSequenceBlk::default()
        };
        let options = SBlastFilterOptions {
            mask_at_hash: false,
            dust_options: None,
            seg_options: None,
            repeat_filter_options: None,
            window_masker_options: None,
            read_quality_options: Some(SReadQualityOptions {
                frac_ambig: 0.25,
                entropy: 16,
            }),
        };
        let mut lcase_mask = blast_mask_loc_new(3).expect("lcase mask");
        assert_eq!(blast_seq_loc_new(&mut lcase_mask.masks[1], 2, 3), 0);

        let mut filter_maskloc = None;
        assert_eq!(
            blast_setup_get_filtering_locations(
                &mut query_blk,
                &query_info,
                crate::program::BLASTN,
                Some(&options),
                &mut filter_maskloc,
                Some(&mut lcase_mask),
            ),
            0
        );
        let mask_loc = filter_maskloc.expect("filter mask loc");
        assert_eq!(mask_loc.masks.len(), 3);
        assert_eq!(
            mask_loc.masks[0].as_ref().map(|node| node.ssl),
            Some(SSeqRange { left: 0, right: 5 })
        );
        assert_eq!(
            mask_loc.masks[1].as_ref().map(|node| node.ssl),
            Some(SSeqRange { left: 0, right: 5 })
        );
        assert!(mask_loc.masks[2].is_none());
        assert!(lcase_mask.masks[1].is_none());
    }

    #[test]
    fn translated_blast_mask_loc_dna_protein_coordinate_conversions() {
        let contexts: Vec<crate::queryinfo::ContextInfo> = (0..crate::util::NUM_FRAMES)
            .map(|context| crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 10,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: crate::util::blast_context_to_frame_blastx(context as u32),
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            })
            .collect();
        let query_info = QueryInfo {
            num_queries: 1,
            contexts,
            max_length: 10,
        };

        let mut dna_masks = blast_mask_loc_new(crate::util::NUM_FRAMES).expect("mask loc");
        assert_eq!(blast_seq_loc_new(&mut dna_masks.masks[0], 0, 8), 0);

        assert_eq!(
            blast_mask_loc_dna_to_protein(Some(&mut dna_masks), &query_info, &[30]),
            0
        );
        let protein_ranges: Vec<Option<SSeqRange>> = dna_masks
            .masks
            .iter()
            .map(|mask| mask.as_ref().map(|node| node.ssl))
            .collect();
        assert_eq!(
            protein_ranges,
            vec![
                Some(SSeqRange { left: 0, right: 2 }),
                Some(SSeqRange { left: 0, right: 2 }),
                Some(SSeqRange { left: 0, right: 2 }),
                Some(SSeqRange { left: 7, right: 9 }),
                Some(SSeqRange { left: 6, right: 9 }),
                Some(SSeqRange { left: 6, right: 9 }),
            ]
        );

        let mut protein_masks = blast_mask_loc_new(crate::util::NUM_FRAMES).expect("mask loc");
        assert_eq!(blast_seq_loc_new(&mut protein_masks.masks[0], 0, 2), 0);
        assert_eq!(blast_seq_loc_new(&mut protein_masks.masks[3], 7, 9), 0);

        assert_eq!(
            blast_mask_loc_protein_to_dna(Some(&mut protein_masks), &query_info, &[30]),
            0
        );
        assert_eq!(
            protein_masks.masks[0].as_ref().map(|node| node.ssl),
            Some(SSeqRange { left: 0, right: 8 })
        );
        assert_eq!(
            protein_masks.masks[3].as_ref().map(|node| node.ssl),
            Some(SSeqRange { left: 0, right: 8 })
        );

        assert_eq!(blast_mask_loc_dna_to_protein(None, &query_info, &[30]), 0);
        assert_eq!(blast_mask_loc_protein_to_dna(None, &query_info, &[30]), 0);
    }

    #[test]
    fn translated_blast_complement_mask_locations_matches_context_offsets() {
        let query_info = QueryInfo {
            num_queries: 2,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: 10,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 11,
                    query_length: 10,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: -1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 22,
                    query_length: 5,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 1,
                    is_valid: false,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 10,
        };

        let mut mask_loc = blast_mask_loc_new(3).expect("mask loc");
        assert_eq!(blast_seq_loc_new(&mut mask_loc.masks[0], 2, 4), 0);
        assert_eq!(
            blast_seq_loc_append(
                &mut mask_loc.masks[0],
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange { left: 7, right: 8 },
                    next: None,
                })),
            ),
            0
        );
        assert_eq!(blast_seq_loc_new(&mut mask_loc.masks[1], 1, 3), 0);
        assert_eq!(
            blast_seq_loc_append(
                &mut mask_loc.masks[1],
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange { left: 6, right: 7 },
                    next: None,
                })),
            ),
            0
        );

        let mut complement = None;
        assert_eq!(
            blast_complement_mask_locations(
                crate::program::BLASTN,
                &query_info,
                Some(&mask_loc),
                Some(&mut complement),
            ),
            0
        );
        let ranges: Vec<SSeqRange> =
            s_blast_seq_loc_list_to_array_of_pointers(complement.as_deref())
                .iter()
                .map(|node| node.ssl)
                .collect();
        assert_eq!(
            ranges,
            vec![
                SSeqRange { left: 0, right: 1 },
                SSeqRange { left: 5, right: 6 },
                SSeqRange { left: 9, right: 9 },
                SSeqRange {
                    left: 11,
                    right: 12
                },
                SSeqRange {
                    left: 15,
                    right: 16
                },
                SSeqRange {
                    left: 20,
                    right: 20
                },
            ]
        );

        let mut unmasked = None;
        assert_eq!(
            blast_complement_mask_locations(
                crate::program::BLASTP,
                &query_info,
                None,
                Some(&mut unmasked),
            ),
            0
        );
        let ranges: Vec<SSeqRange> = s_blast_seq_loc_list_to_array_of_pointers(unmasked.as_deref())
            .iter()
            .map(|node| node.ssl)
            .collect();
        assert_eq!(
            ranges,
            vec![
                SSeqRange { left: 0, right: 9 },
                SSeqRange {
                    left: 11,
                    right: 20
                },
            ]
        );

        assert_eq!(
            blast_complement_mask_locations(crate::program::BLASTN, &query_info, None, None),
            -1
        );
        assert!(blast_is_reverse_strand(true, 1));
        assert!(!blast_is_reverse_strand(true, 0));
        assert!(!blast_is_reverse_strand(false, 1));
    }

    #[test]
    fn translated_blast_setup_mask_query_copies_original_and_masks_contexts() {
        let query_info = QueryInfo {
            num_queries: 1,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 1,
                    query_length: 6,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 8,
                    query_length: 6,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: -1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 6,
        };
        let original: Vec<u8> = (0..16).collect();
        let mut query_blk = BlastSequenceBlk {
            sequence: Some(original.clone()),
            sequence_start: Some(original.clone()),
            ..BlastSequenceBlk::default()
        };
        let mut mask_loc = blast_mask_loc_new(2).expect("mask loc");
        assert_eq!(blast_seq_loc_new(&mut mask_loc.masks[0], 1, 2), 0);
        assert_eq!(blast_seq_loc_new(&mut mask_loc.masks[1], 0, 1), 0);

        blast_setup_mask_query(
            &mut query_blk,
            &query_info,
            &mask_loc,
            crate::program::BLASTN,
        );

        let sequence = query_blk.sequence.as_ref().expect("masked sequence");
        assert_eq!(&sequence[2..=3], &[K_NUCL_MASK; 2]);
        assert_eq!(&sequence[12..=13], &[K_NUCL_MASK; 2]);
        assert_eq!(sequence[1], original[1]);
        assert_eq!(sequence[8], original[8]);
        assert_eq!(
            query_blk.sequence_start_nomask.as_deref(),
            Some(&original[..16])
        );
        assert_eq!(query_blk.sequence_nomask.as_deref(), Some(&original[1..]));
        assert!(query_blk.nomask_allocated);

        let mut empty_blk = BlastSequenceBlk {
            sequence: Some(original.clone()),
            sequence_start: Some(original.clone()),
            ..BlastSequenceBlk::default()
        };
        let empty_mask = blast_mask_loc_new(2).expect("mask loc");
        blast_setup_mask_query(
            &mut empty_blk,
            &query_info,
            &empty_mask,
            crate::program::BLASTN,
        );
        assert_eq!(empty_blk.sequence.as_ref(), Some(&original));
        assert!(empty_blk.sequence_start_nomask.is_none());
        assert!(!empty_blk.nomask_allocated);
    }

    #[test]
    fn translated_blast_seq_loc_reverse_matches_c_coordinate_flip() {
        let mut list = None;
        assert_eq!(blast_seq_loc_new(&mut list, 4, 9), 0);
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange {
                        left: 20,
                        right: 25
                    },
                    next: None
                }))
            ),
            0
        );

        blast_seq_loc_reverse(list.as_deref_mut(), 40);
        let ranges: Vec<SSeqRange> = s_blast_seq_loc_list_to_array_of_pointers(list.as_deref())
            .iter()
            .map(|node| node.ssl)
            .collect();
        assert_eq!(
            ranges,
            vec![
                SSeqRange {
                    left: 30,
                    right: 35
                },
                SSeqRange {
                    left: 14,
                    right: 19
                },
            ]
        );
    }

    #[test]
    fn translated_seq_loc_list_invert_matches_lookup_helper() {
        let mut list = None;
        assert_eq!(blast_seq_loc_new(&mut list, 5, 9), 0);
        assert_eq!(
            blast_seq_loc_append(
                &mut list,
                Some(Box::new(BlastSeqLoc {
                    ssl: SSeqRange {
                        left: 15,
                        right: 20,
                    },
                    next: None,
                })),
            ),
            0
        );

        let inverted = s_seq_loc_list_invert(list.as_deref(), 30);
        let ranges: Vec<SSeqRange> = s_blast_seq_loc_list_to_array_of_pointers(inverted.as_deref())
            .iter()
            .map(|node| node.ssl)
            .collect();
        assert_eq!(
            ranges,
            vec![
                SSeqRange { left: 0, right: 4 },
                SSeqRange {
                    left: 10,
                    right: 14,
                },
                SSeqRange {
                    left: 21,
                    right: 29,
                },
            ]
        );

        let mut starts_near_front = None;
        assert_eq!(blast_seq_loc_new(&mut starts_near_front, 2, 6), 0);
        let inverted = s_seq_loc_list_invert(starts_near_front.as_deref(), 10);
        let ranges: Vec<SSeqRange> = s_blast_seq_loc_list_to_array_of_pointers(inverted.as_deref())
            .iter()
            .map(|node| node.ssl)
            .collect();
        assert_eq!(ranges, vec![SSeqRange { left: 7, right: 9 }]);

        let mut full_coverage = None;
        assert_eq!(blast_seq_loc_new(&mut full_coverage, 0, 9), 0);
        assert!(s_seq_loc_list_invert(full_coverage.as_deref(), 10).is_none());

        let mut adjacent = None;
        assert_eq!(blast_seq_loc_new(&mut adjacent, 0, 4), 0);
        assert_eq!(blast_seq_loc_new(&mut adjacent, 5, 9), 0);
        assert!(s_seq_loc_list_invert(adjacent.as_deref(), 10).is_none());

        assert!(s_seq_loc_list_invert(None, 10).is_none());
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
        let ncbi = crate::encoding::encode_ncbistdaa_sequence(&ascii);
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

    #[test]
    fn test_seg_lifecycle_helpers_match_c_shapes() {
        let mut seq = s_ssequence_new();
        assert_eq!(seq.length, 0);
        assert_eq!(seq.parent_length, 0);
        assert!(!seq.dash);
        assert_eq!(seq.entropy, 0.0);
        seq.seq.extend_from_slice(&[1, 2, 3]);
        seq.composition.extend_from_slice(&[2, 1]);
        seq.state.extend_from_slice(&[2, 1]);
        assert!(s_ssequence_free(Some(seq)).is_none());
        assert!(s_ssequence_free(None).is_none());

        let seg = Some(Box::new(SSeg {
            begin: 0,
            end: 3,
            next: Some(Box::new(SSeg {
                begin: 5,
                end: 9,
                next: None,
            })),
        }));
        assert!(s_seg_free(seg).is_none());
        assert!(s_seg_free(None).is_none());
    }

    #[test]
    fn test_seg_window_helpers_open_shift_entropy_and_close() {
        let mut parent = s_ssequence_new();
        parent.seq = vec![1, 3, 4, 5, 6, 7];
        parent.length = parent.seq.len();
        parent.parent_length = parent.seq.len();

        let mut win = s_open_win(&parent, 1, 3).expect("window");
        assert_eq!(win.seq, vec![3, 4, 5]);
        assert_eq!(win.start, 1);
        assert_eq!(win.parent_length, 6);
        assert_eq!(win.state, vec![1, 1, 1]);
        assert_eq!(win.bogus, 0);

        s_entropy_on(&mut win);
        assert!(win.entropy > 0.0);

        assert!(s_shift_win1(&mut win));
        assert_eq!(win.seq, vec![4, 5, 6]);
        assert_eq!(win.start, 2);
        assert_eq!(win.left, 2);
        assert_eq!(win.right, 5);
        assert!(s_shift_win1(&mut win));
        assert_eq!(win.seq, vec![5, 6, 7]);
        assert!(!s_shift_win1(&mut win));

        let mut sv = vec![3, 2, 1];
        s_decrement_sv(&mut sv, 2);
        assert_eq!(sv, vec![3, 1, 1]);
        s_increment_sv(&mut sv, 1);
        assert_eq!(sv, vec![3, 2, 1]);

        assert!(s_close_win(Some(win)).is_none());
        assert!(s_open_win(&parent, 5, 2).is_none());
    }
}
