//! Rust ports of low-level PHI-BLAST pattern bit helpers.

pub const PHI_BITS_PACKED_PER_WORD: i32 = 30;
pub const PHI_ASCII_SIZE: usize = 256;
pub const PHI_MAX_HIT: usize = 20_000;
pub const PHI_BUF_SIZE: usize = 100;
pub const PHI_MAX_WORD_SIZE: usize = 11;
pub const PHI_MAX_WORDS_IN_PATTERN: usize = 100;
pub const PHI_MAX_PATTERN_LENGTH: usize = PHI_BITS_PACKED_PER_WORD as usize * PHI_MAX_WORD_SIZE;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PatternType {
    #[default]
    OneWord = 0,
    MultiWord = 1,
    VeryLong = 2,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaShortPatternItems {
    pub dna_which_prefix_pos: [u32; PHI_ASCII_SIZE],
    pub dna_which_suffix_pos: [u32; PHI_ASCII_SIZE],
}

impl Default for DnaShortPatternItems {
    fn default() -> Self {
        Self {
            dna_which_prefix_pos: [0; PHI_ASCII_SIZE],
            dna_which_suffix_pos: [0; PHI_ASCII_SIZE],
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShortPatternItems {
    pub match_mask: i32,
    pub which_position: [i32; PHI_ASCII_SIZE],
    pub dna_items: Option<DnaShortPatternItems>,
}

impl Default for ShortPatternItems {
    fn default() -> Self {
        Self {
            match_mask: 0,
            which_position: [0; PHI_ASCII_SIZE],
            dna_items: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtraLongPatternItems {
    pub num_places_in_word: [i32; PHI_MAX_WORDS_IN_PATTERN],
    pub spacing: [i32; PHI_MAX_WORDS_IN_PATTERN],
    pub highest_place: i32,
    pub which_most_specific: i32,
}

impl Default for ExtraLongPatternItems {
    fn default() -> Self {
        Self {
            num_places_in_word: [0; PHI_MAX_WORDS_IN_PATTERN],
            spacing: [0; PHI_MAX_WORDS_IN_PATTERN],
            highest_place: 0,
            which_most_specific: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaLongPatternItems {
    pub dna_prefix_sll: Vec<[u32; PHI_ASCII_SIZE]>,
    pub dna_suffix_sll: Vec<[u32; PHI_ASCII_SIZE]>,
}

impl Default for DnaLongPatternItems {
    fn default() -> Self {
        Self {
            dna_prefix_sll: Vec::new(),
            dna_suffix_sll: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LongPatternItems {
    pub num_words: i32,
    pub match_mask_l: [i32; PHI_BUF_SIZE],
    pub bit_pattern_by_letter: Vec<[i32; PHI_MAX_WORD_SIZE]>,
    pub sll: Vec<[i32; PHI_ASCII_SIZE]>,
    pub input_pattern_masked: Vec<i32>,
    pub dna_items: Option<DnaLongPatternItems>,
    pub extra_long_items: Option<ExtraLongPatternItems>,
}

impl Default for LongPatternItems {
    fn default() -> Self {
        Self {
            num_words: 0,
            match_mask_l: [0; PHI_BUF_SIZE],
            bit_pattern_by_letter: vec![[0; PHI_MAX_WORD_SIZE]; PHI_ASCII_SIZE],
            sll: Vec::new(),
            input_pattern_masked: vec![0; PHI_MAX_PATTERN_LENGTH],
            dna_items: None,
            extra_long_items: None,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct PhiPatternSearchBlk {
    pub flag_pattern_length: PatternType,
    pub pattern_probability: f64,
    pub min_pattern_match_length: i32,
    pub one_word_items: ShortPatternItems,
    pub multi_word_items: Option<LongPatternItems>,
    pub num_patterns_db: i32,
    pub pattern: Option<String>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SphiPatternInfo {
    pub offset: i32,
    pub length: i32,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SphiQueryInfo {
    pub num_patterns: i32,
    pub occurrences: Vec<SphiPatternInfo>,
    pub allocated_size: i32,
    pub probability: f64,
    pub pattern: Option<String>,
}

/// Port of NCBI `SPHIQueryInfoNew` (`pattern.c:478`).
pub fn sphi_query_info_new() -> Option<SphiQueryInfo> {
    const MIN_PHI_LOOKUP_SIZE: usize = 8;
    Some(SphiQueryInfo {
        allocated_size: MIN_PHI_LOOKUP_SIZE as i32,
        occurrences: Vec::with_capacity(MIN_PHI_LOOKUP_SIZE),
        ..Default::default()
    })
}

/// Rust ownership equivalent of NCBI `SPHIQueryInfoFree`.
pub fn sphi_query_info_free(mut query_info: Option<SphiQueryInfo>) -> Option<SphiQueryInfo> {
    if let Some(query_info) = query_info.as_mut() {
        query_info.occurrences.clear();
        query_info.pattern = None;
        query_info.num_patterns = 0;
        query_info.allocated_size = 0;
    }
    None
}

/// Port of NCBI `SPHIQueryInfoCopy` (`pattern.c:507`).
pub fn sphi_query_info_copy(pat_info: Option<&SphiQueryInfo>) -> Option<SphiQueryInfo> {
    let pat_info = pat_info?;
    let mut retval = SphiQueryInfo {
        num_patterns: pat_info.num_patterns,
        occurrences: Vec::with_capacity(pat_info.occurrences.len()),
        allocated_size: pat_info.allocated_size,
        probability: pat_info.probability,
        pattern: None,
    };

    if let Some(pattern) = pat_info.pattern.as_ref() {
        retval.pattern = Some(pattern.clone());
    }

    for occurrence in pat_info.occurrences.iter() {
        retval.occurrences.push(SphiPatternInfo {
            offset: occurrence.offset,
            length: occurrence.length,
        });
    }

    Some(retval)
}

fn phi_alphabet_size(is_dna: bool) -> usize {
    if is_dna {
        crate::encoding::BLASTNA_SIZE
    } else {
        crate::encoding::BLASTAA_SIZE
    }
}

fn phi_all_alphabet_bits(is_dna: bool) -> i32 {
    (1 << phi_alphabet_size(is_dna)) - 1
}

fn phi_order(byte: u8, is_dna: bool) -> usize {
    let index = (byte as usize) & 0x7f;
    if is_dna {
        crate::encoding::IUPACNA_TO_BLASTNA[index] as usize
    } else {
        crate::encoding::AMINOACID_TO_NCBISTDAA[index] as usize
    }
}

fn phi_default_probabilities(is_dna: bool) -> Vec<f64> {
    if is_dna {
        vec![1.0 / crate::encoding::BLASTNA_SIZE as f64; crate::encoding::BLASTNA_SIZE]
    } else {
        crate::stat::protein_std_freq_ncbistdaa().to_vec()
    }
}

fn ncbi2na_unpack_base(byte: i32, base_index: i32) -> usize {
    ((byte >> (2 * base_index)) & 3) as usize
}

/// Port of NCBI `s_FindPrefixAndSuffixPos` (`phi_lookup.c:43`).
pub fn s_find_prefix_and_suffix_pos(
    s: &[i32],
    mask: i32,
    mask2: i32,
    prefix_pos: &mut [u32; PHI_ASCII_SIZE],
    suffix_pos: &mut [u32; PHI_ASCII_SIZE],
) {
    let mask_left_plus_one = (mask << 1) + 1;
    for i in 0..PHI_ASCII_SIZE {
        let a1 = ncbi2na_unpack_base(i as i32, 3);
        let a2 = ncbi2na_unpack_base(i as i32, 2);
        let a3 = ncbi2na_unpack_base(i as i32, 1);
        let a4 = ncbi2na_unpack_base(i as i32, 0);
        let mut tmp =
            ((s.get(a4).copied().unwrap_or(0) >> 1) | mask) & s.get(a3).copied().unwrap_or(0);
        tmp = ((tmp >> 1) | mask) & s.get(a2).copied().unwrap_or(0);
        prefix_pos[i] = (mask2 & ((tmp >> 1) | mask) & s.get(a1).copied().unwrap_or(0)) as u32;

        tmp = ((s.get(a1).copied().unwrap_or(0) << 1) | mask_left_plus_one)
            & s.get(a2).copied().unwrap_or(0);
        tmp = ((tmp << 1) | mask_left_plus_one) & s.get(a3).copied().unwrap_or(0);
        suffix_pos[i] = (((((tmp << 1) | mask_left_plus_one) & s.get(a4).copied().unwrap_or(0))
            << 1)
            | mask_left_plus_one) as u32;
    }
}

/// Port of NCBI `s_InitDNAPattern` (`phi_lookup.c:82`).
pub fn s_init_dna_pattern(pattern_blk: &mut PhiPatternSearchBlk) {
    if pattern_blk.flag_pattern_length != PatternType::OneWord {
        let Some(multiword_items) = pattern_blk.multi_word_items.as_mut() else {
            return;
        };
        let num_words = multiword_items.num_words.max(0) as usize;
        let mut dna_items = multiword_items.dna_items.take().unwrap_or_default();
        dna_items
            .dna_prefix_sll
            .resize(num_words, [0; PHI_ASCII_SIZE]);
        dna_items
            .dna_suffix_sll
            .resize(num_words, [0; PHI_ASCII_SIZE]);
        for word_index in 0..num_words {
            let mask1 = multiword_items.match_mask_l[word_index];
            let composite_mask = mask1 + (mask1 >> 1) + (mask1 >> 2) + (mask1 >> 3);
            let sll = multiword_items
                .sll
                .get(word_index)
                .copied()
                .unwrap_or([0; PHI_ASCII_SIZE]);
            s_find_prefix_and_suffix_pos(
                &sll,
                mask1,
                composite_mask,
                &mut dna_items.dna_prefix_sll[word_index],
                &mut dna_items.dna_suffix_sll[word_index],
            );
        }
        multiword_items.dna_items = Some(dna_items);
    } else {
        let word_items = &mut pattern_blk.one_word_items;
        let mut dna_items = word_items.dna_items.take().unwrap_or_default();
        let match_mask = word_items.match_mask;
        let composite_mask = match_mask + (match_mask >> 1) + (match_mask >> 2) + (match_mask >> 3);
        s_find_prefix_and_suffix_pos(
            &word_items.which_position,
            match_mask,
            composite_mask,
            &mut dna_items.dna_which_prefix_pos,
            &mut dna_items.dna_which_suffix_pos,
        );
        word_items.dna_items = Some(dna_items);
    }
}

/// Port of NCBI `s_ExpandPattern` (`phi_lookup.c:130`).
pub fn s_expand_pattern(
    input_pattern_masked: &mut [i32],
    input_pattern: &mut [u8],
    length: i32,
    max_length: i32,
    all_alphabet_bits: i32,
) -> i32 {
    for i in 0..length.max(0) as usize {
        let this_place_masked = -input_pattern_masked[i];
        if this_place_masked > 0 {
            input_pattern_masked[i] = all_alphabet_bits;
            let mut temp_pattern_mask = vec![0; length.max(0) as usize];
            let mut temp_pattern = vec![0; length.max(0) as usize];
            for j in 0..length.max(0) as usize {
                temp_pattern_mask[j] = input_pattern_masked[j];
                temp_pattern[j] = input_pattern[j];
            }
            let mut rec_return_value_2 = s_expand_pattern(
                input_pattern_masked,
                input_pattern,
                length,
                max_length,
                all_alphabet_bits,
            );
            if rec_return_value_2 == -1 {
                return -1;
            }
            let mut rec_return_value_1 = rec_return_value_2;
            for num_pos in 0..=this_place_masked {
                if num_pos == 1 {
                    continue;
                }
                for k in 0..length.max(0) as usize {
                    if k == i {
                        for _ in 0..num_pos {
                            if rec_return_value_1 >= max_length {
                                return -1;
                            }
                            input_pattern_masked[rec_return_value_1 as usize] = all_alphabet_bits;
                            rec_return_value_1 += 1;
                        }
                    } else {
                        if rec_return_value_1 >= max_length {
                            return -1;
                        }
                        input_pattern_masked[rec_return_value_1 as usize] = temp_pattern_mask[k];
                        input_pattern[rec_return_value_1 as usize] = temp_pattern[k];
                        rec_return_value_1 += 1;
                    }
                    if rec_return_value_1 >= max_length {
                        return -1;
                    }
                }
                let start = rec_return_value_2.max(0) as usize;
                rec_return_value_1 = s_expand_pattern(
                    &mut input_pattern_masked[start..],
                    &mut input_pattern[start..],
                    length + num_pos - 1,
                    max_length - rec_return_value_2,
                    all_alphabet_bits,
                );
                if rec_return_value_1 == -1 {
                    return -1;
                }
                rec_return_value_2 += rec_return_value_1;
                rec_return_value_1 = rec_return_value_2;
            }
            return rec_return_value_1;
        }
    }
    length
}

/// Port of NCBI `s_PackPattern` (`phi_lookup.c:200`).
pub fn s_pack_pattern(input_pattern: &[u8], length: i32) -> i32 {
    let mut return_value = 0;
    for i in 0..length.max(0) as usize {
        if input_pattern.get(i).copied().unwrap_or(0) != 0 {
            return_value += 1 << i;
        }
    }
    return_value
}

/// Port of NCBI `s_PackLongPattern` (`phi_lookup.c:218`).
pub fn s_pack_long_pattern(
    num_places: i32,
    input_pattern: &[u8],
    pattern_blk: &mut PhiPatternSearchBlk,
) {
    let Some(multiword_items) = pattern_blk.multi_word_items.as_mut() else {
        return;
    };
    multiword_items.num_words = (num_places - 1) / PHI_BITS_PACKED_PER_WORD + 1;
    for word_index in 0..multiword_items.num_words.max(0) as usize {
        let mut bit_pattern = 0;
        for i in 0..PHI_BITS_PACKED_PER_WORD as usize {
            if input_pattern
                .get(word_index * PHI_BITS_PACKED_PER_WORD as usize + i)
                .copied()
                .unwrap_or(0)
                != 0
            {
                bit_pattern += 1 << i;
            }
        }
        multiword_items.match_mask_l[word_index] = bit_pattern;
    }
    for char_index in 0..crate::encoding::BLASTAA_SIZE {
        for word_index in 0..multiword_items.num_words.max(0) as usize {
            let mut bit_pattern = 0;
            for i in 0..PHI_BITS_PACKED_PER_WORD as usize {
                let pattern_index = word_index * PHI_BITS_PACKED_PER_WORD as usize + i;
                if ((1 << char_index)
                    & multiword_items
                        .input_pattern_masked
                        .get(pattern_index)
                        .copied()
                        .unwrap_or(0))
                    != 0
                {
                    bit_pattern |= 1 << i;
                }
            }
            multiword_items.bit_pattern_by_letter[char_index][word_index] = bit_pattern;
        }
    }
}

/// Port of NCBI `s_NumOfOne` (`phi_lookup.c:259`).
pub fn s_num_of_one(mut a: i32) -> i32 {
    let mut return_value = 0;
    while a > 0 {
        if a % 2 == 1 {
            return_value += 1;
        }
        a >>= 1;
    }
    return_value
}

/// Port of NCBI `s_PackVeryLongPattern` (`phi_lookup.c:275`).
pub fn s_pack_very_long_pattern(
    input_pattern_masked: &[i32],
    num_places_in_pattern: i32,
    pattern_blk: &mut PhiPatternSearchBlk,
) {
    let Some(multiword_items) = pattern_blk.multi_word_items.as_mut() else {
        return;
    };
    let mut extra_items = ExtraLongPatternItems::default();
    let mut most_specific = 1.0;
    let mut pattern_word_probability = 1.0;
    let mut word_index = 0usize;
    let mut place_in_word = 0usize;
    let mut place_index = 0usize;

    while place_index <= num_places_in_pattern.max(0) as usize {
        if place_index == num_places_in_pattern.max(0) as usize
            || input_pattern_masked.get(place_index).copied().unwrap_or(0) < 0
            || place_in_word == PHI_BITS_PACKED_PER_WORD as usize
        {
            if word_index >= PHI_MAX_WORDS_IN_PATTERN {
                break;
            }
            multiword_items.match_mask_l[word_index] = 1 << (place_in_word as i32 - 1);
            if multiword_items.sll.len() <= word_index {
                multiword_items
                    .sll
                    .resize(word_index + 1, [0; PHI_ASCII_SIZE]);
            }
            for char_index in 0..crate::encoding::BLASTAA_SIZE {
                let mut one_word_mask = 0;
                for place_in_word_2 in 0..place_in_word {
                    let source_index = place_index - place_in_word + place_in_word_2;
                    if ((1 << char_index)
                        & input_pattern_masked.get(source_index).copied().unwrap_or(0))
                        != 0
                    {
                        one_word_mask |= 1 << place_in_word_2;
                    }
                }
                multiword_items.sll[word_index][char_index] = one_word_mask;
            }
            extra_items.num_places_in_word[word_index] = place_in_word as i32;
            if pattern_word_probability < most_specific {
                most_specific = pattern_word_probability;
                extra_items.which_most_specific = word_index as i32;
            }
            if place_index == num_places_in_pattern.max(0) as usize {
                extra_items.spacing[word_index] = 0;
                word_index += 1;
            } else if input_pattern_masked.get(place_index).copied().unwrap_or(0) < 0 {
                extra_items.spacing[word_index] = -input_pattern_masked[place_index];
                word_index += 1;
            } else {
                place_index = place_index.saturating_sub(1);
                extra_items.spacing[word_index] = 0;
                word_index += 1;
            }
            place_in_word = 0;
            pattern_word_probability = 1.0;
        } else {
            pattern_word_probability *= s_num_of_one(input_pattern_masked[place_index]) as f64
                / crate::encoding::BLASTAA_SIZE as f64;
            place_in_word += 1;
        }
        place_index += 1;
    }
    multiword_items.num_words = word_index as i32;
    multiword_items.extra_long_items = Some(extra_items);
}

/// Port of NCBI `s_PatternSearchItemsInit` (`phi_lookup.c:354`).
pub fn s_pattern_search_items_init() -> PhiPatternSearchBlk {
    PhiPatternSearchBlk {
        flag_pattern_length: PatternType::OneWord,
        pattern_probability: 1.0,
        min_pattern_match_length: 0,
        multi_word_items: Some(LongPatternItems::default()),
        ..Default::default()
    }
}

/// Port of NCBI `s_MakePatternUpperCase` (`phi_lookup.c:374`).
pub fn s_make_pattern_upper_case(pattern_in: &str) -> String {
    pattern_in
        .bytes()
        .map(|byte| {
            if byte.is_ascii_lowercase() {
                byte.to_ascii_uppercase()
            } else {
                byte
            }
        })
        .map(char::from)
        .collect()
}

fn parse_number_until(bytes: &[u8], index: &mut usize, stop: &[u8]) -> i32 {
    let start = *index;
    while *index < bytes.len() && !stop.contains(&bytes[*index]) {
        *index += 1;
    }
    std::str::from_utf8(&bytes[start..*index])
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(0)
}

/// Port of NCBI `SPHIPatternSearchBlkNew` (`phi_lookup.c:388`).
///
/// C obtains residue probabilities from `BlastScoreBlk`; Rust accepts them as
/// an optional NCBIstdaa/BLASTNA-indexed slice and otherwise uses the same
/// standard background defaults used by the scoring module.
pub fn sphi_pattern_search_blk_new(
    pattern_in: &str,
    is_dna: bool,
    residue_probabilities: Option<&[f64]>,
) -> Result<PhiPatternSearchBlk, String> {
    const WILDCARD_THRESHOLD: i32 = 30;
    let all_alphabet_bits = phi_all_alphabet_bits(is_dna);
    let probabilities = residue_probabilities
        .map(|p| p.to_vec())
        .unwrap_or_else(|| phi_default_probabilities(is_dna));
    let mut pattern_blk = s_pattern_search_items_init();
    pattern_blk.pattern = Some(s_make_pattern_upper_case(pattern_in));
    let pattern = pattern_blk
        .pattern
        .as_ref()
        .expect("pattern is set")
        .clone();
    let bytes = pattern.as_bytes();
    if bytes.len() >= PHI_MAX_PATTERN_LENGTH {
        return Err(format!(
            "Pattern is too long ({} but only {} supported)",
            bytes.len(),
            PHI_MAX_PATTERN_LENGTH
        ));
    }

    let mut local_pattern = vec![0u8; PHI_MAX_PATTERN_LENGTH];
    let mut pos_index = 0usize;
    let mut current_wildcard_product = 1;
    let mut wildcard_product = 1;
    let mut current_set_mask = 0i32;

    let Some(multiword_items) = pattern_blk.multi_word_items.as_mut() else {
        return Err("pattern search block missing multiword items".to_string());
    };

    let mut char_index = 0usize;
    while char_index < bytes.len() {
        let mut next_char = bytes[char_index];
        if matches!(next_char, b'\0' | b'\r' | b'\n') {
            break;
        }
        if matches!(next_char, b'-' | b'.' | b'>' | b' ' | b'<') {
            char_index += 1;
            continue;
        }

        let char_set_mask: i32;
        let position_probability: f64;
        if next_char != b'[' && next_char != b'{' {
            if next_char == b'X' {
                if bytes.get(char_index + 1).copied() == Some(b'(') {
                    char_index += 1;
                    let mut second_index = char_index;
                    while second_index < bytes.len()
                        && bytes[second_index] != b','
                        && bytes[second_index] != b')'
                    {
                        second_index += 1;
                    }
                    if bytes.get(second_index).copied() == Some(b')') {
                        char_index -= 1;
                        char_set_mask = all_alphabet_bits;
                        position_probability = 1.0;
                    } else {
                        char_index += 1;
                        let mut parse_index = char_index;
                        let min_wildcard = parse_number_until(bytes, &mut parse_index, b",");
                        parse_index += usize::from(parse_index < bytes.len());
                        let max_wildcard_raw = parse_number_until(bytes, &mut parse_index, b")");
                        let max_wildcard = max_wildcard_raw - min_wildcard;
                        current_wildcard_product *= max_wildcard + 1;
                        if current_wildcard_product > wildcard_product {
                            wildcard_product = current_wildcard_product;
                        }
                        pattern_blk.min_pattern_match_length += min_wildcard;
                        for _ in 0..min_wildcard.max(0) {
                            if pos_index >= PHI_MAX_PATTERN_LENGTH {
                                return Err("Pattern too long".to_string());
                            }
                            multiword_items.input_pattern_masked[pos_index] = all_alphabet_bits;
                            pos_index += 1;
                        }
                        if max_wildcard != 0 {
                            if pos_index >= PHI_MAX_PATTERN_LENGTH {
                                return Err("Pattern too long".to_string());
                            }
                            multiword_items.input_pattern_masked[pos_index] = -max_wildcard;
                            pos_index += 1;
                            pattern_blk.pattern_probability *= max_wildcard as f64;
                        }
                        while char_index < bytes.len() && bytes[char_index] != b')' {
                            char_index += 1;
                        }
                        char_index += 1;
                        continue;
                    }
                } else {
                    char_set_mask = all_alphabet_bits;
                    position_probability = 1.0;
                }
            } else if next_char == b'U' && !is_dna {
                char_set_mask = all_alphabet_bits * 2 + 1;
                position_probability = 1.0;
            } else {
                let prev_set_mask = current_set_mask;
                let order = phi_order(next_char, is_dna);
                current_set_mask = 1 << order;
                char_set_mask = current_set_mask;
                if (prev_set_mask & current_set_mask) == 0 {
                    current_wildcard_product = 1;
                }
                position_probability = probabilities.get(order).copied().unwrap_or(0.0);
            }
        } else if next_char == b'[' {
            let mut mask = 0;
            let mut probability = 0.0;
            loop {
                char_index += 1;
                next_char = *bytes
                    .get(char_index)
                    .ok_or_else(|| "pattern description has an unterminated bracket".to_string())?;
                if next_char == b']' {
                    break;
                }
                if !next_char.is_ascii_alphabetic() {
                    return Err(
                        "pattern description has a non-alphabetic character inside a bracket"
                            .to_string(),
                    );
                }
                let order = phi_order(next_char, is_dna);
                mask |= 1 << order;
                probability += probabilities.get(order).copied().unwrap_or(0.0);
            }
            let prev_set_mask = current_set_mask;
            current_set_mask = mask;
            if (prev_set_mask & current_set_mask) == 0 {
                current_wildcard_product = 1;
            }
            char_set_mask = mask;
            position_probability = probability;
        } else {
            let mut mask = all_alphabet_bits;
            let mut probability = 1.0;
            loop {
                char_index += 1;
                next_char = *bytes
                    .get(char_index)
                    .ok_or_else(|| "pattern description has an unterminated brace".to_string())?;
                if next_char == b'}' {
                    break;
                }
                let order = phi_order(next_char, is_dna);
                mask -= mask & (1 << order);
                probability -= probabilities.get(order).copied().unwrap_or(0.0);
            }
            let prev_set_mask = current_set_mask;
            current_set_mask = mask;
            if (prev_set_mask & current_set_mask) == 0 {
                current_wildcard_product = 1;
            }
            char_set_mask = mask;
            position_probability = probability;
        }

        if bytes.get(char_index + 1).copied() == Some(b'(') {
            char_index += 2;
            let num_identical = parse_number_until(bytes, &mut char_index, b")");
            pattern_blk.min_pattern_match_length += num_identical;
            while char_index < bytes.len() && bytes[char_index] != b')' {
                char_index += 1;
            }
            for _ in 0..num_identical.max(0) {
                if pos_index >= PHI_MAX_PATTERN_LENGTH {
                    return Err("Pattern is too long".to_string());
                }
                multiword_items.input_pattern_masked[pos_index] = char_set_mask;
                pos_index += 1;
                pattern_blk.pattern_probability *= position_probability;
            }
        } else {
            if pos_index >= PHI_MAX_PATTERN_LENGTH {
                return Err("Pattern is too long".to_string());
            }
            multiword_items.input_pattern_masked[pos_index] = char_set_mask;
            pos_index += 1;
            pattern_blk.min_pattern_match_length += 1;
            pattern_blk.pattern_probability *= position_probability;
        }
        if pos_index >= PHI_MAX_PATTERN_LENGTH {
            return Err("Pattern is too long".to_string());
        }
        char_index += 1;
    }

    while pos_index > 0 && multiword_items.input_pattern_masked[pos_index - 1] < 0 {
        pos_index -= 1;
    }

    let mut first_non_wild = 0usize;
    while first_non_wild < pos_index
        && multiword_items.input_pattern_masked[first_non_wild] == all_alphabet_bits
    {
        first_non_wild += 1;
    }
    if first_non_wild < pos_index && multiword_items.input_pattern_masked[first_non_wild] < 0 {
        let mut second_index = first_non_wild + 1;
        while second_index < pos_index && multiword_items.input_pattern_masked[second_index] <= 0 {
            second_index += 1;
        }
        let mut write_index = first_non_wild;
        while second_index < pos_index {
            multiword_items.input_pattern_masked[write_index] =
                multiword_items.input_pattern_masked[second_index];
            write_index += 1;
            second_index += 1;
        }
        pos_index = write_index;
    }

    if pos_index == 0 {
        return Err("Pattern is empty".to_string());
    }
    local_pattern[pos_index - 1] = 1;
    if pattern_blk.pattern_probability > 1.0 {
        pattern_blk.pattern_probability = 1.0;
    }

    let temp_input_pattern_masked = multiword_items.input_pattern_masked[..pos_index].to_vec();
    let temp_pos_index = pos_index;
    let expanded_pos_index = s_expand_pattern(
        &mut multiword_items.input_pattern_masked,
        &mut local_pattern,
        pos_index as i32,
        PHI_MAX_PATTERN_LENGTH as i32,
        all_alphabet_bits,
    );
    if expanded_pos_index == -1 || (expanded_pos_index > PHI_BITS_PACKED_PER_WORD && is_dna) {
        pattern_blk.flag_pattern_length = PatternType::VeryLong;
        s_pack_very_long_pattern(
            &temp_input_pattern_masked,
            temp_pos_index as i32,
            &mut pattern_blk,
        );
        let multiword_items = pattern_blk
            .multi_word_items
            .as_mut()
            .expect("multiword items retained");
        for (index, value) in temp_input_pattern_masked.iter().copied().enumerate() {
            multiword_items.input_pattern_masked[index] = value;
        }
        if let Some(extra) = multiword_items.extra_long_items.as_mut() {
            extra.highest_place = temp_pos_index as i32;
        }
        if is_dna {
            s_init_dna_pattern(&mut pattern_blk);
        }
        return Ok(pattern_blk);
    }
    if expanded_pos_index > PHI_BITS_PACKED_PER_WORD {
        pattern_blk.flag_pattern_length = PatternType::MultiWord;
        s_pack_long_pattern(expanded_pos_index, &local_pattern, &mut pattern_blk);
        return Ok(pattern_blk);
    }

    pattern_blk.one_word_items.match_mask = s_pack_pattern(&local_pattern, expanded_pos_index);
    let multiword_items = pattern_blk
        .multi_word_items
        .as_ref()
        .expect("multiword items retained");
    for char_index in 0..phi_alphabet_size(is_dna) {
        let mut this_mask = 0;
        for char_set_mask in 0..expanded_pos_index.max(0) as usize {
            if ((1 << char_index) & multiword_items.input_pattern_masked[char_set_mask]) != 0 {
                this_mask |= 1 << char_set_mask;
            }
        }
        pattern_blk.one_word_items.which_position[char_index] = this_mask;
    }
    if is_dna {
        s_init_dna_pattern(&mut pattern_blk);
    }

    let _warn_too_many_wildcards = wildcard_product > WILDCARD_THRESHOLD;
    Ok(pattern_blk)
}

/// Rust ownership equivalent of NCBI `SPHIPatternSearchBlkFree`.
pub fn sphi_pattern_search_blk_free(
    pattern_blk: Option<PhiPatternSearchBlk>,
) -> Option<PhiPatternSearchBlk> {
    let Some(mut pattern_blk) = pattern_blk else {
        return None;
    };
    if let Some(mut multiword_items) = pattern_blk.multi_word_items.take() {
        multiword_items.extra_long_items.take();
        multiword_items.dna_items.take();
        multiword_items.input_pattern_masked.clear();
        multiword_items.sll.clear();
        multiword_items.bit_pattern_by_letter.clear();
    }
    if pattern_blk.flag_pattern_length != PatternType::VeryLong {
        pattern_blk.one_word_items.dna_items.take();
        pattern_blk.one_word_items.which_position = [0; PHI_ASCII_SIZE];
    }
    pattern_blk.pattern.take();
    None
}

/// Port of NCBI `s_PHIBlastAddPatternHit` (`pattern.c:530`).
pub fn s_phi_blast_add_pattern_hit(
    pattern_info: Option<&mut SphiQueryInfo>,
    offset: i32,
    length: i32,
) -> i16 {
    let Some(pattern_info) = pattern_info else {
        return -1;
    };
    if pattern_info.num_patterns >= pattern_info.allocated_size {
        pattern_info.allocated_size *= 2;
        let additional = (pattern_info.allocated_size as usize)
            .saturating_sub(pattern_info.occurrences.capacity());
        pattern_info.occurrences.reserve(additional);
    }
    pattern_info
        .occurrences
        .push(SphiPatternInfo { offset, length });
    pattern_info.num_patterns += 1;
    0
}

/// Port of NCBI `_PHIGetRightOneBits` (`pattern.c:61`).
pub fn _phi_get_right_one_bits(s: i32, mask: i32) -> (i32, i32) {
    let checking_matches = s & mask;
    let mut left_index = -1;
    let mut right_index = 0;

    while right_index < PHI_BITS_PACKED_PER_WORD {
        if ((checking_matches >> right_index) & 1) == 1 {
            break;
        }
        if ((mask >> right_index) & 1) == 1 {
            left_index = right_index;
        }
        right_index += 1;
    }

    if right_index == PHI_BITS_PACKED_PER_WORD {
        right_index = 0;
    }
    (right_index, left_index)
}

/// Port of NCBI `s_LenOf` (`pattern.c:93`).
pub fn s_len_of(s: i32, mask: i32) -> i32 {
    let (right_one, right_mask_only) = _phi_get_right_one_bits(s, mask);
    right_one - right_mask_only
}

/// Port of NCBI `_PHIPatternWordsLeftShift` (`pattern.c:236`).
pub fn _phi_pattern_words_left_shift(a: &mut [i32], mut b: u8, num_words: i32) {
    let overflow_threshold = 1 << PHI_BITS_PACKED_PER_WORD;
    for value in a.iter_mut().take(num_words.max(0) as usize) {
        let x = (*value << 1) + i32::from(b);
        if x >= overflow_threshold {
            *value = x - overflow_threshold;
            b = 1;
        } else {
            *value = x;
            b = 0;
        }
    }
}

/// Port of NCBI `_PHIPatternWordsBitwiseOr` (`pattern.c:257`).
pub fn _phi_pattern_words_bitwise_or(a: &mut [i32], b: &[i32], num_words: i32) {
    for (left, right) in a.iter_mut().zip(b.iter()).take(num_words.max(0) as usize) {
        *left |= *right;
    }
}

/// Port of NCBI `_PHIPatternWordsBitwiseAnd` (`pattern.c:265`).
pub fn _phi_pattern_words_bitwise_and(
    result: &mut [i32],
    a: &[i32],
    b: &[i32],
    num_words: i32,
) -> i32 {
    let mut return_value = 0;
    for ((out, left), right) in result
        .iter_mut()
        .zip(a.iter())
        .zip(b.iter())
        .take(num_words.max(0) as usize)
    {
        *out = *left & *right;
        if *out != 0 {
            return_value = 1;
        }
    }
    return_value
}

/// Port of NCBI `s_LenOfL` (`pattern.c:286`).
pub fn s_len_of_l(s: &[i32], mask: &[i32], num_words: i32) -> i32 {
    let mut first_one_in_mask = -1;
    for word_index in 0..num_words.max(0) as usize {
        let s_word = s.get(word_index).copied().unwrap_or(0);
        let mask_word = mask.get(word_index).copied().unwrap_or(0);
        for bit_index in 0..PHI_BITS_PACKED_PER_WORD {
            let absolute = word_index as i32 * PHI_BITS_PACKED_PER_WORD + bit_index;
            if ((s_word >> bit_index) & 1) == 1 {
                return absolute - first_one_in_mask;
            }
            if ((mask_word >> bit_index) & 1) == 1 {
                first_one_in_mask = absolute;
            }
        }
    }
    -1
}

/// Port of NCBI `_PHIBlastFindHitsShort` (`pattern.c:104`).
pub fn _phi_blast_find_hits_short(
    hit_array: &mut Vec<i32>,
    seq: &[u8],
    len1: i32,
    pattern_blk: &PhiPatternSearchBlk,
) -> i32 {
    let pattern_items = &pattern_blk.one_word_items;
    let mask = pattern_items.match_mask;
    let mask_shift_plus_1 = (mask << 1) + 1;
    let mut prefix_matched_bit_pattern = 0i32;
    let mut num_matches = 0usize;
    let len = len1.max(0) as usize;

    for (i, &letter) in seq.iter().take(len).enumerate() {
        prefix_matched_bit_pattern = ((prefix_matched_bit_pattern << 1) | mask_shift_plus_1)
            & pattern_items.which_position[letter as usize];
        if (prefix_matched_bit_pattern & mask) != 0 {
            hit_array.push(i as i32);
            hit_array.push(i as i32 - s_len_of(prefix_matched_bit_pattern, mask) + 1);
            num_matches += 2;
            if num_matches == PHI_MAX_HIT {
                break;
            }
        }
    }

    num_matches as i32
}

/// Port of NCBI `s_FindHitsShortDNA` (`pattern.c:154`).
pub fn s_find_hits_short_dna(
    hit_array: &mut Vec<i32>,
    seq: &[u8],
    mut pos: i32,
    len: i32,
    pattern_blk: &PhiPatternSearchBlk,
) -> i32 {
    let pattern_items = &pattern_blk.one_word_items;
    let Some(dna_items) = pattern_items.dna_items.as_ref() else {
        return 0;
    };
    let match_mask = pattern_items.match_mask as u32;
    let mask2 = match_mask * PHI_BITS_PACKED_PER_WORD as u32 + 15;
    let mask_shift_plus_1 = ((pattern_items.match_mask << 1) + 1) as u32;
    let mut twice_num_hits = 0usize;

    let mut prefix_matched_bit_pattern;
    let end;
    let remain;
    let mut seq_offset = 0usize;

    if pos != 0 {
        pos = 4 - pos;
        let Some(&first) = seq.first() else {
            return 0;
        };
        prefix_matched_bit_pattern = ((match_mask * (((1u32 << (pos + 1)) - 1) * 2))
            + ((1u32 << (pos + 1)) - 1))
            & dna_items.dna_which_suffix_pos[first as usize];
        seq_offset = 1;
        end = (len - pos) / 4;
        remain = (len - pos) % 4;
    } else {
        prefix_matched_bit_pattern = mask_shift_plus_1;
        end = len / 4;
        remain = len % 4;
    }

    for i in 0..end.max(0) as usize {
        let Some(&byte) = seq.get(seq_offset + i) else {
            return twice_num_hits as i32;
        };
        let mut tmp = prefix_matched_bit_pattern & dna_items.dna_which_prefix_pos[byte as usize];
        if tmp != 0 {
            for j in 0..4 {
                if (tmp & match_mask) != 0 {
                    let end_pos = i as i32 * 4 + j + pos;
                    hit_array.push(end_pos);
                    hit_array
                        .push(end_pos - s_len_of((tmp & match_mask) as i32, match_mask as i32) + 1);
                    twice_num_hits += 2;
                    if twice_num_hits == PHI_MAX_HIT {
                        return twice_num_hits as i32;
                    }
                }
                tmp <<= 1;
            }
        }
        prefix_matched_bit_pattern = ((prefix_matched_bit_pattern << 4) | mask2)
            & dna_items.dna_which_suffix_pos[byte as usize];
    }

    if let Some(&byte) = seq.get(seq_offset + end.max(0) as usize) {
        let mut tmp = prefix_matched_bit_pattern & dna_items.dna_which_prefix_pos[byte as usize];
        for j in 0..remain.max(0) {
            if (tmp & match_mask) != 0 {
                let end_pos = end * 4 + j + pos;
                hit_array.push(end_pos);
                hit_array
                    .push(end_pos - s_len_of((tmp & match_mask) as i32, match_mask as i32) + 1);
                twice_num_hits += 2;
                if twice_num_hits == PHI_MAX_HIT {
                    return twice_num_hits as i32;
                }
            }
            tmp <<= 1;
        }
    }

    twice_num_hits as i32
}

/// Port of NCBI `s_FindHitsShortHead` (`pattern.c:227`).
pub fn s_find_hits_short_head(
    hit_array: &mut Vec<i32>,
    seq: &[u8],
    start: i32,
    len: i32,
    is_dna: u8,
    pattern_blk: &PhiPatternSearchBlk,
) -> i32 {
    if is_dna != 0 {
        let byte_start = (start.max(0) / 4) as usize;
        return s_find_hits_short_dna(
            hit_array,
            seq.get(byte_start..).unwrap_or(&[]),
            start.rem_euclid(4),
            len,
            pattern_blk,
        );
    }
    let start = start.max(0) as usize;
    _phi_blast_find_hits_short(hit_array, seq.get(start..).unwrap_or(&[]), len, pattern_blk)
}

/// Port of NCBI `s_FindHitsLong` (`pattern.c:315`).
pub fn s_find_hits_long(
    hit_array: &mut Vec<i32>,
    seq: &[u8],
    len1: i32,
    pattern_blk: &PhiPatternSearchBlk,
) -> i32 {
    let Some(pattern_items) = pattern_blk.multi_word_items.as_ref() else {
        return 0;
    };
    let num_words = pattern_items.num_words.max(0);
    let words = num_words as usize;
    if words == 0 {
        return 0;
    }

    let mut match_result = vec![0; words];
    let mut mask = pattern_items.match_mask_l[..words.min(PHI_BUF_SIZE)].to_vec();
    mask.resize(words, 0);
    let mut prefix_matched_bit_pattern = vec![0; words];

    _phi_pattern_words_left_shift(&mut mask, 1, num_words);
    let mut twice_num_hits = 0usize;
    for (i, &letter) in seq.iter().take(len1.max(0) as usize).enumerate() {
        _phi_pattern_words_left_shift(&mut prefix_matched_bit_pattern, 0, num_words);
        _phi_pattern_words_bitwise_or(&mut prefix_matched_bit_pattern, &mask, num_words);
        let letter_bits = pattern_items
            .bit_pattern_by_letter
            .get(letter as usize)
            .map(|row| row.as_slice())
            .unwrap_or(&[]);
        for (left, right) in prefix_matched_bit_pattern
            .iter_mut()
            .zip(letter_bits.iter())
            .take(words)
        {
            *left &= *right;
        }
        if _phi_pattern_words_bitwise_and(
            &mut match_result,
            &prefix_matched_bit_pattern,
            &pattern_items.match_mask_l,
            num_words,
        ) != 0
        {
            hit_array.push(i as i32);
            hit_array.push(
                i as i32 - s_len_of_l(&match_result, &pattern_items.match_mask_l, num_words) + 1,
            );
            twice_num_hits += 2;
        }
    }

    twice_num_hits as i32
}

fn short_items_for_long_word(
    pattern_items: &LongPatternItems,
    word_index: usize,
    is_dna: u8,
) -> ShortPatternItems {
    let mut dna_items = None;
    if is_dna != 0 {
        if let Some(long_dna) = pattern_items.dna_items.as_ref() {
            let prefix = long_dna
                .dna_prefix_sll
                .get(word_index)
                .copied()
                .unwrap_or([0; PHI_ASCII_SIZE]);
            let suffix = long_dna
                .dna_suffix_sll
                .get(word_index)
                .copied()
                .unwrap_or([0; PHI_ASCII_SIZE]);
            dna_items = Some(DnaShortPatternItems {
                dna_which_prefix_pos: prefix,
                dna_which_suffix_pos: suffix,
            });
        }
    }

    ShortPatternItems {
        match_mask: pattern_items
            .match_mask_l
            .get(word_index)
            .copied()
            .unwrap_or(0),
        which_position: pattern_items
            .sll
            .get(word_index)
            .copied()
            .unwrap_or([0; PHI_ASCII_SIZE]),
        dna_items,
    }
}

/// Port of NCBI `s_FindHitsVeryLong` (`pattern.c:368`).
pub fn s_find_hits_very_long(
    hit_array: &mut Vec<i32>,
    seq: &[u8],
    len: i32,
    is_dna: u8,
    pattern_blk: &PhiPatternSearchBlk,
) -> i32 {
    let Some(multiword_items) = pattern_blk.multi_word_items.as_ref() else {
        return 0;
    };
    let Some(extra_items) = multiword_items.extra_long_items.as_ref() else {
        return 0;
    };
    let most_specific_word = extra_items.which_most_specific.max(0) as usize;
    if most_specific_word >= multiword_items.num_words.max(0) as usize {
        return 0;
    }

    let word_blk = PhiPatternSearchBlk {
        one_word_items: short_items_for_long_word(multiword_items, most_specific_word, is_dna),
        ..PhiPatternSearchBlk::default()
    };
    let mut twice_num_hits = s_find_hits_short_head(hit_array, seq, 0, len, is_dna, &word_blk);
    if twice_num_hits < 2 {
        hit_array.clear();
        return 0;
    }

    for word_index in (most_specific_word + 1)..multiword_items.num_words.max(0) as usize {
        let word_blk = PhiPatternSearchBlk {
            one_word_items: short_items_for_long_word(multiword_items, word_index, is_dna),
            ..PhiPatternSearchBlk::default()
        };
        let mut hit_array1 = Vec::new();
        for hit in hit_array
            .chunks_exact(2)
            .take((twice_num_hits / 2) as usize)
        {
            let start = hit[0] + 1;
            let scan_len = (len - hit[0] - 1).min(
                extra_items
                    .spacing
                    .get(word_index - 1)
                    .copied()
                    .unwrap_or(0)
                    + extra_items
                        .num_places_in_word
                        .get(word_index)
                        .copied()
                        .unwrap_or(0),
            );
            let base = hit_array1.len();
            let twice_hits_one_call =
                s_find_hits_short_head(&mut hit_array1, seq, start, scan_len, is_dna, &word_blk);
            for idx in (base..base + twice_hits_one_call.max(0) as usize).step_by(2) {
                hit_array1[idx] += hit[0] + 1;
                hit_array1[idx + 1] = hit[1];
            }
        }
        twice_num_hits = hit_array1.len() as i32;
        if twice_num_hits < 2 {
            hit_array.clear();
            return 0;
        }
        hit_array.clear();
        hit_array.extend_from_slice(&hit_array1);
    }

    let mut word_index = most_specific_word as i32 - 1;
    while word_index >= 0 {
        let word_index_usize = word_index as usize;
        let word_blk = PhiPatternSearchBlk {
            one_word_items: short_items_for_long_word(multiword_items, word_index_usize, is_dna),
            ..PhiPatternSearchBlk::default()
        };
        let mut hit_array1 = Vec::new();
        for hit in hit_array
            .chunks_exact(2)
            .take((twice_num_hits / 2) as usize)
        {
            let mut start = hit[1]
                - extra_items
                    .spacing
                    .get(word_index_usize)
                    .copied()
                    .unwrap_or(0)
                - extra_items
                    .num_places_in_word
                    .get(word_index_usize)
                    .copied()
                    .unwrap_or(0);
            if start < 0 {
                start = 0;
            }
            let base = hit_array1.len();
            let twice_hits_one_call = s_find_hits_short_head(
                &mut hit_array1,
                seq,
                start,
                hit[1] - start,
                is_dna,
                &word_blk,
            );
            for idx in (base..base + twice_hits_one_call.max(0) as usize).step_by(2) {
                hit_array1[idx] = hit[0];
                hit_array1[idx + 1] += start;
            }
        }
        twice_num_hits = hit_array1.len() as i32;
        if twice_num_hits < 2 {
            hit_array.clear();
            return 0;
        }
        hit_array.clear();
        hit_array.extend_from_slice(&hit_array1);
        word_index -= 1;
    }

    twice_num_hits
}

/// Port of NCBI `s_PHIGetShortPattern` (`phi_gapalign.c:320`).
pub fn s_phi_get_short_pattern(
    seq: &[u8],
    len: i32,
    start: &mut i32,
    end: &mut i32,
    pattern_blk: &PhiPatternSearchBlk,
) -> i16 {
    let pattern_items = &pattern_blk.one_word_items;
    let mask = pattern_items.match_mask;
    let mask_shift_plus_1 = (mask << 1) + 1;
    let mut prefix_matched_bit_pattern = 0;

    for &letter in seq.iter().take(len.max(0) as usize) {
        prefix_matched_bit_pattern = ((prefix_matched_bit_pattern << 1) | mask_shift_plus_1)
            & pattern_items.which_position[letter as usize];
    }

    let (right_one, right_mask_only) = _phi_get_right_one_bits(prefix_matched_bit_pattern, mask);
    *start = right_mask_only + 1;
    *end = right_one;
    0
}

/// Port of NCBI `s_PHIGetLongPattern` (`phi_gapalign.c:362`).
pub fn s_phi_get_long_pattern(
    seq: &[u8],
    len: i32,
    start: &mut i32,
    end: &mut i32,
    pattern_blk: &PhiPatternSearchBlk,
) {
    let Some(multiword_items) = pattern_blk.multi_word_items.as_ref() else {
        *start = 0;
        *end = 0;
        return;
    };
    let num_words = multiword_items.num_words.max(0);
    let words = num_words as usize;
    if words == 0 {
        *start = 0;
        *end = 0;
        return;
    }

    let mut mask = multiword_items.match_mask_l[..words.min(PHI_BUF_SIZE)].to_vec();
    mask.resize(words, 0);
    let mut prefix_matched_bit_pattern = vec![0; words];
    _phi_pattern_words_left_shift(&mut mask, 1, num_words);

    for &letter in seq.iter().take(len.max(0) as usize) {
        _phi_pattern_words_left_shift(&mut prefix_matched_bit_pattern, 0, num_words);
        _phi_pattern_words_bitwise_or(&mut prefix_matched_bit_pattern, &mask, num_words);
        let letter_bits = multiword_items
            .bit_pattern_by_letter
            .get(letter as usize)
            .map(|row| row.as_slice())
            .unwrap_or(&[]);
        for (left, right) in prefix_matched_bit_pattern
            .iter_mut()
            .zip(letter_bits.iter())
            .take(words)
        {
            *left &= *right;
        }
    }

    let mut matched = prefix_matched_bit_pattern.clone();
    _phi_pattern_words_bitwise_and(
        &mut matched,
        &prefix_matched_bit_pattern,
        &multiword_items.match_mask_l,
        num_words,
    );

    let mut right_mask_only = -1;
    let mut found = false;
    let mut found_word = 0usize;
    let mut found_bit = 0i32;
    'outer: for word_index in 0..words {
        for bit_index in 0..PHI_BITS_PACKED_PER_WORD {
            if ((matched[word_index] >> bit_index) & 1) == 1 {
                found = true;
                found_word = word_index;
                found_bit = bit_index;
                break 'outer;
            }
            if ((multiword_items.match_mask_l[word_index] >> bit_index) & 1) == 1 {
                right_mask_only = word_index as i32 * PHI_BITS_PACKED_PER_WORD + bit_index;
            }
        }
    }

    *start = right_mask_only + 1;
    if found {
        *end = found_word as i32 * PHI_BITS_PACKED_PER_WORD + found_bit;
    } else {
        *end = words as i32 * PHI_BITS_PACKED_PER_WORD;
    }
}

/// Port of NCBI `s_PHIGetExtraLongPattern` (`phi_gapalign.c:419`).
pub fn s_phi_get_extra_long_pattern(
    seq: &[u8],
    len: i32,
    hit_array: &mut Vec<i32>,
    pattern_blk: &PhiPatternSearchBlk,
) -> i16 {
    hit_array.clear();
    let Some(multiword_items) = pattern_blk.multi_word_items.as_ref() else {
        return -1;
    };
    let Some(extra_items) = multiword_items.extra_long_items.as_ref() else {
        return -1;
    };
    let num_words = multiword_items.num_words.max(0) as usize;
    if num_words == 0 {
        return -1;
    }

    let mut pos: usize;
    hit_array.push(extra_items.num_places_in_word[0]);
    let mut active_len = 1usize;

    for word_index in 1..num_words {
        let word_blk = PhiPatternSearchBlk {
            one_word_items: short_items_for_long_word(multiword_items, word_index, 0),
            ..PhiPatternSearchBlk::default()
        };
        let mut hit_array1 = Vec::new();
        pos = 0;
        let mut j = 0usize;
        while j < active_len {
            let last_offset = hit_array[j + word_index - 1];
            let scan_len = (len - last_offset).min(
                extra_items.spacing[word_index - 1] + extra_items.num_places_in_word[word_index],
            );
            let mut one_word_hit_array = Vec::new();
            let seq_start = last_offset.max(0) as usize;
            let twice_hits_one_word = _phi_blast_find_hits_short(
                &mut one_word_hit_array,
                seq.get(seq_start..).unwrap_or(&[]),
                scan_len.max(0),
                &word_blk,
            );
            for hit_index in (0..twice_hits_one_word.max(0) as usize).step_by(2) {
                for word_index2 in 0..word_index {
                    hit_array1.push(hit_array[j + word_index2]);
                }
                hit_array1
                    .push(hit_array1[pos + word_index - 1] + one_word_hit_array[hit_index] + 1);
                pos += word_index + 1;
            }
            j += word_index;
        }
        active_len = pos;
        hit_array.clear();
        hit_array.extend_from_slice(&hit_array1);
    }

    let mut j = 0usize;
    while j + num_words <= active_len {
        if hit_array[j + num_words - 1] == len {
            let selected = hit_array[j..j + num_words].to_vec();
            hit_array.clear();
            hit_array.extend_from_slice(&selected);
            return 0;
        }
        j += num_words;
    }

    -1
}

/// Port of NCBI `FindPatternHits` (`pattern.c:468`).
pub fn find_pattern_hits(
    hit_array: &mut Vec<i32>,
    seq: &[u8],
    len: i32,
    is_dna: u8,
    pattern_blk: &PhiPatternSearchBlk,
) -> i32 {
    match pattern_blk.flag_pattern_length {
        PatternType::OneWord => s_find_hits_short_head(hit_array, seq, 0, len, is_dna, pattern_blk),
        PatternType::MultiWord => s_find_hits_long(hit_array, seq, len, pattern_blk),
        PatternType::VeryLong => s_find_hits_very_long(hit_array, seq, len, is_dna, pattern_blk),
    }
}

/// Port of NCBI `PHIBlastScanSubject` (`phi_lookup.c:725`).
pub fn phi_blast_scan_subject(
    pattern_blk: &PhiPatternSearchBlk,
    subject: &[u8],
    offset_ptr: &mut i32,
    offset_pairs: &mut Vec<crate::extend::PhiInitialHit>,
    _array_size: i32,
    is_dna: u8,
) -> i32 {
    *offset_ptr = subject.len() as i32;

    let mut hit_array = Vec::with_capacity(PHI_MAX_HIT);
    let twice_num_hits = find_pattern_hits(
        &mut hit_array,
        subject,
        subject.len() as i32,
        is_dna,
        pattern_blk,
    );

    let mut count = 0i32;
    for hit in hit_array
        .chunks_exact(2)
        .take((twice_num_hits / 2).max(0) as usize)
    {
        offset_pairs.push(crate::extend::PhiInitialHit {
            subject_start: hit[1],
            subject_end: hit[0],
        });
        count += 1;
    }

    count
}

/// Port of NCBI internal `s_PHISaveInitialHit` (`phi_extend.c:42`).
pub fn s_phi_save_initial_hit(
    init_hitlist: &mut crate::extend::InitHitList,
    offset_pair: crate::extend::PhiInitialHit,
) -> i16 {
    if crate::extend::blast_save_initial_hit(
        init_hitlist,
        offset_pair.subject_start,
        offset_pair.subject_end,
        None,
    ) {
        0
    } else {
        -1
    }
}

/// Port of NCBI `PHIBlastWordFinder` (`phi_extend.c:53`).
pub fn phi_blast_word_finder(
    subject: &[u8],
    pattern_blk: &PhiPatternSearchBlk,
    max_hits: i32,
    init_hitlist: &mut crate::extend::InitHitList,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
    is_dna: u8,
) -> i16 {
    let mut totalhits = 0i32;
    let mut first_offset = 0i32;
    let last_offset = subject.len() as i32;
    let mut offset_pairs = Vec::new();

    while first_offset < last_offset {
        let base = offset_pairs.len();
        let hits = phi_blast_scan_subject(
            pattern_blk,
            subject,
            &mut first_offset,
            &mut offset_pairs,
            max_hits,
            is_dna,
        );
        totalhits += hits;

        for offset_pair in offset_pairs
            .iter()
            .skip(base)
            .take(hits.max(0) as usize)
            .copied()
        {
            let _ = s_phi_save_initial_hit(init_hitlist, offset_pair);
        }
    }

    crate::diagnostics::blast_ungapped_stats_update(ungapped_stats, totalhits, 0, 0);
    0
}

/// Port of NCBI `PHIGetPatternOccurrences` (`pattern.c:557`).
pub fn phi_get_pattern_occurrences(
    pattern_blk: &PhiPatternSearchBlk,
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    is_dna: u8,
    pattern_info: &mut SphiQueryInfo,
    query_length: i32,
) -> i32 {
    let mut hit_array = Vec::with_capacity(query_sequence.len().saturating_mul(2));
    for loc in locations {
        let from = loc.left.max(0) as usize;
        let to = loc.right.max(loc.left) as usize;
        let Some(sequence) =
            query_sequence.get(from..=to.min(query_sequence.len().saturating_sub(1)))
        else {
            continue;
        };
        hit_array.clear();
        let loc_length = loc.right - loc.left + 1;
        let twice_num_hits =
            find_pattern_hits(&mut hit_array, sequence, loc_length, is_dna, pattern_blk);
        for hit in hit_array
            .chunks_exact(2)
            .take((twice_num_hits / 2).max(0) as usize)
        {
            let offset = hit[1] + loc.left;
            let length = hit[0] - hit[1] + 1;
            if offset == 0 && length == query_length {
                return i32::MAX;
            }
            let _ = s_phi_blast_add_pattern_hit(Some(pattern_info), offset, length);
        }
    }

    pattern_info.num_patterns
}

#[cfg(test)]
mod tests {
    use super::*;

    fn protein_word_pattern(letters: &[u8]) -> ShortPatternItems {
        let mut which_position = [0; PHI_ASCII_SIZE];
        for (index, &letter) in letters.iter().enumerate() {
            which_position[letter as usize] |= 1 << index;
        }
        ShortPatternItems {
            match_mask: 1 << (letters.len() as i32 - 1),
            which_position,
            dna_items: None,
        }
    }

    #[test]
    fn phi_get_right_one_bits_and_len_match_c_rules() {
        assert_eq!(_phi_get_right_one_bits(0b1000, 0b1010), (3, 1));
        assert_eq!(s_len_of(0b1000, 0b1010), 2);
        assert_eq!(_phi_get_right_one_bits(0, 0b1010), (0, 3));
        assert_eq!(s_len_of(0, 0b1010), -3);
    }

    #[test]
    fn phi_pattern_word_bit_operations_match_c_rules() {
        let mut words = vec![(1 << (PHI_BITS_PACKED_PER_WORD - 1)) - 1, 0];
        _phi_pattern_words_left_shift(&mut words, 1, 2);
        assert_eq!(words, vec![(1 << PHI_BITS_PACKED_PER_WORD) - 1, 0]);
        _phi_pattern_words_left_shift(&mut words, 0, 2);
        assert_eq!(words, vec![(1 << PHI_BITS_PACKED_PER_WORD) - 2, 1]);

        let mut left = vec![0b0011, 0b0100];
        _phi_pattern_words_bitwise_or(&mut left, &[0b1100, 0b0010], 2);
        assert_eq!(left, vec![0b1111, 0b0110]);

        let mut result = vec![0, 0];
        assert_eq!(
            _phi_pattern_words_bitwise_and(&mut result, &[0b1010, 0], &[0b1100, 0], 2),
            1
        );
        assert_eq!(result, vec![0b1000, 0]);
    }

    #[test]
    fn phi_len_of_l_matches_multiword_mask_scan() {
        let s = [0, 0b1000];
        let mask = [0b1000, 0b0010];
        assert_eq!(s_len_of_l(&s, &mask, 2), 2);
        assert_eq!(s_len_of_l(&[0], &[0b1], 1), -1);
    }

    #[test]
    fn phi_blast_find_hits_short_finds_one_word_protein_pattern() {
        let pattern_blk = PhiPatternSearchBlk {
            one_word_items: protein_word_pattern(b"AB"),
            ..Default::default()
        };
        let mut hits = Vec::new();

        assert_eq!(
            _phi_blast_find_hits_short(&mut hits, b"ABAB", 4, &pattern_blk),
            4
        );
        assert_eq!(hits, vec![1, 0, 3, 2]);

        hits.clear();
        assert_eq!(
            s_find_hits_short_head(&mut hits, b"XXAB", 2, 2, 0, &pattern_blk),
            2
        );
        assert_eq!(hits, vec![1, 0]);
    }

    #[test]
    fn sphi_pattern_search_blk_new_builds_one_word_pattern() {
        let pattern_blk = sphi_pattern_search_blk_new("a-b", false, None).expect("pattern");
        assert_eq!(pattern_blk.flag_pattern_length, PatternType::OneWord);
        assert_eq!(pattern_blk.pattern.as_deref(), Some("A-B"));
        assert_eq!(pattern_blk.min_pattern_match_length, 2);

        let mut hits = Vec::new();
        let seq = crate::encoding::encode_ncbistdaa_sequence(b"XXAB");
        assert_eq!(
            find_pattern_hits(&mut hits, &seq, seq.len() as i32, 0, &pattern_blk),
            2
        );
        assert_eq!(hits, vec![3, 2]);
        assert!(sphi_pattern_search_blk_free(Some(pattern_blk)).is_none());
    }

    #[test]
    fn sphi_pattern_search_blk_new_builds_long_pattern() {
        let pattern = format!("A{}B", "X".repeat(PHI_BITS_PACKED_PER_WORD as usize - 1));
        let pattern_blk = sphi_pattern_search_blk_new(&pattern, false, None).expect("pattern");
        assert_eq!(pattern_blk.flag_pattern_length, PatternType::MultiWord);
        assert_eq!(
            pattern_blk
                .multi_word_items
                .as_ref()
                .expect("multi")
                .num_words,
            2
        );

        let mut seq_ascii = Vec::from(&b"A"[..]);
        seq_ascii.extend(std::iter::repeat(b'W').take(PHI_BITS_PACKED_PER_WORD as usize - 1));
        seq_ascii.push(b'B');
        let seq = crate::encoding::encode_ncbistdaa_sequence(&seq_ascii);
        let mut hits = Vec::new();
        assert_eq!(
            find_pattern_hits(&mut hits, &seq, seq.len() as i32, 0, &pattern_blk),
            2
        );
        assert_eq!(hits, vec![30, 0]);
    }

    #[test]
    fn sphi_pattern_search_blk_new_builds_very_long_variable_pattern() {
        let pattern_blk = sphi_pattern_search_blk_new("ABX(2,40)CD", false, None).expect("pattern");
        assert_eq!(pattern_blk.flag_pattern_length, PatternType::VeryLong);
        let extra = pattern_blk
            .multi_word_items
            .as_ref()
            .and_then(|items| items.extra_long_items.as_ref())
            .expect("extra long");
        assert!(extra.highest_place > 0);
        assert!(extra.num_places_in_word[0] > 0);
    }

    #[test]
    fn sphi_pattern_search_blk_new_initializes_dna_prefix_suffix() {
        let pattern_blk = sphi_pattern_search_blk_new("AC", true, None).expect("dna pattern");
        assert_eq!(pattern_blk.flag_pattern_length, PatternType::OneWord);
        let dna = pattern_blk
            .one_word_items
            .dna_items
            .as_ref()
            .expect("dna items");
        assert!(dna.dna_which_prefix_pos.iter().any(|&value| value != 0));
        assert!(dna.dna_which_suffix_pos.iter().any(|&value| value != 0));
    }

    #[test]
    fn find_hits_long_matches_multiword_pattern() {
        let mut long_items = LongPatternItems {
            num_words: 2,
            ..Default::default()
        };
        long_items.match_mask_l[1] = 1 << 0;
        long_items.bit_pattern_by_letter[b'A' as usize][0] = 1 << 0;
        long_items.bit_pattern_by_letter[b'X' as usize][0] = (1 << PHI_BITS_PACKED_PER_WORD) - 2;
        long_items.bit_pattern_by_letter[b'B' as usize][1] = 1 << 0;
        let pattern_blk = PhiPatternSearchBlk {
            flag_pattern_length: PatternType::MultiWord,
            multi_word_items: Some(long_items),
            ..Default::default()
        };
        let mut hits = Vec::new();
        let mut seq = Vec::from(&b"A"[..]);
        seq.extend(std::iter::repeat(b'X').take(29));
        seq.push(b'B');

        assert_eq!(
            find_pattern_hits(&mut hits, &seq, seq.len() as i32, 0, &pattern_blk),
            2
        );
        assert_eq!(hits, vec![30, 0]);
    }

    #[test]
    fn phi_get_short_pattern_reports_suffix_match_bounds() {
        let pattern_blk = PhiPatternSearchBlk {
            one_word_items: protein_word_pattern(b"AB"),
            ..Default::default()
        };
        let mut start = -1;
        let mut end = -1;

        assert_eq!(
            s_phi_get_short_pattern(b"AB", 2, &mut start, &mut end, &pattern_blk),
            0
        );
        assert_eq!((start, end), (0, 1));
    }

    #[test]
    fn phi_get_long_pattern_reports_multiword_suffix_bounds() {
        let mut long_items = LongPatternItems {
            num_words: 2,
            ..Default::default()
        };
        long_items.match_mask_l[1] = 1 << 0;
        long_items.bit_pattern_by_letter[b'A' as usize][0] = 1 << 0;
        long_items.bit_pattern_by_letter[b'X' as usize][0] = (1 << PHI_BITS_PACKED_PER_WORD) - 2;
        long_items.bit_pattern_by_letter[b'B' as usize][1] = 1 << 0;
        let pattern_blk = PhiPatternSearchBlk {
            flag_pattern_length: PatternType::MultiWord,
            multi_word_items: Some(long_items),
            ..Default::default()
        };
        let mut seq = Vec::from(&b"A"[..]);
        seq.extend(std::iter::repeat(b'X').take(29));
        seq.push(b'B');
        let mut start = -1;
        let mut end = -1;

        s_phi_get_long_pattern(&seq, seq.len() as i32, &mut start, &mut end, &pattern_blk);

        assert_eq!((start, end), (0, 30));
    }

    #[test]
    fn find_hits_very_long_extends_short_words() {
        let mut long_items = LongPatternItems {
            num_words: 2,
            sll: vec![[0; PHI_ASCII_SIZE]; 2],
            ..Default::default()
        };
        long_items.match_mask_l[0] = 1 << 1;
        long_items.match_mask_l[1] = 1 << 1;
        long_items.sll[0][b'A' as usize] = 1 << 0;
        long_items.sll[0][b'B' as usize] = 1 << 1;
        long_items.sll[1][b'C' as usize] = 1 << 0;
        long_items.sll[1][b'D' as usize] = 1 << 1;
        let mut extra = ExtraLongPatternItems::default();
        extra.num_places_in_word[0] = 2;
        extra.num_places_in_word[1] = 2;
        extra.which_most_specific = 0;
        long_items.extra_long_items = Some(extra);
        let pattern_blk = PhiPatternSearchBlk {
            flag_pattern_length: PatternType::VeryLong,
            multi_word_items: Some(long_items),
            ..Default::default()
        };
        let mut hits = Vec::new();

        assert_eq!(
            find_pattern_hits(&mut hits, b"ABCDXX", 6, 0, &pattern_blk),
            2
        );
        assert_eq!(hits, vec![3, 0]);
    }

    #[test]
    fn phi_get_extra_long_pattern_reports_word_offsets() {
        let mut long_items = LongPatternItems {
            num_words: 2,
            sll: vec![[0; PHI_ASCII_SIZE]; 2],
            ..Default::default()
        };
        long_items.match_mask_l[0] = 1 << 1;
        long_items.match_mask_l[1] = 1 << 1;
        long_items.sll[0][b'A' as usize] = 1 << 0;
        long_items.sll[0][b'B' as usize] = 1 << 1;
        long_items.sll[1][b'C' as usize] = 1 << 0;
        long_items.sll[1][b'D' as usize] = 1 << 1;
        let mut extra = ExtraLongPatternItems::default();
        extra.num_places_in_word[0] = 2;
        extra.num_places_in_word[1] = 2;
        extra.spacing[0] = 0;
        long_items.extra_long_items = Some(extra);
        let pattern_blk = PhiPatternSearchBlk {
            flag_pattern_length: PatternType::VeryLong,
            multi_word_items: Some(long_items),
            ..Default::default()
        };
        let mut offsets = Vec::new();

        assert_eq!(
            s_phi_get_extra_long_pattern(b"ABCD", 4, &mut offsets, &pattern_blk),
            0
        );
        assert_eq!(offsets, vec![2, 4]);
    }

    #[test]
    fn phi_get_pattern_occurrences_adds_hits_and_detects_full_query() {
        let pattern_blk = PhiPatternSearchBlk {
            flag_pattern_length: PatternType::OneWord,
            one_word_items: protein_word_pattern(b"AB"),
            ..Default::default()
        };
        let ranges = [crate::util::SSeqRange { left: 2, right: 5 }];
        let mut info = sphi_query_info_new().expect("query info");

        assert_eq!(
            phi_get_pattern_occurrences(&pattern_blk, b"XXABAB", &ranges, 0, &mut info, 6),
            2
        );
        assert_eq!(
            info.occurrences,
            vec![
                SphiPatternInfo {
                    offset: 2,
                    length: 2
                },
                SphiPatternInfo {
                    offset: 4,
                    length: 2
                }
            ]
        );

        let full = [crate::util::SSeqRange { left: 0, right: 1 }];
        let mut full_info = sphi_query_info_new().expect("query info");
        assert_eq!(
            phi_get_pattern_occurrences(&pattern_blk, b"AB", &full, 0, &mut full_info, 2),
            i32::MAX
        );
    }

    #[test]
    fn sphi_query_info_lifecycle_and_add_hit_match_c_shape() {
        let mut info = sphi_query_info_new().expect("query info");
        assert_eq!(info.allocated_size, 8);
        assert_eq!(info.num_patterns, 0);

        for i in 0..9 {
            assert_eq!(s_phi_blast_add_pattern_hit(Some(&mut info), i, i + 3), 0);
        }

        assert_eq!(info.num_patterns, 9);
        assert_eq!(info.allocated_size, 16);
        assert_eq!(
            info.occurrences[8],
            SphiPatternInfo {
                offset: 8,
                length: 11
            }
        );
        assert_eq!(s_phi_blast_add_pattern_hit(None, 0, 0), -1);
        info.probability = 0.25;
        info.pattern = Some("A-B".to_string());
        let copy = sphi_query_info_copy(Some(&info)).expect("copied query info");
        assert_eq!(copy, info);
        assert_ne!(copy.occurrences.as_ptr(), info.occurrences.as_ptr());
        assert_eq!(sphi_query_info_copy(None), None);
        assert!(sphi_query_info_free(Some(info)).is_none());
    }

    #[test]
    fn phi_blast_scan_subject_returns_subject_pattern_bounds() {
        let pattern_blk = sphi_pattern_search_blk_new("A-B", false, None).expect("pattern");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXABAB");
        let mut offset_ptr = 0;
        let mut pairs = Vec::new();

        assert_eq!(
            phi_blast_scan_subject(&pattern_blk, &subject, &mut offset_ptr, &mut pairs, 32, 0,),
            2
        );

        assert_eq!(offset_ptr, 6);
        assert_eq!(
            pairs,
            vec![
                crate::extend::PhiInitialHit {
                    subject_start: 2,
                    subject_end: 3,
                },
                crate::extend::PhiInitialHit {
                    subject_start: 4,
                    subject_end: 5,
                },
            ]
        );
    }

    #[test]
    fn phi_blast_word_finder_saves_hits_and_updates_stats() {
        let pattern_blk = sphi_pattern_search_blk_new("A-B", false, None).expect("pattern");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXABAB");
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            phi_blast_word_finder(
                &subject,
                &pattern_blk,
                32,
                &mut init_hitlist,
                Some(&mut stats),
                0,
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 2);
        assert_eq!(init_hitlist.hits[0].query_offset, 2);
        assert_eq!(init_hitlist.hits[0].subject_offset, 3);
        assert_eq!(init_hitlist.hits[1].query_offset, 4);
        assert_eq!(init_hitlist.hits[1].subject_offset, 5);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.num_seqs_lookup_hits, 1);
        assert_eq!(stats.init_extends, 0);
        assert_eq!(stats.good_init_extends, 0);
    }
}
