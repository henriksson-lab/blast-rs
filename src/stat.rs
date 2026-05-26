//! Rust equivalent of blast_stat.c — Karlin-Altschul statistics.
//! This is the mathematical core for computing E-values and bit scores.

/// Karlin-Altschul statistical parameters for one context.
#[derive(Debug, Clone, Default)]
pub struct KarlinBlk {
    pub lambda: f64, // Lambda parameter
    pub k: f64,      // K parameter
    pub log_k: f64,  // ln(K)
    pub h: f64,      // H (relative entropy)
    /// NCBI `BlastScoreBlk::round_down` (`blast_stat.c:1868`). When
    /// `true`, odd scores are rounded down to even before e-value and
    /// bit-score computation. Set by `scaled_nucl_gapped_kbp_lookup` for
    /// scoring systems whose table values are only valid at even
    /// scores (e.g. `(2,-3)`, `(2,-5)`, `(2,-7)`, `(3,-4)`).
    pub round_down: bool,
}

#[derive(Debug, Clone, Default)]
pub struct BlastScoreMatrix {
    pub nrows: usize,
    pub ncols: usize,
    pub data: Vec<Vec<i32>>,
}

impl BlastScoreMatrix {
    /// blast-rs: Owned score-matrix constructor; not a direct NCBI C port.
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self {
            nrows,
            ncols,
            data: vec![vec![0; ncols]; nrows],
        }
    }
}

/// Rust-owned `BlastScoreBlk` counterpart for score/matrix/statistics setup.
///
/// The C structure is a large owner of score matrices, Karlin arrays, Gumbel
/// parameters, and option-derived flags. This Rust type keeps those fields as
/// owned values and vectors so parameter setup can pass the same logical inputs
/// C passes through raw pointers.
#[derive(Debug, Clone)]
pub struct BlastScoreBlk {
    pub alphabet_code: i32,
    pub alphabet_size: usize,
    pub protein_alphabet: bool,
    pub matrix: BlastScoreMatrix,
    pub name: Option<String>,
    pub loscore: i32,
    pub hiscore: i32,
    pub penalty: i32,
    pub reward: i32,
    pub read_in_matrix: bool,
    pub scale_factor: f64,
    pub number_of_contexts: i32,
    pub sfp: Vec<Option<ScoreFreq>>,
    pub kbp: Vec<KarlinBlk>,
    pub kbp_std: Vec<KarlinBlk>,
    pub kbp_psi: Vec<KarlinBlk>,
    pub kbp_gap: Vec<KarlinBlk>,
    pub kbp_gap_std: Vec<KarlinBlk>,
    pub kbp_gap_psi: Vec<KarlinBlk>,
    pub kbp_ideal: Option<KarlinBlk>,
    pub gbp: Option<GumbelBlk>,
    pub ambiguous_res: Vec<u8>,
    pub matrix_only_scoring: bool,
    pub round_down: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlastResFreq {
    pub alphabet_code: i32,
    pub prob: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastResComp {
    pub alphabet_code: i32,
    pub comp: Vec<i32>,
}

/// Port of NCBI `BlastScoreBlkNew` (`blast_stat.c:884`).
pub fn blast_score_blk_new(alphabet: i32, number_of_contexts: i32) -> Option<BlastScoreBlk> {
    let alphabet_size = if alphabet == crate::encoding::BLASTNA_SEQ_CODE {
        crate::encoding::BLASTNA_SIZE
    } else {
        crate::encoding::BLASTAA_SIZE
    };
    let protein_alphabet = alphabet == crate::encoding::BLASTAA_SEQ_CODE;
    let context_count = number_of_contexts.max(0) as usize;
    Some(BlastScoreBlk {
        alphabet_code: alphabet,
        alphabet_size,
        protein_alphabet,
        matrix: BlastScoreMatrix::new(alphabet_size, alphabet_size),
        name: None,
        loscore: 0,
        hiscore: 0,
        penalty: 0,
        reward: 0,
        read_in_matrix: false,
        scale_factor: 1.0,
        number_of_contexts,
        sfp: vec![None; context_count],
        kbp: vec![KarlinBlk::default(); context_count],
        kbp_std: vec![KarlinBlk::default(); context_count],
        kbp_psi: vec![KarlinBlk::default(); context_count],
        kbp_gap: Vec::new(),
        kbp_gap_std: vec![KarlinBlk::default(); context_count],
        kbp_gap_psi: vec![KarlinBlk::default(); context_count],
        kbp_ideal: None,
        gbp: std::env::var_os("OLD_FSC")
            .is_none()
            .then(s_blast_gumbel_blk_new),
        ambiguous_res: Vec::new(),
        matrix_only_scoring: false,
        round_down: false,
    })
}

/// Port-shaped equivalent of NCBI `BlastScoreBlkCheck` (`blast_stat.c:853`).
pub fn blast_score_blk_check(sbp: Option<&BlastScoreBlk>) -> i32 {
    let Some(sbp) = sbp else {
        return -1;
    };
    if sbp.kbp.is_empty() || sbp.sfp.is_empty() {
        return 1;
    }
    let contexts = sbp.number_of_contexts.max(0) as usize;
    let found = (0..contexts).any(|index| {
        sbp.sfp.get(index).is_some_and(Option::is_some)
            || sbp.kbp.get(index).is_some_and(KarlinBlk::is_valid)
    });
    if found {
        0
    } else {
        1
    }
}

/// Port of NCBI `BlastScoreBlkFree` (`blast_stat.c:965`).
pub fn blast_score_blk_free(sbp: &mut Option<BlastScoreBlk>) -> Option<BlastScoreBlk> {
    if let Some(sbp) = sbp.as_mut() {
        sbp.name = None;
        sbp.matrix.data.clear();
        sbp.matrix.nrows = 0;
        sbp.matrix.ncols = 0;
        sbp.sfp.clear();
        sbp.kbp.clear();
        sbp.kbp_std.clear();
        sbp.kbp_psi.clear();
        sbp.kbp_gap.clear();
        sbp.kbp_gap_std.clear();
        sbp.kbp_gap_psi.clear();
        sbp.kbp_ideal = None;
        sbp.gbp = None;
        sbp.ambiguous_res.clear();
    }
    *sbp = None;
    None
}

/// Port of NCBI `Blast_ResFreqNew` (`blast_stat.c:1708`).
pub fn blast_res_freq_new(sbp: Option<&BlastScoreBlk>) -> Option<BlastResFreq> {
    let sbp = sbp?;
    Some(BlastResFreq {
        alphabet_code: sbp.alphabet_code,
        prob: vec![0.0; sbp.alphabet_size],
    })
}

/// Port of NCBI `Blast_ResFreqFree` (`blast_stat.c:1693`).
pub fn blast_res_freq_free(rfp: &mut Option<BlastResFreq>) -> Option<BlastResFreq> {
    if let Some(rfp) = rfp.as_mut() {
        rfp.prob.clear();
    }
    *rfp = None;
    None
}

/// Port of NCBI `Blast_ResFreqNormalize` (`blast_stat.c:1835`).
pub fn blast_res_freq_normalize(
    sbp: Option<&BlastScoreBlk>,
    rfp: Option<&mut BlastResFreq>,
    norm: f64,
) -> i16 {
    let (Some(sbp), Some(rfp)) = (sbp, rfp) else {
        return 1;
    };
    if norm == 0.0 {
        return 1;
    }
    let mut sum = 0.0;
    for index in 0..sbp.alphabet_size {
        let p = rfp.prob.get(index).copied().unwrap_or(0.0);
        if p < 0.0 {
            return 1;
        }
        sum += p;
    }
    if sum <= 0.0 {
        return 0;
    }
    for index in 0..sbp.alphabet_size {
        rfp.prob[index] = rfp.prob[index] / sum * norm;
    }
    0
}

/// Port of NCBI `Blast_GetStdAlphabet` (`blast_stat.c:1863`).
pub fn blast_get_std_alphabet(alphabet_code: i32, residues: &mut [u8]) -> i16 {
    let freqs = protein_std_freq_letters();
    if residues.len() < freqs.len() {
        return -2;
    }
    for (index, &(letter, _)) in freqs.iter().enumerate() {
        residues[index] = if alphabet_code == crate::encoding::BLASTAA_SEQ_CODE {
            crate::encoding::AMINOACID_TO_NCBISTDAA[letter.to_ascii_uppercase() as usize]
        } else {
            letter
        };
    }
    freqs.len() as i16
}

/// Port of NCBI `Blast_ResFreqStdComp` (`blast_stat.c:1888`).
pub fn blast_res_freq_std_comp(sbp: Option<&BlastScoreBlk>, rfp: Option<&mut BlastResFreq>) -> i16 {
    let (Some(sbp), Some(rfp)) = (sbp, rfp) else {
        return 1;
    };
    if sbp.protein_alphabet {
        let mut residues = [0u8; 20];
        let retval = blast_get_std_alphabet(sbp.alphabet_code, &mut residues);
        if retval < 1 {
            return retval;
        }
        for (index, &residue) in residues.iter().enumerate() {
            if (residue as usize) < rfp.prob.len() {
                rfp.prob[residue as usize] = protein_std_freq_letters()[index].1;
            }
        }
    } else {
        for index in 0..4.min(rfp.prob.len()) {
            rfp.prob[index] = 25.0;
        }
    }
    blast_res_freq_normalize(Some(sbp), Some(rfp), 1.0)
}

/// NCBI: BlastResCompNew (blast_stat.c:1952).
pub fn blast_res_comp_new(sbp: Option<&BlastScoreBlk>) -> Option<BlastResComp> {
    let sbp = sbp?;
    Some(BlastResComp {
        alphabet_code: sbp.alphabet_code,
        comp: vec![0; sbp.alphabet_size],
    })
}

/// NCBI: BlastResCompDestruct (blast_stat.c:1967).
pub fn blast_res_comp_destruct(rcp: &mut Option<BlastResComp>) -> Option<BlastResComp> {
    if let Some(rcp) = rcp.as_mut() {
        rcp.comp.clear();
    }
    *rcp = None;
    None
}

/// NCBI: BlastResCompStr (blast_stat.c:1984).
pub fn blast_res_comp_str(
    sbp: Option<&BlastScoreBlk>,
    rcp: Option<&mut BlastResComp>,
    string: Option<&[u8]>,
    length: i32,
) -> i16 {
    let (Some(sbp), Some(rcp), Some(string)) = (sbp, rcp, string) else {
        return 1;
    };
    if rcp.alphabet_code != sbp.alphabet_code {
        return 1;
    }
    rcp.comp.fill(0);
    let mask = if sbp.protein_alphabet { 0xff } else { 0x0f };
    for &letter in string.iter().take(length.max(0) as usize) {
        let index = (letter & mask) as usize;
        if index < rcp.comp.len() {
            rcp.comp[index] += 1;
        }
    }
    for &ambiguous in &sbp.ambiguous_res {
        let index = ambiguous as usize;
        if index < rcp.comp.len() {
            rcp.comp[index] = 0;
        }
    }
    0
}

/// Port of NCBI `Blast_ResFreqClr` (`blast_stat.c:2024`).
pub fn blast_res_freq_clr(sbp: Option<&BlastScoreBlk>, rfp: Option<&mut BlastResFreq>) -> i16 {
    let (Some(sbp), Some(rfp)) = (sbp, rfp) else {
        return 1;
    };
    for index in 0..sbp.alphabet_size {
        rfp.prob[index] = 0.0;
    }
    0
}

/// Port of NCBI `Blast_ResFreqResComp` (`blast_stat.c:2044`).
pub fn blast_res_freq_res_comp(
    sbp: Option<&BlastScoreBlk>,
    rfp: Option<&mut BlastResFreq>,
    rcp: Option<&BlastResComp>,
) -> i16 {
    let (Some(sbp), Some(rfp), Some(rcp)) = (sbp, rfp, rcp) else {
        return 1;
    };
    if rfp.alphabet_code != rcp.alphabet_code {
        return 1;
    }
    let sum: f64 = rcp
        .comp
        .iter()
        .take(sbp.alphabet_size)
        .map(|&count| count as f64)
        .sum();
    if sum == 0.0 {
        return blast_res_freq_clr(Some(sbp), Some(rfp));
    }
    for index in 0..sbp.alphabet_size {
        rfp.prob[index] = rcp.comp[index] as f64 / sum;
    }
    0
}

/// Port of NCBI `Blast_ResFreqString` (`blast_stat.c:2078`).
pub fn blast_res_freq_string(
    sbp: Option<&BlastScoreBlk>,
    rfp: Option<&mut BlastResFreq>,
    string: Option<&[u8]>,
    length: i32,
) -> i16 {
    let Some(sbp_ref) = sbp else {
        return 1;
    };
    let Some(mut rcp) = blast_res_comp_new(Some(sbp_ref)) else {
        return 1;
    };
    if blast_res_comp_str(Some(sbp_ref), Some(&mut rcp), string, length) != 0 {
        return 1;
    }
    blast_res_freq_res_comp(Some(sbp_ref), rfp, Some(&rcp))
}

/// blast-rs: Owned deep-copy helper corresponding to `s_BlastScoreBlk_Copy`;
/// not a direct NCBI C port in this source snapshot.
pub fn s_blast_score_blk_copy(src: &BlastScoreBlk) -> BlastScoreBlk {
    let mut matrix = BlastScoreMatrix::new(src.matrix.nrows, src.matrix.ncols);
    for (dst, src_row) in matrix.data.iter_mut().zip(src.matrix.data.iter()) {
        dst.clear();
        dst.extend_from_slice(src_row);
    }

    BlastScoreBlk {
        alphabet_code: src.alphabet_code,
        alphabet_size: src.alphabet_size,
        protein_alphabet: src.protein_alphabet,
        matrix,
        name: src.name.clone(),
        loscore: src.loscore,
        hiscore: src.hiscore,
        penalty: src.penalty,
        reward: src.reward,
        read_in_matrix: src.read_in_matrix,
        scale_factor: src.scale_factor,
        number_of_contexts: src.number_of_contexts,
        sfp: src.sfp.clone(),
        kbp: src.kbp.clone(),
        kbp_std: src.kbp_std.clone(),
        kbp_psi: src.kbp_psi.clone(),
        kbp_gap: src.kbp_gap.clone(),
        kbp_gap_std: src.kbp_gap_std.clone(),
        kbp_gap_psi: src.kbp_gap_psi.clone(),
        kbp_ideal: src.kbp_ideal.clone(),
        gbp: src.gbp.clone(),
        ambiguous_res: src.ambiguous_res.clone(),
        matrix_only_scoring: src.matrix_only_scoring,
        round_down: src.round_down,
    }
}

/// Port of NCBI `BlastScoreBlkNuclMatrixCreate` (`blast_stat.c:1060`).
pub fn blast_score_blk_nucl_matrix_create(sbp: &mut BlastScoreBlk) -> i16 {
    if sbp.alphabet_size != crate::encoding::BLASTNA_SIZE {
        return 1;
    }
    if sbp.matrix.nrows != crate::encoding::BLASTNA_SIZE
        || sbp.matrix.ncols != crate::encoding::BLASTNA_SIZE
    {
        sbp.matrix =
            BlastScoreMatrix::new(crate::encoding::BLASTNA_SIZE, crate::encoding::BLASTNA_SIZE);
    }

    let reward = sbp.reward;
    let penalty = sbp.penalty;
    let mut degeneracy = [0i32; crate::encoding::BLASTNA_SIZE + 1];
    for degen in degeneracy.iter_mut().take(4) {
        *degen = 1;
    }
    for index1 in 4..crate::encoding::BLASTNA_SIZE {
        let mut degen = 0;
        for index2 in 0..4 {
            if (crate::encoding::BLASTNA_TO_NCBI4NA[index1]
                & crate::encoding::BLASTNA_TO_NCBI4NA[index2])
                != 0
            {
                degen += 1;
            }
        }
        degeneracy[index1] = degen;
    }

    for row in &mut sbp.matrix.data {
        row.fill(0);
    }
    for index1 in 0..crate::encoding::BLASTNA_SIZE {
        for index2 in index1..crate::encoding::BLASTNA_SIZE {
            let score = if (crate::encoding::BLASTNA_TO_NCBI4NA[index1]
                & crate::encoding::BLASTNA_TO_NCBI4NA[index2])
                != 0
            {
                crate::math::blast_nint(
                    ((degeneracy[index2] - 1) * penalty + reward) as f64
                        / degeneracy[index2].max(1) as f64,
                ) as i32
            } else {
                penalty
            };
            sbp.matrix.data[index1][index2] = score;
            sbp.matrix.data[index2][index1] = score;
        }
    }
    let strand_sentinel_score = i32::MIN / 2;
    sbp.matrix.data[crate::encoding::BLASTNA_SIZE - 1].fill(strand_sentinel_score);
    for index in 0..crate::encoding::BLASTNA_SIZE {
        sbp.matrix.data[index][crate::encoding::BLASTNA_SIZE - 1] = strand_sentinel_score;
    }
    blast_score_blk_max_score_set(sbp)
}

/// Port-shaped Rust reader for NCBI `BlastScoreBlkNucleotideMatrixRead`.
pub fn blast_score_blk_nucleotide_matrix_read(
    sbp: &mut BlastScoreBlk,
    rows: &[[i32; crate::encoding::BLASTNA_SIZE]; crate::encoding::BLASTNA_SIZE],
) -> i16 {
    if sbp.alphabet_code != crate::encoding::BLASTNA_SEQ_CODE
        && sbp.alphabet_code != crate::encoding::NCBI4NA_SEQ_CODE
    {
        return 1;
    }
    sbp.alphabet_size = crate::encoding::BLASTNA_SIZE;
    sbp.protein_alphabet = false;
    sbp.name = None;
    sbp.matrix =
        BlastScoreMatrix::new(crate::encoding::BLASTNA_SIZE, crate::encoding::BLASTNA_SIZE);
    for (row_index, src) in rows.iter().enumerate() {
        for (col_index, &score) in src.iter().enumerate() {
            if !(BLAST_SCORE_MIN..=BLAST_SCORE_MAX).contains(&score) {
                return 2;
            }
            sbp.matrix.data[row_index][col_index] = score;
        }
    }
    let gap = crate::encoding::BLASTNA_SIZE - 1;
    for index in 0..crate::encoding::BLASTNA_SIZE {
        sbp.matrix.data[gap][index] = sbp.matrix.data[gap][index].min(BLAST_SCORE_MIN);
        sbp.matrix.data[index][gap] = sbp.matrix.data[index][gap].min(BLAST_SCORE_MIN);
    }
    sbp.read_in_matrix = true;
    sbp.matrix_only_scoring = true;
    blast_score_blk_max_score_set(sbp)
}

/// Port-shaped Rust reader for NCBI `BlastScoreBlkProteinMatrixRead`.
pub fn blast_score_blk_protein_matrix_read(sbp: &mut BlastScoreBlk, matrix_name: &str) -> i16 {
    let matrix = if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        &crate::matrix::BLOSUM45
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM50") {
        &crate::matrix::BLOSUM50
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        &crate::matrix::BLOSUM62
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        &crate::matrix::BLOSUM80
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM90") {
        &crate::matrix::BLOSUM90
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        &crate::matrix::PAM30
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        &crate::matrix::PAM70
    } else if matrix_name.eq_ignore_ascii_case("PAM250") {
        &crate::matrix::PAM250
    } else if matrix_name.eq_ignore_ascii_case("IDENTITY") {
        &crate::matrix::IDENTITY
    } else {
        return 1;
    };
    sbp.alphabet_code = crate::encoding::BLASTAA_SEQ_CODE;
    sbp.alphabet_size = crate::encoding::BLASTAA_SIZE;
    sbp.protein_alphabet = true;
    sbp.matrix =
        BlastScoreMatrix::new(crate::encoding::BLASTAA_SIZE, crate::encoding::BLASTAA_SIZE);
    for row_index in 0..crate::encoding::BLASTAA_SIZE {
        for col_index in 0..crate::encoding::BLASTAA_SIZE {
            let score = matrix[row_index][col_index];
            if !(BLAST_SCORE_MIN..=BLAST_SCORE_MAX).contains(&score) {
                return 2;
            }
            sbp.matrix.data[row_index][col_index] = score;
        }
    }
    sbp.name = Some(matrix_name.to_ascii_uppercase());
    sbp.read_in_matrix = true;
    sbp.matrix_only_scoring = false;
    sbp.kbp_ideal = Some(protein_ideal_ungapped_kbp_for_matrix(matrix_name));
    blast_score_blk_max_score_set(sbp)
}

/// NCBI: BlastScoreBlkProteinMatrixLoad (blast_stat.c:1539).
pub fn blast_score_blk_protein_matrix_load(sbp: &mut BlastScoreBlk) -> i16 {
    let matrix_name = sbp.name.as_deref().unwrap_or("BLOSUM62");
    let source = if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        &crate::matrix::BLOSUM45
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM50") {
        &crate::matrix::BLOSUM50
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        &crate::matrix::BLOSUM62
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        &crate::matrix::BLOSUM80
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM90") {
        &crate::matrix::BLOSUM90
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        &crate::matrix::PAM30
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        &crate::matrix::PAM70
    } else if matrix_name.eq_ignore_ascii_case("PAM250") {
        &crate::matrix::PAM250
    } else if matrix_name.eq_ignore_ascii_case("IDENTITY") {
        &crate::matrix::IDENTITY
    } else {
        return 1;
    };

    if sbp.alphabet_size != crate::encoding::BLASTAA_SIZE
        || sbp.matrix.ncols != crate::encoding::BLASTAA_SIZE
        || sbp.matrix.nrows != crate::encoding::BLASTAA_SIZE
    {
        return 1;
    }

    sbp.matrix.data =
        vec![vec![BLAST_SCORE_MIN; crate::encoding::BLASTAA_SIZE]; crate::encoding::BLASTAA_SIZE];
    for i in 0..crate::encoding::BLASTAA_SIZE {
        for j in 0..crate::encoding::BLASTAA_SIZE {
            if i == crate::encoding::NCBISTDAA_U as usize
                || i == crate::encoding::NCBISTDAA_O as usize
                || i == crate::encoding::NCBISTDAA_GAP as usize
                || j == crate::encoding::NCBISTDAA_U as usize
                || j == crate::encoding::NCBISTDAA_O as usize
                || j == crate::encoding::NCBISTDAA_GAP as usize
            {
                continue;
            }
            sbp.matrix.data[i][j] = source[i][j];
        }
    }

    let x = crate::encoding::NCBISTDAA_X as usize;
    let u = crate::encoding::NCBISTDAA_U as usize;
    let o = crate::encoding::NCBISTDAA_O as usize;
    let c = crate::encoding::NCBISTDAA_C as usize;
    for i in 0..crate::encoding::BLASTAA_SIZE {
        sbp.matrix.data[u][i] = sbp.matrix.data[c][i];
        sbp.matrix.data[i][u] = sbp.matrix.data[i][c];
        sbp.matrix.data[o][i] = sbp.matrix.data[x][i];
        sbp.matrix.data[i][o] = sbp.matrix.data[i][x];
    }

    sbp.alphabet_code = crate::encoding::BLASTAA_SEQ_CODE;
    sbp.protein_alphabet = true;
    sbp.read_in_matrix = true;
    sbp.matrix_only_scoring = false;
    sbp.kbp_ideal = Some(protein_ideal_ungapped_kbp_for_matrix(matrix_name));
    blast_score_blk_max_score_set(sbp)
}

/// NCBI: Blast_ScoreBlkMatrixFill (blast_stat.c:1599).
pub fn blast_score_blk_matrix_fill(sbp: &mut BlastScoreBlk) -> i16 {
    let status = if sbp.alphabet_code == crate::encoding::BLASTNA_SEQ_CODE
        || sbp.alphabet_code == crate::encoding::NCBI4NA_SEQ_CODE
    {
        if sbp.read_in_matrix {
            return 1;
        }
        blast_score_blk_nucl_matrix_create(sbp)
    } else {
        blast_score_blk_protein_matrix_load(sbp)
    };
    if status != 0 {
        return status;
    }
    blast_score_blk_max_score_set(sbp)
}

/// NCBI: Blast_ScoreBlkKbpIdealCalc (blast_stat.c:2832).
pub fn blast_score_blk_kbp_ideal_calc(sbp: Option<&mut BlastScoreBlk>) -> i16 {
    let Some(sbp) = sbp else {
        return 1;
    };
    if !sbp.protein_alphabet && sbp.loscore >= sbp.hiscore {
        let status = blast_score_blk_nucl_matrix_create(sbp);
        if status != 0 {
            return status;
        }
    }
    let Some(mut stdrfp) = blast_res_freq_new(Some(sbp)) else {
        return 1;
    };
    if blast_res_freq_std_comp(Some(sbp), Some(&mut stdrfp)) != 0 {
        return 1;
    }
    let Some(mut sfp) = blast_score_freq_new(sbp.loscore, sbp.hiscore) else {
        return 1;
    };
    if blast_score_freq_calc(Some(sbp), Some(&mut sfp), Some(&stdrfp), Some(&stdrfp)) != 0 {
        return 1;
    }
    let mut ideal = blast_karlin_blk_new();
    if blast_karlin_blk_ungapped_calc(Some(&mut ideal), Some(&sfp)) != 0 {
        return 1;
    }
    sbp.kbp_ideal = Some(ideal);
    0
}

/// Port of NCBI `Blast_KarlinBlkCopy` (`blast_stat.c:2871`).
pub fn blast_karlin_blk_copy(kbp_to: Option<&mut KarlinBlk>, kbp_from: Option<&KarlinBlk>) -> i16 {
    let (Some(kbp_to), Some(kbp_from)) = (kbp_to, kbp_from) else {
        return -1;
    };
    *kbp_to = kbp_from.clone();
    0
}

/// Port of NCBI `BlastScoreBlkMaxScoreSet` (`blast_stat.c:1500`).
///
/// NCBI explicitly skips `BLAST_SCORE_MIN`/`BLAST_SCORE_MAX` sentinels
/// (used for gap rows/columns and other "impossible" entries) when
/// computing the loscore/hiscore range. The clamp at the end also
/// restores those sentinels if the loop produced no real scores.
pub fn blast_score_blk_max_score_set(sbp: &mut BlastScoreBlk) -> i16 {
    if sbp.matrix.data.is_empty() {
        return 1;
    }
    let mut loscore = BLAST_SCORE_MAX;
    let mut hiscore = BLAST_SCORE_MIN;
    for row in sbp.matrix.data.iter().take(sbp.alphabet_size) {
        for &score in row.iter().take(sbp.alphabet_size) {
            // NCBI `blast_stat.c:1513`: skip sentinel entries (gap rows,
            // etc.) so they don't pollute the real-score range.
            if score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX {
                continue;
            }
            if loscore > score {
                loscore = score;
            }
            if hiscore < score {
                hiscore = score;
            }
        }
    }
    // NCBI `blast_stat.c:1525`: clamp if no real scores were observed.
    if loscore == BLAST_SCORE_MAX && hiscore == BLAST_SCORE_MIN {
        loscore = BLAST_SCORE_MIN;
        hiscore = BLAST_SCORE_MAX;
    }
    if loscore < BLAST_SCORE_MIN {
        loscore = BLAST_SCORE_MIN;
    }
    if hiscore > BLAST_SCORE_MAX {
        hiscore = BLAST_SCORE_MAX;
    }
    sbp.loscore = loscore;
    sbp.hiscore = hiscore;
    0
}

/// blast-rs: Karlin-block validity helper matching the C predicate shape; not
/// a direct NCBI C port in this source snapshot.
pub fn s_blast_karlin_blk_is_valid(kbp: Option<&KarlinBlk>) -> bool {
    kbp.is_some_and(KarlinBlk::is_valid)
}

/// blast-rs: Shared lookup implementation for valid Karlin blocks; not a
/// direct NCBI C port.
pub fn s_blast_find_valid_karlin_blk<'a>(
    kbp_in: &'a [KarlinBlk],
    query_info: &crate::queryinfo::QueryInfo,
) -> Result<&'a KarlinBlk, i16> {
    for (context_index, context) in query_info.contexts.iter().enumerate() {
        if !context.is_valid {
            continue;
        }
        if let Some(kbp) = kbp_in.get(context_index) {
            if s_blast_karlin_blk_is_valid(Some(kbp)) {
                return Ok(kbp);
            }
        }
    }
    Err(crate::diagnostics::BLASTERR_NOVALIDKARLINALTSCHUL)
}

/// NCBI: s_BlastFindValidKarlinBlk (blast_parameters.c:62).
/// naming: `_c` suffix marks the pointer-shaped compatibility wrapper.
pub fn s_blast_find_valid_karlin_blk_c<'a>(
    kbp_in: &'a [KarlinBlk],
    query_info: Option<&crate::queryinfo::QueryInfo>,
    kbp_ret: Option<&mut Option<&'a KarlinBlk>>,
) -> i16 {
    let (Some(query_info), Some(kbp_ret)) = (query_info, kbp_ret) else {
        return crate::diagnostics::BLASTERR_NOVALIDKARLINALTSCHUL;
    };
    match s_blast_find_valid_karlin_blk(kbp_in, query_info) {
        Ok(kbp) => {
            *kbp_ret = Some(kbp);
            0
        }
        Err(status) => {
            *kbp_ret = None;
            status
        }
    }
}

/// blast-rs: Shared lookup implementation for smallest valid lambda; not a
/// direct NCBI C port.
pub fn s_blast_find_smallest_lambda<'a>(
    kbp_in: &'a [KarlinBlk],
    query_info: &crate::queryinfo::QueryInfo,
) -> Option<(f64, &'a KarlinBlk)> {
    query_info
        .contexts
        .iter()
        .enumerate()
        .filter(|(_, context)| context.is_valid)
        .filter_map(|(context_index, _)| kbp_in.get(context_index))
        .filter(|kbp| s_blast_karlin_blk_is_valid(Some(kbp)))
        .min_by(|left, right| left.lambda.total_cmp(&right.lambda))
        .map(|kbp| (kbp.lambda, kbp))
}

/// NCBI: s_BlastFindSmallestLambda (blast_parameters.c:92).
/// naming: `_c` suffix marks the pointer-shaped compatibility wrapper.
pub fn s_blast_find_smallest_lambda_c<'a>(
    kbp_in: &'a [KarlinBlk],
    query_info: Option<&crate::queryinfo::QueryInfo>,
    kbp_out: Option<&mut Option<&'a KarlinBlk>>,
) -> f64 {
    let Some(query_info) = query_info else {
        if let Some(kbp_out) = kbp_out {
            *kbp_out = None;
        }
        return i32::MAX as f64;
    };
    let result = s_blast_find_smallest_lambda(kbp_in, query_info);
    if let Some(kbp_out) = kbp_out {
        *kbp_out = result.map(|(_, kbp)| kbp);
    }
    result.map(|(lambda, _)| lambda).unwrap_or(i32::MAX as f64)
}

/// NCBI: s_BlastGumbelBlkFree (blast_stat.c:846).
pub fn s_blast_gumbel_blk_free(gbp: &mut Option<GumbelBlk>) -> Option<GumbelBlk> {
    *gbp = None;
    None
}

/// NCBI: s_BlastGumbelBlkNew (blast_stat.c:838).
pub fn s_blast_gumbel_blk_new() -> GumbelBlk {
    GumbelBlk {
        lambda: 0.0,
        a: 0.0,
        b: 0.0,
        alpha: 0.0,
        beta: 0.0,
        sigma: 0.0,
        tau: 0.0,
        db_length: 0,
    }
}

/// blast-rs: Local equivalent of NCBI's `gbp->filled` flag; not a direct NCBI C port.
pub fn gumbel_blk_is_filled(gbp: &GumbelBlk) -> bool {
    gbp.lambda.is_finite() && gbp.lambda > 0.0
}

/// NCBI: Blast_KarlinBlkFree (blast_stat.c:956).
pub fn blast_karlin_blk_free(kbp: &mut Option<KarlinBlk>) -> Option<KarlinBlk> {
    *kbp = None;
    None
}

/// NCBI: BLAST_ScoreSetAmbigRes (blast_stat.c:1012).
pub fn blast_score_set_ambig_res(sbp: Option<&mut BlastScoreBlk>, ambiguous_res: u8) -> i16 {
    let Some(sbp) = sbp else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    let byte = ambiguous_res.to_ascii_uppercase();
    if byte as usize >= 128 {
        return crate::util::BLASTERR_INVALIDPARAM;
    }
    let index = byte as usize;
    let encoded = if sbp.alphabet_code == crate::encoding::BLASTAA_SEQ_CODE {
        crate::encoding::AMINOACID_TO_NCBISTDAA[index]
    } else if sbp.alphabet_code == crate::encoding::BLASTNA_SEQ_CODE {
        crate::encoding::IUPACNA_TO_BLASTNA[index]
    } else if sbp.alphabet_code == crate::encoding::NCBI4NA_SEQ_CODE {
        crate::encoding::IUPACNA_TO_NCBI4NA[index]
    } else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    sbp.ambiguous_res.push(encoded);
    0
}

/// blast-rs: Shared round-down helper for Karlin score conversion; not a
/// direct NCBI C port.
#[inline]
fn round_down_score(score: i32, round_down: bool) -> i32 {
    if round_down {
        score & !1
    } else {
        score
    }
}

impl KarlinBlk {
    /// blast-rs: Karlin-block validity method; not a direct NCBI C port.
    pub fn is_valid(&self) -> bool {
        self.lambda > 0.0 && self.k > 0.0 && self.h > 0.0
    }

    /// Convert a raw score to a bit score.
    /// blast-rs: Formula helper extracted from bit-score calculation; not a
    /// direct NCBI C port:
    /// `(score * Lambda - logK) / NCBIMATH_LN2`. Note that unlike the
    /// E-value path, the bit-score formula does NOT apply `round_down`
    /// even-masking — `Blast_HSPListGetEvalues` (`blast_hits.c:1864-1869`)
    /// applies `score &= ~1`, but `Blast_HSPListGetBitScores` (`:1907`)
    /// only has a commented-out `#if 0` assertion for it.
    pub fn raw_to_bit(&self, raw_score: i32) -> f64 {
        (self.lambda * raw_score as f64 - self.log_k) / crate::math::NCBIMATH_LN2
    }

    /// Convert a raw score and search space to an E-value.
    /// blast-rs: Karlin score-to-E-value method; not a direct NCBI C port:
    /// returns `-1.0` when the Karlin block is degenerate, otherwise
    /// `searchsp * exp(-lambda * S + logK)` (NCBI `:4170` uses `logK`
    /// inside the exponent for numerical stability at large magnitudes).
    /// Applies `round_down` before the formula (`blast_hits.c:1864-1869`).
    pub fn raw_to_evalue(&self, raw_score: i32, search_space: f64) -> f64 {
        if self.lambda < 0.0 || self.k < 0.0 || self.h < 0.0 {
            return -1.0;
        }
        let score = round_down_score(raw_score, self.round_down);
        search_space * (-self.lambda * score as f64 + self.log_k).exp()
    }

    /// Convert an E-value to the minimum raw score needed.
    /// blast-rs: Karlin E-value-to-score method; not a direct NCBI C port. The
    /// e-value is clamped to `K_SMALL_FLOAT` before conversion to avoid
    /// floating-point exceptions on extremely tight cutoffs. Returns
    /// `BLAST_SCORE_MIN` when lambda/K/H are invalid (NCBI `:4054-4057`).
    pub fn evalue_to_raw(&self, evalue: f64, search_space: f64) -> i32 {
        // NCBI checks `Lambda < 0. || K < 0. || H < 0.0`. We don't carry
        // `h` here at the call site, but zero/negative lambda or K are
        // definitive invalid markers; match NCBI's return value.
        if self.lambda < 0.0 || self.k < 0.0 || self.h < 0.0 {
            return BLAST_SCORE_MIN;
        }
        let denom = search_space * self.k;
        if denom <= 0.0 || self.lambda <= 0.0 {
            return BLAST_SCORE_MIN;
        }
        let evalue = evalue.max(K_SMALL_FLOAT);
        let score = -(evalue / denom).ln() / self.lambda;
        score.ceil() as i32
    }
}

/// NCBI: BLAST_KarlinStoE_simple (blast_stat.c:4157).
pub fn blast_karlin_sto_e_simple(score: i32, kbp: Option<&KarlinBlk>, searchsp: i64) -> f64 {
    let Some(kbp) = kbp else {
        return -1.0;
    };
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return -1.0;
    }
    searchsp as f64 * (-kbp.lambda * score as f64 + kbp.log_k).exp()
}

/// NCBI: BlastKarlinEtoS_simple (blast_stat.c:4040).
pub fn blast_karlin_eto_s_simple(evalue: f64, kbp: Option<&KarlinBlk>, searchsp: i64) -> i32 {
    let Some(kbp) = kbp else {
        return BLAST_SCORE_MIN;
    };
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return BLAST_SCORE_MIN;
    }
    let evalue = evalue.max(K_SMALL_FLOAT);
    ((kbp.k * searchsp as f64 / evalue).ln() / kbp.lambda).ceil() as i32
}

/// 1-1 port of `Blast_HSPListGetBitScores` (`blast_hits.c:1907`) for
/// `hspstream::HspList`.
pub fn blast_hsp_list_get_bit_scores(
    hsp_list: Option<&mut crate::hspstream::HspList>,
    _gapped_calculation: bool,
    kbp: &[KarlinBlk],
) -> i32 {
    let Some(hsp_list) = hsp_list else { return 1 };
    for hsp in &mut hsp_list.hsps {
        let Some(ctx_kbp) = kbp.get(hsp.context.max(0) as usize) else {
            return 1;
        };
        hsp.bit_score =
            (hsp.score as f64 * ctx_kbp.lambda - ctx_kbp.log_k) / crate::math::NCBIMATH_LN2;
    }
    0
}

/// `BLAST_GAP_PROB` (`blast_parameters.h:66`): gap probability for
/// ungapped sum-statistics (0.5).
pub const BLAST_GAP_PROB: f64 = 0.5;
/// `BLAST_GAP_PROB_GAPPED` (`blast_parameters.h:67`): gap probability
/// for gapped sum-statistics (1.0).
pub const BLAST_GAP_PROB_GAPPED: f64 = 1.0;
/// Port of NCBI `BLAST_GAP_DECAY_RATE` (`blast_parameters.h:68`):
/// default gap decay rate for ungapped search.
pub const BLAST_GAP_DECAY_RATE: f64 = 0.5;
/// Port of NCBI `BLAST_GAP_DECAY_RATE_GAPPED` (`blast_parameters.h:69`):
/// default gap decay rate for gapped search.
pub const BLAST_GAP_DECAY_RATE_GAPPED: f64 = 0.1;
/// `BLAST_GAP_SIZE` (`blast_parameters.h:70`): default gap size (nt
/// distance) used by HSP linking.
pub const BLAST_GAP_SIZE: i32 = 40;
/// `BLAST_OVERLAP_SIZE` (`blast_parameters.h:71`): default overlap
/// allowed between linked HSPs.
pub const BLAST_OVERLAP_SIZE: i32 = 9;
/// `RESTRICTED_ALIGNMENT_WORST_EVALUE` (`blast_parameters.h:193`):
/// composition-adjusted restricted alignment cutoff.
pub const RESTRICTED_ALIGNMENT_WORST_EVALUE: f64 = 10.0;

/// Port of NCBI `CUTOFF_E_BLASTN` (`blast_parameters.h:76`): default
/// e-value used by the preliminary ungapped cutoff for nucleotide
/// searches.
pub const CUTOFF_E_BLASTN: f64 = 0.05;
/// Port of NCBI `CUTOFF_E_BLASTP` (`blast_parameters.h:77`).
pub const CUTOFF_E_BLASTP: f64 = 1.0e-300;
/// Port of NCBI `CUTOFF_E_BLASTX` (`blast_parameters.h:78`).
pub const CUTOFF_E_BLASTX: f64 = 1.0;
/// Port of NCBI `CUTOFF_E_TBLASTN` (`blast_parameters.h:79`).
pub const CUTOFF_E_TBLASTN: f64 = 1.0;
/// Port of NCBI `CUTOFF_E_TBLASTX` (`blast_parameters.h:80`).
pub const CUTOFF_E_TBLASTX: f64 = 1.0e-300;

// X-dropoff defaults (in bits). NCBI `blast_options.h:122-148`.
/// `BLAST_UNGAPPED_X_DROPOFF_PROT` — ungapped dropoff for protein.
pub const BLAST_UNGAPPED_X_DROPOFF_PROT: i32 = 7;
/// `BLAST_UNGAPPED_X_DROPOFF_NUCL` — ungapped dropoff for nucleotide.
pub const BLAST_UNGAPPED_X_DROPOFF_NUCL: i32 = 20;
/// `BLAST_GAP_X_DROPOFF_PROT` — default preliminary-gapped dropoff for protein.
pub const BLAST_GAP_X_DROPOFF_PROT: i32 = 15;
/// `BLAST_GAP_X_DROPOFF_NUCL` — default preliminary-gapped dropoff for nucleotide.
pub const BLAST_GAP_X_DROPOFF_NUCL: i32 = 30;
/// `BLAST_GAP_X_DROPOFF_GREEDY` — default dropoff for greedy megablast.
pub const BLAST_GAP_X_DROPOFF_GREEDY: i32 = 25;
/// `BLAST_GAP_X_DROPOFF_FINAL_PROT` — final-gapped dropoff for protein.
pub const BLAST_GAP_X_DROPOFF_FINAL_PROT: i32 = 25;
/// `BLAST_GAP_X_DROPOFF_FINAL_NUCL` — final-gapped dropoff for nucleotide.
pub const BLAST_GAP_X_DROPOFF_FINAL_NUCL: i32 = 100;
/// `BLAST_GAP_TRIGGER_PROT` — protein gap-trigger bit threshold (`blast_options.h:137`).
pub const BLAST_GAP_TRIGGER_PROT: f64 = 22.0;
/// `BLAST_GAP_TRIGGER_NUCL` — nucleotide gap-trigger bit threshold (`blast_options.h:140`).
pub const BLAST_GAP_TRIGGER_NUCL: f64 = 27.0;

// Default word sizes (NCBI `blast_options.h:66-73`).
/// `BLAST_WORDSIZE_PROT` — default word size for protein searches (3).
pub const BLAST_WORDSIZE_PROT: i32 = 3;
/// `BLAST_WORDSIZE_NUCL` — default word size for blastn (11).
pub const BLAST_WORDSIZE_NUCL: i32 = 11;
/// `BLAST_WORDSIZE_MEGABLAST` — default word size for contiguous megablast (28).
pub const BLAST_WORDSIZE_MEGABLAST: i32 = 28;
/// `BLAST_WORDSIZE_MAPPER` — default word size for the magicblast mapper (18).
pub const BLAST_WORDSIZE_MAPPER: i32 = 18;

// Default gap costs (NCBI `blast_options.h:84-98`).
/// `BLAST_GAP_OPEN_PROT` — protein gap-open (11).
pub const BLAST_GAP_OPEN_PROT: i32 = 11;
/// `BLAST_GAP_OPEN_NUCL` — blastn gap-open (5).
pub const BLAST_GAP_OPEN_NUCL: i32 = 5;
/// `BLAST_GAP_OPEN_MEGABLAST` — megablast gap-open (0).
pub const BLAST_GAP_OPEN_MEGABLAST: i32 = 0;
/// `BLAST_GAP_EXTN_PROT` — protein gap-extend (1).
pub const BLAST_GAP_EXTN_PROT: i32 = 1;
/// `BLAST_GAP_EXTN_NUCL` — blastn gap-extend (2).
pub const BLAST_GAP_EXTN_NUCL: i32 = 2;
/// `BLAST_GAP_EXTN_MEGABLAST` — megablast gap-extend (0).
pub const BLAST_GAP_EXTN_MEGABLAST: i32 = 0;

// Default match/mismatch scores (NCBI `blast_options.h:151-152`).
/// `BLAST_PENALTY` — default nucleotide mismatch score (-3).
pub const BLAST_PENALTY: i32 = -3;
/// `BLAST_REWARD` — default nucleotide match score (1).
pub const BLAST_REWARD: i32 = 1;

// Default neighboring-word thresholds (NCBI `blast_options.h:104-116`).
/// `BLAST_WORD_THRESHOLD_BLASTP` (11).
pub const BLAST_WORD_THRESHOLD_BLASTP: f64 = 11.0;
/// `BLAST_WORD_THRESHOLD_BLASTP_FAST` — word-size 5 threshold for blastp/x-fast (20).
pub const BLAST_WORD_THRESHOLD_BLASTP_FAST: f64 = 20.0;
/// `BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_6` (21).
pub const BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_6: f64 = 21.0;
/// `BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_7` (20.25).
pub const BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_7: f64 = 20.25;
/// `BLAST_WORD_THRESHOLD_BLASTN` (0 = no threshold).
pub const BLAST_WORD_THRESHOLD_BLASTN: f64 = 0.0;
/// `BLAST_WORD_THRESHOLD_BLASTX` (12).
pub const BLAST_WORD_THRESHOLD_BLASTX: f64 = 12.0;
/// `BLAST_WORD_THRESHOLD_TBLASTN` (13).
pub const BLAST_WORD_THRESHOLD_TBLASTN: f64 = 13.0;
/// `BLAST_WORD_THRESHOLD_TBLASTX` (13).
pub const BLAST_WORD_THRESHOLD_TBLASTX: f64 = 13.0;
/// `BLAST_WORD_THRESHOLD_MEGABLAST` (0).
pub const BLAST_WORD_THRESHOLD_MEGABLAST: f64 = 0.0;

/// `BLAST_SCAN_RANGE_NUCL` — default scan range for blastn (0, `blast_options.h:63`).
pub const BLAST_SCAN_RANGE_NUCL: i32 = 0;
/// `BLAST_GAP_X_DROPOFF_TBLASTX` (0, `blast_options.h:134`).
pub const BLAST_GAP_X_DROPOFF_TBLASTX: i32 = 0;
/// `BLAST_GAP_X_DROPOFF_FINAL_TBLASTX` (0, `blast_options.h:148`).
pub const BLAST_GAP_X_DROPOFF_FINAL_TBLASTX: i32 = 0;

/// `BLAST_DEFAULT_MATRIX` — default scoring matrix for protein searches
/// (`blast_options.h:77`).
pub const BLAST_DEFAULT_MATRIX: &str = "BLOSUM62";

/// `BLAST_PENALTY_MAPPER` (`blast_options.h:154`): mapper nucleotide mismatch.
pub const BLAST_PENALTY_MAPPER: i32 = -4;
/// `BLAST_REWARD_MAPPER` (`blast_options.h:155`): mapper nucleotide match.
pub const BLAST_REWARD_MAPPER: i32 = 1;
/// `BLAST_GAP_OPEN_MAPPER` (`blast_options.h:89`).
pub const BLAST_GAP_OPEN_MAPPER: i32 = 0;
/// `BLAST_GAP_EXTN_MAPPER` (`blast_options.h:98`).
pub const BLAST_GAP_EXTN_MAPPER: i32 = 4;

/// `PSI_INCLUSION_ETHRESH` (`blast_options.h:163`): PSI-BLAST inclusion
/// threshold (e-value below which a hit is used for the next iteration).
pub const PSI_INCLUSION_ETHRESH: f64 = 0.002;
/// `PSI_PSEUDO_COUNT_CONST` (`blast_options.h:164`): PSI-BLAST
/// pseudo-count constant (0 disables).
pub const PSI_PSEUDO_COUNT_CONST: i32 = 0;
/// `DELTA_INCLUSION_ETHRESH` (`blast_options.h:165`): DELTA-BLAST inclusion.
pub const DELTA_INCLUSION_ETHRESH: f64 = 0.05;

/// `BLAST_GENETIC_CODE` (`blast_options.h:168`): the standard genetic
/// code (NCBI table #1).
pub const BLAST_GENETIC_CODE: i32 = 1;
/// `MAX_DB_WORD_COUNT_MAPPER` (`blast_options.h:174`): word cap for
/// the magicblast mapper.
pub const MAX_DB_WORD_COUNT_MAPPER: i32 = 30;

/// `HSP_MAX_WINDOW` — sliding-window size used by
/// `BlastGetOffsetsForGappedAlignment` / `BlastGetStartForGappedAlignment`
/// to pick a high-scoring seed offset near an ungapped HSP
/// (`blast_gapalign_priv.h:120`).
pub const HSP_MAX_WINDOW: usize = 11;

/// NCBI `MININT` (`blast_gapalign.c:58`): `INT4_MIN/2`. Used as the
/// sentinel for impossible DP cells; halving avoids underflow when
/// small mismatches are added to this value during the recurrence.
pub const MININT: i32 = i32::MIN / 2;

/// `BLAST_EXPECT_VALUE` — default e-value threshold
/// (`blast_options.h:158`).
pub const BLAST_EXPECT_VALUE: f64 = 10.0;
/// `BLAST_HITLIST_SIZE` — default number of database sequences to save hits
/// for (`blast_options.h:160`).
pub const BLAST_HITLIST_SIZE: usize = 500;

// Default two-hit window sizes (NCBI `blast_options.h:57-61`).
/// `BLAST_WINDOW_SIZE_PROT` — protein two-hit window (40).
pub const BLAST_WINDOW_SIZE_PROT: i32 = 40;
/// `BLAST_WINDOW_SIZE_NUCL` — blastn (non-MB) window (0 = disabled).
pub const BLAST_WINDOW_SIZE_NUCL: i32 = 0;
/// `BLAST_WINDOW_SIZE_MEGABLAST` — megablast window (0 = disabled).
pub const BLAST_WINDOW_SIZE_MEGABLAST: i32 = 0;
/// `BLAST_WINDOW_SIZE_DISC` — discontiguous-megablast window (40).
pub const BLAST_WINDOW_SIZE_DISC: i32 = 40;

/// NCBI: BLAST_GapDecayDivisor (blast_stat.c:4079).
/// Computes the divisor used by sum-statistics to compensate for the
/// effect of choosing the best among multiple alignments:
/// `(1 - decayrate) * decayrate^(nsegs - 1)`. Typical `decayrate` values
/// are [`BLAST_GAP_DECAY_RATE_GAPPED`] (0.1) and [`BLAST_GAP_DECAY_RATE`] (0.5).
pub fn blast_gap_decay_divisor(decayrate: f64, nsegs: u32) -> f64 {
    // NCBI: `return (1. - decayrate) * BLAST_Powi(decayrate, nsegs - 1);`.
    (1.0 - decayrate) * crate::math::blast_powi(decayrate, (nsegs as i32) - 1)
}

/// NCBI: BLAST_Cutoffs (blast_stat.c:4089).
/// Given a desired
/// e-value `e_in`, a Karlin-Altschul block, a search space size, and an
/// optional decay-rate adjustment, returns a tuple
/// `(cutoff_score, effective_evalue)`.
///
/// * `dodecay=true, 0 < gap_decay_rate < 1` scales the input e-value by
///   `BLAST_GapDecayDivisor(gap_decay_rate, 1)` before converting to a
///   raw score (tightens the cutoff), then divides the recomputed
///   e-value back out on the return trip.
/// * The final cutoff is `max(s_floor, EtoS(e))`; callers pass a floor
///   (typically 1) to match NCBI's behavior where a user-specified
///   minimum score wins if it's larger than the statistics-derived one.
///
/// Returns `(1, e_in)` when the KarlinBlk is uncomputed (lambda/K/H == -1.0),
/// matching NCBI `BLAST_Cutoffs` (`blast_stat.c:4101-4102`) which uses the
/// exact -1.0 sentinel — not `< 0.0` — to detect "computation failed".
/// Note: this differs from `raw_to_evalue` and `evalue_to_raw`, which use
/// `< 0.0` to match their respective NCBI verbatim checks (NCBI itself is
/// inconsistent across these three functions; the Rust port preserves
/// each function's individual NCBI convention).
pub fn blast_cutoffs(
    s_floor: i32,
    e_in: f64,
    kbp: &KarlinBlk,
    searchsp: f64,
    dodecay: bool,
    gap_decay_rate: f64,
) -> (i32, f64) {
    if kbp.lambda == -1.0 || kbp.k == -1.0 || kbp.h == -1.0 {
        return (1, e_in);
    }
    let esave = e_in;
    let mut s = s_floor;
    let mut e = e_in;
    let mut s_changed = false;
    let mut es = 1i32;
    if e > 0.0 {
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            e *= blast_gap_decay_divisor(gap_decay_rate, 1);
        }
        es = kbp.evalue_to_raw(e, searchsp);
    }
    if es > s {
        s_changed = true;
        s = es;
    }
    // Recompute the e-value from the final cutoff when the input e was
    // non-positive or the cutoff didn't change. Mirrors NCBI's
    // `blast_stat.c:4134-4146`.
    let e_out = if esave <= 0.0 || !s_changed {
        let mut recomputed = blast_karlin_sto_e_simple(s, Some(kbp), searchsp as i64);
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            recomputed /= blast_gap_decay_divisor(gap_decay_rate, 1);
        }
        recomputed
    } else {
        esave
    };
    (s, e_out)
}

/// blast-rs: Pointer-style mutable-output adapter for cutoff computation; not
/// a direct NCBI C port.
pub fn blast_cutoffs_in_place(
    score: Option<&mut i32>,
    evalue: Option<&mut f64>,
    kbp: Option<&KarlinBlk>,
    searchsp: i64,
    dodecay: bool,
    gap_decay_rate: f64,
) -> i16 {
    let (Some(score), Some(evalue), Some(kbp)) = (score, evalue, kbp) else {
        return 1;
    };
    if kbp.lambda == -1.0 || kbp.k == -1.0 || kbp.h == -1.0 {
        return 1;
    }

    let mut s = *score;
    let mut e = *evalue;
    let esave = e;
    let mut s_changed = false;
    let mut es = 1;

    if e > 0.0 {
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            e *= blast_gap_decay_divisor(gap_decay_rate, 1);
        }
        es = blast_karlin_eto_s_simple(e, Some(kbp), searchsp);
    }

    if es > s {
        s_changed = true;
        s = es;
        *score = s;
    }

    if esave <= 0.0 || !s_changed {
        e = blast_karlin_sto_e_simple(s, Some(kbp), searchsp);
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            e /= blast_gap_decay_divisor(gap_decay_rate, 1);
        }
        *evalue = e;
    }

    0
}

// NCBI s_BlastSumP interpolation tables (blast_stat.c:4359-4379).
// Retained verbatim so any drift is easy to spot against the C source.

#[rustfmt::skip]
const SUM_P_TAB2: [f64; 19] = [
    0.01669,  0.0249,   0.03683,  0.05390,  0.07794,  0.1111,   0.1559,   0.2146,
    0.2890,   0.3794,   0.4836,   0.5965,   0.7092,   0.8114,   0.8931,   0.9490,
    0.9806,   0.9944,   0.9989,
];

#[rustfmt::skip]
const SUM_P_TAB3: [f64; 38] = [
    0.9806,   0.9944,   0.9989,   0.0001682,0.0002542,0.0003829,0.0005745,0.0008587,
    0.001278, 0.001893, 0.002789, 0.004088, 0.005958, 0.008627, 0.01240,  0.01770,
    0.02505,  0.03514,  0.04880,  0.06704,  0.09103,  0.1220,   0.1612,   0.2097,
    0.2682,   0.3368,   0.4145,   0.4994,   0.5881,   0.6765,   0.7596,   0.8326,
    0.8922,   0.9367,   0.9667,   0.9846,   0.9939,   0.9980,
];

#[rustfmt::skip]
const SUM_P_TAB4: [f64; 55] = [
    2.658e-07,4.064e-07,6.203e-07,9.450e-07,1.437e-06,2.181e-06,3.302e-06,4.990e-06,
    7.524e-06,1.132e-05,1.698e-05,2.541e-05,3.791e-05,5.641e-05,8.368e-05,0.0001237,
    0.0001823,0.0002677,0.0003915,0.0005704,0.0008275,0.001195, 0.001718, 0.002457,
    0.003494, 0.004942, 0.006948, 0.009702, 0.01346,  0.01853,  0.02532,  0.03431,
    0.04607,  0.06128,  0.08068,  0.1051,   0.1352,   0.1719,   0.2157,   0.2669,
    0.3254,   0.3906,   0.4612,   0.5355,   0.6110,   0.6849,   0.7544,   0.8168,
    0.8699,   0.9127,   0.9451,   0.9679,   0.9827,   0.9915,   0.9963,
];

/// NCBI: s_BlastSumP (blast_stat.c:4357).
/// Estimates the Sum P-value by calculation or interpolation. Accuracy:
/// ~2.5 digits throughout the range of `r` (number of segments) and `s`
/// (total score in nats, adjusted by `-r*log(K*N)`).
///
/// For `r = 0` returns `0.0`; for `r = 1` uses the closed-form
/// `1 - exp(-exp(-s))`; for `r` in `2..=4` uses the table-interpolation
/// branches at `blast_stat.c:4394-4404` (either the analytic tail or
/// the `kTable` interpolation); for `r >= 5` delegates to `s_blast_sum_p_calc`
/// (NCBI `s_BlastSumPCalc`, Romberg integration). Always returns `Some`.
pub fn s_blast_sum_p(r: u32, s: f64) -> Option<f64> {
    if r == 1 {
        // NCBI `blast_stat.c:4384`: `return -BLAST_Expm1(-exp(-s));`.
        return Some(-crate::math::expm1(-(-s).exp()));
    }
    if r == 0 {
        return Some(0.0);
    }
    if (2..=4).contains(&r) {
        let r1 = (r - 1) as i32;
        let r_i = r as i32;
        if s >= (r * r + r - 1) as f64 {
            // NCBI `blast_stat.c:4394`: `a = BLAST_LnGammaInt(r+1)`.
            let a = crate::math::blast_ln_gamma_int(r_i + 1);
            return Some(r_i as f64 * (r1 as f64 * s.ln() - s - a - a).exp());
        }
        if s > -2.0 * r as f64 {
            // Table interpolation — verbatim NCBI `blast_stat.c:4397-4403`:
            //   i = (Int4)(a = s+s+(4*r));
            //   a -= i;
            //   i = kTabsize[r2] - i;
            //   return a*kTable[r2][i-1] + (1.-a)*kTable[r2][i];
            let table: &[f64] = match r {
                2 => &SUM_P_TAB2,
                3 => &SUM_P_TAB3,
                4 => &SUM_P_TAB4,
                _ => unreachable!(),
            };
            let mut a = s + s + (4 * r) as f64;
            let mut i = a as i32;
            a -= i as f64;
            let tab_last = (table.len() - 1) as i32; // == NCBI `kTabsize[r2]`.
            i = tab_last - i;
            return Some(a * table[(i - 1) as usize] + (1.0 - a) * table[i as usize]);
        }
        return Some(1.0);
    }
    // r >= 5: delegate to `s_BlastSumPCalc` (Romberg integration).
    Some(s_blast_sum_p_calc(r, s))
}

/// NCBI: SRombergCbackArgs (blast_stat.c:4205).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SRombergCbackArgs {
    pub num_hsps: i32,
    pub num_hsps_minus_2: i32,
    pub adj1: f64,
    pub adj2: f64,
    pub sdvir: f64,
    pub epsilon: f64,
}

/// NCBI: s_OuterIntegralCback (blast_stat.c:4219).
pub fn s_outer_integral_cback(x: f64, callback_args: &SRombergCbackArgs) -> f64 {
    let y = (x - callback_args.sdvir).exp();
    if y == f64::INFINITY {
        return 0.0;
    }
    if callback_args.num_hsps_minus_2 == 0 {
        return (callback_args.adj2 - y).exp();
    }
    if x == 0.0 {
        return 0.0;
    }
    (callback_args.num_hsps_minus_2 as f64 * x.ln() + callback_args.adj2 - y).exp()
}

/// NCBI: s_InnerIntegralCback (blast_stat.c:4240).
pub fn s_inner_integral_cback(s: f64, callback_args: &SRombergCbackArgs) -> f64 {
    let mut outer_args = *callback_args;
    outer_args.adj2 = callback_args.adj1 - s;
    outer_args.sdvir = s / callback_args.num_hsps as f64;
    let mx = if s > 0.0 { outer_args.sdvir + 3.0 } else { 3.0 };
    crate::math::blast_romberg_integrate(
        |x| s_outer_integral_cback(x, &outer_args),
        0.0,
        mx,
        callback_args.epsilon,
        0,
        1,
    )
}

/// NCBI: s_BlastSumPCalc (blast_stat.c:4269).
/// Computes the Sum P-value via double Romberg integration for
/// `r >= 5` (callers with smaller `r` should go through `s_blast_sum_p` instead).
/// Matches the Karlin-Altschul PNAS 1993 formula with the iteratively
/// tightened `itmin` that NCBI uses when the convergence is marginal.
pub fn s_blast_sum_p_calc(r: u32, s: f64) -> f64 {
    let r_i = r as i32;
    if r == 1 {
        if s > 8.0 {
            return (-s).exp();
        }
        // NCBI `blast_stat.c:4271`: `return -BLAST_Expm1(-exp(-s));`.
        return -crate::math::expm1(-(-s).exp());
    }
    if r < 1 {
        return 0.0;
    }

    // Early-out bounds where the integral is essentially 1 ("no
    // significant distinction"). NCBI `blast_stat.c:4286-4305` uses
    // nested `if/else if { if (s <= …) return 1.0; }` — keep the nested
    // layout so the port is line-diffable against the C source.
    let rs = r as f64;
    #[allow(clippy::collapsible_if, clippy::collapsible_else_if)]
    if r < 8 {
        if s <= -2.3 * rs {
            return 1.0;
        }
    } else if r < 15 {
        if s <= -2.5 * rs {
            return 1.0;
        }
    } else if r < 27 {
        if s <= -3.0 * rs {
            return 1.0;
        }
    } else if r < 51 {
        if s <= -3.4 * rs {
            return 1.0;
        }
    } else if r < 101 {
        if s <= -4.0 * rs {
            return 1.0;
        }
    }

    let stddev = rs.sqrt();
    let stddev4 = 4.0 * stddev;
    let r1 = r - 1;

    if r > 100 {
        let est_mean = -(r as f64) * r1 as f64;
        if s <= est_mean - stddev4 {
            return 1.0;
        }
    }

    let logr = rs.ln();
    let mean = rs * (1.0 - logr) - 0.5;
    if s <= mean - stddev4 {
        return 1.0;
    }

    let (t0, itmin0) = if s >= mean {
        (s + 6.0 * stddev, 1)
    } else {
        (mean + 6.0 * stddev, 2)
    };

    // NCBI `blast_stat.c:4338`: `adj1 = num_hsps_minus_2*logr
    //                                   - BLAST_LnGammaInt(r1)
    //                                   - BLAST_LnGammaInt(r)`.
    /// NCBI `kSumpEpsilon` (`blast_stat.c:4276`): Romberg convergence
    /// epsilon for the sum-P double integral.
    const EPSILON: f64 = 0.002;
    let callback_args = SRombergCbackArgs {
        num_hsps: r_i,
        num_hsps_minus_2: r_i - 2,
        adj1: (r_i - 2) as f64 * logr
            - crate::math::blast_ln_gamma_int(r1 as i32)
            - crate::math::blast_ln_gamma_int(r_i),
        adj2: 0.0,
        sdvir: 0.0,
        epsilon: EPSILON,
    };

    // NCBI iteratively tightens `itmin` if the outer integral returns a
    // marginal value (blast_stat.c:4345).
    let mut itmin = itmin0;
    let mut d;
    loop {
        d = crate::math::blast_romberg_integrate(
            |sv| s_inner_integral_cback(sv, &callback_args),
            s,
            t0,
            EPSILON,
            0,
            itmin,
        );
        if d == f64::INFINITY {
            return d;
        }
        if !(s < mean && d < 0.4 && itmin < 4) {
            break;
        }
        itmin += 1;
    }
    d.min(1.0)
}

/// NCBI: BLAST_KarlinPtoE (blast_stat.c:4175).
/// naming: Rust separates the `PtoE` acronym group as `p_to_e`.
pub fn blast_karlin_p_to_e(p: f64) -> f64 {
    if p < 0.0 || p > 1.0 {
        return i32::MIN as f64;
    }
    if p == 1.0 {
        return i32::MAX as f64;
    }
    -crate::math::log1p(-p)
}

/// NCBI: BLAST_KarlinEtoP (blast_stat.c:4189).
/// naming: Rust separates the `EtoP` acronym group as `e_to_p`.
pub fn blast_karlin_e_to_p(evalue: f64) -> f64 {
    -crate::math::expm1(-evalue)
}

const SUM_E_CAP: f64 = i32::MAX as f64;

/// blast-rs: Ergonomic small-gap sum-E helper; not a direct NCBI C port.
/// Computes the e-value of a collection of distinct alignments
/// separated by small gaps. Matches NCBI's formula and cap at
/// `INT4_MAX`. Delegates the P-value step to `s_blast_sum_p`, which handles
/// every `num` via closed-form / table interpolation / Romberg
/// (`s_blast_sum_p_calc`). Returns `Some` on every finite input; the `Option`
/// return is retained for API symmetry with neighbouring helpers.
#[allow(clippy::too_many_arguments)]
pub fn small_gap_sum_e(
    starting_points: i32,
    num: u32,
    mut xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: f64,
    weight_divisor: f64,
) -> Option<f64> {
    let mut sum_e = if num == 1 {
        searchsp_eff * (-xsum).exp()
    } else {
        let pair_search_space = subject_length as f64 * query_length as f64;
        xsum -= pair_search_space.ln() + 2.0 * (num - 1) as f64 * (starting_points as f64).ln();
        // NCBI `blast_stat.c:4452` calls `BLAST_LnFactorial((double)num)`,
        // which goes straight to `s_LnGamma(num+1.0)` (no precomputed-table
        // fast path). Use our matching helper instead of the table-backed
        // `ln_factorial`, since the two paths differ by ~1 ULP for small n.
        xsum -= crate::math::blast_ln_factorial(num as f64);
        let p = s_blast_sum_p(num, xsum)?;
        blast_karlin_p_to_e(p) * (searchsp_eff / pair_search_space)
    };
    if weight_divisor == 0.0 {
        sum_e = SUM_E_CAP;
    } else {
        sum_e /= weight_divisor;
        if sum_e > SUM_E_CAP {
            sum_e = SUM_E_CAP;
        }
    }
    Some(sum_e)
}

/// NCBI: BLAST_SmallGapSumE (blast_stat.c:4418).
/// Port with the C
/// return type and argument shape.
#[allow(clippy::too_many_arguments)]
pub fn blast_small_gap_sum_e(
    starting_points: i32,
    num: i16,
    xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: i64,
    weight_divisor: f64,
) -> f64 {
    small_gap_sum_e(
        starting_points,
        num.max(0) as u32,
        xsum,
        query_length,
        subject_length,
        searchsp_eff as f64,
        weight_divisor,
    )
    .unwrap_or(SUM_E_CAP)
}

/// blast-rs: Ergonomic uneven-gap sum-E helper; not a direct NCBI C port.
/// Used for HSP collections with asymmetric gap widths — e.g. exons
/// separated by introns in translated searches.
#[allow(clippy::too_many_arguments)]
pub fn uneven_gap_sum_e(
    query_start_points: i32,
    subject_start_points: i32,
    num: u32,
    mut xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: f64,
    weight_divisor: f64,
) -> Option<f64> {
    let mut sum_e = if num == 1 {
        searchsp_eff * (-xsum).exp()
    } else {
        let pair_search_space = subject_length as f64 * query_length as f64;
        xsum -= pair_search_space.ln()
            + (num - 1) as f64
                * ((query_start_points as f64).ln() + (subject_start_points as f64).ln());
        // NCBI `blast_stat.c:4511` calls `BLAST_LnFactorial`; see comment
        // in `small_gap_sum_e` — use the helper that matches NCBI exactly.
        xsum -= crate::math::blast_ln_factorial(num as f64);
        let p = s_blast_sum_p(num, xsum)?;
        blast_karlin_p_to_e(p) * (searchsp_eff / pair_search_space)
    };
    if weight_divisor == 0.0 {
        sum_e = SUM_E_CAP;
    } else {
        sum_e /= weight_divisor;
        if sum_e > SUM_E_CAP {
            sum_e = SUM_E_CAP;
        }
    }
    Some(sum_e)
}

/// NCBI: BLAST_UnevenGapSumE (blast_stat.c:4491).
/// Port with the C
/// return type and argument shape.
#[allow(clippy::too_many_arguments)]
pub fn blast_uneven_gap_sum_e(
    query_start_points: i32,
    subject_start_points: i32,
    num: i16,
    xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: i64,
    weight_divisor: f64,
) -> f64 {
    uneven_gap_sum_e(
        query_start_points,
        subject_start_points,
        num.max(0) as u32,
        xsum,
        query_length,
        subject_length,
        searchsp_eff as f64,
        weight_divisor,
    )
    .unwrap_or(SUM_E_CAP)
}

/// blast-rs: Ergonomic large-gap sum-E helper; not a direct NCBI C port.
/// Computes the e-value of a collection of distinct alignments
/// separated by arbitrarily large gaps.
pub fn large_gap_sum_e(
    num: u32,
    mut xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: f64,
    weight_divisor: f64,
) -> Option<f64> {
    let mut sum_e = if num == 1 {
        searchsp_eff * (-xsum).exp()
    } else {
        let q = query_length as f64;
        let s = subject_length as f64;
        // NCBI `blast_stat.c:4555` calls `BLAST_LnFactorial`; see
        // `small_gap_sum_e` for the rationale on bypassing the table.
        xsum -= num as f64 * (s * q).ln() - crate::math::blast_ln_factorial(num as f64);
        let p = s_blast_sum_p(num, xsum)?;
        blast_karlin_p_to_e(p) * (searchsp_eff / (q * s))
    };
    if weight_divisor == 0.0 {
        sum_e = SUM_E_CAP;
    } else {
        sum_e /= weight_divisor;
        if sum_e > SUM_E_CAP {
            sum_e = SUM_E_CAP;
        }
    }
    Some(sum_e)
}

/// NCBI: BLAST_LargeGapSumE (blast_stat.c:4532).
/// Port with the C
/// return type and argument shape.
pub fn blast_large_gap_sum_e(
    num: i16,
    xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: i64,
    weight_divisor: f64,
) -> f64 {
    large_gap_sum_e(
        num.max(0) as u32,
        xsum,
        query_length,
        subject_length,
        searchsp_eff as f64,
        weight_divisor,
    )
    .unwrap_or(SUM_E_CAP)
}

/// Score frequency distribution.
#[derive(Debug, Clone)]
pub struct ScoreFreq {
    pub score_min: i32,
    pub score_max: i32,
    pub obs_min: i32,
    pub obs_max: i32,
    pub score_avg: f64,
    pub sprob: Vec<f64>, // probability for each score value
}

/// NCBI: BlastScoreChk (blast_stat.c:2095).
pub fn blast_score_chk(lo: i32, hi: i32) -> i16 {
    if lo >= 0
        || hi <= 0
        || !(BLAST_SCORE_MIN..=BLAST_SCORE_MAX).contains(&lo)
        || !(BLAST_SCORE_MIN..=BLAST_SCORE_MAX).contains(&hi)
    {
        return 1;
    }
    if hi - lo > BLAST_SCORE_MAX - BLAST_SCORE_MIN {
        return 1;
    }
    0
}

/// NCBI: Blast_ScoreFreqNew (blast_stat.c:2113).
pub fn blast_score_freq_new(score_min: i32, score_max: i32) -> Option<ScoreFreq> {
    if blast_score_chk(score_min, score_max) != 0 {
        return None;
    }
    let range = (score_max - score_min + 1) as usize;
    Some(ScoreFreq {
        score_min,
        score_max,
        obs_min: 0,
        obs_max: 0,
        score_avg: 0.0,
        sprob: vec![0.0; range],
    })
}

/// NCBI: Blast_ScoreFreqFree (blast_stat.c:941).
pub fn blast_score_freq_free(sfp: &mut Option<ScoreFreq>) -> Option<ScoreFreq> {
    if let Some(sfp) = sfp.as_mut() {
        sfp.sprob.clear();
    }
    *sfp = None;
    None
}

impl ScoreFreq {
    #[inline]
    /// blast-rs: Score-index accessor for owned probability storage; not a
    /// direct NCBI C port.
    fn p(&self, score: i32) -> f64 {
        if score < self.score_min || score > self.score_max {
            0.0
        } else {
            self.sprob[(score - self.score_min) as usize]
        }
    }

    #[inline]
    /// blast-rs: Mutable score-index accessor for owned probability storage;
    /// not a direct NCBI C port.
    fn p_mut(&mut self, score: i32) -> Option<&mut f64> {
        if score < self.score_min || score > self.score_max {
            None
        } else {
            Some(&mut self.sprob[(score - self.score_min) as usize])
        }
    }
}

/// NCBI: BlastScoreFreqCalc (blast_stat.c:2151).
pub fn blast_score_freq_calc(
    sbp: Option<&BlastScoreBlk>,
    sfp: Option<&mut ScoreFreq>,
    rfp1: Option<&BlastResFreq>,
    rfp2: Option<&BlastResFreq>,
) -> i16 {
    let (Some(sbp), Some(sfp), Some(rfp1), Some(rfp2)) = (sbp, sfp, rfp1, rfp2) else {
        return 1;
    };
    if sbp.loscore < sfp.score_min || sbp.hiscore > sfp.score_max {
        return 1;
    }
    if rfp1.prob.len() < sbp.alphabet_size
        || rfp2.prob.len() < sbp.alphabet_size
        || sbp.matrix.data.len() < sbp.alphabet_size
    {
        return 1;
    }

    for prob in &mut sfp.sprob {
        *prob = 0.0;
    }
    for index1 in 0..sbp.alphabet_size {
        if sbp.matrix.data[index1].len() < sbp.alphabet_size {
            return 1;
        }
        for index2 in 0..sbp.alphabet_size {
            let score = sbp.matrix.data[index1][index2];
            if score >= sbp.loscore {
                if let Some(slot) = sfp.p_mut(score) {
                    *slot += rfp1.prob[index1] * rfp2.prob[index2];
                }
            }
        }
    }

    let mut score_sum = 0.0;
    let mut obs_min = BLAST_SCORE_MIN;
    let mut obs_max = BLAST_SCORE_MIN;
    for score in sfp.score_min..=sfp.score_max {
        let prob = sfp.p(score);
        if prob > 0.0 {
            score_sum += prob;
            obs_max = score;
            if obs_min == BLAST_SCORE_MIN {
                obs_min = score;
            }
        }
    }
    sfp.obs_min = obs_min;
    sfp.obs_max = obs_max;

    let mut score_avg = 0.0;
    if score_sum > 0.0001 || score_sum < -0.0001 {
        for score in obs_min..=obs_max {
            if let Some(slot) = sfp.p_mut(score) {
                *slot /= score_sum;
                score_avg += score as f64 * *slot;
            }
        }
    }
    sfp.score_avg = score_avg;
    0
}

/// Port of NCBI `Blast_KarlinBlkNew` (`blast_stat.c:2860`).
pub fn blast_karlin_blk_new() -> KarlinBlk {
    KarlinBlk::default()
}

/// blast-rs: Adapter from public `ScoreFreq` to the internal Karlin solver
/// distribution; not a direct NCBI C port.
fn sf_dist_from_score_freq(sfp: &ScoreFreq) -> SfDist {
    SfDist {
        score_min: sfp.score_min,
        score_max: sfp.score_max,
        obs_min: sfp.obs_min,
        obs_max: sfp.obs_max,
        score_avg: sfp.score_avg,
        probs: sfp.sprob.clone(),
    }
}

/// Port of NCBI `Blast_KarlinBlkUngappedCalc` (`blast_stat.c:2699`).
pub fn blast_karlin_blk_ungapped_calc(kbp: Option<&mut KarlinBlk>, sfp: Option<&ScoreFreq>) -> i16 {
    let (Some(kbp), Some(sfp)) = (kbp, sfp) else {
        return 1;
    };
    let dist = sf_dist_from_score_freq(sfp);
    kbp.lambda = compute_lambda(&dist);
    if kbp.lambda < 0.0 {
        kbp.lambda = -1.0;
        kbp.h = -1.0;
        kbp.k = -1.0;
        kbp.log_k = f64::INFINITY;
        return 1;
    }
    kbp.h = blast_karlin_lto_h(&dist, kbp.lambda);
    if kbp.h < 0.0 {
        kbp.lambda = -1.0;
        kbp.h = -1.0;
        kbp.k = -1.0;
        kbp.log_k = f64::INFINITY;
        return 1;
    }
    kbp.k = blast_karlin_lh_to_k(&dist, kbp.lambda, kbp.h);
    if kbp.k < 0.0 {
        kbp.lambda = -1.0;
        kbp.h = -1.0;
        kbp.k = -1.0;
        kbp.log_k = f64::INFINITY;
        return 1;
    }
    kbp.log_k = kbp.k.ln();
    0
}

/// Port of NCBI `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2737`).
pub fn blast_score_blk_kbp_ungapped_calc(
    program: crate::program::ProgramType,
    sbp: Option<&mut BlastScoreBlk>,
    query: Option<&[u8]>,
    query_info: Option<&mut crate::queryinfo::QueryInfo>,
    blast_message: Option<&mut Option<Box<crate::diagnostics::BlastMessage>>>,
) -> i16 {
    let (Some(sbp), Some(query), Some(query_info)) = (sbp, query, query_info) else {
        return 1;
    };

    let status = blast_score_blk_kbp_ideal_calc(Some(sbp));
    if status != 0 {
        return status;
    }
    let Some(mut stdrfp) = blast_res_freq_new(Some(sbp)) else {
        return 1;
    };
    if blast_res_freq_std_comp(Some(sbp), Some(&mut stdrfp)) != 0 {
        return 1;
    }
    let Some(mut rfp) = blast_res_freq_new(Some(sbp)) else {
        return 1;
    };

    let context_count = query_info.contexts.len();
    if sbp.sfp.len() < context_count {
        sbp.sfp.resize_with(context_count, || None);
    }
    if sbp.kbp_std.len() < context_count {
        sbp.kbp_std.resize_with(context_count, KarlinBlk::default);
    }
    if sbp.kbp_psi.len() < context_count {
        sbp.kbp_psi.resize_with(context_count, KarlinBlk::default);
    }
    if sbp.kbp.len() < context_count {
        sbp.kbp.resize_with(context_count, KarlinBlk::default);
    }

    let check_ideal = program == crate::program::BLASTX
        || program == crate::program::TBLASTX
        || program == crate::program::RPS_TBLASTN;
    let mut valid_context = false;
    let mut blast_message = blast_message;

    for context in 0..context_count {
        if !query_info.contexts[context].is_valid {
            continue;
        }
        let query_length = query_info.contexts[context].query_length;
        let context_offset = query_info.contexts[context].query_offset;
        if query_length <= 0 || context_offset < 0 {
            query_info.contexts[context].is_valid = false;
            continue;
        }
        let begin = context_offset as usize;
        let end = begin.saturating_add(query_length as usize);
        if end > query.len() {
            query_info.contexts[context].is_valid = false;
            continue;
        }

        let buffer = &query[begin..end];
        if blast_res_freq_string(Some(sbp), Some(&mut rfp), Some(buffer), query_length) != 0 {
            query_info.contexts[context].is_valid = false;
            continue;
        }
        let Some(mut sfp) = blast_score_freq_new(sbp.loscore, sbp.hiscore) else {
            query_info.contexts[context].is_valid = false;
            continue;
        };
        if blast_score_freq_calc(Some(sbp), Some(&mut sfp), Some(&rfp), Some(&stdrfp)) != 0 {
            query_info.contexts[context].is_valid = false;
            continue;
        }

        let mut kbp = blast_karlin_blk_new();
        let loop_status = blast_karlin_blk_ungapped_calc(Some(&mut kbp), Some(&sfp));
        if loop_status != 0 {
            query_info.contexts[context].is_valid = false;
            sbp.sfp[context] = None;
            sbp.kbp_std[context] = KarlinBlk::default();
            if !crate::program::blast_query_is_translated(program) {
                if let Some(messages) = blast_message.as_mut() {
                    let _ = crate::diagnostics::blast_message_write(
                        *messages,
                        crate::diagnostics::BlastSeverity::Warning,
                        context as i32,
                        crate::diagnostics::K_BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS,
                    );
                }
            }
            continue;
        }
        if check_ideal {
            if let Some(ideal) = sbp.kbp_ideal.as_ref() {
                if kbp.lambda >= ideal.lambda {
                    let _ = blast_karlin_blk_copy(Some(&mut kbp), Some(ideal));
                }
            }
        }

        let mut psi_kbp = blast_karlin_blk_new();
        let loop_status = blast_karlin_blk_ungapped_calc(Some(&mut psi_kbp), Some(&sfp));
        if loop_status != 0 {
            query_info.contexts[context].is_valid = false;
            sbp.sfp[context] = None;
            sbp.kbp_std[context] = KarlinBlk::default();
            sbp.kbp_psi[context] = KarlinBlk::default();
            continue;
        }

        sbp.sfp[context] = Some(sfp);
        sbp.kbp_std[context] = kbp;
        sbp.kbp_psi[context] = psi_kbp;
        valid_context = true;
    }

    if !valid_context {
        if crate::program::blast_query_is_translated(program) {
            if let Some(messages) = blast_message.as_mut() {
                let _ = crate::diagnostics::blast_message_write(
                    *messages,
                    crate::diagnostics::BlastSeverity::Warning,
                    crate::diagnostics::K_BLAST_MESSAGE_NO_CONTEXT,
                    crate::diagnostics::K_BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS,
                );
            }
        }
        return 1;
    }

    sbp.kbp = if crate::program::blast_query_is_pssm(program) {
        sbp.kbp_psi.clone()
    } else {
        sbp.kbp_std.clone()
    };
    0
}

/// Gapped Karlin-Altschul parameters (precomputed tables).
#[derive(Debug, Clone)]
pub struct GappedParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub alpha: f64,
    pub beta: f64,
}

/// Protein matrix statistics row values used in pairwise/XML reporting.
#[derive(Debug, Clone, Copy)]
pub struct ProteinMatrixStats {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub a: f64,
    pub alpha: f64,
    pub sigma: f64,
}

/// Gumbel block parameters for Spouge finite-size correction (FSC) e-value.
/// Port of NCBI Blast_GumbelBlk from blast_stat.h.
#[derive(Debug, Clone)]
pub struct GumbelBlk {
    pub lambda: f64,    // gbp->Lambda (gapped lambda from table)
    pub a: f64,         // gbp->a (alpha from array[6])
    pub b: f64,         // gbp->b = 2*G*(a_un - a)
    pub alpha: f64,     // gbp->Alpha (alpha_v from array[9])
    pub beta: f64,      // gbp->Beta = 2*G*(Alpha_un - Alpha)
    pub sigma: f64,     // gbp->Sigma (from array[10])
    pub tau: f64,       // gbp->Tau = 2*G*(Alpha_un - Sigma)
    pub db_length: i64, // total database length
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatrixStatRow {
    gap_open: i32,
    gap_extend: i32,
    lambda: f64,
    k: f64,
    h: f64,
    alpha: f64,
    beta: f64,
    _gumbel_a: f64,
    alpha_v: f64,
    sigma: f64,
}

impl MatrixStatRow {
    const fn new(
        gap_open: i32,
        gap_extend: i32,
        lambda: f64,
        k: f64,
        h: f64,
        alpha: f64,
        beta: f64,
        _gumbel_a: f64,
        alpha_v: f64,
        sigma: f64,
    ) -> Self {
        Self {
            gap_open,
            gap_extend,
            lambda,
            k,
            h,
            alpha,
            beta,
            _gumbel_a,
            alpha_v,
            sigma,
        }
    }

    /// blast-rs: Convert a matrix table row into public gapped parameters; not
    /// a direct NCBI C port.
    fn gapped_params(self) -> GappedParams {
        GappedParams {
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            lambda: self.lambda,
            k: self.k,
            h: self.h,
            alpha: self.alpha,
            beta: self.beta,
        }
    }

    /// blast-rs: Convert a matrix table row into display statistics; not a
    /// direct NCBI C port.
    fn display_params(self) -> ProteinMatrixStats {
        ProteinMatrixStats {
            lambda: self.lambda,
            k: self.k,
            h: self.h,
            a: self.alpha,
            alpha: self.alpha_v,
            sigma: self.sigma,
        }
    }
}

/// blast-rs: Spouge score-to-E-value helper with gap-decay adjustment; not a
/// direct NCBI C port.
/// Uses per-subject lengths for more accurate e-values than simple Karlin formula.
/// Spouge e-value with NCBI's HSP-linking gap-decay correction
/// applied. Used by translated searches (blastx/tblastn/tblastx) and
/// ungapped paths where `do_sum_stats=TRUE` (set in
/// `blast_options.c:1467`). NCBI applies this divisor in
/// `link_hsps.c:1791` via `link_hsp_params->gap_decay_rate`. The
/// default rate is `BLAST_GAP_DECAY_RATE_GAPPED = 0.1` for protein
/// gapped, giving an effective `/0.9 ≈ ×1.111` scaling vs. the raw
/// per-pair Spouge.
pub fn spouge_evalue_with_gap_decay(
    score: i32,
    kbp: &KarlinBlk,
    gbp: &GumbelBlk,
    query_length: i32,
    subject_length: i32,
) -> f64 {
    let raw = spouge_evalue(score, kbp, gbp, query_length, subject_length);
    raw / blast_gap_decay_divisor(BLAST_GAP_DECAY_RATE_GAPPED, 1)
}

/// blast-rs: Ergonomic Spouge score-to-E-value helper; not a direct NCBI C
/// port. Returns the raw
/// per-pair Spouge finite-size-corrected e-value. Some NCBI call sites
/// (e.g. translated searches with HSP linking via `link_hsps.c:1791`)
/// further divide by `BLAST_GapDecayDivisor(gap_decay_rate, 1)` —
/// callers in those paths should use [`spouge_evalue_with_gap_decay`]
/// instead. Standard blastp / direct traceback paths
/// (`blast_traceback.c:234`, `blast_kappa.c:419`) pass
/// `gap_decay_rate=0` and use the raw value unchanged.
pub fn spouge_evalue(
    score: i32,
    kbp: &KarlinBlk,
    gbp: &GumbelBlk,
    query_length: i32,
    subject_length: i32,
) -> f64 {
    let scale_factor = kbp.lambda / gbp.lambda;
    let db_scale_factor = if gbp.db_length > 0 {
        gbp.db_length as f64 / subject_length as f64
    } else {
        1.0
    };

    let ai_hat = gbp.a * scale_factor;
    let bi_hat = gbp.b;
    let alphai_hat = gbp.alpha * scale_factor;
    let betai_hat = gbp.beta;
    let sigma_hat = gbp.sigma * scale_factor;
    let tau_hat = gbp.tau;

    let aj_hat = ai_hat;
    let bj_hat = bi_hat;
    let alphaj_hat = alphai_hat;
    let betaj_hat = betai_hat;

    const CONST_VAL: f64 = 0.39894228040143267793994605993438; // 1/sqrt(2*PI)

    let y = score as f64;
    let m = query_length as f64;
    let n = subject_length as f64;

    let m_li_y = m - (ai_hat * y + bi_hat);
    let vi_y = (2.0 * alphai_hat / kbp.lambda).max(alphai_hat * y + betai_hat);
    let sqrt_vi_y = vi_y.sqrt();
    let m_f = m_li_y / sqrt_vi_y;
    let p_m_f = crate::math::erfc(-m_f / std::f64::consts::SQRT_2) / 2.0;
    let p1 = m_li_y * p_m_f + sqrt_vi_y * CONST_VAL * (-0.5 * m_f * m_f).exp();

    let n_lj_y = n - (aj_hat * y + bj_hat);
    let vj_y = (2.0 * alphaj_hat / kbp.lambda).max(alphaj_hat * y + betaj_hat);
    let sqrt_vj_y = vj_y.sqrt();
    let n_f = n_lj_y / sqrt_vj_y;
    let p_n_f = crate::math::erfc(-n_f / std::f64::consts::SQRT_2) / 2.0;
    let p2 = n_lj_y * p_n_f + sqrt_vj_y * CONST_VAL * (-0.5 * n_f * n_f).exp();

    let c_y = (2.0 * sigma_hat / kbp.lambda).max(sigma_hat * y + tau_hat);
    let area = p1 * p2 + c_y * p_m_f * p_n_f;

    // NCBI `blast_stat.c:5232` returns the raw computed value:
    //   `e_value = area * k_ * exp(-lambda_ * y_) * db_scale_factor;
    //    ASSERT(e_value >= 0.0); return e_value;`
    // The ASSERT is a no-op in NCBI's release builds, so we mirror that
    // by returning the raw value rather than clamping at zero. In
    // practice `area` is always non-negative for valid alignments, but
    // ULP-level FP errors in `erfc` could conceivably produce tiny
    // negative results — preserving NCBI's behavior keeps byte-identity.
    area * kbp.k * (-kbp.lambda * y).exp() * db_scale_factor
}

/// NCBI: BLAST_SpougeStoE (blast_stat.c:5176).
/// Port with pointer-shaped arguments.
pub fn blast_spouge_sto_e(
    y_: i32,
    kbp: Option<&KarlinBlk>,
    gbp: Option<&GumbelBlk>,
    m_: i32,
    n_: i32,
) -> f64 {
    let (Some(kbp), Some(gbp)) = (kbp, gbp) else {
        return -1.0;
    };
    spouge_evalue(y_, kbp, gbp, m_, n_)
}

/// blast-rs: Ergonomic Spouge E-value-to-score helper; not a direct NCBI C
/// port. Binary-search the
/// raw score `S` such that `BLAST_SpougeStoE(S, kbp, gbp, m, n) <= e0`.
///
/// Used for the per-context cutoff in `BlastHitSavingParametersUpdate`
/// (`blast_parameters.c:940`): the cutoff is computed at `m=query_length` and
/// `n=avg_subject_length`, with `gbp->db_length` set to the full DB. The
/// binary search uses Spouge's full FSC formula, so the cutoff is more
/// conservative than the simple Karlin `evalue_to_raw` formula.
pub fn spouge_etos(e0: f64, kbp: &KarlinBlk, gbp: &GumbelBlk, m: i32, n: i32) -> i32 {
    let db_scale_factor = if gbp.db_length > 0 {
        gbp.db_length as f64
    } else {
        1.0
    };

    // C: `b = MAX((int)(log(db_scale_factor/e0) / kbp->Lambda), 2);`
    let mut b = ((db_scale_factor / e0).ln() / kbp.lambda) as i32;
    if b < 2 {
        b = 2;
    }

    let mut a: i32 = 0;
    let mut e = spouge_evalue(b, kbp, gbp, m, n);

    if e > e0 {
        while e > e0 {
            a = b;
            // C: `b *= 2;` — note the integer overflow possibility on a
            // pathological input is mirrored from C verbatim
            b = b.saturating_mul(2);
            e = spouge_evalue(b, kbp, gbp, m, n);
            if b == i32::MAX {
                break;
            }
        }
    }
    while b - a > 1 {
        let c = (a + b) / 2;
        e = spouge_evalue(c, kbp, gbp, m, n);
        if e > e0 {
            a = c;
        } else {
            b = c;
        }
    }
    let _ = e; // suppress unused on the final iteration
    a
}

/// NCBI: BLAST_SpougeEtoS (blast_stat.c:5236).
/// Port with pointer-shaped arguments.
pub fn blast_spouge_eto_s(
    e0: f64,
    kbp: Option<&KarlinBlk>,
    gbp: Option<&GumbelBlk>,
    m: i32,
    n: i32,
) -> i32 {
    let (Some(kbp), Some(gbp)) = (kbp, gbp) else {
        return BLAST_SCORE_MIN;
    };
    spouge_etos(e0, kbp, gbp, m, n)
}

/// blast-rs: BLOSUM62-specific Gumbel table helper; not a direct NCBI C port.
pub fn protein_gumbel_blk(gap_open: i32, gap_extend: i32, db_length: i64) -> Option<GumbelBlk> {
    matrix_gumbel_blk("BLOSUM62", gap_open, gap_extend, db_length)
}

/// Build Gumbel block for the named protein scoring matrix.
/// blast-rs: Matrix-specific Gumbel table helper; not a direct NCBI C port.
pub fn matrix_gumbel_blk(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    db_length: i64,
) -> Option<GumbelBlk> {
    let rows = matrix_stat_rows(matrix_name)?;
    let ungapped = rows.first()?;
    let row = matrix_stat_row_for_gap_costs(matrix_name, gap_open, gap_extend)?;
    let g = (gap_open + gap_extend) as f64;

    Some(GumbelBlk {
        lambda: row.lambda,
        a: row.alpha,
        b: 2.0 * g * (ungapped.alpha - row.alpha),
        alpha: row.alpha_v,
        beta: 2.0 * g * (ungapped.alpha_v - row.alpha_v),
        sigma: row.sigma,
        tau: 2.0 * g * (ungapped.alpha_v - row.sigma),
        db_length,
    })
}

/// Build Gumbel block for NCBI's IDENTITY matrix.
/// blast-rs: IDENTITY-matrix Gumbel table helper; not a direct NCBI C port.
pub fn identity_gumbel_blk(gap_open: i32, gap_extend: i32, db_length: i64) -> Option<GumbelBlk> {
    matrix_gumbel_blk("IDENTITY", gap_open, gap_extend, db_length)
}

/// NCBI: Blast_KarlinBlkGappedLoadFromTables (blast_stat.c:3577).
pub fn blast_karlin_blk_gapped_load_from_tables(
    kbp: Option<&mut KarlinBlk>,
    gap_open: i32,
    gap_extend: i32,
    matrix_name: Option<&str>,
    standard_only: bool,
) -> i16 {
    let Some(matrix_name) = matrix_name else {
        return -1;
    };
    let Some(rows) = matrix_stat_rows_with_standard_only(matrix_name, standard_only) else {
        return 1;
    };
    let Some(row) = rows
        .iter()
        .copied()
        .find(|row| row.gap_open == gap_open && row.gap_extend == gap_extend)
    else {
        return 2;
    };
    if let Some(kbp) = kbp {
        kbp.lambda = row.lambda;
        kbp.k = row.k;
        kbp.log_k = row.k.ln();
        kbp.h = row.h;
        kbp.round_down = false;
    }
    0
}

/// NCBI: Blast_KarlinBlkGappedCalc (blast_stat.c:3527).
pub fn blast_karlin_blk_gapped_calc(
    kbp: Option<&mut KarlinBlk>,
    gap_open: i32,
    gap_extend: i32,
    matrix_name: Option<&str>,
    error_return: Option<&mut Option<Box<crate::diagnostics::BlastMessage>>>,
) -> i16 {
    let status =
        blast_karlin_blk_gapped_load_from_tables(kbp, gap_open, gap_extend, matrix_name, false);
    if status != 0 {
        report_matrix_table_error(error_return, status, matrix_name, gap_open, gap_extend);
    }
    status
}

/// NCBI: Blast_GumbelBlkLoadFromTables (blast_stat.c:3696).
pub fn blast_gumbel_blk_load_from_tables(
    gbp: Option<&mut GumbelBlk>,
    gap_open: i32,
    gap_extend: i32,
    matrix_name: Option<&str>,
) -> i16 {
    let Some(matrix_name) = matrix_name else {
        return -1;
    };
    if matrix_stat_rows(matrix_name).is_none() {
        return 1;
    }
    let Some(gumbel) = matrix_gumbel_blk(matrix_name, gap_open, gap_extend, 0) else {
        return 2;
    };
    if let Some(gbp) = gbp {
        *gbp = gumbel;
    }
    0
}

/// NCBI: Blast_GumbelBlkCalc (blast_stat.c:3652).
pub fn blast_gumbel_blk_calc(
    gbp: Option<&mut GumbelBlk>,
    gap_open: i32,
    gap_extend: i32,
    matrix_name: Option<&str>,
    error_return: Option<&mut Option<Box<crate::diagnostics::BlastMessage>>>,
) -> i16 {
    let status = blast_gumbel_blk_load_from_tables(gbp, gap_open, gap_extend, matrix_name);
    if status != 0 {
        report_matrix_table_error(error_return, status, matrix_name, gap_open, gap_extend);
    }
    status
}

/// blast-rs: Shared diagnostic formatter for matrix-table lookup errors; not
/// a direct NCBI C port.
fn report_matrix_table_error(
    error_return: Option<&mut Option<Box<crate::diagnostics::BlastMessage>>>,
    status: i16,
    matrix_name: Option<&str>,
    gap_open: i32,
    gap_extend: i32,
) {
    let Some(error_return) = error_return else {
        return;
    };
    let Some(matrix_name) = matrix_name else {
        return;
    };
    match status {
        1 => {
            crate::diagnostics::blast_message_write(
                error_return,
                crate::diagnostics::BlastSeverity::Error,
                crate::diagnostics::K_BLAST_MESSAGE_NO_CONTEXT,
                &format!("{matrix_name} is not a supported matrix"),
            );
            for info in blast_load_matrix_values(false) {
                crate::diagnostics::blast_message_write(
                    error_return,
                    crate::diagnostics::BlastSeverity::Error,
                    crate::diagnostics::K_BLAST_MESSAGE_NO_CONTEXT,
                    &format!("{} is a supported matrix", info.name),
                );
            }
        }
        2 => {
            crate::diagnostics::blast_message_write(
                error_return,
                crate::diagnostics::BlastSeverity::Error,
                crate::diagnostics::K_BLAST_MESSAGE_NO_CONTEXT,
                &format!(
                    "Gap existence and extension values of {gap_open} and {gap_extend} not supported for {matrix_name}"
                ),
            );
            for message in blast_karlin_report_allowed_values(matrix_name) {
                crate::diagnostics::blast_message_write(
                    error_return,
                    crate::diagnostics::BlastSeverity::Error,
                    crate::diagnostics::K_BLAST_MESSAGE_NO_CONTEXT,
                    &message,
                );
            }
        }
        _ => {}
    }
}

/// blast-rs: Uniform nucleotide ungapped-lambda helper; not a direct NCBI C port.
/// Solves: 0.25 * exp(lambda * reward) + 0.75 * exp(lambda * penalty) = 1
/// using Newton's method.
pub fn compute_ungapped_lambda(reward: i32, penalty: i32) -> f64 {
    let r = reward as f64;
    let p = penalty as f64;
    // Initial guess
    let mut lambda = 1.0;
    for _ in 0..100 {
        let er = (lambda * r).exp();
        let ep = (lambda * p).exp();
        let f = 0.25 * er + 0.75 * ep - 1.0;
        let fp = 0.25 * r * er + 0.75 * p * ep;
        if fp.abs() < 1e-30 {
            break;
        }
        let delta = f / fp;
        lambda -= delta;
        if delta.abs() < 1e-12 {
            break;
        }
    }
    lambda
}

/// blast-rs: Uniform nucleotide ungapped-K helper; not a direct NCBI C port.
/// K = (1 - exp(-lambda)) * (1 - exp(-lambda)) / (sum of prob * exp(lambda*score) * score^2 / 2)
/// Simplified for uniform base frequencies.
pub fn compute_ungapped_k(lambda: f64, reward: i32, penalty: i32) -> f64 {
    let r = reward as f64;
    let p = penalty as f64;
    let er = (lambda * r).exp();
    let ep = (lambda * p).exp();
    // H = lambda * (0.25 * r * er + 0.75 * p * ep)
    let h = lambda * (0.25 * r * er + 0.75 * p * ep);
    // K approximation for simple scoring
    // K = H / (lambda * lambda * variance)
    // For nucleotide with uniform freqs, use the standard approximation
    let variance = 0.25 * r * r * er + 0.75 * p * p * ep;
    if variance.abs() < 1e-30 {
        return 0.1;
    }
    h / (lambda * variance)
}

/// blast-rs: Uniform nucleotide ungapped Karlin-block helper; not a direct
/// NCBI C port.
pub fn compute_ungapped_kbp(reward: i32, penalty: i32) -> KarlinBlk {
    let lambda = compute_ungapped_lambda(reward, penalty);
    let k = compute_ungapped_k(lambda, reward, penalty);
    let h = lambda
        * (0.25 * (reward as f64) * (lambda * reward as f64).exp()
            + 0.75 * (penalty as f64) * (lambda * penalty as f64).exp());
    KarlinBlk {
        lambda,
        k,
        log_k: k.ln(),
        h,
        round_down: false,
    }
}

/// Compute length adjustment using the Altschul-Gish formula.
/// This adjusts query and database lengths to account for edge effects.
/// length_adj satisfies: K * (q - length_adj) * (d - N * length_adj) * exp(-lambda * S) = threshold
/// blast-rs: Defaulted length-adjustment helper; not a direct NCBI C port.
/// For callers that have no explicit alpha/beta
/// (typically ungapped or missing-gapped-params paths), NCBI uses
/// `alpha/lambda = 1/H` and `beta = 0`; that convention is applied
/// here so the verbatim bracketing algorithm is exercised.
pub fn compute_length_adjustment(
    query_length: i32,
    db_length: i64,
    num_seqs: i32,
    kbp: &KarlinBlk,
) -> i32 {
    if !kbp.is_valid() || query_length <= 0 || db_length <= 0 {
        return 0;
    }
    if kbp.h <= 0.0 {
        return 0;
    }
    let (adj, _converged) = compute_length_adjustment_exact(
        kbp.k,
        kbp.log_k,
        1.0 / kbp.h,
        0.0,
        query_length,
        db_length,
        num_seqs,
    );
    adj
}

const I2MAX: i32 = i16::MAX as i32;

/// NCBI `blast_stat.c` protein matrix value rows. The first row is the
/// ungapped row; remaining rows are supported gapped costs.
const BLOSUM45_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7, 0.641318, 9.611060, 9.611060,
    ),
    MatrixStatRow::new(
        13, 3, 0.207, 0.049, 0.14, 1.5, -22.0, 0.671128, 35.855900, 35.963900,
    ),
    MatrixStatRow::new(
        12, 3, 0.199, 0.039, 0.11, 1.8, -34.0, 0.691530, 45.693600, 45.851700,
    ),
    MatrixStatRow::new(
        11, 3, 0.190, 0.031, 0.095, 2.0, -38.0, 0.691181, 62.874100, 63.103700,
    ),
    MatrixStatRow::new(
        10, 3, 0.179, 0.023, 0.075, 2.4, -51.0, 0.710529, 88.286800, 88.639100,
    ),
    MatrixStatRow::new(
        16, 2, 0.210, 0.051, 0.14, 1.5, -24.0, 0.666680, 36.279800, 36.452400,
    ),
    MatrixStatRow::new(
        15, 2, 0.203, 0.041, 0.12, 1.7, -31.0, 0.673871, 44.825700, 45.060400,
    ),
    MatrixStatRow::new(
        14, 2, 0.195, 0.032, 0.10, 1.9, -36.0, 0.685753, 60.736200, 61.102300,
    ),
    MatrixStatRow::new(
        13, 2, 0.185, 0.024, 0.084, 2.2, -45.0, 0.698480, 85.148100, 85.689400,
    ),
    MatrixStatRow::new(
        12, 2, 0.171, 0.016, 0.061, 2.8, -65.0, 0.713429, 127.758000, 128.582000,
    ),
    MatrixStatRow::new(
        19, 1, 0.205, 0.040, 0.11, 1.9, -43.0, 0.672302, 53.071400, 53.828200,
    ),
    MatrixStatRow::new(
        18, 1, 0.198, 0.032, 0.10, 2.0, -43.0, 0.682580, 72.342400, 73.403900,
    ),
    MatrixStatRow::new(
        17, 1, 0.189, 0.024, 0.079, 2.4, -57.0, 0.695035, 103.055000, 104.721000,
    ),
    MatrixStatRow::new(
        16, 1, 0.176, 0.016, 0.063, 2.8, -67.0, 0.712966, 170.100000, 173.003000,
    ),
];

const BLOSUM50_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0, 0.609639, 5.388310, 5.388310,
    ),
    MatrixStatRow::new(
        13, 3, 0.212, 0.063, 0.19, 1.1, -16.0, 0.639287, 18.113800, 18.202800,
    ),
    MatrixStatRow::new(
        12, 3, 0.206, 0.055, 0.17, 1.2, -18.0, 0.644715, 22.654600, 22.777700,
    ),
    MatrixStatRow::new(
        11, 3, 0.197, 0.042, 0.14, 1.4, -25.0, 0.656327, 29.861100, 30.045700,
    ),
    MatrixStatRow::new(
        10, 3, 0.186, 0.031, 0.11, 1.7, -34.0, 0.671150, 42.393800, 42.674000,
    ),
    MatrixStatRow::new(
        9, 3, 0.172, 0.022, 0.082, 2.1, -48.0, 0.694326, 66.069600, 66.516400,
    ),
    MatrixStatRow::new(
        16, 2, 0.215, 0.066, 0.20, 1.05, -15.0, 0.633899, 17.951800, 18.092100,
    ),
    MatrixStatRow::new(
        15, 2, 0.210, 0.058, 0.17, 1.2, -20.0, 0.641985, 21.940100, 22.141800,
    ),
    MatrixStatRow::new(
        14, 2, 0.202, 0.045, 0.14, 1.4, -27.0, 0.650682, 28.681200, 28.961900,
    ),
    MatrixStatRow::new(
        13, 2, 0.193, 0.035, 0.12, 1.6, -32.0, 0.660984, 42.059500, 42.471600,
    ),
    MatrixStatRow::new(
        12, 2, 0.181, 0.025, 0.095, 1.9, -41.0, 0.678090, 63.747600, 64.397300,
    ),
    MatrixStatRow::new(
        19, 1, 0.212, 0.057, 0.18, 1.2, -21.0, 0.635714, 26.311200, 26.923300,
    ),
    MatrixStatRow::new(
        18, 1, 0.207, 0.050, 0.15, 1.4, -28.0, 0.643523, 34.903700, 35.734800,
    ),
    MatrixStatRow::new(
        17, 1, 0.198, 0.037, 0.12, 1.6, -33.0, 0.654504, 48.895800, 50.148600,
    ),
    MatrixStatRow::new(
        16, 1, 0.186, 0.025, 0.10, 1.9, -42.0, 0.667750, 76.469100, 78.443000,
    ),
    MatrixStatRow::new(
        15, 1, 0.171, 0.015, 0.063, 2.7, -76.0, 0.694575, 140.053000, 144.160000,
    ),
];

const BLOSUM62_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.964660, 4.964660,
    ),
    MatrixStatRow::new(
        11, 2, 0.297, 0.082, 0.27, 1.1, -10.0, 0.641766, 12.673800, 12.757600,
    ),
    MatrixStatRow::new(
        10, 2, 0.291, 0.075, 0.23, 1.3, -15.0, 0.649362, 16.474000, 16.602600,
    ),
    MatrixStatRow::new(
        9, 2, 0.279, 0.058, 0.19, 1.5, -19.0, 0.659245, 22.751900, 22.950000,
    ),
    MatrixStatRow::new(
        8, 2, 0.264, 0.045, 0.15, 1.8, -26.0, 0.672692, 35.483800, 35.821300,
    ),
    MatrixStatRow::new(
        7, 2, 0.239, 0.027, 0.10, 2.5, -46.0, 0.702056, 61.238300, 61.886000,
    ),
    MatrixStatRow::new(
        6, 2, 0.201, 0.012, 0.061, 3.3, -58.0, 0.740802, 140.417000, 141.882000,
    ),
    MatrixStatRow::new(
        13, 1, 0.292, 0.071, 0.23, 1.2, -11.0, 0.647715, 19.506300, 19.893100,
    ),
    MatrixStatRow::new(
        12, 1, 0.283, 0.059, 0.19, 1.5, -19.0, 0.656391, 27.856200, 28.469900,
    ),
    MatrixStatRow::new(
        11, 1, 0.267, 0.041, 0.14, 1.9, -30.0, 0.669720, 42.602800, 43.636200,
    ),
    MatrixStatRow::new(
        10, 1, 0.243, 0.024, 0.10, 2.5, -44.0, 0.693267, 83.178700, 85.065600,
    ),
    MatrixStatRow::new(
        9, 1, 0.206, 0.010, 0.052, 4.0, -87.0, 0.731887, 210.333000, 214.842000,
    ),
];

const BLOSUM80_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6, 0.564057, 1.918130, 1.918130,
    ),
    MatrixStatRow::new(
        25, 2, 0.342, 0.17, 0.66, 0.52, -1.6, 0.563956, 1.731000, 1.731300,
    ),
    MatrixStatRow::new(
        13, 2, 0.336, 0.15, 0.57, 0.59, -3.0, 0.570979, 2.673470, 2.692300,
    ),
    MatrixStatRow::new(
        9, 2, 0.319, 0.11, 0.42, 0.76, -6.0, 0.587837, 5.576090, 5.667860,
    ),
    MatrixStatRow::new(
        8, 2, 0.308, 0.090, 0.35, 0.89, -9.0, 0.597556, 7.536950, 7.686230,
    ),
    MatrixStatRow::new(
        7, 2, 0.293, 0.070, 0.27, 1.1, -14.0, 0.615254, 11.586600, 11.840400,
    ),
    MatrixStatRow::new(
        6, 2, 0.268, 0.045, 0.19, 1.4, -19.0, 0.644054, 19.958100, 20.441200,
    ),
    MatrixStatRow::new(
        11, 1, 0.314, 0.095, 0.35, 0.90, -9.0, 0.590702, 8.808610, 9.223320,
    ),
    MatrixStatRow::new(
        10, 1, 0.299, 0.071, 0.27, 1.1, -14.0, 0.609620, 13.833800, 14.533400,
    ),
    MatrixStatRow::new(
        9, 1, 0.279, 0.048, 0.20, 1.4, -19.0, 0.623800, 24.252000, 25.490400,
    ),
];

const BLOSUM90_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4, 0.544178, 1.377760, 1.377760,
    ),
    MatrixStatRow::new(
        9, 2, 0.310, 0.12, 0.46, 0.67, -6.0, 0.570267, 4.232290, 4.334170,
    ),
    MatrixStatRow::new(
        8, 2, 0.300, 0.099, 0.39, 0.76, -7.0, 0.581580, 5.797020, 5.961420,
    ),
    MatrixStatRow::new(
        7, 2, 0.283, 0.072, 0.30, 0.93, -11.0, 0.600024, 9.040880, 9.321600,
    ),
    MatrixStatRow::new(
        6, 2, 0.259, 0.048, 0.22, 1.2, -16.0, 0.629344, 16.024400, 16.531600,
    ),
    MatrixStatRow::new(
        11, 1, 0.302, 0.093, 0.39, 0.78, -8.0, 0.576919, 7.143250, 7.619190,
    ),
    MatrixStatRow::new(
        10, 1, 0.290, 0.075, 0.28, 1.04, -15.0, 0.591366, 11.483900, 12.269800,
    ),
    MatrixStatRow::new(
        9, 1, 0.265, 0.044, 0.20, 1.3, -19.0, 0.613013, 21.408300, 22.840900,
    ),
];

const PAM250_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.2252, 0.0868, 0.2223, 0.98, -5.0, 0.660059, 11.754300, 11.754300,
    ),
    MatrixStatRow::new(
        15, 3, 0.205, 0.049, 0.13, 1.6, -23.0, 0.687656, 34.578400, 34.928000,
    ),
    MatrixStatRow::new(
        14, 3, 0.200, 0.043, 0.12, 1.7, -26.0, 0.689768, 43.353000, 43.443800,
    ),
    MatrixStatRow::new(
        13, 3, 0.194, 0.036, 0.10, 1.9, -31.0, 0.697431, 50.948500, 51.081700,
    ),
    MatrixStatRow::new(
        12, 3, 0.186, 0.029, 0.085, 2.2, -41.0, 0.704565, 69.606500, 69.793600,
    ),
    MatrixStatRow::new(
        11, 3, 0.174, 0.020, 0.070, 2.5, -48.0, 0.722438, 98.653500, 98.927100,
    ),
    MatrixStatRow::new(
        17, 2, 0.204, 0.047, 0.12, 1.7, -28.0, 0.684799, 41.583800, 41.735800,
    ),
    MatrixStatRow::new(
        16, 2, 0.198, 0.038, 0.11, 1.8, -29.0, 0.691098, 51.635200, 51.843900,
    ),
    MatrixStatRow::new(
        15, 2, 0.191, 0.031, 0.087, 2.2, -44.0, 0.699051, 67.256700, 67.558500,
    ),
    MatrixStatRow::new(
        14, 2, 0.182, 0.024, 0.073, 2.5, -53.0, 0.714103, 96.315100, 96.756800,
    ),
    MatrixStatRow::new(
        13, 2, 0.171, 0.017, 0.059, 2.9, -64.0, 0.728738, 135.653000, 136.339000,
    ),
    MatrixStatRow::new(
        21, 1, 0.205, 0.045, 0.11, 1.8, -34.0, 0.683265, 48.728200, 49.218800,
    ),
    MatrixStatRow::new(
        20, 1, 0.199, 0.037, 0.10, 1.9, -35.0, 0.689380, 60.832000, 61.514100,
    ),
    MatrixStatRow::new(
        19, 1, 0.192, 0.029, 0.083, 2.3, -52.0, 0.696344, 84.019700, 84.985600,
    ),
    MatrixStatRow::new(
        18, 1, 0.183, 0.021, 0.070, 2.6, -60.0, 0.710525, 113.829000, 115.184000,
    ),
    MatrixStatRow::new(
        17, 1, 0.171, 0.014, 0.052, 3.3, -86.0, 0.727000, 175.071000, 177.196000,
    ),
];

const PAM30_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3, 0.436164, 0.161818, 0.161818,
    ),
    MatrixStatRow::new(
        7, 2, 0.305, 0.15, 0.87, 0.35, -3.0, 0.479087, 1.014010, 1.162730,
    ),
    MatrixStatRow::new(
        6, 2, 0.287, 0.11, 0.68, 0.42, -4.0, 0.499980, 1.688060, 1.951430,
    ),
    MatrixStatRow::new(
        5, 2, 0.264, 0.079, 0.45, 0.59, -7.0, 0.533009, 3.377010, 3.871950,
    ),
    MatrixStatRow::new(
        10, 1, 0.309, 0.15, 0.88, 0.35, -3.0, 0.474741, 1.372050, 1.788770,
    ),
    MatrixStatRow::new(
        9, 1, 0.294, 0.11, 0.61, 0.48, -6.0, 0.492716, 2.463920, 3.186150,
    ),
    MatrixStatRow::new(
        8, 1, 0.270, 0.072, 0.40, 0.68, -10.0, 0.521286, 5.368130, 6.763480,
    ),
    MatrixStatRow::new(
        15, 3, 0.339, 0.28, 1.70, 0.20, -0.5, 0.437688, 0.157089, 0.155299,
    ),
    MatrixStatRow::new(
        14, 2, 0.337, 0.27, 1.62, 0.21, -0.8, 0.440010, 0.206970, 0.198524,
    ),
    MatrixStatRow::new(
        14, 1, 0.333, 0.27, 1.43, 0.23, -1.4, 0.444817, 0.436301, 0.361947,
    ),
    MatrixStatRow::new(
        13, 3, 0.338, 0.27, 1.69, 0.20, -0.5, 0.439086, 0.178973, 0.175436,
    ),
];

const PAM70_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3345, 0.229, 1.029, 0.3250, -0.7, 0.511296, 0.633439, 0.633439,
    ),
    MatrixStatRow::new(
        8, 2, 0.301, 0.12, 0.54, 0.56, -5.0, 0.549019, 2.881650, 3.025710,
    ),
    MatrixStatRow::new(
        7, 2, 0.286, 0.093, 0.43, 0.67, -7.0, 0.565659, 4.534540, 4.785780,
    ),
    MatrixStatRow::new(
        6, 2, 0.264, 0.064, 0.29, 0.90, -12.0, 0.596330, 7.942630, 8.402720,
    ),
    MatrixStatRow::new(
        11, 1, 0.305, 0.12, 0.52, 0.59, -6.0, 0.543514, 3.681400, 4.108020,
    ),
    MatrixStatRow::new(
        10, 1, 0.291, 0.091, 0.41, 0.71, -9.0, 0.560723, 6.002970, 6.716570,
    ),
    MatrixStatRow::new(
        9, 1, 0.270, 0.060, 0.28, 0.97, -14.0, 0.585186, 11.360800, 12.636700,
    ),
    MatrixStatRow::new(
        11, 2, 0.323, 0.186, 0.80, 1.32, -27.0, 0.524062, 1.321301, 1.281671,
    ),
    MatrixStatRow::new(
        12, 3, 0.330, 0.219, 0.93, 0.82, -16.0, 0.516845, 0.818768, 0.811240,
    ),
];

const IDENTITY_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.28768, 0.282, 1.69, 0.1703, -0.3, 0.43828, 0.16804, 0.16804,
    ),
    MatrixStatRow::new(
        15, 2, 0.2835, 0.255, 1.49, 0.19, -1.0, 0.44502, 0.24613, 0.22743,
    ),
];

/// blast-rs: Resolve matrix-statistics rows by matrix name; not a direct NCBI
/// C port.
fn matrix_stat_rows(matrix_name: &str) -> Option<&'static [MatrixStatRow]> {
    if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        Some(BLOSUM45_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM50") {
        Some(BLOSUM50_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        Some(BLOSUM62_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        Some(BLOSUM80_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM90") {
        Some(BLOSUM90_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        Some(PAM30_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        Some(PAM70_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("PAM250") {
        Some(PAM250_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("IDENTITY") {
        Some(IDENTITY_ROWS)
    } else {
        None
    }
}

/// blast-rs: Matrix-statistics row lookup with standard-matrix filtering; not
/// a direct NCBI C port.
fn matrix_stat_rows_with_standard_only(
    matrix_name: &str,
    standard_only: bool,
) -> Option<&'static [MatrixStatRow]> {
    if standard_only && matrix_name.eq_ignore_ascii_case("IDENTITY") {
        None
    } else {
        matrix_stat_rows(matrix_name)
    }
}

/// blast-rs: Matrix-statistics row lookup by gap costs; not a direct NCBI C port.
fn matrix_stat_row_for_gap_costs(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
) -> Option<MatrixStatRow> {
    matrix_stat_rows(matrix_name)?
        .iter()
        .copied()
        .find(|row| row.gap_open == gap_open && row.gap_extend == gap_extend)
}

/// Precomputed Karlin-Altschul parameters for BLOSUM62 with various gap costs.
/// Values from NCBI blast_stat.c blosum62_values[].
/// Format: (gap_open, gap_extend, lambda, K, H, alpha, beta)
pub const BLOSUM62_PARAMS: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (11, 2, 0.297, 0.082, 0.27, 1.1, -10.0),
    (10, 2, 0.291, 0.075, 0.23, 1.3, -15.0),
    (9, 2, 0.279, 0.058, 0.19, 1.5, -19.0),
    (8, 2, 0.264, 0.045, 0.15, 1.8, -26.0),
    (7, 2, 0.239, 0.027, 0.10, 2.5, -46.0),
    (6, 2, 0.201, 0.012, 0.061, 3.3, -58.0),
    (13, 1, 0.292, 0.071, 0.23, 1.2, -11.0),
    (12, 1, 0.283, 0.059, 0.19, 1.5, -19.0),
    (11, 1, 0.267, 0.041, 0.14, 1.9, -30.0),
    (10, 1, 0.243, 0.024, 0.10, 2.5, -44.0),
    (9, 1, 0.206, 0.010, 0.052, 4.0, -87.0),
];

/// Precomputed Karlin-Altschul parameters for NCBI's IDENTITY matrix.
/// Values from `blast_stat.c` `prot_idenity_values[]`.
pub const IDENTITY_PARAMS: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (
        i16::MAX as i32,
        i16::MAX as i32,
        0.28768,
        0.282,
        1.69,
        0.1703,
        -0.3,
    ),
    (15, 2, 0.2835, 0.255, 1.49, 0.19, -1.0),
];

/// blast-rs: BLOSUM62 gapped-parameter lookup helper; not a direct NCBI C port.
pub fn lookup_protein_params(gap_open: i32, gap_extend: i32) -> Option<GappedParams> {
    lookup_matrix_params("BLOSUM62", gap_open, gap_extend)
}

/// blast-rs: Named-matrix gapped-parameter lookup helper; not a direct NCBI C port.
pub fn lookup_matrix_params(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
) -> Option<GappedParams> {
    matrix_stat_row_for_gap_costs(matrix_name, gap_open, gap_extend)
        .map(MatrixStatRow::gapped_params)
}

/// blast-rs: Named-matrix ungapped alpha/beta lookup helper; not a direct NCBI
/// C port.
pub fn lookup_matrix_ungapped_alpha_beta(matrix_name: &str) -> Option<(f64, f64)> {
    matrix_stat_rows(matrix_name)?
        .first()
        .map(|row| (row.alpha, row.beta))
}

/// blast-rs: Named-matrix ungapped display-stat lookup helper; not a direct
/// NCBI C port.
pub fn lookup_matrix_ungapped_display_params(matrix_name: &str) -> Option<ProteinMatrixStats> {
    matrix_stat_rows(matrix_name)?
        .first()
        .copied()
        .map(MatrixStatRow::display_params)
}

/// blast-rs: Named-matrix gapped display-stat lookup helper; not a direct NCBI
/// C port.
pub fn lookup_matrix_display_params(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
) -> Option<ProteinMatrixStats> {
    matrix_stat_row_for_gap_costs(matrix_name, gap_open, gap_extend)
        .map(MatrixStatRow::display_params)
}

/// blast-rs: IDENTITY-matrix gapped-parameter lookup helper; not a direct NCBI
/// C port.
pub fn lookup_identity_params(gap_open: i32, gap_extend: i32) -> Option<GappedParams> {
    lookup_matrix_params("IDENTITY", gap_open, gap_extend)
}

/// blast-rs: BLOSUM62 ungapped Karlin-block table helper; not a direct NCBI C port.
/// Values match NCBI
/// `blosum62_values[0]` (ungapped entry) from `blast_stat.c:259`:
/// Lambda = 0.3176, K = 0.134, H = 0.4012.
pub fn protein_ungapped_kbp() -> KarlinBlk {
    protein_ungapped_kbp_for_matrix("BLOSUM62")
}

/// blast-rs: Named-matrix ungapped Karlin-block table helper; not a direct
/// NCBI C port.
/// Reads from NCBI's standard
/// `blast_stat.c` row 0 table values.
pub fn protein_ungapped_kbp_for_matrix(matrix_name: &str) -> KarlinBlk {
    let row = matrix_stat_rows(matrix_name)
        .and_then(|rows| rows.first().copied())
        .unwrap_or(BLOSUM62_ROWS[0]);
    KarlinBlk {
        lambda: row.lambda,
        k: row.k,
        log_k: row.k.ln(),
        h: row.h,
        round_down: false,
    }
}

/// blast-rs: Protein ideal-Karlin computation helper; not a direct NCBI C port.
/// Computes NCBI's `sbp->kbp_ideal` for a protein matrix.
///
/// `Blast_ScoreBlkKbpIdealCalc` builds a score-frequency distribution from
/// standard residue composition on both axes and runs
/// `Blast_KarlinBlkUngappedCalc`; it does not copy the rounded row-0 matrix
/// table. Translated-query programs use this ideal block when the
/// query-specific lambda is too high.
pub fn protein_ideal_ungapped_kbp_for_matrix(matrix_name: &str) -> KarlinBlk {
    let Some(matrix) = protein_score_matrix_by_name(matrix_name) else {
        return protein_ungapped_kbp_for_matrix(matrix_name);
    };
    protein_ideal_ungapped_kbp_from_matrix(matrix, protein_ungapped_kbp_for_matrix(matrix_name))
}

/// blast-rs: Matrix-backed implementation for ideal protein Karlin blocks; not
/// a direct NCBI C port.
fn protein_ideal_ungapped_kbp_from_matrix(
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    fallback: KarlinBlk,
) -> KarlinBlk {
    let std_freq = protein_std_freq_ncbistdaa();
    let Some((lo, hi)) = protein_observed_score_bounds(matrix, &std_freq) else {
        return fallback;
    };

    let mut sfp = SfDist::new(lo, hi);
    for i in 0..crate::matrix::AA_SIZE {
        for j in 0..crate::matrix::AA_SIZE {
            let s = matrix[i][j];
            if (lo..=hi).contains(&s) {
                *sfp.p_mut(s) += std_freq[i] * std_freq[j];
            }
        }
    }

    let mut obs_min = i32::MAX;
    let mut obs_max = i32::MIN;
    let mut psum = 0.0;
    for s in lo..=hi {
        if sfp.p(s) > 0.0 {
            psum += sfp.p(s);
            obs_min = obs_min.min(s);
            obs_max = obs_max.max(s);
        }
    }
    sfp.obs_min = obs_min;
    sfp.obs_max = obs_max;
    let mut savg = 0.0;
    if psum > 0.0 {
        for s in obs_min..=obs_max {
            *sfp.p_mut(s) /= psum;
            savg += s as f64 * sfp.p(s);
        }
    }
    sfp.score_avg = savg;

    let lambda = compute_lambda(&sfp);
    let h = blast_karlin_lto_h(&sfp, lambda);
    let k = blast_karlin_lh_to_k(&sfp, lambda, h);
    if lambda < 0.0 || h < 0.0 || k < 0.0 {
        return fallback;
    }
    KarlinBlk {
        lambda,
        k,
        log_k: k.ln(),
        h,
        round_down: false,
    }
}

/// blast-rs: Resolve built-in protein matrix data by name; not a direct NCBI C port.
fn protein_score_matrix_by_name(
    matrix_name: &str,
) -> Option<&'static [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]> {
    if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        Some(&crate::matrix::BLOSUM45)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM50") {
        Some(&crate::matrix::BLOSUM50)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        Some(&crate::matrix::BLOSUM62)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        Some(&crate::matrix::BLOSUM80)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM90") {
        Some(&crate::matrix::BLOSUM90)
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        Some(&crate::matrix::PAM30)
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        Some(&crate::matrix::PAM70)
    } else if matrix_name.eq_ignore_ascii_case("PAM250") {
        Some(&crate::matrix::PAM250)
    } else if matrix_name.eq_ignore_ascii_case("IDENTITY") {
        Some(&crate::matrix::IDENTITY)
    } else {
        None
    }
}

/// blast-rs: Compute observed score bounds for a protein matrix and
/// composition; not a direct NCBI C port.
fn protein_observed_score_bounds(
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    std_freq: &[f64; crate::matrix::AA_SIZE],
) -> Option<(i32, i32)> {
    let mut lo = i32::MAX;
    let mut hi = i32::MIN;
    for i in 0..crate::matrix::AA_SIZE {
        if std_freq[i] == 0.0 {
            continue;
        }
        for j in 0..crate::matrix::AA_SIZE {
            if std_freq[j] == 0.0 {
                continue;
            }
            let s = matrix[i][j];
            lo = lo.min(s);
            hi = hi.max(s);
        }
    }
    if lo <= hi {
        Some((lo, hi))
    } else {
        None
    }
}

/// blast-rs: Query-specific ungapped protein Karlin helper; not a direct NCBI C port.
/// Computes the query-specific ungapped Karlin block for a protein query.
/// Mirrors NCBI's `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2737`):
/// `Blast_ResFreqString` over the query → `BlastScoreFreqCalc` →
/// `Blast_KarlinBlkUngappedCalc`. The result drifts slightly from the ideal
/// kbp depending on query composition; the drift is what makes `gap_trigger`
/// match NCBI bit-for-bit at boundary cases (see iter 99).
pub fn query_specific_protein_ungapped_kbp(
    query_aa_ncbistdaa: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
) -> KarlinBlk {
    query_specific_protein_ungapped_kbp_or_fallback(
        query_aa_ncbistdaa,
        matrix,
        protein_ungapped_kbp(),
    )
}

/// blast-rs: Query-specific ungapped protein Karlin helper with matrix fallback;
/// not a direct NCBI C port.
/// Computes query-specific ungapped protein KBP with the selected matrix's
/// standard row as the fallback for unsupported query compositions.
pub fn query_specific_protein_ungapped_kbp_for_matrix(
    query_aa_ncbistdaa: &[u8],
    matrix_name: &str,
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
) -> KarlinBlk {
    query_specific_protein_ungapped_kbp_or_fallback(
        query_aa_ncbistdaa,
        matrix,
        protein_ungapped_kbp_for_matrix(matrix_name),
    )
}

/// blast-rs: Shared query-specific ungapped protein Karlin implementation; not
/// a direct NCBI C port.
fn query_specific_protein_ungapped_kbp_or_fallback(
    query_aa_ncbistdaa: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    table_fallback: KarlinBlk,
) -> KarlinBlk {
    let std_freq = protein_std_freq_ncbistdaa();
    let Some((lo, hi)) = protein_observed_score_bounds(matrix, &std_freq) else {
        return table_fallback;
    };
    // Matrix-based protein setup calls `BLAST_ScoreSetAmbigRes(sbp, 'X')`
    // (`blast_setup.c:384`). `Blast_ResFreqString` therefore clears X only;
    // B/Z/U/*/O/J remain in the query composition exactly as NCBI sees them.
    let ambiguous: [u8; 1] = [crate::encoding::NCBISTDAA_X];
    let contexts = [UngappedKbpContext {
        query_offset: 0,
        query_length: query_aa_ncbistdaa.len() as i32,
        is_valid: true,
    }];
    let mat = |i: usize, j: usize| matrix[i][j];
    let r = ungapped_kbp_calc_with_std(
        query_aa_ncbistdaa,
        &contexts,
        lo,
        hi,
        crate::matrix::AA_SIZE,
        &ambiguous,
        &std_freq,
        &mat,
    );
    let fallback = protein_ideal_ungapped_kbp_from_matrix(matrix, table_fallback);
    r.into_iter().next().flatten().unwrap_or(fallback)
}

/// blast-rs: NCBIstdaa-indexed standard amino-acid background-frequency table;
/// not a direct NCBI C port.
/// NCBIstdaa-indexed standard amino-acid background frequencies.
/// Robinson & Robinson 1991 (matches NCBI's `Blast_ResFreqStdComp` table for
/// `BLASTAA_SEQ_CODE`). Entries for `*`, `-`, ambiguity letters (B/Z/X/U/O/J)
/// are zero — the C path leaves them out of the score-frequency calculation.
pub fn protein_std_freq_ncbistdaa() -> [f64; 28] {
    use crate::matrix::AA_FREQUENCIES;
    // AA_FREQUENCIES is alphabetical (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y).
    // NCBIstdaa indices: -=0 A=1 B=2 C=3 D=4 E=5 F=6 G=7 H=8 I=9 K=10 L=11 M=12
    // N=13 P=14 Q=15 R=16 S=17 T=18 V=19 W=20 X=21 Y=22 Z=23 U=24 *=25 O=26 J=27
    let mut out = [0.0f64; 28];
    for (freq, &code) in AA_FREQUENCIES
        .iter()
        .zip(crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter())
    {
        out[code as usize] = *freq;
    }
    out
}

/// blast-rs: Letter-indexed standard amino-acid background-frequency table;
/// not a direct NCBI C port.
fn protein_std_freq_letters() -> [(u8, f64); 20] {
    [
        (b'A', 78.05),
        (b'C', 19.25),
        (b'D', 53.64),
        (b'E', 62.95),
        (b'F', 38.56),
        (b'G', 73.77),
        (b'H', 21.99),
        (b'I', 51.42),
        (b'K', 57.44),
        (b'L', 90.19),
        (b'M', 22.43),
        (b'N', 44.87),
        (b'P', 52.03),
        (b'Q', 42.64),
        (b'R', 51.29),
        (b'S', 71.20),
        (b'T', 58.41),
        (b'V', 64.41),
        (b'W', 13.30),
        (b'Y', 32.16),
    ]
}

/// blast-rs: Effective-search-space arithmetic helper; not a direct NCBI C port.
/// eff_length = db_length - num_seqs * length_adjustment
/// eff_query_length = query_length - length_adjustment
/// search_space = eff_length * eff_query_length
pub fn compute_search_space(
    query_length: i64,
    db_length: i64,
    num_seqs: i32,
    length_adjustment: i32,
) -> f64 {
    let eff_query = (query_length - length_adjustment as i64).max(1);
    let eff_db = (db_length - num_seqs as i64 * length_adjustment as i64).max(1);
    eff_query as f64 * eff_db as f64
}

const TARGET_HIT_PROB: f64 = 0.98;

/// Rust owner for NCBI `MatrixData` scratch state in `blast_tune.c`.
#[derive(Debug, Clone, Default)]
pub struct MatrixData {
    pub matrix_dim_alloc: i32,
    pub matrix_dim: i32,
    pub hit_probability: f64,
    pub percent_identity: f64,
    pub power_matrix: Vec<f64>,
    pub prod_matrix: Vec<f64>,
}

/// NCBI: s_MatrixDataReset (blast_tune.c:125).
pub fn s_matrix_data_reset(
    m: Option<&mut MatrixData>,
    new_word_size: i32,
    percent_identity: f64,
) -> i16 {
    let Some(m) = m else {
        return -1;
    };

    m.hit_probability = 0.0;
    m.percent_identity = percent_identity;
    m.matrix_dim = new_word_size + 1;

    if m.matrix_dim > m.matrix_dim_alloc {
        let Some(num_cells) = (m.matrix_dim as usize).checked_mul(m.matrix_dim as usize) else {
            m.power_matrix.clear();
            m.prod_matrix.clear();
            return -2;
        };
        m.matrix_dim_alloc = m.matrix_dim;
        m.power_matrix.resize(num_cells, 0.0);
        m.prod_matrix.resize(num_cells, 0.0);
    }
    0
}

/// NCBI: s_SetInitialMatrix (blast_tune.c:178).
fn s_set_initial_matrix(matrix: &mut [f64], matrix_dim: i32, identity: f64) {
    let dim = matrix_dim.max(0) as usize;
    matrix[..dim * dim].fill(0.0);
    for i in 0..dim.saturating_sub(1) {
        let row = i * dim;
        matrix[row] = 1.0 - identity;
        matrix[row + i + 1] = identity;
    }
    if dim > 0 {
        matrix[(dim - 1) * dim + dim - 1] = 1.0;
    }
}

/// NCBI: s_MatrixMultiply (blast_tune.c:210).
fn s_matrix_multiply(a: &[f64], identity: f64, prod: &mut [f64], dim: i32) {
    let dim = dim.max(0) as usize;
    let comp_identity = 1.0 - identity;

    for i in 0..dim {
        let row = i * dim;
        let accum: f64 = a[row..row + dim.saturating_sub(1)].iter().sum();
        prod[row] = comp_identity * accum;
    }

    for i in 0..dim {
        let row = i * dim;
        for j in 1..dim {
            prod[row + j] = identity * a[row + j - 1];
        }
    }

    if dim > 0 {
        for i in 0..dim {
            let row = i * dim;
            prod[row + dim - 1] += a[row + dim - 1];
        }
    }
}

/// NCBI: s_MatrixSquare (blast_tune.c:243).
pub fn s_matrix_square(a: &[f64], prod: &mut [f64], dim: i32) {
    let dim = dim.max(0) as usize;
    let full_entries = dim & !3;

    for i in 0..dim {
        let row = i * dim;
        for j in 0..dim {
            let mut accum = 0.0;
            let mut k = 0;
            while k < full_entries {
                accum += a[row + k] * a[j + k * dim]
                    + a[row + k + 1] * a[j + (k + 1) * dim]
                    + a[row + k + 2] * a[j + (k + 2) * dim]
                    + a[row + k + 3] * a[j + (k + 3) * dim];
                k += 4;
            }
            while k < dim {
                accum += a[row + k] * a[j + k * dim];
                k += 1;
            }
            prod[row + j] = accum;
        }
    }
}

/// NCBI: s_FindHitProbability (blast_tune.c:293).
fn s_find_hit_probability(
    m: &mut MatrixData,
    word_size: i32,
    min_percent_identity: f64,
    min_align_length: i32,
) -> i16 {
    if min_align_length == 0 {
        return -3;
    }
    if s_matrix_data_reset(Some(m), word_size, min_percent_identity) != 0 {
        return -4;
    }

    let dim = m.matrix_dim as usize;
    s_set_initial_matrix(&mut m.power_matrix, m.matrix_dim, min_percent_identity);

    let mut mask = 0x8000_0000u32;
    while (min_align_length as u32 & mask) == 0 {
        mask /= 2;
    }

    mask /= 2;
    let mut num_squares = 0;
    while mask != 0 {
        if num_squares == 0 {
            s_matrix_multiply(
                &m.power_matrix[..dim * dim],
                m.percent_identity,
                &mut m.prod_matrix[..dim * dim],
                m.matrix_dim,
            );
        } else {
            s_matrix_square(
                &m.power_matrix[..dim * dim],
                &mut m.prod_matrix[..dim * dim],
                m.matrix_dim,
            );
        }
        std::mem::swap(&mut m.prod_matrix, &mut m.power_matrix);

        if (min_align_length as u32 & mask) != 0 {
            s_matrix_multiply(
                &m.power_matrix[..dim * dim],
                m.percent_identity,
                &mut m.prod_matrix[..dim * dim],
                m.matrix_dim,
            );
            std::mem::swap(&mut m.prod_matrix, &mut m.power_matrix);
        }

        mask /= 2;
        num_squares += 1;
    }

    m.hit_probability = m.power_matrix[dim - 1];
    0
}

/// NCBI: s_FindWordSize (blast_tune.c:362).
pub fn s_find_word_size(
    m: Option<&mut MatrixData>,
    min_percent_identity: f64,
    min_align_length: i32,
) -> i32 {
    let Some(m) = m else {
        return 0;
    };

    let k_min_w: f64 = 4.0;
    let k_max_w: f64 = 110.0;

    let mut w1 = 28.0;
    if s_find_hit_probability(m, (w1 + 0.5) as i32, min_percent_identity, min_align_length) != 0 {
        return 0;
    }
    let mut p1 = m.hit_probability - TARGET_HIT_PROB;

    let mut w0 = 11.0;
    if s_find_hit_probability(m, (w0 + 0.5) as i32, min_percent_identity, min_align_length) != 0 {
        return 0;
    }
    let mut p0 = m.hit_probability - TARGET_HIT_PROB;

    if p1 > 0.0 {
        while p1 > 0.0 && w1 < k_max_w {
            w0 = w1;
            p0 = p1;
            w1 = (2.0 * w1).min(k_max_w);
            if s_find_hit_probability(m, (w1 + 0.5) as i32, min_percent_identity, min_align_length)
                != 0
            {
                return 0;
            }
            p1 = m.hit_probability - TARGET_HIT_PROB;
        }

        if p1 > 0.0 {
            return (w1 + 0.5) as i32;
        }
    } else if p0 < 0.0 {
        w1 = w0;
        p1 = p0;
        w0 = k_min_w;
        if s_find_hit_probability(m, (w0 + 0.5) as i32, min_percent_identity, min_align_length) != 0
        {
            return 0;
        }
        p0 = m.hit_probability - TARGET_HIT_PROB;

        if p0 < 0.0 {
            return (w0 + 0.5) as i32;
        }
    }

    while (w1 - w0).abs() > 1.0 {
        let w2 = (w0 + w1) / 2.0;
        if s_find_hit_probability(m, (w2 + 0.5) as i32, min_percent_identity, min_align_length) != 0
        {
            return 0;
        }
        let p2 = m.hit_probability - TARGET_HIT_PROB;

        if p2 > 0.0 {
            w0 = w2;
            p0 = p2;
        } else {
            w1 = w2;
            p1 = p2;
        }
    }

    let _ = (p0, p1);
    (w0 + 0.5) as i32
}

/// Port of NCBI `BLAST_FindBestNucleotideWordSize` (`blast_tune.c:468`).
pub fn blast_find_best_nucleotide_word_size(
    min_percent_identity: f64,
    mut min_align_length: i32,
) -> i32 {
    if !(0.6..1.0).contains(&min_percent_identity) {
        return 0;
    }
    if min_align_length > 10_000 {
        min_align_length = 10_000;
    } else if min_align_length < 0 {
        return 0;
    } else if min_align_length < 8 {
        return 4;
    }

    let mut matrix_data = MatrixData::default();
    s_find_word_size(
        Some(&mut matrix_data),
        min_percent_identity,
        min_align_length,
    )
}

// ---------------------------------------------------------------------------
// Gapped Karlin-Altschul parameter lookup (exact C-compatible tables)
// Port of Blast_KarlinBlkNuclGappedCalc from blast_stat.c
// ---------------------------------------------------------------------------

/// Each row: [gap_open, gap_extend, lambda, K, H, alpha, beta, theta]
pub type KbpTableRow = [f64; 8];

const KBPT_1_5: &[KbpTableRow] = &[
    [0.0, 0.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0],
    [3.0, 3.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0],
];
const KBPT_1_4: &[KbpTableRow] = &[
    [0.0, 0.0, 1.383, 0.738, 1.36, 1.02, 0.0, 100.0],
    [1.0, 2.0, 1.36, 0.67, 1.2, 1.1, 0.0, 98.0],
    [0.0, 2.0, 1.26, 0.43, 0.90, 1.4, -1.0, 91.0],
    [2.0, 1.0, 1.35, 0.61, 1.1, 1.2, -1.0, 98.0],
    [1.0, 1.0, 1.22, 0.35, 0.72, 1.7, -3.0, 88.0],
];
const KBPT_2_7: &[KbpTableRow] = &[
    [0.0, 0.0, 0.69, 0.73, 1.34, 0.515, 0.0, 100.0],
    [2.0, 4.0, 0.68, 0.67, 1.2, 0.55, 0.0, 99.0],
    [0.0, 4.0, 0.63, 0.43, 0.90, 0.7, -1.0, 91.0],
    [4.0, 2.0, 0.675, 0.62, 1.1, 0.6, -1.0, 98.0],
    [2.0, 2.0, 0.61, 0.35, 0.72, 1.7, -3.0, 88.0],
];
const KBPT_1_3: &[KbpTableRow] = &[
    [0.0, 0.0, 1.374, 0.711, 1.31, 1.05, 0.0, 100.0],
    [2.0, 2.0, 1.37, 0.70, 1.2, 1.1, 0.0, 99.0],
    [1.0, 2.0, 1.35, 0.64, 1.1, 1.2, -1.0, 98.0],
    [0.0, 2.0, 1.25, 0.42, 0.83, 1.5, -2.0, 91.0],
    [2.0, 1.0, 1.34, 0.60, 1.1, 1.2, -1.0, 97.0],
    [1.0, 1.0, 1.21, 0.34, 0.71, 1.7, -2.0, 88.0],
];
const KBPT_2_5: &[KbpTableRow] = &[
    [0.0, 0.0, 0.675, 0.65, 1.1, 0.6, -1.0, 99.0],
    [2.0, 4.0, 0.67, 0.59, 1.1, 0.6, -1.0, 98.0],
    [0.0, 4.0, 0.62, 0.39, 0.78, 0.8, -2.0, 91.0],
    [4.0, 2.0, 0.67, 0.61, 1.0, 0.65, -2.0, 98.0],
    [2.0, 2.0, 0.56, 0.32, 0.59, 0.95, -4.0, 82.0],
];
const KBPT_1_2: &[KbpTableRow] = &[
    [0.0, 0.0, 1.28, 0.46, 0.85, 1.5, -2.0, 96.0],
    [2.0, 2.0, 1.33, 0.62, 1.1, 1.2, 0.0, 99.0],
    [1.0, 2.0, 1.30, 0.52, 0.93, 1.4, -2.0, 97.0],
    [0.0, 2.0, 1.19, 0.34, 0.66, 1.8, -3.0, 89.0],
    [3.0, 1.0, 1.32, 0.57, 1.0, 1.3, -1.0, 99.0],
    [2.0, 1.0, 1.29, 0.49, 0.92, 1.4, -1.0, 96.0],
    [1.0, 1.0, 1.14, 0.26, 0.52, 2.2, -5.0, 85.0],
];
const KBPT_2_3: &[KbpTableRow] = &[
    [0.0, 0.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0],
    [4.0, 4.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0],
    [2.0, 4.0, 0.615, 0.37, 0.72, 0.85, -3.0, 97.0],
    [0.0, 4.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0],
    [3.0, 3.0, 0.615, 0.37, 0.68, 0.9, -3.0, 97.0],
    [6.0, 2.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0],
    [5.0, 2.0, 0.625, 0.41, 0.78, 0.8, -2.0, 99.0],
    [4.0, 2.0, 0.61, 0.35, 0.68, 0.9, -3.0, 96.0],
    [2.0, 2.0, 0.515, 0.14, 0.33, 1.55, -9.0, 81.0],
];
const KBPT_3_4: &[KbpTableRow] = &[
    [6.0, 3.0, 0.389, 0.25, 0.56, 0.7, -5.0, 95.0],
    [5.0, 3.0, 0.375, 0.21, 0.47, 0.8, -6.0, 92.0],
    [4.0, 3.0, 0.351, 0.14, 0.35, 1.0, -9.0, 86.0],
    [6.0, 2.0, 0.362, 0.16, 0.45, 0.8, -4.0, 88.0],
    [5.0, 2.0, 0.330, 0.092, 0.28, 1.2, -13.0, 81.0],
    [4.0, 2.0, 0.281, 0.046, 0.16, 1.8, -23.0, 69.0],
];
const KBPT_1_1: &[KbpTableRow] = &[
    [3.0, 2.0, 1.09, 0.31, 0.55, 2.0, -2.0, 99.0],
    [2.0, 2.0, 1.07, 0.27, 0.49, 2.2, -3.0, 97.0],
    [1.0, 2.0, 1.02, 0.21, 0.36, 2.8, -6.0, 92.0],
    [0.0, 2.0, 0.80, 0.064, 0.17, 4.8, -16.0, 72.0],
    [4.0, 1.0, 1.08, 0.28, 0.54, 2.0, -2.0, 98.0],
    [3.0, 1.0, 1.06, 0.25, 0.46, 2.3, -4.0, 96.0],
    [2.0, 1.0, 0.99, 0.17, 0.30, 3.3, -10.0, 90.0],
];
const KBPT_4_5: &[KbpTableRow] = &[
    [0.0, 0.0, 0.22, 0.061, 0.22, 1.0, -15.0, 74.0],
    [6.0, 5.0, 0.28, 0.21, 0.47, 0.6, -7.0, 93.0],
    [5.0, 5.0, 0.27, 0.17, 0.39, 0.7, -9.0, 90.0],
    [4.0, 5.0, 0.25, 0.10, 0.31, 0.8, -10.0, 83.0],
    [3.0, 5.0, 0.23, 0.065, 0.25, 0.9, -11.0, 76.0],
];
const KBPT_3_2: &[KbpTableRow] = &[[5.0, 5.0, 0.208, 0.030, 0.072, 2.9, -47.0, 77.0]];
const KBPT_5_4: &[KbpTableRow] = &[
    [10.0, 6.0, 0.163, 0.068, 0.16, 1.0, -19.0, 85.0],
    [8.0, 6.0, 0.146, 0.039, 0.11, 1.3, -29.0, 76.0],
];

struct KbpTableMeta {
    table: &'static [KbpTableRow],
    gap_open_max: i32,
    gap_extend_max: i32,
    round_down: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MatrixInfo {
    pub name: String,
    pub values: &'static [MatrixStatRow],
    pub prefs: Vec<i32>,
    pub max_number_values: i32,
}

/// NCBI: MatrixInfoNew (blast_stat.c:2878).
pub fn matrix_info_new(name: &str) -> Option<MatrixInfo> {
    let values = matrix_stat_rows(name)?;
    Some(MatrixInfo {
        name: name.to_string(),
        values,
        prefs: Vec::new(),
        max_number_values: values.len() as i32,
    })
}

/// NCBI: MatrixInfoDestruct (blast_stat.c:2890).
pub fn matrix_info_destruct(matrix_info: &mut Option<MatrixInfo>) -> Option<MatrixInfo> {
    if let Some(info) = matrix_info.as_mut() {
        info.name.clear();
        info.prefs.clear();
        info.max_number_values = 0;
    }
    *matrix_info = None;
    None
}

const STANDARD_MATRIX_NAMES: &[&str] = &[
    "BLOSUM80", "BLOSUM62", "BLOSUM50", "BLOSUM45", "PAM250", "BLOSUM90", "PAM30", "PAM70",
];

const ALL_MATRIX_NAMES: &[&str] = &[
    "BLOSUM80", "BLOSUM62", "BLOSUM50", "BLOSUM45", "PAM250", "BLOSUM90", "PAM30", "PAM70",
    "IDENTITY",
];

/// NCBI: BlastLoadMatrixValues (blast_stat.c:2952).
/// naming: Returns an owned `Vec<MatrixInfo>` instead of a C linked list.
pub fn blast_load_matrix_values(standard_only: bool) -> Vec<MatrixInfo> {
    let names = if standard_only {
        STANDARD_MATRIX_NAMES
    } else {
        ALL_MATRIX_NAMES
    };
    names
        .iter()
        .filter_map(|name| matrix_info_new(name))
        .collect()
}

/// NCBI: BlastMatrixValuesDestruct (blast_stat.c:2929).
/// naming: Clears the owned `Vec<MatrixInfo>` instead of freeing a C linked list.
pub fn blast_matrix_values_destruct(values: &mut Vec<MatrixInfo>) -> Vec<MatrixInfo> {
    values.clear();
    Vec::new()
}

/// NCBI: BLAST_PrintMatrixMessage (blast_stat.c:3760).
/// naming: Returns a Rust `String` instead of writing to a C message buffer.
pub fn blast_print_matrix_message(matrix_name: &str, standard_only: bool) -> String {
    let mut out = format!("{matrix_name} is not a supported matrix, supported matrices are:\n");
    for info in blast_load_matrix_values(standard_only) {
        out.push_str(&format!("{} \n", info.name));
    }
    out
}

#[derive(Debug, Clone, PartialEq)]
pub struct MatrixValues {
    pub open: Vec<i32>,
    pub extension: Vec<i32>,
    pub lambda: Vec<f64>,
    pub k: Vec<f64>,
    pub h: Vec<f64>,
    pub alpha: Vec<f64>,
    pub beta: Vec<f64>,
    pub pref_flags: Vec<i32>,
}

/// blast-rs: Matrix-name to preferred-row index lookup; not a direct NCBI C port.
fn preferred_matrix_row_index(matrix_name: &str, rows_len: usize) -> Option<usize> {
    let index = if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        7
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM50") {
        9
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        9
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        8
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM90") {
        6
    } else if matrix_name.eq_ignore_ascii_case("PAM250") {
        8
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        5
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        5
    } else if matrix_name.eq_ignore_ascii_case("IDENTITY") {
        1
    } else {
        return None;
    };
    (index < rows_len).then_some(index)
}

/// NCBI: Blast_GetMatrixValues (blast_stat.c:3013).
/// naming: Returns owned parallel vectors instead of C output arrays.
pub fn blast_get_matrix_values(matrix_name: Option<&str>) -> MatrixValues {
    let Some(matrix_name) = matrix_name else {
        return MatrixValues {
            open: Vec::new(),
            extension: Vec::new(),
            lambda: Vec::new(),
            k: Vec::new(),
            h: Vec::new(),
            alpha: Vec::new(),
            beta: Vec::new(),
            pref_flags: Vec::new(),
        };
    };
    let Some(rows) = matrix_stat_rows(matrix_name) else {
        return MatrixValues {
            open: Vec::new(),
            extension: Vec::new(),
            lambda: Vec::new(),
            k: Vec::new(),
            h: Vec::new(),
            alpha: Vec::new(),
            beta: Vec::new(),
            pref_flags: Vec::new(),
        };
    };
    let best_index = preferred_matrix_row_index(matrix_name, rows.len());
    MatrixValues {
        open: rows.iter().map(|row| row.gap_open).collect(),
        extension: rows.iter().map(|row| row.gap_extend).collect(),
        lambda: rows.iter().map(|row| row.lambda).collect(),
        k: rows.iter().map(|row| row.k).collect(),
        h: rows.iter().map(|row| row.h).collect(),
        alpha: rows.iter().map(|row| row.alpha).collect(),
        beta: rows.iter().map(|row| row.beta).collect(),
        pref_flags: rows
            .iter()
            .enumerate()
            .map(|(index, _)| {
                if Some(index) == best_index {
                    BLAST_MATRIX_BEST
                } else {
                    BLAST_MATRIX_NOMINAL
                }
            })
            .collect(),
    }
}

/// NCBI: BLAST_GetAlphaBeta (blast_stat.c:3094).
/// Port for protein
/// matrix rows. Nucleotide callers use [`blast_get_nucl_alpha_beta`].
pub fn blast_get_alpha_beta(
    matrix_name: Option<&str>,
    alpha: &mut f64,
    beta: &mut f64,
    gapped: bool,
    gap_open: i32,
    gap_extend: i32,
    kbp_ungapped: Option<&KarlinBlk>,
) {
    let values = blast_get_matrix_values(matrix_name);
    let num_values = values.open.len();
    if gapped {
        if gap_open == 0 && gap_extend == 0 {
            for index in 1..num_values {
                if values.pref_flags[index] == BLAST_MATRIX_BEST {
                    *alpha = values.alpha[index];
                    *beta = values.beta[index];
                    break;
                }
            }
        } else {
            for index in 1..num_values {
                if values.open[index] == gap_open && values.extension[index] == gap_extend {
                    *alpha = values.alpha[index];
                    *beta = values.beta[index];
                    break;
                }
            }
        }
    } else if num_values > 0 {
        *alpha = values.alpha[0];
        *beta = values.beta[0];
    } else if let Some(kbp) = kbp_ungapped {
        *alpha = kbp.lambda / kbp.h;
        *beta = 0.0;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NuclValuesArray {
    pub normal: Vec<KbpTableRow>,
    pub non_affine: Vec<KbpTableRow>,
    pub gap_open_max: i32,
    pub gap_extend_max: i32,
    pub round_down: bool,
}

/// NCBI `kSmallFloat` (`blast_stat.c:4049`): smallest positive e-value
/// used as a clamp in `BlastKarlinEtoS_simple` to prevent FP exceptions.
pub const K_SMALL_FLOAT: f64 = 1.0e-297;

/// NCBI `NUM_FRAMES` (`blast_def.h:88`): number of translation frames
/// for protein searches (6 = ±3).
pub const NUM_FRAMES: usize = 6;

/// NCBI `CODON_LENGTH` (`blast_def.h:63`): nucleotides per amino acid.
pub const CODON_LENGTH: usize = 3;

/// NCBI `DEFAULT_LONGEST_INTRON` (`blast_def.h:78`): default longest
/// intron length (nt) for linking spliced HSPs in translated searches.
pub const DEFAULT_LONGEST_INTRON: usize = 122;

/// NCBI `COMPRESSION_RATIO` (`blast_def.h:83`): NCBI2na packs 4 bases
/// per byte.
pub const COMPRESSION_RATIO: usize = 4;

/// NCBI `NUM_STRANDS` (`blast_def.h:93`): plus / minus.
pub const NUM_STRANDS: usize = 2;

/// NCBI `GENCODE_STRLEN` (`blast_def.h:98`): fixed 64 codons in the
/// genetic-code table (all unique `NNN` triplets in NCBIstdaa).
pub const GENCODE_STRLEN: usize = 64;

/// NCBI: BLAST_SCORE_MIN (blast_stat.h:121).
/// Minimum allowed
/// score (one-letter comparison). NCBI defines it as `INT2_MIN`
/// (-32768); we keep the value literally so downstream "impossible
/// score" sentinels compare identically.
pub const BLAST_SCORE_MIN: i32 = i16::MIN as i32;

/// NCBI: BLAST_SCORE_MAX (blast_stat.h:122).
/// Maximum allowed
/// score (one-letter comparison) = `INT2_MAX` (32767).
pub const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

/// NCBI `BLAST_MATRIX_NOMINAL` (`blast_stat.h:53`): matrix quality flag
/// for "acceptable values, not recommended".
pub const BLAST_MATRIX_NOMINAL: i32 = 0;
/// NCBI `BLAST_MATRIX_PREFERRED` (`blast_stat.h:54`): preferred values.
pub const BLAST_MATRIX_PREFERRED: i32 = 1;
/// NCBI `BLAST_MATRIX_BEST` (`blast_stat.h:55`): best value (one per matrix).
pub const BLAST_MATRIX_BEST: i32 = 2;

/// `DBSEQ_CHUNK_OVERLAP` (`blast_hits.h:192`): overlap between adjacent
/// chunks when a subject is split for parallel/restarted scanning.
pub const DBSEQ_CHUNK_OVERLAP: usize = 100;

/// `PV_ARRAY_BYTES` (`blast_lookup.h:42`): bytes in the PV-array native word.
pub const PV_ARRAY_BYTES: usize = 4;
/// `PV_ARRAY_BTS` (`blast_lookup.h:43`): bits-to-shift from lookup
/// index to pv-array index (word size 32 → shift 5).
pub const PV_ARRAY_BTS: usize = 5;
/// `PV_ARRAY_MASK` (`blast_lookup.h:44`): mask off low 5 bits to get
/// the bit-within-word offset.
pub const PV_ARRAY_MASK: u32 = 31;

/// `BLAST_SEQSRC_MINGAP` (`blast_seqsrc.h:203`): minimal gap allowed in
/// range list for a sequence source chunk reader.
pub const BLAST_SEQSRC_MINGAP: usize = 1024;
/// `BLAST_SEQSRC_OVERHANG` (`blast_seqsrc.h:204`): extension applied to
/// each new range added to the seq-src range list.
pub const BLAST_SEQSRC_OVERHANG: usize = 1024;
/// `BLAST_SEQSRC_MINLENGTH` (`blast_seqsrc.h:205`): minimum sequence
/// length a seq-src returns.
pub const BLAST_SEQSRC_MINLENGTH: usize = 10;

/// `BLAST_SEQSRC_EXCLUDED` (`blast_seqsrc.h:290`): seq-src status code —
/// sequence excluded due to filtering.
pub const BLAST_SEQSRC_EXCLUDED: i32 = -3;
/// `BLAST_SEQSRC_ERROR` (`blast_seqsrc.h:291`): generic error retrieving sequence.
pub const BLAST_SEQSRC_ERROR: i32 = -2;
/// `BLAST_SEQSRC_EOF` (`blast_seqsrc.h:292`): no more sequences available.
pub const BLAST_SEQSRC_EOF: i32 = -1;
/// `BLAST_SEQSRC_SUCCESS` (`blast_seqsrc.h:293`): successful retrieval.
pub const BLAST_SEQSRC_SUCCESS: i32 = 0;

// Program-type bitmask flags (`blast_program.h:48-64`).
/// `PROTEIN_QUERY_MASK`.
pub const PROTEIN_QUERY_MASK: u32 = 1 << 0;
/// `PROTEIN_SUBJECT_MASK`.
pub const PROTEIN_SUBJECT_MASK: u32 = 1 << 1;
/// `NUCLEOTIDE_QUERY_MASK`.
pub const NUCLEOTIDE_QUERY_MASK: u32 = 1 << 2;
/// `NUCLEOTIDE_SUBJECT_MASK`.
pub const NUCLEOTIDE_SUBJECT_MASK: u32 = 1 << 3;
/// `TRANSLATED_QUERY_MASK`.
pub const TRANSLATED_QUERY_MASK: u32 = 1 << 4;
/// `TRANSLATED_SUBJECT_MASK`.
pub const TRANSLATED_SUBJECT_MASK: u32 = 1 << 5;
/// `PSSM_QUERY_MASK`.
pub const PSSM_QUERY_MASK: u32 = 1 << 6;
/// `PSSM_SUBJECT_MASK`.
pub const PSSM_SUBJECT_MASK: u32 = 1 << 7;
/// `PATTERN_QUERY_MASK`.
pub const PATTERN_QUERY_MASK: u32 = 1 << 8;

/// blast-rs: Nucleotide Karlin table selector; not a direct NCBI C port.
fn nucleotide_kbp_table_meta(reward: i32, penalty: i32) -> Option<KbpTableMeta> {
    match (reward, penalty) {
        (1, -5) => Some(KbpTableMeta {
            table: KBPT_1_5,
            gap_open_max: 3,
            gap_extend_max: 3,
            round_down: false,
        }),
        (1, -4) => Some(KbpTableMeta {
            table: KBPT_1_4,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -7) => Some(KbpTableMeta {
            table: KBPT_2_7,
            gap_open_max: 4,
            gap_extend_max: 4,
            round_down: true,
        }),
        (1, -3) => Some(KbpTableMeta {
            table: KBPT_1_3,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -5) => Some(KbpTableMeta {
            table: KBPT_2_5,
            gap_open_max: 4,
            gap_extend_max: 4,
            round_down: true,
        }),
        (1, -2) => Some(KbpTableMeta {
            table: KBPT_1_2,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -3) => Some(KbpTableMeta {
            table: KBPT_2_3,
            gap_open_max: 6,
            gap_extend_max: 4,
            round_down: true,
        }),
        (3, -4) => Some(KbpTableMeta {
            table: KBPT_3_4,
            gap_open_max: 6,
            gap_extend_max: 3,
            round_down: true,
        }),
        (1, -1) => Some(KbpTableMeta {
            table: KBPT_1_1,
            gap_open_max: 4,
            gap_extend_max: 2,
            round_down: false,
        }),
        (3, -2) => Some(KbpTableMeta {
            table: KBPT_3_2,
            gap_open_max: 5,
            gap_extend_max: 5,
            round_down: false,
        }),
        (4, -5) => Some(KbpTableMeta {
            table: KBPT_4_5,
            gap_open_max: 12,
            gap_extend_max: 8,
            round_down: false,
        }),
        (5, -4) => Some(KbpTableMeta {
            table: KBPT_5_4,
            gap_open_max: 25,
            gap_extend_max: 10,
            round_down: false,
        }),
        _ => None,
    }
}

/// NCBI: s_SplitArrayOf8 (blast_stat.c:3152).
pub fn s_split_array_of8(input: &[KbpTableRow]) -> (i16, Vec<KbpTableRow>, Vec<KbpTableRow>, bool) {
    if input.is_empty() {
        return (1, Vec::new(), Vec::new(), false);
    }
    if input[0][0] == 0.0 && input[0][1] == 0.0 {
        (0, input[1..].to_vec(), vec![input[0]], true)
    } else {
        (0, input.to_vec(), Vec::new(), false)
    }
}

/// NCBI: s_AdjustGapParametersByGcd (blast_stat.c:3186).
pub fn s_adjust_gap_parameters_by_gcd(
    normal: &mut [KbpTableRow],
    non_affine: &mut [KbpTableRow],
    gap_open_max: &mut i32,
    gap_extend_max: &mut i32,
    divisor: i32,
) -> i16 {
    if divisor <= 0 {
        return 1;
    }
    if divisor == 1 {
        return 0;
    }
    for row in normal.iter_mut().chain(non_affine.iter_mut()) {
        row[0] *= divisor as f64;
        row[1] *= divisor as f64;
        row[2] /= divisor as f64;
        row[5] /= divisor as f64;
    }
    *gap_open_max *= divisor;
    *gap_extend_max *= divisor;
    0
}

/// NCBI: s_GetNuclValuesArray (blast_stat.c:3238).
/// naming: Returns owned vectors in `NuclValuesArray` instead of C output arrays.
pub fn s_get_nucl_values_array(reward: i32, penalty: i32) -> Result<NuclValuesArray, i16> {
    let divisor = crate::math::blast_gcd(reward, penalty.abs());
    if divisor <= 0 {
        return Err(-1);
    }
    let nr = reward / divisor;
    let np = penalty / divisor;
    let Some(meta) = nucleotide_kbp_table_meta(nr, np) else {
        return Err(-1);
    };
    let (status, mut normal, mut non_affine, _) = s_split_array_of8(meta.table);
    if status != 0 {
        return Err(status);
    }
    let mut gap_open_max = meta.gap_open_max;
    let mut gap_extend_max = meta.gap_extend_max;
    let status = s_adjust_gap_parameters_by_gcd(
        &mut normal,
        &mut non_affine,
        &mut gap_open_max,
        &mut gap_extend_max,
        divisor,
    );
    if status != 0 {
        return Err(status);
    }
    Ok(NuclValuesArray {
        normal,
        non_affine,
        gap_open_max,
        gap_extend_max,
        round_down: meta.round_down,
    })
}

/// NCBI: BLAST_GetNucleotideGapExistenceExtendParams (blast_stat.c:3402).
pub fn blast_get_nucleotide_gap_existence_extend_params(
    reward: i32,
    penalty: i32,
    gap_existence: &mut i32,
    gap_extension: &mut i32,
) -> i16 {
    let values = match s_get_nucl_values_array(reward, penalty) {
        Ok(values) => values,
        Err(status) => return status,
    };
    if *gap_existence == 0 && *gap_extension == 0 && !values.non_affine.is_empty() {
        return 0;
    }
    let found = values.normal.iter().any(|row| {
        crate::math::blast_nint(row[0]) as i32 == *gap_existence
            && crate::math::blast_nint(row[1]) as i32 == *gap_extension
    });
    if !found && (*gap_existence < values.gap_open_max || *gap_extension < values.gap_extend_max) {
        *gap_existence = values.gap_open_max;
        *gap_extension = values.gap_extend_max;
    }
    0
}

/// NCBI: Blast_KarlinBlkNuclGappedCalc (blast_stat.c:3846).
pub fn blast_karlin_blk_nucl_gapped_calc(
    kbp: Option<&mut KarlinBlk>,
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    kbp_ungap: Option<&KarlinBlk>,
    round_down: Option<&mut bool>,
    error_return: Option<&mut Option<Box<crate::diagnostics::BlastMessage>>>,
) -> i16 {
    let values = match s_get_nucl_values_array(reward, penalty) {
        Ok(values) => values,
        Err(status) => return status,
    };
    if let Some(round_down) = round_down {
        *round_down = values.round_down;
    }
    let (Some(kbp), Some(kbp_ungap)) = (kbp, kbp_ungap) else {
        return 1;
    };

    if gap_open == 0 && gap_extend == 0 && !values.non_affine.is_empty() {
        let row = values.non_affine[0];
        kbp.lambda = row[2];
        kbp.k = row[3];
        kbp.log_k = kbp.k.ln();
        kbp.h = row[4];
        kbp.round_down = values.round_down;
        return 0;
    }

    for row in &values.normal {
        if crate::math::blast_nint(row[0]) as i32 == gap_open
            && crate::math::blast_nint(row[1]) as i32 == gap_extend
        {
            kbp.lambda = row[2];
            kbp.k = row[3];
            kbp.log_k = kbp.k.ln();
            kbp.h = row[4];
            kbp.round_down = values.round_down;
            return 0;
        }
    }

    if gap_open >= values.gap_open_max && gap_extend >= values.gap_extend_max {
        *kbp = kbp_ungap.clone();
        kbp.round_down = values.round_down;
        return 0;
    }

    if let Some(error_return) = error_return {
        let mut message = format!(
            "Gap existence and extension values {} and {} are not supported for substitution scores {} and {}\n",
            gap_open, gap_extend, reward, penalty
        );
        for row in &values.normal {
            message.push_str(&format!(
                "{} and {} are supported existence and extension values\n",
                crate::math::blast_nint(row[0]) as i32,
                crate::math::blast_nint(row[1]) as i32
            ));
        }
        message.push_str(&format!(
            "{} and {} are supported existence and extension values\n",
            values.gap_open_max, values.gap_extend_max
        ));
        message.push_str(&format!(
            "Any values more stringent than {} and {} are supported\n",
            values.gap_open_max, values.gap_extend_max
        ));
        let _ = crate::diagnostics::blast_message_write(
            error_return,
            crate::diagnostics::BlastSeverity::Error,
            crate::diagnostics::K_BLAST_MESSAGE_NO_CONTEXT,
            &message,
        );
        return 1;
    }

    0
}

/// NCBI: BLAST_CheckRewardPenaltyScores (blast_stat.c:3454).
pub fn blast_check_reward_penalty_scores(reward: i32, penalty: i32) -> bool {
    s_get_nucl_values_array(reward, penalty).is_ok()
}

/// NCBI: BLAST_PrintAllowedValues (blast_stat.c:3792).
/// naming: Returns a Rust `String` instead of writing to a C message buffer.
pub fn blast_print_allowed_values(matrix_name: &str, gap_open: i32, gap_extend: i32) -> String {
    let mut out = format!(
        "Gap existence and extension values of {} and {} not supported for {}\nsupported values are:\n",
        gap_open, gap_extend, matrix_name
    );
    if let Some(rows) = matrix_stat_rows(matrix_name) {
        for row in rows {
            out.push_str(&format!("{}, {}\n", row.gap_open, row.gap_extend));
        }
    }
    out
}

/// NCBI: BlastKarlinReportAllowedValues (blast_stat.c:3478).
/// naming: Returns Rust message strings instead of appending to a C message list.
pub fn blast_karlin_report_allowed_values(matrix_name: &str) -> Vec<String> {
    matrix_stat_rows(matrix_name)
        .map(|rows| {
            rows.iter()
                .map(|row| {
                    format!(
                        "Gap existence and extension values of {} and {} are supported",
                        row.gap_open, row.gap_extend
                    )
                })
                .collect()
        })
        .unwrap_or_default()
}

/// blast-rs: Nucleotide gapped-Karlin lookup helper with scaled-score support;
/// not a direct NCBI C port.
/// Looks up gapped KBP params for nucleotide, matching C's
/// Blast_KarlinBlkNuclGappedCalc table semantics.
/// Returns Ok((lambda, k, log_k, h, round_down)) or Err message.
pub fn scaled_nucl_gapped_kbp_lookup(
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    ungapped: &KarlinBlk,
) -> Result<(KarlinBlk, bool), String> {
    let divisor = crate::math::blast_gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    let meta = nucleotide_kbp_table_meta(nr, np)
        .ok_or_else(|| format!("Unsupported scores {} {}", reward, penalty))?;

    let round_down = meta.round_down;

    // Split: first row with gap_open==0 && gap_extend==0 is the linear entry
    let (affine, linear) =
        if !meta.table.is_empty() && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
            (&meta.table[1..], Some(&meta.table[0]))
        } else {
            (meta.table, None)
        };

    // GCD scaling (Rust extension beyond NCBI `Blast_KarlinBlkNuclGappedCalc`,
    // `blast_stat.c:3846`): NCBI's lookup reduces `(reward, penalty)` by
    // their GCD when selecting which table to use but then compares gap
    // costs against the reduced-system table values directly — users
    // must supply gap costs in the reduced units. Rust instead scales
    // the table entries by the divisor so users can pass gap costs in
    // the scaled system (e.g. `(reward=10, penalty=-20, gap_open=30,
    // gap_extend=10)` reduces to KBPT_1_2 with scaled lambda = row[2]/10).
    let (mut go_max, mut ge_max) = (meta.gap_open_max, meta.gap_extend_max);
    let scale = |row: &KbpTableRow, d: i32| -> (f64, f64, f64, f64, f64) {
        let go = row[0] * d as f64;
        let ge = row[1] * d as f64;
        let lam = row[2] / d as f64;
        (go, ge, lam, row[3], row[4])
    };
    if divisor != 1 {
        go_max *= divisor;
        ge_max *= divisor;
    }

    // Linear (non-affine) case
    if gap_open == 0 && gap_extend == 0 {
        if let Some(lin) = linear {
            let (_, _, lam, k, h) = scale(lin, divisor);
            return Ok((
                KarlinBlk {
                    lambda: lam,
                    k,
                    log_k: k.ln(),
                    h,
                    round_down,
                },
                round_down,
            ));
        }
    }

    // Search affine entries
    for row in affine {
        let (go, ge, lam, k, h) = scale(row, divisor);
        if go as i32 == gap_open && ge as i32 == gap_extend {
            return Ok((
                KarlinBlk {
                    lambda: lam,
                    k,
                    log_k: k.ln(),
                    h,
                    round_down,
                },
                round_down,
            ));
        }
    }

    // Fallback: gap costs exceed table maximum → use ungapped params
    if gap_open >= go_max && gap_extend >= ge_max {
        let mut kbp = ungapped.clone();
        kbp.round_down = round_down;
        return Ok((kbp, round_down));
    }

    Err(format!(
        "Unsupported gap costs {} {} for scores {} {}",
        gap_open, gap_extend, reward, penalty
    ))
}

/// blast-rs: Tuple-returning length-adjustment helper; not a direct NCBI C port.
/// Computes the length adjustment for effective search space.
/// Returns (length_adjustment, converged).
pub fn compute_length_adjustment_exact(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i32,
    db_length: i64,
    db_num_seqs: i32,
) -> (i32, bool) {
    let m = query_length as f64;
    let n = db_length as f64;
    let nn = db_num_seqs as f64;

    // Compute ell_max using quadratic formula
    let a = nn;
    let mb = m * nn + n;
    let c = n * m - m.max(n) / k;
    if c < 0.0 {
        return (0, false);
    }
    let ell_max_init = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).sqrt());

    let mut ell_min = 0.0_f64;
    let mut ell_max = ell_max_init;
    let mut ell_next = 0.0_f64;
    let mut converged = false;

    for i in 1..=20 {
        let ell = ell_next;
        let ss = (m - ell) * (n - nn * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;

        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }

        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 {
                ell_max
            } else {
                (ell_min + ell_max) / 2.0
            };
        }
    }

    let mut adj = ell_min as i32;
    if converged {
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - nn * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell {
                adj = ell as i32;
            }
        }
    }
    (adj, converged)
}

/// NCBI: BLAST_ComputeLengthAdjustment (blast_stat.c:5041).
/// Port with pointer-style output and integer convergence status.
pub fn blast_compute_length_adjustment(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i32,
    db_length: i64,
    db_num_seqs: i32,
    length_adjustment: Option<&mut i32>,
) -> i32 {
    let Some(length_adjustment) = length_adjustment else {
        return 1;
    };
    let (adj, converged) = compute_length_adjustment_exact(
        k,
        log_k,
        alpha_d_lambda,
        beta,
        query_length,
        db_length,
        db_num_seqs,
    );
    *length_adjustment = adj;
    if converged {
        0
    } else {
        1
    }
}

/// blast-rs: Tuple-returning nucleotide alpha/beta lookup helper with
/// scaled-score support; not a direct NCBI C port.
pub fn scaled_nucl_alpha_beta(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_lambda: f64,
    ungapped_h: f64,
    gapped: bool,
) -> (f64, f64) {
    if !gapped {
        return (
            ungapped_lambda / ungapped_h,
            s_get_ungapped_beta(reward, penalty),
        );
    }

    // Rust extension beyond NCBI `Blast_GetNuclAlphaBeta`
    // (`blast_stat.c:3965`): scale gap costs and alpha by the GCD
    // divisor so scaled scoring systems like `(10, -20)` match the
    // reduced `(1, -2)` table. NCBI compares against raw table values
    // and forces users to supply reduced-system gap costs.
    let divisor = crate::math::blast_gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    if let Some(meta) = nucleotide_kbp_table_meta(nr, np) {
        let (affine, linear) =
            if !meta.table.is_empty() && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
                (&meta.table[1..], Some(&meta.table[0]))
            } else {
                (meta.table, None)
            };

        if gap_open == 0 && gap_extend == 0 {
            if let Some(lin) = linear {
                let a = lin[5] / divisor as f64;
                let b = lin[6];
                return (a, b);
            }
        }

        for row in affine {
            let go = (row[0] * divisor as f64) as i32;
            let ge = (row[1] * divisor as f64) as i32;
            if go == gap_open && ge == gap_extend {
                let a = row[5] / divisor as f64;
                let b = row[6];
                return (a, b);
            }
        }
    }

    // Fallback: ungapped values. NCBI `Blast_GetNuclAlphaBeta`
    // (`blast_stat.c:4022-4026`): on lookup miss, returns
    // `alpha = Lambda/H` and `beta = s_GetUngappedBeta(reward, penalty)`.
    (
        ungapped_lambda / ungapped_h,
        s_get_ungapped_beta(reward, penalty),
    )
}

/// NCBI: Blast_GetNuclAlphaBeta (blast_stat.c:3965).
pub fn blast_get_nucl_alpha_beta(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    lambda: f64,
    h: f64,
    gapped: bool,
    alpha: &mut f64,
    beta: &mut f64,
) -> i16 {
    let values = match s_get_nucl_values_array(reward, penalty) {
        Ok(values) => values,
        Err(status) => return status,
    };
    let mut found = false;
    if gapped {
        if gap_open == 0 && gap_extend == 0 && !values.non_affine.is_empty() {
            let row = values.non_affine[0];
            *alpha = row[5];
            *beta = row[6];
            found = true;
        } else {
            for row in &values.normal {
                if crate::math::blast_nint(row[0]) as i32 == gap_open
                    && crate::math::blast_nint(row[1]) as i32 == gap_extend
                {
                    *alpha = row[5];
                    *beta = row[6];
                    found = true;
                    break;
                }
            }
        }
    }
    if !found {
        *alpha = lambda / h;
        *beta = s_get_ungapped_beta(reward, penalty);
    }
    0
}

/// NCBI: Blast_FillResidueProbability (blast_stat.c:4577).
pub fn blast_fill_residue_probability(sequence: &[u8], res_prob: &mut [f64]) {
    let mut frequency = vec![0i32; crate::encoding::BLASTAA_SIZE];
    let mut denominator = sequence.len() as i32;
    for &residue in sequence {
        if residue == crate::encoding::NCBISTDAA_X {
            denominator -= 1;
        } else if (residue as usize) < frequency.len() {
            frequency[residue as usize] += 1;
        }
    }
    for index in 0..res_prob.len().min(frequency.len()) {
        res_prob[index] = if frequency[index] == 0 || denominator <= 0 {
            0.0
        } else {
            frequency[index] as f64 / denominator as f64
        };
    }
}

/// NCBI: RPSfindUngappedLambda (blast_stat.c:4610).
/// naming: Rust splits NCBI's `RPSfind` mixed-case prefix as `rps_find`.
pub fn rps_find_ungapped_lambda(matrix_name: Option<&str>) -> f64 {
    let values = blast_get_matrix_values(matrix_name);
    values.lambda.first().copied().unwrap_or(0.0)
}

/// NCBI: RPSFillScores (blast_stat.c:4649).
pub fn rps_fill_scores(
    matrix: &[Vec<i32>],
    query_prob_array: &[f64],
    return_sfp: &mut ScoreFreq,
    alphabet_size: usize,
) {
    let mut min_score = 0;
    let mut max_score = 0;
    for row in matrix {
        for j in 0..alphabet_size.min(row.len()) {
            if j == crate::encoding::NCBISTDAA_X as usize {
                continue;
            }
            let score = row[j];
            if score > BLAST_SCORE_MIN && score < min_score {
                min_score = score;
            }
            if score > max_score {
                max_score = score;
            }
        }
    }

    return_sfp.obs_min = min_score;
    return_sfp.obs_max = max_score;
    return_sfp.score_min = min_score;
    return_sfp.score_max = max_score;
    return_sfp.sprob = vec![0.0; (max_score - min_score + 1).max(0) as usize];
    if matrix.is_empty() {
        return;
    }
    let recip_length = 1.0 / matrix.len() as f64;
    for row in matrix {
        for j in 0..alphabet_size.min(row.len()).min(query_prob_array.len()) {
            if j == crate::encoding::NCBISTDAA_X as usize {
                continue;
            }
            let score = row[j];
            if score >= min_score {
                let index = (score - min_score) as usize;
                if let Some(slot) = return_sfp.sprob.get_mut(index) {
                    *slot += recip_length * query_prob_array[j];
                }
            }
        }
    }
    return_sfp.score_avg = (min_score..=max_score)
        .map(|score| {
            let index = (score - min_score) as usize;
            score as f64 * return_sfp.sprob.get(index).copied().unwrap_or(0.0)
        })
        .sum();
}

/// blast-rs: Internal `ScoreFreq` adapter for the public lambda-NR wrapper;
/// not a direct NCBI C port.
fn blast_karlin_lambda_nr_from_score_freq(sfp: &ScoreFreq, initial_lambda: f64) -> f64 {
    let dist = sf_dist_from_score_freq(sfp);
    if dist.score_avg >= 0.0 {
        return -1.0;
    }
    let low = dist.obs_min;
    let high = dist.obs_max;
    if blast_score_chk(low, high) != 0 {
        return -1.0;
    }
    let mut d = -low;
    for i in 1..=(high - low) {
        if d <= 1 {
            break;
        }
        if dist.p(low + i) != 0.0 {
            d = crate::math::blast_gcd(d, i);
        }
    }
    solve_lambda(&dist, d, low, high, initial_lambda)
}

/// NCBI: Blast_KarlinLambdaNR (blast_stat.c:2567).
/// Port with an explicit initial lambda.
pub fn blast_karlin_lambda_nr(sfp: Option<&ScoreFreq>, initial_lambda: f64) -> f64 {
    let Some(sfp) = sfp else {
        return -1.0;
    };
    blast_karlin_lambda_nr_from_score_freq(sfp, initial_lambda)
}

/// NCBI: RPSRescalePssm (blast_stat.c:4693).
/// naming: Uses owned Rust matrix rows instead of `_PSIAllocateMatrix`.
pub fn rps_rescale_pssm(
    scaling_factor: f64,
    rps_query_length: i32,
    rps_query_seq: Option<&[u8]>,
    db_seq_length: i32,
    pos_matrix: Option<&[Vec<i32>]>,
    sbp: Option<&BlastScoreBlk>,
) -> Option<Vec<Vec<i32>>> {
    let (Some(rps_query_seq), Some(pos_matrix), Some(sbp)) = (rps_query_seq, pos_matrix, sbp)
    else {
        return None;
    };
    if scaling_factor == 0.0 || rps_query_length < 0 || db_seq_length < 0 {
        return None;
    }
    let db_seq_length = db_seq_length as usize;
    if pos_matrix.len() < db_seq_length {
        return None;
    }
    let query_len = (rps_query_length as usize).min(rps_query_seq.len());
    let mut res_prob = vec![0.0; crate::encoding::BLASTAA_SIZE];
    blast_fill_residue_probability(&rps_query_seq[..query_len], &mut res_prob);

    let alphabet_size = sbp.alphabet_size.min(crate::encoding::BLASTAA_SIZE);
    if alphabet_size == 0
        || pos_matrix
            .iter()
            .take(db_seq_length)
            .any(|row| row.len() < alphabet_size)
    {
        return None;
    }
    let mut return_sfp = ScoreFreq {
        score_min: 0,
        score_max: 0,
        obs_min: 0,
        obs_max: 0,
        score_avg: 0.0,
        sprob: Vec::new(),
    };
    rps_fill_scores(
        &pos_matrix[..db_seq_length],
        &res_prob,
        &mut return_sfp,
        alphabet_size,
    );

    let initial_ungapped_lambda = rps_find_ungapped_lambda(sbp.name.as_deref());
    if initial_ungapped_lambda <= 0.0 {
        return None;
    }
    let scaled_initial_ungapped_lambda = initial_ungapped_lambda / scaling_factor;
    let correct_ungapped_lambda =
        blast_karlin_lambda_nr(Some(&return_sfp), scaled_initial_ungapped_lambda);
    if correct_ungapped_lambda == -1.0 {
        return None;
    }
    let final_lambda = correct_ungapped_lambda / scaled_initial_ungapped_lambda;

    let mut return_matrix =
        vec![vec![BLAST_SCORE_MIN; crate::encoding::BLASTAA_SIZE]; db_seq_length];
    for index in 0..db_seq_length {
        for inner_index in 0..alphabet_size {
            let score = pos_matrix[index][inner_index];
            return_matrix[index][inner_index] = if score <= BLAST_SCORE_MIN
                || inner_index == crate::encoding::NCBISTDAA_X as usize
            {
                score
            } else {
                crate::math::blast_nint(score as f64 * final_lambda) as i32
            };
        }
    }
    Some(return_matrix)
}

pub const COMPRESSED_REVERSE_LOOKUP_SIZE: usize = crate::encoding::BLASTAA_SIZE + 1;
pub type CompressedReverseLookup =
    [[i8; COMPRESSED_REVERSE_LOOKUP_SIZE]; COMPRESSED_REVERSE_LOOKUP_SIZE];

pub const S_ALPHABET10: &str = "IJLMV AST BDENZ KQR G FY P H C W";
pub const S_ALPHABET15: &str = "ST IJV LM KR EQZ A G BD P N F Y H C W";

/// Rust-owned counterpart of NCBI `SCompressedAlphabet`.
#[derive(Debug, Clone, Default)]
pub struct SCompressedAlphabet {
    pub compressed_alphabet_size: i32,
    pub compress_table: Vec<u8>,
    pub matrix: Option<BlastScoreMatrix>,
}

/// NCBI: s_BuildCompressedTranslation (blast_stat.c:4736).
pub fn s_build_compressed_translation(
    trans_string: &str,
    table: &mut [u8],
    compressed_alphabet_size: i32,
    rev_table: &mut CompressedReverseLookup,
) {
    for value in table.iter_mut().take(crate::encoding::BLASTAA_SIZE) {
        *value = compressed_alphabet_size as u8;
    }
    for row in rev_table.iter_mut() {
        row.fill(-1);
    }

    let mut compressed_letter = 0usize;
    let mut j = 0usize;
    for byte in trans_string.bytes() {
        if byte.is_ascii_whitespace() {
            compressed_letter += 1;
            j = 0;
        } else if byte.is_ascii_alphabetic() {
            let aa_letter = crate::encoding::AMINOACID_TO_NCBISTDAA[byte as usize] as usize;
            if aa_letter < table.len() {
                table[aa_letter] = compressed_letter as u8;
            }
            if compressed_letter < rev_table.len() && j + 1 < rev_table[compressed_letter].len() {
                rev_table[compressed_letter][j] = aa_letter as i8;
                j += 1;
                rev_table[compressed_letter][j] = -1;
            }
        }
    }
}

/// NCBI: s_GetCompressedProbs (blast_stat.c:4772).
pub fn s_get_compressed_probs(
    sbp: Option<&BlastScoreBlk>,
    compressed_prob: &mut [f64],
    compressed_alphabet_size: i32,
    rev_table: &CompressedReverseLookup,
) -> i16 {
    let Some(sbp) = sbp else {
        return -1;
    };
    let mut rfp = match blast_res_freq_new(Some(sbp)) {
        Some(rfp) => rfp,
        None => return -1,
    };
    if blast_res_freq_std_comp(Some(sbp), Some(&mut rfp)) != 0 {
        return -1;
    }

    for value in compressed_prob
        .iter_mut()
        .take(crate::encoding::BLASTAA_SIZE)
    {
        *value = 0.0;
    }

    for letter in 0..compressed_alphabet_size.max(0) as usize {
        let mut prob_sum = 0.0;
        for &aa in &rev_table[letter] {
            if aa < 0 {
                break;
            }
            prob_sum += rfp.prob.get(aa as usize).copied().unwrap_or(0.0);
        }
        for &aa in &rev_table[letter] {
            if aa < 0 {
                break;
            }
            if let Some(slot) = compressed_prob.get_mut(aa as usize) {
                *slot = if prob_sum == 0.0 {
                    0.0
                } else {
                    rfp.prob[aa as usize] / prob_sum
                };
            }
        }
    }
    0
}

/// NCBI: s_BuildCompressedScoreMatrix (blast_stat.c:4818).
pub fn s_build_compressed_score_matrix(
    sbp: Option<&BlastScoreBlk>,
    new_alphabet: Option<&mut SCompressedAlphabet>,
    mut matrix_scale_factor: f64,
    rev_table: &CompressedReverseLookup,
) -> i16 {
    let (Some(sbp), Some(new_alphabet)) = (sbp, new_alphabet) else {
        return -1;
    };
    let lambda = rps_find_ungapped_lambda(sbp.name.as_deref());
    if lambda <= 0.0 {
        return -1;
    }
    matrix_scale_factor /= lambda;

    let Some(std_freqs) =
        crate::matrix::get_matrix_freq_ratios(sbp.name.as_deref().unwrap_or_default())
    else {
        return -2;
    };

    let compressed_alphabet_size = new_alphabet.compressed_alphabet_size;
    let mut compressed_prob = vec![0.0; crate::encoding::BLASTAA_SIZE];
    if s_get_compressed_probs(
        Some(sbp),
        &mut compressed_prob,
        compressed_alphabet_size,
        rev_table,
    ) < 0
    {
        return -3;
    }

    let mut new_matrix = BlastScoreMatrix::new(
        crate::encoding::BLASTAA_SIZE,
        compressed_alphabet_size.max(0) as usize,
    );
    let min_freq = BLAST_SCORE_MIN as f64 / matrix_scale_factor;
    for q in 0..crate::encoding::BLASTAA_SIZE {
        for s in 0..compressed_alphabet_size.max(0) as usize {
            let mut val = 0.0;
            for &aa in &rev_table[s] {
                if aa < 0 {
                    break;
                }
                let aa = aa as usize;
                val += std_freqs[q][aa] * compressed_prob[aa];
            }
            val = if val < 1e-8 { min_freq } else { val.ln() };
            new_matrix.data[q][s] = crate::math::blast_nint(val * matrix_scale_factor) as i32;
        }
    }
    new_alphabet.matrix = Some(new_matrix);
    0
}

/// Port of NCBI `SCompressedAlphabetNew` (`blast_stat.c:4887`).
pub fn s_compressed_alphabet_new(
    sbp: Option<&BlastScoreBlk>,
    compressed_alphabet_size: i32,
    matrix_scale_factor: f64,
) -> Option<SCompressedAlphabet> {
    if compressed_alphabet_size != 10 && compressed_alphabet_size != 15 {
        return None;
    }
    let alphabet_string = if compressed_alphabet_size == 10 {
        S_ALPHABET10
    } else {
        S_ALPHABET15
    };
    let mut rev_table = [[-1; COMPRESSED_REVERSE_LOOKUP_SIZE]; COMPRESSED_REVERSE_LOOKUP_SIZE];
    let mut new_alphabet = SCompressedAlphabet {
        compressed_alphabet_size,
        compress_table: vec![0; crate::encoding::BLASTAA_SIZE],
        matrix: None,
    };
    s_build_compressed_translation(
        alphabet_string,
        &mut new_alphabet.compress_table,
        compressed_alphabet_size,
        &mut rev_table,
    );
    if s_build_compressed_score_matrix(
        sbp,
        Some(&mut new_alphabet),
        matrix_scale_factor,
        &rev_table,
    ) < 0
    {
        return None;
    }
    Some(new_alphabet)
}

/// Port of NCBI `SCompressedAlphabetFree` (`blast_stat.c:4916`).
pub fn s_compressed_alphabet_free(
    alphabet: &mut Option<SCompressedAlphabet>,
) -> Option<SCompressedAlphabet> {
    if let Some(alphabet) = alphabet.as_mut() {
        alphabet.compress_table.clear();
        if let Some(matrix) = alphabet.matrix.as_mut() {
            matrix.data.clear();
            matrix.nrows = 0;
            matrix.ncols = 0;
        }
        alphabet.matrix = None;
    }
    *alphabet = None;
    None
}

/// NCBI: s_GetUngappedBeta (blast_stat.c:3935).
/// `(1,-1)` and
/// `(2,-3)` scoring systems use `beta = -2`; every other combination uses
/// `beta = 0`.
pub fn s_get_ungapped_beta(reward: i32, penalty: i32) -> f64 {
    if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) {
        -2.0
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Ungapped Karlin-Altschul parameter computation (exact C-compatible)
// Port of Blast_ScoreBlkKbpUngappedCalc and sub-functions from blast_stat.c
// ---------------------------------------------------------------------------

/// Score frequency distribution (internal to KBP computation). Rust
/// analog of NCBI `Blast_ScoreFreq` (`blast_stat.h:128`); `probs`
/// corresponds to NCBI's `sprob` array (shifted by `score_min`).
struct SfDist {
    score_min: i32,
    score_max: i32,
    obs_min: i32,
    obs_max: i32,
    score_avg: f64,
    probs: Vec<f64>, // indexed by [score - score_min]
}

impl SfDist {
    /// blast-rs: Owned score-frequency distribution constructor; not a direct
    /// NCBI C port.
    fn new(lo: i32, hi: i32) -> Self {
        let n = (hi - lo + 1) as usize;
        SfDist {
            score_min: lo,
            score_max: hi,
            obs_min: 0,
            obs_max: 0,
            score_avg: 0.0,
            probs: vec![0.0; n],
        }
    }
    #[inline]
    /// blast-rs: Score-index accessor for internal Karlin probabilities; not a
    /// direct NCBI C port.
    fn p(&self, s: i32) -> f64 {
        if s < self.score_min || s > self.score_max {
            0.0
        } else {
            self.probs[(s - self.score_min) as usize]
        }
    }
    #[inline]
    /// blast-rs: Mutable score-index accessor for internal Karlin probabilities;
    /// not a direct NCBI C port.
    fn p_mut(&mut self, s: i32) -> &mut f64 {
        &mut self.probs[(s - self.score_min) as usize]
    }
}

/// Newton-Raphson solver for Lambda in x = exp(-lambda) space.
/// Port of NlmKarlinLambdaNR from blast_stat.c.
/// NCBI `BLAST_KARLIN_K_SUMLIMIT_DEFAULT` (`blast_stat.c:64`): inner-sum
/// convergence threshold in the K solver.
const BLAST_KARLIN_K_SUMLIMIT_DEFAULT: f64 = 0.0001;
/// NCBI `BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT` (`blast_stat.c:66`):
/// accuracy target for lambda solve.
const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
/// NCBI `BLAST_KARLIN_LAMBDA_ITER_DEFAULT` (`blast_stat.c:68`): number
/// of bisection iterations after the Newton phase (20 + 17 = 37 total).
const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
/// NCBI `BLAST_KARLIN_K_ITER_MAX` (`blast_stat.c:72`): upper limit on
/// iterations for the K solver.
const BLAST_KARLIN_K_ITER_MAX: usize = 100;

/// blast-rs: Internal Karlin lambda Newton/bisection solver; not a direct NCBI
/// C port. Solves the Karlin-Altschul lambda parameter
/// over `x = exp(-lambda)`.
fn solve_lambda(sfp: &SfDist, d: i32, low: i32, high: i32, lambda0: f64) -> f64 {
    let x0 = (-lambda0).exp();
    let mut x = if x0 > 0.0 && x0 < 1.0 { x0 } else { 0.5 };
    let mut a = 0.0_f64;
    let mut b = 1.0_f64;
    let mut f = 4.0_f64;
    let mut is_newton = false;
    // NCBI's `Blast_KarlinLambdaNR` calls `NlmKarlinLambdaNR(... 20, 20 +
    // BLAST_KARLIN_LAMBDA_ITER_DEFAULT, ...)` — passing `itmax=20` and
    // `maxNewton=37`. In NCBI's implementation `itmax` is the loop bound and
    // `maxNewton` is the Newton-fallback threshold; with maxNewton > itmax
    // the threshold never fires within the loop. So the effective behavior
    // is "20 iterations, Newton always tried." Match NCBI exactly.
    let max_iter = 20;
    let max_newton = 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT; // 37, threshold beyond max_iter
    let tolx = BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT;

    for k in 0..max_iter {
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;

        // Horner evaluation of polynomial and derivative.
        let mut g = 0.0_f64;
        f = sfp.p(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + sfp.p(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + sfp.p(0) - 1.0;
        i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + sfp.p(i);
            i += d;
        }

        if f > 0.0 {
            a = x;
        } else if f < 0.0 {
            b = x;
        } else {
            break;
        }
        if b - a < 2.0 * a * (1.0 - b) * tolx {
            x = (a + b) / 2.0;
            break;
        }

        if k >= max_newton || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break;
                }
            }
        }
    }
    -x.ln() / d as f64
}

/// NCBI `BLAST_KARLIN_LAMBDA0_DEFAULT` (`blast_stat.c:70`): initial
/// guess for the Newton-Raphson lambda solver.
const BLAST_KARLIN_LAMBDA0_DEFAULT: f64 = 0.5;

/// blast-rs: Internal lambda computation over `SfDist`; not a direct NCBI C port.
/// naming: Short Rust helper name is scoped to the internal Karlin solver; the
/// public compatibility wrapper is `blast_karlin_lambda_nr`.
/// Uses NCBI `Blast_KarlinLambdaNR` (`blast_stat.c:2567`) score-GCD reduction
/// followed by Newton-Raphson solve on `x = exp(-lambda)`.
fn compute_lambda(sfp: &SfDist) -> f64 {
    if sfp.score_avg >= 0.0 {
        return -1.0;
    }
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    // NCBI `Blast_KarlinLambdaNR` (`blast_stat.c:2581`) gates on
    // `BlastScoreChk(low, high)` before the GCD loop and returns -1 if
    // the score range is invalid (lo>=0, hi<=0, out-of-range, or span
    // exceeding `BLAST_SCORE_RANGE_MAX`).
    if blast_score_chk(low, high) != 0 {
        return -1.0;
    }
    let mut d = -low;
    for i in 1..=(high - low) {
        if d <= 1 {
            break;
        }
        if sfp.p(low + i) != 0.0 {
            d = crate::math::blast_gcd(d, i);
        }
    }
    solve_lambda(sfp, d, low, high, BLAST_KARLIN_LAMBDA0_DEFAULT)
}

/// NCBI: BlastKarlinLtoH (blast_stat.c:2607).
/// Computes H (relative entropy) from score frequencies and Lambda. NCBI gates on
/// `BlastScoreChk(low, high)` before the formula and returns -1 if the
/// score range is invalid.
fn blast_karlin_lto_h(sfp: &SfDist, lambda: f64) -> f64 {
    if lambda < 0.0 {
        return -1.0;
    }
    if blast_score_chk(sfp.obs_min, sfp.obs_max) != 0 {
        return -1.0;
    }
    let etonlam = (-lambda).exp();
    let mut sum = sfp.obs_min as f64 * sfp.p(sfp.obs_min);
    for s in (sfp.obs_min + 1)..=sfp.obs_max {
        sum = s as f64 * sfp.p(s) + etonlam * sum;
    }
    let scale = crate::math::blast_powi(etonlam, sfp.obs_max);
    if scale > 0.0 {
        lambda * sum / scale
    } else {
        lambda * (lambda * sfp.obs_max as f64 + sum.ln()).exp()
    }
}

/// NCBI: BlastKarlinLHtoK (blast_stat.c:2346).
/// naming: Rust separates the `LHtoK` acronym group as `lh_to_k`.
/// Computes K from Lambda and H using the C-compatible DP algorithm.
fn blast_karlin_lh_to_k(sfp: &SfDist, lambda: f64, h: f64) -> f64 {
    if lambda <= 0.0 || h <= 0.0 || sfp.score_avg >= 0.0 {
        return -1.0;
    }
    let olow = sfp.obs_min;
    let ohigh = sfp.obs_max;
    let mut divisor = -olow;
    for i in 1..=(ohigh - olow) {
        if divisor <= 1 {
            break;
        }
        if sfp.p(olow + i) != 0.0 {
            divisor = crate::math::blast_gcd(divisor, i);
        }
    }

    let high = ohigh / divisor;
    let low = olow / divisor;
    let lam_s = lambda * divisor as f64;
    let range = high - low;
    let ftcf = h / lam_s;
    let eml = (-lam_s).exp();

    // Special cases
    if low == -1 && high == 1 {
        let pl = sfp.p(olow);
        let ph = sfp.p(ohigh);
        return (pl - ph) * (pl - ph) / pl;
    }
    if low == -1 || high == 1 {
        let mut f = ftcf;
        if high != 1 {
            let sa = sfp.score_avg / divisor as f64;
            f = (sa * sa) / f;
        }
        return f * (1.0 - eml);
    }

    // DP
    let prob_at = |i: i32| -> f64 {
        let s = olow + i;
        if s >= sfp.score_min && s <= sfp.score_max {
            sfp.p(s)
        } else {
            0.0
        }
    };
    let ru = range as usize;
    let max_iter = BLAST_KARLIN_K_ITER_MAX;
    let mut asp = vec![0.0_f64; max_iter * ru + 1];
    asp[0] = 1.0;
    let mut outer_sum = 0.0_f64;
    let mut low_as = 0i32;
    let mut high_as = 0i32;
    let mut inner_sum = 1.0_f64;
    let mut oldsum;

    // NCBI's `BlastKarlinLHtoK` (`blast_stat.c:2346`) uses an idiomatic C
    // for-loop:
    //   for (iterCounter = 0;
    //        (iterCounter < iterlimit) && (innerSum > sumlimit);
    //        outerSum += innerSum /= ++iterCounter)
    // The increment expression DIVIDES innerSum in-place by the new
    // iterCounter BEFORE adding to outerSum. So the loop CONDITION at
    // iter K (K > 0) compares `innerSum_{K-1} / K > sumlimit` — i.e. the
    // divided value carries into the next iter's check. That makes the
    // effective sumlimit grow as K (sumlimit * K threshold). We mirror
    // this by mutating `inner_sum` after the body.
    let mut iter_counter = 0i32;
    while iter_counter < max_iter as i32 && inner_sum > BLAST_KARLIN_K_SUMLIMIT_DEFAULT {
        let mut first = range;
        let mut last = range;
        low_as += low;
        high_as += high;
        let span = (high_as - low_as) as isize;
        let mut pp = span;
        while pp >= 0 {
            let p1s = pp - first as isize;
            let p1e = pp - last as isize;
            let mut isum = 0.0;
            let mut p1 = p1s;
            let mut p2 = first;
            while p1 >= p1e {
                if p1 >= 0 {
                    isum += asp[p1 as usize] * prob_at(p2);
                }
                p1 -= 1;
                p2 += 1;
            }
            asp[pp as usize] = isum;
            if first > 0 {
                first -= 1;
            }
            if pp <= range as isize {
                last -= 1;
            }
            pp -= 1;
        }
        pp += 1;
        inner_sum = asp[pp as usize];
        let mut i = low_as + 1;
        while i < 0 {
            pp += 1;
            inner_sum = asp[pp as usize] + inner_sum * eml;
            i += 1;
        }
        inner_sum *= eml;
        while i <= high_as {
            pp += 1;
            inner_sum += asp[pp as usize];
            i += 1;
        }
        oldsum = inner_sum;
        let _ = oldsum;
        // Match NCBI's `outerSum += innerSum /= ++iterCounter`:
        // first increment iterCounter, then divide innerSum in place,
        // then add to outerSum.
        iter_counter += 1;
        inner_sum /= iter_counter as f64;
        outer_sum += inner_sum;
    }

    -(-2.0 * outer_sum).exp() / (ftcf * crate::math::expm1(-lam_s))
}

/// Context info for ungapped KBP computation.
pub struct UngappedKbpContext {
    pub query_offset: i32,
    pub query_length: i32,
    pub is_valid: bool,
}

/// blast-rs: Shared owned-data implementation for ungapped KBP contexts; not a
/// direct NCBI C port.
/// Computes ungapped KBP for all contexts. Returns per-context `Option<KarlinBlk>`.
/// Implements the Blast_ScoreBlkKbpUngappedCalc nucleotide path —
/// std composition is fixed at 0.25 for A/C/G/T).
pub fn ungapped_kbp_calc(
    query: &[u8],
    contexts: &[UngappedKbpContext],
    loscore: i32,
    hiscore: i32,
    alphabet_size: usize,
    ambiguous: &[u8],
    matrix: &dyn Fn(usize, usize) -> i32,
) -> Vec<Option<KarlinBlk>> {
    // Standard composition: 0.25 each for A/C/G/T (indices 0-3)
    let mut std_freq = vec![0.0f64; alphabet_size];
    for i in 0..4.min(alphabet_size) {
        std_freq[i] = 0.25;
    }
    ungapped_kbp_calc_with_std(
        query,
        contexts,
        loscore,
        hiscore,
        alphabet_size,
        ambiguous,
        &std_freq,
        matrix,
    )
}

/// Like [`ungapped_kbp_calc`] but takes an explicit `std_freq` array.
/// Used by protein paths (blastx/tblastn/tblastx) where the standard
/// background is BLOSUM62's Robinson & Robinson amino acid frequencies
/// rather than the uniform 0.25 nucleotide composition.
///
/// Mirrors `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2737`): per-context
/// `Blast_ResFreqString` → `BlastScoreFreqCalc` → `Blast_KarlinBlkUngappedCalc`.
/// blast-rs: Shared ungapped KBP implementation with explicit standard
/// composition; not a direct NCBI C port.
pub fn ungapped_kbp_calc_with_std(
    query: &[u8],
    contexts: &[UngappedKbpContext],
    loscore: i32,
    hiscore: i32,
    alphabet_size: usize,
    ambiguous: &[u8],
    std_freq: &[f64],
    matrix: &dyn Fn(usize, usize) -> i32,
) -> Vec<Option<KarlinBlk>> {
    let std_freq = std_freq.to_vec();
    let mut results = Vec::with_capacity(contexts.len());
    for ctx in contexts {
        if !ctx.is_valid || ctx.query_length <= 0 {
            results.push(None);
            continue;
        }
        let off = ctx.query_offset as usize;
        let len = ctx.query_length as usize;
        let buf = &query[off..off + len];

        // Count residue composition. The 4-bit mask was needed only when the
        // caller was nucleotide-only (BLASTNA bytes are 0..15); now the
        // function is also used for protein (NCBIstdaa, bytes 0..27).
        // Removing the mask is a no-op for nucleotide (bytes already <16) and
        // unblocks protein. The `idx < alphabet_size` guard still bounds.
        let mut counts = vec![0i32; alphabet_size];
        for &b in buf {
            let idx = b as usize;
            if idx < alphabet_size {
                counts[idx] += 1;
            }
        }
        for &a in ambiguous {
            if (a as usize) < alphabet_size {
                counts[a as usize] = 0;
            }
        }
        let sum: f64 = counts.iter().map(|&c| c as f64).sum();
        let mut qfreq = vec![0.0f64; alphabet_size];
        if sum > 0.0 {
            for i in 0..alphabet_size {
                qfreq[i] = counts[i] as f64 / sum;
            }
        }

        // Compute score frequency
        let mut sfp = SfDist::new(loscore, hiscore);
        for i in 0..alphabet_size {
            for j in 0..alphabet_size {
                let s = matrix(i, j);
                if s >= loscore && s <= hiscore {
                    *sfp.p_mut(s) += qfreq[i] * std_freq[j];
                }
            }
        }
        // Find obs range and normalize
        let mut obs_min = i32::MAX;
        let mut obs_max = i32::MIN;
        let mut psum = 0.0;
        for s in loscore..=hiscore {
            if sfp.p(s) > 0.0 {
                psum += sfp.p(s);
                if s < obs_min {
                    obs_min = s;
                }
                obs_max = s;
            }
        }
        sfp.obs_min = obs_min;
        sfp.obs_max = obs_max;
        let mut savg = 0.0;
        if psum > 0.0 {
            for s in obs_min..=obs_max {
                *sfp.p_mut(s) /= psum;
                savg += s as f64 * sfp.p(s);
            }
        }
        sfp.score_avg = savg;

        // Compute Lambda, H, K
        let lambda = compute_lambda(&sfp);
        if lambda < 0.0 {
            results.push(None);
            continue;
        }
        let h = blast_karlin_lto_h(&sfp, lambda);
        if h < 0.0 {
            results.push(None);
            continue;
        }
        let k = blast_karlin_lh_to_k(&sfp, lambda, h);
        if k < 0.0 {
            results.push(None);
            continue;
        }
        results.push(Some(KarlinBlk {
            lambda,
            k,
            log_k: k.ln(),
            h,
            round_down: false,
        }));
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blast_score_blk_new_free_and_nucl_matrix_create_match_core_fields() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 2).expect("score block");
        assert_eq!(sbp.alphabet_size, crate::encoding::BLASTNA_SIZE);
        assert!(!sbp.protein_alphabet);
        assert_eq!(sbp.number_of_contexts, 2);
        assert_eq!(sbp.kbp.len(), 2);
        assert_eq!(sbp.kbp_std.len(), 2);
        assert_eq!(sbp.scale_factor, 1.0);
        assert_eq!(blast_score_blk_check(None), -1);
        assert_eq!(blast_score_blk_check(Some(&sbp)), 1);
        sbp.kbp[0] = KarlinBlk {
            lambda: 0.5,
            k: 0.2,
            log_k: 0.2_f64.ln(),
            h: 0.4,
            round_down: false,
        };
        assert_eq!(blast_score_blk_check(Some(&sbp)), 0);

        sbp.reward = 1;
        sbp.penalty = -3;
        assert_eq!(blast_score_blk_nucl_matrix_create(&mut sbp), 0);
        assert_eq!(sbp.matrix.data[0][0], 1);
        assert_eq!(sbp.matrix.data[0][1], -3);
        assert_eq!(
            sbp.matrix.data[crate::encoding::BLASTNA_SIZE - 1][0],
            i32::MIN / 2
        );
        assert_eq!(sbp.hiscore, 1);
        // NCBI `BlastScoreBlkMaxScoreSet` (`blast_stat.c:1513`) skips
        // `BLAST_SCORE_MIN` sentinels (gap rows) when computing loscore,
        // so the result is the lowest REAL score (-3 = penalty), not the
        // sentinel.
        assert_eq!(sbp.loscore, -3);

        let mut owned = Some(sbp);
        assert!(blast_score_blk_free(&mut owned).is_none());
        assert!(owned.is_none());
    }

    #[test]
    fn blast_score_blk_max_score_set_ignores_extra_matrix_cells() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 1).expect("score block");
        sbp.alphabet_size = 2;
        sbp.matrix = BlastScoreMatrix {
            nrows: 3,
            ncols: 3,
            data: vec![vec![5, -4, 99], vec![-3, 7, -88], vec![123, -99, 1000]],
        };

        assert_eq!(blast_score_blk_max_score_set(&mut sbp), 0);
        assert_eq!(sbp.loscore, -4);
        assert_eq!(sbp.hiscore, 7);
    }

    #[test]
    fn blast_score_blk_max_score_set_clamps_all_sentinel_matrix() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 1).expect("score block");
        sbp.alphabet_size = 2;
        sbp.matrix = BlastScoreMatrix {
            nrows: 2,
            ncols: 2,
            data: vec![
                vec![BLAST_SCORE_MIN, BLAST_SCORE_MAX],
                vec![BLAST_SCORE_MAX, BLAST_SCORE_MIN],
            ],
        };

        assert_eq!(blast_score_blk_max_score_set(&mut sbp), 0);
        assert_eq!(sbp.loscore, BLAST_SCORE_MIN);
        assert_eq!(sbp.hiscore, BLAST_SCORE_MAX);
    }

    #[test]
    fn blast_score_blk_protein_matrix_read_sets_bounds_and_name() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("protein score block");
        assert_eq!(blast_score_blk_protein_matrix_read(&mut sbp, "BLOSUM62"), 0);
        assert!(sbp.protein_alphabet);
        assert_eq!(sbp.name.as_deref(), Some("BLOSUM62"));
        assert_eq!(sbp.matrix.data[1][1], crate::matrix::BLOSUM62[1][1]);
        assert_eq!(sbp.loscore, -4);
        assert_eq!(sbp.hiscore, 11);
    }

    #[test]
    fn blast_score_blk_protein_matrix_load_matches_c_special_residue_rules() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("protein score block");
        sbp.name = Some("BLOSUM62".to_string());
        assert_eq!(blast_score_blk_protein_matrix_load(&mut sbp), 0);

        let gap = crate::encoding::NCBISTDAA_GAP as usize;
        let c = crate::encoding::NCBISTDAA_C as usize;
        let x = crate::encoding::NCBISTDAA_X as usize;
        let u = crate::encoding::NCBISTDAA_U as usize;
        let o = crate::encoding::NCBISTDAA_O as usize;
        let a = crate::encoding::NCBISTDAA_A as usize;

        assert_eq!(sbp.matrix.data[gap][a], BLAST_SCORE_MIN);
        assert_eq!(sbp.matrix.data[a][gap], BLAST_SCORE_MIN);
        assert_eq!(sbp.matrix.data[u][a], sbp.matrix.data[c][a]);
        assert_eq!(sbp.matrix.data[a][u], sbp.matrix.data[a][c]);
        assert_eq!(sbp.matrix.data[o][a], sbp.matrix.data[x][a]);
        assert_eq!(sbp.matrix.data[a][o], sbp.matrix.data[a][x]);
        assert_eq!(sbp.hiscore, 11);

        sbp.name = Some("NOT_A_MATRIX".to_string());
        assert_eq!(blast_score_blk_protein_matrix_load(&mut sbp), 1);
    }

    #[test]
    fn s_blast_score_blk_copy_deep_copies_owned_fields() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("protein score block");
        assert_eq!(blast_score_blk_protein_matrix_read(&mut sbp, "PAM30"), 0);
        sbp.kbp[0] = protein_ungapped_kbp_for_matrix("PAM30");
        sbp.ambiguous_res.push(crate::encoding::NCBISTDAA_X);

        let mut copy = s_blast_score_blk_copy(&sbp);
        assert_eq!(copy.name.as_deref(), Some("PAM30"));
        assert_eq!(copy.matrix.data[1][1], sbp.matrix.data[1][1]);
        assert_eq!(copy.ambiguous_res, sbp.ambiguous_res);
        copy.matrix.data[1][1] = -999;
        copy.kbp[0].lambda = 99.0;
        copy.ambiguous_res.push(crate::encoding::NCBISTDAA_B);

        assert_ne!(copy.matrix.data[1][1], sbp.matrix.data[1][1]);
        assert_ne!(copy.kbp[0].lambda, sbp.kbp[0].lambda);
        assert_ne!(copy.ambiguous_res, sbp.ambiguous_res);
    }

    #[test]
    fn s_blast_find_valid_karlin_blk_returns_first_valid_context() {
        let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[10, 10, 10]);
        query_info.contexts[0].is_valid = false;
        let kbps = vec![
            KarlinBlk::default(),
            KarlinBlk {
                lambda: 0.7,
                k: 0.2,
                log_k: 0.2_f64.ln(),
                h: 0.4,
                round_down: false,
            },
            KarlinBlk {
                lambda: 0.5,
                k: 0.3,
                log_k: 0.3_f64.ln(),
                h: 0.6,
                round_down: false,
            },
        ];

        let kbp = s_blast_find_valid_karlin_blk(&kbps, &query_info).expect("valid Karlin block");
        assert_eq!(kbp.lambda, 0.7);
        assert!(s_blast_karlin_blk_is_valid(Some(kbp)));
        let mut out = None;
        assert_eq!(
            s_blast_find_valid_karlin_blk_c(&kbps, Some(&query_info), Some(&mut out)),
            0
        );
        assert_eq!(out.expect("valid Karlin block").lambda, 0.7);

        let invalid = vec![KarlinBlk::default(); 3];
        assert_eq!(
            s_blast_find_valid_karlin_blk(&invalid, &query_info).unwrap_err(),
            crate::diagnostics::BLASTERR_NOVALIDKARLINALTSCHUL
        );
        let mut out = Some(&kbps[0]);
        assert_eq!(
            s_blast_find_valid_karlin_blk_c(&invalid, Some(&query_info), Some(&mut out)),
            crate::diagnostics::BLASTERR_NOVALIDKARLINALTSCHUL
        );
        assert!(out.is_none());
    }

    #[test]
    fn s_blast_find_smallest_lambda_uses_valid_contexts_only() {
        let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[10, 10, 10]);
        query_info.contexts[1].is_valid = false;
        let kbps = vec![
            KarlinBlk {
                lambda: 0.8,
                k: 0.2,
                log_k: 0.2_f64.ln(),
                h: 0.4,
                round_down: false,
            },
            KarlinBlk {
                lambda: 0.1,
                k: 0.2,
                log_k: 0.2_f64.ln(),
                h: 0.4,
                round_down: false,
            },
            KarlinBlk {
                lambda: 0.6,
                k: 0.3,
                log_k: 0.3_f64.ln(),
                h: 0.6,
                round_down: false,
            },
        ];

        let (lambda, kbp) =
            s_blast_find_smallest_lambda(&kbps, &query_info).expect("smallest lambda");
        assert_eq!(lambda, 0.6);
        assert_eq!(kbp.k, 0.3);
        let mut out = None;
        let lambda = s_blast_find_smallest_lambda_c(&kbps, Some(&query_info), Some(&mut out));
        assert_eq!(lambda, 0.6);
        assert_eq!(out.expect("smallest lambda").k, 0.3);
        let lambda = s_blast_find_smallest_lambda_c(&[], Some(&query_info), None);
        assert_eq!(lambda, i32::MAX as f64);
    }

    #[test]
    fn karlin_and_gumbel_free_helpers_clear_slots() {
        let mut kbp = Some(KarlinBlk {
            lambda: 0.5,
            k: 0.2,
            log_k: 0.2_f64.ln(),
            h: 0.4,
            round_down: false,
        });
        assert!(blast_karlin_blk_free(&mut kbp).is_none());
        assert!(kbp.is_none());

        let mut gbp = Some(GumbelBlk {
            lambda: 0.5,
            a: 1.0,
            b: 2.0,
            alpha: 3.0,
            beta: 4.0,
            sigma: 5.0,
            tau: 6.0,
            db_length: 7,
        });
        assert!(s_blast_gumbel_blk_free(&mut gbp).is_none());
        assert!(gbp.is_none());
    }

    #[test]
    fn blast_score_set_ambig_res_encodes_by_score_block_alphabet() {
        let mut protein =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("protein score block");
        assert_eq!(blast_score_set_ambig_res(Some(&mut protein), b'X'), 0);
        assert_eq!(protein.ambiguous_res, vec![crate::encoding::NCBISTDAA_X]);

        let mut blastna =
            blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 1).expect("blastna score block");
        assert_eq!(blast_score_set_ambig_res(Some(&mut blastna), b'N'), 0);
        assert_eq!(blastna.ambiguous_res, vec![14]);

        let mut ncbi4na =
            blast_score_blk_new(crate::encoding::NCBI4NA_SEQ_CODE, 1).expect("ncbi4na score block");
        assert_eq!(blast_score_set_ambig_res(Some(&mut ncbi4na), b'N'), 0);
        assert_eq!(ncbi4na.ambiguous_res, vec![15]);
        assert_eq!(
            blast_score_set_ambig_res(None, b'X'),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn test_blast_gap_decay_divisor_matches_ncbi_formula() {
        // NCBI `BLAST_GapDecayDivisor(decayrate, nsegs)` =
        // `(1 - decayrate) * decayrate^(nsegs - 1)`. Spot-check a few
        // typical values.
        let eps = 1e-12;
        // nsegs=1 → (1-r) * r^0 = 1-r.
        assert!((blast_gap_decay_divisor(0.1, 1) - 0.9).abs() < eps);
        assert!((blast_gap_decay_divisor(0.5, 1) - 0.5).abs() < eps);
        // nsegs=2 → (1-r) * r.
        assert!((blast_gap_decay_divisor(0.1, 2) - 0.09).abs() < eps);
        assert!((blast_gap_decay_divisor(0.5, 2) - 0.25).abs() < eps);
        // nsegs=3 → (1-r) * r^2.
        assert!((blast_gap_decay_divisor(0.1, 3) - 0.009).abs() < eps);
    }

    #[test]
    fn test_blast_cutoffs_without_decay_matches_evalue_to_raw() {
        // With dodecay=false the adjustment is skipped and the cutoff is
        // just `EtoS(e)` (clamped to `s_floor`).
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let searchsp = 1.0e9;
        let e_in = 10.0;
        let expected_s = kbp.evalue_to_raw(e_in, searchsp);
        let (s, _e_out) = blast_cutoffs(1, e_in, &kbp, searchsp, false, 0.1);
        assert_eq!(s, expected_s);
    }

    #[test]
    fn test_blast_cutoffs_with_decay_tightens_cutoff() {
        // With dodecay=true and gap_decay_rate=0.1, e is multiplied by
        // `blast_gap_decay_divisor(0.1, 1) = 0.9` before conversion, giving a
        // slightly higher raw cutoff score than the no-decay path.
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let searchsp = 1.0e9;
        let e_in = 10.0;
        let (s_no_decay, _) = blast_cutoffs(1, e_in, &kbp, searchsp, false, 0.1);
        let (s_decay, _) = blast_cutoffs(1, e_in, &kbp, searchsp, true, 0.1);
        assert!(
            s_decay >= s_no_decay,
            "decay={s_decay} vs no_decay={s_no_decay}"
        );
        // Independent calculation of the expected decayed cutoff.
        let expected = kbp.evalue_to_raw(0.9 * e_in, searchsp);
        assert_eq!(s_decay, expected);
    }

    #[test]
    fn test_blast_cutoffs_respects_floor() {
        // When the user supplies a floor larger than the stats-derived
        // cutoff, the floor wins and the recomputed e-value corresponds
        // to that floor (with the decay adjustment folded back in).
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let searchsp = 1.0e9;
        let e_in = 10.0;
        let floor = 1000;
        let (s, e_out) = blast_cutoffs(floor, e_in, &kbp, searchsp, true, 0.1);
        assert_eq!(s, floor);
        // e_in did not win → recompute e from s=floor, then divide by
        // the decay divisor (0.9) to get the reported e-value.
        let expected_e = kbp.raw_to_evalue(floor, searchsp) / 0.9;
        assert!(
            (e_out - expected_e).abs() < 1e-14,
            "e_out={e_out} expected={expected_e}"
        );
    }

    #[test]
    fn test_blast_cutoffs_invalid_kbp_returns_one() {
        // Lambda/K/H < 0 → NCBI returns immediately with S=1.
        let kbp = KarlinBlk {
            lambda: -1.0,
            k: 0.41,
            log_k: 0.0,
            h: 0.78,
            round_down: false,
        };
        let (s, e_out) = blast_cutoffs(5, 10.0, &kbp, 1.0e9, true, 0.1);
        assert_eq!(s, 1);
        assert_eq!(e_out, 10.0);
    }

    #[test]
    fn translated_karlin_simple_wrappers_match_c_shape() {
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: true,
        };
        let searchsp = 1_000_000_000_i64;
        let evalue = blast_karlin_sto_e_simple(9, Some(&kbp), searchsp);
        let expected = searchsp as f64 * (-kbp.lambda * 9.0 + kbp.log_k).exp();
        assert!((evalue - expected).abs() < expected.abs() * 1e-12);
        assert_ne!(evalue, kbp.raw_to_evalue(9, searchsp as f64));

        let score = blast_karlin_eto_s_simple(evalue, Some(&kbp), searchsp);
        assert_eq!(score, 9);
        assert_eq!(
            blast_karlin_eto_s_simple(
                1.0,
                Some(&KarlinBlk {
                    lambda: -1.0,
                    ..kbp
                }),
                searchsp
            ),
            BLAST_SCORE_MIN
        );
        assert_eq!(blast_karlin_sto_e_simple(9, None, searchsp), -1.0);

        assert_eq!(blast_karlin_p_to_e(1.0), i32::MAX as f64);
        assert_eq!(blast_karlin_p_to_e(-0.1), i32::MIN as f64);
        assert!(blast_karlin_p_to_e(f64::NAN).is_nan());
        assert_eq!(blast_karlin_e_to_p(0.0), 0.0);
        assert!((blast_karlin_e_to_p(std::f64::consts::LN_2) - 0.5).abs() < 1e-12);

        let mut s = 1;
        let mut e = 10.0;
        assert_eq!(
            blast_cutoffs_in_place(Some(&mut s), Some(&mut e), Some(&kbp), searchsp, false, 0.0),
            0
        );
        let (tuple_s, tuple_e) = blast_cutoffs(1, 10.0, &kbp, searchsp as f64, false, 0.0);
        assert_eq!(s, tuple_s);
        assert_eq!(e, tuple_e);
        assert_eq!(
            blast_cutoffs_in_place(None, Some(&mut e), Some(&kbp), searchsp, false, 0.0),
            1
        );
    }

    #[test]
    fn test_sum_p_r1_formula() {
        // r=1 closed form: p = 1 - exp(-exp(-s)).
        // For s=0, exp(0)=1, 1-exp(-1) ≈ 0.632120.
        let p = s_blast_sum_p(1, 0.0).unwrap();
        assert!((p - (1.0 - (-1.0_f64).exp())).abs() < 1e-12);
        // For s=5, small p.
        let p = s_blast_sum_p(1, 5.0).unwrap();
        assert!(p > 0.0 && p < 0.01);
        // For s=-2, close to 1.
        let p = s_blast_sum_p(1, -2.0).unwrap();
        assert!(p > 0.99);
    }

    #[test]
    fn test_sum_p_r2_tail_formula_against_hand_calc() {
        // r=2, s=6 is in the tail regime (s >= r^2 + r - 1 = 5).
        // NCBI formula: r * exp((r-1)*ln(s) - s - 2*ln(r!)).
        // = 2 * exp(ln(6) - 6 - 2*ln(2)) = 2 * (6/4) * exp(-6) = 3*exp(-6).
        let p = s_blast_sum_p(2, 6.0).unwrap();
        let expected = 3.0 * (-6.0_f64).exp();
        assert!(
            (p - expected).abs() < 1e-12,
            "s_blast_sum_p(2, 6) = {p}, expected {expected}"
        );
    }

    #[test]
    fn test_sum_p_r3_tail_formula_against_hand_calc() {
        // r=3, s=11 is in the tail regime (s >= r^2 + r - 1 = 11).
        // NCBI formula: r * exp((r-1)*ln(s) - s - 2*ln(r!)).
        // = 3 * exp(2*ln(11) - 11 - 2*ln(6)) = 3 * (121/36) * exp(-11).
        let p = s_blast_sum_p(3, 11.0).unwrap();
        let expected = 3.0 * (121.0 / 36.0) * (-11.0_f64).exp();
        assert!(
            (p - expected).abs() < 1e-12,
            "s_blast_sum_p(3, 11) = {p}, expected {expected}"
        );
    }

    #[test]
    fn test_sum_p_r2_table_interpolation() {
        // r=2 should produce a P-value in (0,1] for plausible s.
        for s in [-2.0, 0.0, 2.0, 5.0, 10.0, 15.0] {
            let p = s_blast_sum_p(2, s).unwrap();
            assert!(
                (0.0..=1.0).contains(&p),
                "s_blast_sum_p(2, {}) = {} not in [0,1]",
                s,
                p
            );
        }
        // At very negative s, should saturate to 1.
        assert_eq!(s_blast_sum_p(2, -10.0).unwrap(), 1.0);
    }

    #[test]
    fn test_sum_p_r_gte_5_uses_romberg_branch() {
        // r >= 5 goes through `s_blast_sum_p_calc` (Romberg integration). For
        // very negative s it saturates to 1 via the early-out; for
        // moderate s it produces a valid P-value in (0, 1].
        assert_eq!(s_blast_sum_p(5, -20.0).unwrap(), 1.0);
        let p = s_blast_sum_p(10, 5.0).unwrap();
        assert!(
            (0.0..=1.0).contains(&p),
            "s_blast_sum_p(10, 5.0) = {p} not in [0,1]"
        );
        // r=1 via s_blast_sum_p should still match s_blast_sum_p_calc(r=1).
        let s = 3.0;
        assert!((s_blast_sum_p(1, s).unwrap() - s_blast_sum_p_calc(1, s)).abs() < 1e-14);
    }

    #[test]
    fn test_sum_p_calc_r1_matches_closed_form() {
        // s_blast_sum_p_calc(1, s) has two branches (s>8 and s<=8); they are
        // both small-x approximations of the same Karlin-Altschul
        // P-value and agree to ~6 significant figures at the crossover.
        let a = s_blast_sum_p_calc(1, 8.0);
        let b = s_blast_sum_p_calc(1, 8.0 + 1e-12);
        assert!((a - b).abs() < 1e-6, "a={a} b={b}");
    }

    #[test]
    fn translated_sum_p_romberg_callbacks_match_c_shape() {
        let args = SRombergCbackArgs {
            num_hsps: 2,
            num_hsps_minus_2: 0,
            adj1: -crate::math::blast_ln_gamma_int(1) - crate::math::blast_ln_gamma_int(2),
            adj2: -3.0,
            sdvir: 1.5,
            epsilon: 0.002,
        };
        assert_eq!(
            s_outer_integral_cback(0.0, &args),
            (args.adj2 - (-1.5_f64).exp()).exp()
        );

        let mut args3 = args;
        args3.num_hsps = 3;
        args3.num_hsps_minus_2 = 1;
        assert_eq!(s_outer_integral_cback(0.0, &args3), 0.0);
        let inner = s_inner_integral_cback(3.0, &args);
        assert!(inner.is_finite());
        assert!(inner > 0.0);

        let via_sum = s_blast_sum_p_calc(5, 5.0);
        assert!((0.0..=1.0).contains(&via_sum));
    }

    #[test]
    fn test_sum_e_num_gte_5_now_uses_romberg() {
        // Since s_blast_sum_p handles r >= 5 via Romberg, the sum-E wrappers
        // now return Some for num >= 5 too.
        let e = small_gap_sum_e(40, 5, 30.0, 100, 1000, 1.0e9, 1.0);
        assert!(e.is_some(), "expected Some, got None");
    }

    #[test]
    fn test_karlin_p_to_e_basic() {
        // P = 0 → E = 0.
        assert_eq!(blast_karlin_p_to_e(0.0), 0.0);
        // P = 1 → E = INT4_MAX (matches NCBI `BLAST_KarlinPtoE` exactly).
        assert_eq!(blast_karlin_p_to_e(1.0), i32::MAX as f64);
        // P = 0.5 → E = -ln(0.5) = ln(2) ≈ 0.6931.
        assert!((blast_karlin_p_to_e(0.5) - std::f64::consts::LN_2).abs() < 1e-12);
        // Bad input returns sentinel.
        assert_eq!(blast_karlin_p_to_e(-0.1), i32::MIN as f64);
        assert_eq!(blast_karlin_p_to_e(1.5), i32::MIN as f64);
    }

    #[test]
    fn test_sum_e_num1_reduces_to_searchsp_exp() {
        // For num=1 all three sum-E variants collapse to
        // `searchsp_eff * exp(-xsum) / weight_divisor`.
        let xsum: f64 = 10.0;
        let searchsp: f64 = 1.0e9;
        let w: f64 = 0.9;
        let expected = searchsp * (-xsum).exp() / w;
        for got in [
            small_gap_sum_e(40, 1, xsum, 100, 1000, searchsp, w).unwrap(),
            uneven_gap_sum_e(40, 4000, 1, xsum, 100, 1000, searchsp, w).unwrap(),
            large_gap_sum_e(1, xsum, 100, 1000, searchsp, w).unwrap(),
        ] {
            assert!(
                (got - expected).abs() < 1e-6,
                "got={got} expected={expected}"
            );
        }
    }

    #[test]
    fn translated_sum_e_c_wrappers_match_option_helpers() {
        let q = 100;
        let s = 1000;
        let searchsp_i = 1_000_000_000_i64;
        let searchsp = searchsp_i as f64;
        let w = 0.5;
        let xsum = 120.0;
        assert_eq!(
            blast_small_gap_sum_e(40, 2, xsum, q, s, searchsp_i, w),
            small_gap_sum_e(40, 2, xsum, q, s, searchsp, w).unwrap()
        );
        assert_eq!(
            blast_uneven_gap_sum_e(40, 4000, 2, xsum, q, s, searchsp_i, w),
            uneven_gap_sum_e(40, 4000, 2, xsum, q, s, searchsp, w).unwrap()
        );
        assert_eq!(
            blast_large_gap_sum_e(2, xsum, q, s, searchsp_i, w),
            large_gap_sum_e(2, xsum, q, s, searchsp, w).unwrap()
        );
        assert_eq!(
            blast_small_gap_sum_e(40, 1, 0.0, q, s, searchsp_i, 0.0),
            i32::MAX as f64
        );
    }

    #[test]
    fn test_sum_e_caps_at_int4_max() {
        // Zero weight_divisor saturates to INT4_MAX.
        let capped = small_gap_sum_e(40, 1, 0.0, 100, 1000, 1.0e15, 0.0).unwrap();
        assert_eq!(capped, i32::MAX as f64);
    }

    #[test]
    fn test_large_gap_sum_e_num2_matches_hand_calc() {
        // Reconstruct the NCBI formula (`blast_stat.c:4557-4567`) for
        // num=2 step-by-step and verify `large_gap_sum_e` agrees.
        let q: i32 = 100;
        let s: i32 = 1000;
        let searchsp: f64 = 1.0e9;
        let w: f64 = 0.5;
        let raw_xsum: f64 = 120.0;
        let num: u32 = 2;

        // Manual: adjusted = xsum - num*ln(q*s) + ln_fact(num).
        let adjusted = raw_xsum - num as f64 * ((q as f64) * (s as f64)).ln()
            + crate::math::ln_factorial(num as i32);
        let p = s_blast_sum_p(num, adjusted).unwrap();
        let expected = blast_karlin_p_to_e(p) * (searchsp / (q as f64 * s as f64)) / w;

        let got = large_gap_sum_e(num, raw_xsum, q, s, searchsp, w).unwrap();
        assert!(
            (got - expected).abs() < expected.abs() * 1e-12 + 1e-12,
            "got={got} expected={expected}"
        );
    }

    #[test]
    fn test_small_gap_sum_e_num2_matches_hand_calc() {
        // Hand-verify NCBI `BLAST_SmallGapSumE` (`blast_stat.c:4440-4456`)
        // at num=2. Double-adjust xsum, fold in ln_factorial(2), then
        // `blast_karlin_p_to_e(s_blast_sum_p) * (searchsp / (q*s)) / w`.
        let q: i32 = 100;
        let s: i32 = 1000;
        let searchsp: f64 = 1.0e9;
        let w: f64 = 0.9;
        let raw_xsum: f64 = 120.0;
        let num: u32 = 2;
        let starting_points: i32 = 40;

        let pair_search_space = s as f64 * q as f64;
        let adjusted = raw_xsum
            - pair_search_space.ln()
            - 2.0 * (num - 1) as f64 * (starting_points as f64).ln()
            - crate::math::ln_factorial(num as i32);
        let p = s_blast_sum_p(num, adjusted).unwrap();
        let expected = blast_karlin_p_to_e(p) * (searchsp / pair_search_space) / w;

        let got = small_gap_sum_e(starting_points, num, raw_xsum, q, s, searchsp, w).unwrap();
        assert!(
            (got - expected).abs() < expected.abs() * 1e-12 + 1e-12,
            "got={got} expected={expected}"
        );
    }

    #[test]
    fn test_sum_p_r0_is_zero() {
        assert_eq!(s_blast_sum_p(0, 0.0), Some(0.0));
        assert_eq!(s_blast_sum_p(0, 5.0), Some(0.0));
    }

    #[test]
    fn test_evalue_to_raw_roundtrips_raw_to_evalue() {
        // For any positive raw score, `raw_to_evalue` then `evalue_to_raw`
        // should round-trip within ±1 (because `evalue_to_raw` rounds up
        // and `raw_to_evalue` is smooth).
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let search_space = 1.0e9;
        for raw in [10, 25, 50, 100, 200] {
            let e = kbp.raw_to_evalue(raw, search_space);
            let back = kbp.evalue_to_raw(e, search_space);
            assert!(
                (back - raw).abs() <= 1,
                "raw={raw} e={e} back={back} (should be within 1)"
            );
        }
    }

    #[test]
    fn test_round_down_affects_evalue_not_bit_score() {
        // NCBI `blast_hits.c:1864-1869` applies `score &= ~1` before
        // E-value calculation when `sbp->round_down` is set, but
        // `Blast_HSPListGetBitScores` (`blast_hits.c:1907,1927`) does
        // NOT apply the same mask to bit scores — only a commented-out
        // `#if 0` assertion hints at it. Verify Rust matches that split.
        let mut kbp = KarlinBlk {
            lambda: 0.6,
            k: 0.4,
            log_k: 0.4_f64.ln(),
            h: 0.7,
            round_down: false,
        };
        // Without round_down, odd score 9 is used as-is.
        let bit_odd = kbp.raw_to_bit(9);
        let ev_odd = kbp.raw_to_evalue(9, 1.0e9);
        // Even score 8 is a reference point.
        let bit_even = kbp.raw_to_bit(8);
        let ev_even = kbp.raw_to_evalue(8, 1.0e9);
        assert_ne!(bit_odd, bit_even);
        // With round_down, odd score 9 rounds to 8 for E-value only.
        kbp.round_down = true;
        let bit_rd = kbp.raw_to_bit(9);
        let ev_rd = kbp.raw_to_evalue(9, 1.0e9);
        // Bit score IGNORES round_down (per NCBI):
        assert!((bit_rd - bit_odd).abs() < 1e-14);
        // E-value RESPECTS round_down (per NCBI):
        assert!((ev_rd - ev_even).abs() < ev_even.abs() * 1e-14);
        // round_down=true but even input → bit and evalue unchanged.
        let bit_even_rd = kbp.raw_to_bit(8);
        assert!((bit_even_rd - bit_even).abs() < 1e-14);
        let _ = ev_odd;
    }

    #[test]
    fn blast_hsp_list_get_bit_scores_uses_hsp_context_kbp() {
        fn hsp(score: i32, context: i32) -> crate::hspstream::Hsp {
            crate::hspstream::Hsp {
                score,
                num_ident: 0,
                bit_score: -1.0,
                evalue: 0.0,
                query_offset: 0,
                query_end: 10,
                query_gapped_start: 0,
                subject_offset: 0,
                subject_end: 10,
                subject_gapped_start: 0,
                context,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            }
        }

        let kbp = vec![
            KarlinBlk {
                lambda: 0.2,
                k: 0.5,
                log_k: 0.5_f64.ln(),
                h: 0.1,
                round_down: false,
            },
            KarlinBlk {
                lambda: 0.6,
                k: 0.4,
                log_k: 0.4_f64.ln(),
                h: 0.7,
                round_down: true,
            },
        ];
        let mut list = crate::hspstream::HspList::new(3);
        list.add_hsp(hsp(20, 0));
        list.add_hsp(hsp(9, 1));
        list.add_hsp(hsp(7, -1));

        assert_eq!(
            blast_hsp_list_get_bit_scores(Some(&mut list), true, &kbp),
            0
        );
        assert_eq!(list.hsps[0].bit_score, kbp[0].raw_to_bit(20));
        assert_eq!(list.hsps[1].bit_score, kbp[1].raw_to_bit(9));
        assert_eq!(list.hsps[2].bit_score, kbp[0].raw_to_bit(7));

        assert_eq!(blast_hsp_list_get_bit_scores(None, false, &kbp), 1);

        let mut missing_context = crate::hspstream::HspList::new(4);
        missing_context.add_hsp(hsp(11, 2));
        assert_eq!(
            blast_hsp_list_get_bit_scores(Some(&mut missing_context), false, &kbp),
            1
        );
    }

    #[test]
    fn test_scaled_nucl_gapped_kbp_lookup_sets_round_down_for_2_minus_3() {
        // `(reward, penalty) = (2, -3)` scoring system has only
        // even-score entries, so NCBI sets `sbp->round_down = TRUE`.
        // Verify the returned KBP carries the flag.
        let ungapped = KarlinBlk {
            lambda: 0.55,
            k: 0.21,
            log_k: 0.21_f64.ln(),
            h: 0.46,
            round_down: false,
        };
        // Gap costs that exist in KBPT_2_3.
        let (kbp, rd) = scaled_nucl_gapped_kbp_lookup(4, 4, 2, -3, &ungapped).unwrap();
        assert!(rd, "expected round_down=true for (2,-3)");
        assert!(kbp.round_down, "returned KBP should carry round_down=true");
    }

    #[test]
    fn test_evalue_to_raw_pathological_inputs() {
        // Invalid lambda or searchsp returns `BLAST_SCORE_MIN`, matching
        // NCBI `BlastKarlinEtoS_simple` (`blast_stat.c:4054-4057`).
        let kbp = KarlinBlk {
            lambda: 0.0, // invalid: `denom <= 0.0 || lambda <= 0.0` triggers.
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        assert_eq!(kbp.evalue_to_raw(1.0, 1.0e9), BLAST_SCORE_MIN);
        let kbp_ok = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        // Zero search space → denom == 0 → return `BLAST_SCORE_MIN`.
        assert_eq!(kbp_ok.evalue_to_raw(1.0, 0.0), BLAST_SCORE_MIN);
        // Non-positive e-values are clamped to 1e-297 (NCBI
        // `BlastKarlinEtoS_simple`, `blast_stat.c:4059`) rather than
        // early-returning. The resulting raw score is finite and large.
        let clamped = kbp_ok.evalue_to_raw(-1.0, 1.0e9);
        let expected = kbp_ok.evalue_to_raw(1.0e-297, 1.0e9);
        assert_eq!(clamped, expected);
        assert_eq!(kbp_ok.evalue_to_raw(0.0, 1.0e9), expected);
    }

    #[test]
    fn test_karlin_blk_valid() {
        let kbp = KarlinBlk {
            lambda: 0.208,
            k: 0.049,
            log_k: 0.049_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        assert!(kbp.is_valid());
    }

    #[test]
    fn test_raw_to_bit() {
        let kbp = KarlinBlk {
            lambda: 0.208,
            k: 0.049,
            log_k: 0.049_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let bit = kbp.raw_to_bit(50);
        assert!(bit > 0.0);
        // lambda * 50 / ln(2) - ln(K)/ln(2) ≈ 0.208*50/0.693 - ln(0.049)/0.693
        // ≈ 15.0 + 4.35 ≈ 19.35
        assert!((bit - 19.35).abs() < 0.5);
    }

    #[test]
    fn test_compute_ungapped_lambda() {
        // For reward=1, penalty=-3: known lambda ≈ 1.374
        let lambda = compute_ungapped_lambda(1, -3);
        assert!(
            (lambda - 1.374).abs() < 0.01,
            "lambda for 1/-3 should be ~1.374, got {}",
            lambda
        );

        // For reward=2, penalty=-3: different lambda
        let lambda2 = compute_ungapped_lambda(2, -3);
        assert!(
            lambda2 > 0.0 && lambda2 < lambda,
            "lambda for 2/-3 should be positive and < 1/-3 lambda"
        );
    }

    #[test]
    fn test_compute_ungapped_kbp() {
        let kbp = compute_ungapped_kbp(1, -3);
        assert!(kbp.lambda > 0.0);
        assert!(kbp.k > 0.0);
        assert!(kbp.h > 0.0);
        assert!(kbp.is_valid());
    }

    #[test]
    fn test_search_space() {
        let ss = compute_search_space(100, 1000000, 2000, 20);
        assert!(ss > 0.0);
        // (100-20) * (1000000 - 2000*20) = 80 * 960000 = 76800000
        assert_eq!(ss, 76800000.0);
    }

    // --- Ported from NCBI scoreblk_unit_test.cpp ---

    /// Port of GetScoreBlockNucl: verify nucleotide scoring properties.
    /// NCBI checks: alphabet_size=16, loscore=-3, hiscore=1, penalty=-3, reward=1
    #[test]
    fn test_scoreblk_nucl_properties() {
        let m = crate::matrix::nucleotide_matrix(1, -3);
        // Find min/max scores in the matrix
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..16 {
            for j in 0..16 {
                if m[i][j] < lo {
                    lo = m[i][j];
                }
                if m[i][j] > hi {
                    hi = m[i][j];
                }
            }
        }
        assert_eq!(lo, -3, "nucleotide loscore should be -3");
        assert_eq!(hi, 1, "nucleotide hiscore should be 1");
    }

    /// Port of GetScoreBlockProtein: verify BLOSUM62 scoring properties.
    /// NCBI checks: alphabet_size=BLASTAA_SIZE(28), loscore=-4, hiscore=11
    #[test]
    fn test_scoreblk_protein_properties() {
        let m = &crate::matrix::BLOSUM62;
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..crate::matrix::AA_SIZE {
            for j in 0..crate::matrix::AA_SIZE {
                if m[i][j] < lo {
                    lo = m[i][j];
                }
                if m[i][j] > hi {
                    hi = m[i][j];
                }
            }
        }
        assert_eq!(lo, -4, "BLOSUM62 loscore should be -4");
        assert_eq!(hi, 11, "BLOSUM62 hiscore should be 11");
    }

    /// Port of NuclGappedCalc: verify gapped KBP for reward=1, penalty=-2, gap_open=3, gap_extend=1.
    /// NCBI expects: Lambda≈1.32, K≈0.57, round_down=false
    #[test]
    fn test_nucl_gapped_calc_1_2_3_1() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (kbp, round_down) = scaled_nucl_gapped_kbp_lookup(3, 1, 1, -2, &ungapped).unwrap();
        assert!(!round_down);
        assert!(
            (kbp.lambda - 1.32).abs() < 0.01,
            "Lambda should be ~1.32, got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.57).abs() < 0.01,
            "K should be ~0.57, got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: verify alpha/beta for reward=1, penalty=-2, gap_open=3, gap_extend=1.
    /// NCBI expects: alpha≈1.3, beta≈-1.0
    #[test]
    fn test_scaled_nucl_alpha_beta_1_2_3_1() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (alpha, beta) = scaled_nucl_alpha_beta(1, -2, 3, 1, ungapped.lambda, ungapped.h, true);
        assert!(
            (alpha - 1.3).abs() < 0.01,
            "alpha should be ~1.3, got {}",
            alpha
        );
        assert!(
            (beta - (-1.0)).abs() < 0.01,
            "beta should be ~-1.0, got {}",
            beta
        );
    }

    /// Port of NuclGappedCalc: high gap costs fall back to ungapped params.
    /// reward=1, penalty=-2, gap_open=4, gap_extend=2 → copies ungapped Lambda/K.
    #[test]
    fn test_nucl_gapped_fallback_to_ungapped() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (kbp, round_down) = scaled_nucl_gapped_kbp_lookup(4, 2, 1, -2, &ungapped).unwrap();
        assert!(!round_down);
        assert_eq!(
            kbp.lambda, ungapped.lambda,
            "High gap costs should fall back to ungapped Lambda"
        );
        assert_eq!(
            kbp.k, ungapped.k,
            "High gap costs should fall back to ungapped K"
        );
    }

    /// Port of NuclGappedCalc: ungapped alpha/beta when gap costs exceed table maximum.
    /// Alpha = Lambda/H, beta = 0.0
    #[test]
    fn test_scaled_nucl_alpha_beta_ungapped_fallback() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (alpha, beta) = scaled_nucl_alpha_beta(1, -2, 4, 2, ungapped.lambda, ungapped.h, true);
        assert!(
            (alpha - ungapped.lambda / ungapped.h).abs() < 1e-6,
            "alpha should equal Lambda/H for unsupported gap costs"
        );
        assert_eq!(beta, 0.0);
    }

    /// Port of NuclGappedCalc: scaled-up values (reward=10, penalty=-20, gap_open=30, gap_extend=10).
    /// GCD=10, so internally maps to (1, -2, 3, 1). Lambda should be ~0.132 (1.32/10).
    #[test]
    fn test_nucl_gapped_scaled_values() {
        let ungapped = compute_ungapped_kbp(10, -20);
        let (kbp, round_down) = scaled_nucl_gapped_kbp_lookup(30, 10, 10, -20, &ungapped).unwrap();
        assert!(!round_down);
        assert!(
            (kbp.lambda - 0.132).abs() < 0.01,
            "Lambda for 10/-20/30/10 should be ~0.132 (scaled), got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.57).abs() < 0.01,
            "K for 10/-20/30/10 should be ~0.57 (unscaled), got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: reward=2, penalty=-7, gap_open=4, gap_extend=2 should set round_down=true.
    /// NCBI expects: Lambda≈0.675, K≈0.62, round_down=true
    #[test]
    fn test_nucl_gapped_round_down() {
        let ungapped = compute_ungapped_kbp(2, -7);
        let (kbp, round_down) = scaled_nucl_gapped_kbp_lookup(4, 2, 2, -7, &ungapped).unwrap();
        assert!(round_down, "2/-7 should set round_down=true");
        assert!(kbp.round_down, "returned KBP should preserve round_down");
        assert!(
            (kbp.lambda - 0.675).abs() < 0.01,
            "Lambda should be ~0.675, got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.62).abs() < 0.01,
            "K should be ~0.62, got {}",
            kbp.k
        );
    }

    #[test]
    fn test_nucl_gapped_fallback_preserves_round_down() {
        let ungapped = compute_ungapped_kbp(2, -7);
        let (kbp, round_down) = scaled_nucl_gapped_kbp_lookup(99, 99, 2, -7, &ungapped).unwrap();

        assert!(round_down, "2/-7 table should set round_down=true");
        assert!(
            kbp.round_down,
            "ungapped fallback should preserve round_down"
        );
        assert_eq!(kbp.lambda, ungapped.lambda);
        assert_eq!(kbp.k, ungapped.k);
        assert_eq!(kbp.h, ungapped.h);
    }

    /// Port of NuclGappedCalc: invalid gap costs should return Err.
    /// reward=4, penalty=-5, gap_open=3, gap_extend=2 is not in the table.
    #[test]
    fn test_nucl_gapped_invalid_gap_costs() {
        let ungapped = compute_ungapped_kbp(4, -5);
        let result = scaled_nucl_gapped_kbp_lookup(3, 2, 4, -5, &ungapped);
        assert!(
            result.is_err(),
            "gap_open=3, gap_extend=2 for 4/-5 should fail"
        );
    }

    /// Port of NuclGappedCalc: invalid gap costs should return Err.
    /// reward=1, penalty=-2, gap_open=1, gap_extend=3 is not in the table.
    #[test]
    fn test_nucl_gapped_invalid_gap_costs_2() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let result = scaled_nucl_gapped_kbp_lookup(1, 3, 1, -2, &ungapped);
        assert!(
            result.is_err(),
            "gap_open=1, gap_extend=3 for 1/-2 should fail"
        );
    }

    /// Port of NuclGappedCalc: unsupported reward/penalty pair returns Err.
    /// reward=2, penalty=-1 is not supported.
    #[test]
    fn test_nucl_gapped_unsupported_scores() {
        let ungapped = compute_ungapped_kbp(2, -1);
        let result = scaled_nucl_gapped_kbp_lookup(1, 3, 2, -1, &ungapped);
        assert!(result.is_err(), "2/-1 is not a supported scoring pair");
    }

    /// Port of EqualRewardPenaltyLHtoK: when reward==|penalty|, K should be 1/3.
    /// NCBI uses full Blast_ScoreBlkKbpIdealCalc with DP-based K computation.
    /// The full ungapped_kbp_calc function uses the same algorithm.
    #[test]
    fn test_equal_reward_penalty_k() {
        // Use the full KBP computation via ungapped_kbp_calc
        // Query: uniform ACGT, 1000 bases
        let query: Vec<u8> = (0..1000).map(|i| (i % 4) as u8).collect();
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 1000,
            is_valid: true,
        }];
        let m = crate::matrix::nucleotide_matrix(2, -2);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -2,
            2,
            16,
            &(4..16).collect::<Vec<u8>>(),
            &|i, j| m[i][j],
        );
        let kbp = results[0].as_ref().unwrap();
        assert!(
            (kbp.k - 1.0 / 3.0).abs() < 0.02,
            "K for equal reward/penalty should be ~1/3, got {}",
            kbp.k
        );
    }

    /// Port of BlastResFreqStdCompNucleotideTest: nucleotide base frequencies should be 0.25 each.
    #[test]
    fn test_nucleotide_standard_frequencies() {
        let freqs = &crate::matrix::NT_FREQUENCIES;
        // First 4 bases (A, C, G, T) should each be 0.25
        for i in 0..4 {
            assert!(
                (freqs[i] - 0.25).abs() < 0.001,
                "Base {} frequency should be 0.25, got {}",
                i,
                freqs[i]
            );
        }
        // Ambiguity codes should be 0
        for i in 4..freqs.len() {
            assert!(
                freqs[i].abs() < 0.001,
                "Ambiguity code {} frequency should be 0, got {}",
                i,
                freqs[i]
            );
        }
    }

    /// Port of BlastResFreqStdCompProteinTest: verify specific amino acid frequencies.
    /// Compact AA_FREQUENCIES indices (alphabetical): A=0,C=1,D=2,E=3,F=4,G=5,
    /// H=6,I=7,K=8,L=9,M=10,N=11,P=12,Q=13,R=14,S=15,T=16,V=17,W=18,Y=19.
    /// Expected x100000 from NCBI's `Robinson_prob` (`blast_stat.c:1795`)
    /// divided by 1000: F→3856, M→2243, Y→3216, C→1925.
    #[test]
    fn test_protein_standard_frequencies() {
        let freqs = &crate::matrix::AA_FREQUENCIES;
        let check = |idx: usize, expected_x100k: i32, name: &str| {
            let got = (freqs[idx] * 100000.0).round() as i32;
            assert!(
                (got - expected_x100k).abs() <= 5,
                "{} (idx {}) expected ~{}, got {} (raw: {})",
                name,
                idx,
                expected_x100k,
                got,
                freqs[idx]
            );
        };
        check(1, 1925, "C"); // Cys
        check(4, 3856, "F"); // Phe
        check(10, 2243, "M"); // Met
        check(19, 3216, "Y"); // Tyr
    }

    #[test]
    fn translated_residue_frequency_helpers_match_c_shape() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("score block");
        assert_eq!(blast_score_set_ambig_res(Some(&mut sbp), b'X'), 0);
        let mut residues = [0u8; 20];
        assert_eq!(blast_get_std_alphabet(sbp.alphabet_code, &mut residues), 20);
        assert_eq!(residues[0], crate::encoding::NCBISTDAA_A);
        assert_eq!(
            blast_get_std_alphabet(sbp.alphabet_code, &mut residues[..10]),
            -2
        );

        let mut rfp = blast_res_freq_new(Some(&sbp)).expect("res freq");
        assert_eq!(blast_res_freq_std_comp(Some(&sbp), Some(&mut rfp)), 0);
        assert!((rfp.prob.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(rfp.prob[crate::encoding::NCBISTDAA_C as usize] > 0.0);

        let mut custom = BlastResFreq {
            alphabet_code: sbp.alphabet_code,
            prob: vec![2.0; sbp.alphabet_size],
        };
        assert_eq!(
            blast_res_freq_normalize(Some(&sbp), Some(&mut custom), 2.0),
            0
        );
        assert!((custom.prob.iter().sum::<f64>() - 2.0).abs() < 1e-12);
        custom.prob[0] = -1.0;
        assert_eq!(
            blast_res_freq_normalize(Some(&sbp), Some(&mut custom), 1.0),
            1
        );

        let mut rcp = blast_res_comp_new(Some(&sbp)).expect("res comp");
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDX");
        assert_eq!(
            blast_res_comp_str(Some(&sbp), Some(&mut rcp), Some(&query), query.len() as i32),
            0
        );
        assert_eq!(rcp.comp[crate::encoding::NCBISTDAA_A as usize], 1);
        assert_eq!(rcp.comp[crate::encoding::NCBISTDAA_X as usize], 0);

        let mut query_freq = blast_res_freq_new(Some(&sbp)).expect("res freq");
        assert_eq!(
            blast_res_freq_res_comp(Some(&sbp), Some(&mut query_freq), Some(&rcp)),
            0
        );
        assert!((query_freq.prob.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert_eq!(blast_res_freq_clr(Some(&sbp), Some(&mut query_freq)), 0);
        assert_eq!(query_freq.prob.iter().sum::<f64>(), 0.0);
        assert_eq!(
            blast_res_freq_string(
                Some(&sbp),
                Some(&mut query_freq),
                Some(&query),
                query.len() as i32
            ),
            0
        );
        assert!((query_freq.prob.iter().sum::<f64>() - 1.0).abs() < 1e-12);

        let mut huge_counts = BlastResComp {
            alphabet_code: sbp.alphabet_code,
            comp: vec![0; sbp.alphabet_size],
        };
        huge_counts.comp[crate::encoding::NCBISTDAA_A as usize] = i32::MAX;
        huge_counts.comp[crate::encoding::NCBISTDAA_C as usize] = i32::MAX;
        assert_eq!(
            blast_res_freq_res_comp(Some(&sbp), Some(&mut query_freq), Some(&huge_counts)),
            0
        );
        assert!((query_freq.prob.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert_eq!(query_freq.prob[crate::encoding::NCBISTDAA_A as usize], 0.5);
        assert_eq!(query_freq.prob[crate::encoding::NCBISTDAA_C as usize], 0.5);

        let mut rfp_slot = Some(query_freq);
        assert!(blast_res_freq_free(&mut rfp_slot).is_none());
        let mut rcp_slot = Some(rcp);
        assert!(blast_res_comp_destruct(&mut rcp_slot).is_none());
    }

    #[test]
    fn translated_score_freq_new_validates_score_range() {
        assert_eq!(blast_score_chk(-3, 5), 0);
        assert_eq!(blast_score_chk(0, 5), 1);
        assert_eq!(blast_score_chk(-5, 0), 1);
        assert_eq!(blast_score_chk(BLAST_SCORE_MIN - 1, 5), 1);
        assert_eq!(blast_score_chk(-5, BLAST_SCORE_MAX + 1), 1);
        let sfp = blast_score_freq_new(-3, 5).expect("score freq");
        assert_eq!(sfp.score_min, -3);
        assert_eq!(sfp.score_max, 5);
        assert_eq!(sfp.obs_min, 0);
        assert_eq!(sfp.obs_max, 0);
        assert_eq!(sfp.sprob.len(), 9);
        assert!(blast_score_freq_new(0, 5).is_none());
        assert!(blast_score_freq_new(-5, 0).is_none());
    }

    #[test]
    fn translated_score_freq_calc_and_ungapped_karlin_match_c_shape() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 1).expect("score block");
        sbp.reward = 1;
        sbp.penalty = -2;
        assert_eq!(blast_score_blk_nucl_matrix_create(&mut sbp), 0);

        let mut rfp1 = blast_res_freq_new(Some(&sbp)).expect("freq1");
        let mut rfp2 = blast_res_freq_new(Some(&sbp)).expect("freq2");
        assert_eq!(blast_res_freq_std_comp(Some(&sbp), Some(&mut rfp1)), 0);
        assert_eq!(blast_res_freq_std_comp(Some(&sbp), Some(&mut rfp2)), 0);
        let mut sfp = blast_score_freq_new(sbp.loscore, sbp.hiscore).expect("score freq");

        assert_eq!(
            blast_score_freq_calc(Some(&sbp), Some(&mut sfp), Some(&rfp1), Some(&rfp2)),
            0
        );
        assert_eq!(sfp.obs_min, -2);
        assert_eq!(sfp.obs_max, 1);
        assert!((sfp.p(1) - 0.25).abs() < 1e-12);
        assert!((sfp.p(-2) - 0.75).abs() < 1e-12);
        assert!((sfp.score_avg + 1.25).abs() < 1e-12);

        let mut kbp = blast_karlin_blk_new();
        assert_eq!(
            blast_karlin_blk_ungapped_calc(Some(&mut kbp), Some(&sfp)),
            0
        );
        assert!(kbp.is_valid());
        let lambda_residual =
            0.25 * (kbp.lambda * 1.0).exp() + 0.75 * (kbp.lambda * -2.0).exp() - 1.0;
        assert!(lambda_residual.abs() < 1e-10);

        let mut sfp_slot = Some(sfp);
        assert!(blast_score_freq_free(&mut sfp_slot).is_none());
        assert!(sfp_slot.is_none());
    }

    #[test]
    fn blast_karlin_lh_to_k_named_port_keeps_c_special_case() {
        let mut sfp = SfDist::new(-1, 1);
        sfp.obs_min = -1;
        sfp.obs_max = 1;
        *sfp.p_mut(-1) = 0.75;
        *sfp.p_mut(1) = 0.25;
        sfp.score_avg = -0.5;

        let k = blast_karlin_lh_to_k(&sfp, 1.0, 1.0);
        assert!((k - (1.0 / 3.0)).abs() < 1e-15);
        assert_eq!(blast_karlin_lh_to_k(&sfp, 0.0, 1.0), -1.0);
        sfp.score_avg = 0.0;
        assert_eq!(blast_karlin_lh_to_k(&sfp, 1.0, 1.0), -1.0);
    }

    #[test]
    fn translated_score_blk_kbp_ungapped_calc_fills_context_arrays() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFGHIKLMNPQRSTVWY");
        let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("score block");
        sbp.name = Some("BLOSUM62".to_string());
        assert_eq!(blast_score_blk_matrix_fill(&mut sbp), 0);
        assert_eq!(
            blast_score_set_ambig_res(Some(&mut sbp), crate::encoding::NCBISTDAA_X),
            0
        );
        let mut messages = None;

        assert_eq!(
            blast_score_blk_kbp_ungapped_calc(
                crate::program::BLASTP,
                Some(&mut sbp),
                Some(&query),
                Some(&mut query_info),
                Some(&mut messages),
            ),
            0
        );

        assert!(messages.is_none());
        assert!(query_info.contexts[0].is_valid);
        assert!(sbp.kbp_ideal.as_ref().is_some_and(KarlinBlk::is_valid));
        assert!(sbp.sfp[0].is_some());
        assert!(sbp.kbp_std[0].is_valid());
        assert!(sbp.kbp_psi[0].is_valid());
        assert_eq!(sbp.kbp[0].lambda, sbp.kbp_std[0].lambda);
    }

    #[test]
    fn nucleotide_matrix_fill_rejects_unsupported_read_in_matrix() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 1).expect("score block");
        sbp.read_in_matrix = true;
        sbp.name = Some("NUC.4.4".to_string());
        sbp.reward = 1;
        sbp.penalty = -3;

        assert_ne!(blast_score_blk_matrix_fill(&mut sbp), 0);
        assert_eq!(sbp.hiscore, 0);
        assert_eq!(sbp.loscore, 0);
    }

    #[test]
    fn translated_score_blk_kbp_ungapped_calc_warns_for_invalid_plain_query() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"XXXX");
        let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("score block");
        sbp.name = Some("BLOSUM62".to_string());
        assert_eq!(blast_score_blk_matrix_fill(&mut sbp), 0);
        assert_eq!(
            blast_score_set_ambig_res(Some(&mut sbp), crate::encoding::NCBISTDAA_X),
            0
        );
        let mut messages = None;

        assert_eq!(
            blast_score_blk_kbp_ungapped_calc(
                crate::program::BLASTP,
                Some(&mut sbp),
                Some(&query),
                Some(&mut query_info),
                Some(&mut messages),
            ),
            1
        );

        assert!(!query_info.contexts[0].is_valid);
        let message = messages.expect("warning");
        assert_eq!(message.severity, crate::diagnostics::BlastSeverity::Warning);
        assert_eq!(message.context, 0);
    }

    /// Port of EvalueForProteinFSC: E-values should never be negative.
    /// Uses BLOSUM62 with gap_open=11, gap_extend=1.
    #[test]
    fn test_evalue_never_negative() {
        let kbp = lookup_protein_params(11, 1).unwrap();
        let kbp = KarlinBlk {
            lambda: kbp.lambda,
            k: kbp.k,
            log_k: kbp.k.ln(),
            h: kbp.h,
            round_down: false,
        };
        // Test cases from NCBI stat_unit_test.cpp
        let cases = [
            (1201, 294, 422),
            (1204, 294, 416),
            (1179, 294, 418),
            (2332, 1801, 1671),
        ];
        for &(score, len1, len2) in &cases {
            let ss = len1 as f64 * len2 as f64;
            let evalue = kbp.raw_to_evalue(score, ss);
            assert!(
                evalue >= 0.0,
                "E-value should be >= 0 for score={}, lens=({},{}), got {}",
                score,
                len1,
                len2,
                evalue
            );
        }
    }

    /// Verify gapped KBP lookup for all supported nucleotide reward/penalty combos.
    /// Not all combos have a linear (0,0) entry; test the first affine entry from each table.
    #[test]
    fn test_all_nucl_reward_penalty_combos() {
        // (reward, penalty, gap_open, gap_extend) - first affine entry from each KBPT table
        let cases: &[(i32, i32, i32, i32)] = &[
            (1, -5, 3, 3),
            (1, -4, 1, 2),
            (1, -3, 2, 2),
            (1, -2, 2, 2),
            (1, -1, 3, 2),
            (2, -7, 2, 4),
            (2, -5, 2, 4),
            (2, -3, 4, 4),
            (3, -4, 6, 3),
            (3, -2, 5, 5),
            (4, -5, 6, 5),
            (5, -4, 10, 6),
        ];
        for &(reward, penalty, go, ge) in cases {
            let ungapped = compute_ungapped_kbp(reward, penalty);
            let result = scaled_nucl_gapped_kbp_lookup(go, ge, reward, penalty, &ungapped);
            assert!(
                result.is_ok(),
                "Gapped lookup should work for {}/{}/{}/{}",
                reward,
                penalty,
                go,
                ge
            );
            let (kbp, _) = result.unwrap();
            assert!(
                kbp.is_valid(),
                "KBP should be valid for {}/{}/{}/{}",
                reward,
                penalty,
                go,
                ge
            );
        }
    }

    /// Verify gapped KBP for all entries in KBPT_1_3 table via scaled_nucl_gapped_kbp_lookup.
    #[test]
    fn test_nucl_gapped_all_1_3_entries() {
        let ungapped = compute_ungapped_kbp(1, -3);
        // Entries from KBPT_1_3 (the C-compatible table)
        let entries: &[(i32, i32, f64, f64)] = &[
            (2, 2, 1.37, 0.70),
            (1, 2, 1.35, 0.64),
            (0, 2, 1.25, 0.42),
            (2, 1, 1.34, 0.60),
            (1, 1, 1.21, 0.34),
        ];
        for &(go, ge, expected_lambda, expected_k) in entries {
            let result = scaled_nucl_gapped_kbp_lookup(go, ge, 1, -3, &ungapped);
            assert!(
                result.is_ok(),
                "1/-3 gap_open={} gap_extend={} should work",
                go,
                ge
            );
            let (kbp, _) = result.unwrap();
            assert!(
                (kbp.lambda - expected_lambda).abs() < 0.01,
                "Lambda mismatch for 1/-3/{}/{}: expected {}, got {}",
                go,
                ge,
                expected_lambda,
                kbp.lambda
            );
            assert!(
                (kbp.k - expected_k).abs() < 0.01,
                "K mismatch for 1/-3/{}/{}: expected {}, got {}",
                go,
                ge,
                expected_k,
                kbp.k
            );
        }
    }

    /// Verify protein gapped KBP for all BLOSUM62 table entries.
    #[test]
    fn test_protein_gapped_all_blosum62_entries() {
        for &(go, ge, expected_lambda, expected_k, expected_h, _alpha, _beta) in BLOSUM62_PARAMS {
            let params = lookup_protein_params(go, ge);
            assert!(
                params.is_some(),
                "BLOSUM62 gap_open={} gap_extend={} should be in table",
                go,
                ge
            );
            let p = params.unwrap();
            assert!((p.lambda - expected_lambda).abs() < 1e-6);
            assert!((p.k - expected_k).abs() < 1e-6);
            assert!((p.h - expected_h).abs() < 1e-6);
        }
    }

    #[test]
    fn test_matrix_display_params_use_selected_table_rows() {
        let ungapped =
            lookup_matrix_ungapped_display_params("BLOSUM62").expect("BLOSUM62 ungapped row");
        assert!((ungapped.a - 0.7916).abs() < 1e-12);
        assert!((ungapped.alpha - 4.964660).abs() < 1e-12);
        assert!((ungapped.sigma - 4.964660).abs() < 1e-12);

        let gapped = lookup_matrix_display_params("PAM30", 9, 1).expect("PAM30 9/1 row");
        assert!((gapped.lambda - 0.294).abs() < 1e-12);
        assert!((gapped.a - 0.48).abs() < 1e-12);
        assert_ne!(gapped.alpha, 42.602800);
        assert_ne!(gapped.sigma, 43.636200);
    }

    /// Verify length adjustment converges for typical search parameters.
    #[test]
    fn test_length_adjustment_exact_convergence() {
        let kbp = lookup_protein_params(11, 1).unwrap();
        let (adj, converged) = compute_length_adjustment_exact(
            kbp.k,
            kbp.k.ln(),
            kbp.alpha / kbp.lambda,
            kbp.beta,
            300,
            1_000_000,
            5000,
        );
        assert!(converged, "Length adjustment should converge");
        assert!(adj > 0, "Length adjustment should be positive");
        assert!(adj < 300, "Length adjustment should be < query length");
    }

    #[test]
    fn translated_length_adjustment_wrapper_matches_c_status_shape() {
        let kbp = lookup_protein_params(11, 1).unwrap();
        let mut adj = -1;
        let status = blast_compute_length_adjustment(
            kbp.k,
            kbp.k.ln(),
            kbp.alpha / kbp.lambda,
            kbp.beta,
            300,
            1_000_000,
            5000,
            Some(&mut adj),
        );
        let (expected, converged) = compute_length_adjustment_exact(
            kbp.k,
            kbp.k.ln(),
            kbp.alpha / kbp.lambda,
            kbp.beta,
            300,
            1_000_000,
            5000,
        );
        assert_eq!(adj, expected);
        assert_eq!(status, if converged { 0 } else { 1 });

        assert_eq!(
            blast_compute_length_adjustment(
                kbp.k,
                kbp.k.ln(),
                kbp.alpha / kbp.lambda,
                kbp.beta,
                300,
                1_000_000,
                5000,
                None,
            ),
            1
        );

        let mut adj = 99;
        assert_eq!(
            blast_compute_length_adjustment(
                1.0e-300,
                (1.0e-300_f64).ln(),
                kbp.alpha / kbp.lambda,
                kbp.beta,
                10,
                10,
                1,
                Some(&mut adj),
            ),
            1
        );
        assert_eq!(adj, 0);
    }

    /// Verify E-value <-> raw score round-trip consistency.
    #[test]
    fn test_evalue_raw_score_roundtrip() {
        let kbp = KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let search_space = 1e9;
        // For a given raw score, convert to evalue, then back to raw score
        for raw in [20, 50, 100, 200] {
            let evalue = kbp.raw_to_evalue(raw, search_space);
            let recovered = kbp.evalue_to_raw(evalue, search_space);
            assert!(
                (recovered - raw).abs() <= 1,
                "Round-trip failed: raw={} -> evalue={:.2e} -> recovered={}",
                raw,
                evalue,
                recovered
            );
        }
    }

    /// Verify ungapped KBP computation for protein (BLOSUM62).
    #[test]
    fn test_protein_ungapped_kbp_values() {
        let kbp = protein_ungapped_kbp();
        assert!(kbp.is_valid());
        assert!((kbp.lambda - 0.3176).abs() < 0.01);
        assert!((kbp.k - 0.134).abs() < 0.01);
        assert!((kbp.h - 0.401).abs() < 0.01);
    }

    #[test]
    fn test_protein_ideal_ungapped_kbp_is_computed_not_rounded_table_row() {
        let ideal = protein_ideal_ungapped_kbp_for_matrix("BLOSUM62");
        let table = protein_ungapped_kbp_for_matrix("BLOSUM62");

        assert!(ideal.is_valid());
        assert!((ideal.lambda - 0.317605957635731).abs() < 1e-12);
        assert!((ideal.k - 0.133956144488482).abs() < 1e-12);
        assert!((ideal.h - 0.401214524497119).abs() < 1e-12);
        assert_ne!(ideal.log_k, table.log_k);

        let ideal_evalue = ideal.raw_to_evalue(38, 64.0) / blast_gap_decay_divisor(0.5, 1);
        assert_eq!(format_evalue_for_test(ideal_evalue), "9.83e-05");
    }

    #[test]
    fn test_protein_ideal_ungapped_kbp_uses_selected_matrix() {
        let ideal = protein_ideal_ungapped_kbp_for_matrix("BLOSUM45");
        let table = protein_ungapped_kbp_for_matrix("BLOSUM45");

        assert!(ideal.is_valid());
        assert!((ideal.lambda - table.lambda).abs() > 1e-6);
        assert_ne!(ideal.log_k, table.log_k);
    }

    #[test]
    fn test_protein_ideal_ungapped_kbp_uses_full_selected_score_range() {
        let ideal = protein_ideal_ungapped_kbp_for_matrix("PAM30");
        let table = protein_ungapped_kbp_for_matrix("PAM30");

        assert!(ideal.is_valid());
        assert!((ideal.lambda - table.lambda).abs() > 1e-6);
        assert_ne!(ideal.log_k, table.log_k);
    }

    #[test]
    fn test_query_specific_protein_ungapped_kbp_uses_selected_matrix_bounds() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"CWWYC");

        let blosum62 = query_specific_protein_ungapped_kbp(&query, &crate::matrix::BLOSUM62);
        let blosum45 = query_specific_protein_ungapped_kbp(&query, &crate::matrix::BLOSUM45);

        assert!(blosum45.is_valid());
        assert!(blosum62.is_valid());
        assert!((blosum45.lambda - blosum62.lambda).abs() > 1e-6);
    }

    fn format_evalue_for_test(val: f64) -> String {
        if val < 0.001 {
            let s = format!("{:.2e}", val);
            if let Some(e_pos) = s.find('e') {
                let (mantissa, exp_part) = s.split_at(e_pos);
                let exp_str = &exp_part[1..];
                let (sign, digits) = if let Some(rest) = exp_str.strip_prefix('-') {
                    ("-", rest)
                } else if let Some(rest) = exp_str.strip_prefix('+') {
                    ("", rest)
                } else {
                    ("", exp_str)
                };
                if digits.len() < 2 {
                    return format!(
                        "{}e{}{:02}",
                        mantissa,
                        sign,
                        digits.parse::<u32>().unwrap_or(0)
                    );
                }
            }
            s
        } else {
            format!("{val}")
        }
    }

    /// Verify ungapped KBP calc produces valid results for context-based computation.
    #[test]
    fn test_ungapped_kbp_calc_nucleotide() {
        // A simple BLASTNA-encoded query: ACGTACGT = [0,1,2,3,0,1,2,3]
        let query: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 8,
            is_valid: true,
        }];
        let m = crate::matrix::nucleotide_matrix(1, -3);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -3,
            1,
            16,
            &(4..16).collect::<Vec<u8>>(), // ambiguous codes 4-15
            &|i, j| m[i][j],
        );
        assert_eq!(results.len(), 1);
        assert!(results[0].is_some(), "Should produce valid KBP");
        let kbp = results[0].as_ref().unwrap();
        assert!(kbp.is_valid());
        assert!(kbp.lambda > 0.0);
        assert!(kbp.k > 0.0);
        assert!(kbp.h > 0.0);
    }

    /// Verify invalid context (zero-length) produces None.
    #[test]
    fn test_ungapped_kbp_calc_invalid_context() {
        let query: Vec<u8> = vec![0, 1, 2, 3];
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 0,
            is_valid: false,
        }];
        let m = crate::matrix::nucleotide_matrix(1, -3);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -3,
            1,
            16,
            &(4..16).collect::<Vec<u8>>(),
            &|i, j| m[i][j],
        );
        assert!(results[0].is_none());
    }

    /// Verify search space with zero length adjustment.
    #[test]
    fn test_search_space_no_adjustment() {
        let ss = compute_search_space(500, 10_000_000, 1000, 0);
        assert_eq!(ss, 500.0 * 10_000_000.0);
    }

    /// Verify search space clamps to 1 when adjustment exceeds lengths.
    #[test]
    fn test_search_space_clamped() {
        let ss = compute_search_space(10, 100, 5, 50);
        // eff_query = max(10-50, 1) = 1
        // eff_db = max(100 - 5*50, 1) = max(-150, 1) = 1
        assert_eq!(ss, 1.0);
    }

    #[test]
    fn test_matrix_data_reset_reuses_and_resizes_scratch() {
        let mut m = MatrixData::default();
        assert_eq!(s_matrix_data_reset(Some(&mut m), 11, 0.9), 0);
        assert_eq!(m.matrix_dim, 12);
        assert_eq!(m.matrix_dim_alloc, 12);
        assert_eq!(m.power_matrix.len(), 144);
        assert_eq!(m.prod_matrix.len(), 144);
        assert_eq!(m.percent_identity, 0.9);

        m.hit_probability = 0.5;
        assert_eq!(s_matrix_data_reset(Some(&mut m), 4, 0.75), 0);
        assert_eq!(m.matrix_dim, 5);
        assert_eq!(m.matrix_dim_alloc, 12);
        assert_eq!(m.power_matrix.len(), 144);
        assert_eq!(m.hit_probability, 0.0);
        assert_eq!(s_matrix_data_reset(None, 4, 0.75), -1);
    }

    #[test]
    fn test_matrix_square_matches_naive_product() {
        let a = vec![
            1.0, 2.0, 3.0, 4.0, //
            0.5, 1.5, 2.5, 3.5, //
            4.0, 3.0, 2.0, 1.0, //
            2.0, 0.0, 1.0, 3.0,
        ];
        let mut prod = vec![0.0; 16];
        s_matrix_square(&a, &mut prod, 4);

        let mut expected = vec![0.0; 16];
        for i in 0..4 {
            for j in 0..4 {
                expected[i * 4 + j] = (0..4).map(|k| a[i * 4 + k] * a[k * 4 + j]).sum();
            }
        }
        assert_eq!(prod, expected);
    }

    #[test]
    fn test_find_best_nucleotide_word_size_matches_ncbi_sanity_paths() {
        assert_eq!(blast_find_best_nucleotide_word_size(1.0, 100), 0);
        assert_eq!(blast_find_best_nucleotide_word_size(0.59, 100), 0);
        assert_eq!(blast_find_best_nucleotide_word_size(0.9, -1), 0);
        assert_eq!(blast_find_best_nucleotide_word_size(0.9, 7), 4);

        let mut m = MatrixData::default();
        let direct = s_find_word_size(Some(&mut m), 0.9, 100);
        assert_eq!(blast_find_best_nucleotide_word_size(0.9, 100), direct);
        assert!((4..=110).contains(&direct));
    }

    /// Verify length adjustment is 0 for degenerate inputs.
    #[test]
    fn test_length_adjustment_degenerate() {
        let kbp = KarlinBlk {
            lambda: 0.0,
            k: 0.0,
            log_k: f64::NEG_INFINITY,
            h: 0.0,
            round_down: false,
        };
        let adj = compute_length_adjustment(100, 1000000, 100, &kbp);
        assert_eq!(adj, 0);
    }

    #[test]
    fn translated_matrix_and_table_entry_points_match_c_shape() {
        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("score block");
        sbp.name = Some("BLOSUM62".to_string());
        assert_eq!(blast_score_blk_matrix_fill(&mut sbp), 0);
        assert_eq!(sbp.hiscore, 11);
        assert_eq!(blast_score_blk_kbp_ideal_calc(Some(&mut sbp)), 0);
        assert!(sbp.kbp_ideal.as_ref().is_some_and(KarlinBlk::is_valid));
        assert_eq!(blast_score_blk_kbp_ideal_calc(None), 1);

        let src = KarlinBlk {
            lambda: 0.3,
            k: 0.1,
            log_k: 0.1_f64.ln(),
            h: 0.2,
            round_down: true,
        };
        let mut dst = KarlinBlk::default();
        assert_eq!(blast_karlin_blk_copy(Some(&mut dst), Some(&src)), 0);
        assert_eq!(dst.lambda, src.lambda);
        assert_eq!(dst.round_down, src.round_down);
        assert_eq!(blast_karlin_blk_copy(None, Some(&src)), -1);
    }

    #[test]
    fn translated_gapped_and_gumbel_table_loaders_match_c_statuses() {
        let mut kbp = KarlinBlk::default();
        assert_eq!(
            blast_karlin_blk_gapped_load_from_tables(
                Some(&mut kbp),
                11,
                1,
                Some("BLOSUM62"),
                false,
            ),
            0
        );
        assert!((kbp.lambda - 0.267).abs() < 1e-12);
        assert_eq!(
            blast_karlin_blk_gapped_load_from_tables(None, 11, 1, None, false),
            -1
        );
        assert_eq!(
            blast_karlin_blk_gapped_calc(None, 11, 1, Some("NOT_A_MATRIX"), None),
            1
        );
        assert_eq!(
            blast_karlin_blk_gapped_calc(None, 99, 99, Some("BLOSUM62"), None),
            2
        );

        let mut gbp = GumbelBlk {
            lambda: 0.0,
            a: 0.0,
            b: 0.0,
            alpha: 0.0,
            beta: 0.0,
            sigma: 0.0,
            tau: 0.0,
            db_length: 0,
        };
        assert_eq!(
            blast_gumbel_blk_load_from_tables(Some(&mut gbp), 11, 1, Some("BLOSUM62")),
            0
        );
        assert!((gbp.lambda - 0.267).abs() < 1e-12);
        assert_eq!(blast_gumbel_blk_calc(None, 11, 1, None, None), -1);
        assert_eq!(
            blast_gumbel_blk_calc(None, 99, 99, Some("BLOSUM62"), None),
            2
        );
    }

    #[test]
    fn translated_spouge_wrappers_match_c_shape() {
        let mut kbp = KarlinBlk::default();
        assert_eq!(
            blast_karlin_blk_gapped_load_from_tables(
                Some(&mut kbp),
                11,
                1,
                Some("BLOSUM62"),
                false,
            ),
            0
        );
        let mut gbp = GumbelBlk {
            lambda: 0.0,
            a: 0.0,
            b: 0.0,
            alpha: 0.0,
            beta: 0.0,
            sigma: 0.0,
            tau: 0.0,
            db_length: 1_000_000,
        };
        assert_eq!(
            blast_gumbel_blk_load_from_tables(Some(&mut gbp), 11, 1, Some("BLOSUM62")),
            0
        );
        gbp.db_length = 1_000_000;

        let score = 72;
        let m = 250;
        let n = 400;
        let evalue = blast_spouge_sto_e(score, Some(&kbp), Some(&gbp), m, n);
        let expected = spouge_evalue(score, &kbp, &gbp, m, n);
        assert!((evalue - expected).abs() <= expected.abs() * 1e-12);

        let cutoff = blast_spouge_eto_s(evalue, Some(&kbp), Some(&gbp), m, n);
        assert_eq!(cutoff, spouge_etos(evalue, &kbp, &gbp, m, n));
        assert_eq!(blast_spouge_sto_e(score, None, Some(&gbp), m, n), -1.0);
        assert_eq!(
            blast_spouge_eto_s(evalue, Some(&kbp), None, m, n),
            BLAST_SCORE_MIN
        );
    }

    #[test]
    fn translated_matrix_table_wrappers_report_c_style_errors() {
        let all = blast_load_matrix_values(false);
        assert!(all.iter().any(|info| info.name == "IDENTITY"));
        let standard = blast_load_matrix_values(true);
        assert!(!standard.iter().any(|info| info.name == "IDENTITY"));
        let mut destructed = all.clone();
        assert!(blast_matrix_values_destruct(&mut destructed).is_empty());
        assert!(destructed.is_empty());
        assert_eq!(
            blast_karlin_blk_gapped_load_from_tables(None, 15, 2, Some("IDENTITY"), true,),
            1
        );

        let mut messages = None;
        assert_eq!(
            blast_karlin_blk_gapped_calc(None, 11, 1, Some("NOT_A_MATRIX"), Some(&mut messages),),
            1
        );
        let message = messages.expect("matrix error message");
        assert_eq!(message.severity, crate::diagnostics::BlastSeverity::Error);
        assert_eq!(message.message, "NOT_A_MATRIX is not a supported matrix");
        assert!(message
            .next
            .as_ref()
            .is_some_and(|next| next.message == "BLOSUM80 is a supported matrix"));

        let mut messages = None;
        assert_eq!(
            blast_gumbel_blk_calc(None, 99, 99, Some("BLOSUM62"), Some(&mut messages),),
            2
        );
        let message = messages.expect("gap error message");
        assert!(message.message.contains(
            "Gap existence and extension values of 99 and 99 not supported for BLOSUM62"
        ));
        assert!(message.next.as_ref().is_some_and(|next| {
            next.message == "Gap existence and extension values of 32767 and 32767 are supported"
        }));
    }

    #[test]
    fn translated_nucleotide_gap_helpers_and_alpha_beta_match_c_shape() {
        let values = s_get_nucl_values_array(2, -4).expect("values");
        assert!(!values.normal.is_empty());
        assert_eq!(values.gap_open_max, 4);
        assert_eq!(values.gap_extend_max, 4);
        assert!(blast_check_reward_penalty_scores(1, -2));
        assert!(!blast_check_reward_penalty_scores(2, -1));
        assert_eq!(s_get_nucl_values_array(2, -1).unwrap_err(), -1);

        let mut gap_open = 1;
        let mut gap_extend = 1;
        assert_eq!(
            blast_get_nucleotide_gap_existence_extend_params(1, -2, &mut gap_open, &mut gap_extend,),
            0
        );
        assert_eq!((gap_open, gap_extend), (1, 1));

        let ungapped = compute_ungapped_kbp(1, -2);
        let mut alpha = 0.0;
        let mut beta = 0.0;
        assert_eq!(
            blast_get_nucl_alpha_beta(
                1,
                -2,
                3,
                1,
                ungapped.lambda,
                ungapped.h,
                true,
                &mut alpha,
                &mut beta,
            ),
            0
        );
        assert!((alpha - 1.3).abs() < 1e-6);
        assert!((beta + 1.0).abs() < 1e-6);

        assert_eq!(s_get_ungapped_beta(1, -1), -2.0);
        assert_eq!(s_get_ungapped_beta(2, -3), -2.0);
        assert_eq!(s_get_ungapped_beta(1, -2), 0.0);
        assert_eq!(
            blast_get_nucl_alpha_beta(
                2,
                -1,
                1,
                1,
                ungapped.lambda,
                ungapped.h,
                true,
                &mut alpha,
                &mut beta
            ),
            -1
        );

        assert_eq!(
            blast_get_nucl_alpha_beta(
                1,
                -2,
                99,
                99,
                ungapped.lambda,
                ungapped.h,
                true,
                &mut alpha,
                &mut beta,
            ),
            0
        );
        assert!((alpha - ungapped.lambda / ungapped.h).abs() < 1e-12);
        assert_eq!(beta, 0.0);
    }

    #[test]
    fn translated_nucl_gapped_calc_sets_table_linear_fallback_and_error_message() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let mut kbp = KarlinBlk::default();
        let mut round_down = true;
        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                Some(&mut kbp),
                3,
                1,
                1,
                -2,
                Some(&ungapped),
                Some(&mut round_down),
                None,
            ),
            0
        );
        assert!(!round_down);
        assert!((kbp.lambda - 1.32).abs() < 1e-12);
        assert!((kbp.k - 0.57).abs() < 1e-12);

        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                Some(&mut kbp),
                0,
                0,
                1,
                -2,
                Some(&ungapped),
                Some(&mut round_down),
                None,
            ),
            0
        );
        assert!((kbp.lambda - 1.28).abs() < 1e-12);

        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                Some(&mut kbp),
                99,
                99,
                1,
                -2,
                Some(&ungapped),
                Some(&mut round_down),
                None,
            ),
            0
        );
        assert!((kbp.lambda - ungapped.lambda).abs() < 1e-12);

        let before = kbp.clone();
        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                Some(&mut kbp),
                1,
                3,
                1,
                -2,
                Some(&ungapped),
                Some(&mut round_down),
                None,
            ),
            0
        );
        assert_eq!(kbp.lambda, before.lambda);

        let mut messages = None;
        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                Some(&mut kbp),
                1,
                3,
                1,
                -2,
                Some(&ungapped),
                Some(&mut round_down),
                Some(&mut messages),
            ),
            1
        );
        let message = messages.expect("error message");
        assert_eq!(message.severity, crate::diagnostics::BlastSeverity::Error);
        assert!(message
            .message
            .contains("are not supported for substitution scores 1 and -2"));
        assert!(message
            .message
            .contains("Any values more stringent than 2 and 2 are supported"));
    }

    #[test]
    fn translated_residue_probability_and_rps_fill_scores_match_c_shape() {
        let sequence = [
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
            crate::encoding::NCBISTDAA_X,
        ];
        let mut probs = vec![0.0; crate::encoding::BLASTAA_SIZE];
        blast_fill_residue_probability(&sequence, &mut probs);
        assert!((probs[crate::encoding::NCBISTDAA_A as usize] - 2.0 / 3.0).abs() < 1e-12);
        assert!((probs[crate::encoding::NCBISTDAA_C as usize] - 1.0 / 3.0).abs() < 1e-12);
        assert_eq!(probs[crate::encoding::NCBISTDAA_X as usize], 0.0);
        assert!((rps_find_ungapped_lambda(Some("BLOSUM62")) - 0.3176).abs() < 1e-12);
        assert_eq!(rps_find_ungapped_lambda(Some("NOT_A_MATRIX")), 0.0);
        assert_eq!(rps_find_ungapped_lambda(None), 0.0);

        let matrix = vec![vec![-1, 2, BLAST_SCORE_MIN], vec![3, -2, 1]];
        let query_probs = vec![0.25, 0.75, 0.0];
        let mut sfp = blast_score_freq_new(-5, 5).expect("score freq");
        rps_fill_scores(&matrix, &query_probs, &mut sfp, 3);
        assert_eq!(sfp.obs_min, -2);
        assert_eq!(sfp.obs_max, 3);
        assert!((sfp.score_avg - 0.25).abs() < 1e-12);

        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("score block");
        sbp.name = Some("BLOSUM62".to_string());
        let pos_matrix: Vec<Vec<i32>> = crate::matrix::BLOSUM62
            .iter()
            .map(|row| row.to_vec())
            .collect();
        let mut full_sfp = ScoreFreq {
            score_min: 0,
            score_max: 0,
            obs_min: 0,
            obs_max: 0,
            score_avg: 0.0,
            sprob: Vec::new(),
        };
        rps_fill_scores(
            &pos_matrix,
            &probs,
            &mut full_sfp,
            crate::encoding::BLASTAA_SIZE,
        );
        assert!(blast_karlin_lambda_nr(Some(&full_sfp), 0.3) > 0.0);
        assert_eq!(blast_karlin_lambda_nr(None, 0.3), -1.0);
        let rescaled = rps_rescale_pssm(
            2.0,
            sequence.len() as i32,
            Some(&sequence),
            pos_matrix.len() as i32,
            Some(&pos_matrix),
            Some(&sbp),
        )
        .expect("rescaled pssm");
        assert_eq!(rescaled.len(), pos_matrix.len());
        assert_eq!(rescaled[0].len(), crate::encoding::BLASTAA_SIZE);
        assert_eq!(
            rescaled[0][crate::encoding::NCBISTDAA_X as usize],
            pos_matrix[0][crate::encoding::NCBISTDAA_X as usize]
        );
        assert_ne!(
            rescaled[0][crate::encoding::NCBISTDAA_A as usize],
            pos_matrix[0][crate::encoding::NCBISTDAA_A as usize]
        );
        assert!(
            rps_rescale_pssm(0.0, 1, Some(&sequence), 1, Some(&pos_matrix), Some(&sbp)).is_none()
        );
        assert!(
            rps_rescale_pssm(2.0, 1, Some(&sequence), 3, Some(&pos_matrix), Some(&sbp)).is_none()
        );
        let mut short_row_matrix = pos_matrix.clone();
        short_row_matrix[0].truncate(crate::encoding::BLASTAA_SIZE - 1);
        assert!(rps_rescale_pssm(
            2.0,
            1,
            Some(&sequence),
            1,
            Some(&short_row_matrix),
            Some(&sbp)
        )
        .is_none());
        assert!(
            rps_rescale_pssm(2.0, -1, Some(&sequence), 1, Some(&pos_matrix), Some(&sbp)).is_none()
        );
        assert!(rps_rescale_pssm(2.0, 1, None, 1, Some(&pos_matrix), Some(&sbp)).is_none());
        assert!(rps_rescale_pssm(2.0, 1, Some(&sequence), 1, None, Some(&sbp)).is_none());
        assert!(rps_rescale_pssm(2.0, 1, Some(&sequence), 1, Some(&pos_matrix), None).is_none());
        sbp.name = Some("NOT_A_MATRIX".to_string());
        assert!(
            rps_rescale_pssm(2.0, 1, Some(&sequence), 1, Some(&pos_matrix), Some(&sbp)).is_none()
        );
    }

    #[test]
    fn translated_compressed_alphabet_helpers_match_c_shape() {
        let mut table = vec![0u8; crate::encoding::BLASTAA_SIZE];
        let mut rev_table =
            [[-1i8; COMPRESSED_REVERSE_LOOKUP_SIZE]; COMPRESSED_REVERSE_LOOKUP_SIZE];
        s_build_compressed_translation(S_ALPHABET10, &mut table, 10, &mut rev_table);

        assert_eq!(table[crate::encoding::NCBISTDAA_I as usize], 0);
        assert_eq!(table[crate::encoding::NCBISTDAA_J as usize], 0);
        assert_eq!(table[crate::encoding::NCBISTDAA_L as usize], 0);
        assert_eq!(
            table[crate::encoding::aminoacid_to_ncbistdaa_base(b'M') as usize],
            0
        );
        assert_eq!(
            table[crate::encoding::aminoacid_to_ncbistdaa_base(b'V') as usize],
            0
        );
        assert_eq!(table[crate::encoding::NCBISTDAA_A as usize], 1);
        assert_eq!(
            table[crate::encoding::aminoacid_to_ncbistdaa_base(b'W') as usize],
            9
        );
        assert_eq!(rev_table[0][0], crate::encoding::NCBISTDAA_I as i8);
        assert_eq!(rev_table[0][5], -1);

        let mut sbp =
            blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("score block");
        sbp.name = Some("BLOSUM62".to_string());
        assert_eq!(
            blast_score_set_ambig_res(Some(&mut sbp), crate::encoding::NCBISTDAA_X),
            0
        );
        let mut compressed_prob = vec![0.0; crate::encoding::BLASTAA_SIZE];
        assert_eq!(
            s_get_compressed_probs(Some(&sbp), &mut compressed_prob, 10, &rev_table),
            0
        );
        for letter in 0..10 {
            let mut sum = 0.0;
            for &aa in &rev_table[letter] {
                if aa < 0 {
                    break;
                }
                sum += compressed_prob[aa as usize];
            }
            assert!(
                (sum - 1.0).abs() < 1e-12,
                "compressed letter {letter} conditional sum {sum}"
            );
        }
        assert_eq!(
            s_get_compressed_probs(None, &mut compressed_prob, 10, &rev_table),
            -1
        );

        let mut alphabet = SCompressedAlphabet {
            compressed_alphabet_size: 10,
            compress_table: table.clone(),
            matrix: None,
        };
        assert_eq!(
            s_build_compressed_score_matrix(Some(&sbp), Some(&mut alphabet), 1.0, &rev_table),
            0
        );
        let matrix = alphabet.matrix.as_ref().expect("compressed score matrix");
        assert_eq!(matrix.nrows, crate::encoding::BLASTAA_SIZE);
        assert_eq!(matrix.ncols, 10);
        let ratios = crate::matrix::get_matrix_freq_ratios("BLOSUM62").unwrap();
        let compressed_a = table[crate::encoding::NCBISTDAA_A as usize] as usize;
        let mut ratio = 0.0;
        for &aa in &rev_table[compressed_a] {
            if aa < 0 {
                break;
            }
            let aa = aa as usize;
            ratio += ratios[crate::encoding::NCBISTDAA_A as usize][aa] * compressed_prob[aa];
        }
        let lambda = rps_find_ungapped_lambda(Some("BLOSUM62"));
        let expected = crate::math::blast_nint(ratio.ln() / lambda) as i32;
        assert_eq!(
            matrix.data[crate::encoding::NCBISTDAA_A as usize][compressed_a],
            expected
        );
        assert_eq!(
            s_build_compressed_score_matrix(None, Some(&mut alphabet), 1.0, &rev_table),
            -1
        );

        let allocated =
            s_compressed_alphabet_new(Some(&sbp), 15, 1.0).expect("compressed alphabet");
        assert_eq!(allocated.compressed_alphabet_size, 15);
        assert_eq!(
            allocated.compress_table.len(),
            crate::encoding::BLASTAA_SIZE
        );
        assert_eq!(allocated.matrix.as_ref().unwrap().ncols, 15);
        assert!(s_compressed_alphabet_new(Some(&sbp), 11, 1.0).is_none());
        let mut allocated = Some(allocated);
        assert!(s_compressed_alphabet_free(&mut allocated).is_none());
        assert!(allocated.is_none());
    }

    #[test]
    fn translated_allowed_values_helpers_report_matrix_rows() {
        let info = matrix_info_new("BLOSUM62").expect("matrix info");
        assert_eq!(info.name, "BLOSUM62");
        assert!(info.max_number_values > 0);
        let mut info_slot = Some(info.clone());
        assert!(matrix_info_destruct(&mut info_slot).is_none());
        assert!(info_slot.is_none());
        let matrix_message = blast_print_matrix_message("NOT_A_MATRIX", true);
        assert!(matrix_message
            .starts_with("NOT_A_MATRIX is not a supported matrix, supported matrices are:\n"));
        assert!(matrix_message.contains("BLOSUM62 \n"));
        assert!(!matrix_message.contains("IDENTITY \n"));

        let allowed = blast_print_allowed_values("BLOSUM62", 1, 1);
        assert!(allowed.contains("supported values"));
        assert!(allowed.contains("32767, 32767\n") || allowed.contains("32767, 32767,"));
        assert!(allowed.contains("11, 1\n"));
        assert!(!blast_karlin_report_allowed_values("BLOSUM62").is_empty());
    }

    #[test]
    fn translated_matrix_values_and_alpha_beta_match_c_table_shape() {
        let values = blast_get_matrix_values(Some("BLOSUM62"));
        assert_eq!(
            values.open.len(),
            matrix_stat_rows("BLOSUM62").unwrap().len()
        );
        assert_eq!(values.open[9], 11);
        assert_eq!(values.extension[9], 1);
        assert_eq!(values.pref_flags[9], BLAST_MATRIX_BEST);
        assert!(values
            .pref_flags
            .iter()
            .enumerate()
            .all(|(index, flag)| index == 9 || *flag == BLAST_MATRIX_NOMINAL));
        assert!(blast_get_matrix_values(None).open.is_empty());
        assert!(blast_get_matrix_values(Some("NOT_A_MATRIX"))
            .open
            .is_empty());

        let mut alpha = -1.0;
        let mut beta = -1.0;
        blast_get_alpha_beta(Some("BLOSUM62"), &mut alpha, &mut beta, true, 0, 0, None);
        assert!((alpha - values.alpha[9]).abs() < 1e-12);
        assert!((beta - values.beta[9]).abs() < 1e-12);

        blast_get_alpha_beta(Some("BLOSUM62"), &mut alpha, &mut beta, true, 10, 1, None);
        let row = matrix_stat_row_for_gap_costs("BLOSUM62", 10, 1).expect("row");
        assert!((alpha - row.alpha).abs() < 1e-12);
        assert!((beta - row.beta).abs() < 1e-12);

        let fallback = KarlinBlk {
            lambda: 0.5,
            h: 0.25,
            k: 0.1,
            log_k: 0.1_f64.ln(),
            round_down: false,
        };
        blast_get_alpha_beta(
            Some("NOT_A_MATRIX"),
            &mut alpha,
            &mut beta,
            false,
            0,
            0,
            Some(&fallback),
        );
        assert_eq!(alpha, 2.0);
        assert_eq!(beta, 0.0);
    }
}
