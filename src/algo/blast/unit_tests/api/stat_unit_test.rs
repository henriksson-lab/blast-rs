// Upstream source:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/stat_unit_test.cpp

#[cfg(test)]
mod tests {
    use crate::algo::blast::core::blast_stat::{
        blast_gumbel_blk_calc, blast_karlin_blk_free, blast_karlin_blk_gapped_calc,
        blast_karlin_blk_new, blast_spouge_sto_e, GumbelBlk,
    };

    #[test]
    fn evalue_for_protein_fsc() {
        let matrix = "BLOSUM62";
        let gap_open = 11;
        let gap_extend = 1;
        let scores = [1201, 1204, 1179, 2332];
        let len1 = [294, 294, 294, 1801];
        let len2 = [422, 416, 418, 1671];

        let mut kbp = Some(blast_karlin_blk_new());
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
            0,
            blast_karlin_blk_gapped_calc(
                kbp.as_mut(),
                gap_open,
                gap_extend,
                Some(matrix),
                None
            )
        );
        assert_eq!(
            0,
            blast_gumbel_blk_calc(Some(&mut gbp), gap_open, gap_extend, Some(matrix), None)
        );

        for index in 0..scores.len() {
            let evalue = blast_spouge_sto_e(
                scores[index],
                kbp.as_ref(),
                Some(&gbp),
                len1[index],
                len2[index],
            );
            assert!(evalue >= 0.0);
        }

        assert!(blast_karlin_blk_free(&mut kbp).is_none());
    }
}
