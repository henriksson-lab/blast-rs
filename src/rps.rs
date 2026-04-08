//! RPS-BLAST (Reverse Position-Specific BLAST) support.
//! Searches a query sequence against a database of PSSMs (domain profiles).

use crate::pssm::Pssm;
#[cfg(test)]
use crate::matrix::AA_SIZE;

/// An RPS database entry (domain profile).
#[derive(Debug, Clone)]
pub struct RpsProfile {
    pub id: String,
    pub description: String,
    pub pssm: Pssm,
}

/// Search a query against a collection of RPS profiles.
pub fn rps_blast_search(
    query: &[u8], // NCBIstdaa encoded
    profiles: &[RpsProfile],
    evalue_threshold: f64,
    kbp_lambda: f64,
    kbp_k: f64,
) -> Vec<RpsHit> {
    let mut hits = Vec::new();

    for profile in profiles {
        let pssm = &profile.pssm;
        if query.len() < pssm.length { continue; }

        // Scan query with this PSSM
        for qi in 0..=(query.len() - pssm.length) {
            let mut score = 0i32;
            for k in 0..pssm.length {
                score += pssm.score_at(k, query[qi + k]);
            }

            if score > 0 {
                let search_space = (query.len() * pssm.length) as f64;
                let evalue = kbp_k * search_space * (-kbp_lambda * score as f64).exp();
                if evalue <= evalue_threshold {
                    let mut ident = 0;
                    // Count "identities" as positions where the PSSM score is positive
                    for k in 0..pssm.length {
                        if pssm.score_at(k, query[qi + k]) > 0 { ident += 1; }
                    }
                    hits.push(RpsHit {
                        profile_id: profile.id.clone(),
                        query_start: qi as i32 + 1,
                        query_end: (qi + pssm.length) as i32,
                        score,
                        evalue,
                        bit_score: (kbp_lambda * score as f64 - kbp_k.ln()) / std::f64::consts::LN_2,
                        num_ident: ident,
                        align_length: pssm.length as i32,
                    });
                }
            }
        }
    }

    hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    hits
}

/// A hit from RPS-BLAST search.
#[derive(Debug, Clone)]
pub struct RpsHit {
    pub profile_id: String,
    pub query_start: i32,
    pub query_end: i32,
    pub score: i32,
    pub evalue: f64,
    pub bit_score: f64,
    pub num_ident: i32,
    pub align_length: i32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rps_blast() {
        let mut matrix = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 1..21 { matrix[i][i] = 5; }
        for i in 1..21 { for j in 1..21 { if i != j { matrix[i][j] = -2; } } }

        let profile_seq = vec![1u8, 2, 3, 4, 5]; // ABCDE
        let pssm = Pssm::from_sequence(&profile_seq, &matrix);
        let profiles = vec![RpsProfile {
            id: "domain1".to_string(),
            description: "test domain".to_string(),
            pssm,
        }];

        // Query contains the domain
        let query = vec![10u8, 10, 1, 2, 3, 4, 5, 10, 10];
        let hits = rps_blast_search(&query, &profiles, 1.0, 0.3176, 0.134);
        assert!(!hits.is_empty(), "Should find domain");
        assert_eq!(hits[0].profile_id, "domain1");
        assert_eq!(hits[0].query_start, 3); // 1-based
    }
}
