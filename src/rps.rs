//! RPS-BLAST (Reverse Position-Specific BLAST) support.
//! Searches a query sequence against a database of PSSMs (domain profiles).

#[cfg(test)]
use crate::matrix::AA_SIZE;
use crate::pssm::Pssm;

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

    // Total database length = sum of all profile PSSM lengths
    let total_db_len: usize = profiles.iter().map(|p| p.pssm.length).sum();

    for profile in profiles {
        let pssm = &profile.pssm;
        if query.len() < pssm.length {
            continue;
        }

        // Scan query with this PSSM
        for qi in 0..=(query.len() - pssm.length) {
            let mut score = 0i32;
            for k in 0..pssm.length {
                score += pssm.score_at(k, query[qi + k]);
            }

            if score > 0 {
                let search_space = (query.len() * total_db_len) as f64;
                let evalue = kbp_k * search_space * (-kbp_lambda * score as f64).exp();
                if evalue <= evalue_threshold {
                    let mut ident = 0;
                    // Count identities: positions where query matches the PSSM consensus residue
                    for k in 0..pssm.length {
                        if let Some(consensus) = pssm.consensus_at(k) {
                            if query[qi + k] == consensus {
                                ident += 1;
                            }
                        }
                    }
                    hits.push(RpsHit {
                        profile_id: profile.id.clone(),
                        query_start: qi as i32 + 1,
                        query_end: (qi + pssm.length) as i32,
                        score,
                        evalue,
                        bit_score: (kbp_lambda * score as f64 - kbp_k.ln())
                            / std::f64::consts::LN_2,
                        num_ident: ident,
                        align_length: pssm.length as i32,
                    });
                }
            }
        }
    }

    hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
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
    fn test_rps_search_space_calculation() {
        // Verify search space uses sum of profile lengths (total_db_len).
        // search_space = query.len() * sum(profile.pssm.length for all profiles)
        let mut matrix = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 1..21 {
            matrix[i][i] = 5;
        }
        for i in 1..21 {
            for j in 1..21 {
                if i != j {
                    matrix[i][j] = -2;
                }
            }
        }

        let profile1_seq = vec![1u8, 2, 3]; // length 3
        let profile2_seq = vec![4u8, 5, 6, 7, 8]; // length 5
        let pssm1 = Pssm::from_sequence(&profile1_seq, &matrix);
        let pssm2 = Pssm::from_sequence(&profile2_seq, &matrix);

        let profiles = vec![
            RpsProfile {
                id: "dom1".to_string(),
                description: "d1".to_string(),
                pssm: pssm1,
            },
            RpsProfile {
                id: "dom2".to_string(),
                description: "d2".to_string(),
                pssm: pssm2,
            },
        ];

        // Total db length = 3 + 5 = 8
        let total_db_len: usize = profiles.iter().map(|p| p.pssm.length).sum();
        assert_eq!(
            total_db_len, 8,
            "Total DB length should be sum of profile lengths"
        );

        // Query of length 10 -> search_space = 10 * 8 = 80
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let expected_search_space = (query.len() * total_db_len) as f64;
        assert!(
            (expected_search_space - 80.0).abs() < f64::EPSILON,
            "Search space should be query_len * total_db_len = 80"
        );

        // Run the actual search and verify hits use this search space
        // by checking the E-value formula: E = K * search_space * exp(-lambda * score)
        let lambda = 0.3176;
        let k = 0.134;
        let hits = rps_blast_search(&query, &profiles, 10.0, lambda, k);

        // Profile 1 matches at position 0 (query starts with 1,2,3)
        if let Some(hit) = hits.iter().find(|h| h.profile_id == "dom1") {
            let expected_e = k * expected_search_space * (-lambda * hit.score as f64).exp();
            assert!(
                (hit.evalue - expected_e).abs() < 1e-6,
                "E-value should use total DB length in search space. Expected {}, got {}",
                expected_e,
                hit.evalue
            );
        }
    }

    #[test]
    fn test_rps_blast() {
        let mut matrix = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 1..21 {
            matrix[i][i] = 5;
        }
        for i in 1..21 {
            for j in 1..21 {
                if i != j {
                    matrix[i][j] = -2;
                }
            }
        }

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
