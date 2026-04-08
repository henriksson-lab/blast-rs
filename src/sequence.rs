//! Sequence encoding, complement, and manipulation utilities.


/// Encode an IUPAC nucleotide string to BLASTNA.
pub fn encode_nucleotide(iupac: &[u8]) -> Vec<u8> {
    iupac.iter().map(|&b| iupac_to_blastna(b)).collect()
}

/// Convert IUPAC character to BLASTNA encoding.
#[inline]
pub fn iupac_to_blastna(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'R' | b'r' => 4,
        b'Y' | b'y' => 5,
        b'M' | b'm' => 6,
        b'K' | b'k' => 7,
        b'W' | b'w' => 8,
        b'S' | b's' => 9,
        b'B' | b'b' => 10,
        b'D' | b'd' => 11,
        b'H' | b'h' => 12,
        b'V' | b'v' => 13,
        b'N' | b'n' => 14,
        _ => 15,
    }
}

/// Complement a BLASTNA-encoded nucleotide.
#[inline]
pub fn complement_blastna(b: u8) -> u8 {
    const TABLE: [u8; 16] = [3,2,1,0, 5,4,7,6, 8,9,13,12,11,10,14,15];
    if b < 16 { TABLE[b as usize] } else { 15 }
}

/// Reverse complement a BLASTNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_blastna(b)).collect()
}

/// Encode a query for blastn: returns (plus_strand, minus_strand) in BLASTNA.
pub fn encode_query_blastn(iupac: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let plus = encode_nucleotide(iupac);
    let minus = reverse_complement(&plus);
    (plus, minus)
}

/// Build the concatenated query block for blastn (with sentinel bytes).
/// Returns (encoded_block, context_offsets) where context_offsets has
/// (offset_relative_to_sequence, length) for each context.
pub fn build_query_block(queries: &[&[u8]]) -> (Vec<u8>, Vec<(i32, i32)>) {
    let mut encoded = Vec::new();
    let mut offsets = Vec::new();

    for query in queries {
        let (plus, minus) = encode_query_blastn(query);

        // Leading sentinel
        encoded.push(15);
        let start = encoded.len();
        encoded.extend_from_slice(&plus);
        offsets.push(((start - 1) as i32, plus.len() as i32)); // offset relative to sequence_start+1

        // Middle sentinel
        encoded.push(15);
        let start = encoded.len();
        encoded.extend_from_slice(&minus);
        offsets.push(((start - 1) as i32, minus.len() as i32));
    }
    // Trailing sentinel
    encoded.push(15);

    (encoded, offsets)
}

/// BLASTNA to IUPAC character for display.
pub fn blastna_to_iupac(b: u8) -> char {
    match b {
        0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T',
        4 => 'R', 5 => 'Y', 6 => 'M', 7 => 'K',
        8 => 'W', 9 => 'S', 10 => 'B', 11 => 'D',
        12 => 'H', 13 => 'V', 14 => 'N', _ => '-',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_nucleotide() {
        assert_eq!(encode_nucleotide(b"ACGT"), vec![0, 1, 2, 3]);
        assert_eq!(encode_nucleotide(b"acgt"), vec![0, 1, 2, 3]);
        assert_eq!(encode_nucleotide(b"N"), vec![14]);
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement_blastna(0), 3); // A -> T
        assert_eq!(complement_blastna(1), 2); // C -> G
        assert_eq!(complement_blastna(2), 1); // G -> C
        assert_eq!(complement_blastna(3), 0); // T -> A
        assert_eq!(complement_blastna(14), 14); // N -> N
    }

    #[test]
    fn test_reverse_complement() {
        let seq = encode_nucleotide(b"ACGT");
        let rc = reverse_complement(&seq);
        assert_eq!(rc, encode_nucleotide(b"ACGT")); // ACGT is its own RC
    }

    #[test]
    fn test_reverse_complement_asymmetric() {
        let seq = encode_nucleotide(b"AACG");
        let rc = reverse_complement(&seq);
        assert_eq!(rc, encode_nucleotide(b"CGTT"));
    }

    #[test]
    fn test_encode_query_blastn() {
        let (plus, minus) = encode_query_blastn(b"ACG");
        assert_eq!(plus, vec![0, 1, 2]);
        assert_eq!(minus, vec![1, 2, 3]); // CGT
    }
}
