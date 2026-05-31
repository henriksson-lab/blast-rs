/// Encode a Blast-def-line-set ASN.1 BER header matching NCBI format.
pub(crate) fn encode_defline_asn1(title: &str, oid: i32) -> Vec<u8> {
    let mut buf = Vec::new();
    // Blast-def-line-set ::= SEQUENCE OF Blast-def-line.
    buf.extend_from_slice(&[0x30, 0x80]);
    {
        // Blast-def-line ::= SEQUENCE.
        buf.extend_from_slice(&[0x30, 0x80]);
        {
            // title [0] VisibleString.
            buf.extend_from_slice(&[0xa0, 0x80]);
            buf.push(0x1a);
            let title_bytes = title.as_bytes();
            encode_asn1_length(&mut buf, title_bytes.len());
            buf.extend_from_slice(title_bytes);
            buf.extend_from_slice(&[0x00, 0x00]);

            // seqid [1] SET OF Seq-id.
            buf.extend_from_slice(&[0xa1, 0x80]);
            {
                if let Some(first_token) = title.split_whitespace().next() {
                    encode_seqid_chain(&mut buf, first_token);
                }

                // Seq-id ::= CHOICE { general Dbtag }.
                buf.extend_from_slice(&[0x30, 0x80]);
                {
                    buf.extend_from_slice(&[0xaa, 0x80]);
                    {
                        // Dbtag ::= SEQUENCE { db VisibleString, tag Object-id }.
                        buf.extend_from_slice(&[0x30, 0x80]);
                        {
                            // db [0] VisibleString "BL_ORD_ID".
                            buf.extend_from_slice(&[0xa0, 0x80]);
                            buf.push(0x1a);
                            buf.push(9);
                            buf.extend_from_slice(b"BL_ORD_ID");
                            buf.extend_from_slice(&[0x00, 0x00]);

                            // tag [1] Object-id ::= CHOICE { id INTEGER }.
                            buf.extend_from_slice(&[0xa1, 0x80]);
                            {
                                buf.extend_from_slice(&[0xa0, 0x80]);
                                buf.push(0x02);
                                let oid_bytes = encode_asn1_integer(oid);
                                encode_asn1_length(&mut buf, oid_bytes.len());
                                buf.extend_from_slice(&oid_bytes);
                                buf.extend_from_slice(&[0x00, 0x00]);
                            }
                            buf.extend_from_slice(&[0x00, 0x00]);
                        }
                        buf.extend_from_slice(&[0x00, 0x00]);
                    }
                    buf.extend_from_slice(&[0x00, 0x00]);
                }
                buf.extend_from_slice(&[0x00, 0x00]);
            }
            buf.extend_from_slice(&[0x00, 0x00]);
        }
        buf.extend_from_slice(&[0x00, 0x00]);
    }
    buf.extend_from_slice(&[0x00, 0x00]);
    buf
}

fn encode_seqid_chain(buf: &mut Vec<u8>, first_token: &str) {
    let parts: Vec<&str> = first_token.split('|').collect();
    let mut i = 0usize;
    while i < parts.len() {
        match parts[i] {
            "lcl" if i + 1 < parts.len() => {
                encode_local_seqid(buf, parts[i + 1]);
                i += 2;
            }
            "gi" if i + 1 < parts.len() => {
                if let Ok(gi) = parts[i + 1].parse::<i32>() {
                    encode_gi_seqid(buf, gi);
                }
                i += 2;
            }
            "pdb" if i + 1 < parts.len() => {
                let chain = parts.get(i + 2).copied().filter(|chain| !chain.is_empty());
                encode_pdb_seqid(buf, parts[i + 1], chain);
                i += if chain.is_some() { 3 } else { 2 };
            }
            "gnl" if i + 2 < parts.len() => {
                encode_general_seqid(buf, parts[i + 1], parts[i + 2]);
                i += 3;
            }
            db if textseq_tag(db).is_some() => {
                let acc = parts.get(i + 1).copied().unwrap_or("");
                let name = parts.get(i + 2).copied().filter(|name| !name.is_empty());
                if !acc.is_empty() {
                    encode_textseq_seqid(buf, textseq_tag(db).unwrap(), Some(acc), name);
                } else if let Some(name) = name {
                    encode_textseq_seqid(buf, textseq_tag(db).unwrap(), None, Some(name));
                }
                i += if parts.get(i + 2).is_some() { 3 } else { 2 };
            }
            _ => {
                encode_textseq_seqid(buf, 0xa9, Some(first_token), None);
                return;
            }
        }
    }
}

fn textseq_tag(db: &str) -> Option<u8> {
    match db {
        "gb" => Some(0xa4),
        "emb" => Some(0xa5),
        "pir" => Some(0xa6),
        "sp" => Some(0xa7),
        "ref" => Some(0xa9),
        "dbj" => Some(0xac),
        "prf" => Some(0xad),
        "tpg" => Some(0xaf),
        "tpe" => Some(0xb0),
        "tpd" => Some(0xb1),
        _ => None,
    }
}

fn split_accession_version(value: &str) -> (&str, Option<i32>) {
    if let Some(dot) = value.rfind('.') {
        let version = &value[dot + 1..];
        if !version.is_empty() && version.bytes().all(|b| b.is_ascii_digit()) {
            if let Ok(version) = version.parse::<i32>() {
                return (&value[..dot], Some(version));
            }
        }
    }
    (value, None)
}

fn encode_textseq_seqid(buf: &mut Vec<u8>, tag: u8, accession: Option<&str>, name: Option<&str>) {
    let (accession, version) = accession
        .map(split_accession_version)
        .map(|(acc, ver)| (Some(acc), ver))
        .unwrap_or((None, None));
    buf.extend_from_slice(&[tag, 0x80, 0x30, 0x80]);
    if let Some(name) = name {
        encode_visible_string_field(buf, 0xa0, name);
    }
    if let Some(accession) = accession {
        encode_visible_string_field(buf, 0xa1, accession);
    }
    if let Some(version) = version {
        encode_integer_field(buf, 0xa3, version);
    }
    buf.extend_from_slice(&[0x00, 0x00]);
    buf.extend_from_slice(&[0x00, 0x00]);
}

fn encode_gi_seqid(buf: &mut Vec<u8>, gi: i32) {
    buf.extend_from_slice(&[0xab, 0x80, 0x02]);
    let gi_bytes = encode_asn1_integer(gi);
    encode_asn1_length(buf, gi_bytes.len());
    buf.extend_from_slice(&gi_bytes);
    buf.extend_from_slice(&[0x00, 0x00]);
}

fn encode_local_seqid(buf: &mut Vec<u8>, id: &str) {
    buf.extend_from_slice(&[0xa0, 0x80]);
    if let Ok(value) = id.parse::<i32>() {
        encode_integer_field(buf, 0xa0, value);
    } else {
        encode_visible_string_field(buf, 0xa1, id);
    }
    buf.extend_from_slice(&[0x00, 0x00]);
}

fn encode_pdb_seqid(buf: &mut Vec<u8>, mol: &str, chain: Option<&str>) {
    buf.extend_from_slice(&[0xae, 0x80, 0x30, 0x80]);
    encode_visible_string_field(buf, 0xa0, mol);
    if let Some(chain) = chain {
        encode_visible_string_field(buf, 0xa3, chain);
    }
    buf.extend_from_slice(&[0x00, 0x00]);
    buf.extend_from_slice(&[0x00, 0x00]);
}

fn encode_general_seqid(buf: &mut Vec<u8>, db: &str, tag: &str) {
    buf.extend_from_slice(&[0xaa, 0x80, 0x30, 0x80]);
    encode_visible_string_field(buf, 0xa0, db);
    buf.extend_from_slice(&[0xa1, 0x80]);
    if let Ok(id) = tag.parse::<i32>() {
        encode_integer_field(buf, 0xa0, id);
    } else {
        encode_visible_string_field(buf, 0xa1, tag);
    }
    buf.extend_from_slice(&[0x00, 0x00]);
    buf.extend_from_slice(&[0x00, 0x00]);
    buf.extend_from_slice(&[0x00, 0x00]);
}

fn encode_visible_string_field(buf: &mut Vec<u8>, tag: u8, value: &str) {
    buf.extend_from_slice(&[tag, 0x80, 0x1a]);
    encode_asn1_length(buf, value.len());
    buf.extend_from_slice(value.as_bytes());
    buf.extend_from_slice(&[0x00, 0x00]);
}

fn encode_integer_field(buf: &mut Vec<u8>, tag: u8, value: i32) {
    buf.extend_from_slice(&[tag, 0x80, 0x02]);
    let bytes = encode_asn1_integer(value);
    encode_asn1_length(buf, bytes.len());
    buf.extend_from_slice(&bytes);
    buf.extend_from_slice(&[0x00, 0x00]);
}

fn encode_asn1_length(buf: &mut Vec<u8>, len: usize) {
    if len < 128 {
        buf.push(len as u8);
    } else if len < 256 {
        buf.push(0x81);
        buf.push(len as u8);
    } else {
        buf.push(0x82);
        buf.push((len >> 8) as u8);
        buf.push(len as u8);
    }
}

fn encode_asn1_integer(val: i32) -> Vec<u8> {
    if val == 0 {
        vec![0]
    } else if val > 0 && val < 128 {
        vec![val as u8]
    } else if val > 0 && val < 32768 {
        if val < 256 {
            vec![0, val as u8]
        } else {
            vec![(val >> 8) as u8, val as u8]
        }
    } else {
        vec![
            (val >> 24) as u8,
            (val >> 16) as u8,
            (val >> 8) as u8,
            val as u8,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defline_asn1_includes_visible_title_accession_and_ord_id() {
        let encoded = encode_defline_asn1("seq1 sample title", 300);
        assert!(encoded.windows(17).any(|w| w == b"seq1 sample title"));
        assert!(encoded.windows(4).any(|w| w == b"seq1"));
        assert!(encoded.windows(9).any(|w| w == b"BL_ORD_ID"));
        assert!(encoded.windows(2).any(|w| w == [0x01, 0x2c]));
    }

    #[test]
    fn defline_asn1_parses_gi_refseq_chain() {
        let encoded = encode_defline_asn1("gi|12345|ref|NC_000001.11| chromosome", 7);
        assert_eq!(
            crate::db::index::extract_gi_from_header(&encoded),
            Some(12345)
        );
        assert_eq!(
            crate::db::index::extract_blast_seqid_chain_from_asn(&encoded).as_deref(),
            Some("gi|12345|ref|NC_000001.11|")
        );
        assert_eq!(
            crate::db::index::extract_textseq_bare_accession_from_asn(&encoded).as_deref(),
            Some("NC_000001")
        );
    }

    #[test]
    fn defline_asn1_parses_general_id_chain() {
        let encoded = encode_defline_asn1("gnl|TRACE_ASSM|ti222 trace read", 9);
        assert!(encoded.windows(10).any(|w| w == b"TRACE_ASSM"));
        assert!(encoded.windows(5).any(|w| w == b"ti222"));
    }

    #[test]
    fn defline_asn1_parses_local_id() {
        let encoded = encode_defline_asn1("lcl|local_subject local title", 11);
        assert!(encoded.windows(2).any(|w| w == [0xa0, 0x80]));
        assert!(encoded.windows(2).any(|w| w == [0xa1, 0x80]));
        assert!(encoded.windows(13).any(|w| w == b"local_subject"));
    }

    #[test]
    fn defline_asn1_parses_pdb_id() {
        let encoded = encode_defline_asn1("pdb|1ABC|A crystal structure", 12);
        assert!(encoded.windows(4).any(|w| w == [0xae, 0x80, 0x30, 0x80]));
        assert!(encoded.windows(4).any(|w| w == b"1ABC"));
        assert!(encoded.windows(1).any(|w| w == b"A"));
    }
}
