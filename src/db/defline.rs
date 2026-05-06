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
                if let Some(accession) = title.split_whitespace().next() {
                    buf.extend_from_slice(&[0xa9, 0x80]);
                    {
                        buf.extend_from_slice(&[0x30, 0x80]);
                        {
                            buf.extend_from_slice(&[0xa1, 0x80]);
                            buf.push(0x1a);
                            encode_asn1_length(&mut buf, accession.len());
                            buf.extend_from_slice(accession.as_bytes());
                            buf.extend_from_slice(&[0x00, 0x00]);
                        }
                        buf.extend_from_slice(&[0x00, 0x00]);
                    }
                    buf.extend_from_slice(&[0x00, 0x00]);
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
}
