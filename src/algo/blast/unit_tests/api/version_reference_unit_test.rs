// Upstream source:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp

#[cfg(test)]
mod tests {
    use crate::algo::blast::api::version::{
        BlastVersion, Publication, Reference, K_BLAST_MAJOR_VERSION, K_BLAST_MINOR_VERSION,
        K_BLAST_PATCH_VERSION,
    };

    fn do_test_reference(keyword: &str, pubmed_id: &str, publication: Publication) {
        let reference = Reference::get_string(publication);
        assert!(reference.contains(keyword));

        let url = Reference::get_pubmed_url(publication);
        assert!(url.contains(pubmed_id));
    }
}
