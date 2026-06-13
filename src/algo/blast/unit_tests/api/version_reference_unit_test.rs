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

    #[test]
    fn test_version() {
        let version = BlastVersion::new();
        assert_eq!(K_BLAST_MAJOR_VERSION, version.major);
        assert_eq!(K_BLAST_MINOR_VERSION, version.minor);
        assert_eq!(K_BLAST_PATCH_VERSION, version.patch);
        assert_eq!(
            format!(
                "{}.{}.{}+",
                K_BLAST_MAJOR_VERSION, K_BLAST_MINOR_VERSION, K_BLAST_PATCH_VERSION
            ),
            version.print()
        );
    }

    #[test]
    fn test_gapped_blast_reference() {
        do_test_reference(
            "Gapped BLAST and PSI-BLAST",
            "9254694",
            Publication::GappedBlast,
        );
    }

    #[test]
    fn test_mega_blast_reference() {
        do_test_reference(
            "A greedy algorithm for aligning DNA",
            "10890397",
            Publication::MegaBlast,
        );
    }

    #[test]
    fn test_composition_based_statistics_reference() {
        do_test_reference(
            "with composition-based statistics",
            "11452024",
            Publication::CompBasedStats,
        );
    }

    #[test]
    fn test_phi_blast_reference() {
        do_test_reference("using patterns as seeds", "9705509", Publication::PhiBlast);
    }

    #[test]
    fn test_invalid_reference() {
        assert!(Reference::get_string(Publication::MaxPublications).is_empty());
        assert!(Reference::get_pubmed_url(Publication::MaxPublications).is_empty());
    }

    #[test]
    fn print_all_references_html() {
        let publications = [
            Publication::GappedBlast,
            Publication::PhiBlast,
            Publication::MegaBlast,
            Publication::CompBasedStats,
            Publication::CompAdjustedMatrices,
            Publication::IndexedMegablast,
            Publication::DeltaBlast,
        ];

        let mut text = String::new();
        for publication in publications {
            text.push_str(Reference::get_string(publication));
            text.push('\n');
            text.push_str("<a href=");
            text.push_str(Reference::get_pubmed_url(publication));
            text.push_str(">Get abstract</a><br>\n");
        }

        assert_eq!(publications.len(), text.matches("www").count());
    }
}
