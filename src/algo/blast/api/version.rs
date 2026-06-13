//! Rust port of NCBI `version.cpp` and `version.hpp`.

pub const K_BLAST_MAJOR_VERSION: i32 = 2;
pub const K_BLAST_MINOR_VERSION: i32 = 17;
pub const K_BLAST_PATCH_VERSION: i32 = 0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastVersion {
    pub major: i32,
    pub minor: i32,
    pub patch: i32,
}

impl BlastVersion {
    pub const fn new() -> Self {
        Self {
            major: K_BLAST_MAJOR_VERSION,
            minor: K_BLAST_MINOR_VERSION,
            patch: K_BLAST_PATCH_VERSION,
        }
    }

    pub fn print(&self) -> String {
        format!("{}.{}.{}+", self.major, self.minor, self.patch)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Publication {
    GappedBlast = 0,
    PhiBlast,
    MegaBlast,
    CompBasedStats,
    CompAdjustedMatrices,
    IndexedMegablast,
    DeltaBlast,
    MaxPublications,
}

const REFERENCES: [&str; Publication::MaxPublications as usize + 1] = [
    "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of protein database search programs\", Nucleic Acids Res. 25:3389-3402.",
    "Zheng Zhang, Alejandro A. Sch&auml;ffer, Webb Miller, Thomas L. Madden, David J. Lipman, Eugene V. Koonin, and Stephen F. Altschul (1998), \"Protein sequence similarity searches using patterns as seeds\", Nucleic Acids Res. 26:3986-3990.",
    "Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), \"A greedy algorithm for aligning DNA sequences\", J Comput Biol 2000; 7(1-2):203-14.",
    "Alejandro A. Sch&auml;ffer, L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001), \"Improving the accuracy of PSI-BLAST protein database searches with composition-based statistics and other refinements\", Nucleic Acids Res. 29:2994-3005.",
    "Stephen F. Altschul, John C. Wootton, E. Michael Gertz, Richa Agarwala, Aleksandr Morgulis, Alejandro A. Sch&auml;ffer, and Yi-Kuo Yu (2005) \"Protein database searches using compositionally adjusted substitution matrices\", FEBS J. 272:5101-5109.",
    "Aleksandr Morgulis, George Coulouris, Yan Raytselis, Thomas L. Madden, Richa Agarwala, Alejandro A. Sch&auml;ffer (2008), \"Database Indexing for Production MegaBLAST Searches\", Bioinformatics 24:1757-1764.",
    "Grzegorz M. Boratyn, Alejandro A. Schaffer, Richa Agarwala, Stephen F. Altschul, David J. Lipman and Thomas L. Madden (2012) \"Domain enhanced lookup time accelerated BLAST\", Biology Direct 7:12.",
    "",
];

const PUBMED_URLS: [&str; Publication::MaxPublications as usize + 1] = [
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9254694&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9705509&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=10890397&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=11452024&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=16218944&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/pubmed/18567917",
    "//www.ncbi.nlm.nih.gov/pubmed/22510480",
    "",
];

pub struct Reference;

impl Reference {
    pub fn get_string(pub_: Publication) -> &'static str {
        REFERENCES[pub_ as usize]
    }

    pub fn get_html_free_string(pub_: Publication) -> String {
        Self::get_string(pub_).replacen("&auml;", "a", 1)
    }

    pub fn get_pubmed_url(pub_: Publication) -> &'static str {
        PUBMED_URLS[pub_ as usize]
    }
}
