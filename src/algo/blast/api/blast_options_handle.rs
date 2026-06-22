use std::collections::BTreeSet;
use std::sync::Arc;

use crate::algo::blast::api::blast_options_cxx::{CBlastOptions, EApiLocality, EProgram};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ETaskSet {
    NuclNucl,
    ProtProt,
    Mapping,
    All,
}

#[derive(Clone, Debug)]
pub struct CBlastOptionsFactory {
    pub tasks: Vec<String>,
    pub documentation: Vec<(String, String)>,
    pub default_locality: EApiLocality,
}

#[derive(Clone, Debug)]
pub struct CBlastOptionsHandle {
    pub opts: Option<Arc<CBlastOptions>>,
    pub defaults_mode: bool,
    pub locality: EApiLocality,
    pub program: EProgram,
}

impl CBlastOptionsFactory {
    pub fn get_tasks(choice: ETaskSet) -> BTreeSet<String> {
        let mut retval = BTreeSet::new();
        if choice == ETaskSet::NuclNucl || choice == ETaskSet::All {
            retval.insert("blastn".to_string());
            retval.insert("blastn-short".to_string());
            retval.insert("megablast".to_string());
            retval.insert("dc-megablast".to_string());
            retval.insert("vecscreen".to_string());
            retval.insert("rmblastn".to_string());
        }

        if choice == ETaskSet::ProtProt || choice == ETaskSet::All {
            retval.insert("blastp".to_string());
            retval.insert("blastp-short".to_string());
            retval.insert("blastp-fast".to_string());
        }

        if choice == ETaskSet::All {
            retval.insert("psiblast".to_string());
            retval.insert("phiblastp".to_string());
            retval.insert("rpsblast".to_string());
            retval.insert("rpstblastn".to_string());
            retval.insert("blastx".to_string());
            retval.insert("blastx-fast".to_string());
            retval.insert("deltablast".to_string());
            retval.insert("tblastn".to_string());
            retval.insert("tblastn-fast".to_string());
            retval.insert("psitblastn".to_string());
            retval.insert("tblastx".to_string());
            retval.insert("kblastp".to_string());
        }

        if choice == ETaskSet::Mapping || choice == ETaskSet::All {
            retval.insert("mapper".to_string());
            retval.insert("mapr2g".to_string());
            retval.insert("mapr2r".to_string());
            retval.insert("mapg2g".to_string());
        }

        retval
    }
}
