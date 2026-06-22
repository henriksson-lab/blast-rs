use std::sync::Arc;

use crate::algo::blast::api::blast_options_cxx::EProgram;
use crate::algo::blast::api::blast_options_handle::{CBlastOptionsFactory, ETaskSet};
use crate::algo::blast::api::query_data::{
    CSeqId, CSeqLoc, CSeqLocInfo, EBlastProgramType, TMaskedQueryRegions,
};
use crate::algo::blast::core::blast_filter::BlastMaskLoc;
use crate::algo::blast::core::blast_util::blast_number_2_program;

#[derive(Clone, Debug)]
pub struct CFrameFinder {
    pub frame: i32,
}

#[derive(Clone, Debug)]
pub struct CAutomaticGenCodeSingleton {
    pub ref_counter: u32,
    pub genetic_codes: Vec<i32>,
}

#[derive(Clone, Debug)]
pub struct CBlastAppDiagHandler {
    pub saved_messages: Vec<String>,
    pub save: bool,
}

#[derive(Clone, Debug)]
pub struct CConstRefCSeqIdLessThan {
    pub lhs: Option<Arc<CSeqId>>,
    pub rhs: Option<Arc<CSeqId>>,
}

#[derive(Clone, Debug)]
pub struct BlastAuxConversionState {
    pub program: EBlastProgramType,
    pub seq_loc: Option<Arc<CSeqLoc>>,
    pub seq_ids: Vec<Arc<CSeqId>>,
    pub masked_regions: TMaskedQueryRegions,
    pub seq_loc_infos: Vec<Arc<CSeqLocInfo>>,
    pub mask: *const BlastMaskLoc,
}

pub fn throw_if_invalid_task(task: &str) {
    let valid_tasks = CBlastOptionsFactory::get_tasks(ETaskSet::All);
    if !valid_tasks.contains(task) {
        panic!("'{task}' is not a supported task");
    }
}

pub fn eprogram_to_task_name(p: EProgram) -> String {
    match p {
        EProgram::EBlastn => "blastn",
        EProgram::EMegablast => "megablast",
        EProgram::EDiscMegablast => "dc-megablast",
        EProgram::EBlastp => "blastp",
        EProgram::EBlastx => "blastx",
        EProgram::ETblastn => "tblastn",
        EProgram::ETblastx => "tblastx",
        EProgram::ERpsBlast => "rpsblast",
        EProgram::ERpsTblastn => "rpstblastn",
        EProgram::EPsiBlast => "psiblast",
        EProgram::EPsiTblastn => "psitblastn",
        EProgram::EPhiBlastp => "phiblastp",
        EProgram::EPhiBlastn => "phiblastn",
        EProgram::EDeltaBlast => "deltablast",
        EProgram::EVecScreen => "vecscreen",
        EProgram::EMapper => "mapr2g",
        _ => panic!("Invalid EProgram value: {p:?}"),
    }
    .to_string()
}

pub fn eprogram_to_eblast_program_type(p: EProgram) -> EBlastProgramType {
    match p {
        EProgram::EBlastn
        | EProgram::EMegablast
        | EProgram::EDiscMegablast
        | EProgram::EVecScreen => EBlastProgramType::Blastn,
        EProgram::EMapper => EBlastProgramType::Mapping,
        EProgram::EBlastp => EBlastProgramType::Blastp,
        EProgram::EBlastx => EBlastProgramType::Blastx,
        EProgram::ETblastn => EBlastProgramType::Tblastn,
        EProgram::ETblastx => EBlastProgramType::Tblastx,
        EProgram::ERpsBlast => EBlastProgramType::RpsBlast,
        EProgram::ERpsTblastn => EBlastProgramType::RpsTblastn,
        EProgram::EPsiBlast | EProgram::EDeltaBlast => EBlastProgramType::PsiBlast,
        EProgram::EPsiTblastn => EBlastProgramType::PsiTblastn,
        EProgram::EPhiBlastp => EBlastProgramType::PhiBlastp,
        EProgram::EPhiBlastn => EBlastProgramType::PhiBlastn,
        _ => EBlastProgramType::Undefined,
    }
}

pub fn blast_program_name_from_type(program: EBlastProgramType) -> String {
    let (status, program_string) = blast_number_2_program(program as u32, Some(()));
    if status == 0 {
        program_string.unwrap_or_default()
    } else {
        String::new()
    }
}

pub fn program_name_to_enum(program_name: &str) -> EProgram {
    assert!(!program_name.is_empty());

    let lowercase_program_name = program_name.to_ascii_lowercase();
    if lowercase_program_name.starts_with("blastn") {
        EProgram::EBlastn
    } else if lowercase_program_name.starts_with("rmblastn") {
        EProgram::EBlastn
    } else if lowercase_program_name.starts_with("blastp") {
        EProgram::EBlastp
    } else if lowercase_program_name == "blastx" {
        EProgram::EBlastx
    } else if lowercase_program_name == "tblastn" {
        EProgram::ETblastn
    } else if lowercase_program_name == "tblastx" {
        EProgram::ETblastx
    } else if lowercase_program_name == "rpsblast" {
        EProgram::ERpsBlast
    } else if lowercase_program_name == "rpstblastn" {
        EProgram::ERpsTblastn
    } else if lowercase_program_name == "megablast" {
        EProgram::EMegablast
    } else if lowercase_program_name == "psiblast" {
        EProgram::EPsiBlast
    } else if lowercase_program_name == "psitblastn" {
        EProgram::EPsiTblastn
    } else if lowercase_program_name == "dc-megablast" {
        EProgram::EDiscMegablast
    } else if lowercase_program_name == "deltablast" {
        EProgram::EDeltaBlast
    } else if lowercase_program_name == "vecscreen" {
        EProgram::EVecScreen
    } else if lowercase_program_name == "mapper" {
        EProgram::EMapper
    } else if lowercase_program_name == "mapr2g" {
        EProgram::EMapper
    } else if lowercase_program_name == "mapr2r" {
        EProgram::EMapper
    } else if lowercase_program_name == "mapg2g" {
        EProgram::EMapper
    } else {
        panic!("Program type '{program_name}' not supported");
    }
}
