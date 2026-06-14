pub struct CRpsAuxFile;
pub struct CRpsPssmFile;
pub struct CRpsLookupTblFile;
pub struct CRpsFreqsFile;
pub struct CRpsObsrFile;
pub struct CRpsFreqRatiosFile;
pub struct BlastRPSInfo;

pub enum CBlastRPSInfoOpenFlag {
    LookupTableFile = 1,
    PssmFile = 2,
    AuxInfoFile = 4,
    FrequenciesFile = 8,
    ObservationsFile = 16,
    FreqRatiosFile = 32,
    RpsBlast = 7,
    DeltaBlast = 24,
    RpsBlastWithCbs = 39,
}

pub struct CBlastRPSInfo {
    pub aux_file: Option<Box<CRpsAuxFile>>,
    pub pssm_file: Option<Box<CRpsPssmFile>>,
    pub lut_file: Option<Box<CRpsLookupTblFile>>,
    pub freqs_file: Option<Box<CRpsFreqsFile>>,
    pub obsr_file: Option<Box<CRpsObsrFile>>,
    pub freq_ratios_file: Option<Box<CRpsFreqRatiosFile>>,
    pub rps_info: Option<Box<BlastRPSInfo>>,
}
