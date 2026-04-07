use std::env;
use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let ncbi_root = manifest_dir
        .parent()
        .unwrap()
        .join("ncbi-blast-2.17.0+-src/c++");
    let core_src = ncbi_root.join("src/algo/blast/core");
    let comp_src = ncbi_root.join("src/algo/blast/composition_adjustment");
    let scoremat_src = ncbi_root.join("src/util/tables");
    let ncbi_include = ncbi_root.join("include");
    let stubs_dir = manifest_dir.join("stubs");

    // Collect all C source files from the blast core
    let core_c_files: Vec<PathBuf> = std::fs::read_dir(&core_src)
        .expect("Failed to read blast core source directory")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map_or(false, |ext| ext == "c"))
        .collect();

    // Collect composition adjustment C files
    let comp_c_files: Vec<PathBuf> = std::fs::read_dir(&comp_src)
        .expect("Failed to read composition_adjustment source directory")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map_or(false, |ext| ext == "c"))
        .collect();

    // Build the C library
    let mut build = cc::Build::new();
    build
        .warnings(false)
        // Stubs first so they override NCBI toolkit headers
        .include(&stubs_dir)
        // Then real NCBI headers for blast core and util/tables
        .include(&ncbi_include)
        // Define macros
        .define("NCBI_XBLAST_EXPORT", None)
        .define("NCBI_TABLES_EXPORT", None)
        .opt_level(3);

    // Disable OpenMP — threading will be handled by Rust/Rayon
    // build.flag_if_supported("-fopenmp");

    for f in &core_c_files {
        build.file(f);
    }
    for f in &comp_c_files {
        build.file(f);
    }
    // Score matrix data
    build.file(scoremat_src.join("raw_scoremat.c"));

    // MT_LOCK stubs
    build.file(stubs_dir.join("mt_lock_stub.c"));

    build.compile("blast_core");

    // No OpenMP linking needed
    println!("cargo:rustc-link-lib=m");

    // Generate bindings with bindgen
    let bindings = bindgen::Builder::default()
        .header(wrapper_header(&stubs_dir))
        .clang_arg(format!("-I{}", stubs_dir.display()))
        .clang_arg(format!("-I{}", ncbi_include.display()))
        .clang_arg("-DNCBI_XBLAST_EXPORT=")
        .clang_arg("-DNCBI_TABLES_EXPORT=")
        // Allowlist blast core API
        .allowlist_function("Blast.*")
        .allowlist_function("BLAST.*")
        .allowlist_function("BlastMem.*")
        .allowlist_function("ListNode.*")
        .allowlist_function("NCBISM_.*")
        .allowlist_function("PSI.*")
        .allowlist_function("Kappa.*")
        .allowlist_function("s_.*Blast.*")
        .allowlist_function("SBlast.*")
        .allowlist_function("LookupTable.*")
        .allowlist_function("BlastHSPCollector.*")
        .allowlist_function("BlastHSPBestHit.*")
        .allowlist_function("BlastHSPWriter.*")
        .allowlist_function("BlastSeqSrc.*")
        .allowlist_function("_BlastSeqSrcImpl_.*")
        .allowlist_function("BlastSeqBlk.*")
        .allowlist_function("BlastSequenceBlk.*")
        .allowlist_function("BlastSetUp_SeqBlkNew")
        .allowlist_function("GetProgram.*")
        .allowlist_function("NCBI_.*")
        .allowlist_function("NlmKarlin.*")
        .allowlist_function("Nlm_.*")
        // Allowlist blast types
        .allowlist_type("Blast.*")
        .allowlist_type("BLAST.*")
        .allowlist_type("SBlast.*")
        .allowlist_type("BlastSeqSrc.*")
        .allowlist_type("ListNode")
        .allowlist_type("SNCBIPackedScoreMatrix")
        .allowlist_type("SNCBIFullScoreMatrix")
        .allowlist_type("TNCBIScore")
        .allowlist_type("GapEditScript")
        .allowlist_type("GapEditBlock")
        .allowlist_type("GapState.*")
        .allowlist_type("SGreedy.*")
        .allowlist_type("Lookup.*")
        .allowlist_type("EBlast.*")
        .allowlist_type("Int[1248]")
        .allowlist_type("Uint[1248]")
        .allowlist_type("Boolean")
        .allowlist_type("SPHIPatternSearchBlk")
        // Allowlist constants
        .allowlist_var("NCBI_FSM_DIM")
        .allowlist_var("NCBISM_.*")
        .allowlist_var("BLASTNA_.*")
        .allowlist_var("BLAST_.*")
        .allowlist_var("kBlast.*")
        .allowlist_var("kDust.*")
        // Settings
        .derive_debug(true)
        .derive_default(true)
        .size_t_is_usize(true)
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

/// Create a wrapper header that includes the main blast core headers
fn wrapper_header(stubs_dir: &PathBuf) -> String {
    let header_path = stubs_dir.join("wrapper.h");
    // Always read from the checked-in wrapper.h
    if header_path.exists() {
        return header_path.to_string_lossy().to_string();
    }
    {
        std::fs::write(
            &header_path,
            r#"/* Auto-generated wrapper header for bindgen */
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_stat.h>
#include <algo/blast/core/blast_gapalign.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_traceback.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_seqsrc_impl.h>
#include <algo/blast/core/blast_hspstream.h>
#include <algo/blast/core/blast_encoding.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_parameters.h>
#include <algo/blast/core/blast_query_info.h>
#include <algo/blast/core/blast_program.h>
#include <algo/blast/core/blast_extend.h>
#include <algo/blast/core/blast_diagnostics.h>
#include <algo/blast/core/blast_lookup.h>
#include <algo/blast/core/blast_nalookup.h>
#include <algo/blast/core/blast_aalookup.h>
#include <algo/blast/core/blast_hspfilter.h>
#include <algo/blast/core/blast_psi.h>
#include <algo/blast/core/blast_kappa.h>
#include <algo/blast/core/blast_sw.h>
#include <algo/blast/core/blast_seg.h>
#include <algo/blast/core/blast_rps.h>
#include <algo/blast/core/gapinfo.h>
#include <algo/blast/core/gencode_singleton.h>
#include <algo/blast/core/ncbi_math.h>
#include <algo/blast/core/split_query.h>
#include <algo/blast/core/blast_tune.h>
#include <util/tables/raw_scoremat.h>
"#,
        )
        .expect("Failed to write wrapper.h");
    }
    header_path.to_string_lossy().to_string()
}
