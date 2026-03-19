# src

Contains source code for tools used to analyze mCAs. Most tools can be accessed through the `mCAs_WGS` binary.

* `prepare_bcf.cpp`, `prepare_bcf.h`, `run_prepare_bcf.cpp`: tools for efficiently masking CNVs, annotating in allelic depths, calculating/loading ref bias, performing genotype QC by imputation, etc.
* `hmm.cpp`, `hmm.h`, `make_calls.cpp`, `genome_rules.c`, `genome_rules.h`: reimplementation of core MoChA functionality
* `gwas_tools.cpp`, `gwas_tools.h`, `write_CNLOH_GRM.cpp`, `compute_polygenic_drive.cpp`, `compute_cis_fisher_p.cpp`: tools that are helpful for performing GWAS on cis CN-LOH phenotype
* `vcf_extract_common_AC_AN.cpp`, `dragen_vcf_extract_hetADs.cpp`: tools for quickly unpacking WGS vcf files

## `read_depth`

Files in `read_depth` contain source code relevant to the depth profiling pipeline. The important files are summarized below:

* `GenerateRefMasks.cpp`: generates binary GRCh38 masks used in the depth computation
* `CounteReadsSaveDiscordant.cpp`: takes a CRAM file and saves bin level read counts (1kb bins) and all discordant reads in the CRAM file
* `ComputeDepthProfiles.cpp`: creates 20kb and 50kb GC corrected depth profiles and calls germline CNVs to be masked in mCA analysis
* `ComputeRegionDepthsPCadj.cpp`: computes PC adjusted depth profile within the given query regions 