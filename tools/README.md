# tools

Useful scripts for mCA analysis. The important pipeline tools are in: 
* `extract_cram_AD.sh`: extract raw and deduped ADs from either WGS or WES cram files at het sites
  * takes as input an ID, a region, and a bcf file containing het sites for extraction
* `compute_region_depth_PCadj.sh`: extract read depth for a set of regions in a set of individuals

## `RAP` 
Contains scripts used specifically to analyze data stored on UKB-RAP (not portable to other compute environments)
* `extract_DRAGEN_VCF.sh`: quickly extracts het genotypes from individual level VCFs
  * takes an ID as input and runs `../bin/dragen_vcf_extract_hetADs` on the individual level VCF
  * used in `../imputationQC/run_imputationQC.sh`, `../benchmark_mocha/run_benchmark_mocha.sh`, and `../prefilter/run_prefilter_pipeline.sh`
  * relies on precomputed GC-profiles and depth PCs (used in `../bin/computeRegionDepthsPCadj`)
  * can optionally output bin level region depths
  * can only extract depths for individuals in the same WGS_GROUP.batch
* `extract_cram_AD_RAP.sh`: wrapper for running `extract_cram_AD.sh` on RAP
* `compute_region_depth_PCadj_RAP.sh`: wrapper for running `compute_region_depth_PCadj.sh` on RAP