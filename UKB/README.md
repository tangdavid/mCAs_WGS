# Analysis of UK Biobank

These folders contain all scripts used to analyze UK Biobank.

The pipeline we used to call mCAs in UK Biobank followed this sequence:
`read_depth_profiling` → `generate_impute_ref` → `imputationQC` → `cramADs` → `prefilter` → `make_calls`

The MoChA benchmarks are located in `benchmark_mocha`.

Post-calling analyses (roughly in the order they appear in the manuscript) are organized across the following directories:
- `breakpoints` - mCA hotspot analysis
- `focal` - 13q deletion analysis  
- `cancer` - extracting cancer phenotypes
- `transGWAS` - trans-GWAS analysis
- `cisGWAS` - cis-GWAS and burden testing
- `polygenic_drive` - polygenic score analysis
- `plots` - figure generation notebooks


## Table of Contents

- [Analysis of UK Biobank](#analysis-of-uk-biobank)
  - [Table of Contents](#table-of-contents)
  - [Practical note about code portability](#practical-note-about-code-portability)
    - [A note about sample batching in UKB](#a-note-about-sample-batching-in-ukb)
  - [Software dependencies](#software-dependencies)
  - [Depth pipeline](#depth-pipeline)
    - [`read_depth_profiling`](#read_depth_profiling)
  - [Generating calls](#generating-calls)
    - [`generate_impute_ref`](#generate_impute_ref)
    - [`imputationQC`](#imputationqc)
    - [`cramADs`](#cramads)
    - [`prefilter`](#prefilter)
    - [`make_calls`](#make_calls)
    - [`plots`](#plots)
  - [MoChA benchmarking](#mocha-benchmarking)
    - [`benchmark_mocha`](#benchmark_mocha)
  - [Analysis of mCA hotspots](#analysis-of-mca-hotspots)
    - [`breakpoints`](#breakpoints)
    - [`plots`](#plots-1)
  - [Analysis of 13q deletions](#analysis-of-13q-deletions)
    - [`focal`](#focal)
    - [`transGWAS`](#transgwas)
    - [`plots`](#plots-2)
  - [Burden test for cis-acting CN-LOH](#burden-test-for-cis-acting-cn-loh)
    - [`cisGWAS`](#cisgwas)
    - [`plots`](#plots-3)
  - [Common variant associations with CN-LOH](#common-variant-associations-with-cn-loh)
    - [`cisGWAS`](#cisgwas-1)
    - [`polygenic_drive`](#polygenic_drive)
    - [`plots`](#plots-4)

## Practical note about code portability

The scripts used to analyze UK Biobank are **not readily portable** to different compute environments. For instance, we have hard coded in absolute paths specific to our RAP workspace for convenience when running the pipeline. We have provided the code as is for transparency and completeness, providing a reference for other analyses of similarly large cohorts.

### A note about sample batching in UKB

Most of our scripts used to analyze UKB used an output naming convention corresponding to our sample batching scheme which we spell out here for clarity. 

We analyzed UKB in batches of 1000 samples for computational efficiency. Each samples was assigned to a super batch (`WGS_GROUP`) and a sub batch within the super batch (`BATCH`). Many output files end up with the prefix `$WGS_GROUP.batch$BATCH`. There are 6 different WGS_GROUPs (`WGS1`,...,`WGS6`) corresponding to different sequencing centers and rounds of sequencing.

| Variable         | Description                                                                                   |
|------------------|-----------------------------------------------------------------------------------------------|
| `WGS_GROUP`      | Identifier for one of six major groups of WGS samples (e.g., WGS1, WGS2, etc.)                 |
| `BATCH`          | Sub-batch within a WGS group                                              |
| `ID`/`SAMPLE`             | Sample identifier                                                                             |
| `OUT_DIR`        | Output directory for results                                                                  |
| `PFX`            | Prefix for output files, often combining WGS_GROUP and BATCH (e.g. `WGS2.batch58`)                                 |


## Software dependencies
Our analysis of UKB used the following software tools:
```
samtools v1.15.1
bcftools v1.15.1
plink v1.9
plink v2.0
IMPUTE5 v2.0.0
XCFtools v5.0.0
SHAPEIT5 v5.1.1
GCTA v1.94
REGENIE v3.1.1
python 3.8.10
```

## Depth pipeline

### `read_depth_profiling`
* `run_cram_depth_pipeline_RAP.sh`: extracts read depth profiles from UK Biobank WGS CRAM files for use in the mCA calling pipeline
  * processes CRAM files in batches organized by WGS_GROUP and BATCH identifiers
  * uses `../../bin/countReadsSaveDiscordant` to count reads and save discordant reads from CRAM files
  * creates GC-adjusted depth profiles at multiple resolutions (20kb and 50kb)
  * generates outputs used by the main mCA calling pipeline

## Generating calls

### `generate_impute_ref`
* `add_AC_snp_scaffold_RAP_launch.sh`: add AC and AN INFO fields to SNP-array scaffold
* `run_generate_reference_RAP.sh`: subset SHAPEIT reference panel by deCODE vs Sanger and create two xcf reference panels for IMPUTE5 input
  * uses `../../bin/vcf_extract_common_AC_AN` to quickly get common variants from SHAPEIT VCF

### `imputationQC`
* `run_imputeQC.sh`: impute genotypes at all common variants using the deCODE and Sanger reference panels. Set imputed genotypes that disagree with WGS genotypes (from DRAGEN individual level vcfs) to missing.
  * depends on `../../tools/RAP/extract_DRAGEN_VCF.sh` to quickly extract het genotypes from individual level VCFs
  * runs imputation on all chromosomes by calling `impute_phase_chr.sh` for a given batch of ~1000 individuals

### `cramADs`
* `run_extractADs_chrom.sh`, `run_extractADs_chrom_launch.sh`: extract deduped ADs for all chromosomes in a given batch
  * used for benchmarking the prefiltering step and for computing ref bias on a subset of 10k individuals
  * relies on `../../tools/RAP/extract_cram_AD_RAP.sh`
  * works for both WGS and WES
* `run_compute_ref_bias.sh`, `run_compute_ref_bias_launch.sh`: compute empirical ref bias using the 10k individuals with full AD extraction
  * masks out CNV in ref bias computation
  * uses the `../../bin/mCAs_WGS` binary
  * works for both WGS and WES
* `collate_extractedADs.sh`, `merge_ref_bias_files.sh`: scripts to merge the outputs across chromosomes and batches

### `prefilter`
* `run_prefilter_pipeline.sh`, `run_prefilter_pipeline_launch.sh`: perform nomination of candidate chromosomes to follow up on using the uncorrected allelic depths in the individual level VCF files
  * uses `MoChA` for prefiltering
  * saves a tarball with bcfs containing allelic depths for only the nominated chromosomes
  * each batch takes ~10 hours 

### `make_calls`
* `run_make_calls.sh`, `run_make_calls_launch.sh`: runs a HMM on the deduped allelic depths contained in the prefiltered vcf files
  * each batch only takes ~1 hour
  * calls `../../tools/RAP/compute_region_depth_PCadj_RAP.sh` to extract depth information
  * calls `call_isochromosomes.sh` to postprocess isochromosomes called as CN-LOH
  * calls `validate_calls_WES.sh` to validate mCAs with WES allelic imbalance
* `call_isochromosomes.sh`: reassigns erroneously called isodisomies as isochromosomes by checking depth in the p- and q- arms
* `validate_calls_WES.sh`: checks concordance of allelic imbalance in WES and WGS
  * calls `compute_WES_AI.sh` for all (ID,mCA) pairs in the call set
  * calls `compute_expected_replication.py` to compute the expected validation rate
* `run_all_calls.sh`, `run_all_calls_launch.sh`: perform full analysis on the subset of 10k individuals with genome-wide deduped ADs
* `filter_mCA_calls.py`: applies minimal post-hoc filters to clean up the call set
* `merge_AS_bcf.sh`, `collate_calls.sh`: scripts to merge outputs across chromosomes and samples

### `plots`
`fig1_panels.ipynb` contains all code used to generate call set overview panels in notebook format

## MoChA benchmarking

### `benchmark_mocha`
* `run_benchmark_mocha.sh`, `run_benchmark_mocha_launch.sh`: benchmarks the mCA calling pipeline against MoChA (MOsaic CHromosomal Aberrations) method
  * processes the same samples using both the custom HMM pipeline and the established MoChA approach
  * uses imputed common variant genotypes for comparison
  * relies on `../../tools/RAP/extract_DRAGEN_VCF.sh` to extract genotypes from individual level VCFs
  * installs and runs MoChA software for comparison calling
  * generates comparative call sets for validation and method assessment
* `compare_calls.sh`: compares mCA calls between the custom pipeline and MoChA method
  * evaluates concordance and discordance between calling approaches
  * generates metrics for method validation and performance assessment

## Analysis of mCA hotspots

### `breakpoints`
* `extract_discordant_reads.sh`, `run_extract_discordant_reads_from_calls.sh`, `run_extract_discordant_reads_from_calls_launch.sh`: find all discordant reads where both read pairs map within 50kb of the boundaries of an mCA that does not extend to the telomere
  * calls `reconstruct_complex_sv.sh` to match discordant reads to mCA breakpoints
  * bins discordant reads into 1kb bins by the left and right breakpoints to ensure consistency of read orientation
* `get_discordant_breakpoints.sh`: takes in the set of extracted discordant reads and assigns breakpoints to simple deletions and duplications
  * checks orientation to make sure there is agreement between read depth and breakpoint
  * assigns the closest discordant pair to the mCA in the rare case where there are two discordant pairs that map within 50kb of both breakpoints
  * calls `merge_adjacent_bins.sh` to handle cases where breakpoints lie at boundaries of 1kb discordant read bins
* `get_split_breakpoints.sh`: finds exact breakpoints for mCAs with breakpoints localized by discordant reads
  * calls `get_split_reads.sh` to parse cram files for split reads based on supplemental alignment
  * calls `compute_microhomology.sh` to compute microhomology between the left and right alignment of split reads by sequence lookup in the GRCh38 reference

### `plots`
`focal_regions.py` contains code for defining focal regions, computing focal indices, and assigning genes to hotspots. `fig2_panels.ipynb` contains all code used to generate plots

## Analysis of 13q deletions

### `focal`
* `extract_focal_region_depth.sh`, `run_extract_focal_region_depth.sh`, `run_extract_focal_region_depth_launch.sh`: extract depth from a region defined in the `focal_region_list.txt` file (used exclusively for 13q14 depth in the paper analyses)
* `find_focal_breakpoint.sh`, `run_find_focal_breakpoint.sh`, `run_find_focal_breakpoint_launch.sh`: scripts to extract depth within the boundaries of discordant reads
  * extracts 1kb depth profiles within the larger focal region and calls `get_discordant_read_depth.sh` to summarize depth within each of the discordant reads found in the region
  * calls `bin_discordant.sh` to process raw discordant reads extracted from CRAM files
* `run_extract_focal_region_depth.sh`, `resolve_overlapping_discordant.py`: perform variable selection on the discordant read depths in `find_focal_breakpoint.sh` to prevent calling of spurious breakpoints

### `transGWAS`
* `make_13q_phenotype.py`: create file with mCA+CLL phenotypes for GWAS with regenie
  * the 13q phenotype contains a union of calls made from the BAF pipeline (HMM) and the depth pipeline
  * appends CLL read from `stdin`
* `regenie_CLL_step1.sh`, `regenie_CLL_step2.sh`, `regenie_CLL_step1_launch.sh`, `regenie_CLL_step2_launch.sh`: runs regenie steps 1 and 2 on the phenotypes generated in `make_13q_phenotype.py`

### `plots`
`fig3_panels.ipynb` contains all code used to generate plots. This notebook also contains the entirety of the survival analysis code. The first few cells of this notebook also contains the code used to recalibrate z-scores and generate a 13q deletion call set based on depth. 

## Burden test for cis-acting CN-LOH

### `cisGWAS`
* `run_extract_gnomAD_variants.sh`: extracts VEP annotations from gnomAD WES vcfs
  * run on O2 compute cluster rather than RAP
  * relies on PrimateAI-3D scores for missense variant annotation
  * calls `parse_gnomAD_vep` (source code in `parse_gnomAD_vep.cpp`) to keep only moderate and high confidence VEP protein coding variants
  * launched by `run_extract_gnomAD_variants_launch.sh` (RAP) and `run_extract_gnomAD_variants_sbatch.sh` (O2)
* `run_regenie_write_mask.sh`, `run_regenie_write_mask_launch.sh`: writes bgen and vcf burden masks using regenie
  * relies on `extractBgenVariants` to efficiently find protein coding variants from WGS_500k bgen files
  * merges SNP and indels with CNV calls (from HI-CNV WES)
  * writes a vcf with >=20 MAC protein coding variants, burden masks, and SNP-array scaffold for SHAPEIT5 phasing
  * uses `make_burden_masks.sh` to define burden mask annotations
* `make_burden_masks.sh`: defines burden masks used by regenie
  * only uses LoF and missense variants with PrimateAI-3D score > 0.6
  * performs QC based on UKB AF consistency with gnomAD non_ukb_nfe AF  (drop UKB variants that are ten times more common than in gnomAD)
* `run_phase_rare.sh`, `run_phase_rare_launch.sh`: phases the burden masks generated by `run_regenie_write_mask.sh`
  * uses `phase_rare_static` (v5.1.1)
* `annotate_AS.sh`, `run_annotate_AS_launch.sh`: concatenates the chunks from `run_phase_rare.sh` and annotates in allelic shift information for chromosomes affected by CN-LOH
  * relies on AS annotations from the mCA call bcfs
  * will only annotate AS for the scaffold variants (which can then be extended to the burden masks when running the GWAS)
* `run_burden_GWAS.sh`, `run_burden_GWAS_launch.sh`: runs the Fisher and Binomial test for cis associations with CN-LOH using the phased burden masks
* `apply_bh_FDR.sh`: applies BH procedure to control FDR in GWAS
* `run_refine_phase.sh`, `run_refine_phase_launch.sh`: run read based phasing on all `LoF.all_transcripts.missense8.0.001` masks
  * calls `refine_phase.sh`
  * summarizes the read based phasing at the gene level (number of variants made hom by CN-LOH, number of variants removed by CN-LOH, and number of variants that cannot be confidently phased)
  * (additional script to run read based phasing on FDR significant masks with the strongest association in `run_refine_phase_FDR_sig.sh`)
* `refine_phase.sh`: performs read based phasing for a given gene and a given mask
  * outputs phase of variants relative to CN-LOH phase (only considers variants that are overlapped by a CN-LOH)
  * relies on `extract_LoF_variant.sh` to find the specific LoF that a burden mask carrier has
  * performs read based phasing with `read_based_phasing.sh` (and optionally `run_whatshap_1.sh`)
  * outputs the read-based and statistical phase of every variant
* `extract_LoF_variant.sh`: extracts the specific LoF that a burden mask carrier has
* (helper tools) `extract_VAF.sh`, `extract_VAF_mask.sh`, and `extract_VAF_sample.sh`: extracts VAF of variants in burden masks
* `process_depmap_data.sh`: small script to download DepMap data for the genome wide phase analysis

### `plots`
`fig4_panels.ipynb` contains code for generating pathway enrichment, penetrance, and fraction of CN-LOH explained by rare variant plots. `fig4_panels_circos.ipynb` contains code to generate the ideogram and was run locally due to package install issues

## Common variant associations with CN-LOH

### `cisGWAS`
* `run_common_var_GWAS.sh`, `run_common_var_GWAS_launch.sh`: scripts to perform common variant GWAS of CN-LOH
  * performs both the Fisher's test and the binomial test for allelic shift
  * performs analyses in batches and then meta-analyzes the results across all batches to avoid having to form a large bcf with all common variants
  * tests variants included from within cohort imputation

### `polygenic_drive`
* `generate_snp_array_bcf.sh`, `run_generate_snp_array_bcf.sh`: annotates AS tag from mCA calls into the SNP-array scaffold which is then used by `run_compute_polygenic_drive.sh`
* `run_compute_polygenic_drive.sh`: computes the difference in PRS for various blood traits between arms of chromosomes duplicated and removed by CN-LOH
* `estimate_common_var_h2.sh`, `estimate_common_var_h2_launch.sh`: estimates common variant heritability of CN-LOH shift direction using `gcta --reml`
* `annotate_blood_sumstats.sh`: merges blood GWAS summary stats with CN-LOH binomial test summary stats
* `run_h2_simulations.sh`, `simulate_h2_estimation.py`: runs simulations to evaluate the calibration of CN-LOH heritability estimation

### `plots`
`fig5_panels.ipynb` contains code to generate Manhattan plot and polygenic score heatmaps
