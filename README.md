# mCAs_WGS

A pipeline for detecting mosaic chromosomal alterations (mCAs) from whole genome sequencing (WGS) data. This repository includes scripts, binaries, and workflows for processing large-scale sequencing datasets to identify and analyze mCAs.

## Repository Structure
- `bin/` — Compiled binaries and tools
    - Note: this repository ships with some precompiled binaries in `bin/` for convenience; if a binary doesn't run on your system you can compile the project from source (see the "Build & software dependencies" section).
- `src/` — Source code (C/C++ and scripts)
- `pipeline/` — Pipeline scripts for running analyses
- `tools/` — Helper scripts for data processing
- `UKB/` — Custom scripts and project-specific workflows used for the manuscript. These are not easily portable and are included primarily for transparency and reproducibility.

## Getting Started
1. Clone the repository:
   ```sh
   git clone https://github.com/tangdavid/mCAs_WGS.git
   cd mCAs_WGS
   ```
2. Run the pipeline scripts to process your data. For example:
   ```sh
   # Step 1: Compute read-depth profiles from CRAM files
   bash pipeline/run_cram_depth_pipeline.sh -m crams.manifest -o depth_out -n cohortA

   # Step 2: Call mCAs on a chromosome
   bash pipeline/run_call_mCAs_pipeline.sh \
     -m crams.manifest \
     -o results \
     -c chr1 \
     --depth-dir depth_out \
     --ref-bias ref_bias.txt \
     --genome-blacklist blacklist.bed \
     --pc-dir pc_dir \
     --ref-file ref.fa \
     --phased-bcf phased.bcf
   ```
   Adjust the arguments and file paths as needed for your analysis.


### Quick check: prebuilt binaries and recompiling
If a binary doesn't run on your system (`bin/mCAs_WGS`, `bin/countReadsSaveDiscordant`, `bin/computeDepthProfiles`, or `bin/computeRegionDepthsPCadj`), recompile from source using the Makefile. Detailed build instructions are provided below: see [Software dependencies and build](#software-dependencies-and-build).

Common reasons prebuilt binaries fail: incompatible glibc or CPU architecture, missing dynamic libraries, or insufficient permissions. Recompiling on your system will usually resolve these.

## Example with 1000G
1. Ensure that `parallel`, `pigz`, and `awscli` are installed 
    ```sh
    sudo apt --yes install parallel
    sudo apt --yes install pigz
    sudo apt --yes install awscli
    ```
2. Download resources 
    ```sh
    OUT_DIR=/path/to/output
    wget -r -np -nH --cut-dirs=1 -R "index.html*" -P $OUT_DIR https://data.broadinstitute.org/lohlab/mCAs_WGS/
    ```
3. Create the 1000G CRAM manifest 
    ```sh
    
    aws s3 ls --no-sign-request s3://1000genomes-dragen-v4.0.3/data/1kgp_cram/ --recursive \
        | egrep cram$ \
        | awk '
    function basename(file) {
        sub(".*/", "", file)
        split(file, a,".")
        return a[1]
    }
    BEGIN {
        pfx="s3://1000genomes-dragen-v4.0.3"
    }
    {
        print basename($4),pfx"/"$4
    } > $OUT_DIR/1000G.cram_manifest.txt
    ```
4. Generate depth profiles for 1000G individuals
    ```sh
    pipeline/run_cram_depth_pipeline.sh -m $OUT_DIR/1000G.cram_manifest.txt -o $OUT_DIR -n 1000G
    ```
5. Call mCAs
    ```sh
    CHR=... # chose a chromosome to analyze
    REF_FILE=/path/to/GRCh38/fasta # can be downloaded from https://software.broadinstitute.org/software/mocha/
    PHASED_BCF=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

    pipeline/run_call_mCAs_pipeline.sh \
        -m $OUT_DIR/1000G.cram_manifest.txt \
        -o $OUT_DIR \
        -c $CHR \
        --depth-dir $OUT_DIR \
        --ref-bias $OUT_DIR/mCAs_WGS/10k.WGS.ref_bias.txt.gz \
        --genome-blacklist $OUT_DIR/mCAs_WGS/lc_sv_cnv_masks.hg38.bed \
        --pc-dir $OUT_DIR/mCAs_WGS \
        --ref-file $REF_FILE \
        --phased-bcf $PHASED_BCF 
    ```




The depth profiling for 1000G takes ~30 hours when run with 16 cores and the mCA calling takes between 3 and 14 hours when parallelized across chromosomes (using 16 cores for each chromosome).

## License
This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

## Software dependencies and build

The most important libraries required to build the C++ binaries in this repo are:

- HTSlib (libhts) — required for reading CRAM/BCF/VCF/SAM and linked by the C++ tools in `src/`.
- Boost (development headers and libraries) — at least `boost_iostreams` and `boost_program_options` are used/linked by the Makefile.
- libdeflate (development package) — used by HTSlib for high-performance DEFLATE compression; the top-level Makefile links `-ldeflate` if available.

Other required development libraries (for full CRAM support and common features):

- zlib development headers (e.g. `zlib1g-dev`)
- libbz2-dev and liblzma-dev (recommended)
- libcurl development headers (optional but recommended if you need HTTP/S/FTP/S3/GCS access in HTSlib)

Tools used by pipeline scripts (not required to build the C++ code but needed for running pipelines):

- parallel
- pigz
- awscli (only required if you use S3 URLs)

Example (Debian/Ubuntu) packages to install common prerequisites:

```sh
sudo apt update
sudo apt install -y build-essential autoconf automake pkg-config git \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev \
    libboost-iostreams-dev libboost-program-options-dev \
    parallel pigz awscli
# optional but recommended for performance
sudo apt install -y libdeflate-dev
```

Notes about library locations
- You must install HTSlib and Boost on your system, or place them somewhere locally and point the Makefile at them.
- The top-level `Makefile` accepts a single `HTSLD` prefix (recommended) or the pair `HTSLD_INC`/`HTSLD_LIB` to locate HTSlib.

Important: set `HTSLD` (or `HTSLD_INC`/`HTSLD_LIB`) when running `make`.

Examples:

```sh
# Use a single HTSlib prefix (recommended)
HTSLD=/opt/htslib make

# Or set include/lib explicitly
HTSLD_INC=/usr/local/include HTSLD_LIB=/usr/local/lib make
```

- If you installed HTSlib system-wide under a standard prefix (e.g. `/usr/local`), you can also point `HTSLD` to that prefix or omit it after running `sudo make install` for HTSlib.

Building HTSlib (included copy or from source)

You can either install HTSlib system-wide or build it locally and point the top-level Makefile at the local tree. The following builds HTSlib locally. 

```sh
sudo apt --yes install libdeflate-dev # (to get libdeflate)
wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
bunzip2 htslib-1.18.tar.bz2
cd htslib-1.18
./configure --with-libdeflate
make
```

Troubleshooting

- If the linker complains about `-lhts` or `-ldeflate`, set `HTSLD_INC` and `HTSLD_LIB` in the top-level `Makefile` to point to the directories containing the HTSlib headers and libraries.
- To enable libdeflate support in HTSlib, install `libdeflate-dev` before running `./configure` for HTSlib.

## Databases used
Reference files can be found at: https://data.broadinstitute.org/lohlab/mCAs_WGS/. 

In brief, the various reference files are:
* `GRCh38.GCmasks.bin` contains binary masks for GC content and blacklisted regions in GRCh38 used in the depth computation
* `lc_sv_cnv_masks.hg38.bed` contains regions excluded by the low complexity, sv, and cnv filters
* `profiles_{20,50}kb.maxDiff0.08.isMale{0,1}.PCs.txt` and `profiles_{20,50}kb.maxDiff0.08.isMale{0,1}.baselines.txt` contain PC coefficients and bin-level baseline depth measurments computed in UK Biobank
* `10k.WGS.ref_bias.txt.gz` contains the empirical REF-bias for MAF>0.01 variants in UK Biobank


## Contact
For questions or contributions, please open an issue.
