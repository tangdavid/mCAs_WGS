##############################################################################
# Makefile cleaned and documented
# - added clearer variable names and flags

# Primary binary
BFILE := bin/mCAs_WGS
# - added .PHONY and a 'help' target
# - reduced repetition for read-depth binaries
##############################################################################

BOOST_INC ?= /usr/include
BOOST_LIB ?= /usr/lib/x86_64-linux-gnu
LIBDEFLATE_LIB ?= /usr/lib/x86_64-linux-gnu
# Path to HTSlib locations. You can either set HTSLD_INC and HTSLD_LIB separately,
# or set a single HTSLD prefix and the Makefile will derive include/lib from it.
HTSLD ?= htslib-1.18-libdeflate-nocurl
# Path to HTSlib include and lib directories (can be left empty if HTSLD is set)
HTSLD_INC ?=
HTSLD_LIB ?=
# Optional Eigen include path (not required for main build). Set if you have Eigen installed locally.
EIGEN_INC ?=

CXX := g++ -std=c++17
# Compiler flags
CXXFLAGS_REL := -O3 -fopenmp
CXXFLAGS_DBG := -O0 -g
CXXFLAGS_WRN := -Wall -Wextra -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated -Wno-unused-parameter

# Default target-specific flags: release by default

all: CXXFLAGS := $(CXXFLAGS_REL) $(CXXFLAGS_WRN)
all: LDFLAGS := $(CXXFLAGS_REL)
all: $(BFILE) depth

debug: CXXFLAGS := $(CXXFLAGS_DBG) $(CXXFLAGS_WRN)
debug: LDFLAGS := $(CXXFLAGS_DBG)
debug: $(BFILE) depth

# Files and paths
HFILE := $(shell find src -name '*.h')
CFILE := $(shell find src -name '*.cpp' | LC_ALL=C sort)
OFILE := obj/hmm.o obj/make_calls.o obj/prepare_bcf.o obj/run_prepare_bcf.o obj/main.o obj/genome_rules.o obj/compute_cis_fisher_p.o obj/compute_polygenic_drive.o obj/gwas_tools.o
VPATH := $(shell for file in `find src -name '*.cpp' | LC_ALL=C sort`; do echo $$(dirname $$file); done)

# If a single HTSLD prefix is provided, derive include/lib from it when missing.
ifneq ($(strip $(HTSLD)),)
ifeq ($(strip $(HTSLD_INC)),)
HTSLD_INC := $(HTSLD)
endif
ifeq ($(strip $(HTSLD_LIB)),)
HTSLD_LIB := $(HTSLD)
endif
endif

# HTSlib include/lib must be provided. Fail early if not set.
ifeq ($(strip $(HTSLD_INC)),)
$(error HTSLD_INC is not set. Please set HTSLD_INC to the HTSlib include directory (e.g. /usr/local/include or /path/to/htslib/include))
endif

ifeq ($(strip $(HTSLD_LIB)),)
$(error HTSLD_LIB is not set. Please set HTSLD_LIB to the HTSlib library directory (e.g. /usr/local/lib or /path/to/htslib/lib))
endif

# Required HTS flags
HTS_CFLAGS := -I$(HTSLD_INC)
HTS_LDFLAGS := -L$(HTSLD_LIB) -L$(LIBDEFLATE_LIB)
IFLAG := $(HTS_CFLAGS) -I$(BOOST_INC)

# Linker libs (static htslib + boost iostreams + program_options, then dynamic system libs)
LDLIBS := -Wl,-Bstatic -lhts -ldeflate -lboost_iostreams -lboost_program_options -Wl,-Bdynamic -lpthread -lz

##############################################################################
# Primary build
##############################################################################

$(BFILE): $(OFILE)
	$(CXX) $(CXXFLAGS) $^ $(HTS_LDFLAGS) -L$(BOOST_LIB) -o $@ $(LDLIBS)

##############################################################################
# Read-depth helper binaries (kept as small standalone tools)
##############################################################################

READ_DEPTH_PIPELINE := src/read_depth

depth: bin/computeRegionDepthsPCadj bin/computeDepthProfiles bin/countReadsSaveDiscordant bin/generateRefMasks

bin/generateRefMasks: $(READ_DEPTH_PIPELINE)/GenerateRefMasks.cpp
	$(CXX) -O3 -Wall $< -o $@

bin/countReadsSaveDiscordant: $(READ_DEPTH_PIPELINE)/CountReadsSaveDiscordant.cpp
	$(CXX) -O3 -Wall $< -o $@ $(HTS_CFLAGS) $(HTS_LDFLAGS) $(LDLIBS) -L$(BOOST_LIB)

bin/computeDepthProfiles: $(READ_DEPTH_PIPELINE)/ComputeDepthProfiles.cpp
	$(CXX) -O3 -Wall $< -o $@ $(HTS_CFLAGS) $(HTS_LDFLAGS) $(LDLIBS) -L$(BOOST_LIB)

bin/computeRegionDepthsPCadj: $(READ_DEPTH_PIPELINE)/ComputeRegionDepthsPCadj.cpp $(READ_DEPTH_PIPELINE)/RefMasks.hpp $(READ_DEPTH_PIPELINE)/RegionDepthUtils.cpp
	$(CXX) -O3 -Wall $< -o $@ $(HTS_CFLAGS) $(HTS_LDFLAGS) $(LDLIBS) -L$(BOOST_LIB)

##############################################################################
# Standalone tools
##############################################################################

# Build a standalone executable for writing CN-LOH GRM
# This links the write_CNLOH_GRM source with the helper object files it needs
check-eigen:
	@if [ -z "$(strip $(EIGEN_INC))" ]; then \
		echo "EIGEN_INC is not set. Please set EIGEN_INC to the Eigen include directory (e.g. /usr/include/eigen3)"; \
		exit 1; \
	fi

bin/write_CNLOH_GRM: check-eigen src/write_CNLOH_GRM.cpp obj/gwas_tools.o obj/prepare_bcf.o obj/genome_rules.o
	$(CXX) $(CXXFLAGS) $< obj/gwas_tools.o obj/prepare_bcf.o obj/genome_rules.o -I$(EIGEN_INC) $(HTS_CFLAGS) $(HTS_LDFLAGS) -L$(BOOST_LIB) -o $@ $(LDLIBS)

##############################################################################
# Object build rules
##############################################################################

obj/%.o: src/%.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(IFLAG)

obj/genome_rules.o: src/genome_rules.c src/genome_rules.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(IFLAG)

##############################################################################
# Utilities
##############################################################################

.PHONY: all debug clean depth help check-eigen

help:
	@printf "Usage:\n"
	@printf "  make          Build release mCAs_WGS binary (default)\n"
	@printf "  make debug    Build with debug flags\n"
	@printf "  make depth    Build read-depth helper binaries in bin/\n"
	@printf "  make clean    Remove build products\n"
	@printf "\nTargets built by this Makefile: %s\n" "$(BFILE)"

clean:
	rm -f obj/*.o $(BFILE) bin/computeRegionDepthsPCadj bin/computeDepthProfiles bin/countReadsSaveDiscordant bin/generateRefMasks

