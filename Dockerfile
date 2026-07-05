# ============================================================================

# Population-Specific LD Reference Panel - Phase 2 & 3 Docker Image

# Phase 2: LD Block Partitioning

# Phase 3: LD Matrix Calculation & HDF5 Panel Construction

# ============================================================================



FROM continuumio/miniconda3:latest



LABEL maintainer="your_email@example.com"

LABEL description="Docker image for LD Block Partitioning and LD Matrix Calculation"

LABEL version="1.0"



# ============================================================================

# 1. SYSTEM DEPENDENCIES

# ============================================================================

RUN apt-get update && apt-get install -y \

    wget \

    unzip \

    git \

    curl \

    build-essential \

    libssl-dev \

    libcurl4-openssl-dev \

    zlib1g-dev \

    libbz2-dev \

    liblzma-dev \

    parallel \

    screen \

    default-jdk \

    && rm -rf /var/lib/apt/lists/*



# ============================================================================

# 2. CONFIGURE CONDA CHANNELS

# ============================================================================

RUN conda config --add channels defaults && \

    conda config --add channels conda-forge && \

    conda config --add channels bioconda



# ============================================================================

# 3. INSTALL BIOINFORMATICS TOOLS (Phase 2 & 3 specific)

# ============================================================================

RUN conda install -y \

    plink=1.9.0.10c \

    bcftools=1.17 \

    tabix \

    samtools=1.17 \

    r-base=4.2.3 \

    && conda clean --all -y



# ============================================================================

# 4. INSTALL PYTHON PACKAGES

# ============================================================================

RUN pip install --no-cache-dir \

    numpy==1.24.3 \

    pandas==2.0.3 \

    scipy==1.11.2 \

    h5py==3.9.0 \

    ldetect==0.3.3 \

    snpflip==0.1.5 \

    matplotlib==3.7.2 \

    seaborn==0.12.2



# ============================================================================

# 5. INSTALL R PACKAGES (for recombination interpolation)

# ============================================================================

RUN R --vanilla -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')" && \

    R --vanilla -e "BiocManager::install('GenomicRanges', update=FALSE)" && \

    R --vanilla -e "install.packages(c('data.table', 'dplyr'), repos='https://cloud.r-project.org/')"



# ============================================================================

# 6. SET WORKING DIRECTORY

# ============================================================================

WORKDIR /app



# ============================================================================

# 7. CLONE THE REPOSITORY (Phase 2 & 3 scripts)

# ============================================================================

RUN git clone https://github.com/pankaj-iitj/population-specific-reference-panel.git . && \

    rm -rf .git



# ============================================================================

# 8. MAKE SCRIPTS EXECUTABLE

# ============================================================================

RUN chmod +x scripts/*.sh



# ============================================================================

# 9. CREATE OUTPUT DIRECTORIES FOR PHASE 2 & 3

# ============================================================================

RUN mkdir -p \

    /data/input \

    /data/output/phase2_ld_blocks \

    /data/output/phase3_ld_matrices \

    /data/output/reference_panel \

    /data/temp \

    /data/logs



# ============================================================================

# 10. CREATE VOLUMES FOR DATA MOUNTING

# ============================================================================

VOLUME ["/data/input", "/data/output", "/data/reference_genomes", "/data/recombination_maps", "/data/temp"]



# ============================================================================

# 11. SET ENVIRONMENT VARIABLES

# ============================================================================

ENV JAVA_TOOL_OPTIONS="-Xmx4g"

ENV PATH="/opt/conda/bin:$PATH"

ENV PYTHONUNBUFFERED=1

ENV THREADS=8

ENV MEMORY=32



# ============================================================================

# 12. CREATE ENTRYPOINT SCRIPT FOR PHASE 2 & 3

# ============================================================================

RUN mkdir -p /scripts && cat > /scripts/entrypoint.sh << 'EOF'

#!/bin/bash

set -e



# Activate conda environment

source /opt/conda/bin/activate base



# Function to display help

show_help() {

    echo "================================================"

    echo "LD Reference Panel - Phase 2 & 3 Pipeline"

    echo "================================================"

    echo ""

    echo "Usage: docker run ldpanel-phase23:latest <COMMAND> [OPTIONS]"

    echo ""

    echo "Commands:"

    echo "  phase2             Run Phase 2: LD Block Partitioning"

    echo "  phase3             Run Phase 3: LD Matrix Calculation"

    echo "  full               Run both Phase 2 and Phase 3"

    echo "  shell              Start interactive shell"

    echo "  test               Run installation tests"

    echo "  help               Show this help message"

    echo ""

    echo "Example:"

    echo "  docker run -v ~/data:/data ldpanel-phase23:latest phase2 \\"

    echo "    --input-vcf /data/input/variants.vcf.gz \\"

    echo "    --genetic-map /data/recombination_maps/ \\"

    echo "    --output /data/output/phase2_ld_blocks/"

    echo ""

}



# Function to test installation

test_installation() {

    echo "Testing installation..."

    echo ""

    

    echo "✓ Checking Python packages..."

    python3 -c "import numpy, pandas, scipy, h5py, ldetect; print('  All Python packages OK')"

    

    echo "✓ Checking R packages..."

    R --vanilla -q -e "library('GenomicRanges'); cat('  GenomicRanges OK\n')" 2>/dev/null || echo "  Warning: GenomicRanges not fully available"

    

    echo "✓ Checking bioinformatics tools..."

    echo "  PLINK:" $(plink --version | head -1)

    echo "  BCFtools:" $(bcftools --version | head -1)

    echo "  Samtools:" $(samtools --version | head -1)

    echo "  Tabix:" $(tabix -h 2>&1 | grep -i tabix | head -1)

    

    echo ""

    echo "✓ All installations verified successfully!"

}



# Function to run Phase 2

run_phase2() {

    echo "================================================"

    echo "Phase 2: LD Block Partitioning"

    echo "================================================"

    echo ""

    

    # Check if input VCF exists

    if [ ! -f "$INPUT_VCF" ]; then

        echo "ERROR: Input VCF file not found: $INPUT_VCF"

        exit 1

    fi

    

    # Create output directory

    mkdir -p "$OUTPUT_PHASE2"

    

    echo "Input VCF: $INPUT_VCF"

    echo "Genetic Maps: $GENETIC_MAP_DIR"

    echo "Output Directory: $OUTPUT_PHASE2"

    echo ""

    

    # Step 1: Convert VCF to PLINK format

    echo "Step 1: Converting VCF to PLINK binary format..."

    for chr in {1..22}; do

        if bcftools view -r $chr "$INPUT_VCF" 2>/dev/null | grep -q "^#CHROM"; then

            echo "  Processing chromosome $chr..."

            bcftools view -r $chr "$INPUT_VCF" -O z -o "$OUTPUT_PHASE2/chr${chr}.vcf.gz"

            plink --vcf "$OUTPUT_PHASE2/chr${chr}.vcf.gz" \

                  --double-id \

                  --export ped \

                  --make-bed \

                  --out "$OUTPUT_PHASE2/chr${chr}" 2>/dev/null

        fi

    done

    

    # Step 2: Genotype QC

    echo ""

    echo "Step 2: Running genotype quality control..."

    bash /app/scripts/geno_qc.sh

    

    # Step 3: Interpolate recombination rates

    echo ""

    echo "Step 3: Interpolating recombination rates..."

    R --vanilla < /app/scripts/interpolate_recomb.R

    

    # Step 4: Run LDetect for LD block partitioning

    echo ""

    echo "Step 4: Running LDetect for LD block partitioning..."

    bash /app/scripts/block_partition.sh

    

    echo ""

    echo "Phase 2 completed! LD blocks saved to: $OUTPUT_PHASE2"

}



# Function to run Phase 3

run_phase3() {

    echo "================================================"

    echo "Phase 3: LD Matrix Calculation"

    echo "================================================"

    echo ""

    

    # Check if LD blocks exist

    if [ ! -d "$LD_BLOCKS_DIR" ]; then

        echo "ERROR: LD blocks directory not found: $LD_BLOCKS_DIR"

        exit 1

    fi

    

    # Create output directory

    mkdir -p "$OUTPUT_PHASE3"

    

    echo "LD Blocks Directory: $LD_BLOCKS_DIR"

    echo "Output Directory: $OUTPUT_PHASE3"

    echo ""

    

    # Calculate LD matrices

    echo "Calculating LD matrices for each block..."

    bash /app/scripts/calc_LD.sh

    

    # Create HDF5 reference panel

    echo ""

    echo "Creating HDF5 reference panel..."

    python3 /app/scripts/CreateHDF5panel.py

    

    echo ""

    echo "Phase 3 completed! HDF5 reference panel created."

}



# Main script logic

case "${1:-help}" in

    phase2)

        INPUT_VCF="${2:-/data/input/variants.vcf.gz}"

        GENETIC_MAP_DIR="${3:-/data/recombination_maps}"

        OUTPUT_PHASE2="${4:-/data/output/phase2_ld_blocks}"

        run_phase2

        ;;

    phase3)

        INPUT_VCF="${2:-/data/input/variants.vcf.gz}"

        LD_BLOCKS_DIR="${3:-/data/output/phase2_ld_blocks}"

        OUTPUT_PHASE3="${4:-/data/output/phase3_ld_matrices}"

        run_phase3

        ;;

    full)

        INPUT_VCF="${2:-/data/input/variants.vcf.gz}"

        GENETIC_MAP_DIR="${3:-/data/recombination_maps}"

        OUTPUT_PHASE2="${4:-/data/output/phase2_ld_blocks}"

        OUTPUT_PHASE3="${5:-/data/output/phase3_ld_matrices}"

        run_phase2

        run_phase3

        ;;

    test)

        test_installation

        ;;

    shell)

        exec /bin/bash

        ;;

    help)

        show_help

        ;;

    *)

        echo "Unknown command: $1"

        show_help

        exit 1

        ;;

esac

EOF



chmod +x /scripts/entrypoint.sh



# ============================================================================

# 13. HEALTH CHECK

# ============================================================================

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \

    CMD python3 -c "import h5py, ldetect, numpy; exit(0)" || exit 1



# ============================================================================

# 14. SET ENTRYPOINT

# ============================================================================

ENTRYPOINT ["/scripts/entrypoint.sh"]

CMD ["help"]
