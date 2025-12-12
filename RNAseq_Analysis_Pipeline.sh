#!/bin/bash

###############################################################################
# RNA-Seq Analysis Pipeline with Trimming, Alignment, Quantification, and TPM
# 
# Description:
#   This script performs RNA-Seq analysis using:
#   1. Trimmomatic for quality control
#   2. HISAT2 for alignment to reference genome
#   3. featureCounts for gene quantification
#   4. Generates comprehensive gene count matrices
#   5. Calculates TPM for CIBERSORTx and other applications
#
# Usage:
#   ./RNAseq_Analysis_Pipeline.sh [options] <input_directory>
#
# Options:
#   -h, --help      Show this help message and exit
#   --skip-trim     Skip the trimming step (use existing files in paired/)
#
# Input Requirements:
#   - Paired-end FASTQ files named with _R1.fastq.gz and _R2.fastq.gz suffixes
#   - Files should be in the specified input directory
#
# Output:
#   - trimmed_results/: Trimmomatic logs
#   - paired/: Trimmed paired-end files
#   - aligned/: HISAT2 alignment results (BAM files)
#   - quantification/: featureCounts results
#   - final_gene_count_matrix.tsv: Final gene count matrix
#   - gene_lengths.tsv: Gene lengths extracted from annotation
#   - TPM_matrix.tsv: TPM normalized matrix (for CIBERSORTx)
#   - CIBERSORTx_input.tsv: Formatted TPM matrix for CIBERSORTx
#
# Dependencies:
#   - Trimmomatic, HISAT2, samtools, featureCounts, R with tidyverse
#   - Reference genome and annotation files
#
# Version: 2.0
# Author: Assistant
# Date: Aug-2025
###############################################################################

##################################################
### CONFIGURATION - EDIT THESE PARAMETERS FIRST ##
##################################################

# System resources
TOTAL_THREADS=80
THREADS_PER_TRIMMING_JOB=20

# Trimmomatic parameters
ADAPTERS_FILE="/home/user/miniconda3/share/trimmomatic/adapters/TruSeq3-PE.fa"
HEADCROP=7
LEADING=20
TRAILING=20
SLIDINGWINDOW="20:20"
MINLEN=40

# HISAT2 parameters
HISAT2_INDEX="/home/user/Mamun/ref_genome/hisat2_index/grch44_snp_tran"
GENOME_ANNOTATION="/home/user/Mamun/ref_genome/gencode.v44.annotation.gtf"

# featureCounts parameters
FEATURECOUNTS_THREADS=64

# TPM calculation
MIN_TPM_VALUE=0.01  # Minimum TPM value to keep (genes with very low expression)

##################################################
### MAIN SCRIPT - EDIT BELOW THIS LINE CAREFULLY #
##################################################

# Initialize color variables
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
reset=$(tput sgr0)

# Calculate maximum parallel jobs
MAX_PARALLEL_JOBS=$((TOTAL_THREADS / THREADS_PER_TRIMMING_JOB))
[ $MAX_PARALLEL_JOBS -lt 1 ] && MAX_PARALLEL_JOBS=1

usage() {
  echo -e "${green}
$(basename "$0") - RNA-Seq Analysis Pipeline with TPM

${yellow}Usage:${reset}
  $(basename "$0") [options] <input_directory>

${yellow}Options:${reset}
  -h, --help      Show this help message and exit
  --skip-trim     Skip the trimming step (use existing files in paired/)

${yellow}Configuration:${reset}
  Threads: $TOTAL_THREADS total, $THREADS_PER_TRIMMING_JOB per trimming job
  HISAT2 Index: $HISAT2_INDEX

${yellow}Output Folders:${reset}
  trimmed_results/ - Trimmomatic logs
  paired/ - Trimmed paired-end files
  aligned/ - HISAT2 alignment results (BAM files)
  quantification/ - featureCounts results
  
${yellow}Output Files:${reset}
  final_gene_count_matrix.tsv - Raw gene count matrix
  final_gene_count_matrix_nonzero.tsv - Non-zero count matrix
  gene_lengths.tsv - Gene lengths from annotation
  TPM_matrix.tsv - TPM normalized matrix
  CIBERSORTx_input.tsv - Formatted TPM for CIBERSORTx
${reset}"
  exit 0
}

create_directories() {
  echo -e "${yellow}Creating required directories...${reset}"
  mkdir -p trimmed_results paired aligned quantification normalized
}

extract_gene_lengths() {
  echo -e "${yellow}Extracting gene lengths from annotation...${reset}"
  
  if [ ! -f "$GENOME_ANNOTATION" ]; then
    echo -e "${red}Error: Annotation file not found: $GENOME_ANNOTATION${reset}"
    return 1
  fi
  
  # Extract gene lengths from GTF annotation (Gencode format)
  # Format: gene_id, gene_name, length
  echo -e "${yellow}Using Perl for GTF parsing...${reset}"
  
  awk -F'\t' '$3 == "gene"' "$GENOME_ANNOTATION" | \
  perl -ne '
    chomp; 
    @f = split/\t/; 
    $attr = $f[8]; 
    $attr =~ s/;\s*/;/g; 
    %h = (); 
    foreach $pair (split(/;/, $attr)) { 
        if($pair =~ /([^=]+)\s+"([^"]+)"/) { 
            $h{$1} = $2; 
        } 
    } 
    if($h{"gene_id"} && $h{"gene_name"}) { 
        $len = $f[4] - $f[3] + 1; 
        print "$h{\"gene_id\"}\t$h{\"gene_name\"}\t$len\n"; 
    }
  ' | sort -u > "gene_lengths.tsv"
  
  # Check if extraction was successful
  if [ -s "gene_lengths.tsv" ]; then
    echo -e "${green}Gene lengths extracted: $(wc -l < gene_lengths.tsv) genes${reset}"
    echo -e "${yellow}First few genes:${reset}"
    head -5 "gene_lengths.tsv"
  else
    echo -e "${red}Error: Failed to extract gene lengths${reset}"
    echo -e "${yellow}Creating empty gene lengths file for now...${reset}"
    touch "gene_lengths.tsv"
  fi
}

check_r_availability() {
  echo -e "${yellow}Checking for R and Rscript...${reset}"
  
  # Check if Rscript is available
  if command -v Rscript &> /dev/null; then
    echo -e "${green}Rscript found: $(which Rscript)${reset}"
    R_VERSION=$(Rscript --version 2>&1 | head -1)
    echo -e "${green}R version: $R_VERSION${reset}"
    return 0
  else
    echo -e "${red}Rscript not found in PATH${reset}"
    echo -e "${yellow}Checking for R...${reset}"
    
    if command -v R &> /dev/null; then
      echo -e "${green}R found: $(which R)${reset}"
      R_VERSION=$(R --version 2>&1 | head -1)
      echo -e "${green}R version: $R_VERSION${reset}"
      
      # Check if we can use R CMD BATCH instead
      echo -e "${yellow}Will use R CMD BATCH for R scripts${reset}"
      return 1
    else
      echo -e "${red}R not found either. TPM calculation will be skipped.${reset}"
      return 2
    fi
  fi
}

calculate_tpm_matrix() {
  echo -e "${yellow}Calculating TPM matrix...${reset}"
  
  # Check if required files exist
  if [ ! -f "final_gene_count_matrix.tsv" ]; then
    echo -e "${red}Error: Count matrix not found${reset}"
    return 1
  fi
  
  if [ ! -s "gene_lengths.tsv" ]; then
    echo -e "${red}Error: Gene lengths file is empty${reset}"
    return 1
  fi
  
  # First check R availability
  check_r_availability
  r_status=$?
  
  if [ $r_status -eq 2 ]; then
    echo -e "${red}Skipping TPM calculation: R/Rscript not available${reset}"
    echo -e "${yellow}Please install R and tidyverse package, then run the pipeline again${reset}"
    return 1
  fi
  
  # Create R script for TPM calculation
  cat > "calculate_tpm.R" << 'EOF'
#!/usr/bin/env Rscript

# Check for required packages
required_packages <- c("tidyverse", "data.table")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Missing required R packages:", paste(missing_packages, collapse=", "), "\n")
  cat("Please install with: install.packages(c('", paste(missing_packages, collapse="', '"), "'))\n", sep="")
  quit(status=1)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Function to calculate TPM
calculate_tpm <- function(count_matrix, gene_lengths) {
  # Ensure count matrix and gene lengths are data frames
  counts <- as.data.frame(count_matrix)
  
  # Read gene lengths
  if (is.character(gene_lengths)) {
    lengths_df <- read.table(gene_lengths, header = FALSE, sep = "\t", 
                            stringsAsFactors = FALSE, col.names = c("gene_id", "gene_name", "length"))
  } else {
    lengths_df <- gene_lengths
  }
  
  # Match gene IDs between counts and lengths
  # Try to match by gene_id (ENSEMBL ID)
  common_genes <- intersect(rownames(counts), lengths_df$gene_id)
  
  if (length(common_genes) == 0) {
    # Try to match by gene_name
    common_genes <- intersect(rownames(counts), lengths_df$gene_name)
    if (length(common_genes) > 0) {
      cat("Matching by gene names...\n")
      # Create a mapping
      name_to_id <- setNames(lengths_df$gene_id, lengths_df$gene_name)
      # Reorder counts to match lengths
      counts <- counts[common_genes, , drop = FALSE]
      rownames(counts) <- name_to_id[common_genes]
      common_genes <- rownames(counts)
    }
  }
  
  cat("Genes in count matrix:", nrow(counts), "\n")
  cat("Genes with length information:", length(common_genes), "\n")
  
  if (length(common_genes) == 0) {
    stop("No common genes found between count matrix and gene lengths")
  }
  
  # Subset to common genes
  counts_subset <- counts[common_genes, , drop = FALSE]
  lengths_subset <- lengths_df[lengths_df$gene_id %in% common_genes, ]
  
  # Ensure same order
  lengths_subset <- lengths_subset[match(common_genes, lengths_subset$gene_id), ]
  
  # Step 1: Calculate RPK (Reads Per Kilobase)
  # Divide counts by gene length in kilobases
  lengths_kb <- lengths_subset$length / 1000
  rpk <- counts_subset / lengths_kb
  
  # Step 2: Calculate scaling factors per sample
  scaling_factors <- colSums(rpk) / 1e6
  
  # Step 3: Calculate TPM
  tpm_matrix <- t(t(rpk) / scaling_factors)
  
  # Use gene names as rownames
  rownames(tpm_matrix) <- lengths_subset$gene_name
  
  return(tpm_matrix)
}

# Read data with duplicate handling
cat("Reading count matrix...\n")
count_matrix <- tryCatch({
  # First, check for duplicate gene names
  data <- read.table("final_gene_count_matrix.tsv", 
                     header = TRUE, 
                     row.names = NULL, 
                     sep = "\t",
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
  
  # Check for duplicates in the first column
  gene_names <- data[, 1]
  duplicates <- gene_names[duplicated(gene_names)]
  
  if (length(duplicates) > 0) {
    cat("Found", length(duplicates), "duplicate gene names. Aggregating by sum...\n")
    
    # Aggregate duplicates by sum
    aggregated <- aggregate(data[, -1], by = list(gene = gene_names), FUN = sum)
    
    # Set rownames
    rownames(aggregated) <- aggregated$gene
    aggregated$gene <- NULL
    
    aggregated
  } else {
    # No duplicates, proceed normally
    rownames(data) <- data[, 1]
    data[, -1, drop = FALSE]
  }
}, error = function(e) {
  stop("Error reading count matrix: ", e$message)
})

cat("Count matrix dimensions:", dim(count_matrix), "\n")
cat("First few gene names:", paste(head(rownames(count_matrix), 5), collapse=", "), "\n")

# Check if gene lengths file exists
if (!file.exists("gene_lengths.tsv") || file.size("gene_lengths.tsv") == 0) {
  stop("Gene lengths file is missing or empty")
}

# Calculate TPM
cat("Calculating TPM...\n")
tpm_matrix <- tryCatch({
  calculate_tpm(count_matrix, "gene_lengths.tsv")
}, error = function(e) {
  stop("Error calculating TPM: ", e$message)
})

# Remove genes with very low expression (optional)
cat("Filtering lowly expressed genes...\n")
min_tpm <- 0.01  # Default value
tpm_max <- apply(tpm_matrix, 1, max)
tpm_filtered <- tpm_matrix[tpm_max > min_tpm, ]

cat("Genes before filtering:", nrow(tpm_matrix), "\n")
cat("Genes after filtering (TPM >", min_tpm, "):", nrow(tpm_filtered), "\n")

# Save TPM matrix
write.table(tpm_filtered, 
           "TPM_matrix.tsv", 
           sep = "\t", 
           quote = FALSE, 
           col.names = NA)

# Create summary
cat("\nTPM Matrix Summary:\n")
cat("===================\n")
cat("Samples:", ncol(tpm_filtered), "\n")
cat("Genes:", nrow(tpm_filtered), "\n")
cat("File: TPM_matrix.tsv\n")
cat("\nTPM value ranges:\n")
print(summary(as.vector(tpm_filtered)))
max_val <- max(tpm_filtered, na.rm = TRUE)
cat("\nMaximum TPM value:", max_val, "\n")
if (max_val < 50) {
  cat("WARNING: Maximum TPM value < 50. CIBERSORTx may assume log-space data.\n")
}

# Create CIBERSORTx input
cat("\nCreating CIBERSORTx input...\n")
cibersortx_input <- tpm_filtered

# Remove any duplicate gene symbols
duplicate_genes <- rownames(cibersortx_input)[duplicated(rownames(cibersortx_input))]
if (length(duplicate_genes) > 0) {
  cat("Found", length(duplicate_genes), "duplicate gene symbols. Removing duplicates...\n")
  # Keep the first occurrence of each duplicate
  cibersortx_input <- cibersortx_input[!duplicated(rownames(cibersortx_input)), ]
}

# Write CIBERSORTx input
write.table(cibersortx_input,
           "CIBERSORTx_input.tsv",
           sep = "\t",
           quote = FALSE,
           col.names = NA)

cat("\nCIBERSORTx Input Summary:\n")
cat("=========================\n")
cat("File: CIBERSORTx_input.tsv\n")
cat("Unique genes:", nrow(cibersortx_input), "\n")
cat("Samples:", ncol(cibersortx_input), "\n")
cat("\nFirst few genes:\n")
print(head(rownames(cibersortx_input), 5))

cat("\nDone!\n")
EOF
  
  # Run the R script with appropriate method
  echo -e "${yellow}Running TPM calculation...${reset}"
  
  if command -v Rscript &> /dev/null; then
    Rscript calculate_tpm.R
    r_exit=$?
  elif command -v R &> /dev/null; then
    # Use R CMD BATCH as alternative
    echo -e "${yellow}Using R CMD BATCH instead of Rscript...${reset}"
    R CMD BATCH calculate_tpm.R calculate_tpm.Rout
    r_exit=$?
    # Check the output file for errors
    if [ -f "calculate_tpm.Rout" ]; then
      tail -20 calculate_tpm.Rout
    fi
  else
    echo -e "${red}No R interpreter found. Skipping TPM calculation.${reset}"
    r_exit=1
  fi
  
  if [ $r_exit -eq 0 ]; then
    echo -e "${green}TPM calculation completed successfully!${reset}"
  else
    echo -e "${red}Error in TPM calculation${reset}"
    echo -e "${yellow}Check calculate_tpm.R or calculate_tpm.Rout for details${reset}"
  fi
  
  # Clean up R script
  rm -f calculate_tpm.R
  
  # Show file sizes if TPM files were created
  if [ -f "TPM_matrix.tsv" ]; then
    echo -e "${yellow}Generated TPM files:${reset}"
    ls -lh TPM_matrix.tsv CIBERSORTx_input.tsv 2>/dev/null || echo "TPM files not found"
  fi
  
  return $r_exit
}

create_cpm_matrix() {
  echo -e "${yellow}Creating CPM matrix (optional)...${reset}"
  
  if [ ! -f "final_gene_count_matrix.tsv" ]; then
    echo -e "${red}Error: Count matrix not found${reset}"
    return 1
  fi
  
  cat > "calculate_cpm.R" << 'EOF'
#!/usr/bin/env Rscript

# Simple CPM calculation with duplicate handling
calculate_cpm <- function(counts) {
  # Calculate CPM
  cpm_matrix <- t(t(counts) / colSums(counts)) * 1e6
  return(cpm_matrix)
}

# Read count matrix with duplicate handling
counts <- read.table("final_gene_count_matrix.tsv", 
                    header = TRUE, 
                    row.names = NULL, 
                    sep = "\t",
                    check.names = FALSE,
                    stringsAsFactors = FALSE)

# Check for duplicates in the first column
gene_names <- counts[, 1]
duplicates <- gene_names[duplicated(gene_names)]

if (length(duplicates) > 0) {
  cat("Found", length(duplicates), "duplicate gene names. Aggregating by sum...\n")
  
  # Aggregate duplicates by sum
  aggregated <- aggregate(counts[, -1], by = list(gene = gene_names), FUN = sum)
  
  # Set rownames
  rownames(aggregated) <- aggregated$gene
  aggregated$gene <- NULL
  
  counts_matrix <- aggregated
} else {
  # No duplicates
  rownames(counts) <- counts[, 1]
  counts_matrix <- counts[, -1, drop = FALSE]
}

# Calculate CPM
cpm_matrix <- calculate_cpm(counts_matrix)

# Save CPM matrix
write.table(cpm_matrix, 
           "CPM_matrix.tsv", 
           sep = "\t", 
           quote = FALSE, 
           col.names = NA)

cat("CPM matrix created: CPM_matrix.tsv\n")
cat("Dimensions:", dim(cpm_matrix), "\n")
EOF
  
  # Try to run with available R method
  if command -v Rscript &> /dev/null; then
    Rscript calculate_cpm.R
  elif command -v R &> /dev/null; then
    R CMD BATCH calculate_cpm.R calculate_cpm.Rout
    if [ -f "calculate_cpm.Rout" ]; then
      tail -5 calculate_cpm.Rout
    fi
  else
    echo -e "${red}R not available. Skipping CPM calculation.${reset}"
    rm -f calculate_cpm.R
    return 1
  fi
  
  rm -f calculate_cpm.R
  
  if [ -f "CPM_matrix.tsv" ]; then
    echo -e "${green}CPM matrix created${reset}"
  else
    echo -e "${red}Failed to create CPM matrix${reset}"
  fi
}

run_trimmomatic() {
  local R1=$1
  local R2=$2
  local sample_name=$3
  
  R1_pair="paired/${sample_name}_R1_paired.fastq"
  R1_unpair="paired/${sample_name}_R1_unpaired.fastq"
  R2_pair="paired/${sample_name}_R2_paired.fastq"
  R2_unpair="paired/${sample_name}_R2_unpaired.fastq"
  
  echo -e "${yellow}[Trimming] ${sample_name}${reset}"
  
  trimmomatic PE \
    -threads $THREADS_PER_TRIMMING_JOB \
    -phred33 "$R1" "$R2" "$R1_pair" "$R1_unpair" "$R2_pair" "$R2_unpair" \
    HEADCROP:$HEADCROP \
    ILLUMINACLIP:"$ADAPTERS_FILE":2:30:10:3:TRUE \
    LEADING:$LEADING \
    TRAILING:$TRAILING \
    SLIDINGWINDOW:$SLIDINGWINDOW \
    MINLEN:$MINLEN \
    &>> "trimmed_results/${sample_name}.log"
  
  # Remove unpaired files
  rm -f "$R1_unpair" "$R2_unpair"
  
  if [ $? -eq 0 ]; then
    echo -e "${green}Trimming completed for ${sample_name}${reset}"
  else
    echo -e "${red}Error trimming ${sample_name}${reset}"
    return 1
  fi
}

run_hisat2() {
  local sample_name=$1
  local skip_trim=$2
  
  echo -e "${yellow}[HISAT2 Alignment] ${sample_name}${reset}"
  
  # Determine input files
  if [ "$skip_trim" = true ]; then
    # Use original input files (compressed)
    R1="${sample_name}_R1.fastq.gz"
    R2="${sample_name}_R2.fastq.gz"
  else
    # Use trimmed files (uncompressed)
    R1="paired/${sample_name}_R1_paired.fastq"
    R2="paired/${sample_name}_R2_paired.fastq"
  fi
  
  # Check if input files exist
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo -e "${red}Error: Input files not found for ${sample_name}${reset}"
    echo -e "${red}Looking for: $R1 and $R2${reset}"
    return 1
  fi
  
  # Run HISAT2
  hisat2 -x "$HISAT2_INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    -S "aligned/${sample_name}_aligned.sam" \
    --dta \
    --summary-file "aligned/${sample_name}_alignment_summary.txt" \
    -p $((TOTAL_THREADS / MAX_PARALLEL_JOBS)) \
    &>> "aligned/${sample_name}_hisat2.log"
  
  if [ $? -ne 0 ]; then
    echo -e "${red}HISAT2 alignment failed for ${sample_name}${reset}"
    echo -e "${yellow}Check log file: aligned/${sample_name}_hisat2.log${reset}"
    return 1
  fi
  
  echo -e "${green}HISAT2 alignment completed for ${sample_name}${reset}"
}

convert_sam_to_bam() {
  local sample_name=$1
  
  echo -e "${yellow}[SAM to BAM Conversion] ${sample_name}${reset}"
  
  # Convert SAM to sorted BAM and index
  samtools view -@ $((TOTAL_THREADS / MAX_PARALLEL_JOBS)) -bS "aligned/${sample_name}_aligned.sam" | \
  samtools sort -@ $((TOTAL_THREADS / MAX_PARALLEL_JOBS)) -o "aligned/${sample_name}.sorted.bam" -
  
  if [ $? -ne 0 ]; then
    echo -e "${red}SAM to BAM conversion failed for ${sample_name}${reset}"
    return 1
  fi
  
  # Index BAM file
  samtools index -@ $((TOTAL_THREADS / MAX_PARALLEL_JOBS)) "aligned/${sample_name}.sorted.bam"
  
  # Remove SAM file to save space
  rm "aligned/${sample_name}_aligned.sam"
  
  echo -e "${green}BAM conversion completed for ${sample_name}${reset}"
}

run_featurecounts() {
  local sample_name=$1
  
  echo -e "${yellow}[featureCounts] ${sample_name}${reset}"
  
  # Run featureCounts on individual sample
  featureCounts -T $((TOTAL_THREADS / MAX_PARALLEL_JOBS)) \
    -a "$GENOME_ANNOTATION" \
    -o "quantification/${sample_name}_counts.txt" \
    -p --countReadPairs \
    "aligned/${sample_name}.sorted.bam" \
    &>> "quantification/${sample_name}_featurecounts.log"
  
  if [ $? -eq 0 ]; then
    echo -e "${green}featureCounts completed for ${sample_name}${reset}"
  else
    echo -e "${red}Error running featureCounts on ${sample_name}${reset}"
    echo -e "${yellow}Check log file: quantification/${sample_name}_featurecounts.log${reset}"
    return 1
  fi
}

generate_gene_count_matrix() {
  echo -e "${yellow}Generating final gene count matrix...${reset}"
  
  # Get list of all count files
  count_files=(quantification/*_counts.txt)
  
  if [ ${#count_files[@]} -eq 0 ]; then
    echo -e "${red}No count files found for matrix generation${reset}"
    return 1
  fi
  
  # Create gene mapping from annotation
  echo -e "${yellow}Creating gene mapping...${reset}"
  
  # Extract gene_id to gene_name mapping from GTF
  awk -F'\t' '$3 == "gene"' "$GENOME_ANNOTATION" | \
  perl -ne '
    chomp; 
    @f=split/\t/; 
    $attr=$f[8]; 
    $attr=~s/;\s*/;/g; 
    %h=(); 
    foreach $pair (split(/;/,$attr)) { 
        if($pair=~/([^=]+)\s+"([^"]+)"/) { 
            $h{$1}=$2; 
        } 
    } 
    if($h{"gene_id"} && $h{"gene_name"}) { 
        print "$h{\"gene_id\"}\t$h{\"gene_name\"}\n"; 
    }
  ' | sort -u > "quantification/gene_id_to_name.tsv"
  
  if [ ! -s "quantification/gene_id_to_name.tsv" ]; then
    echo -e "${red}Warning: Could not create gene mapping. Using gene IDs only.${reset}"
    touch "quantification/gene_id_to_name.tsv"
  fi
  
  # Process each count file and extract gene counts
  for count_file in "${count_files[@]}"; do
    sample_name=$(basename "$count_file" | sed 's/_counts\.txt//')
    
    # Skip the summary file
    if [[ "$sample_name" == "all_samples" ]]; then
      continue
    fi
    
    # Check if count file exists and has content
    if [ ! -f "$count_file" ] || [ ! -s "$count_file" ]; then
      echo -e "${red}Error: Count file $count_file is empty or missing${reset}"
      continue
    fi
    
    # Extract counts (skip first two header lines)
    tail -n +3 "$count_file" 2>/dev/null | cut -f1,7 > "quantification/${sample_name}_temp.tsv" 2>/dev/null
    
    # Check if extraction worked
    if [ ! -s "quantification/${sample_name}_temp.tsv" ]; then
      echo -e "${red}Error: Could not extract counts from $count_file${reset}"
      continue
    fi
    
    # Sort by gene ID
    sort -k1,1 "quantification/${sample_name}_temp.tsv" > "quantification/${sample_name}_sorted.tsv"
  done
  
  # Check if we have any sorted files
  sorted_files=(quantification/*_sorted.tsv)
  if [ ${#sorted_files[@]} -eq 0 ]; then
    echo -e "${red}Error: No valid count data found${reset}"
    return 1
  fi
  
  # Create a list of all unique gene IDs from all samples
  echo -e "${yellow}Creating list of all genes...${reset}"
  cat quantification/*_sorted.tsv | cut -f1 | sort -u > "quantification/all_gene_ids.tsv"
  
  # Create header
  echo -n "gene_name" > "quantification/header.tsv"
  for count_file in "${count_files[@]}"; do
    sample_name=$(basename "$count_file" | sed 's/_counts\.txt//')
    if [[ "$sample_name" != "all_samples" ]]; then
      echo -ne "\t$sample_name" >> "quantification/header.tsv"
    fi
  done
  echo "" >> "quantification/header.tsv"
  
  # Create matrix using awk
  echo -e "${yellow}Creating matrix using awk...${reset}"
  
  # First, create a temporary file with gene IDs and names
  if [ -s "quantification/gene_id_to_name.tsv" ]; then
    join -t $'\t' -a 1 -e 0 -o 1.1,2.2 \
      <(sort -k1,1 quantification/all_gene_ids.tsv) \
      <(sort -k1,1 quantification/gene_id_to_name.tsv) \
      > "quantification/gene_id_name_mapping.tsv"
  else
    # If no gene mapping, use gene IDs as names
    sort -k1,1 quantification/all_gene_ids.tsv | \
    awk '{print $1 "\t" $1}' > "quantification/gene_id_name_mapping.tsv"
  fi
  
  # Initialize matrix with gene names
  cut -f2 "quantification/gene_id_name_mapping.tsv" > "quantification/matrix_body.tsv"
  
  # Add counts for each sample
  for count_file in "${count_files[@]}"; do
    sample_name=$(basename "$count_file" | sed 's/_counts\.txt//')
    if [[ "$sample_name" != "all_samples" ]] && [ -f "quantification/${sample_name}_sorted.tsv" ]; then
      echo -e "${yellow}Adding counts for $sample_name...${reset}"
      
      # Create temporary file with counts for this sample
      join -t $'\t' -a 1 -e 0 -o 2.2 \
        <(sort -k1,1 quantification/all_gene_ids.tsv) \
        <(sort -k1,1 "quantification/${sample_name}_sorted.tsv") \
        > "quantification/${sample_name}_counts_only.tsv"
      
      # Paste counts to matrix
      paste "quantification/matrix_body.tsv" "quantification/${sample_name}_counts_only.tsv" > "quantification/matrix_body_temp.tsv"
      mv "quantification/matrix_body_temp.tsv" "quantification/matrix_body.tsv"
      
      # Clean up
      rm "quantification/${sample_name}_counts_only.tsv"
    fi
  done
  
  # Combine header and body
  cat "quantification/header.tsv" "quantification/matrix_body.tsv" > "final_gene_count_matrix.tsv"
  
  # Remove zero-count genes
  awk 'NR==1 {print; next} {sum=0; for(i=2;i<=NF;i++) sum+=$i} sum>0' "final_gene_count_matrix.tsv" > "final_gene_count_matrix_nonzero.tsv"
  
  # Cleanup temporary files
  rm -f quantification/*_temp.tsv quantification/*_sorted.tsv \
         quantification/gene_id_to_name.tsv quantification/all_gene_ids.tsv \
         quantification/header.tsv quantification/matrix_body.tsv \
         quantification/gene_id_name_mapping.tsv
  
  echo -e "${green}Final matrices created:${reset}"
  echo -e "${yellow}  - final_gene_count_matrix.tsv (all genes)${reset}"
  echo -e "${yellow}  - final_gene_count_matrix_nonzero.tsv (genes with expression in at least one sample)${reset}"
  
  # Show preview of the matrix
  echo -e "${yellow}Matrix preview:${reset}"
  head -n 5 "final_gene_count_matrix.tsv" | cut -f1-5 | column -t 2>/dev/null || \
    head -n 5 "final_gene_count_matrix.tsv"
  
  # Show some statistics
  if [ -f "final_gene_count_matrix.tsv" ]; then
    total_genes=$(wc -l < "final_gene_count_matrix.tsv" 2>/dev/null || echo "0")
    nonzero_genes=$(wc -l < "final_gene_count_matrix_nonzero.tsv" 2>/dev/null || echo "0")
    echo -e "${yellow}Total genes: $((total_genes - 1))${reset}"
    echo -e "${yellow}Genes with expression: $((nonzero_genes - 1))${reset}"
  fi
}

process_sample() {
  local sample_name=$1
  local skip_trim=$2
  
  # Step 1: Trimming (unless skipped)
  if [ "$skip_trim" = false ]; then
    run_trimmomatic "${sample_name}_R1.fastq.gz" "${sample_name}_R2.fastq.gz" "$sample_name" || return 1
  fi
  
  # Step 2: HISAT2 alignment
  run_hisat2 "$sample_name" "$skip_trim" || return 1
  
  # Step 3: Convert SAM to BAM
  convert_sam_to_bam "$sample_name" || return 1
  
  # Step 4: featureCounts (individual)
  run_featurecounts "$sample_name" || return 1
  
  echo -e "${green}Completed analysis for ${sample_name}${reset}"
}

main() {
  local skip_trim=false
  local input_dir=""
  
  # Parse arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -h|--help) usage ;;
      --skip-trim) 
        skip_trim=true
        shift ;;
      *)
        if [ -d "$1" ]; then
          input_dir="$1"
          shift
        else
          echo -e "${red}Error: Invalid argument $1${reset}"
          usage
          exit 1
        fi
        ;;
    esac
  done

  if [ -z "$input_dir" ]; then
    echo -e "${red}Error: No input directory specified${reset}"
    usage
    exit 1
  fi

  if [ ! -d "$input_dir" ]; then
    echo -e "${red}Error: Input directory $input_dir does not exist${reset}"
    exit 1
  fi

  # Change to input directory
  cd "$input_dir" || exit 1
  
  create_directories
  
  # Extract gene lengths first (needed for TPM)
  extract_gene_lengths
  
  # Find and process samples
  echo -e "${yellow}Looking for FASTQ files in $(pwd)...${reset}"
  
  local samples_found=0
  local samples_to_process=()
  
  for r1_file in *_R1.fastq.gz; do
    if [[ -f "$r1_file" ]]; then
      sample_name=$(basename "$r1_file" | sed 's/_R1\.fastq\.gz//')
      r2_file="${r1_file/_R1.fastq.gz/_R2.fastq.gz}"
      
      if [ ! -f "$r2_file" ]; then
        echo -e "${red}Error: Missing R2 file for sample ${sample_name}${reset}"
        continue
      fi
      
      # Check if this sample has already been processed
      if [ ! -f "aligned/${sample_name}.sorted.bam" ] || [ ! -f "quantification/${sample_name}_counts.txt" ]; then
        samples_to_process+=("$sample_name")
        samples_found=$((samples_found + 1))
      else
        echo -e "${yellow}Sample ${sample_name} already processed, skipping...${reset}"
      fi
    fi
  done
  
  if [ ${#samples_to_process[@]} -eq 0 ]; then
    echo -e "${yellow}All samples already processed. Generating matrices...${reset}"
    if [ ! -f "final_gene_count_matrix.tsv" ]; then
      generate_gene_count_matrix
    fi
    calculate_tpm_matrix
    create_cpm_matrix
    exit 0
  fi
  
  echo -e "${green}Found ${#samples_to_process[@]} samples to process${reset}"
  
  # Process samples in parallel
  local pids=()
  for sample_name in "${samples_to_process[@]}"; do
    echo -e "\n${green}Processing sample: ${sample_name}${reset}"
    process_sample "$sample_name" "$skip_trim" &
    pids+=($!)
    
    # Limit number of parallel jobs
    if [[ $(jobs -r -p | wc -l) -ge $MAX_PARALLEL_JOBS ]]; then
      wait -n
    fi
  done
  
  # Wait for all jobs to complete and check for errors
  local failed=0
  for pid in "${pids[@]}"; do
    wait $pid
    if [ $? -ne 0 ]; then
      failed=1
    fi
  done
  
  if [ $failed -eq 1 ]; then
    echo -e "${red}Some samples failed to process. Check logs above.${reset}"
  fi
  
  # Generate final gene count matrix after processing all samples
  generate_gene_count_matrix
  
  # Calculate TPM matrix (if count matrix was created)
  if [ -f "final_gene_count_matrix.tsv" ] && [ -s "final_gene_count_matrix.tsv" ]; then
    calculate_tpm_matrix
    create_cpm_matrix
  else
    echo -e "${red}Error: Count matrix not created, skipping TPM calculation${reset}"
  fi
  
  echo -e "${green}\nPipeline completed!${reset}"
  echo -e "${yellow}Processed ${#samples_to_process[@]} samples${reset}"
  
  # List output files
  echo -e "${yellow}Output files generated:${reset}"
  [ -f "final_gene_count_matrix.tsv" ] && echo -e "${yellow}  - Raw counts: final_gene_count_matrix.tsv${reset}"
  [ -f "TPM_matrix.tsv" ] && echo -e "${yellow}  - TPM matrix: TPM_matrix.tsv${reset}"
  [ -f "CIBERSORTx_input.tsv" ] && echo -e "${yellow}  - CIBERSORTx input: CIBERSORTx_input.tsv${reset}"
  [ -f "CPM_matrix.tsv" ] && echo -e "${yellow}  - CPM matrix: CPM_matrix.tsv${reset}"
  [ -f "gene_lengths.tsv" ] && echo -e "${yellow}  - Gene lengths: gene_lengths.tsv${reset}"
  
  # Final check for CIBERSORTx requirements
  if [ -f "CIBERSORTx_input.tsv" ]; then
    echo -e "\n${yellow}CIBERSORTx Input Check:${reset}"
    # Simple check using awk instead of R
    max_value=$(awk 'NR==1 {next} {for(i=2;i<=NF;i++) if($i>max) max=$i} END {print max+0}' CIBERSORTx_input.tsv 2>/dev/null)
    if [ ! -z "$max_value" ] && [ "$max_value" != "0" ]; then
      echo -e "${yellow}Maximum value in CIBERSORTx input: $max_value${reset}"
      if (( $(echo "$max_value < 50" | bc -l 2>/dev/null) )) || [ "$max_value" -lt 50 ] 2>/dev/null; then
        echo -e "${red}WARNING: Maximum value < 50. CIBERSORTx may assume log-space data!${reset}"
      else
        echo -e "${green}OK: Maximum value > 50 (non-log space)${reset}"
      fi
    fi
  fi
}

# Run main function
main "$@"
