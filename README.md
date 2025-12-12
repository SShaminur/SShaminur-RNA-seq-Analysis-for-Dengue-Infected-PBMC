# RNA-seq Pipeline for CIBERSORTx Digital Cytometry in PBMC Dengue Virus Analysis


## Overview

A comprehensive, automated RNA-Seq analysis pipeline for processing raw sequencing data through quality control, alignment, quantification, and normalization to produce gene expression matrices suitable for downstream analysis including CIBERSORTx.

This pipeline performs:
- **Quality Control**: Adapter trimming and quality filtering with Trimmomatic
- **Alignment**: Read mapping to reference genome using HISAT2
- **Quantification**: Gene-level read counting with featureCounts
- **Normalization**: Generation of TPM (Transcripts Per Million) and CPM (Counts Per Million) matrices
- **CIBERSORTx Preparation**: Formatted output files ready for cell type deconvolution

## Key Features

- **Parallel Processing**: Efficient utilization of multi-core systems
- **Resume Capability**: Skips already processed samples
- **Comprehensive Output**: Raw counts, normalized matrices, and diagnostic logs
- **Quality Checks**: Automatic validation of input and output
- **Flexible Configuration**: Easy parameter adjustment for different datasets

## Pipeline Workflow

```mermaid
graph LR
    A[Raw FASTQ Files] --> B[Trimmomatic<br/>Quality Control]
    B --> C[HISAT2<br/>Alignment]
    C --> D[SAMtools<br/>BAM Conversion]
    D --> E[featureCounts<br/>Quantification]
    E --> F[Gene Count Matrix]
    F --> G[TPM/CPM Calculation]
    G --> H[CIBERSORTx-ready Files]
```

## Installation

### Prerequisites

- Linux/Unix environment
- Bash shell
- Sufficient disk space and memory for RNA-Seq analysis

### Dependencies

The pipeline requires the following tools installed and available in `PATH`:

1. **Trimmomatic** (v0.39+) - Quality control
2. **HISAT2** (v2.2.1+) - Read alignment
3. **SAMtools** (v1.13+) - BAM file processing
4. **featureCounts** (v2.0.1+) - Gene quantification
5. **R** (v4.0.0+) with **tidyverse** - TPM calculation (optional but recommended)

### Reference Files Required

- HISAT2 genome index
- GTF/GFF3 annotation file
- Trimmomatic adapters file

## Usage

### Basic Usage

```bash
./RNAseq_Analysis_Pipeline.sh /path/to/fastq_directory
```

### With Options

```bash
# Skip trimming step (use existing trimmed files)
./RNAseq_Analysis_Pipeline.sh --skip-trim /path/to/fastq_directory

# Show help
./RNAseq_Analysis_Pipeline.sh -h
```

### Input Requirements

- Paired-end FASTQ files named with `_R1.fastq.gz` and `_R2.fastq.gz` suffixes
- Files should be in the specified input directory

### Configuration

Edit the following parameters at the top of the script:

```bash
# System resources
TOTAL_THREADS=80
THREADS_PER_TRIMMING_JOB=20

# Trimmomatic parameters
ADAPTERS_FILE="/path/to/adapters.fa"
HEADCROP=7
LEADING=20
TRAILING=20

# HISAT2 parameters
HISAT2_INDEX="/path/to/hisat2/index"
GENOME_ANNOTATION="/path/to/annotation.gtf"
```

## Output Structure

```
project_directory/
├── aligned/                 # HISAT2 alignment results
│   ├── *.sorted.bam        # Sorted BAM files
│   ├── *.bai              # BAM index files
│   └── *_alignment_summary.txt
├── paired/                 # Trimmed FASTQ files
├── trimmed_results/        # Trimmomatic logs
├── quantification/         # featureCounts results
│   ├── *_counts.txt       # Individual sample counts
│   └── gene_id_to_name.tsv
├── normalized/            # Normalized matrices (TPM/CPM)
│
├── final_gene_count_matrix.tsv        # Raw gene counts
├── final_gene_count_matrix_nonzero.tsv # Non-zero genes only
├── TPM_matrix.tsv                    # TPM normalized matrix
├── CIBERSORTx_input.tsv              # Formatted TPM for CIBERSORTx
├── CPM_matrix.tsv                    # CPM normalized matrix
└── gene_lengths.tsv                  # Gene lengths from annotation
```

## Output Files Description

| File | Description | Format |
|------|-------------|--------|
| `final_gene_count_matrix.tsv` | Raw read counts per gene | TSV (genes × samples) |
| `TPM_matrix.tsv` | TPM normalized expression matrix | TSV (genes × samples) |
| `CIBERSORTx_input.tsv` | TPM matrix formatted for CIBERSORTx | TSV (genes × samples) |
| `CPM_matrix.tsv` | CPM normalized expression matrix | TSV (genes × samples) |
| `gene_lengths.tsv` | Gene length information | TSV (gene_id, gene_name, length) |

## Quality Control Metrics

The pipeline generates several QC files:
- Trimmomatic logs: `trimmed_results/*.log`
- HISAT2 alignment summaries: `aligned/*_alignment_summary.txt`
- featureCounts logs: `quantification/*_featurecounts.log`

## CIBERSORTx Compatibility

The pipeline prepares TPM matrices specifically formatted for CIBERSORTx requirements:
- Gene symbols as row names
- Samples as column names
- Minimum TPM filtering applied
- Automatic check for log-space vs linear-space data

## Error Handling

The pipeline includes comprehensive error checking:
- Input file validation
- Tool availability verification
- Process failure detection
- Resume capability from failed steps

## Performance Notes

- Designed for high-performance computing environments
- Memory requirements: ~8GB RAM per thread for alignment
- Storage requirements: ~20-50GB per sample (intermediate files)
- Time: ~2-6 hours per sample depending on data size and hardware

## Troubleshooting

### Common Issues

1. **Missing dependencies**: Ensure all required tools are installed and in `PATH`
2. **Insufficient memory**: Adjust `TOTAL_THREADS` and `THREADS_PER_TRIMMING_JOB`
3. **File naming**: Ensure FASTQ files follow the required naming convention
4. **R package errors**: Install required R packages with `install.packages("tidyverse")`

### Log Files

Check the following log files for error details:
- `trimmed_results/*.log` - Trimmomatic issues
- `aligned/*_hisat2.log` - Alignment problems
- `quantification/*_featurecounts.log` - Counting errors

## Citation

If you use this pipeline in your research, please cite the relevant tools:
- Trimmomatic: Bolger et al., 2014
- HISAT2: Kim et al., 2019
- featureCounts: Liao et al., 2014
- CIBERSORTx: Newman et al., 2019

## License

This pipeline is provided for academic use. Users are responsible for ensuring they have appropriate licenses for all dependent software.

## Support

For issues and questions, please check:
- Tool-specific documentation
- GitHub repository issues
- Bioinformatics community forums

