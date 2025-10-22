# FASTX IO Wrapper

A flexible `FastxIOWrapper` API for reading and writing **FASTA/FASTQ** (including gzip/bgzip-compressed files) seamlessly.

---

## 🧬 Overview

This project provides two main components:

### 1. `fastx_io_wrapper.py`
A high-level wrapper for FASTA/FASTQ file I/O supporting:
- Automatic format inference (`.fa`, `.fq`, `.gz`, etc.);
- Context-managed read/write operations;
- Support for single-end and paired-end datasets;
- Quality score simulation utilities (`QualitySimulator`).

Provides:
- `FastxIOWrapper`: unified I/O for FASTA/FASTQ (+gzip/bgzip)
- `PairedFastxIOWrapper`: coordinated read/write for paired-end data
- `QualitySimulator`: generates mock PHRED quality scores and quality strings
- Helper functions for inferring format, trimming, pairing, and converting

### 2. `fasta_splitter.py`
A command-line tool and API for:
- Splitting large FASTA sequences into smaller windows (by size or coverage);
- Generating paired-end reads with configurable gaps or overlaps;
- Simulating FASTQ quality scores for downstream testing;
- Writing FASTA or FASTQ output, including compressed files.

# FASTA_Splitter.py

Implements:
- Sequence window splitting by size or coverage;
- Paired-end read generation with optional overlaps;
- Integration with FastxIOWrapper and QualitySimulator;
- Command-line interface (CLI) for batch operations.
- Handles .fa, .fasta, .fq, .fastq, .gz, .bgz automatically.
- Detects format from file extension or input stream content.

## 🚀 Command-Line Usage

```bash
python3 fasta_splitter.py \
  -i genome.fa \
  -w 150 \
  -c 1.0 \
  -p read1.fq.gz read2.fq.gz \
  -g 0
```

Splits genome.fa into 150-bp reads in paired-end format, achieving roughly 1× coverage, outputting to read1.fq.gz and read2.fq.gz in a compressed format (BGZIP), with fixed qualities.

---


# 🧩 Dependencies

- biopython
- tqdm
