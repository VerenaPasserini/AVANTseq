# AVANTseq: Automated Variant Analysis for Next-gen Targeted Sequencing in Cancer research

**AVANTseq** is a modular, Snakemake-based workflow for high-confidence somatic variant calling from paired end targeted NGS data. It includes:

- A pipeline for generating a custom **Panel of Normals (PoN)**
- A downstream **variant calling pipeline** using GATK Mutect2 with the generated PoN

---

## Overview

This repository contains two main Snakemake workflows:

1. **PoN Pipeline** – Generates a custom Panel of Normals from a set of normal samples.
2. **AVANTseq Pipeline** – Uses Mutect2 to call somatic variants in tumor samples, leveraging the custom PoN.

Both workflows are modular, configurable via YAML, and built for reproducibility and scalability.

---

## 📁 Repository Structure

```plaintext
AVANTseq/
├── rules/
│   ├── trim.smk                # Trim raw fastq reads with atropos
│   ├── align.smk               # Align trimmed reads using bwa mem
│   ├── qc.smk                  # Check sequencing and alignemnt quality
│   ├── pon.smk                 # Generate a custom PoN and merge with an existing one
│   ├── variants.smk            # Call and annotate somatic variants with mutect2 and funcotator
├── config/
│   ├── config.yaml             # Configuration file containing the paths for required files
│   ├── samples_normal.yaml     # Configuration file containing normal samples list
│   ├── samples_tumor.yaml      # Configuration file containing tumor samples list
├── CreatePoN.smk               # Top-level Snakefile to create custom PoN from normal samples
├── AVANTseq.smk                # Top-level Snakefile to call somatic variants form tumor samples using the PoN previously generated
├── LICENSE.txt                 # MIT license file
└── README.md                   # This file
```

## 🔧 Requirements

To run this pipeline, the following tools must be installed and available in your `PATH`:

- [Atropos](https://atropos.readthedocs.io/) – adapter trimming and filtering  
- [BWA](http://bio-bwa.sourceforge.net/) – read alignment  
- [bcftools](http://www.htslib.org/) – VCF/BAM processing and filtering  
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) – quality control of FASTQ files  
- [GATK 4.x](https://gatk.broadinstitute.org/) – variant calling (Mutect2, etc.)  
- [MultiQC](https://multiqc.info/) – summary reports of QC metrics  
- [Snakemake](https://snakemake.readthedocs.io/) – workflow management  
- [samtools](http://www.htslib.org/) – BAM file processing   
- [vt](https://genome.sph.umich.edu/wiki/Vt) – VCF normalization

---

## 1. [Panel of Normals (PoN) Pipeline](docs/PoN.md)
## 2. [AVANTseq Variant Calling Pipeline](docs/AVANTseq.md)

---

## Configuration Files

The configuration YAML file should define all necessary file paths and sample names required by the pipeline.

### Required paths include:

- `work_dir`: Working directory for input and output files
- `ref_bwa`: Reference genome BWA index file (for alignment)
- `ref_fa`: Reference genome FASTA file (for variant calling)
- `bed`: BED file for coverage metrics
- `baits`: Interval list for baited regions
- `targets`: Interval list for target regions
- `pon`: Public Panel of Normals VCF file
- `merged_pon`: Merged custom and public PoN VCF file (generated with the CreatePoN.smk pipeline)
- `germline_resource`: Germline allele frequency resource VCF
- `vcf_exac`: ExAC common variants VCF file
- `data_source`: Funcotator data source directory (for annotation)

Please refer to the `config.yaml` file provided in the `config/` folder for detailed descriptions of each parameter.

## Sample Files

The `samples_normal.yaml` and `samples_tumor.yaml` files contain a simple list of sample names corresponding to paired-end targeted sequencing raw data files.

- **Content:** Only sample names, without file extensions or formats.
- **Data location:** FASTQ files should be stored in the `fastq/` directory inside the `work_dir` specified in the configuration file.
- **File naming convention:** Each sample should have paired-end FASTQ files named as `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`.

Example:
```yaml
samples:
  - "Sample1"
  - "Sample2"
  - "Sample3"
```

## Tip

- Test the workflow with a dry run:

```bash
snakemake -s AVANTseq.smk --dry-run
```
- This workflow assumes input FASTQ files are gzip-compressed (`.fastq.gz`). If your input files are uncompressed (`.fastq`), please update the `trim.smk` rule accordingly by replacing the expected file extensions.  


## License

This project is licensed under the MIT License. See the `LICENSE.txt` file for details.

---

## Contact

For issues, questions, or contributions, please contact:

**Verena Passerini**  
📧 info@verenapasserini.com 
🔗 [github.com/VerenaPasserini/AVANTseq](https://github.com/VerenaPasserini/AVANTseq)
