# AVANTseq: Automated Variant Analysis for Next-gen Targeted Sequencing in Cancer research

**AVANTseq** is a modular, Snakemake-based workflow for high-confidence somatic variant calling from targeted NGS data. It includes:

- A pipeline for generating a custom **Panel of Normals (PoN)**
- A downstream **variant calling pipeline** using GATK Mutect2 with the generated PoN

---

## 🚀 Overview

This repository contains two main Snakemake workflows:

1. **PoN Pipeline** – Generates a custom Panel of Normals from a set of normal samples.
2. **AVANTseq Pipeline** – Uses Mutect2 to call somatic variants in tumor samples, leveraging the custom PoN.

Both workflows are modular, configurable via YAML, and built for reproducibility and scalability.

---

## 📁 Repository Structure

```plaintext
avantseq/
├── rules/
│   ├── trim.smk                # Trim raw fastq reads with atropos
│   ├── align.smk               # Align trimmed reads using bwa mem
│   ├── qc.smk.smk              # Check sequencing and alignemnt quality
│   ├── pon.smk                 # Generate a custom PoN and merge with an existing one
│   ├── variants.smk            # Call and annotate somatic variants with mutect2 and funcotator
├── config/
│   ├── config.yaml             # Configuration file containing the paths for required files
│   ├── samples_normal.yaml     # Configuration file containing normal samples list
│   ├── samples_tumor.yaml      # Configuration file containing tumor samples list
├── CreatePoN.smk               # Top-level Snakefile to create custom PoN from normal samples
├── AVANTseq.smk                # Top-level Snakefile to call somatic variants form tumor samples using the PoN previously generated
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

## 🧪 1. Panel of Normals (PoN) Pipeline

### Purpose

Creates a high-quality PoN VCF file from multiple normal BAM files. This PoN helps filter out recurrent sequencing artifacts and germline variants during somatic variant calling.

### Configuration

Please refer to the sections below for detailed descriptions of:

- The main configuration file `config/config.yaml`, including reference files, target regions, and optional panel of normals.
- The sample list `config/samples_normal.yaml`, which should contain the list of normal samples.

### Run the PoN Pipeline

```bash
snakemake -s CreatePoN.smk --cores 8 
```

### 📤 PoN Outputs

The PoN pipeline generates the following output files, grouped by analysis step. These include quality control metrics, alignment files, and resources required to build a panel of normals for somatic variant calling with Mutect2.

---

#### 🔹 Quality Control and Alignment

- `alignment/bams/{sample}.bam`  
  BAM file of aligned reads for each normal sample, generated using **BWA-MEM**, sorted and with duplicates marked.

- `alignment/bams/{sample}.bam.bai`  
  BAM index file for rapid access using **samtools index**.

- `alignment/qc/{sample}_fastqc.html`  
  Quality control report from **FastQC**, assessing read quality, GC bias, adapter contamination, etc.  
  [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- `alignment/qc/multiqc_report.html`  
  Summary report from **MultiQC**, aggregating FastQC results across all normal samples.  
  [MultiQC documentation](https://multiqc.info/)

---

#### 🔹 Coverage and Hybrid Capture QC

- `variants/qc/CoverageSummary.txt`  
  Global coverage summary across target regions for all normal samples.

- `variants/qc/HsMetrics/{sample}_metrics.txt`  
  Hybrid selection metrics from **Picard CollectHsMetrics**, including bait coverage and on-target rates.  
  [Picard CollectHsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectHsMetrics)

---

#### 🔹 Somatic Variant Calling (for PoN generation)

- `variants/mutect2/{sample}.vcf.gz`  
  Per-sample variant calls from **Mutect2** in tumor-only mode, used for PoN construction.  
  [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

---

#### 🔹 Panel of Normals Output

- `variants/mutect2/custom_pon.vcf.gz`  
  Raw panel of normals VCF file generated from all normal sample calls.

- `merged_pon.vcf.gz`  
  Final merged panel of normals VCF file, used as an input for somatic variant calling in tumor samples.

---

## 🧬 2. AVANTseq Variant Calling Pipeline

### Purpose

Performs somatic variant calling on tumor samples (optionally with matched normals) using Mutect2. The PoN is used to remove recurrent technical artifacts.

### Configuration

Edit the configuration file `config/config.yaml` with the following:

- Path to input BAM files  
- Reference genome path  
- Target regions (BED or interval list)  
- Optional public PoN file path for merging

Edit the configuration file `config/samples_tumor.yaml` with the following:

- List of tumor sample names  

### Run the AVANTseq Pipeline

```bash
snakemake -s AVANTseq.smk --cores 8 
```

### 🗂 Outputs

The AVANTseq pipeline generates the following output files for each sample, grouped by analysis step. 

---

#### 🔹 Quality Control and Alignment

- `alignment/bams/{sample}.bam`  
  Aligned sequencing reads in **BAM** format, generated with **BWA-MEM** and post-processed (sorted, duplicate-marked).

- `alignment/bams/{sample}.bam.bai`  
  Index file for the BAM, created with **samtools index** for rapid access.

- `alignment/qc/{sample}_fastqc.html`  
  Per-sample quality control report from **FastQC**, assessing read quality, GC content, adapter content, etc.  
  [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- `alignment/qc/multiqc_report.html`  
  Aggregated report from **MultiQC**, summarizing FastQC and other QC metrics across all samples.  
  [MultiQC documentation](https://multiqc.info/)

---

#### 🔹 Coverage and Hybrid Capture QC

- `variants/qc/CoverageSummary.txt`  
  Summary statistics on sequencing depth over the target regions, calculated across all samples.

- `variants/qc/HsMetrics/{sample}_metrics.txt`  
  Hybrid selection metrics generated with **Picard CollectHsMetrics**, including on/off-target efficiency, mean target coverage, and bait performance.  
  [Picard CollectHsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectHsMetrics)

---

#### 🔹 Variant Calling (Mutect2)

- `variants/mutect2/{sample}.vcf.gz`  
  Raw somatic variant calls generated by **GATK Mutect2**, identifying SNVs and indels.  
  [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

- `variants/mutect2/{sample}.getpileupsummaries.table`  
  Table summarizing read counts at common germline sites, used for contamination estimation.

- `variants/mutect2/{sample}.calculatecontamination.table`  
  Estimated contamination levels per sample, based on pileup summaries.

- `variants/mutect2/{sample}.segments.table`  
  Copy number segments inferred by **GATK ModelSegments**, used for downstream filtering.

---

#### 🔹 Filtering and Annotation

- `variants/filtered/{sample}_filtered.vcf.gz`  
  Filtered variant calls using **FilterMutectCalls**, removing likely false positives.

- `variants/filtered/{sample}_filtered_norm_dec.vcf.gz`  
  Normalized and decomposed VCF, prepared for annotation using tools like **vt** 
  [Unified representation of genetic variants](https://academic.oup.com/bioinformatics/article/31/13/2202/196142)

- `variants/maf/{sample}.maf`  
  Final annotated variant list in **Mutation Annotation Format (MAF)**, produced with **Funcotator**, suitable for reporting and downstream analysis.  
  [Funcotator documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037593891-Funcotator)


## 🛠 Configuration File

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

## 📝 samples_*.yaml

The `samples_normal.yaml` and `samples_tumor.yaml` files contains a simple list of sample names corresponding to paired-end targeted sequencing data.

- **Content:** Only sample names, without file extensions or formats.
- **Data location:** FASTQ files should be stored in the `fastq/` directory inside the `work_dir` specified in the config file.
- **File naming convention:** Each sample should have paired-end FASTQ files named as `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`.

Example:
```yaml
samples:
  - Sample1
  - Sample2
  - Sample3
```

## 💡 Tip

- Test workflow with a dry run:

  ```bash
  snakemake -s AVANTseq.smk --dry-run
  ```

## 📜 License

This project is licensed under the MIT License. See the `LICENSE.txt` file for details.

---

## 📬 Contact

For issues, questions, or contributions, please contact:

**Your Name**  
📧 your.email@institute.org  
🔗 [github.com/yourusername](https://github.com/yourusername)
