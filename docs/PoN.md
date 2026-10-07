## 1. Panel of Normals (PoN) Pipeline

### Purpose

Creates a high-quality PoN VCF file from multiple normal samples and merges it with a public PoN. This PoN helps filter out recurrent sequencing artifacts and germline variants during somatic variant calling. The following rule graph visualizes the PoN Snakemake workflow structure, highlighting rule dependencies and execution order:

![CreatePoN Snakemake rule graph: trimming, BWA alignment, duplicate marking, QC, Mutect2 on normals, GenomicsDBImport, CreateSomaticPanelOfNormals and merging with a public PoN](dag_createpon.png "Rule graph for CreatePoN.smk")

### Configuration

Please refer to the sections below for detailed descriptions of:

- The main configuration file `config/config.yaml`, including reference files, target regions, the public panel of normals (`pon`) and the output path of the merged PoN (`merged_pon`).
- The sample list `config/samples_normal.yaml`, which should contain the list of normal samples.

### Run the PoN Pipeline

```bash
# run from the repository root
snakemake -s CreatePoN.smk --cores 8
```

### PoN Outputs

The PoN pipeline generates the following output files, grouped by analysis step. These include quality control metrics, alignment files, and resources required to build a panel of normals for somatic variant calling with Mutect2.

---

#### Quality Control and Alignment

- `alignment/bams/{sample}.bam`  
  BAM file of aligned reads for each normal sample, generated using **BWA-MEM**, sorted and with duplicates marked.

- `alignment/bams/{sample}.bam.bai`  
  BAM index file for rapid access using **samtools index**.

- `alignment/qc/fastqc/{sample}_R1_fastqc.html`, `alignment/qc/fastqc/{sample}_R2_fastqc.html`  
  Quality control reports from **FastQC** on the raw reads (R1 and R2), assessing read quality, GC bias, adapter contamination, etc.  
  [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- `alignment/qc/multiqc_report.html`  
  Summary report from **MultiQC**, aggregating FastQC results across all normal samples.  
  [MultiQC documentation](https://multiqc.info/)

---

#### Coverage and Hybrid Capture QC

- `variants/qc/CoverageSummary.txt`  
  Read counts over the covered regions for all normal samples, calculated with **deepTools multiBamSummary**.

- `variants/qc/HsMetrics/{sample}_metrics.txt`  
  Hybrid selection metrics from **GATK (Picard) CollectHsMetrics**, including bait coverage and on-target rates.  
  [Picard CollectHsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectHsMetrics)

---

#### Somatic Variant Calling (for PoN generation)

- `variants/mutect2/pon/{sample}.vcf.gz`  
  Per-sample variant calls from **Mutect2** in tumor-only mode, used for PoN construction.  
  [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

---

#### Panel of Normals Output

- `variants/mutect2/pon/custom_pon.vcf.gz`  
  Custom panel of normals VCF file generated from all normal sample calls with **CreateSomaticPanelOfNormals**.

- `variants/mutect2/pon/pon_sorted.vcf.gz`  
  Sorted custom panel of normals.

- `merged_pon.vcf.gz` (path set by `merged_pon` in `config/config.yaml`)  
  Final panel of normals (custom + public), used by `AVANTseq.smk` for somatic variant calling in tumor samples.

---