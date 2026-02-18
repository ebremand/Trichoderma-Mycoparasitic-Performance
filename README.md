# Trichoseed – Comparative Genomics and Transcriptomics of *Trichoderma atroviride*

This repository contains scripts, data, and analysis pipelines for the genomic and transcriptomic study of six *Trichoderma atroviride* strains, focusing on genome comparison, secondary metabolism, CAZymes, and gene co-expression networks in response to pathogen confrontation.  

---

## Table of Contents
- [Project Overview](#project-overview)
- [Materials and Methods](#materials-and-methods)
  - [DNA Sequencing and Genome Comparison](#dna-sequencing-and-genome-comparison)
  - [RNA Sequencing and WGCNA Analysis](#rna-sequencing-and-wgcna-analysis)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Usage](#usage)
- [References](#references)

---

## Project Overview

The aim of this project is to explore the molecular basis of mycoparasitism in *T. atroviride* strains by combining long-read genomic sequencing with transcriptomic profiling under confrontation conditions with phytopathogenic fungi. We integrated genome assembly, annotation, comparative genomics, and weighted gene co-expression network analysis (WGCNA) to identify modules and candidate genes associated with antagonistic potential.

---

## Materials and Methods

### DNA Sequencing and Genome Comparison

- **Strains:** Six *T. atroviride* strains.  
- **DNA Extraction:** High-molecular-weight genomic DNA was extracted from 7-day-old cultures grown on pectocellulosic membranes placed on PDA. DNA was purified using a modified Möller protocol with RNase treatment immediately after proteinase K digestion and resuspended in 10 mM Tris (pH 8).  
- **Sequencing:** Genomes were sequenced at Novogene (Cambridge, UK) using the PacBio Revio platform with HiFi long-read sequencing.  
- **Assembly:** De novo assembly was performed using Flye v2.9.5. Assembly quality was evaluated with QUAST v5.2.0 and BUSCO v5.7.1 (fungi_odb10 database).  
- **Comparative Analysis:**  
  - Pairwise genomic distances were calculated using Mash v2.3.  
  - UPGMA trees were constructed in R using the `ape` package v5.8.1.  
- **Annotation:**  
  - Gene prediction: Helixer v0.3.4; proteome extraction via gffread v0.12.7.  
  - Functional annotation: eggNOG-mapper v2.1.13 (PFAM domains).  
  - Secondary metabolites: antiSMASH v8.0.4.  
  - Gene Ontology: STRING v12.0.  
  - CAZymes: dbCAN3 v12 (only genes supported by ≥2 databases retained).  
  - Peptidases: BLASTp against MEROPS v12.4.  
  - Effector-like proteins: SignalP v5.0 + EffectorP v3.0; classified into PFAM-associated effectors and SSCPs (<300 aa, >3 cysteines, no PFAM).  
- **Orthologs:** OrthoFinder v2.5.2 with Diamond alignments; shared ortholog groups visualized with UpSetR v1.4.0.  

> **Notes:** MAT and MET data for GH family phylogenetic trees are included in the `supplementary/trees` directory.  

---

### RNA Sequencing and WGCNA Analysis

- **Confrontation Assays:**  
  - Six *T. atroviride* strains vs. three phytopathogens (A. brassicicola, R. solani, G. ultimum) on PDA with pectocellulosic membranes.  
  - Sampling: before physical contact and three days after the first sampling. Growth rates were standardized to ensure comparable distances.  
  - Controls: self-confrontation of *T. atroviride*.  
  - Replicates: three independent replicates per condition.  

- **RNA Extraction & Sequencing:**  
  - RNA extraction: NucleoSpin RNA Set for NucleoZOL + RNA Clean & Concentrator-5 kit.  
  - RNA quality: TapeStation (Agilent).  
  - Sequencing: BGI (China), DNBseq™ paired-end 150 bp (~20M reads per sample).  

- **Read Processing:**  
  - Alignment to MMS1508 transcriptome using Salmon v1.10.3.  
  - Quality checks: FastQC v0.12.1, MultiQC v1.22.2.  
  - TPM normalization for PCA (FactoMineR v2.11).  

- **WGCNA Analysis:**  
  - Filtering: genes with mean counts ≥ 0.5; VST normalization (DESeq2 v1.42.0).  
  - Variance filter: genes with variance < 0.3 removed.  
  - Module detection: WGCNA v1.73, signed hybrid TOM, soft threshold = 8, min module size = 30, merge cut height = 0.20.  
  - Module-trait correlations: related to phenotypic measurements.  
  - Candidate genes: selected based on high module membership (MM) and gene significance (GS).  

- **Differential Expression Analysis:**  
  - Package: AskoR v1.0.0 using edgeR v4.0.3.  
  - Criteria: |log2FC| > 1, FDR < 0.05.  
  - Comparisons: pre- vs post-confrontation, weak vs highly mycoparasitic strains, for each pathogen and combined conditions.  

---

## Repository Structure
Trichoseed/
│
├─ data/
│ ├─ counts.txt # RNA-seq count matrix
│ ├─ parasitisme.txt # Trait table for WGCNA
│ ├─ supplementary/trees/ # GH family phylogenetic trees (MAT & MET)
│ └─ genome_files/ # Genome assemblies & annotations
│
├─ scripts/
│ ├─ WGCNA_analysis.R # Full WGCNA workflow
│ ├─ genome_comparison.R # Comparative genomics & ortholog analysis
│ └─ preprocessing_scripts/ # Data preprocessing scripts
│
├─ results/
│ ├─ Gene2Module.csv # Module assignments
│ ├─ MEs.csv # Module eigengenes
│ ├─ GS-MM_trait.csv # Module membership vs gene significance
│ └─ Genes_trait.csv # Candidate genes
│
└─ README.md


---

## Installation

Required software and R packages:  

- R ≥ 3.6  
- R packages: `WGCNA`, `flashClust`, `data.table`, `DESeq2`, `edgeR`, `cluster`, `FactoMineR`, `UpSetR`, `ape`  
- Python 3.8 (optional for additional preprocessing)  
- Tools: Flye, QUAST, BUSCO, Helixer, gffread, eggNOG-mapper, antiSMASH, dbCAN3, MEROPS, SignalP, EffectorP, OrthoFinder, Salmon  

---

## Usage

1. Preprocess counts and trait files in `data/`.  
2. Run `scripts/WGCNA_analysis.R` for co-expression network construction and module-trait analysis.  
3. Run `scripts/genome_comparison.R` for genome assembly QC, annotation, and ortholog comparison.  
4. Inspect results in `results/`.  

---

## References

- Möller et al., 1992. *Nucleic Acids Research*  
- Kolmogorov et al., 2019. *Nature Biotechnology*  
- Gurevich et al., 2013. *Bioinformatics*  
- Manni et al., 2021. *Molecular Biology and Evolution*  
- Ondov et al., 2016. *Genome Biology*  
- Stiehler et al., 2021. *Bioinformatics*  
- Cantalapiedra et al., 2021. *Molecular Biology and Evolution*  
- Blin et al., 2025. *Nucleic Acids Research*  
- Langfelder & Horvath, 2008. *BMC Bioinformatics*  
- Love et al., 2014. *Genome Biology*  
- Carvalho et al., 2021. *BMC Genomics*  

---

*This repository is maintained by Etienne Bremand, IRHS – Université d’Angers, France.*

