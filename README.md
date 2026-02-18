# Genomic & Transcriptomic Analysis of Trichoderma atroviride

This repository contains scripts and analyses for de novo genome assembly, functional annotation, comparative genomics, RNA-Seq processing, and co-expression network analysis (WGCNA) of six Trichoderma atroviride strains.

---


## Table of Contents
- [Project Overview](#project-overview)
- [Arborescence](#Arborescence)
- [Data Availability](#Data-Availability)
- [Genomic Analysis](#Genomic-Analysis)
- [Transcriptomic Analysis](#Transcriptomic-Analysis)


---

## Project Overview

This repository contains scripts and pipelines for genomic and transcriptomic analyses of six Trichoderma atroviride strains. The main analyses include:
- De novo genome assembly and quality assessment
- Structural and functional annotation
- Comparative genomics and ortholog analysis
- RNA-Seq processing, differential expression, and co-expression network analysis (WGCNA)

---


## Arborescence

Genomic
├── 1-Assembly
│   ├── Busco.sh              # Genome quality assessment with BUSCO
│   ├── Flye.sh               # De novo genome assembly with Flye
│   ├── Mash_Tree.sh          # Genomic distance calculation & UPGMA tree
│   └── Quast.sh              # Assembly quality assessment with QUAST
├── 2-Structural_annotation
│   ├── EDTA.sh               # Transposable element annotation
│   ├── gffread.sh            # Extract proteomes from GFF files
│   └── Helixer.txt           # Gene prediction
├── 3-Functional_annotation
│   ├── Antismash.txt         # Secondary metabolite gene clusters
│   ├── dbCAN3.txt            # CAZyme annotation
│   ├── eggNOG.txt            # Functional annotation (eggNOG)
│   ├── MEROPS.sh             # Peptidase annotation
│   ├── SignalP_EffectorP.sh  # Effector protein prediction
│   └── STRING.txt            # Gene Ontology annotation
├── 4-Genomic_comparison
│   ├── Orthofinder.sh        # Ortholog group inference
│   └── Upset_graph.R         # UpSet plot of shared ortholog groups
Transcription
├── 1-Alignments_Counts
│   ├── Salmon.sh             # Read alignment and quantification with Salmon
│   └── RNAseq_PCA.R          # PCA on RNA-Seq data
├── 2-AskoR
│   └── RNAseq_DE_AskoR.R     # Differential expression analysis using AskoR
└── 3-WGCNA
    └── WGCNA_RNAseq.R        # Weighted gene co-expression network analysis


## Data Availability

All genomic and transcriptomic data used in this analysis are publicly available:
- Genome assemblies: deposited in [NCBI GenBank / ENA] (provide accession numbers)
- Genome sequencing reads (PacBio HiFi): available under [NCBI SRA / ENA]
- RNA-Seq reads: available under [NCBI SRA / ENA]


## Genomic Analysis

The genomic analysis pipeline is organized into four main steps: 
1) Assembly, 2) Structural Annotation, 3) Functional Annotation, 
and 4) Comparative Genomics. Each step uses dedicated scripts 
to automate analyses and generate outputs for downstream interpretation.

-----------------------
1. Genome Assembly

Objective: Assemble PacBio HiFi reads into complete genomes 
and assess their quality.

- Flye.sh
  Runs Flye for de novo genome assembly of long reads.
  Input: raw PacBio HiFi reads.
  Output: draft genome assembly in FASTA format.

- Quast.sh
  Evaluates assembly quality using QUAST.
  Metrics: N50, L50, total length, number of contigs, GC content.
  Input: genome assembly FASTA files.
  Output: summary statistics and detailed reports.

- Busco.sh
  Assesses completeness of assemblies using BUSCO.
  Input: genome assembly FASTA files.
  Output: BUSCO scores (complete, duplicated, fragmented, missing).

- Mash_Tree.sh
  Computes pairwise genomic distances using Mash.
  Input: genome assemblies.
  Output: distance matrix and Newick-format tree.
  Note: The Newick tree can be imported into MEGA or other phylogenetic 
  software to generate a phylogenetic tree figure.

-----------------------
2. Structural Annotation

Objective: Identify gene models, extract proteomes, and annotate transposable elements.

- Helixer.txt
  Performs gene prediction on assembled genomes.
  Input: genome assembly FASTA.
  Output: predicted gene models in GFF format.

- gffread.sh
  Extracts protein sequences (proteomes) from genome assembly and GFF annotation files.
  Input: genome assembly FASTA and GFF files from Helixer prediction.
  Output: FASTA files containing predicted proteins.

- EDTA.sh
  Annotates transposable elements using the EDTA pipeline.
  Input: genome assemblies.
  Output: TE annotations in GFF and summary tables.

-----------------------

3. Functional Annotation
Objective: Assign functional information to predicted genes and proteins. 
Several of these tools were used online and therefore do not have scripts. 
Accordingly, .txt files were created to describe the workflow and the parameters used.

- eggNOG.txt
  Functional annotation using eggNOG-mapper.
  Assigns gene function, and PFAM domains.
  Input: protein sequences.
  Output: annotation tables (TXT/CSV).

- dbCAN3.txt
  Identifies carbohydrate-active enzymes (CAZymes) using dbCAN3.
  Input: protein sequences.
  Output: CAZyme classification tables.

- MEROPS.sh
  Detects and classifies peptidases using MEROPS database.
  Input: protein sequences.
  Output: peptidase annotation tables.

- SignalP_EffectorP.sh
  Predicts secreted proteins and candidate effectors.
  Input: protein sequences.
  Output: list of secreted/effector-like proteins.

- Antismash.txt
  Detects secondary metabolite biosynthetic gene clusters using antiSMASH.
  Input: genome assemblies.
  Output: cluster predictions and summary tables.

- STRING.txt
  Retrieves Gene Ontology (GO) annotations from the STRING database.
  Input: protein sequences.
  Output: GO terms and functional classification tables.

-----------------------

4. Comparative Genomics
Objective: Identify orthologous gene groups across strains 
and visualize shared gene content.

- Orthofinder.sh
  Infers orthologous groups between strains using OrthoFinder.
  Input: predicted proteomes from all strains.
  Output: orthogroup tables, gene trees, and statistics.

- Upset_graph.R
  Generates UpSet plots to visualize shared orthologs across strains.
  Input: orthogroup tables from OrthoFinder.
  Output: publication-ready figures showing overlap between strains.






## Transcriptomic Analysis

The transcriptomic analysis pipeline includes three main steps: 
1) Read alignment and PCA, 2) Differential Expression, 
and 3) Co-expression Network Analysis (WGCNA). Each step 
uses dedicated scripts to process RNA-Seq data and generate 
interpretable outputs.

-----------------------

1. Read Alignment & PCA
Objective: Align RNA-Seq reads, quantify gene expression, 
and explore global patterns in the dataset.

- Salmon.sh
  Aligns RNA-Seq reads to the reference transcriptomes using Salmon.
  Input: paired-end RNA-Seq FASTQ files, reference transcriptome.
  Output: quantification files (TPM and counts) for each sample.

- RNAseq_PCA.R
  Performs Principal Component Analysis on normalized expression matrices.
  Input: TPM matrices generated by Salmon.
  Output: PCA plots to visualize sample clustering and variance.

---------------------------

2. Differential Expression
Objective: Identify genes significantly regulated between conditions.

- RNAseq_DE_AskoR.R
  Performs differential expression analysis using AskoR.
  Input: count matrices from Salmon, sample metadata.
  Output: tables of differentially expressed genes with fold changes 
          and adjusted p-values.

---------------------------

3. Co-expression Network Analysis (WGCNA)

Objective: Identify modules of co-expressed genes.

- WGCNA_RNAseq.R
  Performs weighted gene co-expression network analysis using WGCNA:
    - Filters low-expressed genes (mean count ≥ 0.5)
    - Variance stabilization with DESeq2 VST
    - Constructs signed hybrid networks
    - Detects modules (soft threshold = 8, min module size = 30, 
      merge cut height = 0.2)
  Input: normalized expression matrices (TPM or VST counts).
  Output: module assignments.






