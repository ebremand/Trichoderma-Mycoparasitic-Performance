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

The aim of this project is to explore the molecular basis of mycoparasitism in *T. atroviride* strains by combining long-read genomic sequencing with transcriptomic profiling under confrontation conditions with phytopathogenic fungi. We integrated genome assembly, annotation, comparative genomics, and weighted gene co-expression network analysis (WGCNA) to identify modules and candidate genes associated with antagonistic potential.

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



