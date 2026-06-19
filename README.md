# From TP53 Mutations to Network Hubs

## Overview

This project presents an end-to-end cancer genomics workflow designed to investigate how pathogenic TP53 mutations propagate from genomic alterations to transcriptomic dysregulation and network-level molecular vulnerabilities in breast cancer.

By integrating clinical variant data, RNA-seq analysis, functional enrichment, and protein-protein interaction (PPI) networks, this workflow provides a systems-level view of the biological consequences of TP53 dysfunction.

The project combines computational biology, statistical analysis, data visualization, and network science to identify key genes, pathways, and molecular hubs associated with TP53-mutant tumors.

---

## Objectives

* Identify clinically significant TP53 mutation hotspots using ClinVar data.
* Characterize transcriptomic differences between TP53-mutant and TP53 wild-type breast cancer samples.
* Determine biological pathways significantly affected by TP53 disruption.
* Construct and analyze protein interaction networks.
* Identify central hub genes that may represent potential therapeutic targets.
* Establish a foundation for downstream drug repurposing analyses.

---

## Workflow

### Phase 1: TP53 Variant Analysis

Pathogenic single nucleotide variants (SNVs) were extracted from ClinVar and analyzed to identify recurrent mutation hotspots across the TP53 protein.

#### Key Analyses

* Variant filtering and preprocessing
* Pathogenic SNV extraction
* Mutation frequency analysis
* Protein domain mapping
* Hotspot identification

#### Outputs

* Mutation frequency histogram
* Lollipop plot of recurrent mutation positions
* Domain-wise mutation distribution

---

### Phase 2: Differential Expression Analysis

RNA-seq expression data from TCGA-BRCA cohorts were analyzed using DESeq2 to compare TP53-mutant and TP53 wild-type samples.

#### Key Analyses

* Data normalization
* Principal Component Analysis (PCA)
* Differential expression analysis
* Identification of significantly dysregulated genes

#### Outputs

* PCA visualization
* Volcano plot
* Heatmap of top differentially expressed genes

#### Results

* 4,820 significantly differentially expressed genes identified between TP53-mutant and TP53 wild-type cohorts.

---

### Phase 3: Functional Enrichment Analysis

Differentially expressed genes were subjected to pathway enrichment analysis to uncover biological processes and signaling pathways affected by TP53 mutations.

#### Key Analyses

* KEGG pathway enrichment
* Functional annotation
* Biological pathway interpretation

#### Outputs

* KEGG dot plots
* KEGG bar plots
* Enriched pathway summaries

---

### Phase 4: Protein-Protein Interaction Network Analysis

Significant genes were mapped onto protein interaction networks to identify central molecular regulators and network bottlenecks.

#### Key Analyses

* PPI network construction
* Network topology analysis
* Hub gene identification using CytoHubba
* MCC-based ranking of network nodes

#### Outputs

* Cytoscape network visualizations
* Hub gene rankings
* Cluster analysis

#### Key Findings

Several highly connected hub genes were identified, including:

* H4C1
* H3C1
* H3C13
* CXCL1
* CXCL5
* CXCL10

These genes may represent critical downstream effectors of TP53-associated transcriptomic dysregulation.

---

## Technologies and Tools

### Programming Languages

* Python
* R

### Bioinformatics Tools

* DESeq2
* Cytoscape
* CytoHubba

### Databases

* ClinVar
* TCGA-BRCA

### Python Libraries

* Pandas
* NumPy
* Matplotlib
* Seaborn
* Biopython

### Statistical and Visualization Methods

* Differential Expression Analysis
* Principal Component Analysis (PCA)
* Functional Enrichment Analysis
* Protein-Protein Interaction Networks
* Network Topology Analysis

---

## Repository Structure

```text
├── data/
│   ├── raw/
│   └── processed/
│
├── scripts/
│   ├── variant_analysis/
│   ├── differential_expression/
│   ├── enrichment_analysis/
│   └── network_analysis/
│
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
│
├── notebooks/
│
├── README.md
└── requirements.txt
```

---

## Biological Significance

TP53 is one of the most frequently mutated genes in human cancers. While individual mutations are well documented, understanding how these alterations propagate through molecular networks remains a major challenge.

This project bridges that gap by connecting:

**Genomic Variants → Transcriptomic Changes → Pathway Perturbations → Network Hubs**

providing a comprehensive framework for exploring disease mechanisms and identifying potential therapeutic intervention points.

---

## Future Directions

* Drug repurposing through signature reversal analysis
* Integration of additional omics layers
* Survival and clinical outcome analysis
* Multi-cancer cohort comparisons
* Machine learning-based biomarker discovery

---

## Author

**Ahmad bin Shahid**
