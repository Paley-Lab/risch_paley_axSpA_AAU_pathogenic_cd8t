# HLA-B27 Pathogenic CD8+ T Cell Analysis Scripts

Analysis code to reproduce figures and analyses from the manuscript:

> Risch et al. **Refinement of the Pathogenic T Cell Receptor Motif Suggests Type 17 CD8+ T Cells Initiate HLA-B*27-associated Autoimmunity**. *Journal TBD.*

This repository covers the analyses shown in the Main and Supplementary Figures. 

---

## Table of contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Data availability](#data-availability)
- [Contact](#contact)
- [Citation](#citation)

---

## Overview and Directory Structure

This repository contains the analysis code used to generate the results and figures in the manuscript above. Analyses include:

- Custom read counting and enrichment analysis for the CDR3b saturation mutagenesis experiment (Fig. 1-- directory also contains figures from Supp. Fig. 1)
- Analysis and Criteria 1/2/3 testing in high-throughput TCR sequencing datasets (Fig. 2-- directory also contains figures from Supp. Fig. 2)
- Analysis of the YeiH+CD8+ T cells from peripheral blood of HLA-B27+ healthy controls and axSpA/AAU patients (Fig. 3-- directory also contains figures from Supp. Fig. 3)
- Analysis of the 4 publicly available joint, eye, and blood scRNA-seq and scTCR-seq datasets (Fig. 4-- directory also contains figures from Supp. Fig. 4)
- Analysis of the CD8+ T cells from the gut and blood of HLA-B27+ IBD patients  (Fig. 5-- directory also contains figures from Supp. Fig. 5)
- Plotting of type 17 genes and signatures in the various datasets (Fig. 6-- directory also contains figures from Supp. Fig. 6)

Additionally, the pathogenic gene signature is provided as a .gmt file in this repository. To reproduce the figures in our manuscript, please use the version of the gene signature *without* TCR genes (TRAV21 and TRBV9 removed). 

## Data availability

- **Raw sequencing data:** deposited in dbGaP under accession **[TBD]**.
- **Processed data objects** are available on Zenodo.
  - Full processed high-throughput TCR-seq matrices (Omniscope OS-T and TIRTL-seq) for Figure 2 are available on Zenodo.
  - For the single-cell datasets, metadata matrices needed to reproduce our filtering and annotations from the raw data are available on Zenodo.

## Contact

For questions about the code, please contact:

- Isabel Risch, irisch@wustl.edu or isabelrisch7@gmail.com — analysis lead
- Michael Paley, paleym@wustl.edu — corresponding author

## Citation

If you use this code, please cite:

```
[Citation TBD]
```
