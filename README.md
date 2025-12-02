# bio312-vangl1-ntd
Analysis for VANGL1 variant–conservation paper (BIO 312): alignment, conservation scoring, enrichment tests, figures

# VANGL1 Neural Tube Defect (NTD) Variant Analysis  
BIO 312 – Molecular Evolution & Bioinformatics  
Dalisha Mehta

## Overview
This project investigates evolutionary conservation, domain structure, and variant distribution in the VANGL1 gene, a planar cell polarity protein implicated in neural tube defects (NTDs). Using multiple sequence alignments, constraint analyses, phylogenetics, and structural modeling, the study evaluates whether NTD-associated variants disproportionately affect conserved or functionally important regions of the protein.

All analyses were performed on AWS EC2 using Python, pandas, matplotlib, Biopython, and custom scripts.

---

## Project Objectives
1. **Characterize VANGL1 domain architecture** using curated transmembrane and C-terminal annotations.  
2. **Quantify evolutionary conservation** across the full protein and within functional regions.  
3. **Compare constraint at NTD-associated vs. benign variants** using gnomAD data and Rate4Site scores.  
4. **Assess enrichment** of NTD variants in high-constraint residues and the PDZ-binding motif.  
5. **Visualize VANGL1 structure** using AlphaFold-predicted 3D coordinates.  
6. **Reconstruct phylogeny** of VANGL1 and VANGL2 orthologs to place sequence patterns in evolutionary context.

## 1.1  Start the lab, make sure your instance is running on EC2 and log in via ssh.


## 1.2 Clone Lab 13

On the command line, clone the Lab 13 repository.

```bash
git clone git@github.com:Bio312/lab13-$MYGIT.git
```

You may be asked to enter your GitHub username and password. Git will now clone a copy of today's lab into a folder called lab13-$MYGIT (where \$MYGIT is your GitHub username). Go there:

```bash
cd lab13-$MYGIT
```