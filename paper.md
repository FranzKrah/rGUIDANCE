---
title: '``rGUIDANCE`` – alignment confidence score computation in R'
authors:
- affiliation: 1, 2
  name: Franz-Sebastian Krah
  orcid: 0000-0001-7866-7508
- affiliation: 2
  name: Christoph Heibl
  orcid: 0000-0002-7655-3299
date: "06 August 2019"
output:
  pdf_document: default
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- multiple sequence alignment
- phylogeny
- genome
- guide tree
- R toolbox
- evolution
affiliations:
- index: 1
  name: Plant Biodiversity Research Group, Department of Ecology & Ecosystem Management,
    Technical University of Munich, 85354 Freising, Germany
- index: 2
  name: Bavarian Forest National Park, 94481 Grafenau, Germany
---


# Summary

**R** has become the toolbox of many researchers in ecology and evolutionary biology, both fields increasingly integrating molecular phylogenetic information. This modular toolbox allows to implement pipelines from sequence retrieval and alignment, to phylogeny estimation and high-end visualization of results. One essential step within this pipeline is currently missing from the **R** package ecosystem: detection of unreliably aligned regions within multiple sequence alignment (MSA). Although alignments are fundamental to phylogenetically informed analyses, they often contain extended unreliable regions. The alignment confidence score GUIDANCE demonstrated high accuracy in detecting such regions of low quality. <br />
Here we introduce the **R** package ````rGUIDANCE````, which fully implements the alignment confidence algorithms GUIDANCE, HoT and GUIDANCE2. We will demonstrate the core functionality of ``rGUIDANCE`` and how it can be easily integrated into a phylogeny inference pipeline. <br />
``rGUIDANCE`` is a free and open-source **R** package, available via GitHub at https://github.com/FranzKrah/rGUIDANCE\\

# Introduction
The free software environment for statistical computing and graphics, **R**, has developed as an indispensable toolbox in ecology and evolution and in data science in general [@Muenchen2014]. A myriad of introductory books to statistics, ecological and evolutionary analysis are based on **R** (e.g., https://www.r-project.org/doc/bib/R-books.html, UseR! series). Within ecology and evolution, a growing body of research relies on phylogenetic information [@Cavender-BaresKozak2009], e.g., phylogenetic diversity or ancestral character reconstruction. Once users have gained basic knowledge, the modular **R** toolbox allows to code pipelines, e.g., for phylogeny inference, meeting individual’s needs. Such a pipeline may contain the following steps and **R** packages (only some are listed): (1) sequence retrieval, either manual, semi-automated (*ape*) or automated (*phylotaR*), (2) sequence handling (*ape*, *seqinr*), (3) multiple sequence alignments (MSA) via interface functions to MSA programs (*ape*, *ips*) (4) phylogeny estimation based on MSA via interface functions to tree inference programs (*ape*, *phyclust*, *ips*) (5) divergence time estimation (*ape*); (6) comparative phylogenetic methods such as ancestral character estimation (*ape*, *phytools*, *geiger*) or phylogenetic community ecology such as phylogenetic diversity (*picante*, *vegan*, *MicEco*) and (7) phylogeny visualization (*ape*, *phytools* or *ggtree*). 

Although it is possible to implement such a pipeline in the **R** framework, an essential step is currently not implemented: MSA confidence estimation, which would be situated between step (3) and (4). Currently a typical workflow as outlined above includes a single MSA, which is assumed to be correct. However, benchmark studies show that alignment accuracy is often low, depending on the MSA program and settings used [@NuinWang2006; @ThompsonLinard2011]. Numerous programs have been developed to filter MSA (e.g., trimAI, GBlocks) from ambiguously aligned regions, however simple filtering of alignments was shown to worsen phylogeny estimation in many cases [@TanMuffato2015]. One promising option is to include alignment uncertainty in the phylogeny inference to down-weight sites of low reliability. Such an option is available via a combination of GUIDANCE [@PennPrivman10] and RAxML [@Stamatakis14]. The GUIDANCE column score can be passed as individual weights to each column of the alignment to RAxML via the –a flag (for a practical example: @KrahBassler18). However, GUIDANCE is currently not available in **R** and thus integration of alignment reliability is hampered. Here we thus introduce the **R** package ``rGUIDANCE`` and its core functionality. ``rGUIDANCE`` fully implements GUIDANCE, HoT and GUIDANCE2, which were shown to detect unreliable areas with high accuracy [@PennPrivman10; @SelaAshkenazy15; @LandanGraur08]. ``rGUIDANCE`` further implements basic MSA comparison tools using **Rcpp** to facilitate further development of the **R** toolbox in ecology and evolution.

# Core package functionality
Here, we present the core functions of the ``rGUIDANCE`` package. ``rGUIDANCE`` makes use of several **R** packages that allow interaction with sequence data, e.g., *ape*, *ips*, *adephylo*. Further, ``rGUIDANCE`` uses **Rcpp** to implement MSA comparison tools, which are the basis for GUIDANCE computations. This includes the Cmatrix computation as well as the sum-of-pairs score (SP) [@PennPrivman10]. At the core of ``rGUIDANCE`` are the functions ``guidance``, ``HoT`` and ``guidance2``, which implement algorithms of the same name. Please note that the accuracy of GUIDANCE and ``rGUIDANCE`` scores in identifying alignment errors is almost identical (Fig. 1). Here we present how ``rGUIDANCE`` fits into the **R** ecosystem to build a basic pipeline for phylogeny inference integrating GUIDANCE column confidence scores (Fig. 2).


```
####################### Pipeline ####################### 
## [1] Use phylotaR to retrieve sequence clusters
# for detailed code see Vignette of R package

## [2] Load sequences for the first sequence cluster
cl0 <- read.FASTA("PATHTO/cl0.fas")

## [3] Use GUIDANCE to calculate column score (CS)
g <- guidance(cl0, ncore = 12)
msa0 <- g@msa
sc <- scores(g, "column_raxml", na.rm = FALSE)

## [4] Use alignment and CS for phylogeny inference
tr.w <- raxml(msa0, m = "GTRGAMMA", f = "a", 
              N = 100, p = 1234, x = 1234,
              exec = "PATHTO/raxmlHPC-PTHREADS-AVX", 
              threads = 12, weights = sc$column_raxml,
              outgroup = "Helvella_aestivalis")

## [5] Divergence time estimation
tr.w.ult <- chronos(tr.w$bipartitions)

## [6] Phylogeny visualization
plot(tr.w.ult)
```
For more detailed and up-to-date examples and tutorials, see the ``rGUIDANCE`` GitHub page and vignettes therein. 

# Figures
![Fig. 1 Accuracy of GUIDANCE scores computed with the original GUIDANCE implementation and rGUIDANCE, demonstrating equally high accuracy. Receiver operating characteristic (ROC) curves for GUIDANCE scores (red), HoT scores (blue) and GUIDANCE2 (green) of aligned residue pairs relative to the BAliBASE benchmark database (@ThompsonLinard2011). For detailed explanation of ROC curves and the BALIBASE benchmark dataset, see (@PennPrivman10). We followed the methods described therein. In short: We applied both implementations to each BAliBASE data set (128 datasets), using the MAFFT alignment program, generating GUIDANCE residue pair scores for each pair of aligned residues in the base MSA. We then used the BAliBASE reference MSAs in order to assess the predictive power of the residue pair scores to identify alignment errors. Each aligned residue pair in the base MSA was marked as correct/incorrect by comparing it with the reference MSA (BALiBASE). A receiver operating characteristic (ROC) analysis (@Fawcett2006) was then applied (R package ROCR @SingSander2005) to evaluate the accuracy of the GUIDANCE confidence measure.](paper_fig1_accuracy.png)
\pagebreak


![Fig. 2  The estimated 28S gene tree of *Helvella* received substantially greater overall bootstrap support and has 39% more internal nodes resolved with a bootstrap of 70 or higher, when alignment uncertainty as represented by the GUIDANCE column score is taken into account (red), compared to the estimated topology when alignment uncertainty is ignored (blue). Photo by F.-S. Krah.](paper_fig2_example.png)
\pagebreak

# Conclusions
With the **R** package ``rGUIDANCE`` we hope to enrich the **R** toolbox for ecological and evolutionary research. ``rGUIDANCE`` provides implementations of well-performing MSA reliability score programs, GUDIANCE, HoT and GUIDANCE2. Alignment column scores can be computed and easily integrated into a phylogeny pipeline as was exemplarily demonstrated. The **R** package further provides further functions such as sum-of-pairs scores from MSA comparisons. These functions facilitate the modular development of further MSA reliability scores. Finally, we hope that more phylogeny-based analyses will integrate alignment uncertainty in the phylogeny inference step and thus decrease bias in ecological and evolutionary studies.

# References
