---
title: 'scmorph: a Python package for analysing single-cell morphological profiles'
tags:
  - Python
  - biology
  - drug discovery
  - morphological profiling
  - single-cell
authors:
  - name: Jesko Wagner
    orcid: 0000-0001-9805-7192
    affiliation: 1
  - name: Hugh Warden
    orcid: 0000-0002-4308-7316
    affiliation: 1
  - name: Ava Khamseh
    orcid: 0000-0001-5203-2205
    corresponding: true
    affiliation: "1, 2"
  - name: Sjoerd&#160;Viktor&#160;Beentjes
    orcid: 0000-0002-7998-4262
    corresponding: true
    affiliation: "1, 3"
affiliations:
  - name: MRC Human Genetics Unit, Institute of Genetics and Cancer, University of Edinburgh, Edinburgh, UK
    index: 1
  - name: School of Informatics, University of Edinburgh, Edinburgh, UK
    index: 2
  - name: School of Mathematics, University of Edinburgh, Edinburgh, UK
    index: 3
date: 17 October 2024
bibliography: paper_minimal.bib
---

# Summary
`scmorph` is a Python package to analyse single-cell data from morphological
profiling experiments which generate large tabular data. `scmorph` combines
domain-specific methods such as single-cell hit calling and batch correction
with the versatile and scalable [`scverse`](https://scverse.org) tools to offer
feature selection, dimensionality reduction and more. Overall, `scmorph` brings
together a host of single-cell morphological profiling methods, making it
applicable for a wide range of experimental designs and workflows.


# Statement of need
Morphological profiling has become an essential tool in biology and drug
discovery, but there is a lack of open-source software for analysing single-cell
morphological data. Existing solutions are commercial, do not scale to large
datasets, or do not offer single-cell specific methods [@OmtaEtAl2016;
@SerranoEtAl2025]. `scmorph` offers a comprehensive set of methods for analysing
single-cell morphological data, which do not require averaging of features
across cells. By integrating with the growing `scverse` of single-cell tools,
`scmorph` also opens up advanced processing capabilities including access to
deep learning tools [@WolfEtAl2018].

Briefly, `scmorph` provides five modules to analyze morphological profiles:

- Reading and writing (IO). `scmorph` allows reading data from `csv`, `sql`,
  `sqlite`, and `h5ad` files, including from the popular CellProfiler software
  [@StirlingEtAl2021]. Once converted, `scmorph` works with AnnData objects
  stored as h5ad, which track processing steps and can easily be written to disk
  [@VirshupEtAl2024].
- Quality control. `scmorph` integrates two levels of unsupervised quality
  control: image-level and single-cell level. Image-level correction is
  performed with a kNN-based outlier detection method, whereas single-cell
  profiles that are outliers are detected via `pyod` [@LiEtAl2022].
- Preprocessing. Provided functions perform feature selection, compute PCA
  coordinates, and optionally aggregate data. For the first time in the field of
  morphological profiling, `scmorph` integrates scone as batch correction
  function, which retains interpretability of features [@ColeEtAl2019].
  Additionally, the integrated feature selection methods can remove features
  associated with known confounders or with high correlation structures, as is
  common in morphological profiling experiments [@KruskalWallis1952;
  @LinHan2021].
- Plotting. `scmorph` uses scanpy for easy plotting of PCA and UMAP coordinates,
  either in 2D or as cumulative densities, which can be useful for identifying
  technical artifacts such as batch effects [@WolfEtAl2018]. It also provides
  methods for plotting features per experimental group, such as plates.
- Downstream analysis. For experiments focused on profiling non-dynamic
  responses, such as a small molecule library, `scmorph` integrates functions to
  perform hit calling from single-cell profiles using the Kolmogorov–Smirnov
  statistic of single-cells to controls in PCA space. For dynamic systems such
  as differentiating cells, `scmorph` incorporates differential trajectory
  inference modelling via  `slingshot` and `condiments` through the `rpy2`
  translation layer [@StreetEtAl2018; @RouxdeBezieuxEtAl2024].

![Overview of `scmorph` functionality. `scmorph` processes profiles generated
with software such as CellProfiler to facilitate downstream analysis by
performing batch correction, image- and single-cell QC, feature selection, hit
calling and trajectory inference. All methods are built with single-cell
analysis in mind and do not require subsampling.](scmorph_overview.png)

In contrast to the commonly used `pycytominer` package [@SerranoEtAl2025] and
`SPACe` [@StossiEtAl2024], `scmorph` offers (i) interpretable batch correction
techniques compatible with single-cell profiles, (ii) enhanced feature selection
with an adapted Chatterjee correlation coefficient or Kruskal-Wallis test
[@LinHan2021; @KruskalWallis1952], (iii) lineage trajectory inference
[@StreetEtAl2018; @BezieuxEtAl2021], and (iv) the option to analyse
multi-nucleated cells. Compared to `pycytominer`, `scmorph` also performs
single-cell based hit calling. And unlike `SPACe`, `scmorph` is agnostic to the
segmentation and feature extraction methods used upstream and therefore
compatible with CellProfiler. `scmorph` also benefits from improvements of
`AnnData` and `scanpy`, such as enabling out-of-core processing crucial to big
data analysis [@VirshupEtAl2024; @WolfEtAl2018].

Already, `scmorph` has been used to quality control morphological profiling
experiments involving differentiating liver cells [@GrahamEtAl2025]. `scmorph`
is also involved in three projects involving small compound and microRNA
perturbations in the domains of drug discovery and fundamental research,
spanning datasets of >20M cells. Going forward, we envision that `scmorph` will
enable analysis of complex and large morphological profiling experiments.

# Acknowledgements
JW and HW are funded by an MRC Unit Award.

# References