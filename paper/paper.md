---
title: 'scmorph: a python package for analysing single-cell morphological profiles'
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
profiling experiments, which generates large tabular data. To facilitate feature
selection, dimensionality reduction, aggregation, and hit calling, `scmorph` is
built on the versatile and scalable [`scverse`](https://scverse.org) tools.
`scmorph` brings together a host of statistically robust methods, making it
applicable for a wide range of experimental designs and workflows.


# Statement of need
Morphological profiling has become an essential tool in biology and drug
discovery, but there is a lack of open-source software tools for analysing
single-cell morphological data. Existing tools are commercial, do not scale to
large datasets, or do not offer single-cell specific tools [@OmtaEtAl2016;
@SerranoEtAl2024]. `scmorph` complements existing tools by providing a
comprehensive set of methods for analysing single-cell morphological data, none
of which require image-level averaging of features. By integrating with the
growing `scverse` of single-cell tools, `scmorph` also opens up advanced
processing capabilities including access to deep learning tools [@WolfEtAl2018].

Briefly, `scmorph` provides five modules to analyze morphological profiles:

- Reading and writing (IO). To enable easy reading in of data generated with the
  popular CellProfiler software, `scmorph` allows reading data from `csv`,
  `sql`, `sqlite`, and `h5ad` files. Once converted, `scmorph` works with
  AnnData objects stored as h5ad, which track processing steps and can easily be
  written out [@VirshupEtAl2024].
- Quality control. `scmorph` can train classifiers to remove low-quality images
  using a subset of labeled images. To further reduce erronous profiles, scmorph
  can also perform single-cell-level quality control through unsupervised
  outlier detection and removal. Of note, the integrated batch correction does
  not change the scale of features, so interpretability is retained
  [@ColeEtAl2019].
- Preprocessing. Provided functions help reduce batch effects, perform feature
  selection, compute PCA coordinates, and optionally aggregate data.
- Plotting. `scmorph` uses scanpy for easy plotting of PCA and UMAP coordinates,
  either in 2D or as cumulative densities, which can be useful for identifying
  technical artifacts such as batch effects [@WolfEtAl2018].
- Trajectory inference. `scmorph` integrates an interface to the trajectory
  inference tools `slingshot` and `condiments` through the `rpy2` translation
  layer [@StreetEtAl2018; @RouxdeBezieuxEtAl2024].

![Overview of `scmorph` functionality. `scmorph` processes profiles generated
with software such as CellProfiler to reduce batch effects, clean images, cells
and features, and perform downstream analyses such as trajectory inference and
hit calling. All methods are built with single-cell analysis in mind and do not
require subsampling.](scmorph_overview.png)

In contrast to the commonly used `pycytominer` package [@SerranoEtAl2024] and
the new `SPACe` software [@StossiEtAl2024], `scmorph` offers (i) batch
correction techniques compatible with single-cell profiles, (ii) enhanced
feature selection with an adapted Chatterjee correlation coefficient or
Kruskal-Wallis test [@LinHan2021; @KruskalWallis1952], (iii) lineage trajectory
inference [@StreetEtAl2018; @BezieuxEtAl2021], and (iv) the option to analyse
multi-nucleated profiles. Compared to `pycytominer`, `scmorph` also performs
single-cell based hit calling. And unlike `SPACe`, `scmorph` is agnostic to the
segmentation and feature extraction methods used and therefore compatible with
CellProfiler. `scmorph` also benefits from improvements of `AnnData` and
`scanpy`, which `scmorph` is based on, which, going forward, will enable
out-of-core computation crucial to big data analysis [@VirshupEtAl2024;
@WolfEtAl2018].

Already, `scmorph` has been used to quality control morphological profiling
experiments involving differentiating liver cells [@GrahamEtAl2023]. `scmorph`
is also involved in three projects involving small compound and microRNA
perturbations in the domains of drug discovery and fundamental research,
spanning datasets of >20M cells. Going forward, we envision that `scmorph` will
enable analysis of complex and large morphological profiling experiments.

# Acknowledgements
JW and HW are funded by an MRC Unit Award.

# References
