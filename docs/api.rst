scmorph API
=============

.. module:: scmorph

Import ``scmorph`` as:

.. code-block:: python

    import scmorph as sm


Reading and writing data: ``io``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: scmorph.io


``scmorph`` can read data from a variety of sources, including data exported by CellProfiler.
Once loaded in, all data is treated as an :doc:`AnnData <anndata:index>` object.
This has the advantage of being a fast, standard format that can be used with many
existing single-cell tools, such as :doc:`scanpy <scanpy:index>`.

.. note::
   If you would like to learn more about the ``h5ad`` file format, please see
   :doc:`anndata <anndata:index>`, which is used to read and write these files.

.. note::
    scmorph only processes continuous, non-radial features, i.e. features like number of nearest neighbors (discrete),
    X/Y coordinates (discrete and unfinformative) and BoundingBox (rotation-sensitive) are discarded.
    You may see a warning message about this: consider this an information rather than as an error.

.. autosummary::
    :toctree: generated/

    read
    read_cellprofiler_csv
    read_cellprofiler_batches
    read_sql
    make_AnnData
    split_feature_names

Preprocessing: ``pp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: scmorph.pp

Preprocessing tools that do not produce output, but modify the data to prepare it for downstream analysis.

Basic Preprocessing
-------------------

.. autosummary::
    :toctree: generated/

    drop_na
    scale
    scale_by_batch

Batch Effects
-------------------

Tools to remove batch effects from single-cell morphological data.

.. autosummary::
    :toctree: generated/

    remove_batch_effects

Feature Selection
-------------------

Tools to reduce number of features based on correlations.

.. autosummary::
    :toctree: generated/

    select_features
    corr

Aggregation
-------------------

Tools to compare aggregate profiles.
Additionally, different distance metrics are available.
For a simple aggregation, use ``aggregate``. For a statistically robust distance
metric, use ``aggregate_mahalanobis``.

.. autosummary::
    :toctree: generated/

    aggregate
    aggregate_ttest
    tstat_distance
    aggregate_pc
    aggregate_mahalanobis

Dimensionality-reduction
----------------------------

Tools to perform dimensionality-reduction.

.. autosummary::
    :toctree: generated/

    pca
    neighbors
    umap

Quality Control: ``qc``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: scmorph.qc

Tools to filter cells and images based on quality control metrics and morphological profiles.
For cells, unsupervised filtering is done using :doc:`pyod <pyod:index>` through ``filter_outliers``.
For images, semi-supervised filtering is done using machine-learning methods trained on
image-level data and a subset of labels with ``qc_images``.

While the former can be performed on any dataset, it is likely not as accurate and
may remove underrepresented cell types.

.. autosummary::
    :toctree: generated/

    filter_outliers
    read_image_qc
    qc_images

Visualization: ``pl``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: scmorph.pl

Tools to plot data, often from dimensionality-reduction techniques.
Most of these functions are wrappers around :doc:`scanpy <scanpy:index>` functions.

.. autosummary::
    :toctree: generated/

    pca
    umap
    cumulative_density
    ridge_plot

Datasets: ``datasets``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: scmorph.datasets

Datasets that are included with ``scmorph`` for testing and demonstration purposes.
Currently, this includes various versions of the data in :cite:p:`Rohban2017`.

.. autosummary::
    :toctree: generated/

    rohban2017
    rohban2017_minimal
    rohban2017_imageQC
    rohban2017_minimal_csv
