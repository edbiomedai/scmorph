scmorph API
=============

.. module:: scmorph

Import scmorph as:
``import scmorph as sm``

Reading and writing data: `io`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To read and write data from CellProfiler and to `h5ad` files, scmorph offers some convenience functions.

.. note::
    scmorph only processes continuous, non-radial features, i.e. features like number of nearest neighbors (discrete),
    X/Y coordinates (discrete and unfinformative) and BoundingBox (rotation-sensitive) are discarded.
    You may see a warning message about this: consider this an information rather than as an error.

.. autosummary::

    read
    read_cellprofiler
    read_cellprofiler_batches

Preprocessing: `pp`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Preprocessing tools that do not produce output, but modify the data to prepare it for downstream analysis.

Basic Preprocessing
-------------------

.. autosummary::

    pp.drop_na
    pp.scale
    pp.scale_by_batch

Batch Effects
-------------------

Tools to remove batch effects from single-cell morphological data.

.. autosummary::

    pp.remove_batch_effects

Feature Selection
-------------------

Tools to reduce number of features based on correlations.

.. autosummary::

    pp.select_features
    pp.corr_features

Aggregation
-------------------

Tools to compare aggregate profiles. Different distance metrics are available.

.. autosummary::

    pp.aggregate_ttest
    pp.tstat_distance
    pp.aggregate_pca
    pp.aggregate_mahalanobis
