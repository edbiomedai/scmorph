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

Tools to reduce number of features based on correlation or confounder association.

.. autosummary::
    :toctree: generated/

    select_features
    corr
    kruskal_test
    kruskal_filter

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
