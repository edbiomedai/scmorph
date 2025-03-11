Quality Control: ``qc``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: scmorph.qc

Tools to filter cells and images based on quality control metrics and morphological profiles.
For cells, unsupervised filtering is done using :doc:`pyod <pyod:index>` through ``filter_outliers``.
For images, unsupervised filtering is done with ``qc_images_by_dissimilarity``.

Note that performing unsupervised QC may remove underrepresented cell types or removes images with
suitable quality and should therefore be tailored to your analysis needs.

.. autosummary::
    :toctree: generated/

    filter_outliers
    qc_images_by_dissimilarity
    count_cells_per_group
