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
