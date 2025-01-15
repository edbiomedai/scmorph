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
