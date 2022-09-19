
Usage
=====

This section outlines how scmorph can be used to process CellPainting datasets.
Please note that since scmorph is under active development some of this functionality may change in future releases.
Many of the functions shown below have additional arguments to fit your needs.
If you are missing any functionality, please open a GitHub issue!

.. note:

    If you have previously used Scanpy, many of scmorph's functions will be familiar to you.
    scmorph builds on Scanpy and adds functionality for morphological datasets.
    That said, the underlying data structure is exactly the same.

    If you have never heard of Scanpy: don't worry! The below guide should cover most of your questions.

To use this package, import it:

.. code-block:: python

   import scmorph as sm

Load data
---------

The following command loads a sample morphological dataset created with CellProfiler from a csv file.

.. note::
    SQlite databases are not yet supported.
    Instead, ``scmorph`` reads data as csv files, with cells in rows and features in columns (as produced by CellProfiler)
    There should be some columns containing metadata, followed by features in the remaining columns.

.. code-block:: python

   adata = sm.datasets.multi_plate_experiment() # read data example data

QC
--------

There are three broad ways to QC your data:


#. Supervised image-based QC: you can provide an additional csv file containing per-image labels whether the image is of sufficient quality
#. Unsupervised single-cell QC: tries to automatically detect and remove outlier cells, based on their morphological profiles
#. Supervised single-cell QC: removes e.g. poorly segmented cells using a subset of labeled cells

Supervised image-based QC
^^^^^^^^^^^^^^^^^^^^^^^^^

In a best-case scenario, we have labeled a subset of images with information about their quality.
For example, in our test data set we labeled images based on whether they contained large speckles in the DAPI channel.

.. note::
    Note that here we labeled images on if they contained large particles in the DAPI channel.
    However, this method can be used with any image-level annotations.
    Depending on how small the artifacts are and how many images are labeled the QC may work better or worse.

.. code-block:: python

   img_qc = sm.datasets.image_qc_data() # read image QC data from csv file
   sm.qc.qc_images(adata, img_qc) # perform QC

Unsupervised single-cell QC
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this task, we will try to remove outlier cells, which are presumably poorly measured cells, based on their morphological profiles.
Since manually labeling single-cells is a tedious task, here we provide an unsupervised approach to outlier removal based on `PyOD's <http://jmlr.org/papers/v20/19-011.html>_` implementation of EODF.

.. note::
    This method has two caveats:

    #. It needs an educated guess about the fraction of cells that are of poor quality.
    #. It also requires those poor-quality cells to be at the edges of their feature distribution. In other words, if for example poorly segmented cells are phenotypically similar to real cells, this method will not work.


.. code-block:: python

   fraction_outliers = 0.05 # 5% of cells are expected to be outliers
   adata = sm.qc.filter_outliers(adata, fraction_outliers)

Supervised single-cell QC
^^^^^^^^^^^^^^^^^^^^^^^^^

For a tutorial on how to perform single-cell supervised QC, please see: :doc:` Building a single-cell QC model <tutorials/outlier_removal>`.

Dimensionality reduction
------------------------

To get an idea of if cells cluster distinctly, we can perform a principle component analysis. This works best if we scale our data first.

.. code-block:: python

   sm.pp.scale(adata) # scale data
   sm.pp.pca(adata) # compute PCA
   sm.pl.pca(adata) # plot PCA

For a non-linear transformation of the data, we can just as easily visualize cells in a UMAP:

.. code-block:: python

   sm.pp.neighors(adata)
   sm.pp.umap(adata)
   sm.pl.umap(adata)

Saving data
^^^^^^^^^^^

Saving processed data is also straightforward:

.. code-block:: python

   out_path = "output_file.h5ad"
   adata.write(out_path)
