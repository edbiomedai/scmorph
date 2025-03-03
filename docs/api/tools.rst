Analysis tools: ``tl``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: scmorph.tl


``scmorph`` can perform trajectory analysis using the `Slingshot <https://bioconductor.org/packages/release/bioc/html/slingshot.html>`_ algorithm.
Additionally, it can compute the single-cell distance between treatments and control cell, which can be used for hit calling.

.. autosummary::
    :toctree: generated/

    get_ks
    slingshot
    test_common_trajectory
    test_differential_differentiation
    test_differential_progression
