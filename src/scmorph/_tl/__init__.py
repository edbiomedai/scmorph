import scanpy as sc

from ._trajectories import (
    slingshot,
    test_common_trajectory,
    test_differential_differentiation,
    test_differential_progression,
)

leiden = sc.tl.leiden
leiden.__doc__ = (
    "| Copied from :ref:`scanpy.tl.leiden`. [Wolf18]_" + sc.tl.leiden.__doc__
)
