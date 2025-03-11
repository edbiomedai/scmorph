import os
from pathlib import Path
from textwrap import dedent
from typing import TypedDict, cast

import anndata as ad
import numpy as np
import pytest

import scmorph as sm

######################################################
# These functions are reused from scanpy
# (c) 2025 scverse®
# (c) 2017 F. Alexander Wolf, P. Angerer, Theis Lab
# # Licensed under the BSD 3-Clause License


class CompareResult(TypedDict):
    rms: float
    expected: str
    actual: str
    diff: str
    tol: int


@pytest.fixture
def check_same_image(add_nunit_attachment):
    from urllib.parse import quote

    from matplotlib.testing.compare import compare_images

    def check_same_image(
        expected: Path | os.PathLike,
        actual: Path | os.PathLike,
        *,
        tol: int,
        basename: str = "",
    ) -> None:
        __tracebackhide__ = True

        def fmt_descr(descr):
            return f"{descr} ({basename})" if basename else descr

        result = cast(
            CompareResult | None,
            compare_images(str(expected), str(actual), tol=tol, in_decorator=True),
        )
        if result is None:
            return

        add_nunit_attachment(result["expected"], fmt_descr("Expected"))
        add_nunit_attachment(result["actual"], fmt_descr("Result"))
        add_nunit_attachment(result["diff"], fmt_descr("Difference"))

        result_urls = {
            k: f"file://{quote(v)}" if isinstance(v, str) else v for k, v in result.items()
        }
        msg = dedent(
            """\
            Image files did not match.
            RMS Value:  {rms}
            Expected:   {expected}
            Actual:     {actual}
            Difference: {diff}
            Tolerance:  {tol}
            """
        ).format_map(result_urls)
        raise AssertionError(msg)

    return check_same_image


@pytest.fixture
def image_comparer(check_same_image):
    from matplotlib import pyplot as plt

    def save_and_compare(*path_parts: Path | os.PathLike, tol: int):
        __tracebackhide__ = True

        base_pth = Path(*path_parts)

        if not base_pth.is_dir():
            base_pth.mkdir()
        expected_pth = base_pth / "expected.png"
        actual_pth = base_pth / "actual.png"
        plt.savefig(actual_pth, dpi=40)
        plt.close()
        if not expected_pth.is_file():
            msg = f"No expected output found at {expected_pth}."
            raise OSError(msg)
        check_same_image(expected_pth, actual_pth, tol=tol)

    return save_and_compare


###################### End of reused functions ######################


@pytest.fixture
def adata_fixed_values():
    n_treated = 100
    n_control = 100
    n_features = 2
    np.random.seed(2025)
    X = np.vstack(
        [
            np.random.normal(0, 1, (n_treated, n_features)),
            np.random.normal(1, 1, (n_control, n_features)),
        ]
    )
    obs = np.array(["treated"] * n_treated + ["DMSO"] * n_control)
    adata = ad.AnnData(X, obs=obs)
    adata.obs.columns = ["Treatment"]
    adata.obs["Well"] = np.repeat(range(10), 20)
    return adata


@pytest.fixture
def adata_no_na():
    adata = sm.datasets.rohban2017_minimal()
    sm.pp.drop_na(adata)
    return adata


@pytest.fixture
def pca_result(adata_no_na):
    adata = adata_no_na
    sm.pp.scale(adata)
    sm.pp.pca(adata, n_comps=2)
    return adata


@pytest.fixture
def umap_result(pca_result):
    pca_result = pca_result[:100].copy()
    sm.pp.neighbors(pca_result)
    sm.pp.umap(pca_result)
    return pca_result
