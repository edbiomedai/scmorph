from functools import partial
from pathlib import Path

import pytest

import scmorph as sm

HERE: Path = Path(__file__).parent
ROOT = HERE / "_images"


@pytest.fixture
def rohban():
    return sm.datasets.rohban2017_minimal()


def test_plot_pca(pca_result, image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)
    sm.pl.pca(pca_result, color="Image_Metadata_Plate", show=False)
    save_and_compare_images("pca")


def test_plot_umap(umap_result, image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)
    sm.pl.umap(umap_result, color="Image_Metadata_Plate", show=False)
    save_and_compare_images("umap")


def test_cumulative_density(pca_result, image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)
    sm.pl.cumulative_density(
        pca_result, x=[0, 1], color="Image_Metadata_Plate", layer="pca", show=False
    )
    save_and_compare_images("cumulative_density_pca")

    sm.pl.cumulative_density(
        pca_result, x="Nuclei_AreaShape_Area", color="Image_Metadata_Plate", show=False
    )
    save_and_compare_images("cumulative_density_feature")


def test_ridge_plot(rohban, image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=5)
    sm.pl.ridge_plot(rohban, "Nuclei_AreaShape_Area", "Image_Metadata_Plate", show=False)
    save_and_compare_images("ridge_plot")
