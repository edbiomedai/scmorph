import numpy as np
from anndata import AnnData
from scanpy.pp import subsample


def filter_outliers(
    adata: AnnData,
    outliers: float = 0.05,
    fraction: float | None = None,
    n_obs: int | None = None,
    detect_only: bool = False,
    n_cores: int = 1,
) -> AnnData:
    """

    Filter outlier observations from an AnnData object.

    Note
    ----------
    The `outliers` argument determines how many cells will be classified as outlier cells.
    Since it is an arbitrary threshold this will depend on your dataset and downstream analysis.
    We encourage you to try different values and see which one works best for your dataset.

    Parameters
    ----------
    adata
            Annotated data matrix.

    outliers
            Expected fraction of outlier cells.

    fraction: float
            During training, subsample to this `fraction` of the number of observations.

    n_obs
            During training, subsample to this number of observations.
            We recommend 10,000 or fewer, as this results in faster training with adequate accuracy.

    detect_only
            Whether to only detect outliers but not filter them.

    n_cores
            Number of cores to use for parallelization. -1 for all cores.

    Returns
    -------
    adata
    """
    from pyod.models.ecod import ECOD

    model = ECOD(contamination=outliers, n_jobs=n_cores)

    # sampling
    if (fraction is not None and fraction < 1) or (n_obs is not None and n_obs < adata.n_obs):
        adata_ss = subsample(adata, fraction=fraction, n_obs=n_obs, copy=True)

    else:
        adata_ss = adata

    model = model.fit(adata_ss.X)

    chunk_size = min(len(adata), 10000)
    predictions = np.empty(len(adata), dtype=float)
    counter = 0

    for chunk, _, _ in adata.chunked_X(chunk_size=chunk_size):
        predictions[counter : counter + chunk_size] = model.predict(chunk)
        counter += chunk.shape[0]

    adata.obs["outlier"] = predictions

    if not detect_only:
        adata = adata[adata.obs["outlier"] == 0]

    return adata
