from anndata import AnnData


def calculate_qc_metrics(adata: AnnData) -> AnnData:
    """
    Calculate QC metrics with Scanpy and sensible defaults

    Parameters
    ----------
    adata
            Annotated data matrix

    Returns
    -------
    adata
    """
    import scanpy as sc

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[],
        expr_type="value",
        var_type="feature",
        log1p=False,
        inplace=True,
    )
    return adata
