import warnings
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from matplotlib.axes import Axes
from seaborn.axisgrid import FacetGrid

__all__ = ["pca", "umap", "cumulative_density", "ridge_plot"]


def pca(adata: AnnData, annotate_var_explained: bool = True, **kwargs) -> Axes | list[Axes] | None:
    """
    Principal component analysis :cite:p:`Pedregosa2011`.

    This function wraps the higher-level function :func:`~scanpy.pl.pca` (:cite:p:`Wolf18`).

    Parameters
    ----------
    adata
        AnnData object

    annotate_var_explained
        Annotate the variance explained by each component in the plot axis labels.

    kwargs
        Additional arguments passed to :func:`~scanpy.pl.pca`.
    """
    return sc.pl.pca(
        adata,
        annotate_var_explained=annotate_var_explained,
        **kwargs,
    )


def umap(adata: AnnData, **kwargs: Any) -> Axes | list[Axes] | None:
    """
    Uniform Manifold Approximation and Projection :cite:p:`McInnes2018`.

    This function wraps the higher-level function :func:`~scanpy.pl.pca` (:cite:p:`Wolf18`).

    Parameters
    ----------
    adata
        AnnData object

    kwargs
        Additional arguments passed to :func:`~scanpy.pl.umap`.

    Returns
    -------
    Axes | list[Axes] | None
        UMAP plot
    """
    return sc.pl.umap(adata, **kwargs)


def cumulative_density(
    adata: AnnData,
    x: int | str | list[int | str],
    layer: str = "X",
    color: str | None = None,
    n_col: int = 3,
    xlim: tuple[float, float] | None = None,
    xlabel: str | None = None,
    show: bool = True,
    **kwargs: Any,
) -> FacetGrid:
    """
    Plot cumulative densities of variables in AnnData

    Parameters
    ----------
    adata
        AnnData object

    x
        Name or index of variable(s) to plot

    layer
        Where to find values for the variable. Useful if you want to plot "pca" or "umap" values.


    color
        Variable in "obs" to color by

    n_col
        Number of columns to facet by

    xlim
        Limits of x-axis

    xlabel
        Label for x-axis

    show
        Show the plot

    kwargs
        Other arguments passed to :class:`~seaborn.FacetGrid`

    Returns
    -------
    plt.Figure
        Plots of cumulative densities of variables in AnnData
    """
    import numpy as np
    import seaborn as sns

    if not isinstance(x, list):
        x = [x]
    n_col = min(n_col, len(x))
    if layer == "X":
        x_vals = adata[:, x].to_df()  # type: ignore[index]
    else:
        adlayer = f"X_{layer}"
        if isinstance(x[0], int):
            x_vals = adata.obsm[adlayer][:, x]  # type: ignore[index]
        else:
            x_vals = adata.obsm[adlayer].loc[:, x]  # type: ignore[index]

    if layer == "X":
        col_name = "var"
    elif layer == "pca":
        col_name = "PC"
    elif layer == "umap":
        col_name = "UMAP"
    else:
        col_name = layer[2:]

    df = pd.DataFrame(x_vals)
    if color is not None:
        df = pd.concat(
            [adata.obs[color].reset_index(drop=True), df.reset_index(drop=True)],
            axis=1,
        )

    df = pd.melt(
        df,
        id_vars=[color] if color else None,
        var_name=col_name,
        value_name="value",
    )

    if col_name == "PC":
        var = np.round(adata.uns["pca"]["variance_ratio"] * 100, 1).astype(str)
        var = [f"{i}%" for i in var]
        df.loc[:, col_name] = [f"{i+1}, ({var[i]})" for i in df.loc[:, col_name]]
    if col_name == "UMAP":
        df.loc[:, col_name] = df.loc[:, col_name] + 1

    fg = sns.displot(df, x="value", kind="ecdf", col=col_name, hue=color, col_wrap=n_col, **kwargs)

    if color is not None:
        sns.move_legend(fg, "center right", bbox_to_anchor=(1, 0.5))

    fg.set(ylim=(0, 1))
    if xlabel:
        fg.set(xlabel=xlabel)
    if xlim:
        fg.set(xlim=xlim)

    if show:
        plt.show()

    return fg


def ridge_plot(adata: AnnData, x: str, y: str, layer: str = "X", n_col: int = 1, **kwargs: Any) -> FacetGrid:
    """
    Plot features as ridge plot.

    Helps to distinguish the distribution of a feature across different categories, such as plates.

    Parameters
    ----------
    adata
        Annotated data matrix.

    x
        Name of column containing feature values.

    y
        Name of column containing category values.

    layer
        Where to find values for the variable. Useful if you want to plot "pca" or "umap" values.

    n_col
        How many columns to plot over.

    kwargs
        Other arguments passed to seaborn.FacetGrid.

    Returns
    -------
    A ridge plot with categories split out and colored by.
    """
    warnings.filterwarnings(action="ignore", category=UserWarning, module=r".*seaborn")
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    if layer == "X":
        x_vals = adata[:, x].X.flatten()  # type: ignore[index]
    else:
        adlayer = f"X_{layer}"
        x_vals = adata.obsm[adlayer][:, adata.var_names.get_loc(x)]  # type: ignore[index]

    df = pd.DataFrame({x: x_vals, y: adata.obs[y].values})

    # Initialize the FacetGrid object
    g = sns.FacetGrid(df, col=y, hue=y, aspect=10, height=0.5, sharey=False, col_wrap=n_col, **kwargs)

    # Draw the densities in a few steps
    g.map(
        sns.kdeplot,
        x,
        bw_adjust=0.5,
        clip_on=False,
        fill=True,
        alpha=0.9,
        linewidth=1.5,
    )
    g.map(sns.kdeplot, x, clip_on=False, color="w", lw=1.5, bw_adjust=0.5)

    # passing color=None to refline() uses the hue mapping
    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):  # type: ignore
        ax = plt.gca()
        ax.text(
            0,
            0.2,
            label,
            fontweight="bold",
            color=color,
            ha="right",
            va="center",
            transform=ax.transAxes,
        )

    g.map(label, x)

    # Set the subplots to overlap
    g.figure.subplots_adjust(hspace=-0.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    return g
