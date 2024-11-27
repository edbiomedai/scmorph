import functools
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from anndata import AnnData

__all__ = ["pca", "umap", "cumulative_density", "ridge_plot"]

pca = functools.partial(sc.pl.pca, annotate_var_explained=True)
pca.__doc__ = "| Copied from :ref:`scanpy.pl.pca`, but with annotate_var_explained=True by default.\n" + str(
    sc.pl.pca.__doc__
)
umap = sc.pl.umap
umap.__doc__ = "| Copied from :ref:`scanpy.pl.umap`.\n" + str(umap.__doc__)


def cumulative_density(
    adata: AnnData,
    x: int | str | list[int] | list[str],
    layer: str = "X",
    color: str | None = None,
    n_col: int = 3,
    xlim: tuple[float, float] | None = None,
    xlabel: str | None = None,
    **kwargs: Any,
) -> plt.Figure:
    """
    Plot cumulative densities of variables in AnnData

    Parameters
    ----------
    adata :class:`~anndata.AnnData`
        AnnData object

    x : Union[int, str, List[int], List[str]]
        Name or index of variable(s) to plot

    layer : str
        Where to find values for the variable. Useful if you want to plot "pca" or "umap" values.
        (default: "X")

    color : str
        Variable in "obs" to color by (default: None)

    n_col : int
        Number of columns to facet by (default: 3)

    xlim : Tuple[float, float]
        Limits of x-axis (default: None)

    xlabel : str
        Label for x-axis (default: None)

    kwargs : Any
        Other arguments passed to seaborn.FacetGrid

    Returns
    -------
    plt.Figure
        Plots of cumulative densities of variables in AnnData
    """
    import numpy as np
    import seaborn as sns

    if not isinstance(x, list):
        x = [x]  # type: ignore[assignment]
    n_col = min(n_col, len(x))  # type: ignore[arg-type]
    if layer == "X":
        x_vals = adata[:, x].to_df()
    else:
        adlayer = f"X_{layer}"
        if isinstance(x[0], int):  # type: ignore[index]
            x_vals = adata.obsm[adlayer][:, x]
        else:
            x_vals = adata.obsm[adlayer].loc[:, x]

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
        id_vars=color,
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

    sns.move_legend(fg, "center right", bbox_to_anchor=(1, 0.5))
    fg.set(ylim=(0, 1))
    if xlabel:
        fg.set(xlabel=xlabel)
    if xlim:
        fg.set(xlim=xlim)
    return fg


# Plot individual features
def ridge_plot(df: pd.DataFrame, x: str, y: str, n_col: int = 1, **kwargs: Any) -> plt.Figure:
    """
    Plot features as ridge plot

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Long DataFrame containing features and categories to include in ridge plot

    x : str
        Name of column containing feature values

    y : str
        Name of column containing category values

    n_col : int
        How many columns to plot over (default: 1)

    kwargs : Any
        Other arguments passed to seaborn.FacetGrid

    Returns
    -------
    :class:`matplotlib.pyplot.Figure`
        Plots of cumulative densities of variables in AnnData
    """
    import warnings

    import seaborn as sns

    warnings.filterwarnings(action="ignore", category=UserWarning, module=r".*seaborn")
    # Adapted from https://seaborn.pydata.org/examples/kde_ridgeplot.html
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    # Initialize the FacetGrid object
    g = sns.FacetGrid(df, col=y, hue=y, aspect=10, height=0.5, sharey=False, col_wrap=n_col, **kwargs)

    # Draw the densities in a few steps
    g.map(
        sns.kdeplot,
        "x",
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

    g.map(label, "x")

    # Set the subplots to overlap
    g.figure.subplots_adjust(hspace=-0.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    return g
