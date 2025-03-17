import warnings
from abc import ABC, abstractmethod

import anndata as ad
import numpy as np
import numpy.typing
import pandas as pd
from anndata import ImplicitModificationWarning
from scipy.spatial.distance import mahalanobis
from scipy.stats import kstest
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

import scmorph as sm


class Iterator(ABC):
    """
    A base iterator class to iterate over groups of data.

    Attributes
    ----------
    name
        The name of the iterator.
    index
        The current index of the iterator.
    groups
        A list of groups to iterate over.
    copy
        Whether to return a copy of the data.
    """

    def __init__(self):
        self.name = "Iterator"
        self.index = 0
        self.groups = []
        self.copy = False

    def __iter__(self) -> "Iterator":
        return self

    @abstractmethod
    def __next__(self) -> tuple:
        pass

    def __len__(self) -> int:
        return len(self.groups)


class AnnDataIterator(Iterator):
    """
    A class to iterate over AnnData objects.

    Parameters
    ----------
    adata
        The AnnData object to iterate over.
    id_col
        The column(s) to group by, by default ["PlateID"].
    copy
        Whether to return a copy of the data, by default False.
    as_plate
        Whether to return the data as a Plate object, by default False.
    """

    def __init__(
        self,
        adata: ad.AnnData,
        id_col: str | list[str] = "PlateID",
        copy: bool = False,
        as_plate: bool = False,
    ):
        super().__init__()
        self.name = "AnnDataIterator"
        self.adata = adata
        self.id_col = id_col
        self.copy = copy
        self.groups = list(adata.obs.groupby(id_col, observed=True))
        self.as_plate = as_plate

    def __next__(self) -> tuple:
        if self.index >= len(self.groups):
            raise StopIteration
        group_id, group_obs = self.groups[self.index]
        group_itx = group_obs.index
        self.index += 1
        ret_adata = self.adata[group_itx]
        if self.copy:
            ret_adata = ret_adata.copy()
        if self.as_plate:
            ret_adata = Plate(ret_adata, neg_wells=None, neg_control=None)
        return group_id, ret_adata


class PlateCollectionIterator(Iterator):
    """
    A class to iterate over AnnData collections.

    Parameters
    ----------
    adata
        The list of AnnData objects to iterate over.
    copy
        Whether to return a copy of the data, by default False.
    """

    def __init__(self, adata, copy=False):
        super().__init__()
        self.name = "PlateCollectionIterator"
        self.adata = adata
        self.copy = copy
        self.groups = [i.id for i in adata]
        self.index = 0

    def __next__(self) -> tuple:
        if self.index >= len(self.groups):
            raise StopIteration
        group_id = self.groups[self.index]
        group_itx = self.index
        self.index += 1
        if self.copy:
            return group_id, self.adata[group_itx].copy()
        else:
            return group_id, self.adata[group_itx]


class Plate:
    """Class to store one plate.

    Parameters
    ----------
    adata
        The AnnData object containing the plate data.
    treatment_col
        The column name for treatments, by default "Treatment".
    neg_control
        The negative control treatment, by default "DMSO".
        Either `neg_control` or `neg_wells` must be provided.
    neg_wells
        The negative control wells, by default None.
    plate_name
        The name of the plate, by default "Plate".
    """

    def __init__(
        self,
        adata: ad.AnnData,
        treatment_col: str | None = "Treatment",
        neg_wells: list | None = None,
        neg_control: str | None = "DMSO",
        plate_name: str = "Plate",
    ):
        self.adata = adata
        self.id = plate_name
        self.treatment_col = treatment_col
        self.neg_control = neg_control
        if (neg_wells is None) and (self.neg_control is None):
            self.neg_wells = None
        elif (neg_wells is None) and (self.neg_control is not None):
            assert self.treatment_col in adata.obs.columns, (
                f"{treatment_col} not found in adata.obs columns"
            )
            self.neg_wells = None
        else:
            self.neg_wells = neg_wells

    def filter_treatments(
        self,
        treatments: str | list[str],
        treatment_col: str = "Treatment",
        copy: bool = False,
    ) -> "Plate":
        """
        Filters the data by treatments.

        Parameters
        ----------
        treatments
            The treatment(s) to filter by.
        treatment_col
            The column name for treatments, by default "Treatment".
        copy
            Whether to return a copy of the data, by default False.

        Returns
        -------
        Plate
            A Plate object with the filtered data.
        """
        if isinstance(treatments, str):
            treatments = [treatments]
        out = self.adata[self.adata.obs[treatment_col].isin(treatments)]
        if copy:
            out = out.copy()
        return Plate(
            out,
            treatment_col=treatment_col,
            neg_wells=self.neg_wells,
            neg_control=self.neg_control,
            plate_name=self.id,
        )

    def get_controls(self, copy: bool = False, well_key="Well") -> "Plate":
        """
        Gets the control data.

        Parameters
        ----------
        copy
            Whether to return a copy of the data, by default False.
        well_key
            The column name for wells, by default "Well".

        Returns
        -------
        Plate
            A Plate object with the control data.
        """
        ctrls = None
        if self.neg_wells is not None:
            ctrls = self.filter_treatments(self.neg_wells, treatment_col=well_key, copy=copy)
        elif self.neg_control is not None:
            ctrls = self.filter_treatments(
                self.neg_control, treatment_col=self.treatment_col, copy=copy
            )

        if ctrls is not None:
            if len(ctrls) == 0:
                raise ValueError(
                    f"No matching controls found for {self.id}. Check `neg_wells`,"
                    + " `neg_control` and `treatment_col` arguments"
                )
            return ctrls
        raise ValueError("No negative control provided")

    def embed(self, n_pcs: int = 10) -> None:
        """
        Embeds the data using PCA.

        Parameters
        ----------
        n_pcs
            The number of principal components to retain, by default 10.
        """
        if "X_pca" not in self.adata.obsm.keys():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ImplicitModificationWarning)
                sm.pp.scale(self.adata)
                sm.pp.pca(self.adata)
                n_pcs = min(n_pcs, self.adata.obsm["X_pca"].shape[1])
                self.adata.obsm["X_pca"] = self.adata.obsm["X_pca"][:, :n_pcs]

    def get_pca_data(self) -> np.ndarray:
        """
        Gets the PCA data.

        Returns
        -------
        np.ndarray
            The PCA data.
        """
        return self.adata.obsm["X_pca"]

    def iter_group(
        self, col: str = "Well", copy: bool = False, as_plate: bool = False
    ) -> AnnDataIterator:
        """
        Iterates over groups of data.

        Parameters
        ----------
        col
            The column to group by, by default "Well".
        copy
            Whether to return a copy of the data, by default False.
        as_plate
            Whether to return the data as a Plate object, by default False.

        Returns
        -------
        AnnDataIterator
            An iterator over the groups.
        """
        iterator = AnnDataIterator(self.adata, id_col=col, copy=copy, as_plate=as_plate)
        return iterator

    def __repr__(self) -> str:
        return f"Plate {self.id} ({self.adata.shape[0]} x {self.adata.shape[1]})"

    def __str__(self) -> str:
        return self.__repr__()

    def __len__(self) -> int:
        return self.adata.shape[0]


class PlateCollection:
    """Class to store multiple plates.

    Parameters
    ----------
    adatas
        The AnnData object(s) to create the collection from
    treatment_col
        The column for treatments, by default "Treatment".
    neg_control
        The negative control treatment, by default "DMSO".
        Either `neg_control` or `neg_wells` must be provided.
    neg_wells
        The negative control wells, by default None.
    """

    def __init__(
        self,
        adatas: ad.AnnData | list | None,
        treatment_col: str = "Treatment",
        neg_control: str = "DMSO",
        neg_wells: list | None = None,
    ):
        self.plates = []
        if adatas is not None:
            self.add_plate(
                adatas,
                treatment_col=treatment_col,
                neg_control=neg_control,
                neg_wells=neg_wells,
            )

    @classmethod
    def from_adata(
        cls,
        adata: ad.AnnData,
        batch_key: str = "PlateID",
        treatment_col: str | None = None,
        neg_control: str | None = None,
        neg_wells: list | None = None,
    ) -> "PlateCollection":
        """
        Creates a PlateCollection instance from an AnnData object.

        Parameters
        ----------
        adata
            The AnnData object to create the collection from.
        batch_key
            The column to group by, by default "PlateID".
        treatment_col
            The column for treatments, by default `None`.
        neg_control
            The negative control treatment, by default `None`.
        neg_wells
            The negative control wells, by default None.

        Returns
        -------
        PlateCollection
            A new PlateCollection instance.
        """
        instance = cls(adatas=None)
        iterator = AnnDataIterator(adata, batch_key)
        for plate_id, plate_adata in iterator:
            instance.add_plate(
                plate_adata,
                treatment_col=treatment_col,
                neg_control=neg_control,
                neg_wells=neg_wells,
                plate_name=plate_id,
            )
        return instance

    def add_plate(
        self,
        adata: ad.AnnData | Plate | list,
        treatment_col: str | None = "Treatment",
        neg_control: str | None = "DMSO",
        neg_wells: list | None = None,
        plate_name: str = "Plate",
    ) -> None:
        """
        Adds a plate to the collection.

        Parameters
        ----------
        adata
            The plate(s) to add.
        treatment_col
            The column for treatments, by default "Treatment".
        neg_control
            The negative control treatment, by default "DMSO".
            Either `neg_control` or `neg_wells` must be provided.
        neg_wells
            The negative control wells, by default None.
        plate_name
            The name of the plate, by default "Plate".
        """
        if adata is not None:
            if isinstance(adata, ad.AnnData):
                self.plates.append(
                    Plate(
                        adata,
                        treatment_col=treatment_col,
                        neg_control=neg_control,
                        neg_wells=neg_wells,
                        plate_name=plate_name,
                    )
                )
            elif isinstance(adata, Plate):
                self.plates.append(adata)
            elif isinstance(adata, list):
                for plate in adata:
                    self.add_plate(
                        plate,
                        treatment_col=treatment_col,
                        neg_control=neg_control,
                        neg_wells=neg_wells,
                        plate_name=plate_name,
                    )

    def embed(self, n_pcs: int = 10) -> None:
        """
        Embeds the data using PCA for each plate in the collection.

        Parameters
        ----------
        n_pcs
            The number of principal components to retain, by default 10.
        """
        for plate in self.plates:
            plate.embed(n_pcs)

    def __len__(self) -> int:
        return len(self.plates)

    def __next__(self):
        for plate in self.plates:
            yield plate.id, plate

    def __iter__(self):
        return self.__next__()

    def __repr__(self) -> str:
        header = f"PlateCollection with {len(self.plates)} plates"
        return "\n".join([header] + [str(plate) for plate in self.plates])


class MahalanobisKSTest:
    """
    A class to perform Mahalanobis distance and KS test.

    Parameters
    ----------
    ref
        The reference data.
    """

    def __init__(self, ref: np.ndarray):
        self.ref = ref
        self.compute_ref_stats()
        self.ref_dists = self.mahalanobis(self.ref)

    def compute_ref_stats(self) -> None:
        """Computes the mean and covariance matrix of the reference data."""
        self.centr = np.mean(self.ref, axis=0)
        self.cov = np.cov(self.ref, rowvar=False)
        self.covinv = np.linalg.inv(self.cov)

    def mahalanobis(self, X: np.ndarray) -> np.ndarray:
        """
        Computes the Mahalanobis distance for each observation in X.

        Parameters
        ----------
        X
            The data to compute the Mahalanobis distance for.

        Returns
        -------
        np.ndarray
            The Mahalanobis distances.
        """
        return np.apply_along_axis(lambda x: mahalanobis(self.centr, x, self.covinv), axis=1, arr=X)

    def test_ks(self, X: np.ndarray) -> tuple:
        """
        Performs the KS test comparing the Mahalanobis distances of X to the reference data.

        Parameters
        ----------
        X
            The data to test.

        Returns
        -------
        tuple
            The KS statistic and p-value.
        """
        stat = kstest(self.ref_dists, self.mahalanobis(X), alternative="greater", method="asymp")
        return stat.statistic, stat.pvalue


def get_ks_per_group(ref_data: Plate, plate_data: Plate, id_key: str = "Treatment") -> pd.DataFrame:
    """
    Performs the KS test for each group in the plate data against the reference data.

    Parameters
    ----------
    ref_data
        The reference plate data.
    plate_data
        The plate data to test.
    id_key
        The column name to group by, by default "Treatment".

    Returns
    -------
    A DataFrame containing the KS test results for each group.
    """
    ref_group = ref_data.adata.obs[id_key].unique()[0]
    ref_pca = ref_data.get_pca_data()
    tester = MahalanobisKSTest(ref_pca)

    res = []
    for tar_group, tar_data in plate_data.iter_group(id_key, as_plate=True):
        if ref_group == tar_group:
            continue
        tar_pca = tar_data.get_pca_data()
        ks, pval = tester.test_ks(tar_pca)

        res.append(
            [
                ref_data.id,
                ref_group,
                tar_group,
                ks,
                pval,
            ]
        )

    res = pd.DataFrame(
        res,
        columns=[
            "plate",
            "control",
            "treatment",
            "ks_stat",
            "ks_pval",
        ],
    )

    return res


def get_ks_fdr_per_well(
    plates: PlateCollection,
    well_key: str = "Well",
    fdr_threshold: numpy.typing.ArrayLike = (0.05, 0.1),
    progress: bool = False,
) -> pd.DataFrame:
    """
    Performs pairwise KS tests for controls across wells and applies FDR correction.

    Parameters
    ----------
    plates
        The collection of plates to test.
    well_key
        The column name for wells, by default "Well".
    fdr_threshold
        The false discovery rate threshold, by default [0.05, 0.1].
    progress
        Whether to display a progress bar, by default False.

    Returns
    -------
    A DataFrame containing the KS test results with FDR correction.
    """
    ks_over_groups = []
    if progress:
        plates = tqdm(plates)

    for _, plate in plates:
        ctrls = plate.get_controls(well_key=well_key)
        for _, ref_data in ctrls.iter_group(well_key, as_plate=True):
            ref_data.id = plate.id
            res = get_ks_per_group(ref_data, ctrls, well_key)
            ks_over_groups.append(res)

    # Aggregate results from all pairwise comparisons across plates
    ks_over_groups = pd.concat(ks_over_groups)
    ks_over_groups.attrs["ks_threshold"] = {}

    fdr_threshold = np.array(fdr_threshold)
    q_thresh = np.quantile(ks_over_groups["ks_pval"], fdr_threshold, method="inverted_cdf")

    for p, cur_thresh in zip(fdr_threshold, q_thresh, strict=False):
        ks_over_groups.attrs["ks_threshold"][p] = cur_thresh
        ks_over_groups[f"is_significant_{p}"] = ks_over_groups["ks_pval"] < cur_thresh
    return ks_over_groups


def get_ks_per_treatment(
    plates: PlateCollection,
    treatment_key: str = "Treatment",
    well_key="Well",
    progress: bool = False,
) -> pd.DataFrame:
    """
    Performs the KS test for each treatment group in the plate data.

    Parameters
    ----------
    plates
        The collection of plates to test.
    treatment_key
        The column name for treatments, by default "Treatment".
    well_key
        The column name for wells, by default "Well".
    progress
        Whether to display a progress bar, by default False.

    Returns
    -------
    A DataFrame containing the KS test results for each treatment group.
    """
    ks_over_groups = []
    if progress:
        plates = tqdm(plates)

    for _, plate in plates:
        ctrls = plate.get_controls(well_key=well_key)
        res = get_ks_per_group(ctrls, plate, treatment_key)
        ks_over_groups.append(res)

    # Aggregate results from all pairwise comparisons across plates
    ks_over_groups = pd.concat(ks_over_groups)
    ks_over_groups["ks_qval"] = fdrcorrection(ks_over_groups["ks_pval"])[1]
    return ks_over_groups


def get_ks(
    adata: ad.AnnData,
    treatment_key: str = "Treatment",
    well_key: str = "Well",
    control: str = "DMSO",
    control_wells: list | None = None,
    batch_key: str | None = None,
    progress: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Performs the KS test for all treatments and establishes a KS cutoff based on control wells FDR.

    Parameters
    ----------
    adata
        The AnnData object containing the data.
    treatment_key
        The column name for treatments.
    well_key
        The column name for wells, by default "Well". Used for building background distribution.
    control
        The negative control treatment, by default "DMSO".
        Either `neg_control` or `neg_wells` must be provided.
    control_wells
        The negative control wells.
    batch_key
        Additional grouping column, usually empty (single plate) or plate identifier.
    progress
        Whether to display a progress bar.

    Returns
    -------
    A tuple containing the reference KS results and the treatment KS results.
    """
    if batch_key is not None:
        plates = PlateCollection.from_adata(
            adata,
            batch_key=batch_key,
            treatment_col=treatment_key,
            neg_control=control,
            neg_wells=control_wells,
        )
    else:
        plates = PlateCollection(
            adata,
            treatment_col=treatment_key,
            neg_control=control,
            neg_wells=control_wells,
        )

    # Step 1: embed all plates
    plates.embed()

    # Step 2: establish KS cutoff based on control wells FDR
    if progress:
        print("Building negative control p-value distribution")
    ref_ks = get_ks_fdr_per_well(plates, well_key=well_key, progress=progress)

    # Step 3: perform KS test for all treatments
    if progress:
        print("Computing treatment p-values")
    treat_ks = get_ks_per_treatment(
        plates, treatment_key=treatment_key, well_key=well_key, progress=progress
    )

    p_thresholds = ref_ks.attrs["ks_threshold"].items()
    for q, p_ref_thresh in p_thresholds:
        p_thresh = treat_ks.loc[treat_ks["ks_pval"] < p_ref_thresh, "ks_pval"].max()
        treat_ks[f"is_significant_{q}"] = treat_ks["ks_pval"] <= p_thresh

    for df in [ref_ks, treat_ks]:
        if df["plate"].nunique() == 1:
            df.drop(columns="plate", inplace=True)

    return ref_ks, treat_ks
