from __future__ import annotations

import collections
import itertools
import logging

import networkx as nx
import numpy as np
from scipy.sparse.linalg import lsqr
from scipy.sparse import csr_matrix
import bottleneck as bn
from job_pool import JobPool

from .. import helpers
from .base import ProteinGroupColumns
from .sum_and_ibaq import _get_intensities, get_silac_channels
from .peptide_count import _unique_peptide_counts_per_experiment
from . import fastlfq

# imports for typing
from .. import results
from .. import precursor_quant


logger = logging.getLogger(__name__)


class LFQIntensityColumns(ProteinGroupColumns):
    min_peptide_ratios_lfq: int
    stabilize_large_ratios_lfq: bool
    fast_lfq: bool

    def __init__(
        self,
        min_peptide_ratios_lfq: int,
        stabilize_large_ratios_lfq: bool,
        fast_lfq: bool = False,
        fast_lfq_min_neighbors: int = 3,
        fast_lfq_avg_neighbors: int = 6,
        fast_lfq_min_samples: int = 10,
        num_threads: int = 1,
        protein_group_fdr_threshold: float = 0.01,
    ) -> None:
        self.min_peptide_ratios_lfq = min_peptide_ratios_lfq
        self.stabilize_large_ratios_lfq = stabilize_large_ratios_lfq
        self.fast_lfq = fast_lfq
        self.fast_lfq_min_neighbors = fast_lfq_min_neighbors
        self.fast_lfq_avg_neighbors = fast_lfq_avg_neighbors
        self.fast_lfq_min_samples = fast_lfq_min_samples
        self.num_threads = num_threads

        # only used for reporting number of protein groups at the given threshold
        self.protein_group_fdr_threshold = protein_group_fdr_threshold

    def is_valid(self, protein_group_results: results.ProteinGroupResults):
        return (
            len(protein_group_results.experiments) > 1
            and protein_group_results.num_tmt_channels <= 0
        )

    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        silac_channels = get_silac_channels(protein_group_results.num_silac_channels)
        for experiment in protein_group_results.experiments:
            if protein_group_results.num_silac_channels > 0:
                for silac_channel in silac_channels:
                    protein_group_results.append_header(
                        "LFQ Intensity " + silac_channel + " " + experiment
                    )
            else:
                protein_group_results.append_header("LFQ Intensity " + experiment)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()

        logger.info("Doing quantification: MaxLFQ intensity")
        silac_channels = get_silac_channels(protein_group_results.num_silac_channels)
        num_silac_channels = len(silac_channels)

        fast_lfq_graph = None
        if self.fast_lfq:
            peptides_by_sample = collections.defaultdict(set)
            for pgr in protein_group_results:
                for pq in pgr.precursorQuants:
                    peptides_by_sample[pq.experiment].add(pq.peptide)
            fast_lfq_graph = fastlfq.build_graph(peptides_by_sample)
            fast_lfq_graph = fastlfq.prune_graph(
                fast_lfq_graph,
                min_neighbors=self.fast_lfq_min_neighbors,
                avg_neighbors=self.fast_lfq_avg_neighbors,
            )
            # fastlfq.export_to_vis_js(fast_lfq_graph)

        if self.num_threads > 1:
            processingPool = JobPool(
                processes=self.num_threads,
                maxtasksperchild=10,
                max_jobs_queued=self.num_threads * 3,
                write_progress_to_logger=True,
                print_progress_every=100,
                total_jobs=len(protein_group_results),
            )

        allIntensities = list()
        for i, pgr in enumerate(protein_group_results):
            args = [
                pgr.precursorQuants,
                experiment_to_idx_map,
                post_err_prob_cutoff,
                self.min_peptide_ratios_lfq,
                self.stabilize_large_ratios_lfq,
                fast_lfq_graph,
                self.fast_lfq_min_samples,
                num_silac_channels,
            ]
            if self.num_threads > 1:
                processingPool.applyAsync(_get_lfq_intensities, args)
            else:
                if i % 100 == 0:
                    logger.info(f"Processing protein {i}/{len(protein_group_results)}")
                allIntensities.append(_get_lfq_intensities(*args))

        if self.num_threads > 1:
            allIntensities = processingPool.checkPool(printProgressEvery=100)

        protein_group_counts = np.zeros(
            len(experiment_to_idx_map) * max(1, num_silac_channels), dtype="int"
        )
        for i, (pgr, intensities) in enumerate(
            zip(protein_group_results, allIntensities)
        ):
            pgr.extend(intensities)

            if pgr.qValue < self.protein_group_fdr_threshold:
                protein_group_counts += np.array(
                    [1 if intensity > 0 else 0 for intensity in intensities]
                )

        logger.info(
            f"#Protein groups quantified ({self.protein_group_fdr_threshold*100:g}% protein group-level FDR, LFQ):"
        )
        for experiment, num_protein_groups in zip(
            experiment_to_idx_map.keys(),
            helpers.chunks(protein_group_counts, max(1, num_silac_channels)),
        ):
            if num_silac_channels > 0:
                for silac_idx, silac_channel in enumerate(silac_channels):
                    logger.info(
                        f"    {experiment} {silac_channel}: {num_protein_groups[silac_idx]}"
                    )
            else:
                logger.info(f"    {experiment}: {num_protein_groups[0]}")


def _get_lfq_intensities(
    precursor_list: list[precursor_quant.PrecursorQuant],
    experiment_to_idx_map: dict[str, int],
    post_err_prob_cutoff: float,
    min_peptide_ratios_lfq: int = 2,
    stabilize_large_ratios_lfq: bool = False,
    fast_lfq_graph: nx.Graph | None = None,
    fast_lfq_min_samples: int = 10,
    num_silac_channels: int = 0,
) -> list[float]:
    # in case of SILAC, we output LFQ intensites for each of the channels
    # only, not for the experiment (=sum over channels) itself, which is
    # what MQ does as well
    num_experiments = len(experiment_to_idx_map) * max(1, num_silac_channels)

    peptide_intensities, total_intensity = _get_peptide_intensities(
        precursor_list,
        experiment_to_idx_map,
        post_err_prob_cutoff,
        num_silac_channels,
        num_experiments,
    )
    if len(peptide_intensities) == 0:
        return [0.0] * num_experiments

    log_median_peptide_ratios = _get_log_median_peptide_ratios(
        peptide_intensities, min_peptide_ratios_lfq, fast_lfq_graph, fast_lfq_min_samples
    )
    if len(log_median_peptide_ratios) == 0:
        return _get_summed_intensity_with_min_peptides(
            precursor_list,
            experiment_to_idx_map,
            post_err_prob_cutoff,
            num_silac_channels,
            min_peptide_ratios_lfq,
        )

    if stabilize_large_ratios_lfq:
        log_median_peptide_ratios = _apply_large_ratio_stabilization(
            log_median_peptide_ratios,
            precursor_list,
            experiment_to_idx_map,
            post_err_prob_cutoff,
            num_silac_channels,
        )

    matrix, vector = _build_linear_system(log_median_peptide_ratios, num_experiments)
    if len(vector) == 0:
        return [0.0] * num_experiments

    intensities = _solve_linear_system(matrix, vector)
    intensities = _scale_equal_sum(intensities, total_intensity)

    return intensities


def _get_peptide_intensities(
    precursor_list: list[precursor_quant.PrecursorQuant],
    experiment_to_idx_map: dict[str, int],
    post_err_prob_cutoff: float,
    num_silac_channels: int,
    num_experiments: int,
) -> tuple[dict[tuple[str, int], list[float]], float]:
    """
    Collects all precursor intensities per experiment
    """

    def filter_missing_and_unidentified(p: precursor_quant.PrecursorQuant):
        return p.intensity > 0.0 and (
            helpers.is_mbr(p.post_err_prob) or p.post_err_prob <= post_err_prob_cutoff
        )

    precursor_list = filter(filter_missing_and_unidentified, precursor_list)

    # for each (peptide, charge, experiment, fraction) tuple, sort the
    # lowest PEP (= most confident PSM) on top
    def order_by_pep(p: precursor_quant.PrecursorQuant):
        return (
            p.peptide,
            p.charge,
            p.experiment,
            p.fraction,
            -1 * p.intensity,
            p.post_err_prob,
        )

    precursor_list = sorted(precursor_list, key=order_by_pep)

    peptide_intensities = collections.defaultdict(lambda: [0.0] * num_experiments)
    total_intensity = 0.0
    prev_exp_frac = (None, None)
    prev_precursor = (None, None)
    for precursor in precursor_list:
        exp_idx = experiment_to_idx_map[precursor.experiment]
        # for each (peptide, charge, experiment, fraction) tuple, only
        # use the intensity of the PSM with the lowest PEP (= most
        # confident PSM)
        curr_precursor = (precursor.peptide, precursor.charge)
        curr_exp_frac = (precursor.experiment, precursor.fraction)
        if prev_exp_frac != curr_exp_frac or prev_precursor != curr_precursor:
            if num_silac_channels > 0:
                for silac_idx, silac_intensity in enumerate(precursor.silac_intensities):
                    silac_exp_idx = exp_idx * num_silac_channels + silac_idx
                    peptide_intensities[curr_precursor][silac_exp_idx] += silac_intensity
                    total_intensity += silac_intensity
            else:
                peptide_intensities[curr_precursor][exp_idx] += precursor.intensity
                total_intensity += precursor.intensity
            prev_exp_frac = curr_exp_frac
            prev_precursor = curr_precursor
    return peptide_intensities, total_intensity


def _get_summed_intensity_with_min_peptides(
    precursor_list: list[precursor_quant.PrecursorQuant],
    experiment_to_idx_map: dict[str, int],
    post_err_prob_cutoff: float,
    num_silac_channels: int,
    min_peptide_ratios_lfq: int,
) -> dict[tuple[int, int], float]:
    """
    Return summed intensity for proteins that only occur in a single
    sample or where multiple peptides are detected for a protein, but
    each peptide was only identified in a single sample
    (see also https://github.com/kusterlab/picked_group_fdr/issues/16)
    """
    summed_intensities = _get_intensities(
        precursor_list,
        experiment_to_idx_map,
        post_err_prob_cutoff,
        num_silac_channels,
        remove_silac_summed_intensity_columns=True,
    )

    peptide_counts = _unique_peptide_counts_per_experiment(
        precursor_list, experiment_to_idx_map, post_err_prob_cutoff
    )
    peptide_counts = np.repeat(peptide_counts, max([1, num_silac_channels]))
    return [
        summed_intensity if peptide_count >= min_peptide_ratios_lfq else 0.0
        for summed_intensity, peptide_count in zip(
            summed_intensities,
            peptide_counts,
        )
    ]


def _apply_large_ratio_stabilization(
    log_median_peptide_ratios: dict[tuple[int, int], float],
    precursor_list: list[precursor_quant.PrecursorQuant],
    experiment_to_idx_map: dict[str, int],
    post_err_prob_cutoff: float,
    num_silac_channels: int,
) -> dict[tuple[int, int], float]:
    """Stabilizes log median peptide ratios based on peptide count differences.

    Applies ratio stabilization to account for situations where large differences
    in peptide counts between experiments could lead to unreliable ratio estimates.
    When one experiment has substantially more peptide identifications than another,
    the median peptide ratio is either replaced or blended with the intensity-based
    ratio.

    Args:
        log_median_peptide_ratios: Dictionary mapping experiment index pairs (i, j)
            to their log-transformed median peptide ratios.
        precursor_list: List of PrecursorQuant objects containing peptide MS/MS
            identifications with intensities.
        experiment_to_idx_map: Dictionary mapping experiment names to their indices.
        post_err_prob_cutoff: Posterior error probability threshold for filtering
            identifications. Only PSMs with PEP <= cutoff are included.
        num_silac_channels: Number of SILAC heavy/light channels. 0 for label-free
            quantification. When > 0, intensities are processed per channel.

    Returns:
        Dictionary with the same structure as log_median_peptide_ratios but with
        stabilized log ratios. Experiments with zero peptide counts are skipped,
        and ratios not present in the input dictionary are not added.

    Notes:
        Stabilization logic is applied based on the peptide count ratio between
        experiments:
        - Ratio > 5: Replace log median ratio with log(summed_intensity_1 /
          summed_intensity_2)
        - 2.5 < Ratio <= 5: Use weighted average of intensity-based log ratio
          and original log median ratio, where weight increases with the
          peptide count ratio
        - Ratio <= 2.5: Keep original log median ratio unchanged

        The peptide count ratio is computed as max(pc1, pc2) / min(pc1, pc2).
    """
    summed_intensities = _get_intensities(
        precursor_list,
        experiment_to_idx_map,
        post_err_prob_cutoff,
        num_silac_channels,
        remove_silac_summed_intensity_columns=True,
    )

    peptide_counts = _unique_peptide_counts_per_experiment(
        precursor_list, experiment_to_idx_map, post_err_prob_cutoff
    )
    peptide_counts = np.repeat(peptide_counts, max([1, num_silac_channels]))
    for (i, j), (si1, si2), (pc1, pc2) in zip(
        itertools.combinations(range(len(summed_intensities)), 2),
        itertools.combinations(summed_intensities, 2),
        itertools.combinations(peptide_counts, 2),
    ):
        if pc1 == 0 or pc2 == 0 or (i, j) not in log_median_peptide_ratios:
            continue

        peptide_count_ratio = _get_max_ratio(pc1, pc2)
        if peptide_count_ratio > 5:
            log_median_peptide_ratios[(i, j)] = np.log(si1 / si2)
        elif peptide_count_ratio > 2.5:
            w = (peptide_count_ratio - 2.5) / 2.5
            log_median_peptide_ratios[(i, j)] = (
                w * np.log(si1 / si2) + (1 - w) * log_median_peptide_ratios[(i, j)]
            )

    return log_median_peptide_ratios


def _get_log_median_peptide_ratios(
    peptide_intensities: dict[tuple[str, int], list[float]],
    min_peptide_ratios_lfq: int,
    fast_lfq_graph: nx.Graph | None = None,
    fast_lfq_min_samples: int = 10,
) -> dict[tuple[int, int], float]:
    """
    :param peptide_intensities: rows are peptides, columns are experiments
    :param min_peptide_ratios_lfq: minimum valid ratios needed to perform LFQ
    returns: dictionary of (sample_i, sample_j) -> log(median(ratios))
    """
    intensity_matrix = np.array(list(peptide_intensities.values()))

    nonzeros_per_column = np.count_nonzero(intensity_matrix, axis=0)
    valid_columns = np.argwhere(nonzeros_per_column >= min_peptide_ratios_lfq).flatten()

    valid_counts_matrix = (intensity_matrix[:, valid_columns] > 0).astype(int)
    valid_vals = np.dot(valid_counts_matrix.T, valid_counts_matrix)

    # replace zeroes by NaNs to automatically filter ratios with one missing value
    intensity_matrix[intensity_matrix == 0] = np.nan
    cols = [(idx, i, intensity_matrix[:, i]) for idx, i in enumerate(valid_columns)]
    peptide_ratios, experiment_pairs = list(), list()
    # this for loop is the reason for quadratic runtime increase for more samples
    # (30 seconds per protein for 2400 samples and 600 peptides)
    for (idx_i, i, col_i), (idx_j, j, col_j) in itertools.combinations(cols, 2):
        if (
            fast_lfq_graph is not None
            and len(valid_columns) >= fast_lfq_min_samples
            and not fast_lfq_graph.has_edge(i, j)
        ):
            continue

        if valid_vals[idx_i, idx_j] < min_peptide_ratios_lfq:
            continue

        ratios = col_i / col_j
        # vectorization of computing medians is 25% faster than computing
        # them in a loop here but requires more memory
        peptide_ratios.append(bn.nanmedian(ratios))
        experiment_pairs.append((i, j))

    return dict(zip(experiment_pairs, np.log(peptide_ratios)))


def _get_max_ratio(pc1: float, pc2: float) -> float:
    assert pc1 > 0 and pc2 > 0
    if pc1 < pc2:
        return pc2 / pc1
    else:
        return pc1 / pc2


def _build_linear_system(
    log_median_peptide_ratios: dict[tuple[int, int], float], num_experiments: int
) -> tuple[np.ndarray, np.ndarray]:
    num_ratios = len(log_median_peptide_ratios)
    rows, cols, vals = list(), list(), list()
    seen_columns = set()
    for idx, (i, j) in enumerate(log_median_peptide_ratios.keys()):
        rows.append(idx)
        cols.append(i)
        vals.append(1.0)
        seen_columns.add(i)

        rows.append(idx)
        cols.append(j)
        vals.append(-1.0)
        seen_columns.add(j)

    # anchor the (log) average and experiments without ratios to 0,
    # otherwise the matrix is ill-conditioned and least squares sometimes
    # does not converge
    num_non_zero_columns = len(seen_columns)
    zero_columns = [x for x in range(num_experiments) if x not in seen_columns]
    num_zero_columns = len(zero_columns)

    rows.extend([num_ratios] * num_non_zero_columns)
    cols.extend(list(seen_columns))
    vals.extend([1.0] * num_non_zero_columns)

    rows.extend(range(num_ratios + 1, num_ratios + 1 + num_zero_columns))
    cols.extend(zero_columns)
    vals.extend([1.0] * num_zero_columns)

    matrix = csr_matrix((vals, (rows, cols)), shape=(max(rows) + 1, max(cols) + 1))
    vector = np.array(
        list(log_median_peptide_ratios.values()) + [0.0] * (num_zero_columns + 1)
    )

    # np.savetxt("matrix.csv", matrix, delimiter="\t")
    # np.savetxt("vector.csv", vector, delimiter="\t")

    return matrix, vector


def _solve_linear_system(matrix: csr_matrix, vector: np.ndarray) -> np.ndarray:
    intensities = np.exp(lsqr(matrix, vector)[0])

    zero_columns = np.array(matrix.astype(bool).sum(axis=0))[0] == 1
    intensities[zero_columns] = 0.0
    return intensities


def _scale_equal_sum(intensities: np.ndarray, total_intensity: float) -> list[float]:
    """
    scale intensities such that the total sum stays the same as it was
    before

    :param intensities: array of intensities per sample
    :param total_intensity: sum of intensities prior to method application
    :returns: array of intensities scaled such that sum equals
              total_intensity
    """
    intensity_sum = np.sum(intensities)
    if intensity_sum > 0:
        return ((total_intensity / intensity_sum) * intensities).tolist()
    else:
        return intensities.tolist()
