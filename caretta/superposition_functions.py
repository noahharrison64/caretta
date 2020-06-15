import numba as nb
import numpy as np
from caretta import dynamic_time_warping as dtw, score_functions
from geometricus import moment_utility, utility
import prody as pd

GAP_OPEN = 0
GAP_EXTEND = 0

"""
Provides pairwise superposition functions to use in Caretta
Each takes coords_1, coords_2, parameters, score_function as input
    parameters is a dict with gap_open_penalty, gap_extend_penalty, and other specific parameters as keys

returns score, superposed_coords_1, superposed_coords_2
"""


def dtw_svd_superpose_function(coords_1, coords_2, parameters: dict, score_function=score_functions.get_caretta_score):
    """
    Assumes coords_1 and coords_2 are already in a well-superposed state,
    runs DTW alignment and then superposes with Kabsch on the aligning positions
    """
    score_matrix = score_functions.make_score_matrix(coords_1, coords_2,
                                                     score_function,
                                                     normalized=False)
    _, coords_1, coords_2, common_coords_1, common_coords_2 = _align_and_superpose(coords_1, coords_2, score_matrix, parameters["gap_open_penalty"],
                                                                                   parameters["gap_extend_penalty"])
    return score_functions.get_total_score(common_coords_1, common_coords_2, score_function, False), coords_1, coords_2


def signal_superpose_function(coords_1, coords_2, parameters, score_function=score_functions.get_signal_score):
    """
    Makes initial superposition of coordinates using DTW alignment of overlapping signals
    A signal is a vector of euclidean distances of first (or last) coordinate to all others in a 30-residue stretch
    """
    score_first, c1_first, c2_first = _signal_superpose_index(0, coords_1, coords_2, score_function, parameters["gap_open_penalty"],
                                                              parameters["gap_extend_penalty"])
    score_last, c1_last, c2_last = _signal_superpose_index(-1, coords_1, coords_2, score_function, parameters["gap_open_penalty"],
                                                           parameters["gap_extend_penalty"])
    if score_first > score_last:
        return score_first, c1_first, c2_first
    else:
        return score_last, c1_last, c2_last


def signal_svd_superpose_function(coords_1, coords_2, parameters, score_function=score_functions.get_caretta_score):
    """
    Uses signal_superpose followed by dtw_svd_superpose
    """
    _, coords_1, coords_2 = signal_superpose_function(coords_1, coords_2, parameters)
    return dtw_svd_superpose_function(coords_1, coords_2, parameters, score_function)


@nb.njit
def get_moments_kmer(coords, kmer_size=30, scale=True):
    res = np.zeros((coords.shape[0], moment_utility.NUM_MOMENTS))
    half = kmer_size // 2
    for i in nb.prange(res.shape[0]):
        res[i, :] = moment_utility.get_second_order_moments(coords[max(0, i - half): min(coords.shape[0], i + half)])
    if scale:
        return np.log1p(res)
    else:
        return res


def get_moments_radius(coords, radius=4, scale=True):
    indices = []
    kd_tree = pd.KDTree(coords)
    for i in range(coords.shape[0]):
        kd_tree.search(center=coords[i], radius=radius)
        indices.append(kd_tree.getIndices())
    moments = np.zeros((len(indices), moment_utility.NUM_MOMENTS))
    for i in range(len(indices)):
        moments[i] = moment_utility.get_second_order_moments(coords[indices[i]])
    if scale:
        return np.log1p(moments)
    else:
        return moments


def moment_superpose_function(coords_1, coords_2, parameters, score_function=score_functions.get_caretta_score):
    """
    Uses 4 rotation/translation invariant moments for each 5-mer to run DTW
    """
    if parameters["split_type"] == "kmer":
        moments_1, moments_2 = get_moments_kmer(coords_1, parameters["split_size"], parameters["scale"]), get_moments_kmer(coords_2,
                                                                                                                           parameters["split_size"],
                                                                                                                           parameters["scale"])
    else:
        moments_1, moments_2 = get_moments_radius(coords_1, parameters["split_size"], parameters["scale"]), get_moments_radius(coords_2, parameters[
            "split_size"], parameters["scale"])
    score_matrix = score_functions.make_score_matrix(moments_1, moments_2, score_function, normalized=False)
    score, coords_1, coords_2, _, _ = _align_and_superpose(coords_1, coords_2, score_matrix, parameters["gap_open_penalty"],
                                                           parameters["gap_extend_penalty"])
    return score, coords_1, coords_2


def moment_superpose_function_both(coords_1, coords_2, parameters, score_function=score_functions.get_caretta_score):
    """
    Uses 4 rotation/translation invariant moments for each 5-mer to run DTW
    """
    moments_1_k, moments_2_k = get_moments_kmer(coords_1, parameters["kmer_size"]), get_moments_kmer(coords_2, parameters["kmer_size"])
    moments_1_r, moments_2_r = get_moments_radius(coords_1, parameters["radius"]), get_moments_radius(coords_2, parameters["radius"])

    moments_1 = moments_1_k + moments_1_r
    moments_2 = moments_2_k + moments_2_r
    score_matrix = score_functions.make_score_matrix(moments_1, moments_2, score_function, normalized=False)
    score, coords_1, coords_2, _, _ = _align_and_superpose(coords_1, coords_2, score_matrix, parameters["gap_open_penalty"],
                                                           parameters["gap_extend_penalty"])
    return score, coords_1, coords_2


def moment_svd_superpose_function(coords_1, coords_2, parameters, score_function=score_functions.get_caretta_score):
    """
    Uses moment_superpose followed by dtw_svd_superpose
    """
    _, coords_1, coords_2 = moment_superpose_function(coords_1, coords_2, parameters)
    return dtw_svd_superpose_function(coords_1, coords_2, parameters, score_function)


@nb.njit
def _align_and_superpose(coords_1, coords_2, score_matrix, gap_open_penalty, gap_extend_penalty):
    """
    Runs DTW on a score matrix and Kabsch superposition on resulting alignment
    """
    dtw_aln_array_1, dtw_aln_array_2, score = dtw.dtw_align(score_matrix, gap_open_penalty, gap_extend_penalty)
    pos_1, pos_2 = score_functions.get_common_positions(dtw_aln_array_1, dtw_aln_array_2)
    common_coords_1, common_coords_2 = coords_1[pos_1], coords_2[pos_2]
    coords_1, coords_2, common_coords_2 = paired_svd_superpose_with_subset(coords_1, coords_2, common_coords_1, common_coords_2)
    return score, coords_1, coords_2, common_coords_1, common_coords_2


@nb.njit
def paired_svd_superpose(coords_1: np.ndarray, coords_2: np.ndarray):
    """
    Superpose paired coordinates on each other using Kabsch superposition (SVD)

    Parameters
    ----------
    coords_1
        numpy array of coordinate data for the first protein; shape = (n, 3)
    coords_2
        numpy array of corresponding coordinate data for the second protein; shape = (n, 3)

    Returns
    -------
    rotation matrix, translation matrix for optimal superposition
    """
    centroid_1, centroid_2 = utility.nb_mean_axis_0(coords_1), utility.nb_mean_axis_0(coords_2)
    coords_1_c, coords_2_c = coords_1 - centroid_1, coords_2 - centroid_2
    correlation_matrix = np.dot(coords_2_c.T, coords_1_c)
    u, s, v = np.linalg.svd(correlation_matrix)
    reflect = np.linalg.det(u) * np.linalg.det(v) < 0
    if reflect:
        s[-1] = -s[-1]
        u[:, -1] = -u[:, -1]
    rotation_matrix = np.dot(u, v)
    translation_matrix = centroid_1 - np.dot(centroid_2, rotation_matrix)
    return rotation_matrix.astype(np.float64), translation_matrix.astype(np.float64)


@nb.njit
def paired_svd_superpose_with_subset(coords_1, coords_2, common_coords_1, common_coords_2):
    """
    Superpose two sets of un-aligned coordinates using smaller subsets of aligned coordinates

    Parameters
    ----------
    coords_1
    coords_2
    common_coords_1
    common_coords_2

    Returns
    -------
    superposed coord_1, superposed coords_2, superposed common_coords_2
    """
    rot, tran = paired_svd_superpose(common_coords_1, common_coords_2)
    coords_1 = coords_1 - utility.nb_mean_axis_0(common_coords_1)
    coords_2 = np.dot(coords_2 - utility.nb_mean_axis_0(common_coords_2), rot)
    common_coords_2_rot = apply_rotran(common_coords_2, rot, tran)
    return coords_1, coords_2, common_coords_2_rot


def _signal_superpose_index(index, coords_1, coords_2, score_function, gap_open_penalty=0., gap_extend_penalty=0., size=30, overlap=1):
    """
    Makes initial superposition using DTW alignment of overlapping signals
    A signal is a vector of euclidean distances of first (or last) coordinate to all others in a 30-residue stretch
    """

    def _make_signal_index(coords, idx):
        centroid = coords[idx]
        distances = np.zeros(coords.shape[0])
        for c in range(coords.shape[0]):
            distances[c] = np.sqrt(np.sum((coords[c] - centroid) ** 2, axis=-1))
        return distances

    signals_1 = np.zeros(((coords_1.shape[0] - size) // overlap, size))
    signals_2 = np.zeros(((coords_2.shape[0] - size) // overlap, size))
    middles_1 = np.zeros((signals_1.shape[0], coords_1.shape[1]))
    middles_2 = np.zeros((signals_2.shape[0], coords_2.shape[1]))
    if index == -1:
        index = size - 1
    for x, i in enumerate(range(0, signals_1.shape[0] * overlap, overlap)):
        signals_1[x] = _make_signal_index(coords_1[i: i + size], index)
        middles_1[x] = coords_1[i + index]
    for x, i in enumerate(range(0, signals_2.shape[0] * overlap, overlap)):
        signals_2[x] = _make_signal_index(coords_2[i: i + size], index)
        middles_2[x] = coords_2[i + index]
    score_matrix = score_functions.make_score_matrix(signals_1, signals_2, score_function, normalized=False)
    dtw_1, dtw_2, score = dtw.dtw_align(score_matrix, gap_open_penalty, gap_extend_penalty)
    pos_1, pos_2 = score_functions.get_common_positions(dtw_1, dtw_2)
    aln_coords_1 = np.zeros((len(pos_1), coords_1.shape[1]))
    aln_coords_2 = np.zeros((len(pos_2), coords_2.shape[1]))
    for i, (p1, p2) in enumerate(zip(pos_1, pos_2)):
        aln_coords_1[i] = middles_1[p1]
        aln_coords_2[i] = middles_2[p2]
    coords_1, coords_2, _ = paired_svd_superpose_with_subset(coords_1, coords_2, aln_coords_1, aln_coords_2)
    return score, coords_1, coords_2


@nb.njit
# @numba_cc.export('svd_superimpose', '(f64[:], f64[:])')
def svd_superimpose(coords_1: np.ndarray, coords_2: np.ndarray):
    """
    Superimpose paired coordinates on each other using svd

    Parameters
    ----------
    coords_1
        numpy array of coordinate data for the first protein; shape = (n, 3)
    coords_2
        numpy array of corresponding coordinate data for the second protein; shape = (n, 3)

    Returns
    -------
    rotation matrix, translation matrix for optimal superposition
    """
    centroid_1, centroid_2 = utility.nb_mean_axis_0(coords_1), utility.nb_mean_axis_0(coords_2)
    coords_1_c, coords_2_c = coords_1 - centroid_1, coords_2 - centroid_2
    correlation_matrix = np.dot(coords_2_c.T, coords_1_c)
    u, s, v = np.linalg.svd(correlation_matrix)
    reflect = np.linalg.det(u) * np.linalg.det(v) < 0
    if reflect:
        s[-1] = -s[-1]
        u[:, -1] = -u[:, -1]
    rotation_matrix = np.dot(u, v)
    translation_matrix = centroid_1 - np.dot(centroid_2, rotation_matrix)
    return rotation_matrix.astype(np.float64), translation_matrix.astype(np.float64)


@nb.njit
# @numba_cc.export('apply_rotran', '(f64[:], f64[:], f64[:])')
def apply_rotran(coords: np.ndarray, rotation_matrix: np.ndarray, translation_matrix: np.ndarray) -> np.ndarray:
    """
    Applies a rotation and translation matrix onto coordinates

    Parameters
    ----------
    coords
    rotation_matrix
    translation_matrix

    Returns
    -------
    transformed coordinates
    """
    return np.dot(coords, rotation_matrix) + translation_matrix


# @numba_cc.export('superpose_with_pos', '(f64[:], f64[:], f64[:], f64[:])')
@nb.njit
def superpose_with_pos(coords_1, coords_2, common_coords_1, common_coords_2):
    """
    Superpose two sets of un-aligned coordinates using smaller subsets of aligned coordinates

    Parameters
    ----------
    coords_1
    coords_2
    common_coords_1
    common_coords_2

    Returns
    -------
    superposed coord_1, superposed coords_2, superposed common_coords_2
    """
    rot, tran = svd_superimpose(common_coords_1, common_coords_2)
    coords_1 = coords_1 - utility.nb_mean_axis_0(common_coords_1)
    coords_2 = np.dot(coords_2 - utility.nb_mean_axis_0(common_coords_2), rot)
    common_coords_2_rot = apply_rotran(common_coords_2, rot, tran)
    return coords_1, coords_2, common_coords_2_rot