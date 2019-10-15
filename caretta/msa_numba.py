import typing
from dataclasses import dataclass, field
from pathlib import Path

import numba as nb
import numpy as np
import prody as pd

from caretta import feature_extraction
from caretta import neighbor_joining as nj
from caretta import psa_numba as psa
from caretta import rmsd_calculations, helper
from setup import numba_cc


@nb.njit
@numba_cc.export('get_common_coordinates', '(f64[:], f64[:], i64[:], i64[:], i64)')
def get_common_coordinates(coords_1, coords_2, aln_1, aln_2, gap=-1):
    assert aln_1.shape == aln_2.shape
    pos_1, pos_2 = helper.get_common_positions(aln_1, aln_2, gap)
    return coords_1[pos_1], coords_2[pos_2]


@nb.njit(parallel=True)
@numba_cc.export('make_pairwise_dtw_score_matrix', '(f64[:], i8[:], i64[:], f64, f64, f64, f64, f64)')
def make_pairwise_dtw_score_matrix(coords_array, secondary_array, lengths_array, gamma,
                                   gap_open_penalty: float, gap_extend_penalty: float,
                                   gap_open_sec, gap_extend_sec):
    pairwise_matrix = np.zeros((coords_array.shape[0], coords_array.shape[0]))
    for i in nb.prange(pairwise_matrix.shape[0] - 1):
        for j in range(i + 1, pairwise_matrix.shape[1]):
            dtw_aln_1, dtw_aln_2, score = psa.get_pairwise_alignment(coords_array[i, :lengths_array[i]], coords_array[j, :lengths_array[j]],
                                                                     secondary_array[i, :lengths_array[i]], secondary_array[j, :lengths_array[j]],
                                                                     gamma,
                                                                     gap_open_sec=gap_open_sec,
                                                                     gap_extend_sec=gap_extend_sec,
                                                                     gap_open_penalty=gap_open_penalty,
                                                                     gap_extend_penalty=gap_extend_penalty)
            common_coords_1, common_coords_2 = get_common_coordinates(coords_array[i, :lengths_array[i]],
                                                                      coords_array[j, :lengths_array[j]],
                                                                      dtw_aln_1, dtw_aln_2)
            rotation_matrix, translation_matrix = rmsd_calculations.svd_superimpose(common_coords_1[:, :3], common_coords_2[:, :3])
            common_coords_2[:, :3] = rmsd_calculations.apply_rotran(common_coords_2[:, :3], rotation_matrix, translation_matrix)
            score = rmsd_calculations.get_caretta_score(common_coords_1, common_coords_2, gamma, True)
            pairwise_matrix[i, j] = -score
    pairwise_matrix += pairwise_matrix.T
    return pairwise_matrix


@nb.njit
@numba_cc.export('_get_alignment_data', '(f64[:], f64[:], i8[:], i8[:], f64, f64, f64, f64, f64)')
def _get_alignment_data(coords_1, coords_2, secondary_1, secondary_2, gamma,
                        gap_open_sec, gap_extend_sec,
                        gap_open_penalty: float, gap_extend_penalty: float):
    dtw_aln_1, dtw_aln_2, _ = psa.get_pairwise_alignment(
        coords_1, coords_2,
        secondary_1, secondary_2,
        gamma,
        gap_open_sec=gap_open_sec,
        gap_extend_sec=gap_extend_sec,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty)
    pos_1, pos_2 = helper.get_common_positions(dtw_aln_1, dtw_aln_2, -1)
    coords_1[:, :3], coords_2[:, :3], _ = rmsd_calculations.superpose_with_pos(coords_1[:, :3], coords_2[:, :3],
                                                                               coords_1[pos_1][:, :3], coords_2[pos_2][:, :3])
    aln_coords_1 = helper.get_aligned_data(dtw_aln_1, coords_1, -1)
    aln_coords_2 = helper.get_aligned_data(dtw_aln_2, coords_2, -1)
    aln_sec_1 = helper.get_aligned_string_data(dtw_aln_1, secondary_1, -1)
    aln_sec_2 = helper.get_aligned_string_data(dtw_aln_2, secondary_2, -1)
    return aln_coords_1, aln_coords_2, aln_sec_1, aln_sec_2, dtw_aln_1, dtw_aln_2


@nb.njit
@numba_cc.export('get_mean_coords_extra', 'f64[:](f64[:], f64[:])')
def get_mean_coords_extra(aln_coords_1: np.ndarray, aln_coords_2: np.ndarray) -> np.ndarray:
    """
    Mean of two coordinate sets (of the same shape)

    Parameters
    ----------
    aln_coords_1
    aln_coords_2

    Returns
    -------
    mean_coords
    """
    mean_coords = np.zeros(aln_coords_1.shape)
    for i in range(aln_coords_1.shape[0]):
        mean_coords[i, :-1] = np.array([np.nanmean(np.array([aln_coords_1[i, x], aln_coords_2[i, x]])) for x in range(aln_coords_1.shape[1] - 1)])
        if not np.isnan(aln_coords_1[i, 0]):
            mean_coords[i, -1] += aln_coords_1[i, -1]
        if not np.isnan(aln_coords_2[i, 0]):
            mean_coords[i, -1] += aln_coords_2[i, -1]
    return mean_coords


@nb.njit
@numba_cc.export('get_mean_secondary', 'i8[:](i8[:], i8[:])')
def get_mean_secondary(aln_sec_1: np.ndarray, aln_sec_2: np.ndarray, gap=0) -> np.ndarray:
    """
    Mean of two coordinate sets (of the same shape)

    Parameters
    ----------
    aln_sec_1
    aln_sec_2
    gap

    Returns
    -------
    mean_sec
    """
    mean_sec = np.zeros(aln_sec_1.shape, dtype=aln_sec_1.dtype)
    for i in range(aln_sec_1.shape[0]):
        if aln_sec_1[i] == aln_sec_2[i]:
            mean_sec[i] = aln_sec_1[i]
        else:
            if aln_sec_1[i] != gap:
                mean_sec[i] = aln_sec_1[i]
            elif aln_sec_2[i] != gap:
                mean_sec[i] = aln_sec_2[i]
    return mean_sec


@dataclass(eq=False)
class Structure:
    name: str
    pdb_file: typing.Union[str, Path, None]
    sequence: typing.Union[str, None]
    secondary: np.ndarray = field(repr=False)
    features: typing.Union[np.ndarray, None] = field(repr=False)
    coords: np.ndarray = field(repr=False)

    @classmethod
    def from_pdb_file(cls, pdb_file: typing.Union[str, Path], dssp_dir="caretta_tmp",
                      extract_all_features=True, force_overwrite=False):
        pdb_name = helper.get_file_parts(pdb_file)[1]
        pdb = pd.parsePDB(str(pdb_file))
        alpha_indices = helper.get_alpha_indices(pdb)
        sequence = pdb[alpha_indices].getSequence()
        coordinates = pdb[alpha_indices].getCoords().astype(np.float64)
        only_dssp = (not extract_all_features)
        features = feature_extraction.get_features(str(pdb_file), str(dssp_dir), only_dssp=only_dssp, force_overwrite=force_overwrite)
        return cls(pdb_name, pdb_file, sequence, helper.secondary_to_array(features["secondary"]), features, coordinates)


# class Structure:
#     def __init__(self, name, pdb_file, sequence, secondary, features, coords, add_column=True, consensus_weight=1.):
#         """
#         Makes a Structure object
#
#         Parameters
#         ----------
#         name
#         sequence
#         secondary
#         features
#             dict of features
#         coords
#             array of alpha carbon coordinates
#         add_column
#             adds a column for consensus weight if True
#         """
#         self.name = name
#         self.pdb_file = pdb_file
#         self.sequence = sequence
#         self.secondary = secondary
#         self.features = features
#         if add_column:
#             add = np.zeros((coords.shape[0], 1)) + consensus_weight
#             self.coords = np.hstack((coords, add))
#         else:
#             self.coords = coords


class StructureMultiple:
    """
    Class for multiple structure alignment
    """

    def __init__(self, pdb_files, dssp_dir, consensus_weight=1., num_threads=20, extract_all_features=True, force_overwrite=False):
        pdbs = [pd.parsePDB(filename) for filename in pdb_files]
        alpha_indices = [helper.get_alpha_indices(pdb) for pdb in pdbs]
        sequences = [pdbs[i][alpha_indices[i]].getSequence() for i in range(len(pdbs))]
        coordinates = [np.hstack((pdbs[i][alpha_indices[i]].getCoords().astype(np.float64), np.zeros((len(sequences[i]), 1)) + consensus_weight))
                       for i in range(len(pdbs))]
        only_dssp = (not extract_all_features)
        features = feature_extraction.get_features_multiple(pdb_files, str(dssp_dir),
                                                            num_threads=num_threads, only_dssp=only_dssp, force_overwrite=force_overwrite)
        self.structures = []
        for i in range(len(pdbs)):
            pdb_name = helper.get_file_parts(pdb_files[i])[1]
            self.structures.append(Structure(pdb_name,
                                             pdb_files[i],
                                             sequences[i],
                                             helper.secondary_to_array(features[i]["secondary"]),
                                             features[i],
                                             coordinates[i]))
        self.lengths_array = np.array([len(s.sequence) for s in self.structures])
        self.max_length = np.max(self.lengths_array)
        self.coords_array = np.zeros((len(self.structures), self.max_length, self.structures[0].coords.shape[1]))
        self.secondary_array = np.zeros((len(self.structures), self.max_length))
        for i in range(len(self.structures)):
            self.coords_array[i, :self.lengths_array[i]] = self.structures[i].coords
            self.secondary_array[i, :self.lengths_array[i]] = self.structures[i].secondary
        self.final_structures = []
        self.tree = None
        self.branch_lengths = None
        self.alignment = None

    def align(self, gamma, gap_open_sec, gap_extend_sec, gap_open_penalty, gap_extend_penalty, pw_matrix=None) -> dict:
        print("Aligning...")
        if len(self.structures) == 2:
            dtw_1, dtw_2, _ = psa.get_pairwise_alignment(self.coords_array[0, :self.lengths_array[0]],
                                                         self.coords_array[1, :self.lengths_array[1]],
                                                         self.secondary_array[0, :self.lengths_array[0]],
                                                         self.secondary_array[1, :self.lengths_array[1]],
                                                         gamma,
                                                         gap_open_sec=gap_open_sec,
                                                         gap_extend_sec=gap_extend_sec,
                                                         gap_open_penalty=gap_open_penalty,
                                                         gap_extend_penalty=gap_extend_penalty)
            self.alignment = {self.structures[0].name: "".join([self.structures[0].sequence[i] if i != -1 else '-' for i in dtw_1]),
                              self.structures[1].name: "".join([self.structures[1].sequence[i] if i != -1 else '-' for i in dtw_2])}
            return self.alignment

        if pw_matrix is None:
            pw_matrix = make_pairwise_dtw_score_matrix(self.coords_array,
                                                       self.secondary_array,
                                                       self.lengths_array,
                                                       gamma,
                                                       gap_open_sec, gap_extend_sec,
                                                       gap_open_penalty, gap_extend_penalty)

        print("Pairwise score matrix calculation done")

        tree, branch_lengths = nj.neighbor_joining(pw_matrix)
        self.tree = tree
        self.branch_lengths = branch_lengths
        self.final_structures = [s for s in self.structures]
        msa_alignments = {s.name: {s.name: s.sequence} for s in self.structures}

        print("Neighbor joining tree constructed")

        def make_intermediate_node(n1, n2, n_int):
            name_1, name_2 = self.final_structures[n1].name, self.final_structures[n2].name
            name_int = f"int-{n_int}"
            n1_coords = self.final_structures[n1].coords
            n1_coords[:, -1] *= len(msa_alignments[name_2])
            n1_coords[:, -1] /= 2.
            n2_coords = self.final_structures[n2].coords
            n2_coords[:, -1] *= len(msa_alignments[name_1])
            n2_coords[:, -1] /= 2.
            aln_coords_1, aln_coords_2, aln_sec_1, aln_sec_2, dtw_aln_1, dtw_aln_2 = _get_alignment_data(n1_coords,
                                                                                                         n2_coords,
                                                                                                         self.final_structures[
                                                                                                             n1].secondary,
                                                                                                         self.final_structures[
                                                                                                             n2].secondary,
                                                                                                         gamma,
                                                                                                         gap_open_sec=gap_open_sec,
                                                                                                         gap_extend_sec=gap_extend_sec,
                                                                                                         gap_open_penalty=gap_open_penalty,
                                                                                                         gap_extend_penalty=gap_extend_penalty)
            aln_coords_1[:, -1] *= 2. / len(msa_alignments[name_2])
            aln_coords_2[:, -1] *= 2. / len(msa_alignments[name_1])
            msa_alignments[name_1] = {name: "".join([sequence[i] if i != -1 else '-' for i in dtw_aln_1]) for name, sequence in
                                      msa_alignments[name_1].items()}
            msa_alignments[name_2] = {name: "".join([sequence[i] if i != -1 else '-' for i in dtw_aln_2]) for name, sequence in
                                      msa_alignments[name_2].items()}
            msa_alignments[name_int] = {**msa_alignments[name_1], **msa_alignments[name_2]}

            mean_coords = get_mean_coords_extra(aln_coords_1, aln_coords_2)
            mean_sec = get_mean_secondary(aln_sec_1, aln_sec_2, 0)
            self.final_structures.append(Structure(name_int, None, None, mean_sec, None, mean_coords))

        for x in range(0, self.tree.shape[0] - 1, 2):
            node_1, node_2, node_int = self.tree[x, 0], self.tree[x + 1, 0], self.tree[x, 1]
            assert self.tree[x + 1, 1] == node_int
            make_intermediate_node(node_1, node_2, node_int)

        node_1, node_2 = self.tree[-1, 0], self.tree[-1, 1]
        make_intermediate_node(node_1, node_2, "final")
        alignment = {**msa_alignments[self.final_structures[node_1].name], **msa_alignments[self.final_structures[node_2].name]}
        self.alignment = alignment
        return alignment

    def superpose(self, alignments: dict = None):
        """
        Superpose structures according to alignment
        """
        if alignments is None:
            alignments = self.alignment
        reference_index = 0
        reference_key = self.structures[reference_index].name
        core_indices = np.array([i for i in range(len(alignments[reference_key])) if '-' not in [alignments[n][i] for n in alignments]])
        aln_ref = helper.aligned_string_to_array(alignments[reference_key])
        ref_coords = self.structures[reference_index].coords[np.array([aln_ref[c] for c in core_indices])][:, :3]
        ref_centroid = rmsd_calculations.nb_mean_axis_0(ref_coords)
        ref_coords -= ref_centroid
        for i in range(len(self.structures)):
            if i == reference_index:
                self.structures[i].coords[:, :3] -= ref_centroid
            else:
                aln_c = helper.aligned_string_to_array(alignments[self.structures[i].name])
                common_coords_2 = self.structures[i].coords[np.array([aln_c[c] for c in core_indices])][:, :3]
                rotation_matrix, translation_matrix = rmsd_calculations.svd_superimpose(ref_coords, common_coords_2)
                self.structures[i].coords[:, :3] = rmsd_calculations.apply_rotran(self.structures[i].coords[:, :3], rotation_matrix,
                                                                                  translation_matrix)

    def make_pairwise_rmsd_coverage_matrix(self, alignments: dict = None, superpose_first: bool = True):
        """
        Find RMSDs and coverages of the alignment of each pair of sequences

        Parameters
        ----------
        alignments
        superpose_first
            if True then superposes all structures to first structure first

        Returns
        -------
        RMSD matrix, coverage matrix
        """
        if alignments is None:
            alignments = self.alignment
        num = len(self.structures)
        pairwise_rmsd_matrix = np.zeros((num, num))
        pairwise_rmsd_matrix[:] = np.nan
        pairwise_coverage = np.zeros((num, num))
        pairwise_coverage[:] = np.nan
        if superpose_first:
            self.superpose(alignments)
        for i in range(num):
            for j in range(i + 1, num):
                name_1, name_2 = self.structures[i].name, self.structures[j].name
                if isinstance(alignments[name_1], str):
                    aln_1 = helper.aligned_string_to_array(alignments[name_1])
                    aln_2 = helper.aligned_string_to_array(alignments[name_2])
                else:
                    aln_1 = alignments[name_1]
                    aln_2 = alignments[name_2]
                common_coords_1, common_coords_2 = get_common_coordinates(self.structures[i].coords[:, :3],
                                                                          self.structures[j].coords[:, :3], aln_1, aln_2)
                assert common_coords_1.shape[0] > 0
                if not superpose_first:
                    rot, tran = rmsd_calculations.svd_superimpose(common_coords_1, common_coords_2)
                    common_coords_2 = rmsd_calculations.apply_rotran(common_coords_2, rot, tran)
                pairwise_rmsd_matrix[i, j] = pairwise_rmsd_matrix[j, i] = rmsd_calculations.get_rmsd(common_coords_1, common_coords_2)
                pairwise_coverage[i, j] = pairwise_coverage[j, i] = common_coords_1.shape[0] / len(aln_1)
        return pairwise_rmsd_matrix, pairwise_coverage

    def get_aligned_features(self, alignments: dict = None):
        """
        Get dict of aligned features
        """
        if alignments is None:
            alignments = self.alignment
        feature_names = list(self.structures[0].features.keys())
        aligned_features = {}
        alignment_length = len(alignments[self.structures[0].name])
        for feature_name in feature_names:
            if feature_name == "secondary":
                continue
            aligned_features[feature_name] = np.zeros((len(self.structures), alignment_length))
            aligned_features[feature_name][:] = np.nan
            for p in range(len(self.structures)):
                farray = self.structures[p].features[feature_name]
                if "gnm" in feature_name or "anm" in feature_name:
                    farray = farray / np.nansum(farray ** 2) ** 0.5
                indices = [i for i in range(alignment_length) if alignments[self.structures[p].name][i] != '-']
                aligned_features[feature_name][p, indices] = farray
        return aligned_features

    def write_alignment(self, filename, alignments: dict = None):
        """
        Writes alignment to a fasta file
        """
        if alignments is None:
            alignments = self.alignment
        with open(filename, "w") as f:
            for key in alignments:
                f.write(f">{key}\n{alignments[key]}\n")

    def write_superposed_pdbs(self, output_pdb_folder, alignments: dict = None):
        """
        Superposes PDBs according to alignment and writes transformed PDBs to files
        (View with Pymol)

        Parameters
        ----------
        alignments
        output_pdb_folder
        """
        if alignments is None:
            alignments = self.alignment
        output_pdb_folder = Path(output_pdb_folder)
        if not output_pdb_folder.exists():
            output_pdb_folder.mkdir()
        reference_name = self.structures[0].name
        reference_pdb = pd.parsePDB(self.structures[0].pdb_file)
        core_indices = np.array([i for i in range(len(alignments[reference_name])) if '-' not in [alignments[n][i] for n in alignments]])
        aln_ref = helper.aligned_string_to_array(alignments[reference_name])
        ref_coords_core = reference_pdb[helper.get_alpha_indices(reference_pdb)].getCoords().astype(np.float64)[
            np.array([aln_ref[c] for c in core_indices])]
        ref_centroid = rmsd_calculations.nb_mean_axis_0(ref_coords_core)
        transformation = pd.Transformation(np.eye(3), -ref_centroid)
        reference_pdb = pd.applyTransformation(transformation, reference_pdb)
        pd.writePDB(str(output_pdb_folder / f"{reference_name}.pdb"), reference_pdb)
        for i in range(1, len(self.structures)):
            name = self.structures[i].name
            pdb = pd.parsePDB(self.structures[i].pdb_file)
            aln_name = helper.aligned_string_to_array(alignments[name])
            common_coords_2 = pdb[helper.get_alpha_indices(pdb)].getCoords().astype(np.float64)[np.array([aln_name[c] for c in core_indices])]
            rotation_matrix, translation_matrix = rmsd_calculations.svd_superimpose(ref_coords_core, common_coords_2)
            transformation = pd.Transformation(rotation_matrix, translation_matrix)
            pdb = pd.applyTransformation(transformation, pdb)
            pd.writePDB(str(output_pdb_folder / f"{name}.pdb"), pdb)
