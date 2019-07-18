from warp_aligner import multiple_structure_alignment as msa
from warp_aligner import helper
from pathlib import Path
import prody as pd
import numpy as np
import typing


def get_pdbs(pdb_dir) -> dict:
    """
    Get all the proteins in a directory as pd AtomGroup objects

    Parameters
    ----------
    pdb_dir

    Returns
    -------
    dictionary {name: pd.AtomGroup}
    """
    pdb_dir = Path(pdb_dir)
    pdbs = {}
    for pdb_file in pdb_dir.glob("[!.]*"):
        pdbs[pdb_file.stem] = pd.parsePDB(pdb_file)
    return pdbs


def get_sequences(pdbs: typing.Dict[str, pd.AtomGroup]) -> typing.Dict[str, str]:
    """
    Get pdb sequences
    Parameters
    ----------
    pdbs
        dict of {name: pd.AtomGroup}

    Returns
    -------
    dict of {name: sequence}
    """
    return {n: pdbs[n][helper.get_alpha_indices(pdbs[n])].getSequence() for n in pdbs}


def get_sequence_alignment(sequences: typing.Dict[str, str], directory, name="ref") -> typing.Dict[str, str]:
    """
    Use clustal-omega to make a multiple sequence alignment

    Parameters
    ----------
    sequences
        dict of {name: sequence}
    directory
        dir to save sequence and alignment file
    name
        name to give sequence/alignment file
    Returns
    -------
    dict of {name: aln_sequence}
    """
    directory = Path(directory)
    sequence_file = directory / f"{name}.fasta"
    with open(sequence_file, "w") as f:
        for n in sequences:
            f.write(f">{n}\n{sequences[n]}\n")
    aln_sequence_file = directory / f"{name}_aln.fasta"
    helper.clustal_msa_from_sequences(sequence_file, aln_sequence_file)
    return helper.get_sequences_from_fasta(aln_sequence_file)


def get_msa_class(pdb_dir) -> msa.StructureMultiple:
    """
    Get msa_class (msa.StructureMultiple) from the PDB files in a directory.
    """
    pdbs = get_pdbs(pdb_dir)
    names = list(pdbs.keys())
    sequences = get_sequences(pdbs)
    coordinates = [pdbs[n][helper.get_beta_indices(pdbs[n])].getCoords().astype(np.float64) for n in names]
    structures = msa.make_structures(names, [sequences[n] for n in names], coordinates)
    structures_multiple = msa.StructureMultiple(structures)
    return structures_multiple


def get_msa(pdb_dir, sequence_dir, name="ref", gap_open_penalty: float = 0., gap_extend_penalty: float = 0.):
    """
    Example usage of above functions to make a structure-guided multiple sequence alignment from a directory of PDB files
    """
    msa_class = get_msa_class(pdb_dir)
    aln_sequences = get_sequence_alignment({s.name: s.sequence for s in msa_class.structures}, sequence_dir, name)
    return msa_class.align(aln_sequences, gap_open_penalty, gap_extend_penalty)
