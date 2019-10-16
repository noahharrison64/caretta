from pathlib import Path

import fire

from caretta import msa_numba


def align_cli(input_pdb,
              dssp_dir="caretta_tmp", num_threads=20, extract_all_features=True,
              gap_open_penalty=1., gap_extend_penalty=0.01, consensus_weight=1.,
              write_fasta=True, output_fasta_filename=Path("./result.fasta"),
              write_pdb=True, output_pdb_folder=Path("./result_pdb/"),
              write_features=True, output_feature_filename=Path("./result_features.pkl"),
              write_class=True, output_class_filename=Path("./result_class.pkl"),
              overwrite_dssp=False):
    """
    Caretta aligns protein structures and returns a sequence alignment, a set of aligned feature matrices, superposed PDB files, and
    a class with intermediate structures made during progressive alignment.
    Parameters
    ----------
    input_pdb
        Can be \n
        A folder with input protein files \n
        A file which lists PDB filenames on each line \n
        A file which lists PDB IDs on each line \n
    dssp_dir
        Folder to store temp DSSP files (default caretta_tmp)
    num_threads
        Number of threads to use for feature extraction
    extract_all_features
        True => obtains all features (default True) \n
        False => only DSSP features (faster)
    gap_open_penalty
        default 1
    gap_extend_penalty
        default 0.01
    consensus_weight
        default 1
    write_fasta
        True => writes alignment as fasta file (default True)
    output_fasta_filename
        Fasta file of alignment (default result.fasta)
    write_pdb
        True => writes all protein PDB files superposed by alignment (default True)
    output_pdb_folder
        Folder to write superposed PDB files (default result_pdb)
    write_features
        True => writes aligned features a s a dictionary of numpy arrays into a pickle file (default True)
    output_feature_filename
        Pickle file to write aligned features (default result_features.pkl)
    write_class
        True => writes StructureMultiple class with intermediate structures and tree to pickle file (default True)
    output_class_filename
        Pickle file to write StructureMultiple class (default result_class.pkl)
    overwrite_dssp
        Forces DSSP to rerun (default False)
    """
    msa_numba.StructureMultiple.align_from_pdb_files(input_pdb,
                                                     dssp_dir, num_threads, extract_all_features,
                                                     gap_open_penalty, gap_extend_penalty, consensus_weight,
                                                     write_fasta, output_fasta_filename,
                                                     write_pdb, output_pdb_folder,
                                                     write_features, output_feature_filename,
                                                     write_class, output_class_filename,
                                                     overwrite_dssp)


if __name__ == '__main__':
    fire.Fire(align_cli)