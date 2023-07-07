    
from prody import LOGGER
from caretta import multiple_alignment
from functools import partialmethod
import numba as nb
import MDAnalysis as mda
from typing import Union, List
from pathlib import Path

def run_caretta_alignment(
        input_files: List[Path],
        output_folder: Path,
        nthreads, 
        return_paths: bool,
        verbose: bool = False,
    ) -> Union[List[mda.Universe], List[Path]]:
        
        LOGGER.verbosity = 'warning'
        nb.set_num_threads(nthreads)
        multiple_alignment.trigger_numba_compilation()       

        output = multiple_alignment.align_from_structure_files(
            input_files=input_files,
            output_folder=output_folder,
            write_pdb=True,
            num_threads=nthreads,
            align_cofactors=True,
            verbose=verbose
                )
        out_files =  [list(output[1].pdb_folder.glob(f'*{in_file.stem}*'))[0] for in_file in input_files]
        if return_paths:
           return out_files
        return [mda.Universe(out_file) for out_file in out_files]