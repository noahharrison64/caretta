    
from prody import LOGGER
from functools import partialmethod
import tqdm
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
        
        tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)
        from caretta import multiple_alignment
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
        
        tqdm.__init__ = partialmethod(tqdm.__init__, disable=False)
        if return_paths:
            return [
                f for f in output[1].pdb_folder.glob('*.pdb') if f.stem in [in_f.stem for in_f in input_files]
            ]
        else:
            return [
                mda.Universe(f) for f in output[1].pdb_folder.glob('*.pdb') if f.stem in [in_f.stem for in_f in input_files]
            ]