from pdga.PDGA_only20aa import PDGA
from pdga import get_ap_distancefn, get_ap_fingerprintfn
import sys, json
import numpy as np
from utils.log_util import logger
from multiprocessing import Pool
from utils.file_util import FileUtil


seed = 0
task_name = 'cpp_peptides'
postfix = '_3chars_seq'
orig_filename = 'FGF'
filename = f'{orig_filename}{postfix}'
query_sequences_file = f'projects/{task_name}/orig_data/{filename}.txt'
query_sequences = FileUtil.read_raw_text(query_sequences_file)

total_similar_num = 10_000_000
similar_num_per_seq = total_similar_num // len(query_sequences) + 1
logger.info(f'similar_num_per_seq {similar_num_per_seq}')
# query_smiles = "C[C@H]1CCC[C@@H]2[C@H](CC(/C(C)=C/C3=CSC(C)=N3)OC(C[C@H](O)C(C)(C)C([C@H](C)[C@H]1O)=O)=O)O2"
# query_name = "EpothilonA"

use_map4 = False
if use_map4:
    from pdga.fingerprints_map4 import get_map4_distancefn, get_map4_fingerprintfn
    # To use MAP4 as fingerprint:
    fingerprintfn = get_map4_fingerprintfn(return_strings=True, radius=1)
    distancefn = get_map4_distancefn(use_string=True, tanimoto=True)
else:
    # To use the RDKit Atom pair as fingerprint:
    fingerprintfn = get_ap_fingerprintfn()
    distancefn = get_ap_distancefn()


def run_one_peptide(peptied_num=0):
    query_seq = query_sequences[peptied_num]
    query_name = f'{task_name}_{peptied_num}'
    ga = PDGA(pop_size=200, mut_rate=0.5, gen_gap=0.5, query=query_seq, sim_treshold=0.6,
                porpouse="linear", folder=f"results/{task_name}/{query_name}_seed{seed}",
                fingerprintfn=fingerprintfn,
                distancefn=distancefn,
                query_name=query_name, peptied_num=peptied_num, similar_num = similar_num_per_seq,
                is_peptide_sequence=True,
                verbose=True, seed=seed)

    ga.write_param()
    logger.info(f'Starts peptied_num {peptied_num} ......')
    ga.run()


if __name__ == "__main__":
    # run_one_peptide()
    with Pool(24) as pool:
        for i in pool.imap_unordered(run_one_peptide, range(len(query_sequences))):
            pass
