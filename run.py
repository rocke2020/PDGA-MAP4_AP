from pdga.PDGA_only20aa import PDGA
from pdga import get_map4_distancefn, get_map4_fingerprintfn, get_ap_distancefn, get_ap_fingerprintfn
import sys, json
import numpy as np
from util.log_util import logger
from multiprocessing import Pool

seed = 0
query_sequences_file = 'projects/wound_heal/orig_data/orig_seq_3char.json'
with open(query_sequences_file, 'r', encoding='utf-8') as f:
    query_sequences = json.load(f)


total_similar_num = 1000_000
current_similar_num = total_similar_num // len(query_sequences) + 1
logger.info(f'current_similar_num {current_similar_num}')
# query_smiles = "C[C@H]1CCC[C@@H]2[C@H](CC(/C(C)=C/C3=CSC(C)=N3)OC(C[C@H](O)C(C)(C)C([C@H](C)[C@H]1O)=O)=O)O2"
# query_name = "EpothilonA"

# To use MAP4 as fingerprint:
fingerprintfn = get_map4_fingerprintfn(return_strings=True, radius=1)
distancefn = get_map4_distancefn(use_string=True, tanimoto=True)

# To use the RDKit Atom pair as fingerprint:
#     fingerprintfn = get_ap_fingerprintfn()
#     distancefn = get_ap_distancefn()

def run_one_peptide(peptied_num=0):
    query = query_sequences[peptied_num]
    query_name = f'anti_inflammation_{peptied_num}'        
    # if you want to avoid the possibility of metylation of the amide bond set methyl to False
    ga = PDGA(pop_size=200, mut_rate=0.5, gen_gap=0.5, query=query, sim_treshold=0.6, 
                porpouse="linear", folder=f"results/anti_inflammation/{query_name}_seed{seed}", fingerprintfn=fingerprintfn, 
                distancefn=distancefn, 
                query_name=query_name, peptied_num=peptied_num, similar_num = current_similar_num, 
                is_peptide_sequence=True, methyl=True, 
                verbose=True, seed=seed)

    ga.write_param()
    logger.info(f'Starts peptied_num {peptied_num} ......')
    ga.run()


if __name__ == "__main__":
    # run_one_peptide()
    with Pool(23) as pool:
        for i in pool.imap_unordered(run_one_peptide, range(len(query_sequences))):
            pass
