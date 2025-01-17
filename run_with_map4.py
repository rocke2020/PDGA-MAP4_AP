from genetic_algorithm.PDGA_only20aa import PDGA
from genetic_algorithm.fingerprints import get_ap_distancefn, get_ap_fingerprintfn
import sys, os, shutil
from utils.log_util import logger
from multiprocessing import Pool
from utils.file_util import FileUtil
sys.path.append(os.path.abspath('.'))


use_map4 = False
if use_map4:
    from genetic_algorithm.fingerprints_map4 import get_map4_distancefn, get_map4_fingerprintfn
    # To use MAP4 as fingerprint:
    fingerprintfn = get_map4_fingerprintfn(return_strings=True, radius=1)
    distancefn = get_map4_distancefn(use_string=True, tanimoto=True)
else:
    # To use the RDKit Atom pair as fingerprint:
    fingerprintfn = get_ap_fingerprintfn()
    distancefn = get_ap_distancefn()


## draft codes
seed = 1
task_name = 'cpp_peptides'
postfix = '_3chars_seq'
orig_filename = 'FGF'
filename = f'{orig_filename}{postfix}'
query_sequences_file = f'projects/{task_name}/input_data/{filename}.txt'
query_sequences = FileUtil.read_raw_text(query_sequences_file)

single_char_seqs = FileUtil.read_raw_text(f'projects/{task_name}/input_data/{orig_filename}.txt')
# To track the codes to check mistakes on input files.
logger.info('query_sequences_file %s', query_sequences_file)
logger.info(f'seed {seed}, seqs num {len(single_char_seqs)}, single_char_seqs[:5] {single_char_seqs[:5]}')

total_target_num = 6_000_000
targt_num_per_seq = total_target_num // len(query_sequences) + 1
logger.info(f'targt_num_per_seq {targt_num_per_seq}')


def run_one_peptide(peptied_num=0):
    query_seq = query_sequences[peptied_num]
    query_name = f'{task_name}_{peptied_num}'
    ga = PDGA(pop_size=200, mut_rate=0.5, gen_gap=0.5, query=query_seq, sim_treshold=0.6,
                porpouse="linear", folder=f"results/{task_name}/{query_name}_seed{seed}",
                fingerprintfn=fingerprintfn,
                distancefn=distancefn,
                query_name=query_name, peptied_num=peptied_num, similar_num = targt_num_per_seq,
                is_peptide_sequence=True,
                verbose=False, seed=seed)
    ga.write_param()
    logger.info(f'Starts peptied_num {peptied_num} ......')
    ga.run()


if __name__ == "__main__":
    # run_one_peptide()
    clean_all_files_in_task_result_dir = 1
    if clean_all_files_in_task_result_dir:
        shutil.rmtree(f"./results/{task_name}", ignore_errors=False)
    with Pool(24) as pool:
        for i in pool.imap_unordered(run_one_peptide, range(len(query_sequences))):
            pass
