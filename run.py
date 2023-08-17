from multiprocessing import Pool
import sys, os, shutil
sys.path.append(os.path.abspath('.'))
from utils.arg_util import ArgparseUtil, log_args
from utils.log_util import logger
from utils.file_util import FileUtil
from genetic_algorithm.PDGA_only20aa import PDGA
from genetic_algorithm.fingerprints import get_ap_distancefn, get_ap_fingerprintfn
from peptide_utils.aa_utils import convert_to_3_chars, basic_aa_3chars_lower_to_1chars
from pathlib import Path
from pandas import DataFrame


# Use the RDKit Atom pair as fingerprint:
fingerprintfn = get_ap_fingerprintfn()
distancefn = get_ap_distancefn()
postfix = '_3chars_seq'


def convert_to_3chars_seqs(input_file):
    """  """
    if input_file.is_file():
        sequences = FileUtil.read_raw_text(input_file)
        logger.info(f'input seqs num {len(sequences)}, sequences top 5: {sequences[:5]}')
        query_sequences = convert_to_3_chars(sequences, input_file.parent / f'{input_file.stem}{postfix}.txt')
    else:
        logger.info('%s is not existent as a file', input_file)
        raise Exception()
    return query_sequences


def run_one_peptide(peptied_num=0):
    query_seq = query_sequences[peptied_num]
    query_name = f'peptied_num_{peptied_num}'
    ga = PDGA(pop_size=200, mut_rate=0.5, gen_gap=0.5, query=query_seq, sim_treshold=0.6,
                porpouse="linear", folder=str(result_dir / f"{query_name}_seed{args.seed}"),
                fingerprintfn=fingerprintfn,
                distancefn=distancefn,
                query_name=query_name, peptied_num=peptied_num, similar_num = targt_num_per_seq,
                is_peptide_sequence=True,
                verbose=False, seed=args.seed)
    ga.write_param()
    logger.info(f'Starts peptied_num {peptied_num} ......')
    ga.run()


def merge():
    """
    """
    result_files_count = 0
    sequences = []
    for _dir in result_dir.iterdir():
        for result_file in _dir.iterdir():
            if result_file.name.endswith('results'):
                result_files_count += 1
                with open(result_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        items = line.split()
                        chars = items[1].split('-')
                        _chars = [basic_aa_3chars_lower_to_1chars[char] for char in chars]
                        sequences.append(''.join(_chars))
    sequences = list(set(sequences))
    logger.info('result_files_count %s, unique len(sequences) %s', result_files_count, len(sequences))
    df = DataFrame({'Sequence': sequences})
    df['len'] = df['Sequence'].map(len)
    logger.info('%s', df['len'].describe())
    df = df[(df['len'] >= args.min_len) & (df['len'] <= args.max_len)]
    sequences = df['Sequence']
    logger.info('Valid length seq num %s', len(sequences))
    out_file = result_dir / f'merged_result_seed{args.seed}.txt'
    FileUtil.write_raw_text(sequences, out_file)


if __name__ == "__main__":
    args = ArgparseUtil().genetic_algorithm_generator()
    log_args(args, logger)
    input_data_dir = Path(args.input_data_dir)
    input_file = input_data_dir / args.input_filename
    query_sequences = convert_to_3chars_seqs(input_file)
    targt_num_per_seq = args.total_target_num // len(query_sequences) + 1
    logger.info(f'targt_num_per_seq {targt_num_per_seq}')

    result_dir = Path(f"{args.out_dir}/{args.task_name}/{input_file.stem}")
    clean_all_files_in_task_result_dir = 1
    if clean_all_files_in_task_result_dir and result_dir.is_dir():
        shutil.rmtree(result_dir, ignore_errors=False)
    with Pool(24) as pool:
        for i in pool.imap_unordered(run_one_peptide, range(len(query_sequences))):
            pass
    merge()