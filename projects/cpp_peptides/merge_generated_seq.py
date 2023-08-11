from itertools import count
import re, random
from pathlib import Path
import pickle
import json
import os, sys, shutil
from pandas import DataFrame
sys.path.append(os.path.abspath('.'))
from utils.peptide_util import basic_aa_3chars_to_1chars
from utils.file_util import FileUtil
from utils.log_util import logger


task_name = 'cpp_peptides'  # anti_inflammation wound_heal
result_dir = Path(f'results/{task_name}')
seed = 0
out_dir = Path(f'/mnt/sda/bio_drug_corpus/AIpep/cell_penetrating/ga_outputs')
out_file = out_dir / f'merged_result_seed{seed}.txt'
out_dir.mkdir(exist_ok=1, parents=1)
min_len = 5
max_len = 10

def merge_and_copy():
    """
    443 valid length among orig 1000 target where 693 unique ones.
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
                        _chars = [basic_aa_3chars_to_1chars[char] for char in chars]
                        sequences.append(''.join(_chars))
    sequences = list(set(sequences))
    print(result_files_count)
    print(len(sequences))
    df = DataFrame({'Sequence': sequences})
    df['len'] = df['Sequence'].map(len)
    logger.info('%s', df['len'].describe())
    df = df[(df['len'] >= min_len) & (df['len'] <= max_len)]
    sequences = df['Sequence']
    logger.info('Valid length seq num %s', len(sequences))
    FileUtil.write_raw_text(sequences, out_file)
    shutil.copyfile(result_dir / f'{task_name}_0_seed{seed}/param.txt', out_dir / f'param_seed{seed}.txt')


if __name__ == "__main__":
    merge_and_copy()