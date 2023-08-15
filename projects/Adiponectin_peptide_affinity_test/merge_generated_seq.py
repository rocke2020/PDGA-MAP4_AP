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
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
ic.lineWrapWidth = 120

# Notice to modify the task name which is the input!!
task_name = 'Adiponectin_peptide_affinity_test'  # anti_inflammation wound_heal
result_dir = Path(f'results/{task_name}')

# Output
root_dir = Path(f'/mnt/sda/bio_drug_corpus/peptides/Adiponectin peptide affinity test/ga_outputs')
orig_filename = 'GST'
out_dir = root_dir / orig_filename
out_file = out_dir / f'merged_genetic_algorithm_result.txt'
out_dir.mkdir(exist_ok=1, parents=1)
min_len = 5
max_len = 10


def merge_and_copy():
    """
    seed 1, Valid length seq num 2,589,320 from 5M
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
    ic(result_files_count)
    ic(len(sequences))
    df = DataFrame({'Sequence': sequences})
    df['len'] = df['Sequence'].map(len)
    logger.info('\ndf[len].describe()\n%s', df['len'].describe())
    df = df[(df['len'] >= min_len) & (df['len'] <= max_len)]
    sequences = df['Sequence']
    logger.info('Valid length seq num %s', len(sequences))
    FileUtil.write_raw_text(sequences, out_file)
    for file in result_dir.iterdir():
        if file.is_file(): continue
        for _file in file.iterdir():
            if _file.name == 'param.txt':
                shutil.copyfile(_file, out_dir / f'param.txt')
                return


if __name__ == "__main__":
    merge_and_copy()