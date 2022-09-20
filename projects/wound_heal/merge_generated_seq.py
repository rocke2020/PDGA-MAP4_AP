from itertools import count
import re, random
from pathlib import Path
import pickle
import json
import os, sys, shutil
import os, sys
sys.path.append(os.path.abspath('.'))
from util.peptide_util import basic_aa_3chars_to_1chars

result_dir = Path('results/anti_inflammation_v0.1')
seed = 1
out_dir = Path('/mnt/sda/bio_drug_corpus/AIpep/anti_inflammation/ga_outputs')
out_file = f'/mnt/sda/bio_drug_corpus/AIpep/anti_inflammation/ga_outputs/merged_result_seed{seed}.json'


def merge_and_copy():
    count = 0
    sequences = []
    for _dir in result_dir.iterdir():
        for result_file in _dir.iterdir():
            if result_file.name.endswith('results'):
                count += 1
                with open(result_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        items = line.split()
                        chars = items[1].split('-')
                        _chars = [basic_aa_3chars_to_1chars[char] for char in chars]
                        sequences.append(''.join(_chars))
    sequences = list(set(sequences))
    print(count)
    print(len(sequences))
    with open(out_file, 'w', encoding='utf-8') as f:
        json.dump(sequences, f, ensure_ascii=False, indent=4)
                
    shutil.copyfile('results/anti_inflammation_v0.1/anti_inflammation_0_seed1/param.txt', out_dir / f'param_seed{seed}.txt')
            


if __name__ == "__main__":
    merge_and_copy()
    pass