import re, random
from pathlib import Path
import pickle
import json
import os, sys, shutil
sys.path.append(os.path.abspath('.'))
from util.peptide_util import basic_aa_1chars_to_3chars


file = '/home/qcdong/codes/MultiPep/APD3_data/woundheal.txt'
sequences = set()
with open(file, 'r', encoding='utf-8') as f:
    for i, line in enumerate(f):
        if i == 0: continue
        if not line.strip():
            break
        if line.startswith('AP'): continue
        items = line.split('\t')
        sequences.add(items[1].strip())
print(len(sequences))
with open('projects/wound_heal/orig_data/orig_seq_1char.txt', 'w', encoding='utf-8') as f:
    for seq in sequences:
        f.write(f'{seq}\n')


def convert_to_3_chars(sequences):
    new_seqs = []
    for seq in sequences:
        chars = []
        for char in seq:
            chars.append(basic_aa_1chars_to_3chars[char])
        new_seqs.append('-'.join(chars))
    with open('projects/wound_heal/orig_data/orig_seq_3char.json', 'w', encoding='utf-8') as f:
        json.dump(new_seqs, f, ensure_ascii=False, indent=4)


convert_to_3_chars(sequences)