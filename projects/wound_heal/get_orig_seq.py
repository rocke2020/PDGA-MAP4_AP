import re, random
from pathlib import Path
import pickle
import json
import os, sys, shutil
sys.path.append(os.path.abspath('.'))
from utils.peptide_util import convert_to_3_chars


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


convert_to_3_chars(sequences, 'projects/wound_heal/orig_data/orig_seq_3char.json')