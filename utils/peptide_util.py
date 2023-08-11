from collections import defaultdict
from pandas import Series, DataFrame
import random
import pickle, json
import os, sys
sys.path.append(os.path.abspath('.'))


# 20 natrual
basic_aminoacids = ["A","C","D","E","F","G","H","I","L","M","N","P","K","Q","R","S","T","V","W","Y"]
multi_chars_to_single_char_dict = {
    'Arg': 'R', 'His': 'H', 'Lys': 'K', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Thr': 'T', 'Asn': 'N',
    'Gln': 'Q', 'Cys': 'C', 'Sec': 'U', 'Gly': 'G', 'Pro': 'P', 'Ala': 'A', 'Ile': 'I', 'Leu': 'L',
    'Met': 'M', 'Phe': 'F', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Dap': '1', 'Dab': '2',
    'BOrn': '3', 'BLys': '4', 'Hyp': 'Z', 'Orn': 'O', 'bAla': '!', 'Gaba': '?', 'dDap': '5',
    'dDab': '6',
    'dBOrn': '7', 'dBLys': '8', 'dArg': 'r', 'dHis': 'h', 'dLys': 'k', 'dAsp': 'd', 'dGlu': 'e',
    'dSer': 's',
    'dThr': 't', 'dAsn': 'n', 'dGln': 'q', 'dCys': 'c', 'dSec': 'u', 'dGly': 'g', 'dPro': 'p',
    'dAla': 'a',
    'dIle': 'i', 'dLeu': 'l', 'dMet': 'm', 'dPhe': 'f', 'dTrp': 'w', 'dTyr': 'y', 'dVal': 'v',
    'dHyp': 'z', 'dOrn': 'o', 'a5a': '=', 'a6a': '%', 'a7a': '$', 'a8a': '@', 'a9a': '#',
    'Cys1': 'Ä', 'Cys2': 'Ö', 'Cys3': 'Ü', 'dCys1': 'ä', 'dCys2': 'ö', 'dCys3': 'ü',
    'Ac': '&', 'NH2': '+', 'met': '-', 'cy': 'X'}
basic_aa_1chars_to_3chars = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "K": "Lys",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr"
}
basic_aa_3chars_to_1chars = {
    "Ala": "A",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Pro": "P",
    "Lys": "K",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Val": "V",
    "Trp": "W",
    "Tyr": "Y"
}


## Non-containing methionine sequences were preferred. Methionine，简写M，Met, 甲硫氨酸, for diversity, not used in codes.
# often_AMP_aminoacids = basic_aminoacids.copy()
# often_AMP_aminoacids.remove('M')
C_TERMINUS = "cTerminus"
N_TERMINUS = "nTerminus"
LEAST_SEQ_LENGTH = 4
NOT_HEMOLYTIC_KEY = 'isNotHemolytic'
ACTIVITY_KEY = 'activity'


def randomChoice(l):
    return l[random.randint(0, len(l) - 1)]


def _sample(model, n):
    sampled_seq = model.sample(n)
    sequences = []
    for s in sampled_seq:
        sequences.append(model.vocabulary.tensor_to_seq(s))
    return sequences


def novelty(seqs, list_):
    novel_seq = []
    for s in seqs:
        if s not in list_:
            novel_seq.append(s)
    return novel_seq, (len(novel_seq)/len(seqs))*100


def uniqueness(seqs):
    unique_seqs = defaultdict(int)
    for s in seqs:
        unique_seqs[s] += 1
    return unique_seqs, (len(unique_seqs)/len(seqs))*100


def is_natural(seq):
    try:
        seq = seq.upper()
        for aa in seq:
            if aa not in basic_aminoacids:
                return False
        return True
    except:
        return False


def float_ignore_plus_minus(mynumber):
    try:
        return sum(map(float,mynumber.split("±")))
    except:
        return float("inf")


def check_active(unit, concentration):
    if ((unit == "µM" and concentration < 10) or (unit == "nM" and concentration < 10000)
        or (unit == "µg/ml" and concentration < 32)):
        return True


def check_inactive(unit, concentration):
    if ((unit == "µM" and concentration > 10) or (unit == "nM" and concentration > 10000)
        or (unit == "µg/ml" and concentration > 32)):
        return True


def convert_2d_values_to_row_col(input_s:Series):
    """
    Args: each item in input series is 2d list [[2, 3], [1, 3]]
        [[[2, 3], [1, 3]], ...]
    Returns:
        dataFrame
            0  1
        0  2  3
        0  1  3
    """
    if len(input_s) == 0:
        return DataFrame()
    s = input_s.apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    out = s.apply(lambda x: Series(x))
    return out


def get_terminus_names(dbaasp_peptide):
    n_terminus = dbaasp_peptide.get(N_TERMINUS)
    if n_terminus:
        n_name = n_terminus['name']
    else:
        n_name = 'nan'
    c_terminus = dbaasp_peptide.get(C_TERMINUS)
    if c_terminus:
        c_name = c_terminus['name']
    else:
        c_name = 'nan'
    return n_name, c_name


def is_valid_terminus(n_name, c_name):
    if (n_name == 'nan' or n_name == 'ACT') and (c_name == 'nan' or c_name == 'AMD'):
        return True


def convert_to_3_chars(sequences, out_file):
    """  """
    new_seqs = []
    for seq in sequences:
        chars = []
        for char in seq:
            chars.append(basic_aa_1chars_to_3chars[char])
        new_seqs.append('-'.join(chars))
    with open(out_file, 'w', encoding='utf-8') as f:
        for seq in new_seqs:
            f.write(seq+'\n')


if __name__ == "__main__":
    # basic_aa_3chars_to_1chars = {}
    # for k, v in basic_aa_1chars_to_3chars.items():
    #         basic_aa_3chars_to_1chars[v] = k
    # with open('1.json', 'w', encoding='utf-8') as f:
    #     json.dump(basic_aa_3chars_to_1chars, f, ensure_ascii=False, indent=4)
    print(len(basic_aminoacids))
