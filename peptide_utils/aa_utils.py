amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "L", "M", "N", "P", "K", "Q", "R", "S", "T", "V", "W", "Y"]

basic_aa_1chars_to_3chars_lower = {
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

basic_aa_3chars_lower_to_1chars = {
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

basic_aa_1chars_to_3chars = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "K": "LYS",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR"
}

basic_aa_3chars_to_1chars = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "LYS": "K",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y"
}


def convert_to_3_chars(sequences, out_file, all_upper_3chars=False):
    """  """
    new_seqs = []
    if all_upper_3chars:
        convert_dict = basic_aa_1chars_to_3chars
    else:
        convert_dict = basic_aa_1chars_to_3chars_lower
    for seq in sequences:
        chars = []
        for char in seq:
            chars.append(convert_dict[char])
        new_seqs.append('-'.join(chars))
    with open(out_file, 'w', encoding='utf-8') as f:
        for seq in new_seqs:
            f.write(seq+'\n')
    return new_seqs
