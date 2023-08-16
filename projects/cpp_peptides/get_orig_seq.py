from pathlib import Path
import os, sys, shutil
sys.path.append(os.path.abspath('.'))
from utils.peptide_util import convert_to_3_chars
from utils.file_util import FileUtil


orig_data_dir = Path('projects/cpp_peptides/orig_data')
postfix = '_3chars_seq'

def create_3chars_seq():
    """  """
    for file in orig_data_dir.glob('*.txt'):
        if file.stem.endswith(postfix): continue
        sequences = FileUtil.read_raw_text(file)
        convert_to_3_chars(sequences, orig_data_dir / f'{file.stem}{postfix}.txt')


create_3chars_seq()