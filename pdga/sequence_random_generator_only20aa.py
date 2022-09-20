import random 
from util.peptide_util import basic_aminoacids

class SequenceGenerator:
    def __init__(self, verbose=False):
            # variables for random generation of dendrimers
            self.AA4rndm = basic_aminoacids
            self.min_aa_no = 10
            self.max_aa_no = 20

            self.verbose = verbose
    
    def generate(self):
        """Generates random implicit sequences of max "max_gen_no" generation dendrimers
            with max "max_aa_no" AA in each generation, picking from AA4random, B4random
            (probability of position to be empty intrinsic in these lists). 
        Returns:
            string -- implicit sequence of a random dendrimer
        """
        new_random_seq = []
        aa_count = 0
        max_num = random.randint(self.min_aa_no, self.max_aa_no)
        while aa_count < max_num:
            new_random_seq.append(random.choice(self.AA4rndm))
            aa_count += 1
        new_random_seq = ''.join(new_random_seq)
        return new_random_seq
