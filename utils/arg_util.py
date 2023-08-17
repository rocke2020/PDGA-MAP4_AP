import argparse, os
from datetime import datetime


DATE_TIME = "%Y_%m_%d %H:%M:%S"


class ArgparseUtil(object):
    """
    参数解析工具类
    """
    def __init__(self):
        """ Basic args """
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("--seed", default=2, type=int)

    def af_batch_runner(self):
        """ task args """
        self.parser.add_argument('--gpu_device_id', default=1, type=int, help='the GPU NO.')
        self.parser.add_argument("--input_root_dir", type=str, default='/mnt/sdc/af_input/tasks', help="")
        self.parser.add_argument("--out_root_dir", type=str, default='/mnt/sdc/af_out/tasks', help="")
        self.parser.add_argument("--overwrite_predicted_pdb", type=int, default=1, help="0 is false, 1 is true")
        self.parser.add_argument("--model_preset", type=str, default='multimer', help="")
        self.parser.add_argument("--task_name", type=str, default='', help="")
        self.parser.add_argument("--sub_task_name", type=str, default='', help="")
        self.parser.add_argument("--run_relax", type=int, default=1, help="0 is false, 1 is true")
        self.parser.add_argument("--run_peptide_msa", type=int, default=0, help="0 is false, 1 is true")
        self.parser.add_argument("--reverse", type=int, default=0, help="0 is false, 1 is true")
        self.parser.add_argument("--part_num", type=int, default=1, help="NB: starts from 1")
        self.parser.add_argument("--total_parts", type=int, default=2, help="total parts num of input files")
        self.parser.add_argument('--collect_results', default=0, type=int,
                                 help='0 false, 1 true, not run alpha fold, but just collect results.')
        self.parser.add_argument('--merge_results', default=0, type=int,
                                 help='0 false, 1 true, not run alpha fold, but merge_results from different IPs.')
        args = self.parser.parse_args()
        return args

    def pepetide_generator(self):
        """  """
        self.parser.add_argument('--gpu_device_id', default=1, type=int, help='the GPU NO.')
        self.parser.add_argument("--data_type", type=str, default='acne', help="")
        self.parser.add_argument("--target_sample_num", type=int, default=10000, help="")
        self.parser.add_argument("--sample_batch_size", type=int, default=10000, help="")
        self.parser.add_argument("--min_seq_len", type=int, default=6, help="")
        self.parser.add_argument("--max_seq_len", type=int, default=15, help="")
        self.parser.add_argument("--out_dir", type=str, default='output/generated', help="")
        args = self.parser.parse_args()
        return args

    def genetic_algorithm_generator(self):
        """  """
        self.parser.add_argument("--task_name", type=str, default='', help="")
        self.parser.add_argument("--input_data_dir", type=str, default='', help="")
        self.parser.add_argument("--input_filename", type=str, default='', 
                                 help="Only support .txt file, one seq per line")
        self.parser.add_argument("--out_dir", type=str, default='', help="")
        self.parser.add_argument("--total_target_num", type=int, default=10_000, help="")
        self.parser.add_argument("--min_len", type=int, default=5, help="")
        self.parser.add_argument("--max_len", type=int, default=10, help="")
        args = self.parser.parse_args()
        return args


def save_args(args, output_dir='.', with_time_at_filename=False):
    os.makedirs(output_dir, exist_ok=True)

    t0 = datetime.now().strftime(DATE_TIME)
    if with_time_at_filename:
        out_file = os.path.join(output_dir, f"args-{t0}.txt")
    else:
        out_file = os.path.join(output_dir, f"args.txt")
    with open(out_file, "w", encoding='utf-8') as f:
        f.write(f'{t0}\n')
        for arg, value in vars(args).items():
            f.write(f"{arg}: {value}\n")


def log_args(args, logger):
    for arg, value in vars(args).items():
        logger.info(f"{arg}: {value}")
