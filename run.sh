nohup python run.py \
--task_name cpp_peptides \
--input_data_dir /home/qcdong/codes/PDGA-MAP4_AP/projects/cpp_peptides/input_data \
--input_filename FGF.txt \
--out_dir /home/qcdong/codes/PDGA-MAP4_AP/results \
--total_target_num 10000 \
> genetic_algorithm_generator.log 2>&1 &