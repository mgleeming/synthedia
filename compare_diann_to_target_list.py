import os, sys
import pandas as pd

def compare_diann_to_target_list(in_dir):

    peptide_table = os.path.join(in_dir, [_ for _ in os.listdir(in_dir) if 'peptide_table.tsv' in _][0])
    peptide_table = pd.read_csv(peptide_table, sep = '\t')
    print(peptide_table)

    diann_hits_file = os.path.join(in_dir, 'out', 'out.tsv')
    diann_table = pd.read_csv(diann_hits_file, sep = '\t')
    print(diann_table)

    return

compare_diann_to_target_list('/home/mleeming/Code/synthedia/data/run_length_20')
