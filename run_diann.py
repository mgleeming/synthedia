import os, sys, argparse

DIA_NN_PATH = '/usr/diann/1.8/diann-1.8'
out_dir = os.path.join(input_dir,'out')

try:
    os.makedirs(os.path.join(input_dir,'out'))
except:
    pass

files = os.listdir(input_dir)

mzml_file = [_ for _ in files if '.mzml' in _.lower()]
fasta_file = [_ for _ in files if '.fasta' in _.lower()]

assert len(mzml_file) == 1
assert len(fasta_file) == 1

mzml_file = os.path.join(input_dir, mzml_file[0])
fasta_file = os.path.join(input_dir, fasta_file[0])

cmd = '%s ' %DIA_NN_PATH
cmd += '--f %s ' % mzml_file
cmd += '--fasta %s ' % fasta_file
cmd += '--out %s ' % os.path.join(out_dir, 'out.tsv')
cmd += '--out-lib %s ' % os.path.join(out_dir, 'out.lib.tsv')
cmd += '--lib --predictor  --threads 8 --verbose 3 --qvalue 0.01 --matrices --gen-spec-lib --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --smart-profiling --pg-level 1 --prosit --predictor'




diann.exe --f "D:\MaxQuant Searches\ML\DIA_NN\all\output_120.mzML" --lib "D:\MaxQuant Searches\ML\DIA_NN\all\out_monday-lib.predicted.speclib" --threads 20 --verbose 3 --out "D:\MaxQuant Searches\ML\DIA_NN\all\out_monday.tsv" --qvalue 0.01 --matrices  --out-lib "D:\MaxQuant Searches\ML\DIA_NN\all\out_monday-lib.tsv" --gen-spec-lib --prosit --var-mods 1 --window 6 --mass-acc 22.5447 --mass-acc-ms1 22.5447 --use-quant --double-search --individual-mass-acc --individual-windows --reanalyse --pg-level 0
os.system(cmd)

