import os, sys, argparse
DIA_NN_PATH = '/usr/diann/1.8/diann-1.8'

parser = argparse.ArgumentParser(
    description = 'Run DIA-NN'
)

parser.add_argument( '--dir', required = True, type = str)

options =  parser.parse_args()

out_dir = os.path.join(options.dir,'out')

try:
    os.makedirs(os.path.join(options.dir,'out'))
except:
    pass

files = os.listdir(options.dir)

mzml_file = [_ for _ in files if '.mzml' in _.lower()]
fasta_file = [_ for _ in files if '.fasta' in _.lower()]

assert len(mzml_file) == 1
assert len(fasta_file) == 1

mzml_file = os.path.join(options.dir, mzml_file[0])
fasta_file = os.path.join(options.dir, fasta_file[0])


cmd = '%s ' %DIA_NN_PATH
cmd += '--f %s ' % mzml_file
cmd += '--fasta %s ' % fasta_file
cmd += '--out %s ' % os.path.join(out_dir, 'out.tsv')
cmd += '--out-lib %s ' % os.path.join(out_dir, 'out.lib.tsv')
cmd += '--lib --predictor  --threads 8 --verbose 3 --qvalue 0.01 --matrices --gen-spec-lib --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --smart-profiling --pg-level 1 --prosit --predictor'

print(cmd)

os.system(cmd)
