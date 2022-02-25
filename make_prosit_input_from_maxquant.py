import os, argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description = 'Generate Prosit input from MaxQuant txt directory'
)

parser.add_argument( '--mq_txt_dir', required = False, type = str)
options =  parser.parse_args()

MQ_TXT_DIR = options.mq_txt_dir

evidence = pd.read_csv(os.path.join(MQ_TXT_DIR, 'evidence.txt'), sep = '\t')
evidence = evidence[evidence['Modifications'] == 'Unmodified']

evidence = evidence[['Sequence', 'Charge']]
evidence['collision_energy'] = 30
evidence = evidence.drop_duplicates()

evidence = evidence.rename(columns = {
    'Sequence': 'modified_sequence',
    'Charge': 'precursor_charge',
})
evidence.to_csv('prosit.csv', index = False)
