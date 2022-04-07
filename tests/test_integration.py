import pytest, os, yaml, copy, shutil

import pandas as pd
import numpy as np

from synthedia import assembly
from synthedia.mzml import MZMLReader

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_RESOURCES = os.path.join(ROOT_DIR, 'resources')
TEST_OUTPUTS = os.path.join(ROOT_DIR, 'test_outputs')

def create_test_dir():
    try:
        shutil.rmtree(TEST_OUTPUTS)
    except:
        pass
    try:
        os.makedirs(TEST_OUTPUTS)
    except:
        pass

PARAM_TEMPLATE = yaml.load(open(os.path.join(TEST_RESOURCES, 'default_simulation_args.yaml')), Loader=yaml.FullLoader)
PARAM_TEMPLATE['out_dir'] = TEST_OUTPUTS

class Options():
    def __init__(self, d):
        for k,v in d.items(): setattr(self,k,v)

def read_mzml_file(mzml_file):
    n_ms1, n_ms2 = 0, 0
    for (rt, lvl, mzs, ints) in MZMLReader(mzml_file):
        if lvl == 1:
            n_ms1 += 1
        if lvl == 2:
            n_ms2 += 2
    return n_ms1, n_ms2

def plot_MS1_EICe(mzml_file):
    ll, hl = get_precursor()
    rts, intensities, zero_intensities = [], [], []
    for (rt, lvl, mzs, ints) in MZMLReader(mzml_file):
        if lvl != 1: continue
        mask = np.where((mzs > ll) & (mzs < hl))
        rts.append(rt)
        intensities.append(ints[mask].sum())

        # check for 0 intensity everywhere else
        zero_mask = np.zeros(mzs.shape, dtype='bool')
        zero_mask[:] = True
        zero_mask[mask] = False
        zero_intensities.append(ints[zero_mask].sum())

    assert sum(intensities) > 0
    assert sum(zero_intensities) == 0

def get_precursor():
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    mz = peptide_table['m/z'].iloc[0]
    ll = mz - 1
    hl = mz + 5
    return ll, hl

def update_param(new_params):

    params = copy.deepcopy(PARAM_TEMPLATE)
    for param, value in new_params.items():
        params[param] = value

    return Options(params)

def check_files():
    assert os.path.isfile(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML')), 'mzML file not produced'
    assert os.path.getsize(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML')) > 0, 'mzML file has 0 size'
    assert os.path.isfile(os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv')), 'Peptide table not produced'
    assert os.path.getsize(os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv')) > 0, 'Peptide table has 0 size'

def test_maxquant_default():
    create_test_dir()
    options = update_param({'mq_txt_dir': TEST_RESOURCES})
    assembly.assemble(options)
    check_files()
    n_ms1, n_ms2 = read_mzml_file(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    assert n_ms1 > 0
    assert n_ms2 > 0

def test_prosit_defaults():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv')})
    assembly.assemble(options)
    check_files()
    n_ms1, n_ms2 = read_mzml_file(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    assert n_ms1 > 0
    assert n_ms2 > 0

def test_acquisition_schema():
    create_test_dir()
    options = update_param({'mq_txt_dir': TEST_RESOURCES, 'acquisition_schema' : os.path.join(TEST_RESOURCES, 'acquisition_schema_example.csv')})
    assembly.assemble(options)
    check_files()
    n_ms1, n_ms2 = read_mzml_file(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    assert n_ms1 > 0
    assert n_ms2 > 0

def test_acquisition_schema_no_MS1():
    create_test_dir()
    options = update_param({'mq_txt_dir': TEST_RESOURCES, 'acquisition_schema' : os.path.join(TEST_RESOURCES, 'acquisition_schema_example_no_MS1.csv')})
    assembly.assemble(options)
    check_files()
    n_ms1, n_ms2 = read_mzml_file(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    assert n_ms1 == 0
    assert n_ms2 > 0

def test_no_peptides_in_MS_range():
    with pytest.raises(Exception):
        create_test_dir()
        options = update_param({'mq_txt_dir': TEST_RESOURCES, 'ms2_max_mz': 200, 'ms1_max_mz': 361})
        assembly.assemble(options)

def test_invalid_input():
    with pytest.raises(Exception):
        create_test_dir()
        options = update_param({})
        assembly.assemble(options)

def test_invalid_rt_peak_model():
    with pytest.raises(Exception):
        create_test_dir()
        options = update_param({'mq_txt_dir': TEST_RESOURCES, 'rt_peak_model': 'aaaa'})
        assembly.assemble(options)

def test_invalid_mz_peak_model():
    with pytest.raises(Exception):
        create_test_dir()
        options = update_param({'mq_txt_dir': TEST_RESOURCES, 'mz_peak_model': 'aaaa'})
        assembly.assemble(options)

def test_MS1_EICs():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid': True})
    assembly.assemble(options)
    plot_MS1_EICe(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
