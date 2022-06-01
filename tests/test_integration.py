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

def get_TIC_for_ms_lvl(mzml_file, target_lvl = None):
    rts, intensities, zero_intensities = [], [], []
    for (rt, lvl, mzs, ints) in MZMLReader(mzml_file):
        if target_lvl:
            if lvl != target_lvl: continue
        if ints.sum() > 0:
            rts.append(rt)
            intensities.append(ints.sum())
    return rts, intensities

def get_EIC_for_ms_lvl(mzml_file, target_lvl = None):
    rts, intensity_sums, intensity_maxes= [],[],[]

    if target_lvl == 1:
        ll, hl = get_precursor_10ppm()
    if target_lvl == 2:
        ll, hl = get_most_abundant_fragment()
    for (rt, lvl, mzs, ints) in MZMLReader(mzml_file):
        if lvl != target_lvl: continue

        mask = np.where((mzs > ll) & (mzs < hl))
        mzs = mzs[mask]
        ints = ints[mask]

        if ints.sum() > 0:
            rts.append(rt)
            intensity_maxes.append(ints.max())
            intensity_sums.append(ints.sum())

    return rts, intensity_sums, intensity_maxes

def get_precursor_10ppm():
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    mz = peptide_table['m/z'].iloc[0]
    diff = mz / 1000000 * 10
    ll = mz - diff
    hl = mz + diff + 5
    return ll, hl

def get_most_abundant_fragment():
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    mz = peptide_table['Synthetic max fragment m/z'].iloc[0]
    diff = mz / 1000000 * 50
    ll = mz - diff
    hl = mz + diff
    return ll, hl

def get_precursor():
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    mz = peptide_table['m/z'].iloc[0]
    ll = mz - 1
    hl = mz + 5
    return ll, hl

def get_n_chromatographic_points_for_ms_level(lvl):
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    npoints = peptide_table['MS%s chromatographic points group_0_sample_0'% lvl].iloc[0]
    return npoints

def get_target_peak_area_and_height():
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    area = peptide_table['Total precursor abundance group_0_sample_0'].iloc[0]
    height = peptide_table['Precursor max peak height group_0_sample_0'].iloc[0]
    return area, height

def get_fragment_peak_area_and_height():
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    area = peptide_table['Most abundant fragment total intensity group_0_sample_0'].iloc[0]
    height = peptide_table['Most abundant fragment max height group_0_sample_0'].iloc[0]
    return area, height

def get_target_peak_retention_boundaries():
    peptide_table = pd.read_csv( os.path.join(TEST_OUTPUTS, 'output_peptide_table.tsv') , sep = '\t')
    start = peptide_table['Synthetic RT Start'].iloc[0]
    end = peptide_table['Synthetic RT End'].iloc[0]
    return start/60, end/ 60

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
#
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

def test_prosit_centroid_ms1():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True})
    assembly.assemble(options)
    check_files()
    n_ms1, n_ms2 = read_mzml_file(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    assert n_ms1 > 0
    assert n_ms2 > 0

def test_prosit_centroid_ms2():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms2': True})
    assembly.assemble(options)
    check_files()
    n_ms1, n_ms2 = read_mzml_file(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    assert n_ms1 > 0
    assert n_ms2 > 0

def test_prosit_centroid_both():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'centroid_ms2': True})
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
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True})
    assembly.assemble(options)
    plot_MS1_EICe(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))

def test_n_MS_points():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True})
    assembly.assemble(options)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 2)
    n_pionts = get_n_chromatographic_points_for_ms_level(2)
    assert n_pionts == len(rts)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 1)
    n_pionts = get_n_chromatographic_points_for_ms_level(1)
    assert n_pionts == len(rts)

def test_n_MS_points_narrow():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'rt_peak_fwhm': 0.3})
    assembly.assemble(options)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 2)
    n_pionts = get_n_chromatographic_points_for_ms_level(2)
    assert n_pionts == len(rts)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 1)
    n_pionts = get_n_chromatographic_points_for_ms_level(1)
    assert n_pionts == len(rts)

def test_n_MS_points_wide():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'rt_peak_fwhm': 10})
    assembly.assemble(options)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 2)
    n_pionts = get_n_chromatographic_points_for_ms_level(2)
    assert n_pionts == len(rts)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 1)
    n_pionts = get_n_chromatographic_points_for_ms_level(1)
    assert n_pionts == len(rts)

def test_peptide_total_intensity():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'rt_peak_fwhm': 10})
    assembly.assemble(options)
    rts, intensity_sums, intensity_maxes = get_EIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 1)
    area, height = get_target_peak_area_and_height()

    assert abs (100 - sum(intensity_sums) /  area * 100) < 0.01
    assert abs (100 - max(intensity_maxes) / height * 100) < 0.01


def test_retention_lengths():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True})
    assembly.assemble(options)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    start, end = get_target_peak_retention_boundaries()
    assert np.isclose(rts[0], start)
    assert np.isclose(rts[-1], end)

def test_retention_lengths_narrow():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'rt_peak_fwhm': 1})
    assembly.assemble(options)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    start, end = get_target_peak_retention_boundaries()
    assert np.isclose(rts[0], start)
    assert np.isclose(rts[-1], end)

def test_retention_lengths_wide():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'rt_peak_fwhm': 10})
    assembly.assemble(options)
    rts, intensities = get_TIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'))
    start, end = get_target_peak_retention_boundaries()
    assert np.isclose(rts[0], start)
    assert np.isclose(rts[-1], end)

def test_missing_in_group():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'prob_missing_in_group': 100})
    assembly.assemble(options)
    start, end = get_target_peak_retention_boundaries()
    assert np.isclose(start, 0)
    assert np.isclose(end, 0)

def test_missing_in_sample():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'prob_missing_in_group': 100, 'prob_missing_in_sample': 0})
    assembly.assemble(options)
    start, end = get_target_peak_retention_boundaries()
    assert np.isclose(start, 0)
    assert np.isclose(end, 0)

def test_max_fragment_intensity():
    create_test_dir()
    options = update_param({'prosit': os.path.join(TEST_RESOURCES, 'myPrositLib.csv'), 'centroid_ms1': True, 'centroid_ms2': True})
    assembly.assemble(options)
    rts, intensity_sums, intensity_maxes = get_EIC_for_ms_lvl(os.path.join(TEST_OUTPUTS, 'output_group_0_sample_0.mzML'), target_lvl = 2)
    area, height = get_fragment_peak_area_and_height()
    assert abs (100 - sum(intensity_sums) /  area * 100) < 0.01
    assert abs (100 - max(intensity_maxes) / height * 100) < 0.01

