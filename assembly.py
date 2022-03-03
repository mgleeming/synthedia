import os, sys, time, copy, pickle, argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyopenms import *

from pyteomics import mass, fasta
from numba import jit

parser = argparse.ArgumentParser(
    description = 'Generate DIA data from DDA MaxQuant output.')

parser.add_argument( '--mq_txt_dir', required = False, type = str,
                    help = 'Path to MaxQuat "txt" directory')
parser.add_argument( '--prosit', required = False, type = str,
                    help = 'Path to prosit prediction library')
parser.add_argument( '--use_existing_peptide_file', required = False, type = str,
                    help = 'Use existing peptide definition file')

parser.add_argument( '--ms1_min_mz', required = False, type = float, default = 350,
                    help = 'Minimum m/z at MS1 level')
parser.add_argument( '--ms1_max_mz', required = False, type = float, default = 1600,
                    help = 'Maximum m/z at MS1 level')

parser.add_argument( '--ms2_min_mz', required = False, type = float, default = 100,
                    help = 'Minimum m/z at MS2 level')
parser.add_argument( '--ms2_max_mz', required = False, type = float, default = 2000,
                    help = 'Maximum m/z at MS2 level')

parser.add_argument( '--ms1_stdev', required = False, type = float, default = 0.02,
                    help = 'Mass spectral peak standard deviation at MS1 level' )
parser.add_argument( '--ms2_stdev', required = False, type = float, default = 0.05,
                    help = 'Mass spectral peak standard deviation at MS2 level')

parser.add_argument( '--rt_stdev', required = False, type = float, default = 2,
                    help = 'Chromatographic peak standard deviation')

parser.add_argument( '--ms1_point_diff', required = False, type = float, default = 0.002,
                    help = 'MS1 m/z difference between adjacent point')
parser.add_argument( '--ms2_point_diff', required = False, type = float, default = 0.01,
                    help = 'MS2 m/z difference between adjacent point')

parser.add_argument( '--ms1_scan_duration', required = False, type = float, default = 0.37,
                    help = 'Time in seconds taken to record an MS1 scan.')
parser.add_argument( '--ms2_scan_duration', required = False, type = float, default = 0.037,
                    help = 'Time in seconds taken to record an MS2 scan.')

parser.add_argument( '--original_run_length', required = False, type = float, default = 120,
                    help = 'Length in minutes of original data file. If not given, this will be determined by taking the difference between the minimum and maximum peptide retention times.')

parser.add_argument( '--new_run_length', required = False, type = float, default = 12,
                    help = 'Length in minutes of new data file.')
parser.add_argument( '--ms_clip_window', required = False, type = float, default = 0.5,
                    help = 'm/z window surrounding an MS peak that should be considered when simulating peak intensities. For high resolution data, this normally does not need to be changed.')
parser.add_argument( '--min_peak_fraction', required = False, type = float, default = 0.01,
                    help = 'Peptide elution profiles are simulated as gaussian peaks. This value sets the minimum gaussian curve intensitiy for a peptide to be simulated.')
parser.add_argument( '--mq_pep_threshold', required = False, type = float, default = 0.001,
                    help = 'For MaxQuant input data, use only peptides with a Posterior Error Probability (PEP) less than this value')
parser.add_argument( '--output_label', required = False, type = str, default = 'output',
                    help = 'Prefix for output files')
parser.add_argument( '--isolation_window', required = False, type = int, default = 30,
                    help = 'Length of DIA window in m/z')
parser.add_argument( '--write_empty_spectra', action = 'store_true',
                    help = 'Write empty mass sepctra to the output data file')
parser.add_argument( '--run_diann', action = 'store_true',
                    help = 'Run DIA-NN on the output data file')
parser.add_argument( '--diann_path', required = False, type = str, default = '/usr/diann/1.8/diann-1.8',
                    help = 'Path to DIA-NN')
parser.add_argument( '--out_dir', required = False, type = str, default = os.getcwd(),
                    help = 'Output directory where results should be written')

parser.add_argument( '--write_protein_fasta', action = 'store_true',
                    help = 'Write FASTA file with protein sequences for simulated peptides. If given, a FASTA file must be supplied with the --fasta options')
parser.add_argument( '--fasta', required = False, type = str,
                    help = 'Path to FASTA file from which protein sequences should be taken')
parser.add_argument( '--decoys', required = False, type = int, default = 200,
                    help = 'Write additional non-target protein sequences to output FASTA file')

options =  parser.parse_args()

if not any([options.mq_txt_dir, options.prosit]):
    print('Either an MaxQuant output directory or Prosit file is required')
    print('Exiting')
    sys.exit()

if options.write_protein_fasta == True and options.fasta == None:
    print('Synthedia was asked to write a FASTA file but no input FASTA was given')
    print('Exiting')
    sys.exit()

USE_EXISTING_PEPTIDE_FILE = options.use_existing_peptide_file
WRITE_EMPTY_SPECTRA = options.write_empty_spectra

DIA_NN_PATH = options.diann_path
RUN_DIANN = options.run_diann

MQ_TXT_DIR = options.mq_txt_dir
PROSIT_FILE = options.prosit

# constants
PROTON = 1.007276
IAA = 57.02092

FASTA_FILE = options.fasta
WRITE_PROTEIN_FASTA = options.write_protein_fasta

PEPTIDE_PEP_THRESHOLD = options.mq_pep_threshold
MS1_SCAN_RANGE = [options.ms1_min_mz, options.ms1_max_mz]
MS2_SCAN_RANGE = [options.ms2_min_mz, options.ms2_max_mz]

MS1_SCAN_DURATION = options.ms1_scan_duration
MS2_SCAN_DURATION = options.ms2_scan_duration

ORIGINAL_RUN_LENGTH = options.original_run_length * 60
NEW_RUN_LENGTH = options.new_run_length * 60

MIN_PEAK_FRACTION = options.min_peak_fraction

MS1_STDEV = options.ms1_stdev
MS2_STDEV = options.ms2_stdev
RT_STDEV = options.rt_stdev

MS1_MZS = np.arange(MS1_SCAN_RANGE[0], MS1_SCAN_RANGE[1], options.ms1_point_diff)
MS1_INTS = np.zeros(len(MS1_MZS))

MS2_MZS = np.arange(MS2_SCAN_RANGE[0], MS2_SCAN_RANGE[1], options.ms2_point_diff)
MS2_INTS = np.zeros(len(MS2_MZS))

WINDOW = options.ms_clip_window
ISOLATION_WINDOW = options.isolation_window

DECOY_NUMBER = options.decoys
out_dir = os.path.join( options.out_dir, 'run_length_%s' %str(int(options.new_run_length)))

try:
    os.makedirs(out_dir)
except:
    pass

output_label = os.path.join(out_dir, options.output_label + '_' + str(int(options.new_run_length)))
print('Writing outputs to %s' %output_label)

@jit(nopython=True)
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

class MS_run(object):

    def __init__(self):
        self.consumer = PlainMSDataWritingConsumer( '%s.mzML' %output_label)
        return

    def write_spec(self, spec):
        spec_to_write = MSSpectrum()
        spec_to_write.setRT(spec.rt)

        spec_to_write.setMSLevel(spec.order)

        if spec.order == 2:
            p = Precursor()
            p.setMZ((spec.isolation_hl + spec.isolation_ll) / 2)
            p.setIsolationWindowLowerOffset(spec.isolation_ll)
            p.setIsolationWindowUpperOffset(spec.isolation_hl)
            spec_to_write.setPrecursors( [p] )

        mask = np.where(spec.ints > 0)
        ints = spec.ints[mask]
        mzs = spec.mzs[mask]

        if WRITE_EMPTY_SPECTRA == False:
            if len(ints) > 0:
                spec_to_write.set_peaks([mzs, ints])

                self.consumer.consumeSpectrum(spec_to_write)
        else:
            spec_to_write.set_peaks([mzs, ints])
            self.consumer.consumeSpectrum(spec_to_write)

        del spec_to_write
        return

    def close(self):
        del self.consumer
        return

class Spectrum(object):
    def __init__(self, rt, order, isolation_range):
        self.rt = rt
        self.order = order

        if isolation_range:
            self.isolation_ll = isolation_range[0]
            self.isolation_hl = isolation_range[1]

        return

    def make_spectrum(self):
        if self.order == 1:
            self.mzs = MS1_MZS
            self.ints = copy.deepcopy(MS1_INTS)
        else:
            self.mzs = MS2_MZS
            self.ints = copy.deepcopy(MS2_INTS)
        return

    def add_peaks(self, peptide_scaled_rt, peaks):

        # scaling factor for point on chromatogram
        intensity_scale_factor = gaussian(self.rt, peptide_scaled_rt, RT_STDEV)

        if intensity_scale_factor < MIN_PEAK_FRACTION: return

        if self.order == 1:
            stdev = MS1_STDEV
        else:
            stdev = MS2_STDEV

        for peak in peaks:

            mz_mask = np.where((self.mzs > peak[0] - WINDOW) & (self.mzs < peak[0] + WINDOW))

            peak_ints = gaussian(self.mzs[mz_mask], peak[0], stdev)
            factor = peak[1] * intensity_scale_factor
            peak_ints *= factor
            int_mask = np.where(peak_ints < MIN_PEAK_FRACTION)
            peak_ints[int_mask] = 0
            self.ints[mz_mask] += peak_ints

        return

    def clear(self):
        del self.mzs
        del self.ints
        return

class Peptide(object):

    def __init__(self, evidence_entry = None, msms_entry = None, prosit_entry = None):

        if (evidence_entry is not None) and (msms_entry is not None):
            self.populate_from_mq(evidence_entry, msms_entry)
        elif prosit_entry is not None:
            self.populate_from_prosit(prosit_entry)
        else:
            print('Insufficient data to construct peptides. Exiting')
            sys.exit()
        return

    def populate_from_prosit(self, prosit_entry):
        self.sequence = prosit_entry['StrippedPeptide'].iloc[0]
        self.charge = prosit_entry['PrecursorCharge'].iloc[0]

        self.mass = (prosit_entry['PrecursorMz'].iloc[0] * self.charge) - (self.charge * PROTON)
        self.mz = prosit_entry['PrecursorMz'].iloc[0]

        self.rt = prosit_entry['iRT'].iloc[0] * 60
        self.protein = 'None'

        self.intensity = 1000000

        self.ms2_mzs = [float(_) for _ in prosit_entry['FragmentMz'].to_list()]
        self.ms2_ints = [float(_)*self.intensity for _ in prosit_entry['RelativeIntensity'].to_list()]
        self.ms2_peaks = list(zip(self.ms2_mzs, self.ms2_ints))

        self.ms1_isotopes = self.get_ms1_isotope_pattern()
        self.scaled_rt = self.scale_retention_times()
        return

    def populate_from_mq(self, evidence_entry, msms_entry):
        self.sequence = evidence_entry['Sequence']
        self.charge = evidence_entry['Charge']
        self.mass = evidence_entry['Mass']
        self.mz = evidence_entry['m/z']
        self.rt = evidence_entry['Retention time'] * 60
        self.intensity = evidence_entry['Intensity']
        self.protein = evidence_entry['Proteins']

        # MS2 fragment intensities are substantiall lower than precursor intensity (which is integrated)
        # Should these be scaled?
        self.ms2_mzs = [float(_) for _ in msms_entry['Masses'].iloc[0].split(';')]
        self.ms2_ints = [float(_) for _ in msms_entry['Intensities'].iloc[0].split(';')]
        self.ms2_peaks = list(zip(self.ms2_mzs, self.ms2_ints))

        assert len(self.ms2_mzs) == len(self.ms2_ints)

        self.ms1_isotopes = self.get_ms1_isotope_pattern()
        self.scaled_rt = self.scale_retention_times()

        self.evidence_entry = evidence_entry
        self.msms_entry = msms_entry
        return

    def get_ms1_isotope_pattern(self):
        isotopes = []
        for isotope in mass.mass.isotopologues( sequence = self.sequence, report_abundance = True, overall_threshold = 0.01, elements_with_isotopes = ['C']):
            calc_mz = (isotope[0].mass() + (IAA * self.sequence.count('C')) +  (PROTON * self.charge)) / self.charge
            isotopes.append( [calc_mz, isotope[1] * self.intensity])
        return isotopes

    def scale_retention_times(self):
        return self.rt / ORIGINAL_RUN_LENGTH * NEW_RUN_LENGTH

def make_spectra(run_length):
    '''
    Constructs a list of dicts that define the parameters of mass spectra to be simulated
    '''
    run_template = [
        {'order': 1, 'length': MS1_SCAN_DURATION, 'isolation_range': None},
    ]

    for i in range(MS1_SCAN_RANGE[0], MS1_SCAN_RANGE[1], ISOLATION_WINDOW):
        run_template.append({
            'order': 2, 'length': MS2_SCAN_DURATION, 'isolation_range': [i, i + ISOLATION_WINDOW]
        })

    spectra = []
    total_run_time = 0
    while total_run_time < run_length:
        for entry in run_template:
            spectra.append(
                Spectrum( total_run_time, entry['order'], entry['isolation_range'])
            )
            total_run_time += entry['length']

    return spectra


def read_peptides_from_prosit():

    # read inputs
    prosit = pd.read_csv(PROSIT_FILE, sep = ',')

    peptides = []
    l = len(prosit.groupby(['StrippedPeptide','PrecursorCharge']))
    for i, (index, precursor) in enumerate(prosit.groupby(['StrippedPeptide','PrecursorCharge'])):
        if i % 100 == 0:
            print('\t Constructing peptide %s of %s' %(i, l))

        peptides.append( Peptide(prosit_entry = precursor))
        break
    print('\tFinished constructing %s peptides' %(len(peptides)))
    return peptides

def read_peptides_from_mq():

    # read inputs
    msms = pd.read_csv(os.path.join(MQ_TXT_DIR, 'msms.txt'), sep = '\t')
    evidence = pd.read_csv(os.path.join(MQ_TXT_DIR, 'evidence.txt'), sep = '\t')

    evidence = evidence.sort_values(by = ['PEP'])
    for filterTerm in ['REV_', 'CON_']:
        evidence = evidence.loc[-evidence['Proteins'].str.contains(filterTerm, na=False)]

    evidence = evidence[evidence['PEP'] < PEPTIDE_PEP_THRESHOLD]

    # restrict to unmodified peptides for the moment
    # NB - carbamidomethyl cys is retained
    evidence = evidence[evidence['Modifications'] == 'Unmodified']

    peptides = []
    for i, evidence_row in evidence.iterrows():

        if i % 100 == 0:
            print('\t Constructing peptide %s of %s' %(i, len(evidence)))

        # filter non quantified peptides
        if np.isnan(evidence_row['Intensity']): continue

        # find matching msms entry - this cintains mz2 fragments
        msms_entry = msms[msms['id'] == evidence_row['Best MS/MS']]

        # safety
        assert evidence_row['Sequence'] == msms_entry['Sequence'].iloc[0]

        peptides.append(
            Peptide(evidence_entry = evidence_row, msms_entry = msms_entry)
        )
    print('\tFinished constructing %s peptides' %(len(peptides)))
    return peptides

#@profile
def populate_spectra(peptides, spectra):

    run = MS_run()
    t1 = time.time()
    for si, s in enumerate(spectra):

        if si % 1000 == 0:
            print('\t Writing spectrum %s of %s' %(si, len(spectra)))

        # make spec numpy arrays on the fly to sav mem
        s.make_spectrum()

        peptide_subset = [p for p in peptides if abs(p.scaled_rt - s.rt) < 30]

        for p in peptide_subset:
            if s.order == 1:
                # adds peptides MS1 isotope envelopes
                s.add_peaks(p.scaled_rt, p.ms1_isotopes)
            elif s.order == 2:
                if (p.mz > s.isolation_ll) and (p.mz < s.isolation_hl):
                    s.add_peaks(p.scaled_rt, p.ms2_peaks)

        # write final spec to file
        run.write_spec(s)

        # this deletes spec arrays that aren't needed anymore
        s.clear()

    # close consumer
    run.close()
    return

def write_peptide_target_table(peptides, spectra):

    def get_peptide_elution_window(p, ms1_rts):

        ints = gaussian(ms1_rts, p.scaled_rt, RT_STDEV)
        ints *= p.intensity
        mask = np.where(ints > MIN_PEAK_FRACTION)
        peak_rts = ms1_rts[mask]

        return min(peak_rts), max(peak_rts)

    ms1_rts = np.asarray([s.rt for s in spectra])

    of1 = open('%s_peptide_table.tsv' %output_label,'wt')
    of1.write('%s\n' %'\t'.join(['Protein', 'Sequence', 'Intensity', 'm/z', 'Charge', 'Mass', 'Experimental RT', 'Synthetic RT', 'Synthetic RT Start', 'Synthetic RT End', 'Synthetic m/z 0']))
    for p in peptides:

        rt_min, rt_max = get_peptide_elution_window(p, ms1_rts)

        of1.write('%s\n' %'\t'.join([str(_) for _ in [
            p.protein,
            p.sequence,
            p.intensity,
            p.mz,
            p.charge,
            p.mass,
            p.rt,
            p.scaled_rt,
            '%.3f' %rt_min,
            '%.3f' %rt_max,
            p.ms1_isotopes[0][0]
        ]]))

    of1.close()
    return

def write_target_protein_fasta(peptides):

    if not FASTA_FILE:
        print('No fasta file supplied - no proteins written')
        return

    decoy_counter = 0
    sequences = []
    with fasta.read(FASTA_FILE, use_index = True) as db:
        for i, entry in enumerate(db):
            for p in peptides:
                if p.sequence in entry.sequence:
                    sequences.append([ entry.description, entry.sequence ])
            else:
                if decoy_counter < DECOY_NUMBER:
                    sequences.append([ entry.description, entry.sequence ])
                    decoy_counter += 1

    print('Writing %s sequences including %s decoys to fasta' %(len(sequences), decoy_counter))
    fasta.write(sequences, output = '%s.fasta' %output_label, file_mode = 'w')
    return

def run_diann(input_dir):

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

    os.system(cmd)
    return

def main(options):

    print('Started Synthedia')
    t1 = time.time()

    print('Preparing spectral template')
    spectra = make_spectra(NEW_RUN_LENGTH)

    print('Constructing peptide models')
    if (USE_EXISTING_PEPTIDE_FILE == True) and (os.path.isfile(os.path.join(options.out_dir, 'peptides.pickle')) == True):
        print('Using existing peptide file')
        with open( os.path.join(options.out_dir, 'peptides.pickle') , 'rb') as handle:
            peptides = pickle.load(handle)
    else:
        if MQ_TXT_DIR:
            peptides = read_peptides_from_mq()
        elif PROSIT_FILE:
            peptides = read_peptides_from_prosit()
        with open( os.path.join(options.out_dir, 'peptides.pickle') , 'wb') as handle:
            pickle.dump(peptides, handle, protocol=pickle.HIGHEST_PROTOCOL)

    for p in peptides:
        p.scaled_rt = p.scale_retention_times()

    print('Writing peptide target table')
    write_peptide_target_table(peptides, spectra)

    print('Writing peptides to spectra')
    populate_spectra(peptides, spectra)

    if WRITE_PROTEIN_FASTA:
        print('Writing protein fasta file')
        write_target_protein_fasta(peptides)

    if RUN_DIANN:
        print('Running DIA-NN')
        run_diann(out_dir)

    print('Done!')
    print('Total execution time: %s' %(time.time() - t1))
    return

if __name__ == '__main__':
    options =  parser.parse_args()
    main(options)
