import os, sys, time, copy, pickle, argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyopenms import *

from pyteomics import mass, fasta
from numba import jit

parser = argparse.ArgumentParser(
    description = 'Generate DIA data from DDA MaxQuant output.'
)

parser.add_argument( '--mq_txt_dir', required = False, type = str, default = '.')
parser.add_argument( '--use_pickled_peptides', required = False, type = str,)
parser.add_argument( '--ms1_min_mz', required = False, type = float, default = 350)
parser.add_argument( '--ms1_max_mz', required = False, type = float, default = 1600)
parser.add_argument( '--ms2_min_mz', required = False, type = float, default = 100)
parser.add_argument( '--ms2_max_mz', required = False, type = float, default = 2000)
parser.add_argument( '--ms1_stdev', required = False, type = float, default = 0.02)
parser.add_argument( '--ms2_stdev', required = False, type = float, default = 0.05)
parser.add_argument( '--rt_stdev', required = False, type = float, default = 2)
parser.add_argument( '--ms1_point_diff', required = False, type = float, default = 0.002)
parser.add_argument( '--ms2_point_diff', required = False, type = float, default = 0.01)
parser.add_argument( '--ms1_scan_duration', required = False, type = float, default = 0.37)
parser.add_argument( '--ms2_scan_duration', required = False, type = float, default = 0.037)
parser.add_argument( '--original_run_length', required = False, type = float, default = 120)
parser.add_argument( '--new_run_length', required = False, type = float, default = 12)
parser.add_argument( '--ms_clip_window', required = False, type = float, default = 0.5)
parser.add_argument( '--min_peak_fraction', required = False, type = float, default = 0.01)
parser.add_argument( '--output_label', required = False, type = str, default = 'output')
parser.add_argument( '--isolation_window', required = False, type = int, default = 30)

options =  parser.parse_args()

MQ_TXT_DIR = options.mq_txt_dir

# constants
PROTON = 1.007276
IAA = 57.02092

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

out_dir = 'run_length_%s' %str(options.new_run_length)
os.makedirs(out_dir)
output_label = os.path.join(out_dir, options.output_label + '_' + str(options.new_run_length))

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

        if len(ints) > 0:
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
    def __init__(self, evidence_entry, msms_entry):

        self.sequence = evidence_entry['Sequence']
        self.charge = evidence_entry['Charge']
        self.mass = evidence_entry['Mass']
        self.mz = evidence_entry['m/z']
        self.rt = evidence_entry['Retention time'] * 60
        self.intensity = evidence_entry['Intensity']

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
            isotopes.append(
                [calc_mz, isotope[1] * self.intensity]
            )
        return isotopes

    def scale_retention_times(self):
        return self.rt / ORIGINAL_RUN_LENGTH * NEW_RUN_LENGTH

def make_spectra(run_length):

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


def read_peptides():

    # read inputs
    msms = pd.read_csv(os.path.join(MQ_TXT_DIR, 'msms.txt'), sep = '\t')
    evidence = pd.read_csv(os.path.join(MQ_TXT_DIR, 'evidence.txt'), sep = '\t')
    proteins = pd.read_csv(os.path.join(MQ_TXT_DIR, 'proteinGroups.txt'), sep = '\t')

    for filterTerm in ['REV_', 'CON_']:
        proteins = proteins.loc[-proteins['Protein IDs'].str.contains(filterTerm, na=False)]

    proteins = proteins.sort_values(by=['Razor + unique peptides'], ascending=False)

    selected_proteins = proteins.head(100)['Protein IDs'].to_list()

    # restrict to unmodified peptides for the moment
    # NB - carbamidomethyl cys is retained
    evidence = evidence[evidence['Modifications'] == 'Unmodified']

    peptides = []
    for i, evidence_row in evidence.iterrows():

        if evidence_row['Leading razor protein'] not in selected_proteins: continue

        if i % 100 == 0:
            print('\t Constructing peptide %s of %s' %(i, len(evidence)))

        # filter non quantified peptides
        if np.isnan(evidence_row['Intensity']): continue

        # find matching msms entry - this cintains mz2 fragments
        msms_entry = msms[msms['id'] == evidence_row['Best MS/MS']]

        # safety
        assert evidence_row['Sequence'] == msms_entry['Sequence'].iloc[0]

        peptides.append(
            Peptide(evidence_row, msms_entry)
        )

    print('\tFinished constructing %s peptides' %(len(peptides)))
    return peptides

#@profile
def populate_spectra(peptides, spectra):

    run = MS_run()
    t1 = time.time()
    for si, s in enumerate(spectra):

        if si % 100 == 0:
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
        ints *= p.evidence_entry['Intensity']
        mask = np.where(ints > MIN_PEAK_FRACTION)
        peak_rts = ms1_rts[mask]

        return min(peak_rts), max(peak_rts)

    ms1_rts = np.asarray([s.rt for s in spectra])

    of1 = open('%s_peptide_table.tsv' %output_label,'wt')
    of1.write('%s\n' %'\t'.join(['Protein', 'Sequence', 'Intensity', 'm/z', 'Charge', 'Mass', 'Experimental RT', 'Synthetic RT', 'Synthetic RT Start', 'Synthetic RT End', 'Synthetic m/z 0']))
    for p in peptides:

        rt_min, rt_max = get_peptide_elution_window(p, ms1_rts)

        of1.write('%s\n' %'\t'.join([str(_) for _ in [
            p.evidence_entry['Proteins'],
            p.sequence,
            p.evidence_entry['Intensity'],
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

    selected_proteins = []
    for p in peptides:
        for prot in p.evidence_entry['Proteins'].split(';'):
            selected_proteins.append(prot)
    selected_proteins = list(set(selected_proteins))

    sequences = []
    with fasta.read('UNIPROT_HUMAN_REVIEWED_2020.fasta') as db:
        for i, entry in enumerate(db):
            for s in selected_proteins:
                if s in entry.description:
                    sequences.append([ entry.description, entry.sequence ])
    fasta.write(sequences, output = '%s.fasta' %output_label, file_mode = 'wt')
    return

def main():

    t1 = time.time()

    print('Preparing spectral template')
    spectra = make_spectra(NEW_RUN_LENGTH)

    print('Constructing peptide models')
    get_peptides = False
    if get_peptides:
        peptides = read_peptides()
        with open('peptides.pickle', 'wb') as handle:
            pickle.dump(peptides, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open('peptides.pickle', 'rb') as handle:
            peptides = pickle.load(handle)

    for p in peptides:
        p.scaled_rt = p.scale_retention_times()

    print('Writing peptide target table')
    write_peptide_target_table(peptides, spectra)

    print('Writing protein fasta file')
    write_target_protein_fasta(peptides)

    print('Writing peptides to spectra')
    populate_spectra(peptides, spectra)

    print('Done!')
    print('Total execution time: %s' %(time.time() - t1))
    return

if __name__ == '__main__':
    main()
