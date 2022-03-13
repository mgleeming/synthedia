import os, sys, time, copy, pickle, multiprocessing
import random, math, datetime
import pandas as pd
import numpy as np
from pyopenms import *
from pyteomics import mass, fasta
from numba import jit

from . import plotting

# TODO:jj
# 1) Acquisition schema import - Done
# 2) Spectrum centroiding - Donej
# 3) Automated determination of peak standard deviations - Done
# 4) Randomness in elution profiles - Done
# 6) Decoys?
# 7) Plots? - Done
# 8) Quantification between multiple files - Done
# 9) Peak tailing
# 10) Probability sample in missing
# 11) Probability group in missing
# 12) Probability missing increases as abundance decreases
# 13) add contaminant wall ions
# 14) Chemical noise
# 15) Missing value plots

# constants
PROTON = 1.007276
IAA = 57.02092

@jit(nopython=True)
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

class MZMLReader():
    def __init__(self, file):
        self.od_exp = OnDiscMSExperiment()
        self.od_exp.openFile(file)
        self.num_spectra = self.od_exp.getNrSpectra()

    def __iter__(self):
        self.current_spec = 0
        return self

    def __next__(self):
        if self.current_spec < self.num_spectra:
            s = self.od_exp.getSpectrum(self.current_spec)
            t = s.getRT() / 60
            lvl = s.getMSLevel()
            mzs, ints = s.get_peaks()
            self.current_spec += 1
            return t, lvl, mzs, ints
        else:
            self.close()
            raise StopIteration

    def close(self):
        del self.od_exp
        return

class MZMLWriter():

    def __init__(self, out_file):
        self.consumer = PlainMSDataWritingConsumer(out_file)
        self.n_spec_written = 0
        return

    def write_spec(self, spec):
        spec_to_write = MSSpectrum()

        mask = np.where(spec.ints > 0)
        ints = spec.ints[mask]
        mzs = spec.mzs[mask]

        if options.write_empty_spectra == False:
            if len(ints) > 0:
                spec_to_write.set_peaks([mzs, ints])
            else:
                del spec_to_write
                return
        else:
            spec_to_write.set_peaks([mzs, ints])

        if options.centroid:
            centroided_spectrum = MSSpectrum()
            PeakPickerHiRes().pick(spec_to_write, centroided_spectrum)
            spec_to_write = centroided_spectrum

        spec_to_write.setRT(spec.rt)
        spec_to_write.setMSLevel(spec.order)

        if spec.order == 2:
            p = Precursor()
            p.setMZ((spec.isolation_hl + spec.isolation_ll) / 2)
            p.setIsolationWindowLowerOffset(spec.isolation_ll)
            p.setIsolationWindowUpperOffset(spec.isolation_hl)
            spec_to_write.setPrecursors( [p] )

        self.consumer.consumeSpectrum(spec_to_write)
        self.n_spec_written += 1
        del spec_to_write
        return

    def close(self):
        del self.consumer
        return

class Spectrum():
    def __init__(self, rt, order, isolation_range):
        self.rt = rt
        self.order = order

        if isolation_range:
            self.isolation_ll = isolation_range[0]
            self.isolation_hl = isolation_range[1]

        return

    def make_spectrum(self, MS1_MZS, MS1_INTS, MS2_MZS, MS2_INTS):
        # much faster than creating a new array for every scan
        if self.order == 1:
            self.mzs = MS1_MZS
            self.ints = copy.deepcopy(MS1_INTS)
        else:
            self.mzs = MS2_MZS
            self.ints = copy.deepcopy(MS2_INTS)
        return

    def add_peaks(self, peptide_scaled_rt, peaks, abundance_offset):

        # scaling factor for point on chromatogram
        intensity_scale_factor = gaussian(self.rt, peptide_scaled_rt, options.rt_stdev)

        # apply chromatographic instability if needed
        if options.rt_instability > 0:
            instability_factor = 1 - ( random.randint(0,options.rt_instability * 100) / 10000 )
            intensity_scale_factor *= instability_factor

        # make sure peak decays to 0
        if intensity_scale_factor < options.min_peak_fraction: return

        if self.order == 1:
            stdev = options.ms1_stdev
        else:
            stdev = options.ms2_stdev

        for peak in peaks:

            # apply offset to peak intensity
            log_int = math.log2(peak[1])
            adjusted_log2_int = log_int + abundance_offset
            adjusetd_raw_int = 2 ** adjusted_log2_int

            # calculating the gaussian for the full m/z range is slow
            # subset data to only a small region around the peak to speed calculation
            mz_mask = np.where((self.mzs > peak[0] - options.ms_clip_window) & (self.mzs < peak[0] + options.ms_clip_window))
            peak_ints = gaussian(self.mzs[mz_mask], peak[0], stdev)

            # scale gaussian intensities by chromatogram scaling factor
            factor = adjusetd_raw_int * intensity_scale_factor
            peak_ints *= factor

            # remove low intensity points
            int_mask = np.where(peak_ints < options.min_peak_fraction)
            peak_ints[int_mask] = 0

            # add new data to full spectrum intensity
            self.ints[mz_mask] += peak_ints

        return

    def clear(self):
        # save memory
        del self.mzs
        del self.ints
        return

class Peptide():

    def __init__(self, evidence_entry = None, msms_entry = None, prosit_entry = None):

        if (evidence_entry is not None) and (msms_entry is not None):
            self.populate_from_mq(evidence_entry, msms_entry)
        elif prosit_entry is not None:
            self.populate_from_prosit(prosit_entry)
        else:
            print('Insufficient data to construct peptides. Exiting')
            sys.exit()

        self.ms1_isotopes = None
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
        return

    def populate_from_mq(self, evidence_entry, msms_entry):
        self.sequence = evidence_entry['Sequence']
        self.charge = evidence_entry['Charge']
        self.mass = evidence_entry['Mass']
        self.mz = evidence_entry['m/z']
        self.rt = evidence_entry['Retention time'] * 60
        self.intensity = evidence_entry['Intensity']
        self.protein = evidence_entry['Proteins']

        # MS2 fragment intensities are substantially lower than precursor intensity (which is integrated)
        # Should these be scaled?
        self.ms2_mzs = [float(_) for _ in msms_entry['Masses'].iloc[0].split(';')]
        self.ms2_ints = [float(_) for _ in msms_entry['Intensities'].iloc[0].split(';')]
        self.ms2_peaks = list(zip(self.ms2_mzs, self.ms2_ints))

        assert len(self.ms2_mzs) == len(self.ms2_ints)

        self.evidence_entry = evidence_entry
        self.msms_entry = msms_entry
        return

    def get_ms1_isotope_pattern(self):
        isotopes = []
        for isotope in mass.mass.isotopologues( sequence = self.sequence, report_abundance = True, overall_threshold = 0.01, elements_with_isotopes = ['C']):
            calc_mz = (isotope[0].mass() + (IAA * self.sequence.count('C')) +  (PROTON * self.charge)) / self.charge
            isotopes.append( [calc_mz, isotope[1] * self.intensity])
        self.ms1_isotopes = isotopes
        return

    def scale_retention_times(self, options):
        self.scaled_rt = self.rt / options.original_run_length * options.new_run_length
        return

    def set_abundances(self, group, sample, abundance_offset):

        if not hasattr(self, 'abundances'):
            self.abundances = []
            self.offsets = []

        # first sample of group - add new list
        if sample == 0:
            self.abundances.append([])
            self.offsets.append([])

        log_int = math.log2(self.intensity)
        adjusted_log2_int = log_int + abundance_offset
        adjusetd_raw_int = 2 ** adjusted_log2_int

        self.abundances[group].append(adjusetd_raw_int)
        self.offsets[group].append(abundance_offset)
        return

def make_spectra(options):
    '''
    Constructs a list of dicts that define the parameters of mass spectra to be simulated
    '''
    run_template = []
    if options.acquisition_schema is not None:
        try:
            df = pd.read_csv(options.acquisition_schema)
            for index, row in df.iterrows():
                run_template.append({
                    'order': int(row['ms_level']), 'length': float(row['scan_duration_in_seconds']),
                    'isolation_range': [float(row['isolation_window_lower_mz']), float(row['isolation_window_upper_mz'])]
                })
        except:
            print('Error parsing acquisition schema file')
            print('Exiting')
            sys.exit()
    else:
        run_template.append({
            'order': 1, 'length': options.ms1_scan_duration, 'isolation_range': None
        })
        for i in range(options.ms1_min_mz, options.ms1_max_mz, options.isolation_window):
            run_template.append({
                'order': 2, 'length': options.ms2_scan_duration, 'isolation_range': [i, i + options.isolation_window]
            })

    if (options.all == True) or (options.schema == True):
        plotting.plot_acquisition_schema(options, run_template)

    sys.exit()

    spectra = []
    total_run_time = 0
    while total_run_time < options.new_run_length:
        for entry in run_template:
            spectra.append(
                Spectrum( total_run_time, entry['order'], entry['isolation_range'])
            )
            total_run_time += entry['length']

    return spectra

def read_peptides_from_prosit(options):

    # read inputs
    prosit = pd.read_csv(options.prosit, sep = ',')

    peptides = []
    l = len(prosit.groupby(['StrippedPeptide','PrecursorCharge']))
    for i, (index, precursor) in enumerate(prosit.groupby(['StrippedPeptide','PrecursorCharge'])):
        if i % 100 == 0:
            print('\t Constructing peptide %s of %s' %(i, l))

        # check if precursor out of bounds
        if (float(precursor['PrecursorMz']) < options.ms1_min_mz) or (float(precursor['PrecursorMz']) > options.ms1_max_mz):
            continue

        peptides.append( Peptide(prosit_entry = precursor))
        break

    print('\tFinished constructing %s peptides' %(len(peptides)))
    return peptides

def read_peptides_from_mq(options):

    # read inputs
    msms = pd.read_csv(os.path.join(options.mq_txt_dir, 'msms.txt'), sep = '\t')
    evidence = pd.read_csv(os.path.join(options.mq_txt_dir, 'evidence.txt'), sep = '\t')

    evidence = evidence.sort_values(by = ['PEP'])
    for filterTerm in ['REV_', 'CON_']:
        evidence = evidence.loc[-evidence['Proteins'].str.contains(filterTerm, na=False)]

    evidence = evidence[evidence['PEP'] < options.mq_pep_threshold]

    # restrict to unmodified peptides for the moment
    # NB - carbamidomethyl cys is retained
    evidence = evidence[evidence['Modifications'] == 'Unmodified']

    counter = 0
    peptides = []
    for i, evidence_row in evidence.iterrows():

        counter += 1
        if counter % 100 == 0:
            print('\t Constructing peptide %s of %s' %(counter, len(evidence)))

        # check if precursor out of bounds
        if (float(evidence_row['m/z']) < options.ms1_min_mz) or (float(evidence_row['m/z']) > options.ms1_max_mz):
            continue

        # filter non quantified peptides
        if np.isnan(evidence_row['Intensity']): continue

        # find matching msms entry - this cintains mz2 fragments
        msms_entry = msms[msms['id'] == evidence_row['Best MS/MS']]

        # safety
        assert evidence_row['Sequence'] == msms_entry['Sequence'].iloc[0]

        peptides.append(
            Peptide(evidence_entry = evidence_row, msms_entry = msms_entry)
        )

        if len(peptides) == 50:
            break

    print('\tFinished constructing %s peptides' %(len(peptides)))
    return peptides

def generate_protein_abundances_for_mzml_files(options, peptides):

    # read mq input
    evidence = pd.read_csv(os.path.join(options.mq_txt_dir, 'evidence.txt'), sep = '\t')

    # determine protein abundances per treatment group
    group_abundance_dict = {p:[0] + [np.random.normal(loc = 0, scale = options.between_group_stdev) for _ in range(options.n_groups - 1)] for p in list(set(evidence['Proteins'].to_list()))}

    for p in peptides:
        group_abundances = group_abundance_dict.get(p.protein)
        for groupi, group_abundance in enumerate(group_abundances):
            sample_abundances = [np.random.normal(loc = group_abundance, scale = options.within_group_stdev) for _ in range(options.samples_per_group)]
            for samplei, sample_abundance in enumerate(sample_abundances):
                p.set_abundances(groupi, samplei, sample_abundance)

    return [p for p in peptides if p.abundances != None]

def calculate_scaled_retention_times(options, peptides):

    for p in peptides:
        p.scale_retention_times(options)

    return peptides

def populate_spectra(options, peptides, spectra, groupi, samplei):

    MS1_MZS = np.arange(options.ms1_min_mz, options.ms1_max_mz, options.ms1_point_diff)
    MS1_INTS = np.zeros(len(MS1_MZS))

    MS2_MZS = np.arange(options.ms2_min_mz, options.ms2_max_mz, options.ms2_point_diff)
    MS2_INTS = np.zeros(len(MS2_MZS))

    run = MZMLWriter(
        os.path.join(
            options.out_dir, '%s_group_%s_sample_%s.mzML' %(
                options.output_label, groupi, samplei
            )
        )
    )

    for spectrumi, spectrum in enumerate(spectra):

        if spectrumi % 1000 == 0:
            print('\tGroup %s, Sample %s - Writing spectrum %s of %s' %(groupi, samplei, spectrumi, len(spectra)))

        # make spec numpy arrays on the fly to sav mem
        spectrum.make_spectrum(MS1_MZS, MS1_INTS, MS2_MZS, MS2_INTS)

        peptide_subset = [p for p in peptides if abs(p.scaled_rt - spectrum.rt) < 30]

        for p in peptide_subset:

            if spectrum.order == 1:
                # adds peptides MS1 isotope envelopes
                abundance_offset = p.offsets[groupi][samplei]
                spectrum.add_peaks(p.scaled_rt, p.ms1_isotopes, abundance_offset)

            elif spectrum.order == 2:
                if (p.mz > spectrum.isolation_ll) and (p.mz < spectrum.isolation_hl):
                    abundance_offset = p.offsets[groupi][samplei]
                    spectrum.add_peaks(p.scaled_rt, p.ms2_peaks, abundance_offset)

        # write final spec to file
        run.write_spec(spectrum)

        # this deletes m/z and intensity arrays that aren't needed anymore
        spectrum.clear()

    # close consumer
    run.close()

    return

def write_peptide_target_table(options, peptides, spectra):

    def get_peptide_elution_window(p, ms1_rts):

        ints = gaussian(ms1_rts, p.scaled_rt, options.rt_stdev)
        ints *= p.intensity
        mask = np.where(ints > options.min_peak_fraction)
        peak_rts = ms1_rts[mask]

        return min(peak_rts), max(peak_rts)

    ms1_rts = np.asarray([s.rt for s in spectra])

    of1 = open( os.path.join(options.out_dir, '%s_peptide_table.tsv' %options.output_label),'wt')

    to_write = [
        'Protein',
        'Sequence',
        'Intensity',
        'm/z',
        'Charge',
        'Mass',
        'Experimental RT',
        'Synthetic RT',
        'Synthetic RT Start',
        'Synthetic RT End',
        'Synthetic m/z 0'
    ]
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Intensity group_%s_sample_%s' %(group, sample))
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Offset group_%s_sample_%s' %(group, sample))

    of1.write('%s\n' %'\t'.join([str(_) for _ in to_write]))

    for p in peptides:
        rt_min, rt_max = get_peptide_elution_window(p, ms1_rts)
        to_write = [
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
        ]
        for group in p.abundances:
            for sample in group:
                to_write.append(sample)
        for group in p.offsets:
            for sample in group:
                to_write.append(sample)

        of1.write('%s\n' %'\t'.join([str(_) for _ in to_write]))

    of1.close()
    return

def write_target_protein_fasta(options, peptides):

    if not options.fasta:
        print('No fasta file supplied - no proteins written')
        return

    decoy_counter = 0
    sequences = []
    with fasta.read(options.fasta, use_index = True) as db:
        for i, entry in enumerate(db):
            for p in peptides:
                if p.sequence in entry.sequence:
                    sequences.append([ entry.description, entry.sequence ])
            else:
                if decoy_counter < options.decoys:
                    sequences.append([ entry.description, entry.sequence ])
                    decoy_counter += 1

    print('Writing %s sequences including %s decoys to fasta' %(len(sequences), decoy_counter))
    fasta.write(
        sequences, file_mode = 'w',
        output = os.path.join(options.out_dir, '%s.fasta' %options.output_label)
    )
    return

def run_diann(input_dir):

    options.diann_path = '/usr/diann/1.8/diann-1.8'
    out_dir = os.path.join(input_dir,'out')

    try:
        os.makedirs(os.path.join(input_dir,'out'))
    except:
        pass

    files = os.listdir(input_dir)

    mzml_file = [_ for _ in files if '.mzml' in _.lower()]
    options.fasta = [_ for _ in files if '.fasta' in _.lower()]

    assert len(mzml_file) == 1
    assert len(options.fasta) == 1

    mzml_file = os.path.join(input_dir, mzml_file[0])
    options.fasta = os.path.join(input_dir, options.fasta[0])

    cmd = '%s ' %options.diann_path
    cmd += '--f %s ' % mzml_file
    cmd += '--fasta %s ' % options.fasta
    cmd += '--out %s ' % os.path.join(out_dir, 'out.tsv')
    cmd += '--out-lib %s ' % os.path.join(out_dir, 'out.lib.tsv')
    cmd += '--lib --predictor  --threads 8 --verbose 3 --qvalue 0.01 --matrices --gen-spec-lib --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --smart-profiling --pg-level 1 --prosit --predictor'

    os.system(cmd)
    return

def calculate_peak_parameters(options):

    # resolution = m / dm
    # dm = m / resolution

    # calculate peak FWHM
    dm_ms1 = options.resolution_at / float(options.ms1_resolution)
    dm_ms2 = options.resolution_at / float(options.ms2_resolution)

    # calculate peak standard deviations
    # sigma = FWHM / (2 * sqrt(2 * loge(2)))
    factor = 2.35482004503095
    options.ms1_stdev = dm_ms1 / factor
    options.ms2_stdev = dm_ms2 / factor

    # want at least n points in the peak region above the FWHM
    options.ms1_point_diff = dm_ms1 / options.n_points_gt_fwhm
    options.ms2_point_diff = dm_ms2 / options.n_points_gt_fwhm

    # chromatographic peak stdev
    options.rt_stdev = options.rt_peak_fwhm / factor

    return options


def simulate_isotope_patterns(*peptide_subset):
    for p in peptide_subset:
        p.get_ms1_isotope_pattern()
    return peptide_subset

def assemble(options):

    start = datetime.datetime.now()
    print('Started Synthedia %s' % start)
    print('Executing with %s processors' %options.num_processors)

    if not any([options.mq_txt_dir, options.prosit]):
        print('Either an MaxQuant output directory or Prosit file is required')
        print('Exiting')
        return

    if options.write_protein_fasta == True and options.fasta == None:
        print('Synthedia was asked to write a FASTA file but no input FASTA was given')
        print('Exiting')
        return

    try:
        os.makedirs(options.out_dir)
    except:
        pass

    print('Calculating peak parameters')
    options = calculate_peak_parameters(options)

    options.original_run_length = options.original_run_length * 60
    options.new_run_length = options.new_run_length * 60

    print('Writing outputs to %s' %options.out_dir)

    print('Preparing spectral template')
    spectra = make_spectra(options)

    print('Constructing peptide models')
    if options.use_existing_peptide_file:
        print('Using existing peptide file')

        if not os.path.isfile(options.use_existing_peptide_file):
            print('The specified peptide file was not found')
            print('Exiting')
            return

        with open( options.use_existing_peptide_file , 'rb') as handle:
            peptides = pickle.load(handle)

    else:
        if options.mq_txt_dir:
            peptides = read_peptides_from_mq(options)
        elif options.prosit:
            peptides = read_peptides_from_prosit(options)

        if options.num_processors == 1:
            for p in peptides:
                p.get_ms1_isotope_pattern()
        else:
            pool = multiprocessing.Pool(processes = options.num_processors)

            # split work into equal sized lists
            arg_sets = [[] for _ in range(options.num_processors)]
            counter = 0
            for p in peptides:
                arg_sets[counter].append(p)
                counter += 1
                if counter == options.num_processors:
                    counter = 0

            # send work to procs and collect results
            peptides = []
            print('Simulating isotope patterns')
            for _ in pool.starmap(simulate_isotope_patterns, arg_sets):
                peptides.extend(list(_))

            pool.close()
            pool.join()

        with open( os.path.join(options.out_dir, 'peptides.pickle') , 'wb') as handle:
            pickle.dump(peptides, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print('Generating protein abundance profiles')
    peptides = generate_protein_abundances_for_mzml_files(options, peptides)

    print('Scaling retention times')
    peptides = calculate_scaled_retention_times(options, peptides)

    if len(peptides) == 0:
        print('Error - no peptides to write')
        print('Exiting')
        return

    if options.num_processors == 1:

        for groupi in range(options.n_groups):
            for samplei in range(options.samples_per_group):
                print('Writing peptides to spectra')
                populate_spectra(options, peptides, spectra, groupi, samplei)
    else:
        print('Writing peptides to spectra')
        arg_sets = []
        for groupi in range(options.n_groups):
            for samplei in range(options.samples_per_group):
                arg_sets.append([ options, peptides, spectra, groupi, samplei ])

        pool = multiprocessing.Pool(processes = options.num_processors)
        pool.starmap(populate_spectra, arg_sets)
        pool.close()
        pool.join()

    print('Writing peptide target table')
    write_peptide_target_table(options, peptides, spectra)

    if options.write_protein_fasta:
        print('Writing protein fasta file')
        write_target_protein_fasta(options, peptides)

    if options.run_diann:
        print('Running DIA-NN')
        run_diann(OUT_DIR)

    if (options.all == True) or (options.tic == True):
        print('Plotting TIC')
        plotting.plot_tic(options)

    end = datetime.datetime.now()
    print('Done!')
    print('Total execution time: %s' %(end - start))
    return