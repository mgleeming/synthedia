import os, sys, time, copy, pickle, multiprocessing
import random, math, datetime
import pandas as pd
import numpy as np
from numba import jit

from . import plotting
from .peptide import SyntheticPeptide
from .mzml import Spectrum, MZMLWriter, MZMLReader

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
# 13) Add contaminant wall ions
# 14) Chemical noise
# 15) Missing value plots

@jit(nopython=True)
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

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

        peptides.append( SyntheticPeptide(prosit_entry = precursor))
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
            SyntheticPeptide(evidence_entry = evidence_row, msms_entry = msms_entry)
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
            sample_abundances = [
                np.random.normal(
                    loc = group_abundance, scale = options.within_group_stdev
                ) for _ in range(options.samples_per_group)
            ]
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
                spectrum.add_peaks(options, p.scaled_rt, p.ms1_isotopes, abundance_offset)

            elif spectrum.order == 2:
                if (p.mz > spectrum.isolation_ll) and (p.mz < spectrum.isolation_hl):
                    abundance_offset = p.offsets[groupi][samplei]
                    spectrum.add_peaks(options, p.scaled_rt, p.ms2_peaks, abundance_offset)

        # write final spec to file
        run.write_spec(options, spectrum)

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
