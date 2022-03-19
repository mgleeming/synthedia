import os, sys, time, copy, pickle, multiprocessing
import random, math, datetime, logging
import pandas as pd
import numpy as np

from . import plotting
from .peptide import SyntheticPeptide, calculate_scaled_retention_times
from .mzml import Spectrum, MZMLWriter, MZMLReader
from .peak_models import PeakModels

# TODO
# 1) Acquisition schema import - Done
# 2) Spectrum centroiding - Donej
# 3) Automated determination of peak standard deviations - Done
# 4) Randomness in elution profiles - Done
# 6) Decoys? - Done
# 7) Plots? - Done
# 8) Quantification between multiple files - Done
# 9) Peak tailing - Done
# 10) Probability sample in missing - Done
# 11) Probability group in missing - Done
# 12) Probability missing increases as abundance decreases
# 13) Add contaminant wall ions
# 14) Chemical noise
# 15) Missing value plots
# 16) Make sure different peptides of the same charge state get same tailing factors, abundance ratios, dropout probabilities - Done
# 17) Might need change RT peak simulation limits to a fixed value cutoff rather than a percentage of max intensity - Done
# 18) Handle case where multipe raw files are given in MQ input - Done

def make_spectra(options):

    logger = logging.getLogger("assembly_logger")
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
            logger.error('Error parsing acquisition schema file')
            logger.error('Exiting')
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

def generate_group_and_sample_abundances(options):

    peptide_abundance_offsets_between_groups = [0] + [
        np.random.normal(
            loc = 0, scale = options.between_group_stdev
        ) for _ in range(options.n_groups - 1)
    ]

    sample_abundance_offsets = [
        [
            np.random.normal(
                loc = group_abundance_offset, scale = options.within_group_stdev
            ) for _ in range(options.samples_per_group)
        ] for group_abundance_offset in peptide_abundance_offsets_between_groups
    ]

    return peptide_abundance_offsets_between_groups, sample_abundance_offsets

def generate_group_and_sample_probabilities(options):

    found_in_group = [1 if random.randint(0,100) >= options.prob_missing_in_group else 0 for _ in range(options.n_groups)]
    found_in_sample = []
    for group in found_in_group:

        if group == 1:
            found_in_sample.append([
                1 if random.randint(0,100) >= options.prob_missing_in_sample else 0 for _ in range(options.samples_per_group)
            ])
        else:
            found_in_sample.append([
                0 for _ in range(options.samples_per_group)
            ])

    return found_in_group, found_in_sample

def read_peptides_from_prosit(options):

    logger = logging.getLogger("assembly_logger")

    # read inputs
    prosit = pd.read_csv(options.prosit, sep = ',')

    # iRT values can be negative
    prosit['Retention time'] = prosit['iRT']
    min_rt = min(prosit['iRT'].tolist())
    if min_rt < 0:
        prosit['Retention time'] = prosit['Retention time'] + abs(min_rt)

    # add rt_buffer to start of run
    prosit['Retention time'] = prosit['Retention time'] + options.rt_buffer

    counter = 0
    peptides = []
    for _, peptide in prosit.groupby(['ModifiedPeptide']):

        # simulate abundances
        peptide_abundance_offsets_between_groups, sample_abundance_offsets = generate_group_and_sample_abundances(options)
        found_in_group, found_in_sample = generate_group_and_sample_probabilities(options)

        for __, precursor in peptide.groupby(['PrecursorCharge']):

            counter += 1
            if counter % 100 == 0:
                logger.info('Constructing peptide %s of %s' %(counter, len(prosit)))

            # check if precursor out of bounds
            if (float(precursor['PrecursorMz'].iloc[0]) < options.ms1_min_mz) or (float(precursor['PrecursorMz'].iloc[0]) > options.ms1_max_mz):
                continue

            peptides.append(
                SyntheticPeptide(
                    options,
                    prosit_entry = precursor,
                    peptide_abundance_offsets_between_groups = peptide_abundance_offsets_between_groups,
                    sample_abundance_offsets = sample_abundance_offsets,
                    found_in_group = found_in_group,
                    found_in_sample = found_in_sample,
                )
            )

    logger.info('Finished constructing %s peptides' %(len(peptides)))
    return peptides

def read_decoys_from_msp(options, peptides):

    logger = logging.getLogger("assembly_logger")

    peptide_intensities = [p.intensity for p in peptides]

    # read inputs
    lipids = []
    new_lipid = []
    with open(options.decoy_msp_file,'r') as if1:
        for l in if1:
            if l.strip() == '':
                lipids.append(new_lipid)
                new_lipid = []
            else:
                new_lipid.append(l.strip())

    lipids = [_ for _ in lipids if 'PRECURSORTYPE: [M+H]+' in _]

    counter = 0
    peptides = []

    rt_steps = np.arange(options.rt_buffer,options.new_run_length,0.001)

    for l in lipids:

        # simulate abundances
        peptide_abundance_offsets_between_groups, sample_abundance_offsets = generate_group_and_sample_abundances(options)
        found_in_group, found_in_sample = generate_group_and_sample_probabilities(options)

        fragments = False
        lipid_dict = {'fragments': []}
        for item in l:
            if 'Num Peaks:' in item:
                fragments = True
                continue
            if not fragments:
                try:
                    k,v = item.split(': ')

                    if k == 'RETENTIONTIME':
                        v = rt_steps[random.randint(0,len(rt_steps))]

                    lipid_dict[k] = v
                except:
                    continue
            if fragments:
                mz, intensity = item.split('\t')
                lipid_dict['fragments'].append([float(mz), float(intensity)])

        peptides.append(
            SyntheticPeptide(
                options,
                msp_entry = [lipid_dict, min(peptide_intensities), max(peptide_intensities)],
                peptide_abundance_offsets_between_groups = peptide_abundance_offsets_between_groups,
                sample_abundance_offsets = sample_abundance_offsets,
                found_in_group = found_in_group,
                found_in_sample = found_in_sample,
            )
        )

        if len(peptides) == options.num_decoys:
            break

    logger.info('\tFinished constructing %s decoys' %(len(peptides)))
    return peptides

def read_peptides_from_mq(options):

    logger = logging.getLogger("assembly_logger")

    # read inputs
    msms = pd.read_csv(os.path.join(options.mq_txt_dir, 'msms.txt'), sep = '\t')
    evidence = pd.read_csv(os.path.join(options.mq_txt_dir, 'evidence.txt'), sep = '\t')

    raw_files = list(set(evidence['Raw file'].tolist()))

    if len(raw_files) > 1:
        evidence = evidence[evidence['Raw file'] == raw_files[0]]

    # add rt_buffer to start of run
    evidence['Retention time'] = evidence['Retention time'] + options.rt_buffer

    evidence = evidence.sort_values(by = ['PEP'])

    # filter unwanted proteins
    for filterTerm in options.filterTerm:
        evidence = evidence.loc[-evidence['Proteins'].str.contains(filterTerm, na=False)]

    evidence = evidence[evidence['PEP'] < options.mq_pep_threshold]

    # restrict to unmodified peptides for the moment
    # NB - carbamidomethyl cys is retained
    evidence = evidence[evidence['Modifications'] == 'Unmodified']

    counter = 0
    peptides = []

    # group peptides from sample so all charge states of the same peptide can be given the same abundances
    for __, peptide in evidence.groupby(['Modified sequence']):

        # simulate abundances
        peptide_abundance_offsets_between_groups, sample_abundance_offsets = generate_group_and_sample_abundances(options)
        found_in_group, found_in_sample = generate_group_and_sample_probabilities(options)

        for ___, evidence_row in peptide.iterrows():

            counter += 1
            if counter % 100 == 0:
                logger.info('Constructing peptide %s of %s' %(counter, len(evidence)))

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
                SyntheticPeptide(
                    options,
                    evidence_entry = evidence_row,
                    msms_entry = msms_entry,
                    peptide_abundance_offsets_between_groups = peptide_abundance_offsets_between_groups,
                    sample_abundance_offsets = sample_abundance_offsets,
                    found_in_group = found_in_group,
                    found_in_sample = found_in_sample,
                )
            )

    logger.info('Finished constructing %s peptides' %(len(peptides)))
    return peptides

def populate_spectra(options, peptides, spectra, groupi, samplei):
    logger = logging.getLogger("assembly_logger")

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
            logger.info('\tGroup %s, Sample %s - Writing spectrum %s of %s' %(groupi, samplei, spectrumi, len(spectra)))

        # make spec numpy arrays on the fly to sav mem
        spectrum.make_spectrum(MS1_MZS, MS1_INTS, MS2_MZS, MS2_INTS)

        peptide_subset = [p for p in peptides if abs(p.scaled_rt - spectrum.rt) < 30]

        for p in peptide_subset:

            # skip missing peptides
            if (p.found_in_sample[groupi][samplei] == 0):
                continue

            if spectrum.order == 1:
                spectrum.add_peaks(options, p, groupi, samplei)
            elif spectrum.order == 2:
                if (p.mz > spectrum.isolation_ll) and (p.mz < spectrum.isolation_hl):
                    spectrum.add_peaks(options, p, groupi, samplei)

        # write final spec to file
        run.write_spec(options, spectrum)

        # this deletes m/z and intensity arrays that aren't needed anymore
        spectrum.clear()

    # close consumer
    run.close()

    return

def write_peptide_target_table(options, peptides, spectra):

    ms1_rts = np.asarray([s.rt for s in spectra if s.order == 1])

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
    # precursor abundances in file
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Intensity group_%s_sample_%s' %(group, sample))
    # precursor offsets
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Offset group_%s_sample_%s' %(group, sample))
    # precursor found in group
    for group in range(options.n_groups):
        to_write.append('Found in group_%s' %(group))
    # precursor found in sample
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Found in group_%s_sample_%s' %(group, sample))


    of1.write('%s\n' %'\t'.join([str(_) for _ in to_write]))

    print(ms1_rts)
    for p in peptides:
        rt_min, rt_max = p.calculate_retention_length(options, ms1_rts)
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

        # precursor abundances in file
        for group in p.abundances:
            for sample in group:
                to_write.append(sample)

        # precursor offsets
        for group in p.offsets:
            for sample in group:
                to_write.append(sample)

        # precursor found in group
        for group in range(options.n_groups):
            to_write.append(p.found_in_group[group])

        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append(p.found_in_sample[group][sample])

        of1.write('%s\n' %'\t'.join([str(_) for _ in to_write]))

    of1.close()
    return

def write_target_protein_fasta(options, peptides):

    logger = logging.getLogger("assembly_logger")
    if not options.fasta:
        logger.warning('No fasta file supplied - no proteins written')
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

    logger.info('Writing %s sequences including %s decoys to fasta' %(len(sequences), decoy_counter))
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

def get_rt_range_from_input_data(options):

    if options.mq_txt_dir:
        evidence = pd.read_csv(os.path.join(options.mq_txt_dir, 'evidence.txt'), sep = '\t')
        rts = evidence['Retention time'].tolist()

    elif options.prosit:
        prosit = pd.read_csv(options.prosit, sep = ',')
        rts = prosit['iRT'].tolist()

    # iRT values can be negative
    rt_range = abs(max(rts) - min(rts))

    return rt_range

def get_extra_parameters(options):

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

    # determine peak models to use
    pm = PeakModels()
    options.rt_peak_model = pm.get_rt_peak_model(options)
    options.mz_peak_model = pm.get_mz_peak_model(options)

    if int(options.original_run_length) == 0:
        # get default retention time range
        options.original_run_length = get_rt_range_from_input_data(options)

    if int(options.new_run_length) == 0:
        # fall back to original value
        options.new_run_length = options.original_run_length

    # add buffer so peaks are not simulated on the boundaries of the acquisition
    options.original_run_length += options.rt_buffer * 2

    return options

def simulate_isotope_patterns(*peptide_subset):
    for p in peptide_subset:
        p.get_ms1_isotope_pattern()
    return peptide_subset

def configure_logging(options):
    logger = logging.getLogger("assembly_logger")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    fh = logging.FileHandler(os.path.join(options.out_dir, 'assembly.log'))
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

def assemble(options):

    try:
        os.makedirs(options.out_dir)
    except:
        pass

    logger = configure_logging(options)

    start = datetime.datetime.now()
    logger.info('Started Synthedia %s' % start)
    logger.info('Executing with %s processors' %options.num_processors)

    logger.info('Config args:')
    for k,v in options.__dict__.items():
        logger.info('\t%s: %s' %(k,v))

    if not any([options.mq_txt_dir, options.prosit]):
        logger.error('Either an MaxQuant output directory or Prosit file is required')
        logger.error('Exiting')
        return

    if options.write_protein_fasta == True and options.fasta == None:
        logger.error('Synthedia was asked to write a FASTA file but no input FASTA was given')
        logger.error('Exiting')
        return

    logger.info('Calculating peak parameters')
    options = get_extra_parameters(options)

    options.original_run_length = options.original_run_length * 60
    options.new_run_length = options.new_run_length * 60

    logger.info('Writing outputs to %s' %options.out_dir)

    logger.info('Preparing spectral template')
    spectra = make_spectra(options)

    logger.info('Constructing peptide models')
    if options.use_existing_peptide_file:
        logger.info('Using existing peptide file')

        if not os.path.isfile(options.use_existing_peptide_file):
            logger.info('The specified peptide file was not found')
            logger.info('Exiting')
            return

        with open( options.use_existing_peptide_file , 'rb') as handle:
            peptides = pickle.load(handle)

    else:
        if options.mq_txt_dir:
            peptides = read_peptides_from_mq(options)
        elif options.prosit:
            peptides = read_peptides_from_prosit(options)

        if options.decoy_msp_file:
            logger.info('Reading decoy file')
            decoys = read_decoys_from_msp(options, peptides)
            peptides = peptides + decoys

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
            logger.info('Simulating isotope patterns')
            for _ in pool.starmap(simulate_isotope_patterns, arg_sets):
                peptides.extend(list(_))

            pool.close()
            pool.join()

        with open( os.path.join(options.out_dir, 'peptides.pickle') , 'wb') as handle:
            pickle.dump(peptides, handle, protocol=pickle.HIGHEST_PROTOCOL)


    # TODO
    # fix this
    # some problem with entities being modeled outside rt boundaries
    spec_rts = [s.rt for s in spectra]
    new_peptides = []
    for p in peptides:
        if p.scaled_rt > min(spec_rts) and p.scaled_rt < max(spec_rts):
            new_peptides.append(p)
    peptides = new_peptides

    if options.rescale_rt:
        logger.info('Scaling retention times')
        peptides = calculate_scaled_retention_times(options, peptides)

    if len(peptides) == 0:
        logger.error('Error - no peptides to write')
        logger.error('Exiting')
        return

    if options.num_processors == 1:

        for groupi in range(options.n_groups):
            for samplei in range(options.samples_per_group):
                logger.info('Writing peptides to spectra')
                populate_spectra(options, peptides, spectra, groupi, samplei)
    else:
        logger.info('Writing peptides to spectra')
        arg_sets = []
        for groupi in range(options.n_groups):
            for samplei in range(options.samples_per_group):
                arg_sets.append([ options, peptides, spectra, groupi, samplei ])

        pool = multiprocessing.Pool(processes = options.num_processors)
        pool.starmap(populate_spectra, arg_sets)
        pool.close()
        pool.join()

    logger.info('Writing peptide target table')
    write_peptide_target_table(options, peptides, spectra)

    if options.write_protein_fasta:
        logger.info('Writing protein fasta file')
        write_target_protein_fasta(options, peptides)

    if options.run_diann:
        logger.info('Running DIA-NN')
        run_diann(OUT_DIR)

    if (options.all == True) or (options.tic == True):
        logger.info('Plotting TIC')
        plotting.plot_tic(options)

    end = datetime.datetime.now()
    logger.info('Done!')
    logger.info('Total execution time: %s' %(end - start))
    return
