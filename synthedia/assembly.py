import os, sys, time, copy, pickle, yaml
import random, math, datetime, logging
import pandas as pd
import numpy as np

from . import plotting
from .peptide import SyntheticPeptide, calculate_scaled_retention_times, calculate_retention_lengths, calculate_feature_windows
from .mzml import Spectrum, MZMLWriter, MZMLReader
from .peak_models import PeakModels
from .io import InputReader, AcquisitionSchema

class NoPeptidesToSimulateError(Exception):
    pass

class AcquisitionSchemaError(Exception):
    pass

class IncorrectInputError(Exception):
    pass

def populate_spectra(options, peptides, spectra, groupi, samplei):
    logger = logging.getLogger("assembly_logger")

    MS1_MZS = np.arange(options.ms1_min_mz, options.ms1_max_mz, options.ms1_point_diff)
    MS1_INTS = np.zeros(len(MS1_MZS), dtype = int)

    MS2_MZS = np.arange(options.ms2_min_mz, options.ms2_max_mz, options.ms2_point_diff)
    MS2_INTS = np.zeros(len(MS2_MZS), dtype = int)

    run = MZMLWriter(
        os.path.join(
            options.out_dir, '%s_group_%s_sample_%s.mzML' %(
                options.output_label, groupi, samplei
            )
        ),
        len(spectra)
    )

    min_rt = min([p.min_scaled_peak_rt_list[groupi][samplei] for p in peptides])
    max_rt = max([p.max_scaled_peak_rt_list[groupi][samplei] for p in peptides])

    return_peptides = []

    peptides.sort(key=lambda p: p.scaled_rt_lists[groupi][samplei], reverse = True)

    peptide_subset = []
    for spectrumi, spectrum in enumerate(spectra):

        if (spectrum.rt < min_rt) or (spectrum.rt > max_rt):
            continue

        if spectrumi % 1000 == 0:
            logger.info('\tGroup %s, Sample %s - Writing spectrum %s of %s' %(
                groupi, samplei, spectrumi, len(spectra))
            )

        # make spec numpy arrays on the fly to save memory
        spectrum.make_spectrum(MS1_MZS, MS1_INTS, MS2_MZS, MS2_INTS)

        while True:
            if len(peptides) == 0: break
            if spectrum.rt < peptides[-1].min_scaled_peak_rt_list[groupi][samplei]: break
            if spectrum.rt >= peptides[-1].min_scaled_peak_rt_list[groupi][samplei]:
                peptide_subset.append(peptides.pop())

        while True:
            if len(peptide_subset) == 0: break
            if spectrum.rt <= peptide_subset[0].max_scaled_peak_rt_list[groupi][samplei]:
                break
            if spectrum.rt > peptide_subset[0].max_scaled_peak_rt_list[groupi][samplei]:
                return_peptides.append(peptide_subset.pop(0))

        for p in peptide_subset:

            if spectrum.rt > p.max_scaled_peak_rt_list[groupi][samplei]:
                continue

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

    while len(peptide_subset) != 0:
        return_peptides.append(peptide_subset.pop())

    # close consumer
    run.close()

    return return_peptides

def write_peptide_target_table(options, peptides):

    of1 = open( os.path.join(options.out_dir, '%s_peptide_table.tsv' %options.output_label),'wt')

    to_write = [
        'Protein',
        'Sequence',
        'Intensity',
        'm/z',
        'Charge',
        'Mass',
        'Experimental RT',
        'Synthetic theoretical m/z 0',
        'Synthetic max fragment m/z'
    ]

    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Synthetic RT group_%s_sample_%s' %(group, sample))
            to_write.append('Synthetic RT start group_%s_sample_%s' %(group, sample))
            to_write.append('Synthetic RT end group_%s_sample_%s' %(group, sample))

    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Synthetic m/z group_%s_sample_%s' %(group, sample))

    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Synthetic m/z ppm error group_%s_sample_%s' %(group, sample))

     # precursor measured abundances in file
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Total precursor abundance group_%s_sample_%s' %(group, sample))

    # precursor measured height
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Precursor max peak height group_%s_sample_%s' %(group, sample))

    # total abundance of most abundant fragment ion
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Most abundant fragment total intensity group_%s_sample_%s' %(group, sample))

    # max height of most abundant fragment ion
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('Most abundant fragment max height group_%s_sample_%s' %(group, sample))

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
    # MS1 points per peak
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('MS1 chromatographic points group_%s_sample_%s' %(group, sample))
    # MS2 points per peak
    for group in range(options.n_groups):
        for sample in range(options.samples_per_group):
            to_write.append('MS2 chromatographic points group_%s_sample_%s' %(group, sample))

    of1.write('%s\n' %'\t'.join([str(_) for _ in to_write]))

    for p in peptides:
        to_write = [
            p.protein,
            p.sequence,
            p.intensity,
            p.mz,
            p.charge,
            p.mass,
            p.rt,
            '%.5f' %p.ms1_isotopes[0].base_mz, # exact mz i.e. no error
            '%.5f' %p.max_fragment_mz
        ]

        # precursor simulated retention times
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append('%.2f'%p.scaled_rt_lists[group][sample])
                to_write.append('%.2f'%p.get_min_peak_rt(group, sample))
                to_write.append('%.2f'%p.get_max_peak_rt(group, sample))

        # true m/z values written to mzml
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append('%.5f'%p.ms1_isotopes[0].peak_mz_list[group][sample])

        # true m/z values written to mzml
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append('%.3f'%p.ms1_isotopes[0].peak_mz_error_list[group][sample])

        # precursor measured abundances in file
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append(p.get_total_precursor_intensity(group, sample))

        # precursor measured height
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append(p.get_max_precursor_intensity(group, sample))

        # max intensity fragment total abundance
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append(p.get_total_fragment_intensity(group, sample))

        # max intensity fragment max height
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append(p.get_max_fragment_intensity(group, sample))

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

        # points per MS1 peak
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append(p.get_points_per_peak(group, sample, 1))

        # points per MS2 peak
        for group in range(options.n_groups):
            for sample in range(options.samples_per_group):
                to_write.append(p.get_points_per_peak(group, sample, 2))

        of1.write('%s\n' %'\t'.join([str(_) for _ in to_write]))

    of1.close()
    return

def get_original_run_length(rts):
    if min(rts) >= 0:
        return max(rts)
    if min(rts) < 0:
        return max(rts) + abs(min(rts))

def get_rt_range_from_input_data(options):

    logger = logging.getLogger("assembly_logger")

    if options.mq_txt_dir:
        if not os.path.isfile(os.path.join(options.mq_txt_dir, 'msms.txt')):
            msg = 'The specified MaxQuant msms.txt file does not exist'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)
        if not os.path.isfile(os.path.join(options.mq_txt_dir, 'evidence.txt')):
            msg = 'The specified MaxQuant evidence.txt file does not exist'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)
        evidence = pd.read_csv(os.path.join(options.mq_txt_dir, 'evidence.txt'), sep = '\t')
        rts = evidence['Retention time'].tolist()
    elif options.prosit:
        if not os.path.isfile(os.path.join(options.prosit)):
            msg = 'The specified Prosit library file does not exist'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)
        prosit = pd.read_csv(options.prosit, sep = ',')
        rts = prosit['iRT'].tolist()
    elif options.use_existing_peptide_file:
        if not os.path.isfile(options.use_existing_peptide_file):
            msg = 'The specified peptide file was not found'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)

        with open( options.use_existing_peptide_file , 'rb') as handle:
            peptides = pickle.load(handle)
        rts = [p.rt / 60 for p in peptides]

    return get_original_run_length(rts)

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

    return options

def simulate_isotope_patterns(*peptide_subset):
    for p in peptide_subset:
        p.get_ms1_isotope_pattern()
    return peptide_subset

def configure_logging(options):
    logger = logging.getLogger("assembly_logger")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    if not options.silent:
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

    # save init args
    with open(os.path.join(options.out_dir, '%s_simulation_args.yaml'%options.output_label), 'w') as f:
        yaml.dump(options.__dict__, f)

    logger = configure_logging(options)

    start = datetime.datetime.now()
    logger.info('Started Synthedia %s' % start)

    logger.info('Config args:')
    for k,v in options.__dict__.items():
        logger.info('\t%s: %s' %(k,v))

    if not any([options.mq_txt_dir, options.prosit, options.use_existing_peptide_file]):
        msg = 'Either an MaxQuant output directory, Prosit library, or peptide file from a previous simulation is required'
        logger.error(msg)
        logger.error('Exiting')
        raise IncorrectInputError(msg)

    logger.info('Calculating peak parameters')
    options = get_extra_parameters(options)

    options.original_run_length = options.original_run_length * 60
    options.new_run_length = options.new_run_length * 60

    logger.info('Writing outputs to %s' %options.out_dir)

    logger.info('Preparing spectral template')
    acquisition_schema = AcquisitionSchema(options)
    spectra = acquisition_schema.get_spectra()

    logger.info('Constructing peptide models')
    input_peptides = InputReader(options)
    peptides = input_peptides.get_peptides()

    if options.use_existing_peptide_file:
       logger.info('Scaling retention times')
       peptides = calculate_scaled_retention_times(options, peptides)

    logger.info('Calculating feature windows')
    calculate_feature_windows(options, peptides, spectra)

    logger.info('Calculating retention lengths')
    peptides = calculate_retention_lengths(options, peptides, spectra)

    for groupi in range(options.n_groups):
        for samplei in range(options.samples_per_group):
            logger.info('Writing peptides to spectra')
            peptides = populate_spectra(options, peptides, spectra, groupi, samplei)

    logger.info('Writing peptide target table')
    write_peptide_target_table(options, peptides)

    if (options.all == True) or (options.tic == True):
        logger.info('Plotting TIC')
        plotting.plot_tic(options)

    end = datetime.datetime.now()
    logger.info('Done!')
    logger.info('Total execution time: %s' %(end - start))
    return


