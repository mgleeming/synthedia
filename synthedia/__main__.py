import os, sys, argparse, yaml
import multiprocessing
from synthedia import assembly

def main(args = None):

    parser = argparse.ArgumentParser(
        description = 'Generate synthetic DIA LC-MS/MS bottom up proteomics data with known composition.'
    )

    io_args = parser.add_argument_group("Input/Output")
    io_args.add_argument( '--mq_txt_dir', required = False, type = str,
                        help = 'Path to MaxQuat "txt" directory.')
    io_args.add_argument( '--prosit', required = False, type = str,
                        help = 'Path to prosit prediction library.')
    io_args.add_argument( '--prosit_peptide_abundance_mean', required = False, type = float, default = 22,
                        help = 'Mean log2 abundance used to simulate peptide abundances for prosit input types. Not used for MaxQuat input types.')
    io_args.add_argument( '--prosit_peptide_abundance_stdev', required = False, type = float, default = 3,
                        help = 'Standard deviation of gaussian used to simulate peptide abundances for prosit input types. Not used for MaxQuant input types.')
    io_args.add_argument( '--acquisition_schema', required = False, type = str,
                        help = 'Path to file defining MS2 acquisition schema.')
    io_args.add_argument( '--use_existing_peptide_file', required = False, type = str,
                        help = 'Path to an existin peptide file which will be used.')
    io_args.add_argument( '--out_dir', required = False, type = str, default = os.path.join(os.getcwd(), 'output'),
                        help = 'Output directory where results should be written.')
    io_args.add_argument( '--output_label', required = False, type = str, default = 'output',
                        help = 'Prefix for output files.')
    io_args.add_argument( '--config', required = False, type = str, default = None,
                        help = 'Path to *.yaml config file.')
    io_args.add_argument( '--silent', action= 'store_true',
                        help = 'Do not print logging output to terminal')

    filtering_args = parser.add_argument_group("Filtering")
    filtering_args.add_argument( '--mq_pep_threshold', required = False, type = float, default = 0.001,
                        help = 'For MaxQuant input data, use only peptides with a Posterior Error Probability (PEP) less than this value')
    filtering_args.add_argument('--filterTerm', required = False, type = str, action = 'append', default = ['CON_', 'REV'],
                        help = 'Terms used to filter input maxquant lists to remove unwanted targets. For example contaminant protein can be removed by specifying "--filterTerm CON_". Multiple filters can be applied. For example "--filterTerm CON_ --filterTerm REV_". Filters are used only for MaxQuant input types (no effect for Prosit) and are applied to the "Proteins" column of the evidence.txt table.')

    processing_args = parser.add_argument_group("Processing")
    processing_args.add_argument( '--num_processors', required = False, type = int, default = multiprocessing.cpu_count() ,
                        help = 'Number of cores to use in constructing mzML files. Defaults to all available cores')

    instrument_args = parser.add_argument_group("Instrument Parameters")
    instrument_args.add_argument( '--ms1_min_mz', required = False, type = int, default = 350,
                        help = 'Minimum m/z at MS1 level.')
    instrument_args.add_argument( '--ms1_max_mz', required = False, type = int, default = 1600,
                        help = 'Maximum m/z at MS1 level.')
    instrument_args.add_argument( '--ms2_min_mz', required = False, type = int, default = 100,
                        help = 'Minimum m/z at MS2 level.')
    instrument_args.add_argument( '--ms2_max_mz', required = False, type = int, default = 2000,
                        help = 'Maximum m/z at MS2 level.')
    instrument_args.add_argument( '--ms1_resolution', required = False, type = int, default = 120000,
                        help = 'Mass spectral resolution at MS1 level.')
    instrument_args.add_argument( '--ms2_resolution', required = False, type = int, default = 15000,
                        help = 'Mass spectral resolution at MS2 level.')
    instrument_args.add_argument( '--ms1_scan_duration', required = False, type = float, default = 0.37,
                        help = 'Time in seconds taken to record an MS1 scan.')
    instrument_args.add_argument( '--ms2_scan_duration', required = False, type = float, default = 0.037,
                        help = 'Time in seconds taken to record an MS2 scan.')
    instrument_args.add_argument( '--isolation_window', required = False, type = int, default = 30,
                        help = 'Length of DIA window in m/z.')
    instrument_args.add_argument( '--resolution_at', required = False, type = int, default = 200,
                        help = 'm/z value at which resolution is defined.')
    instrument_args.add_argument( '--n_points_gt_fwhm', required = False, type = int, default = 3,
                        help = 'Number of MS data points greater than the peak FWHM. Increasing this number means each mass spectral peak will be described by more data points but will also slow processing time and increase file size.')
    instrument_args.add_argument( '--esi_instability', required = False, type = float, default = 20,
                        help = 'Simulates imperfection in chromatographic peaks by applying a randomly intensity scaling factor to adjacent scans. A value of 0 indicates no randomness. A value of 100 indicates high spray instability.')
    instrument_args.add_argument( '--ms1_ppm_error_mean', required = False, type = float, default = 0,
                        help = 'The mean value of a Gaussian distribution from which PPM errors for MS1 precursors will be drawn. This value can be negative. Setting both ms1_ppm_error_mean, and ms1_ppm_error_stdev to 0 equates to perfect mass accuracy.')
    instrument_args.add_argument( '--ms1_ppm_error_stdev', required = False, type = float, default = 0,
                        help = 'The standard deviation of a Gaussian distribution from which PPM errors for MS1 precursors will be drawn. Setting both ms1_ppm_error_mean and ms1_ppm_error_stdev to 0 equates to perfect mass accuracy.')
    instrument_args.add_argument( '--ms2_ppm_error_mean', required = False, type = float, default = 0,
                        help = 'The mean value of a Gaussian distribution from which PPM errors for MS2 fragments will be drawn. This value can be negative. Setting both ms1_ppm_error_mean, and ms1_ppm_error_stdev to 0 equates to perfect mass accuracy.')
    instrument_args.add_argument( '--ms2_ppm_error_stdev', required = False, type = float, default = 0,
                        help = 'The standard deviation of a Gaussian distribution from which PPM errors for MS2 fragments will be drawn. Setting both ms1_ppm_error_mean and ms1_ppm_error_stdev to 0 equates to perfect mass accuracy.')

    chromatography_args = parser.add_argument_group("Chromatography")
    chromatography_args.add_argument( '--rt_peak_fwhm', required = False, type = float, default = 4,
                        help = 'Chromatographic peak full with at half maximum intehsity in seconds.')
    chromatography_args.add_argument( '--original_run_length', required = False, type = float, default = 0,
                        help = 'Length in minutes of original data file. If not given, this will be determined by taking the difference between the minimum and maximum peptide retention times. If set to "0", the retention time range will be automatically detected from the input data.')
    chromatography_args.add_argument( '--new_run_length', required = False, type = float, default = 0,
                        help = 'Length in minutes of new data file. If set to "0", the retention time range of the input data will be used.')
    chromatography_args.add_argument( '--rt_buffer', required = False, type = float, default = 5,
                        help = 'Time (in minutes) that should be appended to the beginning and end of the retention time range of a set of input peptides. This helps ensure that peptides at the boundaries of the elution window are simulated completely')
    chromatography_args.add_argument( '--rt_instability', required = False, type = int, default = 15,
                        help = 'Introduces an instability in retentention time values for the same peptide when a multi-group or multi-sample simulation is conducted. The value is the maximum number of seconds by which peptide retention times will differ')

    simulation_args = parser.add_argument_group("Simulation")
    simulation_args.add_argument( '--ms1_min_peak_intensity', required = False, type = float, default = 100,
                        help = 'Peptide elution profiles are simulated as gaussian peaks. This value sets the minimum gaussian curve intensitiy for a peptide to be simulated in MS1 spectra.')
    simulation_args.add_argument( '--ms2_min_peak_intensity', required = False, type = float, default = 10,
                        help = 'Peptide elution profiles are simulated as gaussian peaks. This value sets the minimum gaussian curve intensitiy for a peptide to be simulated in MS2 spectra.')
    simulation_args.add_argument( '--centroid_ms1', action = 'store_true',
                        help = 'If given, simulated MS1 mass spectra will be centroided. Otherwise, profile data will be written.')
    simulation_args.add_argument( '--centroid_ms2', action = 'store_true',
                        help = 'If given, simulated MS2 mass spectra will be centroided. Otherwise, profile data will be written.')
    simulation_args.add_argument( '--write_empty_spectra', action = 'store_true',
                        help = 'Write empty mass sepctra to the output data file')
    simulation_args.add_argument( '--mz_peak_model', required = False, type = str, default = 'gaussian',
                        help = 'The model used to simulate mass spectral peaks. Can be "gaussian", "exponentially_modified_gaussian" or "cauchy"')
    simulation_args.add_argument( '--rt_peak_model', required = False, type = str, default = 'exponentially_modified_gaussian',
                        help = 'The model used to simulate chromatographic peaks. Can be "gaussian", "exponentially_modified_gaussian" or "cauchy"')
    simulation_args.add_argument( '--mz_emg_k', required = False, type = float, default = 2,
                        help = 'Shape factor for exponentially modified gaussian in the mass spectral domain. Must be greater than 0. Increasing K results in more heavily tailed mass spectral peaks. This parameter is inactive unless --mz_peak_model is not set to exponentially_modified_gaussian.')
    simulation_args.add_argument( '--rt_emg_k', required = False, type = float, default = 1,
                        help = 'Shape factor for exponentially modified gaussian in the retention time domain. Must be greater than 0. Increasing K results in more heavily tailed chromatographic peaks. This parameter is inactive unless --rt_peak_model is not set to exponentially_modified_gaussian.')
    simulation_args.add_argument( '--prob_missing_in_sample', required = False, type = float, default = 0,
                        help = 'Probability (0-100) that a peptide is missing in any given sample')
    simulation_args.add_argument( '--prob_missing_in_group', required = False, type = float, default = 0,
                        help = 'Probability (0-100) that a peptide is missing in an entire group')

    plotting_args = parser.add_argument_group("Plotting")
    plotting_args.add_argument( '--tic', action = 'store_true',
                        help = 'Plot TIC for the generated mzML file.')
    plotting_args.add_argument( '--schema', action = 'store_true',
                        help = 'Plot acquisition schema.')
    plotting_args.add_argument( '--all', action = 'store_true',
                        help = 'Plot all graphics.')

    grouping_and_quant_args = parser.add_argument_group("Grouping and quantitation")
    grouping_and_quant_args.add_argument( '--n_groups', required = False, type = int, default = 1,
                        help = 'Number of treatment groups to simulate.')
    grouping_and_quant_args.add_argument( '--samples_per_group', required = False, type = int, default = 1,
                        help = 'Number of individual samples to simulate per treatment group.')
    grouping_and_quant_args.add_argument( '--between_group_stdev', required = False, type = float, default = 1.0,
                        help = 'Standard deviation of a normal distribution from which group means will be drawn.')
    grouping_and_quant_args.add_argument( '--within_group_stdev', required = False, type = float, default = 0.2,
                        help = 'Standard deviation of a normal distribution from which within group samples will be drawn.')

    decoy_args = parser.add_argument_group("Decoys")
    decoy_args.add_argument( '--decoy_msp_file', required = False, type = str,
                        help = 'Path to MSP file. Note - must include retention times.')
    decoy_args.add_argument( '--num_decoys', required = False, type = int, default = 500,
                        help = 'Number of decoy peaks to simulate')
    decoy_args.add_argument( '--simulate_top_n_decoy_fragments', required = False, type = int, default = 15,
                        help = 'Simulate n most intense fragments of the decoy compound.')
    decoy_args.add_argument( '--decoy_abundance_mean', required = False, type = float, default = 22,
                        help = 'Mean log2 abundance used to simulate decoy ion abundances.')
    decoy_args.add_argument( '--decoy_abundance_stdev', required = False, type = float, default = 3,
                        help = 'Standard deviation of gaussian used to simulate decoy ion abundances.')


    options =  parser.parse_args()

    try:
        os.makedirs(options.out_dir)
    except:
        pass

    # want heirarchy to be
    # user specified cli args - highest
    # config file args
    # argparse defaults - lowest

    # want to be able to specify args in config file but also overwrite these with command line inputs
    if options.config: # read config file if given

        # argparse defaults plus user cli args
        cmd_line_args = options.__dict__

        # load config file
        config_file_args = yaml.load(open(cmd_line_args['config']), Loader=yaml.FullLoader)

        # create new options dict - config_file takes priority
        options = dict(list(cmd_line_args.items()) + list(config_file_args.items()))

        # need to overwrite config file args if specified on the cmd line
        for arg in sys.argv:
            if arg[0:2] == '--':
                options[arg.replace('--','')] = cmd_line_args[arg.replace('--','')]

        options = Options(options)

    assembly.assemble(options)

class Options():
    def __init__(self, d):
        for k,v in d.items(): setattr(self,k,v)

if __name__ == '__main__':
    sys.exit(main())
