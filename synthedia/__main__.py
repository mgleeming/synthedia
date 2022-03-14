import os, sys, argparse, yaml
import multiprocessing
from synthedia import assembly

def main(args = None):

    print(sys.argv)
    parser = argparse.ArgumentParser(
        description = 'Generate DIA data from DDA MaxQuant output.'
    )

    io_args = parser.add_argument_group("Input/Output")
    io_args.add_argument( '--mq_txt_dir', required = False, type = str,
                        help = 'Path to MaxQuat "txt" directory.')
    io_args.add_argument( '--prosit', required = False, type = str,
                        help = 'Path to prosit prediction library.')
    io_args.add_argument( '--acquisition_schema', required = False, type = str,
                        help = 'Path to file defining MS2 acquisition schema.')
    io_args.add_argument( '--use_existing_peptide_file', required = False, type = str,
                        help = 'Path to an existin peptide file which will be used.')
    io_args.add_argument( '--write_protein_fasta', action = 'store_true',
                        help = 'Write FASTA file with protein sequences for simulated peptides. If given, a FASTA file must be supplied with the --fasta options.')
    io_args.add_argument( '--fasta', required = False, type = str,
                        help = 'Path to FASTA file from which protein sequences should be taken.')
    io_args.add_argument( '--out_dir', required = False, type = str, default = os.path.join(os.getcwd(), 'output'),
                        help = 'Output directory where results should be written.')
    io_args.add_argument( '--output_label', required = False, type = str, default = 'output',
                        help = 'Prefix for output files.')
    io_args.add_argument( '--config', required = False, type = str, default = None,
                        help = 'Path to *.yaml config file.')


    processing_args = parser.add_argument_group("Processing")
    processing_args.add_argument( '--num_processors', required = False, type = int, default = multiprocessing.cpu_count() ,
                        help = 'Number of cores to use in constructing mzML files. Defaults to all available cores')


    instrument_args = parser.add_argument_group("Instrument Parameters")
    instrument_args.add_argument( '--ms1_min_mz', required = False, type = float, default = 350,
                        help = 'Minimum m/z at MS1 level.')
    instrument_args.add_argument( '--ms1_max_mz', required = False, type = float, default = 1600,
                        help = 'Maximum m/z at MS1 level.')
    instrument_args.add_argument( '--ms2_min_mz', required = False, type = float, default = 100,
                        help = 'Minimum m/z at MS2 level.')
    instrument_args.add_argument( '--ms2_max_mz', required = False, type = float, default = 2000,
                        help = 'Maximum m/z at MS2 level.')
    instrument_args.add_argument( '--ms1_resolution', required = False, type = float, default = 120000,
                        help = 'Mass spectral resolution at MS1 level.')
    instrument_args.add_argument( '--ms2_resolution', required = False, type = float, default = 15000,
                        help = 'Mass spectral resolution at MS2 level.')
    instrument_args.add_argument( '--ms1_scan_duration', required = False, type = float, default = 0.37,
                        help = 'Time in seconds taken to record an MS1 scan.')
    instrument_args.add_argument( '--ms2_scan_duration', required = False, type = float, default = 0.037,
                        help = 'Time in seconds taken to record an MS2 scan.')
    instrument_args.add_argument( '--isolation_window', required = False, type = int, default = 30,
                        help = 'Length of DIA window in m/z.')
    instrument_args.add_argument( '--resolution_at', required = False, type = float, default = 200,
                        help = 'm/z value at which resolution is defined.')
    instrument_args.add_argument( '--n_points_gt_fwhm', required = False, type = int, default = 3,
                        help = 'Number of MS data points greater than the peak FWHM. Increasing this number means each mass spectral peak will be described by more data points but will also slow processing time and increase file size.')


    chromatography_args = parser.add_argument_group("Chromatography")
    chromatography_args.add_argument( '--rt_peak_fwhm', required = False, type = float, default = 7,
                        help = 'Chromatographic peak full with at half maximum intehsity in seconds.')
    chromatography_args.add_argument( '--rt_instability', required = False, type = float, default = 20,
                        help = 'Simulates imperfection in chromatographic peaks by applying a randomly intensity scaling factor to adjacent scans. A value of 0 indicates no randomness. A value of 100 indicates high spray instability.')
    chromatography_args.add_argument( '--original_run_length', required = False, type = float, default = 120,
                        help = 'Length in minutes of original data file. If not given, this will be determined by taking the difference between the minimum and maximum peptide retention times.')
    chromatography_args.add_argument( '--new_run_length', required = False, type = float, default = 12,
                        help = 'Length in minutes of new data file.')


    simulation_args = parser.add_argument_group('Simulation')
    simulation_args.add_argument( '--ms_clip_window', required = False, type = float, default = 0.15,
                        help = 'm/z window surrounding an MS peak that should be considered when simulating peak intensities. For high resolution data, this normally does not need to be changed.')
    simulation_args.add_argument( '--min_peak_fraction', required = False, type = float, default = 0.01,
                        help = 'Peptide elution profiles are simulated as gaussian peaks. This value sets the minimum gaussian curve intensitiy for a peptide to be simulated.')
    simulation_args.add_argument( '--centroid', action = 'store_true',
                        help = 'If given, simulated mass spectra will be centroided. Otherwise, profile data will be written.')
    simulation_args.add_argument( '--write_empty_spectra', action = 'store_true',
                        help = 'Write empty mass sepctra to the output data file')


    filtering_args = parser.add_argument_group('Filtering')
    filtering_args.add_argument( '--mq_pep_threshold', required = False, type = float, default = 0.001,
                        help = 'For MaxQuant input data, use only peptides with a Posterior Error Probability (PEP) less than this value')


    diann_args = parser.add_argument_group('Analysis')
    diann_args.add_argument( '--run_diann', action = 'store_true',
                        help = 'Run DIA-NN on the output data file.')
    diann_args.add_argument( '--diann_path', required = False, type = str, default = '/usr/diann/1.8/diann-1.8',
                        help = 'Path to DIA-NN.')


    plotting_args = parser.add_argument_group('Plotting')
    plotting_args.add_argument( '--tic', action = 'store_true',
                        help = 'Plot TIC for the generated mzML file.')
    plotting_args.add_argument( '--schema', action = 'store_true',
                        help = 'Plot acquisition schema')
    plotting_args.add_argument( '--all', action = 'store_true',
                        help = 'Plot all graphics')

    grouping_and_quant_args = parser.add_argument_group('Grouping and quantitation')
    grouping_and_quant_args.add_argument( '--n_groups', required = False, type = int, default = 1,
                        help = 'Number of treatment groups to simulate.')
    grouping_and_quant_args.add_argument( '--samples_per_group', required = False, type = int, default = 1,
                        help = 'Number of individual samples to simulate per treatment group.')
    grouping_and_quant_args.add_argument( '--between_group_stdev', required = False, type = float, default = 1.0,
                        help = 'Standard deviation of a normal distribution from which group means will be drawn.')
    grouping_and_quant_args.add_argument( '--within_group_stdev', required = False, type = float, default = 0.2,
                        help = 'Standard deviation of a normal distribution from which within group samples will be drawn.')

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

    # save init args
    with open(os.path.join(options.out_dir, '%s_simulation_args.yaml'%options.output_label), 'w') as f:
        yaml.dump(options.__dict__, f)

    assembly.assemble(options)

class Options():
    def __init__(self, d):
        for k,v in d.items(): setattr(self,k,v)

if __name__ == '__main__':
    sys.exit(main())
