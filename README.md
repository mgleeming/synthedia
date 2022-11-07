# Synthedia

[![Synthedia Integration](https://github.com/mgleeming/synthedia/actions/workflows/integration.yaml/badge.svg?branch=main)](https://github.com/mgleeming/synthedia/actions/workflows/integration.yaml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![Python](https://img.shields.io/badge/python-3.8-yellow.svg)
![Python](https://img.shields.io/badge/python-3.9-yellow.svg)

# Table of contents
**[1. Introduction](#introduction)**<br>
**[2. Getting started](#getting-started)**<br>
**[3. Installation](#installation)**<br>
**[4. Input data types](#input-data-types)**<br>
**[5. DIA acquisition strategies](#dia-acquisition-strategies)**<br>
**[6. Decoy signals](#decoy-signals)**<br>
**[7. Compatibility of Synthedia mzML files with other software](#compatibility-of-synthedia-mzml-files-with-other-softeare)**<br>
**[8. Notes about the resulting data](#notes-about-the-resulting-data)**<br>
**[9. Usage](#usage)**<br>
**[10. Synthedia parameter reference](#synthedia-parameter-reference)**<br>
**[11. Viewing mzML files](#viewing-mzml-files)**<br>
**[12. Support](#support)**<br>

## Introduction
A wide range of software is available to identify and quantify peptides in bottom-up proteomics data acquired by LC-MS/MS and, frequently, analysis of the same input data file with two different packages yields different results. These differences likely originate from differences in algorithms for preprocessing data, matching MSn spectra to peptides, quantifying peak areas, matching assignments between samples as well as differences in the myriad parameters that must be set by the user to initiate an analysis. Thorough analysis of these variables is, however, complicated since the ‘true’ sample composition of real data files is almost never completely known.

Synthedia creates synthetic LC-MS/MS data that mimics real data but with a composition that is exactly known. Currently, synthedia support the creation of Data-Independent Acquisition (DIA) style data wherein fixed, large m/z windows are sequentially isolated for fragmentation. We have focused on creating DIA data to date since the complexity of the analysis preformed by processing tools is substantial and the impact of different acquisition methodologies on the eventual outcome is somewhat more difficult to predict.

Synthedia can be configured to produce synthetic DIA data that mimics a wide range of acquisition strategies. For example, data can be simulated that models the same set of peptides eluted over progressively shorter gradients which would allow an assessment of how processing software copes with increasing complexity of data. Data can be modelled with different mass spectral resolutions, chromatographic peak tailing, DIA window strategies among others. Synthedia can also help to assess quantitation by creating multiple sets of data files which simulate replicates in a multi-group comparison. These modules include abilities to set within and between-group variability as well as simulate missing data in both a random and group-wise manner.

## Getting started

Synthedia can be used by downloading and installing the package or through our web server which can be accessed [here](https://synthedia.org/). Please note that producing synthetic data is a time-consuming task and capacity on our server is limited. In most cases, it will be faster to install synthedia on your own computer as this can take advantage of multiple processing cores.

## Installation

### Install via pip

Create a virtual environment (optional):
```
virtualenv venv
source venv/bin/activate
pip install --upgrade pip
```

Note that the [pyOpenMS](https://pyopenms.readthedocs.io/en/latest/index.html) dependency of synthedia requires python 3.7, 3.8 or 3.9. If you have multiple python versions installed, you can construct virtual environments with different versions by:
```
virtualenv venv --python=/path/to/python3.X
```
On linux systems, the path is usually ```/usr/bin/pythonX.X```

Install synthedia using ```pip```
```
pip install synthedia
```


### Install from source
Clone the repo:
```
git clone https://github.com/mgleeming/synthedia.git
cd synthedia
```
Create a virtual environment (optional):
```
virtualenv venv
source venv/bin/activate
pip install --upgrade pip
```

Note that the [pyOpenMS](https://pyopenms.readthedocs.io/en/latest/index.html) dependency of synthedia requires python 3.7, 3.8 or 3.9. If you have multiple python versions installed, you can construct virtual environments with different versions by:
```
virtualenv venv --python=/path/to/python3.X
```
On linux systems, the path is usually ```/usr/bin/pythonX.X```

Install synthedia:
```
pip install .
```

Synthedia can then be invoked from the command line by:
```
synthedia [COMMAND LINE ARGUMENTS]
```

Optionally, you can run the Synthedia integration test to check that everything is working:
```
pip install pytest
pytest -v
```

## Input data types

To create synthetic data, synthedia requires information about peptide fragmentation patterns and some metric as to relative retention time. These can be supplied by either a:
  - [Prosit](https://www.proteomicsdb.org/prosit/)-predicted spectral library
  - MaxQuant 'txt' directory

#### Prosit libraries

[Prosit](https://www.nature.com/articles/s41592-019-0426-7) is a machine learning-based application that predicts MS/MS fragmentation patterns for input peptide sequences. Input sequences can be arbitrary and do not necessarily need to originate from any specific organism or protein. Prosit will generate an output file that contains predicted abundances for peptide sequence ions.

To create a Synthedia-compatible Prosit spectral library:
  1. Go to https://www.proteomicsdb.org/prosit/
  2. Navigate to 'Spectral Library'
  3. Upload your target peptide list as described in the Prosit documentation
  4. **IMPORTANT**: under the 'Output format' header (just prior to submitting the Prosit task), ensure that **Generic Text** is selected

Once processing is complete, the resulting file can be used with Synthedia

When Prosit inputs are given, Synthedia models peptide precursor abundances using a Gaussian distribution of Log2 instnsities which approximates the disributions typically observed upon anlysis of real data. The parameters of this distribution (mean and standard deviation) can be specified with the ```--prosit_peptide_abundance_mean``` and ```--prosit_peptide_abundance_stdev``` parameters which have defaults of 22 and 3 respectively.


#### MaxQuant 'txt' directories

As an alternative to Prosit, Synthedia can read and simulate DIA data based upon the MaxQuant processing results of a file acquired using a Data-Dependent Acquisition (DDA) strategy. In this case, peptide fragment ions are generated based on the matched ions for a PSM reported in the MaxQuant ```msms.txt``` file from the ```Masses``` and ```Intensities``` columns. Note: these are only those fragment ions that MaxQuant assigns as matching a given peptide - they may not necessarily provide a 'full' sequence coverage and may not be correctly assigned in some cases.

Synthedia offers options to filter reverse and contaminant peptides as well as filter PSMs by Posterior Error Probability (PEP) values.


## DIA acquisition strategies

DIA-LC-MS/MS data can be acquired in many ways. The default invocation of synthedia creates non-overlapping, 30 Th windows between m/z 350 and m/z 1600. To simulate data using different DIA acquisition strategies, a file defining the acquisition schema can be supplied. An example acquisition schema file and blank template are provided in the ```templates``` directory.

## Decoy signals

LC-MS/MS analysis of bottom-up proteomics samples can be complicated by the presence of ions derived from non-peptide sample contaminants. To mimic this, decoy ions can be simulated together with peptide signals by specifying a decoy database in NIST '.msp' format. See the section 'NIST Text Format of Individual Spectra' in [this document](https://chemdata.nist.gov/mass-spc/ms-search/docs/Ver20Man_11.pdf) for details on the .msp format.

Custom .msp files can be specified or pre-prepared files can be dowloaded from [MS-DIAL](http://prime.psc.riken.jp/compms/msdial/main.html#MSP). The Synthedia web server uses the "All public MS/MS (13,303 unique compounds)" file from [MS-DIAL](http://prime.psc.riken.jp/compms/msdial/main.html#MSP).

The number of decoy peaks to simulate, as well as the maximum number of fragments to simulate per decoy, can be specified throught the command line arguments ```num_decoys``` and ```simulate_top_n_decoy_fragments```.

## Compatibility of Synthedia mzML files with other softeare

The mzML files generated with Synthedia have been tested with:

- DIA-NN
- EncyclopeDIA
- DIAUmppire
- Skyline
- OpenMS/ToppView
- MSConvert

The mzML files are known to not be compatible with MaxQuant (at least as at MaxQuant Version 2.0.3.0)

## Notes about the resulting data

### Peptide vs protein abundances
In a real experiment, all peptides from an up-regulated protein should be observed with increased abundance compared to a control group. Synthedia models ion abundances at the peptide level only. This means that, in the case of a two-group simulation, peptides from the same protein may have very different (even opposing) directions of abundance change between groups. This is primarily because the Prosit input type (which is preferred) does not contain mappings between peptides and input proteins and we wished to maintain the ability to simulated arbitrary peptide data.

As a result of this, if mzML files generated with synthedia are analysed with DIA analysis software, the main comparison in abundances should be made at the peptide level.

### Synthetic vs real data

While we have endeavoured to provide a range of options that allow for simulation of a broad range of chromatographic and mass spectral variables, many experimental processes are not modelled which will cause deviations between simulated and real data. As such, users are warned that simulations with Synthedia do not completely re-create experimental complexity and should be used as an investigational tool only.

As an example, Synthedia offers the capability to simulate the same set of precursors as if acquired on different length chromatographic gradients. This means that data acquired on a lengthy gradients could be reconstructed to approximate acquisition on shorter gradients. In these cases, Synthedia simply models spectra containing signals from many peptides as a simple superposition of their individual signals. In reality however, ion supression effects would result in data that may look substantially different from that which would be produced if a comparable experimental workflow were to be executed.

### Post-translational modifications

Synthedia currently has very limited support for simulating peptides with post-translational modifications. Currently, all cysteine residues are assumed to be modified by carbamidomethylation. No other post-translational modifications are supported.

## Usage

A basic invokation of synthedia using a MaxQuant input directory and accepting all other parameters as default is below:
```
synthedia --mq_txt_dir /path/to/MQ/results/txt
```
This will create a directory at the current working directory called 'output' containing the processing results.

Note that the default invocation will produce mzML files with profile peaks at both MS1 and MS2 levels. To create centroided data (which is significantly faster):
```
synthedia --mq_txt_dir /path/to/MQ/results/txt --centroid_ms1 --centroid_ms2
```

To create profile MS1 and centroid MS2 data:
```
synthedia --mq_txt_dir /path/to/MQ/results/txt --centroid_ms2
```

To plot diagnostic graphics and provide a custom path to an output directory:
```
synthedia --mq_txt_dir /path/to/MQ/results/txt --all --out_dir /path/to/out/directory
```

To provide a custom acquisition schema:
```
python synthedia --mq_txt_dir /path/to/MQ/results/txt --acquisition_schema /path/to/acquisition_schema.csv
```

To simulate MaxQuant reuslts on a shorter chromatographic gradient:
```
synthedia --mq_txt_dir /path/to/MQ/results/txt --new_run_length 10
```

To simulate a two group analysis with three replicates per group:
```
synthedia --mq_txt_dir /path/to/MQ/results/txt --n_groups 2 --samples_per_group 3
```

Rather than specifying arguments on the commmand line, processing parameters can by placed into a ```.yaml``` file and then supplied as:
```
synthedia --config /path/to/params.yaml
```

Additional or updated parameters can be specified:
```
synthedia --config /path/to/params.yaml  --ms1_resolution 100000
```
In this case the ```ms1_resolution``` parameter value given on the command line is used even if a different value is given in ```/path/to/params.yaml```. That is, the heirarchy is command line parameter > config.yaml > synthedia default.

To quickly see the impact of changing input parameters, a single peptide ion can be modelled from the command line and no additional input files are required (i.e. the requirement for a Prosit or MaxQuant input is removed):
```
synthedia --preview --preview_sequence PEPTIDE --preview_charge 2
```

For example, to model a broad chromatographic peak, the mean peak full-width at half-maximum intensity (FWHM) can be increased:
```
synthedia --preview --preview_sequence PEPTIDE --preview_charge 2 --rt_peak_fwhm_distribution_mean 10
```
A series of plots are generated in the directiry specified in ```--out_dir``` to visualise the change.

## Synthedia parameter reference
    usage: synthedia [-h] [--mq_txt_dir MQ_TXT_DIR] [--prosit PROSIT]
                 [--prosit_peptide_abundance_model PROSIT_PEPTIDE_ABUNDANCE_MODEL]
                 [--prosit_peptide_abundance_mean PROSIT_PEPTIDE_ABUNDANCE_MEAN]
                 [--prosit_peptide_abundance_stdev PROSIT_PEPTIDE_ABUNDANCE_STDEV]
                 [--prosit_peptide_abundance_emg_k PROSIT_PEPTIDE_ABUNDANCE_EMG_K]
                 [--acquisition_schema ACQUISITION_SCHEMA]
                 [--use_existing_peptide_file USE_EXISTING_PEPTIDE_FILE] [--out_dir OUT_DIR]
                 [--output_label OUTPUT_LABEL] [--config CONFIG] [--silent] [--write_params]
                 [--mq_pep_threshold MQ_PEP_THRESHOLD] [--filterTerm FILTERTERM]
                 [--num_processors NUM_PROCESSORS] [--ms1_min_mz MS1_MIN_MZ]
                 [--ms1_max_mz MS1_MAX_MZ] [--ms2_min_mz MS2_MIN_MZ] [--ms2_max_mz MS2_MAX_MZ]
                 [--ms1_resolution MS1_RESOLUTION] [--ms2_resolution MS2_RESOLUTION]
                 [--ms1_scan_duration MS1_SCAN_DURATION] [--ms2_scan_duration MS2_SCAN_DURATION]
                 [--isolation_window ISOLATION_WINDOW] [--resolution_at RESOLUTION_AT]
                 [--n_points_gt_fwhm N_POINTS_GT_FWHM] [--esi_instability ESI_INSTABILITY]
                 [--ms1_ppm_error_mean MS1_PPM_ERROR_MEAN]
                 [--ms1_ppm_error_stdev MS1_PPM_ERROR_STDEV]
                 [--ms2_ppm_error_mean MS2_PPM_ERROR_MEAN]
                 [--ms2_ppm_error_stdev MS2_PPM_ERROR_STDEV]
                 [--rt_peak_fwhm_distribution_model RT_PEAK_FWHM_DISTRIBUTION_MODEL]
                 [--rt_peak_fwhm_distribution_mean RT_PEAK_FWHM_DISTRIBUTION_MEAN]
                 [--rt_peak_fwhm_distribution_stdev RT_PEAK_FWHM_DISTRIBUTION_STDEV]
                 [--rt_peak_fwhm_distribution_emg_k RT_PEAK_FWHM_DISTRIBUTION_EMG_K]
                 [--min_rt_peak_fwhm MIN_RT_PEAK_FWHM] [--original_run_length ORIGINAL_RUN_LENGTH]
                 [--new_run_length NEW_RUN_LENGTH] [--rt_buffer RT_BUFFER]
                 [--rt_instability RT_INSTABILITY] [--ms1_min_peak_intensity MS1_MIN_PEAK_INTENSITY]
                 [--ms2_min_peak_intensity MS2_MIN_PEAK_INTENSITY] [--centroid_ms1] [--centroid_ms2]
                 [--write_empty_spectra] [--mz_peak_model MZ_PEAK_MODEL]
                 [--rt_peak_model RT_PEAK_MODEL] [--mz_emg_k MZ_EMG_K] [--rt_emg_k RT_EMG_K]
                 [--prob_missing_in_sample PROB_MISSING_IN_SAMPLE]
                 [--prob_missing_in_group PROB_MISSING_IN_GROUP] [--no_isotopes] [--tic] [--schema]
                 [--all] [--n_groups N_GROUPS] [--samples_per_group SAMPLES_PER_GROUP]
                 [--between_group_stdev BETWEEN_GROUP_STDEV]
                 [--within_group_stdev WITHIN_GROUP_STDEV] [--decoy_msp_file DECOY_MSP_FILE]
                 [--num_decoys NUM_DECOYS]
                 [--simulate_top_n_decoy_fragments SIMULATE_TOP_N_DECOY_FRAGMENTS]
                 [--decoy_abundance_mean DECOY_ABUNDANCE_MEAN]
                 [--decoy_abundance_stdev DECOY_ABUNDANCE_STDEV] [--preview]
                 [--preview_sequence PREVIEW_SEQUENCE] [--preview_charge PREVIEW_CHARGE]
                 [--preview_abundance PREVIEW_ABUNDANCE]

      Generate synthetic DIA LC-MS/MS bottom up proteomics data with known composition.

      optional arguments:
        -h, --help            show this help message and exit

      Input/Output:
        --mq_txt_dir MQ_TXT_DIR
                              Path to MaxQuat "txt" directory.
        --prosit PROSIT       Path to prosit prediction library.
        --prosit_peptide_abundance_model PROSIT_PEPTIDE_ABUNDANCE_MODEL
                              The model from which log2 abundance values are drawn to generate abundances
                              for prosit precursors. Can be "gaussian", "exponentially_modified_gaussian"
                              or "cauchy". Not used for MaxQuant input types.
        --prosit_peptide_abundance_mean PROSIT_PEPTIDE_ABUNDANCE_MEAN
                              Mean log2 abundance used to simulate peptide abundances for prosit input
                              types. Not used for MaxQuat input types.
        --prosit_peptide_abundance_stdev PROSIT_PEPTIDE_ABUNDANCE_STDEV
                              Standard deviation of gaussian used to simulate peptide abundances for prosit
                              input types. Not used for MaxQuant input types.
        --prosit_peptide_abundance_emg_k PROSIT_PEPTIDE_ABUNDANCE_EMG_K
                              Shape factor for exponentially modified gaussian from which peptide
                              abundances are drawn for prosit input types. Must be greater than 0.
                              Increasing this value results in a more heavily tailed distribution resulting
                              in more peptides higher intensities. Not used for MaxQuant input types
        --acquisition_schema ACQUISITION_SCHEMA
                              Path to file defining MS2 acquisition schema.
        --use_existing_peptide_file USE_EXISTING_PEPTIDE_FILE
                              Path to an existin peptide file which will be used.
        --out_dir OUT_DIR     Output directory where results should be written.
        --output_label OUTPUT_LABEL
                              Prefix for output files.
        --config CONFIG       Path to *.yaml config file.
        --silent              Do not print logging output to terminal
        --write_params        Write parameters to yaml file and exit.

      Filtering:
        --mq_pep_threshold MQ_PEP_THRESHOLD
                              For MaxQuant input data, use only peptides with a Posterior Error Probability
                              (PEP) less than this value
        --filterTerm FILTERTERM
                              Terms used to filter input maxquant lists to remove unwanted targets. For
                              example contaminant protein can be removed by specifying "--filterTerm CON_".
                              Multiple filters can be applied. For example "--filterTerm CON_ --filterTerm
                              REV_". Filters are used only for MaxQuant input types (no effect for Prosit)
                              and are applied to the "Proteins" column of the evidence.txt table.

      Processing:
        --num_processors NUM_PROCESSORS
                              Number of cores to use in constructing mzML files. Defaults to all available
                              cores

      Instrument Parameters:
        --ms1_min_mz MS1_MIN_MZ
                              Minimum m/z at MS1 level.
        --ms1_max_mz MS1_MAX_MZ
                              Maximum m/z at MS1 level.
        --ms2_min_mz MS2_MIN_MZ
                              Minimum m/z at MS2 level.
        --ms2_max_mz MS2_MAX_MZ
                              Maximum m/z at MS2 level.
        --ms1_resolution MS1_RESOLUTION
                              Mass spectral resolution at MS1 level.
        --ms2_resolution MS2_RESOLUTION
                              Mass spectral resolution at MS2 level.
        --ms1_scan_duration MS1_SCAN_DURATION
                              Time in seconds taken to record an MS1 scan.
        --ms2_scan_duration MS2_SCAN_DURATION
                              Time in seconds taken to record an MS2 scan.
        --isolation_window ISOLATION_WINDOW
                              Length of DIA window in m/z.
        --resolution_at RESOLUTION_AT
                              m/z value at which resolution is defined.
        --n_points_gt_fwhm N_POINTS_GT_FWHM
                              Number of MS data points greater than the peak FWHM. Increasing this number
                              means each mass spectral peak will be described by more data points but will
                              also slow processing time and increase file size.
        --esi_instability ESI_INSTABILITY
                              Simulates imperfection in chromatographic peaks by applying a randomly
                              intensity scaling factor to adjacent scans. A value of 0 indicates no
                              randomness. A value of 100 indicates high spray instability.
        --ms1_ppm_error_mean MS1_PPM_ERROR_MEAN
                              The mean value of a Gaussian distribution from which PPM errors for MS1
                              precursors will be drawn. This value can be negative. Setting both
                              ms1_ppm_error_mean, and ms1_ppm_error_stdev to 0 equates to perfect mass
                              accuracy.
        --ms1_ppm_error_stdev MS1_PPM_ERROR_STDEV
                              The standard deviation of a Gaussian distribution from which PPM errors for
                              MS1 precursors will be drawn. Setting both ms1_ppm_error_mean and
                              ms1_ppm_error_stdev to 0 equates to perfect mass accuracy.
        --ms2_ppm_error_mean MS2_PPM_ERROR_MEAN
                              The mean value of a Gaussian distribution from which PPM errors for MS2
                              fragments will be drawn. This value can be negative. Setting both
                              ms1_ppm_error_mean, and ms1_ppm_error_stdev to 0 equates to perfect mass
                              accuracy.
        --ms2_ppm_error_stdev MS2_PPM_ERROR_STDEV
                              The standard deviation of a Gaussian distribution from which PPM errors for
                              MS2 fragments will be drawn. Setting both ms1_ppm_error_mean and
                              ms1_ppm_error_stdev to 0 equates to perfect mass accuracy.

      Chromatography:
        --rt_peak_fwhm_distribution_model RT_PEAK_FWHM_DISTRIBUTION_MODEL
                              The model used to simulate chromatographic peaks. Can be "gaussian",
                              "exponentially_modified_gaussian" or "cauchy"
        --rt_peak_fwhm_distribution_mean RT_PEAK_FWHM_DISTRIBUTION_MEAN
                              Chromatographic peak full with at half maximum intehsity in seconds. If
                              --rt_peak_fwhm_stdev is non-zero, this value is the mean of a Gaussian
                              distribution from which FWHMs for individual peptides are drawn
        --rt_peak_fwhm_distribution_stdev RT_PEAK_FWHM_DISTRIBUTION_STDEV
                              Standard deviation of a Gaussian distribution from which chromatographic peak
                              full with at half maximum intehsity are drawn for indifivual peptides.
        --rt_peak_fwhm_distribution_emg_k RT_PEAK_FWHM_DISTRIBUTION_EMG_K
                              Shape factor for exponentially modified gaussian from which peak FWHMs are
                              drawn. Must be greater than 0. Increasing this value results in a more
                              heavily tailed distribution resulting in more peptides with broad
                              chromatographic profiles
        --min_rt_peak_fwhm MIN_RT_PEAK_FWHM
                              Minimum chromatographic peak full with at half maximum intehsity in seconds.
        --original_run_length ORIGINAL_RUN_LENGTH
                              Length in minutes of original data file. If not given, this will be
                              determined by taking the difference between the minimum and maximum peptide
                              retention times. If set to "0", the retention time range will be
                              automatically detected from the input data.
        --new_run_length NEW_RUN_LENGTH
                              Length in minutes of new data file. If set to "0", the retention time range
                              of the input data will be used.
        --rt_buffer RT_BUFFER
                              Time (in minutes) that should be appended to the beginning and end of the
                              retention time range of a set of input peptides. This helps ensure that
                              peptides at the boundaries of the elution window are simulated completely
        --rt_instability RT_INSTABILITY
                              Introduces an instability in retentention time values for the same peptide
                              when a multi-group or multi-sample simulation is conducted. The value is the
                              maximum number of seconds by which peptide retention times will differ

      Simulation:
        --ms1_min_peak_intensity MS1_MIN_PEAK_INTENSITY
                              Peptide elution profiles are simulated as gaussian peaks. This value sets the
                              minimum gaussian curve intensitiy for a peptide to be simulated in MS1
                              spectra.
        --ms2_min_peak_intensity MS2_MIN_PEAK_INTENSITY
                              Peptide elution profiles are simulated as gaussian peaks. This value sets the
                              minimum gaussian curve intensitiy for a peptide to be simulated in MS2
                              spectra.
        --centroid_ms1        If given, simulated MS1 mass spectra will be centroided. Otherwise, profile
                              data will be written.
        --centroid_ms2        If given, simulated MS2 mass spectra will be centroided. Otherwise, profile
                              data will be written.
        --write_empty_spectra
                              Write empty mass sepctra to the output data file
        --mz_peak_model MZ_PEAK_MODEL
                              The model used to simulate mass spectral peaks. Can be "gaussian",
                              "exponentially_modified_gaussian" or "cauchy"
        --rt_peak_model RT_PEAK_MODEL
                              The model used to simulate chromatographic peaks. Can be "gaussian",
                              "exponentially_modified_gaussian" or "cauchy"
        --mz_emg_k MZ_EMG_K   Shape factor for exponentially modified gaussian in the mass spectral domain.
                              Must be greater than 0. Increasing K results in more heavily tailed mass
                              spectral peaks. This parameter is inactive unless --mz_peak_model is not set
                              to exponentially_modified_gaussian.
        --rt_emg_k RT_EMG_K   Shape factor for exponentially modified gaussian in the retention time
                              domain. Must be greater than 0. Increasing K results in more heavily tailed
                              chromatographic peaks. This parameter is inactive unless --rt_peak_model is
                              not set to exponentially_modified_gaussian.
        --prob_missing_in_sample PROB_MISSING_IN_SAMPLE
                              Probability (0-100) that a peptide is missing in any given sample
        --prob_missing_in_group PROB_MISSING_IN_GROUP
                              Probability (0-100) that a peptide is missing in an entire group
        --no_isotopes         Disable simulation of non-monoisotopic ions. Can be used to increase
                              simulation speed when optimising parameters.

      Plotting:
        --tic                 Plot TIC for the generated mzML file.
        --schema              Plot acquisition schema.
        --all                 Plot all graphics.

      Grouping and quantitation:
        --n_groups N_GROUPS   Number of treatment groups to simulate.
        --samples_per_group SAMPLES_PER_GROUP
                              Number of individual samples to simulate per treatment group.
        --between_group_stdev BETWEEN_GROUP_STDEV
                              Standard deviation of a normal distribution from which group means will be
                              drawn.
        --within_group_stdev WITHIN_GROUP_STDEV
                              Standard deviation of a normal distribution from which within group samples
                              will be drawn.

      Decoys:
        --decoy_msp_file DECOY_MSP_FILE
                              Path to MSP file. Note - must include retention times.
        --num_decoys NUM_DECOYS
                              Number of decoy peaks to simulate
        --simulate_top_n_decoy_fragments SIMULATE_TOP_N_DECOY_FRAGMENTS
                              Simulate n most intense fragments of the decoy compound.
        --decoy_abundance_mean DECOY_ABUNDANCE_MEAN
                              Mean log2 abundance used to simulate decoy ion abundances.
        --decoy_abundance_stdev DECOY_ABUNDANCE_STDEV
                              Standard deviation of gaussian used to simulate decoy ion abundances.

      Preview:
        --preview             Simulate a single peptide without a Prosit/MaxQuant file. Useful for testing
                              parameter sets.
        --preview_sequence PREVIEW_SEQUENCE
                              Primary sequence of peptide to preview.
        --preview_charge PREVIEW_CHARGE
                              Charge of intact peptide to preview.
        --preview_abundance PREVIEW_ABUNDANCE
                              Abundance of peptide to simulate



## Viewing mzML files

The mzML files produced can be viewed in many different freely available software such as the [TOPPView](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_TOPPView.html) package which is part of [OpenMS](https://www.openms.de/).

![image](https://user-images.githubusercontent.com/16992583/160973872-4009c5cf-57c6-49a5-a868-edb21fa5e190.png)

## Support
- If you notice bugs or have suggestions for improvements, please contact us [here](https://github.com/mgleeming/synthedia/issues).
- Or better still, fork and submit a pull request :)
