import os, sys, logging, random, pickle, multiprocessing
import pandas as pd
import numpy as np
from .peptide import SyntheticPeptide
from .mzml import Spectrum
from itertools import repeat

class NoPeptidesToSimulateError(Exception):
    pass

class IncorrectInputError(Exception):
    pass

class AcquisitionSchemaError(Exception):
    pass

class AcquisitionSchema (object):

    def __init__(self, options):

        if options.acquisition_schema:
            acquisition_schema_file = self.read_acquisition_schema_file(options)
            self.schema = self.make_schema_from_file(options, acquisition_schema_file)
        else:
            self.schema = self.make_default_schema(options)

        self.spectra = self.make_spectra(options)
        return

    def get_spectra(self):
        return self.spectra

    def read_acquisition_schema_file(self, options):
        logger = logging.getLogger("assembly_logger")
        logger.info('\tReading acquisition schema file: %s' %(
            options.acquisition_schema
        ))
        if not os.path.isfile(options.acquisition_schema):
            msg = 'The specified acquisition schema file was not found'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)
        try:
            df = pd.read_csv(options.acquisition_schema)
        except:
            msg = 'Error parsing acquisition schema file'
            logger.error(msg)
            logger.error('Exiting')
            raise IncorrectInputError(msg)

        expected_contents = [
            'ms_level','scan_duration_in_seconds',
            'isolation_window_lower_mz','isolation_window_upper_mz'
        ]
        for ec in expected_contents:
            if ec not in list(df):
                msg = 'Expected column %s not found in acquisition schema! Found columns : %s' %(
                    ec, ', '.join(list(df))
                )
                logger.error(msg)
                logger.error('Exiting')
                raise AcquisitionSchemaError(msg)

        return df

    def make_schema_from_file(self, options, acquisition_schema_file):
        schema = []
        logger = logging.getLogger("assembly_logger")
        logger.info('\tCreating acquisition schema from file: %s' %(
            options.acquisition_schema
        ))

        def validate(options, row):

            try:
                assert row['ms_level'] in [1,2]
            except AssertionError:
                msg = 'ms_level must be 1 or 2! Found %s' %row['ms_level']
                logger.error(msg)
                logger.error('Exiting')
                raise AcquisitionSchemaError(msg)

            try:
                assert row['scan_duration_in_seconds'] > 0
            except AssertionError:
                msg = 'Scan duration must be greater than 0! Found %s' %(
                    row['scan_duration_in_seconds']
                )
                logger.error(msg)
                logger.error('Exiting')
                raise AcquisitionSchemaError(msg)

            if int(row['ms_level']) == 2:
                try:
                    assert float(row['isolation_window_lower_mz']) > options.ms1_min_mz
                except AssertionError:
                    msg = 'isolation_window_lower_mz must be greater than ms1_min_mz'
                    logger.error(msg)
                    logger.error('Exiting')
                    raise AcquisitionSchemaError(msg)

                try:
                    assert float(row['isolation_window_upper_mz']) < options.ms1_max_mz
                except AssertionError:
                    msg = 'isolation_window_upper_mz must be lower than ms1_max_mz'
                    logger.error(msg)
                    logger.error('Exiting')
                    raise AcquisitionSchemaError(msg)

            return

        for index, row in acquisition_schema_file.iterrows():

            validate(options, row)

            schema.append({
                'order': int(row['ms_level']),
                'length': float(row['scan_duration_in_seconds']),
                'isolation_range': [
                    float(row['isolation_window_lower_mz']),
                    float(row['isolation_window_upper_mz'])
                ]
            })
        return schema

    def make_default_schema(self, options):
        logger = logging.getLogger("assembly_logger")
        logger.info('\tCreating default acquisition schema')
        schema = [{
            'order': 1, 'length': options.ms1_scan_duration, 'isolation_range': None
        }]
        for i in range(options.ms1_min_mz, options.ms1_max_mz, options.isolation_window):
            schema.append({
                'order': 2,
                'length': options.ms2_scan_duration,
                'isolation_range': [i, i + options.isolation_window]
            })
        return schema

    def make_spectra(self, options):
        spectra = []
        total_run_time = 0
        synthedia_id_counter = 0
        while total_run_time < options.new_run_length + 2 * options.rt_buffer * 60:
            for entry in self.schema:
                spectra.append(
                    Spectrum(
                        synthedia_id_counter,
                        total_run_time, entry['order'],
                        entry['isolation_range'],
                        options
                    )
                )
                total_run_time += entry['length']
                synthedia_id_counter += 1
        return spectra

def simulate_isotope_patterns(options, peptide_subset):
    for p in peptide_subset:
        p.get_ms1_isotope_pattern(options)
    return peptide_subset

class InputReader (object):

    def __init__(self, options):
        logger = logging.getLogger("assembly_logger")


        self.peptides = []

        if options.use_existing_peptide_file:
            self.peptides = self.read_pickle(options)

            # need to check that num_gorups and sample per group is the same as that used for the pickel file
            # -- these figures are used to generate the sample and group abundance offsets

            n_groups_pickle = len(self.peptides[0].peptide_abundance_offsets_between_groups)
            n_samples_pickle = len(self.peptides[0].sample_abundance_offsets[0])

            try:
                assert n_groups_pickle == options.n_groups
                assert n_samples_pickle == options.samples_per_group
            except AssertionError:
                msg = '''The supplied peptide pickle file has different number of groups and/or samples per group compared to the requested parameters.\n \
                    Peptide file n_groups = %s, requested n_groups = %s \n \
                    Peptide file samples_per_group = %s, requested samples_per_group = %s''' %(
                    n_groups_pickle,
                    options.n_groups,
                    n_samples_pickle,
                    options.samples_per_group
                )
                logger.error(msg)
                logger.error('Exiting')
                raise IncorrectInputError(msg)

        else:
            if options.mq_txt_dir:
                logger.info('Reading MaxQuant directory')
                self.peptides = self.read_peptides_from_mq(options)
            elif options.prosit:
                logger.info('Reading Prosit directory')
                self.peptides = self.read_peptides_from_prosit(options)

            if options.decoy_msp_file:
                if options.num_decoys > 0:
                    logger.info('Reading decoy file')
                    decoys = self.read_decoys_from_msp(options, self.peptides)
                    self.peptides = self.peptides + decoys

            logger.info('Simulating isotope patterns')
            if options.num_processors == 1:
                for p in self.peptides:
                    p.get_ms1_isotope_pattern(options)
            else:
                pool = multiprocessing.Pool(processes = options.num_processors)

                # split work into equal sized lists
                arg_sets = [[] for _ in range(options.num_processors)]
                counter = 0
                for p in self.peptides:
                    arg_sets[counter].append(p)
                    counter += 1
                    if counter == options.num_processors:
                        counter = 0

                # send work to procs and collect results
                self.peptides = []
                for _ in pool.starmap(simulate_isotope_patterns, zip(repeat(options), arg_sets)):
                    self.peptides.extend(list(_))

                pool.close()
                pool.join()

        if len(self.peptides) == 0:
            msg = 'No peptides or decoys created - Nothing to simulate!'
            logger.error(msg)
            logger.error('Exiting')
            raise NoPeptidesToSimulateError(msg)

        self.write_pickle(options)
        return

    def get_peptides(self):
        return self.peptides

    def write_pickle(self, options):
        with open( os.path.join(options.out_dir, 'peptides.pickle') , 'wb') as handle:
            pickle.dump(self.peptides, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return

    def read_pickle(self, options):

        if not os.path.isfile(options.use_existing_peptide_file):
            msg = 'The specified peptide file was not found'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)

        with open( options.use_existing_peptide_file , 'rb') as handle:
            peptides = pickle.load(handle)

        return peptides

    def read_peptides_from_mq(self, options):

        logger = logging.getLogger("assembly_logger")

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

    def read_peptides_from_prosit(self, options):

        logger = logging.getLogger("assembly_logger")

        if not os.path.isfile(os.path.join(options.prosit)):
            msg = 'The specified Prosit library file does not exist'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)

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

        n_groups = len(list(prosit.groupby(['ModifiedPeptide', 'PrecursorCharge'])))

        for _, peptide in prosit.groupby(['ModifiedPeptide']):

            # simulate abundances
            peptide_abundance_offsets_between_groups, sample_abundance_offsets = generate_group_and_sample_abundances(options)
            found_in_group, found_in_sample = generate_group_and_sample_probabilities(options)

            for __, precursor in peptide.groupby(['PrecursorCharge']):

                counter += 1
                if counter % 100 == 0:
                    logger.info('Constructing peptide %s of %s' %(counter, n_groups))

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

    def read_decoys_from_msp(self, options, peptides):

        logger = logging.getLogger("assembly_logger")

        peptide_intensities = [p.intensity for p in peptides]

        if not os.path.isfile(os.path.join(options.decoy_msp_file)):
            msg = 'The specified Decoy MSP file does not exist'
            logger.info(msg)
            logger.info('Exiting')
            raise IncorrectInputError(msg)

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

                    if float(mz) < options.ms2_min_mz: continue
                    if float(mz) > options.ms2_max_mz: continue

                    lipid_dict['fragments'].append([float(mz), float(intensity)])

            # check for missing fields that are required later
            if len(lipid_dict['fragments']) == 0: continue
            if 'RETENTIONTIME' not in lipid_dict.keys(): continue
            if 'NAME' not in lipid_dict.keys(): continue
            if 'PRECURSORMZ' not in lipid_dict.keys(): continue
            if 'FORMULA' not in lipid_dict.keys(): continue

            # some formulae have charges etc - skip these
            if not lipid_dict['FORMULA'].isalnum(): continue

            # check if precursor out of bounds
            if (float(lipid_dict['PRECURSORMZ']) < options.ms1_min_mz) or (float(lipid_dict['PRECURSORMZ']) > options.ms1_max_mz):
                continue

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

    sample_abundance_offsets[0][0] = 0
    return peptide_abundance_offsets_between_groups, sample_abundance_offsets

def generate_group_and_sample_probabilities(options):

    if options.prob_missing_in_group == 100:
        found_in_group = [0 for _ in range(options.n_groups)]
    else:
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
