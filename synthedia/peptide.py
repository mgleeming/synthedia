import math, random, sys, copy
import numpy as np
from pyteomics import mass, fasta

# constants
PROTON = 1.007276
IAA = 57.02092
ONE_MILLION = 1000000

class Peak():
    def __init__(self, options, mz, intensity, ms_level = 1):
        self.ms_level = ms_level
        self.base_mz = mz
        self.intensity = intensity
        self.make_peak_intensity_lists(options)
        self.make_mz_lists(options)
        self.make_profile_intensity_lists(options)
        return

    def make_profile_intensity_lists(self, options):
        self.profile_intensity_list = []
        for groupi in range(options.n_groups):
            self.profile_intensity_list.append([])
        return

    def make_mz_lists(self, options):
        self.peak_mz_list = []
        self.peak_mz_error_list = []
        for groupi in range(options.n_groups):
            self.peak_mz_list.append([])
            self.peak_mz_error_list.append([])
            for samplei in range(options.samples_per_group):
                if self.ms_level == 1:
                    error = np.random.normal(
                        loc = options.ms1_ppm_error_mean,
                        scale = options.ms1_ppm_error_stdev
                    )
                else:
                    error = np.random.normal(
                        loc = options.ms2_ppm_error_mean,
                        scale = options.ms2_ppm_error_stdev
                    )
                offset = error * self.base_mz / ONE_MILLION
                mzi = self.base_mz + offset
                self.peak_mz_list[-1].append(mzi)
                self.peak_mz_error_list[-1].append(error)
        return

    def make_peak_intensity_lists(self, options):
        self.total_ms1_intensity = []
        self.max_ms1_intensity = []

        self.total_fragment_intensity = []
        self.max_fragment_intensity = []

        for groupi in range(options.n_groups):
            self.total_ms1_intensity.append([])
            self.max_ms1_intensity.append([])

            self.total_fragment_intensity.append([])
            self.max_fragment_intensity.append([])
            for samplei in range(options.samples_per_group):
                self.total_ms1_intensity[-1].append([])
                self.max_ms1_intensity[-1].append([])

                self.total_fragment_intensity[-1].append([])
                self.max_fragment_intensity[-1].append([])
        return

    def update_ms1_intensities_to_report(self, groupi, samplei, peak_ints):
        self.total_ms1_intensity[groupi][samplei].append(peak_ints.sum())
        self.max_ms1_intensity[groupi][samplei].append(peak_ints.max())
        return

    def update_fragment_intensities_to_report(self, groupi, samplei, peak_ints):
        self.total_fragment_intensity[groupi][samplei].append(peak_ints.sum())
        self.max_fragment_intensity[groupi][samplei].append(peak_ints.max())
        return

    def get_total_ms1_intensity_to_report(self, groupi, samplei):
        return sum(self.total_ms1_intensity[groupi][samplei])

    def get_max_ms1_intensity_to_report(self, groupi, samplei):
        return max(self.max_ms1_intensity[groupi][samplei])

    def get_total_fragment_intensity_to_report(self, groupi, samplei):
        return sum(self.total_fragment_intensity[groupi][samplei])

    def get_max_fragment_intensity_to_report(self, groupi, samplei):
        return max(self.max_fragment_intensity[groupi][samplei])

    def get_limits(self, options, mzs, groupi, samplei):
        try:
            return self.lower_limit, self.higher_limit, self.indicies
        except AttributeError:
            mz_mask = np.where(
                (mzs > self.peak_mz_list[groupi][samplei] - options.ms_clip_window)
                &
                (mzs < self.peak_mz_list[groupi][samplei] + options.ms_clip_window)
            )
            try:
                self.lower_limit = mz_mask[0].min()
                self.higher_limit = mz_mask[0].max()
            except:
                print('error getting limits - exiting')
                sys.exit()
            self.indicies = mz_mask[0]
            return self.lower_limit, self.higher_limit, self.indicies

    def set_peak_intensities(self, options, peak_intensities, groupi, samplei):
        self.profile_intensity_list[groupi].append(peak_intensities)
        return

    def get_peak_intensities(self, groupi, samplei):
        return copy.deepcopy(self.profile_intensity_list[groupi][samplei])

class SyntheticPeptide():
    def __init__(self,
            options,
            evidence_entry = None,
            msms_entry = None,
            prosit_entry = None,
            msp_entry = None,
            peptide_abundance_offsets_between_groups = None,
            sample_abundance_offsets = None,
            found_in_group = None,
            found_in_sample = None
        ):

        self.decoy = False
        if (evidence_entry is not None) and (msms_entry is not None):
            self.populate_from_mq(evidence_entry, msms_entry)
        elif prosit_entry is not None:
            self.populate_from_prosit(options, prosit_entry)
        elif msp_entry is not None:
            self.decoy = True
            self.populate_from_msp(options, *msp_entry)
        else:
            print('Insufficient data to construct peptides. Exiting')
            sys.exit()

        self.ms1_isotopes = None
        self.peptide_abundance_offsets_between_groups = peptide_abundance_offsets_between_groups
        self.sample_abundance_offsets = sample_abundance_offsets

        self.found_in_sample = found_in_sample
        self.found_in_group = found_in_group

        self.scale_retention_times(options)
        self.set_abundances()
        self.configure_ms2_peaks(options)
        self.make_points_per_peak_dict(options)

        self.min_peak_rt_list = []
        self.max_peak_rt_list = []
        for groupi in range(options.n_groups):
            self.min_peak_rt_list.append([])
            self.max_peak_rt_list.append([])
            for samplei in range(options.samples_per_group):
                self.min_peak_rt_list[groupi].append(None)
                self.max_peak_rt_list[groupi].append(0)
        return

    def get_min_peak_rt(self, group, sample):
        if self.min_peak_rt_list[group][sample] == None:
            return 0
        else:
            return self.min_peak_rt_list[group][sample]

    def get_max_peak_rt(self, group, sample):
        return self.max_peak_rt_list[group][sample]

    def update_peptide_retention_times(self, rt, groupi, samplei):
        if self.min_peak_rt_list[groupi][samplei] == None:
            self.min_peak_rt_list[groupi][samplei] = rt
        self.max_peak_rt_list[groupi][samplei] = rt
        return

    def make_points_per_peak_dict(self, options):
        self.points_per_peak_dict = {}
        for groupi in range(options.n_groups):
            for samplei in range(options.samples_per_group):
                for ms_level in [1,2]:
                    self.points_per_peak_dict['%s_%s_%s'%(groupi, samplei, ms_level)] = 0
        return

    def increment_points_per_peak_dict(self, groupi, samplei, ms_level):
        self.points_per_peak_dict['%s_%s_%s'%(groupi, samplei, ms_level)] += 1
        return

    def get_points_per_peak(self, groupi, samplei, ms_level):
        return self.points_per_peak_dict['%s_%s_%s'%(groupi, samplei, ms_level)]

    def get_total_precursor_intensity(self, groupi, samplei):
        try:
            return sum([
                peak.get_total_ms1_intensity_to_report(
                    groupi, samplei
                ) for peak in self.ms1_isotopes
            ])
        except ValueError:
            # occurs if no MS1 intensity
            # i.e. no MS1 spectra
            return 0

    def get_max_precursor_intensity(self, groupi, samplei):
        try:
            return max([
                peak.get_max_ms1_intensity_to_report(
                    groupi, samplei
                ) for peak in self.ms1_isotopes
            ])
        except ValueError:
            # occurs if no MS1 intensity
            # i.e. no MS1 spectra
            return 0

    def get_total_fragment_intensity(self, groupi, samplei):
        try:
            return self.ms2_peaks[self.max_fragment_index].get_total_fragment_intensity_to_report(
                groupi, samplei
            )
        except:
            return 0

    def get_max_fragment_intensity(self, groupi, samplei):
        try:
            return self.ms2_peaks[self.max_fragment_index].get_max_fragment_intensity_to_report(
                groupi, samplei
            )
        except:
            return 0

    def configure_ms2_peaks(self, options):
        self.ms2_peaks = [
            Peak(options, p[0], p[1], ms_level = 2) for p in self.ms2_peaks if all([p[0] > options.ms2_min_mz, p[0] < options.ms2_max_mz])
        ]

        self.max_fragment_index = np.argmax([_.intensity for _ in self.ms2_peaks])
        self.max_fragment_mz = self.ms2_peaks[self.max_fragment_index].base_mz
        return

    def populate_from_msp(self, options, msp_entry, peptide_min, peptide_max):
        self.sequence = msp_entry['NAME']
        self.formula = msp_entry['FORMULA']
        self.ms2_peaks = msp_entry['fragments']

        if len(self.ms2_peaks) > options.simulate_top_n_decoy_fragments:
            self.ms2_peaks.sort(key=lambda x: x[1], reverse = True)
            self.ms2_peaks = self.ms2_peaks[0:options.simulate_top_n_decoy_fragments]

        self.rt = float(msp_entry['RETENTIONTIME'])

        self.intensity = 2 ** np.random.normal(
            loc = options.decoy_abundance_mean,
            scale = options.decoy_abundance_stdev
        )

        self.charge = 1
        self.mz = float(msp_entry['PRECURSORMZ'])
        self.protein = 'DECOY'
        self.mass = self.mz - PROTON
        return

    def populate_from_prosit(self, options, prosit_entry):
        self.sequence = prosit_entry['StrippedPeptide'].iloc[0]
        self.charge = prosit_entry['PrecursorCharge'].iloc[0]

        self.mass = (prosit_entry['PrecursorMz'].iloc[0] * self.charge) - (self.charge * PROTON)
        self.mz = prosit_entry['PrecursorMz'].iloc[0]

        self.rt = prosit_entry['Retention time'].iloc[0] * 60
        self.protein = 'None'

        self.intensity = 2 ** np.random.normal(
            loc = options.prosit_peptide_abundance_mean,
            scale = options.prosit_peptide_abundance_stdev
        )

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

    def get_ms1_isotope_pattern(self, options):
        isotopes = []
        if self.decoy:
            for isotope in mass.mass.isotopologues( formula = self.formula, report_abundance = True, overall_threshold = 0.01, elements_with_isotopes = ['C']):
                calc_mz = (isotope[0].mass() + (IAA * self.sequence.count('C')) +  (PROTON * self.charge)) / self.charge
                isotopes.append(
                    Peak(
                        options, calc_mz, isotope[1] * self.intensity
                    )
                )
        else:
            for isotope in mass.mass.isotopologues( sequence = self.sequence, report_abundance = True, overall_threshold = 0.01, elements_with_isotopes = ['C']):
                calc_mz = (isotope[0].mass() + (IAA * self.sequence.count('C')) +  (PROTON * self.charge)) / self.charge
                isotopes.append(
                    Peak(
                        options, calc_mz, isotope[1] * self.intensity
                    )
                )
        self.ms1_isotopes = isotopes
        return

    def scale_retention_times(self, options):

        self.scaled_base_rt = self.rt / options.original_run_length * options.new_run_length

        self.scaled_rt_lists = []

        for group in range(options.n_groups):
            self.scaled_rt_lists.append([])
            for sample in range(options.samples_per_group):
                sample_rt = self.rt - random.randint(0, options.rt_instability * 10000) / 10000
                self.scaled_rt_lists[-1].append( sample_rt / options.original_run_length * options.new_run_length )
        return

    def set_abundances(self):
        self.abundances = []
        self.offsets = []
        for groupi, group_offset in enumerate(self.peptide_abundance_offsets_between_groups):
            self.abundances.append([])
            self.offsets.append([])
            for samplei, sample_offset in enumerate(self.sample_abundance_offsets[groupi]):

                if self.found_in_sample[groupi][samplei] == 0:
                    self.abundances[groupi].append(0)
                    self.offsets[groupi].append(0)
                else:
                    log_int = math.log2(self.intensity)
                    adjusted_log2_int = log_int + sample_offset
                    adjusetd_raw_int = 2 ** adjusted_log2_int
                    self.abundances[groupi].append(adjusetd_raw_int)
                    self.offsets[groupi].append(sample_offset)

        return

    def calculate_retention_length(self, options, all_ms_rts, all_ids):

        self.min_scaled_peak_rt_list = []
        self.max_scaled_peak_rt_list = []
        self.intensity_scale_factor_list = []

        for group in range(options.n_groups):

            self.min_scaled_peak_rt_list.append([])
            self.max_scaled_peak_rt_list.append([])
            self.intensity_scale_factor_list.append([])

            for sample in range(options.samples_per_group):

                # get rts within clip window
                rt_mask = np.where(
                    (all_ms_rts > self.scaled_rt_lists[group][sample] - options.rt_clip_window)
                    &
                    (all_ms_rts < self.scaled_rt_lists[group][sample] + options.rt_clip_window)
                )

                # subset of mass spectral rts in window
                ms_rts = all_ms_rts[rt_mask]

                # equivalent spectral index subset
                ids_subset = all_ids[rt_mask]

                # base peak model
                model_ints = options.rt_peak_model(ms_rts, **{
                    'mu': self.scaled_rt_lists[group][sample], 'sig': options.rt_stdev, 'emg_k': options.rt_emg_k
                })

                # multiply by peptide intensity
                ints = self.intensity * model_ints

                # mask of above threshold points
                mask = np.where(ints > options.ms2_min_peak_intensity)

                # rts of spectra above threshold
                peak_rts = ms_rts[mask]


                # equivalent mask of spectral indicies
                ids_subset_2 = ids_subset[mask]

                # base model intensities within threshold window
                # -- used to derive intensity scale factors
                model_ints = model_ints[mask]

                self.min_scaled_peak_rt_list[-1].append(min(peak_rts))
                self.max_scaled_peak_rt_list[-1].append(max(peak_rts))

                self.intensity_scale_factor_list[-1].append({ids_subset_2[i]:model_ints[i] for i in range(len(model_ints))})

        return

def calculate_retention_lengths(options, peptides, spectra):

    ms_rts = np.asarray([s.rt for s in spectra])
    ids = np.asarray([s.synthedia_id for s in spectra])

    assert len(ms_rts) == len(ids)

    for p in peptides:
        p.calculate_retention_length(options, ms_rts, ids)

    return peptides

def calculate_scaled_retention_times(options, peptides):

    for p in peptides:
        p.scale_retention_times(options)

    return peptides

def calculate_feature_windows(options, peptides, spectra):

    # min intensity of any peak
    min_intensity = min([
        options.ms1_min_peak_intensity, options.ms2_min_peak_intensity
    ])
    intensity_scale_factor = 1

    # need to find most abundant peak of any peptide and any ms level

    # get peptide with most abundant base intensity
    max_int_peptide_index = np.argmax([p.intensity for p in peptides])
    max_int_peptide = peptides[max_int_peptide_index]

    # find sample with most positive abundance offset between groups/replicates
    max_abundance_offset = max([max(_) for _ in max_int_peptide.offsets])

    # create list of all ms1 and ms2 peaks for this peptide
    all_peaks = max_int_peptide.ms1_isotopes + max_int_peptide.ms2_peaks

    # get most abundant peak
    max_int_peak_index = np.argmax([peak.intensity for peak in all_peaks])
    max_int_peak = all_peaks[max_int_peak_index]

    # calculate most abundant final intensity
    log_int = math.log2(max_int_peak.intensity)
    adjusted_log2_int = log_int + max_abundance_offset
    adjusetd_raw_int = 2 ** adjusted_log2_int
    peak_int = adjusetd_raw_int * intensity_scale_factor

    # calculate rt clip window
    ms_rts = np.asarray([s.rt for s in spectra])

    ints = options.rt_peak_model(ms_rts, **{
        'mu': max_int_peptide.scaled_base_rt, 'sig': options.rt_stdev, 'emg_k': options.rt_emg_k
    })

    ints *= peak_int
    mask = np.where(ints > min_intensity)
    peak_rts = ms_rts[mask]

    options.rt_clip_window = max(peak_rts) - min(peak_rts)

    point_diff = min([options.ms1_point_diff, options.ms2_point_diff])
    min_mz = min([options.ms1_min_mz, options.ms2_min_mz])
    max_mz = max([options.ms1_max_mz, options.ms2_max_mz])
    mzs = np.arange(min_mz, max_mz, point_diff)
    max_stdev = max([options.ms1_stdev, options.ms2_stdev])

    peak_ints = options.mz_peak_model(mzs, **{
        'mu': max_int_peak.base_mz, 'sig': max_stdev, 'emg_k': options.mz_emg_k
    })

    mask = np.where(peak_ints > min_intensity)
    peak_mzs = mzs[mask]

    options.ms_clip_window = max(peak_mzs) - min(peak_mzs)
    return
