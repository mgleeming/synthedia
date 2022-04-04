import math, random
import numpy as np
from pyteomics import mass, fasta

# constants
PROTON = 1.007276
IAA = 57.02092

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
            self.populate_from_prosit(prosit_entry)
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
        return

    def populate_from_msp(self, options, msp_entry, peptide_min, peptide_max):
        self.sequence = msp_entry['NAME']
        self.formula = msp_entry['FORMULA']
        self.ms2_peaks = msp_entry['fragments']
        if len(self.ms2_peaks) > 15:
            self.ms2_peaks.sort(key=lambda x: x[1], reverse = True)
            self.ms2_peaks = self.ms2_peaks[0:options.simulate_top_n_decoy_fragments]
        self.rt = float(msp_entry['RETENTIONTIME'])
        self.intensity = random.randint(peptide_min, peptide_max)
        self.charge = 1
        self.mz = float(msp_entry['PRECURSORMZ'])
        self.protein = 'DECOY'
        self.mass = self.mz - PROTON
        return

    def populate_from_prosit(self, prosit_entry):
        self.sequence = prosit_entry['StrippedPeptide'].iloc[0]
        self.charge = prosit_entry['PrecursorCharge'].iloc[0]

        self.mass = (prosit_entry['PrecursorMz'].iloc[0] * self.charge) - (self.charge * PROTON)
        self.mz = prosit_entry['PrecursorMz'].iloc[0]

        self.rt = prosit_entry['Retention time'].iloc[0] * 60
        self.protein = 'None'

        self.intensity = 100000000

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
        if self.decoy:
            for isotope in mass.mass.isotopologues( formula = self.formula, report_abundance = True, overall_threshold = 0.01, elements_with_isotopes = ['C']):
                calc_mz = (isotope[0].mass() + (IAA * self.sequence.count('C')) +  (PROTON * self.charge)) / self.charge
                isotopes.append( [calc_mz, isotope[1] * self.intensity])
        else:
            for isotope in mass.mass.isotopologues( sequence = self.sequence, report_abundance = True, overall_threshold = 0.01, elements_with_isotopes = ['C']):
                calc_mz = (isotope[0].mass() + (IAA * self.sequence.count('C')) +  (PROTON * self.charge)) / self.charge
                isotopes.append( [calc_mz, isotope[1] * self.intensity])
        self.ms1_isotopes = isotopes
        return

    def scale_retention_times(self, options):
        self.scaled_rt = self.rt / options.original_run_length * options.new_run_length
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

    def calculate_retention_length(self, options, ms_rts):
        rt_mask = np.where(
            (ms_rts > self.scaled_rt - (options.rt_clip_window * 60))
            &
            (ms_rts < self.scaled_rt + (options.rt_clip_window * 60))
        )
        ms_rts = ms_rts[rt_mask]
        ints = options.rt_peak_model(ms_rts, **{
            'mu': self.scaled_rt, 'sig': options.rt_stdev, 'emg_k': options.rt_emg_k
        })
        ints *= self.intensity
        mask = np.where(ints > options.ms2_min_peak_intensity)
        peak_rts = ms_rts[mask]
        self.min_scaled_peak_rt = min(peak_rts)
        self.max_scaled_peak_rt = max(peak_rts)
        return min(peak_rts), max(peak_rts)

def calculate_retention_lengths(options, peptides, spectra):

    ms_rts = np.asarray([s.rt for s in spectra])

    for p in peptides:
        p.calculate_retention_length(options, ms_rts)

    return peptides

def calculate_scaled_retention_times(options, peptides):

    for p in peptides:
        p.scale_retention_times(options)

    return peptides

