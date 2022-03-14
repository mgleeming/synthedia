import math
from pyteomics import mass, fasta

# constants
PROTON = 1.007276
IAA = 57.02092

class SyntheticPeptide():

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


