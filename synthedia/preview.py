import os, random, time
import pandas as pd
from pyteomics import mass, fasta

class PreviewPeptide (object):

    def __init__(self, options):

        name = 'preview_prosit.dat'

        self.out_file = os.path.join(options.out_dir, name)

        self.generate_prosit_table(options)
        self.save_dataframe()
        return

    def generate_prosit_table(self, options):
        self.prosit_df_items = []

        intact_mz = mass.fast_mass(
            options.preview_sequence, charge = options.preview_charge
        )

        for fragment in self.fragments(options.preview_sequence):
            fragment_mz, ion_type, fragment_charge, fragment_position = fragment

            entry = {
                'RelativeIntensity' : 1,
                'FragmentMz' : fragment_mz,
                'ModifiedPeptide' : '_%s_'%options.preview_sequence,
                'LabeledPeptide' : options.preview_sequence,
                'StrippedPeptide' : options.preview_sequence,
                'PrecursorCharge' : options.preview_charge,
                'PrecursorMz' : intact_mz,
                'iRT' : 30,
                'proteotypicity' : 1,
                'FragmentNumber' : fragment_position,
                'FragmentType' : ion_type,
                'FragmentCharge' : fragment_charge,
                'FragmentLossType' : 'noloss'
            }
            self.prosit_df_items.append(entry)

        self.df = pd.DataFrame(self.prosit_df_items)
        return

    def get_dataframe_path(self):
        return self.out_file

    def save_dataframe(self):
        self.df.to_csv(self.out_file, index = False)
        return

    def fragments(self, peptide, types=('b', 'y'), maxcharge=1):
        for i in range(1, len(peptide)-1):
            for ion_type in types:
                for charge in range(1, maxcharge+1):
                    if ion_type[0] in 'abc':
                        mz =  mass.fast_mass(
                                peptide[:i], ion_type=ion_type, charge=charge)
                        fragment_number = len(peptide[:i])
                    else:
                        mz = mass.fast_mass(
                                peptide[i:], ion_type=ion_type, charge=charge)
                        fragment_number = len(peptide[i:])
                    yield mz, ion_type, charge, fragment_number

