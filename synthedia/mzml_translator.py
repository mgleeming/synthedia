import os, copy, random, math, sys
from pyopenms import *

import time
class MZMLReader():
    def __init__(self):
        self.od_exp = OnDiscMSExperiment()
        self.od_exp.openFile('/home/mleeming/share/synthedia_tests/32_MB_small_subset/220329_MB-4007_08_filter.mzML')
        self.num_spectra = self.od_exp.getNrSpectra()

    def __iter__(self):
        self.current_spec = 0
        return self

    def __next__(self):
        if self.current_spec < self.num_spectra:
            s = self.od_exp.getSpectrum(self.current_spec)
            self.current_spec += 1
            return Spectrum(s)
        else:
            self.close()
            raise StopIteration

    def close(self):
        del self.od_exp
        return

class MZMLWriter():

    def __init__(self):#, 'test.mzML', 67677):
        self.n_spectra = 67677
        self.consumer = PlainMSDataWritingConsumer('/home/mleeming/share/synthedia_tests/33_MB_target_channel_only/220329_MB-4007_08_filter.mzML')
        self.consumer.setExpectedSize(self.n_spectra,0)

        print('set instrument\n\n')
        experiment = ExperimentalSettings()
        instrument = Instrument()

        # set ionisation source
        source = IonSource()
        source.setOrder(1)
        source.setIonizationMethod(20)
        instrument.setIonSources([source])

        # set analyser
        # see pyopenms.pyopenms_6.__AnalyzerType.__dict__
        # for mapping of integers to analyser types
        analysers = []
        for analyser_id in [1,13]:
            analyser = MassAnalyzer()
            analyser.setType(analyser_id)
            analyser.setOrder(2)
            analysers.append(analyser)
        instrument.setMassAnalyzers(analysers)

        # set detector
        detector = IonDetector()
        detector.setOrder(3)
        detector.setType(20)
        instrument.setIonDetectors([detector])

        experiment.setInstrument(instrument)
        self.consumer.setExperimentalSettings(experiment)
        self.n_spec_written = 0

        self.ghost_ms1_data = {}
        self.ghost_ms2_data = {}

#        self.get_template_spectra()
        return

    def get_template_spectra(self):
        od_exp = OnDiscMSExperiment()
        od_exp.openFile( '/home/mleeming/Code/synthedia/data/220329_MB-4007_08.mzML' )
        self.ms1_template = None
        self.ms2_template = None
        print('getting template spectra')
        for k in range(od_exp.getNrSpectra()):
            s = od_exp.getSpectrum(k)
            lvl = s.getMSLevel()
            if all([self.ms1_template, self.ms2_template]):
                break
            if lvl == 1:
                self.ms1_template = s
            if lvl == 2:
                self.ms2_template = s
        del od_exp
        return

    def write_spec(self, spec):

        spec_to_write = MSSpectrum()
#        if spec.order == 1:
#            spec_to_write = copy.deepcopy(self.ms1_template)
#        if spec.order == 2:
#            spec_to_write = copy.deepcopy(self.ms2_template)

        ints = spec.ints
        mzs = spec.mzs

        spec_to_write.set_peaks([mzs, ints])

        spec_to_write.setType(spec.spectrum.getType())

        spec_to_write.setRT(spec.rt)
        spec_to_write.setMSLevel(spec.order)

        if spec.order == 2:
            p = Precursor()
            p.setMZ((spec.isolation_hl + spec.isolation_ll) / 2)
            p.setIsolationWindowLowerOffset(spec.lower_offset)
            p.setIsolationWindowUpperOffset(spec.upper_offset)

            m =(spec.isolation_hl + spec.isolation_ll) / 2
            lo = spec.lower_offset
            ho = spec.upper_offset

            if (m > 595) and (m < 601):
                pass
            else:
                return

            print('%.2f\t%.2f\t%.2f\t%.2f' %(m - lo, m, m+ho, lo+ho))

            p.setActivationMethods(set([0]))
            p.setActivationEnergy(30)

            spec_to_write.setPrecursors( [p] )

        # Add scan metadata
        # - some software will crash reading the file without this
        instrument_settings = InstrumentSettings()
        instrument_settings.setPolarity(1)

        scan_window = ScanWindow()

        if spec.order == 1:
            scan_window.begin = 350
            scan_window.end = 1400
        elif spec.order == 2:
            scan_window.begin =200
            scan_window.end = 2000

        instrument_settings.setScanWindows([scan_window])
        spec_to_write.setInstrumentSettings( instrument_settings)

        spec_to_write.setMetaValue('lowest observed m/z', mzs[0])
        spec_to_write.setMetaValue('highest observed m/z', mzs[-1])
        spec_to_write.setMetaValue('base peak intensity', float(max(ints)))
        spec_to_write.setMetaValue('base peak m/z', mzs[np.argmax(ints)])
        spec_to_write.setMetaValue('total ion current', float(sum(ints)))

        spec_to_write.setMetaValue('filter string', spec.spectrum.getMetaValue(b'filter string'))

        spec_to_write.setNativeID(
            'controllerType=0 controllerNumber=1 scan=%s' % self.n_spec_written
        )
        spec_to_write.updateRanges()

        self.consumer.consumeSpectrum(spec_to_write)
        self.n_spec_written += 1
        del spec_to_write
        return

    def close(self):
        del self.consumer
        return

class Spectrum():
    def __init__(self, s):

        self.spectrum = s
        self.rt = s.getRT()
        self.order = s.getMSLevel()
        self.mzs, self.ints = s.get_peaks()

        if self.order == 2:
            self.lower_offset = s.getPrecursors()[0].getIsolationWindowLowerOffset()
            self.upper_offset = s.getPrecursors()[0].getIsolationWindowUpperOffset()
            self.isolation_ll = s.getPrecursors()[0].getMZ() - self.lower_offset
            self.isolation_hl = s.getPrecursors()[0].getMZ() + self.upper_offset

        return

    def clear(self):
        # save memory
        del self.mzs
        del self.ints
        del self.indicies
        return


writer = MZMLWriter()
reader = MZMLReader()


for i , r in enumerate(reader):
    writer.write_spec(r)
writer.close()
