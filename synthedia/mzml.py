import os, copy, random, math, sys
from pyopenms import *
from .peak_models import *

import time
class MZMLReader():
    def __init__(self, file):
        self.od_exp = OnDiscMSExperiment()
        self.od_exp.openFile(file)
        self.num_spectra = self.od_exp.getNrSpectra()


    def __iter__(self):
        self.current_spec = 0
        return self

    def __next__(self):
        if self.current_spec < self.num_spectra:
            s = self.od_exp.getSpectrum(self.current_spec)
            t = s.getRT() / 60
            lvl = s.getMSLevel()
            mzs, ints = s.get_peaks()
            self.current_spec += 1
            return t, lvl, mzs, ints
        else:
            self.close()
            raise StopIteration

    def close(self):
        del self.od_exp
        return

class Server():

    def __init__(self):
        file = '/home/mleeming/share/synthedia_tests/220329_MB-4007_08.mzML'
        self.od_exp = OnDiscMSExperiment()
        self.od_exp.openFile(file)
        self.num_spectra = self.od_exp.getNrSpectra()
        self.current_spec = 0
        return

    def get_next_ms1(self):

        while True:
            self.current_spec += 1
            s = self.od_exp.getSpectrum(self.current_spec)
            lvl = s.getMSLevel()
            if lvl == 1:
                break
        mzs, ints = s.get_peaks()
        return mzs, ints

    def get_next_ms2(self):

        while True:
            self.current_spec += 1
            s = self.od_exp.getSpectrum(self.current_spec)
            lvl = s.getMSLevel()
            if lvl == 2:
                break
        mzs, ints = s.get_peaks()
        return mzs, ints

class TemplateReader():

    def __init__(self):
        file = '/home/mleeming/share/synthedia_tests/oDIA_HEK200ng_500-740_T60_2.mzML'
        print('started open template')
        self.od_exp = OnDiscMSExperiment()
        self.od_exp.openFile(file)
        self.num_spectra = self.od_exp.getNrSpectra()
        print('finished open template')

    def __iter__(self):
        self.current_spec = 28000
        return self

    def __next__(self):
        if self.current_spec < self.num_spectra:
            s = self.od_exp.getSpectrum(self.current_spec)
            self.current_spec += 1
            return TemplateSpectrum(s)
        else:
            self.close()
            raise StopIteration

    def close(self):
        del self.od_exp
        return


class MZMLWriter():

    def __init__(self, out_file, n_spectra):
        self.n_spectra = n_spectra
        self.consumer = PlainMSDataWritingConsumer(out_file)
        self.consumer.setExpectedSize(n_spectra,0)

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
#        self.server = Server()
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
            print(lvl, self.ms1_template, self.ms2_template)
            if lvl == 1:
                self.ms1_template = s
            if lvl == 2:
                self.ms2_template = s
                print(s.getPrecursors())
        del od_exp
        return

    def get_ghost_data(self, options, spec):

        if spec.order == 1:
            if len(self.ghost_ms1_data) > 0:
                return self.ghost_ms1_data
            else:
                self.ghost_ms1_data = spec.prepare_ghosts(options)
                return self.ghost_ms1_data
        if spec.order == 2:
            if len(self.ghost_ms2_data) > 0:
                return self.ghost_ms2_data
            else:
                self.ghost_ms2_data = spec.prepare_ghosts(options)
                return self.ghost_ms2_data

    def write_spec_new(self, options, spec):

        #spec_to_write = MSSpectrum()
        spec_to_write = spec.spectrum
#        if spec.order == 1:
#            spec_to_write = copy.deepcopy(self.ms1_template)
#        if spec.order == 2:
#            spec_to_write = copy.deepcopy(self.ms2_template)

        indicies = sorted(list(set(spec.indicies)))

#        print(len(indicies))
        ints = spec.ints[indicies]
        mzs = spec.mzs[indicies]

        if options.write_empty_spectra == False:
            if len(ints) > 0:
                spec_to_write.set_peaks([mzs, ints])
            else:
                del spec_to_write
                return
        else:
            spec_to_write.set_peaks([mzs, ints])

        spec_to_write.setType(spec.spectrum.getType())
#        spec_to_write.setType(1)

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

            p.setActivationMethods(set([0]))
            p.setActivationEnergy(30)

            spec_to_write.setPrecursors( [p] )

        # Add scan metadata
        # - some software will crash reading the file without this
#        instrument_settings = InstrumentSettings()
#        instrument_settings.setPolarity(1)

#        scan_window = ScanWindow()

#        string = spec.spectrum.getMetaValue(b'filter string')

#        scan_window.begin = int(float(string.split('[')[1].split('-')[0]))
#        scan_window.end = int(float(string.split('[')[1].split(']')[0].split('-')[1]))

#        instrument_settings.setScanWindows([scan_window])
#        spec_to_write.setInstrumentSettings( instrument_settings)

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

    def write_spec(self, options, spec):

        spec_to_write = MSSpectrum()
#        if spec.order == 1:
#            bmz, bint = self.server.get_next_ms1()
#        else:
#            bmz, bint = self.server.get_next_ms2()
        bmz, bint = spec.bmz, spec.bint

#        if spec.order == 1:
#            spec_to_write = copy.deepcopy(self.ms1_template)
#        if spec.order == 2:
#            spec_to_write = copy.deepcopy(self.ms2_template)

#        if len(spec.indicies) == 0:
#            n_ghosts = 20
#            ghost_data = self.get_ghost_data(options, spec)
#            spec.make_ghost_peak(options, n_ghosts, ghost_data)

        indicies = sorted(list(set(spec.indicies)))

#        print(len(indicies))
        ints = spec.ints[indicies]
        mzs = spec.mzs[indicies]
#        ints = spec.ints
#        mzs = spec.mzs

        ints = np.append(ints, bint)
        mzs = np.append(mzs, bmz)
#
        mask = np.argsort(mzs)
#
        mzs = mzs[mask]
        ints = ints[mask]

#        mask = np.where(ints < 1)
#        ints = ints[mask]
#        mzs = mzs[mask]

        if options.write_empty_spectra == False:
            if len(ints) > 0:
                spec_to_write.set_peaks([mzs, ints])
            else:
                del spec_to_write
                return
        else:
            spec_to_write.set_peaks([mzs, ints])

        if options.centroid:
            centroided_spectrum = MSSpectrum()
            PeakPickerHiRes(
                min_required_elements = 1,
                wrote_log_messages = 'true'
            ).pick(spec_to_write, centroided_spectrum)
            spec_to_write = centroided_spectrum
            spec_to_write.setType(1)
        else:
            spec_to_write.setType(2)

#        spec_to_write.setRT(spec.rt)
        spec_to_write.setRT(spec.rt)


        spec_to_write.setMSLevel(spec.order)

        if spec.order == 2:
            p = Precursor()
            p.setMZ((spec.isolation_hl + spec.isolation_ll) / 2)
            p.setIsolationWindowLowerOffset(spec.lower_offset)
            p.setIsolationWindowUpperOffset(spec.upper_offset)

            p.setActivationMethods(set([0]))
            p.setActivationEnergy(30)

            spec_to_write.setPrecursors( [p] )

        # Add scan metadata
        # - some software will crash reading the file without this
        instrument_settings = InstrumentSettings()
        instrument_settings.setPolarity(1)

        scan_window = ScanWindow()

        if spec.order == 1:
            scan_window.begin = options.ms1_min_mz
            scan_window.end = options.ms1_max_mz
        elif spec.order == 2:
            scan_window.begin = options.ms2_min_mz
            scan_window.end = options.ms2_max_mz

        instrument_settings.setScanWindows([scan_window])
        spec_to_write.setInstrumentSettings( instrument_settings)

        spec_to_write.setMetaValue('lowest observed m/z', mzs[0])
        spec_to_write.setMetaValue('highest observed m/z', mzs[-1])
        spec_to_write.setMetaValue('base peak intensity', float(max(ints)))
        spec_to_write.setMetaValue('base peak m/z', mzs[np.argmax(ints)])
        spec_to_write.setMetaValue('total ion current', float(sum(ints)))

#        spec_to_write.setMetaValue('preset scan configuration', 0)

        # MS1 filter string template
        # "FTMS + p NSI Full ms [350.0000-1400.0000]"
        # MS2 filter string template
        # "FTMS + c NSI Full ms2 422.6500@hcd30.00 [200.0000-2000.0000]"
        if spec.order == 1:
            spec_to_write.setMetaValue('filter string',
                "FTMS + p NSI Full ms [%s.0000-%s.0000]" %(int(options.ms1_min_mz), int(options.ms1_max_mz))
            )
        elif spec.order == 2:
            spec_to_write.setMetaValue('filter string',
                "FTMS + c NSI Full ms2 %.4f@cid30.00 [%s.0000-%s.0000]" %(
                    (spec.isolation_hl + spec.isolation_ll) / 2, int(options.ms2_min_mz), int(options.ms2_max_mz)
                )
            )

        if spec.order == 2:
            if len(spec_to_write.get_peaks()[0]) == 0:
                print(mzs)
                print('delete')
                del spec_to_write
                return

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
    def __init__(self, rt, order, isolation_range):
        self.rt = rt
        self.order = order

        if isolation_range:
            self.isolation_ll = isolation_range[0]
            self.isolation_hl = isolation_range[1]

            self.lower_offset = (isolation_range[1] - isolation_range[0]) / 2
            self.upper_offset = (isolation_range[1] - isolation_range[0]) / 2

        return

    def make_spectrum(self, MS1_MZS, MS1_INTS, MS2_MZS, MS2_INTS):
        # much faster than creating a new array for every scan
        if self.order == 1:
            self.mzs = MS1_MZS
            self.ints = np.zeros(len(MS1_INTS))
        else:
            self.mzs = MS2_MZS
            self.ints = np.zeros(len(MS2_INTS))

        self.indicies = []
        return

    def add_peaks(self, options, p, groupi, samplei):

        peptide_scaled_rt = p.scaled_rt
        abundance_offset = p.offsets[groupi][samplei]

        # scaling factor for point on chromatogram
        intensity_scale_factor = options.rt_peak_model(self.rt, **{
            'mu': peptide_scaled_rt, 'sig': options.rt_stdev, 'emg_k': options.rt_emg_k
        })

        # apply chromatographic instability if needed
        if options.rt_instability > 0:
            instability_factor = 1 - ( random.randint(0,options.rt_instability * 100) / 10000 )
            intensity_scale_factor *= instability_factor

        if self.order == 1:
            peaks = p.ms1_isotopes
            stdev = options.ms1_stdev
            min_peak_intensity = options.ms1_min_peak_intensity
        else:
            peaks = p.ms2_peaks
            stdev = options.ms2_stdev
            min_peak_intensity = options.ms2_min_peak_intensity

        for peaki, peak in enumerate(peaks):

            lower_limit, higher_limit, indicies = peak.get_limits(options, self.mzs)
            self.indicies.extend(indicies)

            # apply offset to peak intensity
            log_int = math.log2(peak.intensity)
            adjusted_log2_int = log_int + abundance_offset
            adjusetd_raw_int = 2 ** adjusted_log2_int

            if not hasattr(peak, 'peak_intensities'):

                # calculating the peak for the full m/z range is slow
                # subset data to only a small region around the peak to speed calculation
                peak_ints = options.mz_peak_model(self.mzs[lower_limit:higher_limit], **{
                    'mu': peak.mz, 'sig': stdev, 'emg_k': options.mz_emg_k
                })

#                if self.order == 2:
#                    for fragment_isotope_i in range(1, 3):
#                        isotope_offset = fragment_isotope_i * 1.007276
#                        peak_ints += options.mz_peak_model(self.mzs[lower_limit:higher_limit], **{
#                            'mu': peak.mz + isotope_offset, 'sig': stdev, 'emg_k': options.mz_emg_k
#                        }) * 1/ (fragment_isotope_i + 1) ** 2
                peak.set_peak_intensities(options, peak_ints)

            peak_ints = peak.get_peak_intensities()

            # scale peak intensities by chromatogram scaling factor
            factor = adjusetd_raw_int * intensity_scale_factor
            peak_ints *= factor

            # in the case of MS1 spectra - increment only for one isotope
            # in the case of MS2 spectra - increment only for one fragment
            if peaki == 0:
                # track the number of data points per peak
                if max(peak_ints) > 0:
                    p.increment_points_per_peak_dict(groupi, samplei, self.order)

            # remove low intensity points
            int_mask = np.where(peak_ints < min_peak_intensity)
            peak_ints[int_mask] = 0

            # track peak intensities
            peak.update_intensities_to_report(groupi, samplei, peak_ints)

            # add new data to full spectrum intensity
            self.ints[lower_limit:higher_limit] += peak_ints

        return

    def make_ghost_peak(self, options, n_ghosts, ghost_data):

        ghost_intensities = ghost_data['ghost_intensities']
        ghost_indicies = ghost_data['ghost_indicies']

        ghost_max_int = 100
        for i in range(n_ghosts):
            index = random.randint(0,len(self.mzs) - len(ghost_indicies) - 2)
            mask = ghost_indicies + index
            self.ints[mask] = ghost_intensities * random.randint(10, 500000)
            self.indicies.extend(mask)

        return

    def prepare_ghosts(self, options):
        if self.order == 1:
            stdev = options.ms1_stdev
        else:
            stdev = options.ms2_stdev
        peak_ints = options.mz_peak_model(self.mzs, **{
            'mu': 500, 'sig': stdev, 'emg_k': options.mz_emg_k
        })
        ghost_indicies = np.where( peak_ints > 0.001)
        ghost_intensities = peak_ints[ghost_indicies]
        ghost_intensities = ghost_intensities / max(ghost_intensities)
        ghost_indicies = ghost_indicies[0] - min(ghost_indicies[0])
        return {'ghost_intensities': ghost_intensities, 'ghost_indicies': ghost_indicies}

    def clear(self):
        # save memory
        del self.mzs
        del self.ints
        del self.indicies
        return

class TemplateSpectrum():
    def __init__(self, s):

        self.spectrum = s
        self.rt = s.getRT()
        self.order = s.getMSLevel()
        self.bmz, self.bint = s.get_peaks()
        self.indicies = []
        if self.order == 2:
            self.lower_offset = s.getPrecursors()[0].getIsolationWindowLowerOffset()
            self.upper_offset = s.getPrecursors()[0].getIsolationWindowUpperOffset()
            self.isolation_ll = s.getPrecursors()[0].getMZ() - self.lower_offset
            self.isolation_hl = s.getPrecursors()[0].getMZ() + self.upper_offset

        return

    def make_spectrum(self, MS1_MZS, MS1_INTS, MS2_MZS, MS2_INTS):
        # much faster than creating a new array for every scan
        if self.order == 1:
            self.mzs = MS1_MZS
            self.ints = np.zeros(len(MS1_INTS))
        else:
            self.mzs = MS2_MZS
            self.ints = np.zeros(len(MS2_INTS))

        self.indicies = []
        return

    def add_peaks(self, options, p, groupi, samplei):

        peptide_scaled_rt = p.scaled_rt
        abundance_offset = p.offsets[groupi][samplei]

        # scaling factor for point on chromatogram
        intensity_scale_factor = options.rt_peak_model(self.rt, **{
            'mu': peptide_scaled_rt, 'sig': options.rt_stdev, 'emg_k': options.rt_emg_k
        })

        # apply chromatographic instability if needed
        if options.rt_instability > 0:
            instability_factor = 1 - ( random.randint(0,options.rt_instability * 100) / 10000 )
            intensity_scale_factor *= instability_factor

        if self.order == 1:
            peaks = p.ms1_isotopes
            stdev = options.ms1_stdev
            min_peak_intensity = options.ms1_min_peak_intensity
        else:
            peaks = p.ms2_peaks
            stdev = options.ms2_stdev
            min_peak_intensity = options.ms2_min_peak_intensity

        for peaki, peak in enumerate(peaks):

            lower_limit, higher_limit, indicies = peak.get_limits(options, self.mzs)
            self.indicies.extend(indicies)

            # apply offset to peak intensity
            log_int = math.log2(peak.intensity)
            adjusted_log2_int = log_int + abundance_offset
            adjusetd_raw_int = 2 ** adjusted_log2_int

            if not hasattr(peak, 'peak_intensities'):

                # calculating the peak for the full m/z range is slow
                # subset data to only a small region around the peak to speed calculation
                peak_ints = options.mz_peak_model(self.mzs[lower_limit:higher_limit], **{
                    'mu': peak.mz, 'sig': stdev, 'emg_k': options.mz_emg_k
                })

#                if self.order == 2:
#                    for fragment_isotope_i in range(1, 3):
#                        isotope_offset = fragment_isotope_i * 1.007276
#                        peak_ints += options.mz_peak_model(self.mzs[lower_limit:higher_limit], **{
#                            'mu': peak.mz + isotope_offset, 'sig': stdev, 'emg_k': options.mz_emg_k
#                        }) * 1/ (fragment_isotope_i + 1) ** 2
                peak.set_peak_intensities(options, peak_ints)

            peak_ints = peak.get_peak_intensities()

            # scale peak intensities by chromatogram scaling factor
            factor = adjusetd_raw_int * intensity_scale_factor
            peak_ints *= factor

            # in the case of MS1 spectra - increment only for one isotope
            # in the case of MS2 spectra - increment only for one fragment
            if peaki == 0:
                # track the number of data points per peak
                if max(peak_ints) > 0:
                    p.increment_points_per_peak_dict(groupi, samplei, self.order)

            # remove low intensity points
            int_mask = np.where(peak_ints < min_peak_intensity)
            peak_ints[int_mask] = 0

            # track peak intensities
            peak.update_intensities_to_report(groupi, samplei, peak_ints)

            # add new data to full spectrum intensity
            self.ints[lower_limit:higher_limit] += peak_ints

        return

    def make_ghost_peak(self, options, n_ghosts, ghost_data):

        ghost_intensities = ghost_data['ghost_intensities']
        ghost_indicies = ghost_data['ghost_indicies']

        ghost_max_int = 100
        for i in range(n_ghosts):
            index = random.randint(0,len(self.mzs) - len(ghost_indicies) - 2)
            mask = ghost_indicies + index
            self.ints[mask] = ghost_intensities * random.randint(10, 500000)
            self.indicies.extend(mask)

        return

    def prepare_ghosts(self, options):
        if self.order == 1:
            stdev = options.ms1_stdev
        else:
            stdev = options.ms2_stdev
        peak_ints = options.mz_peak_model(self.mzs, **{
            'mu': 500, 'sig': stdev, 'emg_k': options.mz_emg_k
        })
        ghost_indicies = np.where( peak_ints > 0.001)
        ghost_intensities = peak_ints[ghost_indicies]
        ghost_intensities = ghost_intensities / max(ghost_intensities)
        ghost_indicies = ghost_indicies[0] - min(ghost_indicies[0])
        return {'ghost_intensities': ghost_intensities, 'ghost_indicies': ghost_indicies}

    def clear(self):
        # save memory
        del self.mzs
        del self.ints
        del self.indicies
        return
