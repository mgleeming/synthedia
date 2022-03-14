import os, copy, random, math
from pyopenms import *
from numba import jit

@jit(nopython=True)
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

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

class MZMLWriter():

    def __init__(self, out_file):
        self.consumer = PlainMSDataWritingConsumer(out_file)
        self.n_spec_written = 0
        return

    def write_spec(self, options, spec):
        spec_to_write = MSSpectrum()

        mask = np.where(spec.ints > 0)
        ints = spec.ints[mask]
        mzs = spec.mzs[mask]

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
            PeakPickerHiRes().pick(spec_to_write, centroided_spectrum)
            spec_to_write = centroided_spectrum

        spec_to_write.setRT(spec.rt)
        spec_to_write.setMSLevel(spec.order)

        if spec.order == 2:
            p = Precursor()
            p.setMZ((spec.isolation_hl + spec.isolation_ll) / 2)
            p.setIsolationWindowLowerOffset(spec.isolation_ll)
            p.setIsolationWindowUpperOffset(spec.isolation_hl)
            spec_to_write.setPrecursors( [p] )

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

        return

    def make_spectrum(self, MS1_MZS, MS1_INTS, MS2_MZS, MS2_INTS):
        # much faster than creating a new array for every scan
        if self.order == 1:
            self.mzs = MS1_MZS
            self.ints = copy.deepcopy(MS1_INTS)
        else:
            self.mzs = MS2_MZS
            self.ints = copy.deepcopy(MS2_INTS)
        return

    def add_peaks(self, options, peptide_scaled_rt, peaks, abundance_offset):

        # scaling factor for point on chromatogram
        intensity_scale_factor = gaussian(self.rt, peptide_scaled_rt, options.rt_stdev)

        # apply chromatographic instability if needed
        if options.rt_instability > 0:
            instability_factor = 1 - ( random.randint(0,options.rt_instability * 100) / 10000 )
            intensity_scale_factor *= instability_factor

        # make sure peak decays to 0
        if intensity_scale_factor < options.min_peak_fraction: return

        if self.order == 1:
            stdev = options.ms1_stdev
        else:
            stdev = options.ms2_stdev

        for peak in peaks:

            # apply offset to peak intensity
            log_int = math.log2(peak[1])
            adjusted_log2_int = log_int + abundance_offset
            adjusetd_raw_int = 2 ** adjusted_log2_int

            # calculating the gaussian for the full m/z range is slow
            # subset data to only a small region around the peak to speed calculation
            mz_mask = np.where((self.mzs > peak[0] - options.ms_clip_window) & (self.mzs < peak[0] + options.ms_clip_window))
            peak_ints = gaussian(self.mzs[mz_mask], peak[0], stdev)

            # scale gaussian intensities by chromatogram scaling factor
            factor = adjusetd_raw_int * intensity_scale_factor
            peak_ints *= factor

            # remove low intensity points
            int_mask = np.where(peak_ints < options.min_peak_fraction)
            peak_ints[int_mask] = 0

            # add new data to full spectrum intensity
            self.ints[mz_mask] += peak_ints

        return

    def clear(self):
        # save memory
        del self.mzs
        del self.ints
        return


