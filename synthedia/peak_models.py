import sys
import numpy as np
from scipy.stats import exponnorm, cauchy, norm

class InvalidPeakModelError(Exception):
    pass

class PeakModels(object):

    @staticmethod
    def gaussian(x, **kwargs):
        if x is not None:
            return norm.pdf(x, kwargs['mu'], kwargs['sig'])
        else:
            return norm.rvs(kwargs['mu'], kwargs['sig'])


    @staticmethod
    def exponentially_modified_gaussian(x, **kwargs):
        if kwargs['emg_k'] <= 0:
            return PeakModels.gaussian(x, **kwargs)

        if x is not None:
            return exponnorm.pdf(x, kwargs['emg_k'], kwargs['mu'], kwargs['sig'])
        else:
            return exponnorm.rvs(kwargs['emg_k'], kwargs['mu'], kwargs['sig'])

    @staticmethod
    def cauchy(x, **kwargs):
        if x is not None:
            return cauchy.pdf(x, kwargs['mu'], kwargs['sig'])
        else:
            return cauchy.rvs(kwargs['mu'], kwargs['sig'])

    def get_rt_peak_model(self, options):
        try:
            rt_peak_model = getattr(self, options.rt_peak_model)
        except AttributeError:
            peak_models = []
            for method in dir(self):
                if ('__' in method) or ('get_' in method): continue
                peak_models.append(method)
            raise InvalidPeakModelError(
                'RT peak model selection %s is invalid. Valid options are: %s' % (
                    options.rt_peak_model, ', '.join(peak_models))
            )
        return rt_peak_model

    def get_mz_peak_model(self, options):
        try:
            mz_peak_model = getattr(self, options.mz_peak_model)
        except AttributeError:
            peak_models = []
            for method in dir(self):
                if ('__' in method) or ('get_' in method): continue
                peak_models.append(method)
            raise InvalidPeakModelError(
                'MZ peak model selection %s is invalid. Valid options are: %s' % (
                    options.mz_peak_model, ', '.join(peak_models))
            )
        return mz_peak_model

    def get_rt_peak_fwhm_distribution_model(self, options):
        try:
            rt_peak_fwhm_distribution_model = getattr(self, options.rt_peak_fwhm_distribution_model)
        except AttributeError:
            peak_models = []
            for method in dir(self):
                if ('__' in method) or ('get_' in method): continue
                peak_models.append(method)
            raise InvalidPeakModelError(
                'RT peak FWHM distribution model selection %s is invalid. Valid options are: %s' % (
                    options.rt_peak_fwhm_distribution_model, ', '.join(peak_models))
            )
        return rt_peak_fwhm_distribution_model

    def get_prosit_peptide_abundance_model(self, options):
        try:
            prosit_abundance_model = getattr(self, options.prosit_peptide_abundance_model)
        except AttributeError:
            peak_models = []
            for method in dir(self):
                if ('__' in method) or ('get_' in method): continue
                peak_models.append(method)
            raise InvalidPeakModelError(
                'Prosit abundance model selection %s is invalid. Valid options are: %s' % (
                    options.prosit_peptide_abundance_model, ', '.join(peak_models))
            )
        return prosit_abundance_model

