import sys
import numpy as np
from scipy.stats import exponnorm, cauchy, norm

class PeakModels(object):

    @staticmethod
    def gaussian(x, **kwargs):
        return norm.pdf(x, kwargs['mu'], kwargs['sig'])

    @staticmethod
    def exponentially_modified_gaussian(x, **kwargs):
        if kwargs['emg_k'] <= 0:
            return gaussian(**kwargs)
        return exponnorm.pdf(x, kwargs['emg_k'], kwargs['mu'], kwargs['sig'])

    @staticmethod
    def cauchy(x, **kwargs):
        return cauchy.pdf(x, kwargs['mu'], kwargs['sig'])

    def get_rt_peak_model(self, options):

        try:
            rt_peak_model = getattr(self, options.rt_peak_model)
        except AttributeError:
            print('RT peak model selection %s is invalid. Valid options are:' %options.rt_peak_model)
            for method in dir(self):
                if ('__' in method) or ('get_' in method): continue
                print('\t%s'%method)
            print('Exiting.')
            sys.exit()

        return rt_peak_model

    def get_mz_peak_model(self, options):

        try:
            mz_peak_model = getattr(self, options.mz_peak_model)
        except AttributeError:
            print('MZ peak model selection %s is invalid. Valid options are:' %options.mz_peak_model)
            for method in dir(self):
                if ('__' in method) or ('get_' in method): continue
                print('\t%s'%method)
            print('Exiting.')
            sys.exit()

        return mz_peak_model


