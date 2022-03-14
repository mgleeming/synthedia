import numpy as np
from numba import jit
from scipy.stats import exponnorm, cauchy

@jit(nopython=True)
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def emg(x, k, mu, sig):
    return exponnorm.pdf(x,k,mu,sig)

def cauchy(x, mu, sig):
    return cauchy.pdf(x,mu,sig)
