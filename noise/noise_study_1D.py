from sim import noise_study, extrapolation_study
from D1_v2 import pcc_dmera_1d
from compressor import compress
import numpy as np


def noise_estimate_pcc(n, D, p, samples, verbose=False):
    results = []
    for i in range(samples):
        if verbose:
            print("{}/{}".format(i+1,samples))
        randsite = np.random.randint(2**n)
        if verbose:
            print("randsite = {}".format(randsite))
        circ = compress(pcc_dmera_1d(n, D, [randsite, (randsite+1)%(2**n)]))
        results.append(np.sum(abs(np.linalg.eigvals(noise_study(circ, p)))))
    return results


def extrapolation_estimate_pcc(n, D, p, samples, verbose=False):
    results = []
    for i in range(samples):
        if verbose:
            print("{}/{}".format(i+1,samples))
        randsite = np.random.randint(2**n)
        if verbose:
            print("randsite = {}".format(randsite))
        circ = compress(pcc_dmera_1d(n, D, [randsite, (randsite+1)%(2**n)]))
        results.append(np.sum(abs(np.linalg.eigvals(extrapolation_study(circ, p)))))
    return results


def noise_estimate_pcc_naive(n, D, p, samples, verbose=False):
    results = []
    for i in range(samples):
        if verbose:
            print("{}/{}".format(i+1,samples))
        randsite = np.random.randint(2**n)
        if verbose:
            print("randsite = {}".format(randsite))
        circ = compress(pcc_dmera_1d(n, D, [randsite, (randsite+1)%(2**n)]))
        n_e = 0
        for c in circ:
            for g in c:
                n_e += len(g)
        results.append(n_e * p)
    return results
