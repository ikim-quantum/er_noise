from noise_study_1D import noise_estimate_pcc, noise_estimate_pcc_naive, extrapolation_estimate_pcc
import numpy as np


n = 8
ps = [0.0001*2**i for i in range(10)]
D = 5
samples = 500

for p in ps:
    print("p={}".format(p))
    results = noise_estimate_pcc(n, D, p, samples, verbose=True)
    results_ext = extrapolation_estimate_pcc(n, D, p, samples, verbose=True)
    results_naive = noise_estimate_pcc_naive(n, D, p, samples)
    print("naive = {}+-{}".format(np.mean(results_naive), np.std(results_naive)))
    print("precise = {}+-{}".format(np.mean(results), np.std(results)))
    print("Extrapolated = {}+-{}".format(np.mean(results_ext), np.std(results_ext)))
