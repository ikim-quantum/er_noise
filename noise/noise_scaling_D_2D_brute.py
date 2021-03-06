from noise_study_2D_brute import noise_estimate_pcc, noise_estimate_pcc_naive, extrapolation_estimate_pcc
import numpy as np


ns = [2]
p = 0.001
Ds = [1, 2, 3, 4, 5]

for n in ns:
    print("n={}".format(n))
    for D in Ds:
        print("D={}".format(D))
        results = noise_estimate_pcc(n, D, p, verbose=True)
        results_ext = extrapolation_estimate_pcc(n, D, p, verbose=False)
        results_naive = noise_estimate_pcc_naive(n, D, p)
        print("naive = {}+-{}".format(np.mean(results_naive), np.std(results_naive)))
        print("precise = {}+-{}".format(np.mean(results), np.std(results)))
        print("Extrapolated = {}+-{}".format(np.mean(results_ext), np.std(results_ext)))
