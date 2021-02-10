from noise_study import noise_estimate_pcc, noise_estimate_pcc_naive
import numpy as np


n = 8
p = 0.001
D = 2
samples = 1000

results = noise_estimate_pcc(n, D, p, samples)
results_naive = noise_estimate_pcc_naive(n, D, p, samples)
print("precise = {}+-{}".format(np.mean(results), np.std(results)))
print("naive = {}+-{}".format(np.mean(results_naive), np.std(results_naive)))
