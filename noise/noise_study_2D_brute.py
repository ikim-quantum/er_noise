from sim import noise_study, extrapolation_study
from D2_v2 import pcc_dmera_2d, convert_2d_to_1d
from compressor import compress, qubits
import numpy as np


def noise_estimate_pcc(n, D, p, verbose=False):
    results = []
    for x in range(2**n):
        for y in range(2**n):
            if verbose:
                print("{},{}".format(x, y))
        
            supp_x = [(x, y), ((x+1)%(2**n), y)]
            supp_y = [(x, y), (x, (y+1)%(2**n))]
            
            supp_x_con = [convert_2d_to_1d(c,n) for c in supp_x]
            supp_y_con = [convert_2d_to_1d(c,n) for c in supp_y]
        
            circ = compress(pcc_dmera_2d(n, D, supp_x_con))
            results.append(np.sum(abs(np.linalg.eigvals(noise_study(circ, p)))))
            n_q = len(qubits(circ))

            circ = compress(pcc_dmera_2d(n, D, supp_y_con))
            results.append(np.sum(abs(np.linalg.eigvals(noise_study(circ, p)))))
            n_q = len(qubits(circ))
        
            if verbose:
                print("Number of qubits={}".format(n_q))

    return results


def extrapolation_estimate_pcc(n, D, p, verbose=False):
    results = []
    for x in range(2**n):
        for y in range(2**n):
            if verbose:
                print("{},{}".format(x, y))

        
            supp_x = [(x, y), ((x+1)%(2**n), y)]
            supp_y = [(x, y), (x, (y+1)%(2**n))]
            
            supp_x_con = [convert_2d_to_1d(c,n) for c in supp_x]
            supp_y_con = [convert_2d_to_1d(c,n) for c in supp_y]
        
            circ = compress(pcc_dmera_2d(n, D, supp_x_con))
            results.append(np.sum(abs(np.linalg.eigvals(extrapolation_study(circ, p)))))
            n_q = len(qubits(circ))

            circ = compress(pcc_dmera_2d(n, D, supp_y_con))
            results.append(np.sum(abs(np.linalg.eigvals(extrapolation_study(circ, p)))))
            n_q = len(qubits(circ))
        
            if verbose:
                print("Number of qubits={}".format(n_q))
    return results


def noise_estimate_pcc_naive(n, D, p, verbose=False):
    results = []
    
    for x in range(2**n):
        for y in range(2**n):
            if verbose:
                print("{},{}".format(x, y))

        
            supp_x = [(x, y), ((x+1)%(2**n), y)]
            supp_y = [(x, y), (x, (y+1)%(2**n))]
            
            supp_x_con = [convert_2d_to_1d(c,n) for c in supp_x]
            supp_y_con = [convert_2d_to_1d(c,n) for c in supp_y]

            circ = compress(pcc_dmera_2d(n, D, supp_x_con))
            n_e = 0
            for c in circ:
                for g in c:
                    n_e += len(g)
            results.append(n_e * p)
            
            circ = compress(pcc_dmera_2d(n, D, supp_y_con))
            n_e = 0
            for c in circ:
                for g in c:
                    n_e += len(g)
            results.append(n_e * p)

    return results
