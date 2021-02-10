from sim import noise_study, extrapolation_study
from D2_v2 import pcc_dmera_2d, convert_2d_to_1d
from compressor import compress, qubits
import numpy as np


def noise_estimate_pcc(n, D, p, samples, verbose=False):
    results = []
    for i in range(samples):
        if verbose:
            print("{}/{}".format(i+1,samples))

        randx = np.random.randint(2**n)
        randy = np.random.randint(2**n)
        supp = []
        
        if np.random.rand()<0.5:
            supp = [(randx, randy), ((randx+1)%(2**n), randy)]
        else:
            supp = [(randx, randy), (randx, (randy+1)%(2**n))]
        supp_con = [convert_2d_to_1d(c,n) for c in supp]
        
        if verbose:
            print("randsite = ({},{})".format(randx, randy))
        
        circ = compress(pcc_dmera_2d(n, D, supp_con))
        n_q = len(qubits(circ))
        
        if verbose:
            print("Number of qubits={}".format(n_q))
        
        results.append(np.sum(abs(np.linalg.eigvals(noise_study(circ, p)))))
    return results


def extrapolation_estimate_pcc(n, D, p, samples, verbose=False):
    results = []
    for i in range(samples):
        if verbose:
            print("{}/{}".format(i+1,samples))

        randx = np.random.randint(2**n)
        randy = np.random.randint(2**n)
        supp = []
        
        if np.random.rand()<0.5:
            supp = [(randx, randy), ((randx+1)%(2**n), randy)]
        else:
            supp = [(randx, randy), (randx, (randy+1)%(2**n))]
        supp_con = [convert_2d_to_1d(c,n) for c in supp]
        
        if verbose:
            print("randsite = {}".format(randsite))
        
        circ = compress(pcc_dmera_2d(n, D, supp_con))
        results.append(np.sum(abs(np.linalg.eigvals(extrapolation_study(circ, p)))))
    return results


def noise_estimate_pcc_naive(n, D, p, samples, verbose=False):
    results = []
    for i in range(samples):
        if verbose:
            print("{}/{}".format(i+1,samples))
        
        randx = np.random.randint(2**n)
        randy = np.random.randint(2**n)
        supp = []
        
        if np.random.rand()<0.5:
            supp = [(randx, randy), ((randx+1)%(2**n), randy)]
        else:
            supp = [(randx, randy), (randx, (randy+1)%(2**n))]
        supp_con = [convert_2d_to_1d(c,n) for c in supp]
    
        if verbose:
            print("randsite = {}".format(randsite))

        circ = compress(pcc_dmera_2d(n, D, supp_con))
        n_e = 0
        for c in circ:
            for g in c:
                n_e += len(g)
        results.append(n_e * p)
    return results
