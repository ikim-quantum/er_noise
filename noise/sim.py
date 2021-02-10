# Author : Isaac H. Kim
# Last updated: 12/23/2020
from compressor import qubits
from scipy import sparse
import numpy as np
import random


# Permutes the subsystems according to perm
def syspermute(rho, perm):
    n=len(rho)
    qubits = int(np.log2(n))
    perm = [qubits-i-1 for i in reversed(perm)]
    perm = perm + [qubits+i for i in perm]
    return np.transpose(rho.reshape([2] * qubits *2, order='F'),perm).reshape([n,n], order='F')


# Applies a two-qubit gate
def twoQ(U,rho,control,target):
    dim = len(rho)
    perm = [n for n in range(int(np.log2(dim)))]
    perm.remove(control)
    perm.remove(target)
    perm = [control] + [target] + perm
    iperm=list(np.argsort(list(perm)))
    rho_permuted = syspermute(rho, iperm)
    U_e = sparse.kron(U, sparse.eye(dim/4))
    rho_pc = U_e @ rho_permuted @ U_e.getH()
    return syspermute(rho_pc,perm)


# Creates a random n x n unitary matrix.
def randU(n):
    X = (np.random.randn(n,n) + 1j * np.random.randn(n,n))/np.sqrt(2)
    Q, R = np.linalg.qr(X)
    R = np.diag(np.diag(R) / abs(np.diag(R)))
    return Q @ R


# Apply a depolarizing noise on the indexed qubit with noise rate p.
def depolarizing(rho,index,p):
    dim = len(rho)
    dim_former = pow(2,index)
    dim_latter = dim / 2 / dim_former

    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])

    X_e =  sparse.kron(sparse.kron(sparse.eye(dim_former), X), sparse.eye(dim_latter))
    Y_e =  sparse.kron(sparse.kron(sparse.eye(dim_former), Y), sparse.eye(dim_latter))
    Z_e =  sparse.kron(sparse.kron(sparse.eye(dim_former), Z), sparse.eye(dim_latter))

    return (1-p) * rho + (p/3) * (X_e @ rho @ X_e + Y_e @ rho @ Y_e + Z_e @ rho @ Z_e)


def amp_damping(rho,index,p):
    dim = len(rho)
    dim_former = pow(2,index)
    dim_latter = dim / 2 / dim_former
    K0 = np.array([[1,0], [0,np.sqrt(1-p)]])
    K1 = np.array([[0,np.sqrt(p)], [0, 0]])
    return K0 @ rho @ K0.T.conj() + K1 @ rho @ K1.T.conj()


# Apply partial trace on the indexed qubit
def partial_trace(rho, index):
    dim = len(rho)
    dim_former = pow(2,index)
    dim_latter = dim / 2 / dim_former
    
    unitvectors=np.identity(2)
    v1 = sparse.kron(sparse.kron(sparse.eye(dim_former), unitvectors[0,:]), sparse.eye(dim_latter))
    v2 = sparse.kron(sparse.kron(sparse.eye(dim_former), unitvectors[1,:]), sparse.eye(dim_latter))
    return v1 @ rho @ v1.T.conj() + v2 @ rho @ v2.T.conj() 


# Reset the indexed qubit
def reset(rho, index):
    dim = len(rho)
    dim_former = pow(2,index)
    dim_latter = dim / 2 / dim_former

    unitvectors = np.identity(2)
    p1 = sparse.kron(sparse.kron(sparse.eye(dim_former), np.outer(unitvectors[0,:].T, unitvectors[0,:])), sparse.eye(dim_latter))
    p2 = sparse.kron(sparse.kron(sparse.eye(dim_former), np.outer(unitvectors[0,:].T, unitvectors[1,:])), sparse.eye(dim_latter))
    return p1 @ rho @ p1.T + p2 @ rho @ p2.T


def sim(circ, verbose=False):
    """
    Simulate circuit
    Input format:
    circuit is a list of depth-1 circuits, which themselves are lists.
    If the element is a single integer, that's a reset. Otherwise, we have
    a set of gates. We will apply a random circuit.
    """
    n_q = len(qubits(circ))
    dim = 2**n_q
    tempvec = np.zeros((dim,1))
    tempvec[0] = 1.0
    rho = np.outer(tempvec, tempvec)
    for c in circ:
        for g in c:
            if len(g)==1:
                if verbose:
                    print('Reset {}'.format(g[0]))
                rho = reset(rho, g[0])
            else:
                U = randU(4)
                rho = twoQ(U, rho, g[0], g[1])
    return rho


def noise_study(circ, p, verbose=False):
    qs = qubits(circ)
    n_q = len(qs)
    dim = 2**n_q
    tempvec = np.zeros((dim,1))
    tempvec[0] = 1.0
    rho = np.outer(tempvec, tempvec)
    rho_e = np.outer(tempvec, tempvec)
    for q in qs:
        rho_e = depolarizing(rho_e, q, p)

    for c in circ:
        for g in c:
            if len(g)==1:
                if verbose:
                    print('Reset {}'.format(g[0]))
                rho = reset(rho, g[0])
                rho_e = depolarizing(reset(rho_e, g[0]), g[0], p)
            else:
                if verbose:
                    print('Gate {}'.format(g))
                U = randU(4)
                rho = twoQ(U, rho, g[0], g[1])
                rho_e = depolarizing(depolarizing(twoQ(U,rho_e,g[0],g[1]),g[0],p),g[1],p)

    return rho_e - rho


def extrapolation_study(circ, p, verbose=False):
    qs = qubits(circ)
    n_q = len(qs)
    dim = 2**n_q
    tempvec = np.zeros((dim,1))
    tempvec[0] = 1.0
    rho = np.outer(tempvec, tempvec)
    rho_e = np.outer(tempvec, tempvec)
    rho_2e = np.outer(tempvec, tempvec)
    for q in qs:
        rho_e = depolarizing(rho_e, q, p)
        rho_2e = depolarizing(rho_2e, q, 2*p)

    for c in circ:
        for g in c:
            if len(g)==1:
                if verbose:
                    print('Reset {}'.format(g[0]))
                rho = reset(rho, g[0])
                rho_e = depolarizing(reset(rho_e, g[0]), g[0], p)
                rho_2e = depolarizing(reset(rho_2e, g[0]), g[0], 2*p)
            else:
                if verbose:
                    print('Gate {}'.format(g))
                U = randU(4)
                rho = twoQ(U, rho, g[0], g[1])
                rho_e = depolarizing(depolarizing(twoQ(U,rho_e,g[0],g[1]),g[0],p),g[1],p)
                rho_2e = depolarizing(depolarizing(twoQ(U,rho_2e,g[0],g[1]),g[0],2*p),g[1],2*p)
                
    return rho + rho_2e - 2*rho_e
