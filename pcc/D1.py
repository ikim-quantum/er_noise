###################################
# Past causal cone calculator: 1D #
# Isaac H. Kim 2020/12/24         #
###################################
import numpy as np
from compressor import pcc, qubits, circ_subtract, rearrange, compress, optimal_width_freeze, depth_width_reduced_freeze


def sites_1d(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(int): A list of sites that appear up to the s'th scale.
    """
    return [i*(1<<(n-s)) for i in range(2**s)]


def nn_e(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(tuple): A list of nearest neighbors at scale s. (even)
    """
    mysites = sites_1d(n, s)
    return [(mysites[2*i], mysites[(2*i+1)%(2**s)]) for i in range(2**(s-1))]


def nn_o(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(tuple): A list of nearest neighbors at scale s. (even)
    """
    mysites = sites_1d(n, s)
    return [(mysites[(2*i+1)%(2**s)], mysites[(2*i+2)%(2**s)]) for i in range(2**(s-1))]


def dmera_1d(n, D):
    """
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale

    Returns:
        list(list(tuple)): A list containing lists. Each (sub-)list
                           contains tuples, each of which represent
                           two-qubit gates.
    """
    circ = []
    for s in range(n):
        for j in range(D):
            if j%2==0:
                circ.append(nn_e(n,s+1))
            else:
                circ.append(nn_o(n,s+1))
    
    return circ


def pcc_dmera_1d(n, D, supp):
    """
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers

    Returns:
        list(list(tuple)): Past causal cone of supp
    """
    return pcc(dmera_1d(n,D), supp)


def width_pcc_dmera_1d_nocompression(n, D, supp):
    """
    Optimal width of the circuit for the pcc after compression
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers    

    Returns:
        int: Optimal width
    """
    c = pcc(dmera_1d(n,D), supp)
    width_max = 0
    for s in c:
        mywidth = np.sum([len(g) for g in s])
        if mywidth > width_max:
            width_max = mywidth
    return width_max
    

def width_pcc_dmera_1d(n, D, supp):
    """
    Optimal width of the circuit for the pcc after compression
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers    

    Returns:
        int: Optimal width
    """
    return optimal_width_freeze(pcc_dmera_1d(n,D,supp), supp)


def depth_pcc_dmera_1d(n, D, supp):
    """
    Depth of the width-reduced past causal cone
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers    

    Returns:
        int: Optimal width
    """
    return depth_width_reduced_freeze(pcc_dmera_1d(n,D,supp), supp)
