###################################
# Past causal cone calculator: 1D #
# Isaac H. Kim 2020/12/24         #
###################################
import numpy as np
from compressor import pcc, qubits, circ_subtract, rearrange, compress, optimal_width_freeze, depth_width_reduced_freeze


def sites_2d(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(int): A list of sites that appear up to the s'th scale.
    """
    myslice = [i*(2**(n-s)) for i in range(2**s)]
    return {(x,y): (myslice[x], myslice[y]) for x in range(2**s) for y in range(2**s)}


def nn_xe(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(tuple): A list of nearest neighbors at scale s. (even)
    """
    sites = sites_2d(n,s)
    return [[sites[2*x,y], sites[(2*x+1)%(2**s),y]] for x in range(2**(s)//2) for y in range(2**s)]


def nn_xo(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(tuple): A list of nearest neighbors at scale s. (even)
    """
    sites = sites_2d(n,s)
    return [[sites[(2*x+1)%(2**s),y], sites[(2*x+2)%(2**s),y]] for x in range(2**(s)//2) for y in range(2**s)]


def nn_ye(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(tuple): A list of nearest neighbors at scale s. (even)
    """
    sites = sites_2d(n,s)
    return [[sites[x,2*y], sites[x,(2*y+1)%(2**s)]] for x in range(2**s) for y in range(2**(s)//2)]


def nn_yo(n, s):
    """
    Args:
        n(int): Number of scales
        s(int): The index of the scale. The top is s=1, the second is s=2, etc.

    Returns:
        list(tuple): A list of nearest neighbors at scale s. (even)
    """
    sites = sites_2d(n,s)
    return [[sites[x,(2*y+1)%(2**s)], sites[x,(2*y+2)%(2**s)]] for x in range(2**s) for y in range(2**(s)//2)]


def dmera_2d_preconversion(n, D):
    """
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale

    Returns:
        list(list(tuple(tuple))): A list containing lists. Each (sub-)list
                           contains tuples of tuples, each of which represent
                           two-qubit gates. The output format is *not* consistent
                           with the input format for pcc. This is for testing 
                           purposes only.
    """
    circ = []
    for s in range(n):
        for j in range(D):
            if j%4==0:
                circ.append(nn_xe(n,s+1))
            elif j%4==1:
                circ.append(nn_ye(n,s+1))
            elif j%4==2:
                circ.append(nn_xo(n,s+1))
            else:
                circ.append(nn_yo(n,s+1))

    return circ


def convert_2d_to_1d(coordinate, n):
    """
    Convert a 2D coordinate to a coordinate on a line.

    Args:
        coordinate(tuple): 2D coordinate
        n(int): Number of scales

    Returns:
        int: 1D coordinate
    """
    return coordinate[0] + 2**n * coordinate[1]


def dmera_2d(n, D):
    """
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale

    Returns:
        list(list(tuple)): A list containing lists. Each (sub-)list
                           contains tuples of tuples, each of which represent
                           two-qubit gates. The qubits are represented by
                           integers instead of tuples. This is the one that
                           can be plugged into the pcc function.
    """
    dmera_raw = dmera_2d_preconversion(n, D)
    dmera_out = []
    for c in dmera_raw:
        c_new = []
        for g in c:
            q0 = convert_2d_to_1d(g[0], n)
            q1 = convert_2d_to_1d(g[1], n)
            c_new.append((q0, q1))
        dmera_out.append(c_new)
    return dmera_out


def pcc_dmera_2d(n, D, supp):
    """
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers

    Returns:
        list(list(tuple)): Past causal cone of supp
    """
    return pcc(dmera_2d(n,D),supp)


def width_pcc_dmera_2d(n, D, supp):
    """
    Optimal width of the circuit for the pcc after compression
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers    

    Returns:
        int: Optimal width
    """
    supp_con = [convert_2d_to_1d(c,n) for c in supp]
    return optimal_width_freeze(pcc_dmera_2d(n,D,supp_con),supp_con)


def width_pcc_dmera_2d_nocompression(n, D, supp):
    """
    Uncompressed width of the circuit for the pcc after compression
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers    

    Returns:
        int: Uncompresed width
    """
    supp_con = [convert_2d_to_1d(c,n) for c in supp]
    c = pcc_dmera_2d(n,D,supp_con)
    width_max = 0
    for s in c:
        mywidth = np.sum([len(g) for g in s])
        if mywidth > width_max:
            width_max = mywidth
    return width_max


def depth_pcc_dmera_2d(n, D, supp):
    """
    Depth of the width-reduced past causal cone
    Args:
        n(int): Number of scales
        D(int): Number of cycles per scale
        supp(list): List of integers    

    Returns:
        int: Optimal width
    """
    supp_con = [convert_2d_to_1d(c,n) for c in supp]
    return depth_width_reduced_freeze(pcc_dmera_2d(n,D,supp_con),supp_con)
