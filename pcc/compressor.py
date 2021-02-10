################p###################
# Circuit compressor              #
# Isaac H. Kim 2020/12/25         #
###################################
import numpy as np
import copy


def ancestor_width(circ, supp, verbose=False):
    """
    Args:
        circ(list(list(tuple))): Circuit
        supp(list): List of integers

    Returns:
        int: Width of the past causal cone of supp
    """
    circ_rev= circ[::-1]
    supp_coded = 0
    for s in supp:
        supp_coded |= (1<<s)

    for unitcirc in circ_rev:
        for gate in unitcirc:
            if verbose:
                print("gate={}".format(gate))
            if (1<<gate[0]) & supp_coded:
                if not ((1<<gate[1]) & supp_coded):
                    supp_coded |= (1<<gate[1])
            elif (1<<gate[1]) & supp_coded:
                supp_coded |= (1<<gate[0])
    return bin(supp_coded).count('1')

    
def pcc(circ, supp, verbose=False):
    """
    Args:
        circ(list(list(tuple))): Circuit
        supp(list): List of integers

    Returns:
        list(list(tuple)): Past causal cone of supp
    """
    circ_rev= circ[::-1]
    circ_reduced = []
    supp_coded = 0
    for s in supp:
        supp_coded |= (1<<s)

    for unitcirc in circ_rev:
        circ_temp = []
        for gate in unitcirc:
            if verbose:
                print("gate={}".format(gate))
            if (1<<gate[0]) & supp_coded:
                circ_temp.append(gate)
                if not ((1<<gate[1]) & supp_coded):
                    supp_coded |= (1<<gate[1])
            elif (1<<gate[1]) & supp_coded:
                circ_temp.append(gate)
                supp_coded |= (1<<gate[0])
        circ_reduced.append(circ_temp)
    return circ_reduced[::-1]


def pcc_slower1(circ, supp, verbose=False):
    """
    Args:
        circ(list(list(tuple))): Circuit
        supp(list): List of integers

    Returns:
        list(list(tuple)): Past causal cone of supp
    """
    circ_rev= circ[::-1]
    circ_reduced = []
    supp_temp = copy.deepcopy(supp)
    for unitcirc in circ_rev:
        circ_temp = []
        for gate in unitcirc:
            if verbose:
                print("gate={}".format(gate))
            if gate[0] in supp_temp:
                circ_temp.append(gate)
                if not (gate[1] in supp_temp):
                    supp_temp.append(gate[1])
            elif gate[1] in supp_temp:
                circ_temp.append(gate)
                supp_temp.append(gate[0])
        circ_reduced.append(circ_temp)
    return circ_reduced[::-1]


def pcc_slower2(circ, supp, verbose=False):
    """
    Args:
        circ(list(list(tuple))): Circuit
        supp(list): List of integers

    Returns:
        list(list(tuple)): Past causal cone of supp
    """
    circ_rev= circ[::-1]
    circ_reduced = []
    supp_temp = copy.deepcopy(supp)
    for unitcirc in circ_rev:
        circ_temp = []
        for gate in unitcirc:
            if verbose:
                print("gate={}".format(gate))
            if gate[0] in supp_temp or gate[1] in supp_temp:
                circ_temp.append(gate)
                supp_temp += list(gate)
        circ_reduced.append(circ_temp)
        supp_temp = list(set(supp_temp))
    return circ_reduced[::-1]


def qubits(circ):
    """
    Args:
        circ(list(list(tuple))): Circuit

    Returns:
        list: The list of qubits appearing in the circuit.
    """
    qs = []
    for c in circ:
        for g in c:
            qs += list(g)
        qs = list(set(qs))
    return qs


def circ_subtract(circ1, circ2):
    """
    Let circ1 and circ2 be lists of depth-1 circuits.
    This function subtracts the gates in the k'th depth-1 circuit of circ2
    from that of circ1. Note that this function may lead to misbehaviors
    if circ2 is not a sub-circuit of circ1.

    Args:
        circ1, circ2(list(list(tuple))): Circuits

    Returns:
        list(list(tuple)): The remaining circuit
    """
    if len(circ1) != len(circ2):
        raise ValueError("The inputs must have identical lengths.")

    circ_out = []
    for c1, c2 in zip(circ1, circ2):
        circ_out.append([x for x in c1 if x not in c2])
    return circ_out


def rearrange_freeze_old(circ, frozen):
    """
    Rearrange the circuit for width reduction
    Args:
        circ(list(list(tuple))): Circuit
        frozen(list): Qubits that shouldn't be reset.
    Returns:
        list: Rearranged circuit
    """
    circ_temp = circ
    circ_out = []
    while len(qubits(circ_temp)) > 0:
        qs = qubits(circ_temp)
        n_q_min = len(qs)
        q_min_idx = qs[0]
        for q in qs:
            pcc_temp = pcc(circ_temp, [q])
            if len(qubits(pcc_temp)) < n_q_min:
                n_q_min = len(qubits(pcc_temp))
                q_min_idx = q
        mypcc = pcc(circ_temp, [q_min_idx])
        circ_temp = circ_subtract(circ_temp, mypcc)
        circ_out += [c for c in mypcc if len(c)>0]
#        if q_min_idx in frozen:
#            print("{} is frozen".format(q_min_idx))
#            print("Frozen list={}".format(frozen))
        if not (q_min_idx in frozen):
            circ_out += [q_min_idx]
    return circ_out


def rearrange_freeze(circ, frozen):
    """
    Rearrange the circuit for width reduction
    Args:
        circ(list(list(tuple))): Circuit
        frozen(list): Qubits that shouldn't be reset.
    Returns:
        list: Rearranged circuit
    """
    circ_temp = circ
    circ_out = []
    while len(qubits(circ_temp)) > 0:
        qs = qubits(circ_temp)
        n_q_min = len(qs)
        q_min_idx = qs[0]
        for q in qs:
            width_temp = ancestor_width(circ_temp, [q])
            if width_temp < n_q_min:
                n_q_min = width_temp
                q_min_idx = q
        mypcc = pcc(circ_temp, [q_min_idx])
        circ_temp = circ_subtract(circ_temp, mypcc)
        circ_out += [c for c in mypcc if len(c)>0]
#        if q_min_idx in frozen:
#            print("{} is frozen".format(q_min_idx))
#            print("Frozen list={}".format(frozen))
        if not (q_min_idx in frozen):
            circ_out += [q_min_idx]
    return circ_out


def rearrange(circ):
    """
    Rearrange the circuit for width reduction
    Args:
        circ(list(list(tuple))): Circuit

    Returns:
        list: Rearranged circuit
    """
    circ_temp = circ
    circ_out = []
    while len(qubits(circ_temp)) > 0:
        qs = qubits(circ_temp)
        n_q_min = len(qs)
        q_min_idx = qs[0]
        for q in qs:
            pcc_temp = pcc(circ_temp, [q])
            if len(qubits(pcc_temp)) < n_q_min:
                n_q_min = len(qubits(pcc_temp))
                q_min_idx = q
        mypcc = pcc(circ_temp, [q_min_idx])
        circ_temp = circ_subtract(circ_temp, mypcc)
        circ_out += [c for c in mypcc if len(c)>0]
        circ_out += [q_min_idx]
    return circ_out


def compress(circ):
    """
    Compress the circuit for width reduction
    Args:
        circ(list(list(tuple))): Circuit

    Returns:
        list: Compressed circuit
    """
    circ_compressed = []
    q_max = max(qubits(circ))+1
    circ_r = rearrange(circ)
    q_available = [i for i in range(q_max)][::-1]
    q_using = []
    assignment = {i:None for i in range(q_max)}
    for c in circ_r:
        c_compressed = []
        if type(c)==list:
            for g in c:
                if assignment[g[0]] == None:
                    q_new = q_available.pop()
                    assignment[g[0]]=q_new
                    q_using.append(q_new)
                if assignment[g[1]] == None:
                    q_new = q_available.pop()
                    assignment[g[1]]=q_new
                    q_using.append(q_new)
                g_new = (assignment[g[0]], assignment[g[1]])
                c_compressed.append(g_new)
        elif type(c)==int:
            # Remove the qubits in q_using assigned to c.
            q_using.remove(assignment[c])
            q_available.append(assignment[c])
            # Following is added to show the reset
            c_compressed.append([assignment[c]])
        if len(c_compressed)>0:
            circ_compressed.append(c_compressed)
    return circ_compressed


def compress_freeze(circ, frozen):
    """
    Compress the circuit for width reduction
    Args:
        circ(list(list(tuple))): Circuit
        frozen(list): Frozen qubits. (They shouldn't be reset.)
    Returns:
        list: Compressed circuit
    """
    circ_compressed = []
    q_max = 0
    if qubits(circ):
        q_max = max(qubits(circ))+1
    circ_r = rearrange_freeze(circ, frozen)
    q_available = [i for i in range(q_max)][::-1]
    q_using = []
    assignment = {i:None for i in range(q_max)}
    for c in circ_r:
        c_compressed = []
        if type(c)==list:
            for g in c:
                if assignment[g[0]] == None:
                    q_new = q_available.pop()
                    assignment[g[0]]=q_new
                    q_using.append(q_new)
                if assignment[g[1]] == None:
                    q_new = q_available.pop()
                    assignment[g[1]]=q_new
                    q_using.append(q_new)
                g_new = (assignment[g[0]], assignment[g[1]])
                c_compressed.append(g_new)
        elif type(c)==int:
            # Remove the qubits in q_using assigned to c.
            q_using.remove(assignment[c])
            q_available.append(assignment[c])
            # Following is added to show the reset
            c_compressed.append([assignment[c]])
        if len(c_compressed)>0:
            circ_compressed.append(c_compressed)
    return circ_compressed
        

def optimal_width(circ):
    """
    Optimal width of the circuit after compression
    Args:
        circ(list(list(tuple))): Circuit

    Returns:
        int: Optimal width
    """
    return len(qubits(compress(circ)))


def depth_width_reduced(circ):
    """
    Depth of the width-reduced circuit
    Args:
        circ(list(list(tuple))): Circuit

    Returns:
        int: Optimal width
    """
    return len(compress(circ))


def optimal_width_freeze(circ, frozen):
    """
    Optimal width of the circuit after compression
    Args:
        circ(list(list(tuple))): Circuit
        frozen(list): Frozen qubits
    Returns:
        int: Optimal width
    """
    return len(qubits(compress_freeze(circ, frozen)))


def depth_width_reduced_freeze(circ, frozen):
    """
    Depth of the width-reduced circuit
    Args:
        circ(list(list(tuple))): Circuit
        frozen(list): Frozen aubits
    Returns:
        int: Optimal width
    """
    return len(compress_freeze(circ, frozen))
