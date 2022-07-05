import numpy as np

def compAdd(contribution, compSum):
    """Add a single contribution to a compensated summation.
    
    Add a single contribution to the sum.
    The compensated sum is obtained by summing the final values of s, cs and ccs
    after this method has been called for all contributions.

    :param float contribution: contribution to add to the sum
    :param arr(float) compSum: [3: s, cs, ccs] target for output
    """
    s  = compSum[0]
    cs = compSum[1]

    t = s + contribution
    if np.abs(s) >= np.abs(contribution):
        c = (s - t) + contribution
    else:
        c = (contribution - t) + s
    compSum[0] = t

    t2 = cs + c
    if np.abs(cs) >= np.abs(c):
        cc = (cs - t2) + c
    else:
        cc = (c - t2) + cs
    compSum[1] = t2
    compSum[2] += cc
