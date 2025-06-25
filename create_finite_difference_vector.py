import numpy as np

def diff_vec(V):
    npoints = len(V)
    diffV = np.zeros(npoints)

    # Difference at the first node
    diffV[0] = V[1] - V[0]

    # Difference for the interior nodes
    for ii in range(1, npoints - 1):
        diffV[ii] = (V[ii + 1] - V[ii - 1]) / 2

    # Difference at the last node
    diffV[-1] = V[-1] - V[-2]

    return diffV