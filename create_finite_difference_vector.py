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


        geoms = {
    "1:1":  {"rR": 0.8171466, "B": 0.8171466},
    "1:2":  {"rR": 1.0193330, "B": 0.5096665},
    "1:3":  {"rR": 1.1621323, "B": 0.3873774},
    "1:4":  {"rR": 1.2761753, "B": 0.3190438},
    "1:5":  {"rR": 1.3726637, "B": 0.2745327},
    "1:6":  {"rR": 1.4571150, "B": 0.2428525},
    "1:7":  {"rR": 1.5327045, "B": 0.2189578},
    "1:8":  {"rR": 1.6014473, "B": 0.2001809},
    "1:9":  {"rR": 1.6647121, "B": 0.1849680},
    "1:10": {"rR": 1.7234761, "B": 0.1723476},
}
        
