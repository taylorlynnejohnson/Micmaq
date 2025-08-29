import numpy as np

def exergy_stored_in_phase(PhaseFileName='TCF_Fig22-223.npz', MinUsefulTemp=400, ChargingTemp=550):
  
    # Load phase data
    data = np.load(PhaseFileName)
    Ts = data['Ts']
    r = data['r']
    time = data['time']
    B = data['B']
    dr = r[1] - r[0]
    NPoints = len(r)
    T_ref = data['T_ref']

    exergy = 0
    volume = 0

    for inode in range(NPoints):
        dvol = 2 * np.pi * B * r[inode] * dr
        c_s, k_s, mu_s, rho_s = material_library('rock', Ts[inode])
        if Ts[inode] >= MinUsefulTemp:
            Texergetic = max(0, (Ts[inode] - T_ref))
            volume += dvol
            deltaExergy = c_s * Texergetic * dvol * rho_s
            exergy += deltaExergy

    exergy_MWhr = exergy / (3600E6)

    volCheck = (np.pi * r[-1]**2 - np.pi * r[0]**2) * B
    exergyCapacity = volume * rho_s * c_s * (ChargingTemp - MinUsefulTemp)
    exergyCapacity_MWhr = exergyCapacity / (3600E6)

    endTime = time

    return exergy_MWhr, exergyCapacity_MWhr, endTime
