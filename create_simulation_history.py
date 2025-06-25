import numpy as np

def create_simulation_history(chargeTemp, ambientTemp, chargeHours, holdHours, dischargeHours, chargeFlowRate, dischargeFlowRate):
     
    chargePhase = np.tile([chargeFlowRate, 1, chargeTemp, ambientTemp], (chargeHours, 1))

    holdPhase = np.tile([0.00, 1, ambientTemp, ambientTemp], (holdHours, 1))
    
    dischargePhase = np.tile([dischargeFlowRate, 1, ambientTemp, ambientTemp], (dischargeHours, 1))

    # concatenate all phases to form the simHist matrix
    simHist = np.vstack((chargePhase, holdPhase, dischargePhase))

    return simHist

