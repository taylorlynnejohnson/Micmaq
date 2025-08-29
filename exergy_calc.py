import numpy as np
import pandas as pd
import material_library_fn

def exergy(fileName, MinUsefulTemp=20):
 
    # Load simulation data
    data = np.load(fileName)
    timeHistory = data['timeHistory']
    simHist = data['simHist']
    numPhases = simHist.shape[0]

    timesSim = timeHistory[:, 0]
    A1_Sim_T_g = timeHistory[:, 1]

    # Simulation energy calculations
    Energy_In_WellSim = 0
    Energy_In_A1_Sim = 0
    Energy_Out_A1_Sim = 0
    t_beginPhase = 0

    for iPhase in range(numPhases):
        M_dot = simHist[iPhase, 0]
        phaseDuration = simHist[iPhase, 1] * 3600
        T_g0 = simHist[iPhase, 2]
        T_s0 = simHist[iPhase, 3]
        t_endPhase = t_beginPhase + phaseDuration

        if M_dot > 0:
            c_g, k_g, mu_g, rho_g = material_library_fn('air', T_g0)
            Energy_In_WellSim += M_dot * c_g * (T_g0 - T_s0) * phaseDuration

            indStart = np.searchsorted(timesSim, t_beginPhase)
            indEnd = np.searchsorted(timesSim, t_endPhase)
            for ind in range(indStart, indEnd):
                T_g = A1_Sim_T_g[ind]
                c_g, k_g, mu_g, rho_g = material_library_fn('air', T_g)
                TempDiff = T_g - T_ref
                if T_g > MinUsefulTemp:
                    Energy_In_A1_Sim += M_dot * c_g * TempDiff * TraceTimeIncrement
        elif M_dot < 0:
            c_g, k_g, mu_g, rho_g = material_library_fn('air', T_g0)
            indStart = np.searchsorted(timesSim, t_beginPhase)
            indEnd = np.searchsorted(timesSim, t_endPhase)
            for ind in range(indStart, indEnd):
                T_g = A1_Sim_T_g[ind]
                c_g, k_g, mu_g, rho_g = material_library_fn('air', T_g)
                TempDiff = T_g - T_ref
                if T_g > MinUsefulTemp:
                    Energy_Out_A1_Sim += M_dot * c_g * TempDiff * TraceTimeIncrement

        t_beginPhase = t_endPhase

    Energy_In_WellSim /= (3600 * 1000)  # kWh
    Energy_In_A1_Sim /= (3600 * 1000)  # kWh
    Energy_Out_A1_Sim /= (3600 * 1000)  # kWh


    return Energy_In_A1_Sim, Energy_Out_A1_Sim

