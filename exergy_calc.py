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

    # Read experimental data
    ExpData = pd.read_csv('12-15-2022.txt', delimiter='\t')
    times_Exp = ExpData.iloc[1:, 0].values * 60  # in seconds
    deltaTimeExp = times_Exp[1] - times_Exp[0]  # in seconds
    A1_Exp_T = ExpData.iloc[1:, 6].values
    InjectWell_T = ExpData.iloc[1:, 27].values

    # Experimental energy calculations
    t_beginPhase = 0
    Energy_In_A1_Exp = 0
    Energy_Out_A1_Exp = 0

    for iPhase in range(numPhases):
        M_dot = simHist[iPhase, 0]
        phaseDuration = simHist[iPhase, 1] * 3600
        T_s0 = ExpData.iloc[0, 8]
        t_endPhase = t_beginPhase + phaseDuration

        indStart = np.searchsorted(times_Exp, t_beginPhase)
        indEnd = np.searchsorted(times_Exp, t_endPhase)

        if M_dot > 0:
            for ind in range(indStart, indEnd):
                T_g = A1_Exp_T[ind]
                c_g, k_g, mu_g, rho_g = material_library('air', T_g)
                TempDiff = T_g - T_ref
                if T_g > MinUsefulTemp:
                    Energy_In_A1_Exp += M_dot * c_g * TempDiff * deltaTimeExp
        elif M_dot < 0:
            for ind in range(indStart, indEnd):
                T_g = A1_Exp_T[ind]
                c_g, k_g, mu_g, rho_g = material_library('air', T_g)
                TempDiff = T_g - T_ref
                if T_g > MinUsefulTemp:
                    Energy_Out_A1_Exp += M_dot * c_g * TempDiff * deltaTimeExp

        t_beginPhase = t_endPhase

    Energy_In_A1_Exp /= (3600 * 1000)  # kWh
    Energy_Out_A1_Exp /= (3600 * 1000)  # kWh

    return Energy_In_A1_Sim, Energy_Out_A1_Sim, Energy_In_A1_Exp, Energy_Out_A1_Exp

def material_library(material, temperature):
    # Example material library function for 'air'
    if material == 'air':
        c = 1005  # Specific heat capacity [J/kg-K]
        k = 0.03  # Thermal conductivity [W/m-K]
        mu = 1.8e-5  # Dynamic viscosity [Pa-s]
        rho = 1.2  # Density [kg/m^3]
    return c, k, mu, rho

