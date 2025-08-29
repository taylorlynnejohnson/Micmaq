# radial_thermocline.py

import numpy as np
from material_library_fn import material_library_fn

def radial_thermocline(Ts, Tg, HTF_type, r, numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
                      h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner, M_dot, time, dt,
                      Tg0, T_ref, NodesToTraceInTime, TraceTimeIncrement, timeHistory):
    
    numHistNodes = len(NodesToTraceInTime)
    timeHistRow = len(timeHistory)
    numHistOutputs = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
    timeHistory = np.vstack([timeHistory, np.zeros((numHistOutputs, 1 + 2 * numHistNodes))])
    energy_in_left = 0
    energy_out_right = 0
    pumping_energy = 0
    NPoints = len(Ts)
    dr = r[1] - r[0]   
    Area = 2 * np.pi * r * B  # Cross-sectional area at each radial position
    vol_s = dr * B * 2 * np.pi * r  # Volume of solid material in each radial segment
    # rock material props
    c_s, k_s, mu_s, rho_s = material_library_fn('rock', 20)
    lastHistOutput = time

    for iStep in range(numTimeSteps):
        time += dt
        if iStep % 100000 == 0:
            print(f"Step {iStep}/{numTimeSteps}: Time = {time:.2f}s, M_dot = {M_dot:.4f} kg/s")
        
        # temperature gradients
        Tg_r = np.gradient(Tg, dr)  # Radial gradient of gas temperature
        Ts_r = np.gradient(Ts, dr)  # Radial gradient of solid temperature
        r_Ts_r = r * Ts_r  
        dr_r_Ts_r = np.gradient(r_Ts_r, dr)  
        dr_r_Ts_r_over_r = dr_r_Ts_r / r  
        # material properties for the gas HTF
        c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
        v_s = M_dot / (rho_g * Area)  # air velocity
        v_local = v_s / eps  # local air velocity in the void fraction
        Re = rho_g * np.abs(v_s) * Dp / mu_g  # Reynolds number for packed bed
        Pr = c_g * mu_g / k_g  # Prandtl number

        # Nusselt number calc
        Nu_Gunn = ((7 - (10 * eps) + 5 * (eps ** 2)) * (1 + 0.7 * (Re ** 0.2) * (Pr ** 0.33)) +
                   (1.33 - (2.4 * eps) + 1.2 * (eps ** 2)) * (Re ** 0.7) * (Pr ** 0.33))
        h_sg = (k_g / Dp) * Nu_Gunn  # heat transfer coefficient between gas and solid
        H_sg = 6 * h_sg * (1 - eps) / (alpha * Dp)  # total heat transfer coefficient
        Bi = h_sg * Dp / (2 * k_s)  # Biot number 
        f_v = 610.0 / Re + 13.0 / (Re ** 0.11)  # friction factor for packed bed
        dp_dr = rho_g * (v_s ** 2) * (1 - eps) / (2 * Dp * (eps ** 3)) * f_v  # pressure drop per unit length

        # update gas temperature using convective heat transfer and advection
        k2 = H_sg / (c_g * rho_g * eps)
        Tg += (k2 * (Ts - Tg) - v_local * Tg_r) * dt

        # apply BCs for gas temperature
        if M_dot > 0:
            Tg[0] = Tg0  # gas BC at left end of domain (inlet temperature)
        if M_dot < 0:
            Tg[-1] = Tg0  # gas BC at right end of domain (inlet temperature)

        # update solid temperature using convective heat transfer from HTF
        k1 = H_sg / (c_s * rho_s * (1 - eps))
        Ts += k1 * (Tg - Ts) * dt  # convective heat transfer from gas to solid
        Ts += alpha_gravel * dr_r_Ts_r_over_r * dt  # radial thermal diffusion in packed bed

        # convective heat transfer to environment (top and bottom surfaces)
        # 1D model 
        k_3 = (h_surf_upper + h_surf_lower) / (rho_s * c_s * (1 - eps) * B)
        Ts -= k_3 * (Ts - T_ref) * dt

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
        

        # heat loss from outer cylindrical boundary
        beta_outer = h_surf_outer / (rho_s * c_s * (1 - eps) * dr)
        Ts[-1] -= beta_outer * (Ts[-1] - T_ref) * dt

        # heat loss from inner cylindrical boundary
        beta_inner = h_surf_inner * r[0] / (2 * B * rho_s * c_s * (1 - eps) * dr)
        Ts[0] -= beta_inner * (Ts[0] - T_ref) * dt

        # energy input and output
        power_in_left = M_dot * (Tg[0] - T_ref) * c_g[0]
        energy_in_left += power_in_left * dt
        power_out_right = M_dot * (Tg[-1] - T_ref) * c_g[-1]
        energy_out_right += power_out_right * dt
        pumping_power = np.sum(Area * np.abs(v_s) * dp_dr) * dr
        pumping_energy += pumping_power * dt

        if (time - lastHistOutput) >= TraceTimeIncrement:
            lastHistOutput = time
            timeHistRow += 1
            if timeHistRow >= timeHistory.shape[0]:
                # Dynamically add a new row if timeHistRow exceeds current size
                timeHistory = np.vstack([timeHistory, np.zeros((1, timeHistory.shape[1]))])
            timeHistory[timeHistRow, 0] = time
            for iNode in range(numHistNodes):
                timeHistory[timeHistRow, 2 * iNode + 1: 2 * iNode + 3] = [Tg[NodesToTraceInTime[iNode]], Ts[NodesToTraceInTime[iNode]]]

    energy_in = energy_in_left - energy_out_right
    stored_energy = (np.sum((Ts - T_ref) * c_s * vol_s * rho_s * (1 - eps)) +
                     np.sum((Tg - T_ref) * c_g * vol_s * rho_g * eps))

    return Ts, Tg, timeHistory, time, energy_in, stored_energy, pumping_energy

