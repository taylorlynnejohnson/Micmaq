# one_D_radial_thermocline.py

import numpy as np
from material_library_fn import material_library_fn
from one_D_model import calculate_heat_loss, outer_area

def radial_thermocline(
    Ts, Tg, HTF_type, r,
    numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
    h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
    M_dot, time, dt,
    Tg0, T_ref,
    NodesToTraceInTime, TraceTimeIncrement,
    timeHistory
):
    numHistNodes = len(NodesToTraceInTime)
    timeHistRow  = len(timeHistory)
    # Pre‐extend history to avoid repeated vstack
    numHistOutputs = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
    timeHistory   = np.vstack([timeHistory,
                               np.zeros((numHistOutputs, 1 + 2 * numHistNodes))])

    energy_in_left  = 0.0
    energy_out_right = 0.0
    pumping_energy   = 0.0

    NPoints = len(Ts)
    dr      = r[1] - r[0]

    # Precompute segment data
    Area   = 2.0 * np.pi * r * B                   # m² at each node
    vol_s  = dr * B * 2.0 * np.pi * r               # m³ in each radial segment
    c_s, k_s, mu_s, rho_s = material_library_fn('rock', 20)

    lastHistOutput = time

    dp_time_integral = 0.0
    sim_time_total   = 0.0

    for iStep in range(numTimeSteps):
        time += dt
        if iStep % 100000 == 0:
            print(f"Step {iStep}/{numTimeSteps}: t = {time:.1f}s, M_dot = {M_dot:.4f}")

        # 1) Compute gradients
        Tg_r = np.gradient(Tg, dr)
        Ts_r = np.gradient(Ts, dr)
        rTs_r = r * Ts_r
        dr_rTs_r = np.gradient(rTs_r, dr)
        dr_rTs_r_over_r = dr_rTs_r / r

        # 2) HTF & convection params
        c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
        v_s     = M_dot / (rho_g * Area)
        v_local = v_s / eps
        Re = rho_g * np.abs(v_s) * Dp / mu_g
        Pr = c_g * mu_g / k_g

        # Gunn Nusselt
        Nu = ((7 - 10*eps + 5*eps**2)*(1 + 0.7*Re**0.2 * Pr**0.33)
              + (1.33 - 2.4*eps + 1.2*eps**2)*Re**0.7 * Pr**0.33)
        h_sg = (k_g / Dp) * Nu
        H_sg = 6.0 * h_sg * (1 - eps) / (alpha * Dp)

        # Pressure drop
        fv    = 610.0/Re + 13.0/(Re**0.11)
        dp_dr = rho_g*(v_s**2)*(1-eps)/(2*Dp*eps**3)*fv


        # accumulate dp
        dp_inst = np.sum(dp_dr) * dr
        dp_time_integral += dp_inst * dt
        sim_time_total   += dt


        # 3) Gas temperature update (advection + conv)
        Tg += (H_sg/(rho_g*c_g*eps)*(Ts - Tg) - v_local*Tg_r)*dt
        if M_dot > 0:
            Tg[0] = Tg0
        elif M_dot < 0:
            Tg[-1] = Tg0

        # 4) Solid temperature update (conv + diffusion)
        k1 = H_sg/(rho_s*c_s*(1-eps))
        Ts += k1*(Tg - Ts)*dt
        Ts += alpha_gravel * dr_rTs_r_over_r * dt

        # 5) External heat loss via multi‐layer model at true boundaries
        # Outer cylinder
        Q_outer    = calculate_heat_loss(Ts[-1], T_ref)
        mcp_outer  = rho_s * vol_s[-1] * c_s
        Ts[-1]    -= (Q_outer * dt) / mcp_outer

        # Inner hole surface
        Q_inner   = calculate_heat_loss(Ts[0], T_ref)
        mcp_inner = rho_s * vol_s[0] * c_s
        Ts[0]    -= (Q_inner * dt) / mcp_inner

        # Top & bottom faces (flat annulus)
        A_flat    = np.pi*(r[-1]**2 - r[0]**2)
        Q_full    = calculate_heat_loss(np.mean(Ts), T_ref)
        Q_tb      = Q_full * (A_flat / outer_area)
        mcp_total = np.sum(rho_s * vol_s * c_s)
        dT_tb     = (Q_tb * dt) / mcp_total
        Ts       -= dT_tb

        # 6) Track energy flows
        energy_in_left  += M_dot * (Tg[0] - T_ref) * c_g[0] * dt
        energy_out_right+= M_dot * (Tg[-1] - T_ref)* c_g[-1] * dt
        pumping_energy  += np.sum(Area*np.abs(v_s)*dp_dr)*dr * dt

        # 7) History snapshots
        if (time - lastHistOutput) >= TraceTimeIncrement:
            lastHistOutput = time
            timeHistRow   += 1
            if timeHistRow >= timeHistory.shape[0]:
                timeHistory = np.vstack([timeHistory,
                                         np.zeros((1, timeHistory.shape[1]))])
            timeHistory[timeHistRow, 0] = time
            for jNode in range(numHistNodes):
                Tgas_j = Tg[NodesToTraceInTime[jNode]]
                Tsol_j = Ts[NodesToTraceInTime[jNode]]
                timeHistory[timeHistRow, 2*jNode+1:2*jNode+3] = [Tgas_j, Tsol_j]

    # 8) Final stored energy and return
    energy_in    = energy_in_left - energy_out_right
    stored_energy= (
        np.sum((Ts - T_ref)*c_s*rho_s*vol_s*(1-eps))
      + np.sum((Tg - T_ref)*c_g*rho_g*vol_s*eps)
    )
    dp_avg = dp_time_integral / sim_time_total

    return Ts, Tg, timeHistory, time, energy_in, stored_energy, pumping_energy, dp_avg


# # one_D_radial_thermocline.py
# import numpy as np
# from material_library_fn import material_library_fn
# from one_D_model import calculate_heat_loss

# def radial_thermocline(Ts, Tg, HTF_type, r, numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
#                       h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner, M_dot, time, dt,
#                       Tg0, T_ref, NodesToTraceInTime, TraceTimeIncrement, timeHistory):
    
#     numHistNodes = len(NodesToTraceInTime)
#     timeHistRow = len(timeHistory)
#     numHistOutputs = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
#     timeHistory = np.vstack([timeHistory, np.zeros((numHistOutputs, 1 + 2 * numHistNodes))])
#     energy_in_left = 0
#     energy_out_right = 0
#     pumping_energy = 0
#     NPoints = len(Ts)
#     dr = r[1] - r[0]   
#     Area = 2 * np.pi * r * B  # cross-sectional area at each radial position
#     vol_s = dr * B * 2 * np.pi * r  # v of solid material in each radial segment
#     # rock material props
#     c_s, k_s, mu_s, rho_s = material_library_fn('rock', 20)
#     lastHistOutput = time

#     for iStep in range(numTimeSteps):
#         time += dt
#         if iStep % 100000 == 0:
#             print(f"Step {iStep}/{numTimeSteps}: Time = {time:.2f}s, M_dot = {M_dot:.4f} kg/s")
            
#         # temperature gradients
#         Tg_r = np.gradient(Tg, dr)  # rad gradient of gas temperature
#         Ts_r = np.gradient(Ts, dr)  # ^ of solid temperature
#         r_Ts_r = r * Ts_r  
#         dr_r_Ts_r = np.gradient(r_Ts_r, dr)  
#         dr_r_Ts_r_over_r = dr_r_Ts_r / r  

#         # material properties for the gas HTF
#         c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
#         v_s = M_dot / (rho_g * Area)  # air velocity
#         v_local = v_s / eps  # local air velocity in the void fraction
#         Re = rho_g * np.abs(v_s) * Dp / mu_g  # Reynolds number for packed bed
#         Pr = c_g * mu_g / k_g  # Prandtl number

#         # Nusselt number calc
#         Nu_Gunn = ((7 - (10 * eps) + 5 * (eps ** 2)) * (1 + 0.7 * (Re ** 0.2) * (Pr ** 0.33)) +
#                    (1.33 - (2.4 * eps) + 1.2 * (eps ** 2)) * (Re ** 0.7) * (Pr ** 0.33))
#         h_sg = (k_g / Dp) * Nu_Gunn  # heat transfer coefficient between gas and solid
#         H_sg = 6 * h_sg * (1 - eps) / (alpha * Dp)  # total heat transfer coefficient
#         Bi = h_sg * Dp / (2 * k_s)  # Biot number 
#         f_v = 610.0 / Re + 13.0 / (Re ** 0.11)  # friction factor for packed bed
#         dp_dr = rho_g * (v_s ** 2) * (1 - eps) / (2 * Dp * (eps ** 3)) * f_v  # pressure drop per unit length

#         # update gas temperature using convective heat transfer and advection
#         k2 = H_sg / (c_g * rho_g * eps)
#         Tg += (k2 * (Ts - Tg) - v_local * Tg_r) * dt

#         # apply BCs for gas temperature
#         if M_dot > 0:
#             Tg[0] = Tg0  # gas BC at left end of domain (inlet temperature)
#         if M_dot < 0:
#             Tg[-1] = Tg0  # gas BC at right end of domain (inlet temperature)

#         # update solid temperature using convective heat transfer from HTF
#         k1 = H_sg / (c_s * rho_s * (1 - eps))
#         Ts += k1 * (Tg - Ts) * dt  # convective heat transfer from gas to solid
#         Ts += alpha_gravel * dr_r_Ts_r_over_r * dt  # radial thermal diffusion in packed bed

#         # apply 1D heat loss model
#         avg_temp_surface = np.mean(Ts)  # avg solid temperature
#         heat_loss_rate = calculate_heat_loss(avg_temp_surface, T_ref)
#         Ts -= (heat_loss_rate * dt) / (c_s * rho_s * vol_s)  # update solid temperature with heat loss

#         # energy input and output
#         power_in_left = M_dot * (Tg[0] - T_ref) * c_g[0]
#         energy_in_left += power_in_left * dt
#         power_out_right = M_dot * (Tg[-1] - T_ref) * c_g[-1]
#         energy_out_right += power_out_right * dt
#         pumping_power = np.sum(Area * np.abs(v_s) * dp_dr) * dr
#         pumping_energy += pumping_power * dt

#         if (time - lastHistOutput) >= TraceTimeIncrement:
#             lastHistOutput = time
#             timeHistRow += 1
#             if timeHistRow >= timeHistory.shape[0]:
#                 timeHistory = np.vstack([timeHistory, np.zeros((1, timeHistory.shape[1]))])
#             timeHistory[timeHistRow, 0] = time
#             for iNode in range(numHistNodes):
#                 timeHistory[timeHistRow, 2 * iNode + 1: 2 * iNode + 3] = [Tg[NodesToTraceInTime[iNode]], Ts[NodesToTraceInTime[iNode]]]

#     energy_in = energy_in_left - energy_out_right
#     stored_energy = (np.sum((Ts - T_ref) * c_s * vol_s * rho_s * (1 - eps)) +
#                      np.sum((Tg - T_ref) * c_g * vol_s * rho_g * eps))

#     return Ts, Tg, timeHistory, time, energy_in, stored_energy, pumping_energy