
import numpy as np
from material_library_fn import material_library_fn
from one_D_model import calculate_axial_annular_roof_loss, calculate_bottom_annular_loss, outer_area

# constants
stefan_boltzmann = 5.67e-8  # W/(m²·K⁴)

# Material thermal conductivities (W/(m·K))
k_super_wool   = 0.19   # Super‑wool insulation
k_sand         = 0.114  # Sand
k_fiberglass   = 0.04   # Fiberglass blanket
k_g_steel      = 50     # Galvanized steel tank

# Layer thicknesses (m)
# super wool thickness optimized in main  
thickness_sand        = 0.0254    # Sand
thickness_fiberglass  = 0.0047625 # Fiberglass
thickness_tank        = 0.0047625 # Steel tank

# ground pad under tank
k_concrete      = 1.2    # W/m·K
thickness_ground= 0.5    # m

# surface properties
emissivity      = 0.05   # –
h_convection    = 5.0    # W/(m²·K)

# Geometry
inner_diameter  = 2*0.2286   # central hole r (m)

def radial_thermocline(
    Ts, Tg, HTF_type, r,
    numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
    h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
    M_dot, time, dt,
    Tg0, T_ref,
    NodesToTraceInTime, TraceTimeIncrement,
    timeHistory,
    MinUsefulTemp,            # mi temperature (C) for exergy 
    optimized_insul_thickness,
    rR
):
    
    numHistNodes    = len(NodesToTraceInTime)
    timeHistRow     = len(timeHistory)
    numHistOutputs  = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
    timeHistory     = np.vstack([
        timeHistory,
        np.zeros((numHistOutputs, 1 + 2 * numHistNodes))
    ])
    MinUsefulTemp += 273.15 # C to K
    charging_exergy = 0
    discharging_exergy = 0
    energy_in_left    = 0.0
    energy_out_right  = 0.0
    pumping_energy    = 0.0
    energy_charge_hist  = np.zeros(numTimeSteps)
    energy_discharge_hist = np.zeros(numTimeSteps)
    stored_energy_hist  = np.zeros(numTimeSteps)
    charge_exergy_hist  = np.zeros(numTimeSteps)
    discharge_exergy_hist = np.zeros(numTimeSteps)

    energy_charge_hist = []
    energy_discharge_hist = []
    stored_energy_hist = []
    charge_exergy_hist = []
    discharge_exergy_hist = []
    T_gas_inner = []

    NPoints = len(Ts)
    dr      = r[1] - r[0]

    dt_max = dr*dr/(2*alpha_gravel)
    dt      = min(dt, dt_max*0.5)

    Area   = 2.0 * np.pi * r * B
    vol_s  = dr * B * 2.0 * np.pi * r
    c_s, k_s, mu_s, rho_s = material_library_fn('rock', 20)

    dp_time_integral = 0.0
    sim_time_total   = 0.0
    peak_pump_power  = 0.0
    stored_energy = 0.0

    lastHistOutput = time

    for iStep in range(numTimeSteps):
        Ts_old = Ts.copy()

        time += dt

        Tg_r     = np.gradient(Tg, dr)
        Ts_r     = np.gradient(Ts, dr)
        rTs_r    = r * Ts_r
        dr_rTs_r = np.gradient(rTs_r, dr)
        dr_rTs_r_over_r = dr_rTs_r / r

        c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
        v_s     = M_dot / (rho_g * Area) # air velocity
        v_local = v_s / eps # local air velocity in the void fraction
        Re      = rho_g * np.abs(v_s) * Dp / mu_g # Re for the packed bed
        Pr      = c_g * mu_g / k_g # prandtl number

        Nu = ((7 - 10*eps + 5*eps**2)*(1 + 0.7*Re**0.2 * Pr**0.33)
              + (1.33 - 2.4*eps + 1.2*eps**2)*Re**0.7 * Pr**0.33) # Nusselt number
        h_sg = (k_g / Dp) * Nu # heat transfer coefficient between the gas and solid
        H_sg = 6.0 * h_sg * (1 - eps) / (alpha * Dp) # total heat transfer coefficient
        fv    = 610.0/Re + 13.0/(Re**0.11) # if M_dot != 0 else 0.0 # friction factor for the packed bed
        dp_dr = rho_g * v_s**2 * (1-eps)/(2*Dp*eps**3) * fv # if M_dot != 0 else 0.0 # pressure drop per unit length 

        # if M_dot == 0:
        #     v_local = 0.0
        #     H_sg     = 0.0

        # update gas temperature w convective heat transfer and advection
        k2 = H_sg / (c_g * rho_g * eps)
        # Tg_fwd = np.empty_like(Tg)
        # Tg_fwd[:-1] = (Tg[1:]   - Tg[:-1])/dr
        # Tg_fwd[-1]  = Tg_fwd[-2]

        # Tg_bwd = np.empty_like(Tg)
        # Tg_bwd[1:]  = (Tg[1:]   - Tg[:-1])/dr
        # Tg_bwd[0]   = Tg_bwd[1]

        # Tg_r_up = np.where(v_local > 0, Tg_bwd, Tg_fwd)

        # Tg += (k2 * (Ts - Tg) - v_local * Tg_r_up) * dt
        Tg += (k2 * (Ts - Tg) - v_local * Tg_r) * dt

        if M_dot > 0:
            Tg[0] = Tg0 # gas BC at left end of domain i
        elif M_dot < 0:
            Tg[-1] = Tg0 # gas BC at right end of domain o
        
        #Tg = np.maximum(Tg, T_ref)

        k1 = H_sg / (rho_s * c_s * (1 - eps))
        #Ts += k1 * (Tg - Ts) * dt # convective heat transfer from gas to solid
        Ts = (Ts + k1*dt*Tg)/(1 + k1*dt)

        psi_min = 0.2
        N_ramp  = 3                      
        psi     = np.ones_like(r)           

        # proportions = np.linspace(0.0, 1.0, N_ramp)
        # psi[:N_ramp] = psi_min + proportions*(1.0 - psi_min)

        # now apply the scaled diffusion:
        Ts += psi * alpha_gravel * dr_rTs_r_over_r * dt
      
        # ############################
        Ts_pre_axial = Ts.copy()
  
        Rcond_roof_pa = (
            optimized_insul_thickness/k_super_wool
            + thickness_sand      /k_sand
            + thickness_fiberglass/k_fiberglass
            + thickness_tank      /k_g_steel
        )
        Rcond_bot_pa  = thickness_tank / k_g_steel
        Rground_pa    = thickness_ground / k_concrete
        Rbot_pa       = Rcond_bot_pa + Rground_pa
 
        T_amb_roof = T_ref
        T_amb_bot  = T_ref + 15.0
 
        A_ring = 2*np.pi * r * dr                    # m2
        m_i    = vol_s * rho_s * (1 - eps)           # kg

        beta_roof = np.zeros_like(Ts)
        beta_bot  = np.zeros_like(Ts)

        # — radiative resistance at local Ts
        Rrad_roof_pa = 1.0/(
            emissivity
        * stefan_boltzmann
        * (Ts**2 + T_ref**2)
        * (Ts + T_ref)
        )
 
        Rroof_pa = Rcond_roof_pa + 1.0/(1.0/Rrad_roof_pa + 1.0/h_convection)
 
        beta_roof = A_ring/(m_i*c_s) * (1.0/Rroof_pa)
        beta_bot  = A_ring/(m_i*c_s) * (1.0/Rbot_pa)
 
        beta_tot = beta_roof + beta_bot
        T_eq     = (beta_roof*T_amb_roof + beta_bot*T_amb_bot)/beta_tot
 
        Ts = (Ts_pre_axial + beta_tot*dt*T_eq) / (1.0 + beta_tot*dt)


        ###################################

        # heat loss from outer cylinder
        beta_outer = h_surf_outer / (rho_s * c_s * (1 - eps) * dr)
        Ts[-1] -= beta_outer * (Ts[-1] - T_ref) * dt
        # heat loss from inner cylindrical boundary
        if M_dot != 0:
            beta_inner = h_surf_inner * r[0] / (2 * B * rho_s * c_s * (1 - eps) * dr)
            Ts[0]    -= beta_inner * (Ts[0] - T_ref) * dt

                  
        if iStep % 100000 == 0:
            cfl_max   = np.max(np.abs(v_local)) * dt / dr
            k1dt_max  = np.max(k1) * dt
            dTs    = Ts - Ts_old
            idx       = np.argmax(dTs)
            r_max     = r[idx]
            dTs_max   = dTs[idx]

            print(
                f"Step {iStep}, t={time:.1f}s → "
                f"CFL_max={cfl_max:.2f}, k1·dt_max={k1dt_max:.2f}, "
                f"ΔTs_max={dTs_max:.2f} K @ r={r_max:.3f} m"
            )


        if iStep % 100000 == 0:
            print(f'convective HT g>s : {np.mean(k1 * (Tg - Ts))}, th diffusion: {np.mean(alpha_gravel * dr_rTs_r_over_r * dt)}, outer loss: {np.mean(beta_outer * (Ts[-1] - T_ref) * dt)}, inner loss: {np.mean(beta_outer * (Ts[0] - T_ref) * dt)}')

        # energy input and output
        power_in_left = M_dot * (Tg[0] - T_ref) * c_g[0]
        energy_in_left += power_in_left * dt
        power_out_right = M_dot * (Tg[-1] - T_ref) * c_g[-1]
        energy_out_right += power_out_right * dt

        # pumping power
        pumping_power = np.sum(Area * np.abs(v_s) * dp_dr) * dr
        rho_in = rho_g[0]
        Q_inst = abs(M_dot) / rho_in         # m³/s
        dp_inst = pumping_power / Q_inst
        pumping_energy += pumping_power * dt
        dp_time_integral += dp_inst * dt
        peak_pump_power = max(peak_pump_power, pumping_power)

        T_gas_inner.append(Tg[0])

        # exergy integration 
        if M_dot > 0:
            Psi = 1 if Tg[0] >= MinUsefulTemp else 0
            charging_exergy = Psi * M_dot * c_g[0] * (Tg[0] - T_ref) * dt
            discharge_exergy_hist.append(0)
            charge_exergy_hist.append(charging_exergy / 3.6e6)

        elif M_dot < 0:
            Psi = 1 if Tg[0] >= MinUsefulTemp else 0
            discharging_exergy = Psi* (-M_dot) * c_g[0] * (Tg[0] - T_ref) * dt # change to -1?
            charge_exergy_hist.append(0)
            discharge_exergy_hist.append(discharging_exergy / 3.6e6)
        else:
            discharge_exergy_hist.append(0)
            charge_exergy_hist.append(0)

        # energy integration 
        P_in = M_dot * c_g[0] * (Tg[0] - T_ref)
        if M_dot > 0:
            energy_charge = P_in * dt
            energy_discharge_hist.append(0)
            energy_charge_hist.append(energy_charge / 3.6e6)

        elif M_dot < 0:
            energy_discharge = (-P_in) * dt
            energy_charge_hist.append(0)
            energy_discharge_hist.append(energy_discharge / 3.6e6)
        else:
            energy_charge_hist.append(0)
            energy_discharge_hist.append(0)


        if iStep % 100000 == 0:
              print(f"Step {iStep}/{numTimeSteps}: t = {time:.1f}s, M_dot = {M_dot:.4f}")
              print(f"c_g[0]: {c_g[0]:.3f}, Tg[0]: {Tg[0]-273.15:.3f}, T_ref: {T_ref:.3f}")
              print(f"peak pump power: {peak_pump_power:.3f}, pumping_energy: {pumping_energy:.3f}, charging exergy: {charging_exergy:.2f}, discharging exergy: {discharging_exergy:.2f}")
              
        energy_in_left   += M_dot     * (Tg[0]  - T_ref) * c_g[0]  * dt
        energy_out_right += M_dot     * (Tg[-1] - T_ref) * c_g[-1] * dt

        energy_in       = energy_in_left - energy_out_right
        stored_energy   = (
            np.sum((Ts-T_ref)*c_s*rho_s*vol_s*(1-eps))
            + np.sum((Tg-T_ref)*c_g*rho_g*vol_s*eps)
            )
        stored_energy_hist.append(stored_energy / 3.6e6)

        if (time - lastHistOutput) >= TraceTimeIncrement:
            lastHistOutput = time
            timeHistRow   += 1
            if timeHistRow >= timeHistory.shape[0]:
                timeHistory = np.vstack([
                    timeHistory,
                    np.zeros((1, timeHistory.shape[1]))
                ])
            timeHistory[timeHistRow,0] = time
            for jNode in range(numHistNodes):
                timeHistory[timeHistRow,2*jNode+1:2*jNode+3] = [
                    Tg[NodesToTraceInTime[jNode]],
                    Ts[NodesToTraceInTime[jNode]]
                ]


    dp_avg          = dp_time_integral / sim_time_total
    peak_pump_kW    = peak_pump_power / 1e3
    pumping_kWh     = pumping_energy   / 3.6e6

    return (
        Ts, Tg, timeHistory, time,
        energy_in, stored_energy_hist,
        pumping_kWh, dp_avg,
        peak_pump_kW,
        charge_exergy_hist, discharge_exergy_hist, energy_charge_hist, energy_discharge_hist, T_gas_inner
    )

# import numpy as np
# from material_library_fn import material_library_fn
# from one_D_model import calculate_axial_annular_roof_loss, calculate_bottom_annular_loss, outer_area

# def radial_thermocline(
#     Ts, Tg, HTF_type, r,
#     numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
#     h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
#     M_dot, time, dt,
#     Tg0, T_ref,
#     NodesToTraceInTime, TraceTimeIncrement,
#     timeHistory,
#     MinUsefulTemp,            # mi temperature (C) for exergy 
#     optimized_insul_thickness,
#     rR
# ):
    
#     numHistNodes    = len(NodesToTraceInTime)
#     timeHistRow     = len(timeHistory)
#     numHistOutputs  = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
#     timeHistory     = np.vstack([
#         timeHistory,
#         np.zeros((numHistOutputs, 1 + 2 * numHistNodes))
#     ])
#     MinUsefulTemp += 273.15 # C to K
#     energy_in_left    = 0.0
#     energy_out_right  = 0.0
#     pumping_energy    = 0.0
#     charging_exergy   = 0.0    
#     discharging_exergy= 0.0    

#     NPoints = len(Ts)
#     dr      = r[1] - r[0]

#     Area   = 2.0 * np.pi * r * B
#     vol_s  = dr * B * 2.0 * np.pi * r
#     c_s, k_s, mu_s, rho_s = material_library_fn('rock', 20)

#     dp_time_integral = 0.0
#     sim_time_total   = 0.0
#     peak_pump_power  = 0.0
#     stored_energy = 0.0

#     lastHistOutput = time

#     for iStep in range(numTimeSteps):
#         Ts_old = Ts.copy()

#         time += dt

#         #Tg_r     = np.gradient(Tg, dr)
#         Ts_r     = np.gradient(Ts, dr)
#         rTs_r    = r * Ts_r
#         dr_rTs_r = np.gradient(rTs_r, dr)
#         dr_rTs_r_over_r = dr_rTs_r / r

#         c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
#         v_s     = M_dot / (rho_g * Area) # air velocity
#         v_local = v_s / eps # local air velocity in the void fraction
#         Re      = rho_g * np.abs(v_s) * Dp / mu_g # Re for the packed bed
#         Pr      = c_g * mu_g / k_g # prandtl number

#         Nu = ((7 - 10*eps + 5*eps**2)*(1 + 0.7*Re**0.2 * Pr**0.33)
#               + (1.33 - 2.4*eps + 1.2*eps**2)*Re**0.7 * Pr**0.33) # Nusselt number
#         h_sg = (k_g / Dp) * Nu # heat transfer coefficient between the gas and solid
#         H_sg = 6.0 * h_sg * (1 - eps) / (alpha * Dp) # total heat transfer coefficient
#         fv    = 610.0/Re + 13.0/(Re**0.11) if M_dot != 0 else 0.0 # friction factor for the packed bed
#         dp_dr = rho_g * v_s**2 * (1-eps)/(2*Dp*eps**3) * fv if M_dot != 0 else 0.0 # pressure drop per unit length 

#         if M_dot == 0:
#             v_local = 0.0
#             H_sg     = 0.0

#         k2 = H_sg / (c_g * rho_g * eps)
#         Tg_fwd = np.empty_like(Tg)
#         Tg_fwd[:-1] = (Tg[1:]   - Tg[:-1])/dr
#         Tg_fwd[-1]  = Tg_fwd[-2]

#         Tg_bwd = np.empty_like(Tg)
#         Tg_bwd[1:]  = (Tg[1:]   - Tg[:-1])/dr
#         Tg_bwd[0]   = Tg_bwd[1]

#         Tg_r_up = np.where(v_local > 0, Tg_bwd, Tg_fwd)

#         Tg += (k2 * (Ts - Tg) - v_local * Tg_r_up) * dt
#         #Tg += (k2 * (Ts - Tg) - v_local * Tg_r) * dt

#         if M_dot > 0:
#             Tg[0] = Tg0 # gas BC at left end of domain i
#         elif M_dot < 0:
#             Tg[-1] = Tg0 # gas BC at right end of domain o
        
#         Tg = np.maximum(Tg, T_ref)

#         k1 = H_sg / (rho_s * c_s * (1 - eps))
#         #Ts += k1 * (Tg - Ts) * dt # convective heat transfer from gas to solid
#         Ts = (Ts + k1*dt*Tg)/(1 + k1*dt)
#         Ts += alpha_gravel * dr_rTs_r_over_r * dt # radial thermal diffusion in packed bed
      
#         Q_roof = calculate_axial_annular_roof_loss(temp_surface=Ts, temp_ambient=T_ref, thickness_super_wool=optimized_insul_thickness, outer_diameter=2*rR)
#         Q_bottom = calculate_bottom_annular_loss(temp_surface=Ts, temp_ambient=T_ref, outer_diameter=2*rR)
#         # remove axial heat loss from packed bed
#         mass_gravel = np.sum(vol_s * rho_s * (1 - eps))      # kg
#         # temp drop = (Q_loss * dt) / (mass * c_s)
#         deltaT_roof = (Q_roof * dt) / (mass_gravel * c_s)
#         deltaT_bot = (Q_bottom * dt) / (mass_gravel * c_s)

#         Ts -= deltaT_roof
#         Ts -= deltaT_bot

#         # k3_equiv = (Q_roof + Q_bottom) / ((T_bed_avg - T_ref) * ms·cs)
#         # Ts -= k3_equiv*(Ts - T_ref)*dt

#         # if iStep % 100000 == 0:
#         #     print(f'deltaT_roof: {deltaT_roof},Q_roof: {Q_roof},deltaT_bot: {deltaT_bot},Q_bot: {Q_bottom}, dt: {dt}, mass_gravel: {mass_gravel}, c_s: {c_s}')

#         # heat loss from outer cylinder
#         beta_outer = h_surf_outer / (rho_s * c_s * (1 - eps) * dr)
#         Ts[-1] -= beta_outer * (Ts[-1] - T_ref) * dt
#         # heat loss from inner cylindrical boundary
#         beta_inner = h_surf_inner * r[0] / (2 * B * rho_s * c_s * (1 - eps) * dr)
#         Ts[0] -= beta_inner * (Ts[0] - T_ref) * dt

                  
#         if iStep % 100000 == 0:
#             cfl_max   = np.max(np.abs(v_local)) * dt / dr
#             k1dt_max  = np.max(k1) * dt
#             dTs    = Ts - Ts_old
#             idx       = np.argmax(dTs)
#             r_max     = r[idx]
#             dTs_max   = dTs[idx]

#             print(
#                 f"Step {iStep}, t={time:.1f}s → "
#                 f"CFL_max={cfl_max:.2f}, k1·dt_max={k1dt_max:.2f}, "
#                 f"ΔTs_max={dTs_max:.2f} K @ r={r_max:.3f} m"
#             )


#         if iStep % 100000 == 0:
#             print(f'convective HT g>s : {np.mean(k1 * (Tg - Ts))}, th diffusion: {np.mean(alpha_gravel * dr_rTs_r_over_r * dt)}, deltaT_roof: {np.mean(deltaT_roof)}, deltaT_bot: {np.mean(deltaT_bot)}, outer loss: {np.mean(beta_outer * (Ts[-1] - T_ref) * dt)}, inner loss: {np.mean(beta_outer * (Ts[0] - T_ref) * dt)}')

#         # energy input and output
#         power_in_left = M_dot * (Tg[0] - T_ref) * c_g[0]
#         energy_in_left += power_in_left * dt
#         power_out_right = M_dot * (Tg[-1] - T_ref) * c_g[-1]
#         energy_out_right += power_out_right * dt

#         # pumping power
#         pumping_power = np.sum(Area * np.abs(v_s) * dp_dr) * dr
#         rho_in = rho_g[0]
#         Q_inst = abs(M_dot) / rho_in         # m³/s
#         dp_inst = pumping_power / Q_inst
#         pumping_energy += pumping_power * dt
#         dp_time_integral += dp_inst * dt
#         peak_pump_power = max(peak_pump_power, pumping_power)

#         # exergy integration 
#         if M_dot > 0:
#             Psi = 1 if Tg[0] >= MinUsefulTemp else 0
#             charging_exergy += Psi * M_dot * c_g[0] * (Tg[0] - T_ref) * dt
#         elif M_dot < 0:
#             if iStep % 100000 == 0:
#               print(f"c_g[0]: {c_g[0]:.3f}, Tg[0]: {Tg[0]:.3f}, T_ref: {T_ref:.3f}")
#             Psi = 1 if Tg[0] >= MinUsefulTemp else 0
#             discharging_exergy += Psi* (-M_dot) * c_g[0] * (Tg[0] - T_ref) * dt # change to -1?

#         if iStep % 100000 == 0:
#               print(f"Step {iStep}/{numTimeSteps}: t = {time:.1f}s, M_dot = {M_dot:.4f}")
#               print(f"peak pump power: {peak_pump_power:.3f}, pumping_energy: {pumping_energy:.3f}, charging exergy: {charging_exergy:.2f}, discharging exergy: {discharging_exergy:.2f}")
              
 
#         energy_in_left   += M_dot     * (Tg[0]  - T_ref) * c_g[0]  * dt
#         energy_out_right += M_dot     * (Tg[-1] - T_ref) * c_g[-1] * dt

#         if (time - lastHistOutput) >= TraceTimeIncrement:
#             lastHistOutput = time
#             timeHistRow   += 1
#             if timeHistRow >= timeHistory.shape[0]:
#                 timeHistory = np.vstack([
#                     timeHistory,
#                     np.zeros((1, timeHistory.shape[1]))
#                 ])
#             timeHistory[timeHistRow,0] = time
#             for jNode in range(numHistNodes):
#                 timeHistory[timeHistRow,2*jNode+1:2*jNode+3] = [
#                     Tg[NodesToTraceInTime[jNode]],
#                     Ts[NodesToTraceInTime[jNode]]
#                 ]

        
#     energy_in       = energy_in_left - energy_out_right
#     stored_energy   = (
#         np.sum((Ts-T_ref)*c_s*rho_s*vol_s*(1-eps))
#       + np.sum((Tg-T_ref)*c_g*rho_g*vol_s*eps)
#     )
#     print(f'stored_energy: {stored_energy:.2f}')


#     dp_avg          = dp_time_integral / sim_time_total
#     peak_pump_kW    = peak_pump_power / 1e3
#     pumping_kWh     = pumping_energy   / 3.6e6
#     charge_ex_kWh   = charging_exergy   / 3.6e6
#     discharge_ex_kWh= discharging_exergy/ 3.6e6

#     return (
#         Ts, Tg, timeHistory, time,
#         energy_in, stored_energy,
#         pumping_kWh, dp_avg,
#         peak_pump_kW,
#         charge_ex_kWh, discharge_ex_kWh
#     )


# ###############################
# import numpy as np
# from material_library_fn import material_library_fn
# from one_D_model import calculate_axial_annular_roof_loss, calculate_bottom_annular_loss, outer_area

# def radial_thermocline(
#     Ts, Tg, HTF_type, r,
#     numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
#     h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
#     M_dot, time, dt,
#     Tg0, T_ref,
#     NodesToTraceInTime, TraceTimeIncrement,
#     timeHistory,
#     MinUsefulTemp,            # mi temperature (°C) for exergy counting
#     optimized_insul_thickness,
#     rR
# ):
    
#     numHistNodes    = len(NodesToTraceInTime)
#     timeHistRow     = len(timeHistory)
#     numHistOutputs  = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
#     timeHistory     = np.vstack([
#         timeHistory,
#         np.zeros((numHistOutputs, 1 + 2 * numHistNodes))
#     ])
#     MinUsefulTemp += 273.15 # C to K
#     energy_in_left    = 0.0
#     energy_out_right  = 0.0
#     pumping_energy    = 0.0
#     charging_exergy   = 0.0    
#     discharging_exergy= 0.0    

#     NPoints = len(Ts)
#     dr      = r[1] - r[0]

#     Area   = 2.0 * np.pi * r * B
#     vol_s  = dr * B * 2.0 * np.pi * r
#     c_s, k_s, mu_s, rho_s = material_library_fn('rock', 20)

#     dp_time_integral = 0.0
#     sim_time_total   = 0.0
#     peak_pump_power  = 0.0
#     stored_energy = 0.0

#     lastHistOutput = time

#     for iStep in range(numTimeSteps):
#         time += dt
        
#         # radial temperature gradients
#         Tg_r     = np.gradient(Tg, dr)
#         Ts_r     = np.gradient(Ts, dr)
#         rTs_r    = r * Ts_r
#         dr_rTs_r = np.gradient(rTs_r, dr)
#         dr_rTs_r_over_r = dr_rTs_r / r

#         # fluid & convective material properties for the gas HTF
#         c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
#         v_s     = M_dot / (rho_g * Area) # air velocity
#         v_local = v_s / eps # local air velocity in the void fraction
#         Re      = rho_g * np.abs(v_s) * Dp / mu_g # Re for the packed bed
#         Pr      = c_g * mu_g / k_g # prandtl number

#         Nu = ((7 - 10*eps + 5*eps**2)*(1 + 0.7*Re**0.2 * Pr**0.33)
#               + (1.33 - 2.4*eps + 1.2*eps**2)*Re**0.7 * Pr**0.33) # Nusselt number
#         h_sg = (k_g / Dp) * Nu # heat transfer coefficient between the gas and solid
#         H_sg = 6.0 * h_sg * (1 - eps) / (alpha * Dp) # total heat transfer coefficient
#         fv    = 610.0/Re + 13.0/(Re**0.11) # friction factor for the packed bed
#         dp_dr = rho_g * v_s**2 * (1-eps)/(2*Dp*eps**3) * fv # pressure drop per unit length 

#         # update gas temp using convective heat transfer and advection
#         k2 = H_sg / (c_g * rho_g * eps)
#         Tg += (k2 * (Ts - Tg) - v_local * Tg_r) * dt
#         if M_dot > 0:
#             Tg[0] = Tg0 # gas BC at left end of domain (inlet temperature)
#         elif M_dot < 0:
#             Tg[-1] = Tg0 # gas BC at right end of domain (inlet temperature)

#         # Update solid temperature (convection + diffusion)
#         k1 = H_sg / (rho_s * c_s * (1 - eps))
#         Ts += k1 * (Tg - Ts) * dt # convective heat transfer from gas to solid
#         Ts += alpha_gravel * dr_rTs_r_over_r * dt # radial thermal diffusion in packed bed

#         # Axial heat loss - conductive/convective heat transfer to environment (top and bottom surfaces)
#         # precompute once per time‐step:
#         areas   = 2*np.pi * r * dr                    # m² of each ring
#         vol_s_i  = vol_s                              # m³ of gravel in each ring
#         mass_i   = vol_s_i * rho_s * (1 - eps)        # kg of gravel in each ring
#         A_total  = areas.sum()                        # total annular area
#         dt_s     = dt

#         # compute local roof losses
#         # - each ring “sees” its local Ts[i]
#         # - total Q_roof is distributed proportionally by area_i/A_total
#         Q_roof_total = calculate_axial_annular_roof_loss(
#             temp_surface=Ts,           # *vector* of length NPoints
#             temp_ambient=T_ref,
#             thickness_super_wool=optimized_insul_thickness,
#             outer_diameter=2*rR
#         )
#         # Q_roof_total now a vector of W for each ring

#         # temperature drop in each ring
#         # ΔT_i = (Q_roof_i * dt)/(m_i * c_s)
#         deltaT_roof = Q_roof_total * dt_s / (mass_i * c_s)

#         # same for bottom
#         Q_bot_total = calculate_bottom_annular_loss(
#             temp_surface=Ts,
#             temp_ambient=T_ref,
#             outer_diameter=2*rR
#         )
#         deltaT_bot = Q_bot_total * dt_s / (mass_i * c_s)

#         # then subtract locally
#         Ts -= deltaT_roof
#         Ts -= deltaT_bot

#         # if iStep % 100000 == 0:
#         #     print(f'deltaT_roof: {deltaT_roof},deltaT_bot: {deltaT_bot}')


#         # Radial loss through outer and inner cylinder 
#         # heat loss from outer cylindrical boundary
#         beta_outer = h_surf_outer / (rho_s * c_s * (1 - eps) * dr)
#         Ts[-1] -= beta_outer * (Ts[-1] - T_ref) * dt
#         # heat loss from inner cylindrical boundary
#         beta_inner = h_surf_inner * r[0] / (2 * B * rho_s * c_s * (1 - eps) * dr)
#         Ts[0] -= beta_inner * (Ts[0] - T_ref) * dt

#         # energy input and output
#         power_in_left = M_dot * (Tg[0] - T_ref) * c_g[0]
#         energy_in_left += power_in_left * dt
#         power_out_right = M_dot * (Tg[-1] - T_ref) * c_g[-1]
#         energy_out_right += power_out_right * dt

#         # pumping power
#         pumping_power = np.sum(Area * np.abs(v_s) * dp_dr) * dr
#         rho_in = rho_g[0]
#         Q_inst = abs(M_dot) / rho_in         # m³/s
#         dp_inst = pumping_power / Q_inst
#         pumping_energy += pumping_power * dt
#         dp_time_integral += dp_inst * dt
#         peak_pump_power = max(peak_pump_power, pumping_power)

#         # exergy integration 
#         if M_dot > 0:
            
#             Psi = 1 if Tg[0] >= MinUsefulTemp else 0
#             charging_exergy += Psi * M_dot * c_g[0] * (Tg[0] - T_ref) * dt
#         elif M_dot < 0:
#             Psi = 1 if Tg[0] >= MinUsefulTemp else 0
#             discharging_exergy += Psi* (-M_dot) * c_g[0] * (Tg[0] - T_ref) * dt # change to -1?

#         if iStep % 100000 == 0:
#               print(f"Step {iStep}/{numTimeSteps}: t = {time:.1f}s, M_dot = {M_dot:.4f}")
#               print(f"peak pump power: {peak_pump_power:.3f}, pumping_energy: {pumping_energy:.3f}, charging exergy: {charging_exergy:.2f}, discharging exergy: {discharging_exergy:.2f}")
              
 
#         energy_in_left   += M_dot     * (Tg[0]  - T_ref) * c_g[0]  * dt
#         energy_out_right += M_dot     * (Tg[-1] - T_ref) * c_g[-1] * dt

#         if (time - lastHistOutput) >= TraceTimeIncrement:
#             lastHistOutput = time
#             timeHistRow   += 1
#             if timeHistRow >= timeHistory.shape[0]:
#                 timeHistory = np.vstack([
#                     timeHistory,
#                     np.zeros((1, timeHistory.shape[1]))
#                 ])
#             timeHistory[timeHistRow,0] = time
#             for jNode in range(numHistNodes):
#                 timeHistory[timeHistRow,2*jNode+1:2*jNode+3] = [
#                     Tg[NodesToTraceInTime[jNode]],
#                     Ts[NodesToTraceInTime[jNode]]
#                 ]

        
#     energy_in       = energy_in_left - energy_out_right
#     stored_energy   = (
#         np.sum((Ts-T_ref)*c_s*rho_s*vol_s*(1-eps))
#       + np.sum((Tg-T_ref)*c_g*rho_g*vol_s*eps)
#     )
#     print(f'stored_energy: {stored_energy:.2f}')


#     dp_avg          = dp_time_integral / sim_time_total
#     peak_pump_kW    = peak_pump_power / 1e3
#     pumping_kWh     = pumping_energy   / 3.6e6
#     charge_ex_kWh   = charging_exergy   / 3.6e6
#     discharge_ex_kWh= discharging_exergy/ 3.6e6

#     return (
#         Ts, Tg, timeHistory, time,
#         energy_in, stored_energy,
#         pumping_kWh, dp_avg,
#         peak_pump_kW,
#         charge_ex_kWh, discharge_ex_kWh
#     )


# ########################

# import numpy as np
# from material_library_fn import material_library_fn
# from one_D_model import calculate_axial_annular_roof_loss, calculate_bottom_annular_loss, outer_area

# # constants
# stefan_boltzmann = 5.67e-8  # W/(m²·K⁴)

# # Material thermal conductivities (W/(m·K))
# k_super_wool   = 0.19   # Super‑wool insulation
# k_sand         = 0.114  # Sand
# k_fiberglass   = 0.04   # Fiberglass blanket
# k_g_steel      = 50     # Galvanized steel tank

# # Layer thicknesses (m)
# # super wool thickness optimized in main  
# thickness_sand        = 0.0254    # Sand
# thickness_fiberglass  = 0.0047625 # Fiberglass
# thickness_tank        = 0.0047625 # Steel tank

# # ground pad under tank
# k_concrete      = 1.2    # W/m·K
# thickness_ground= 0.5    # m

# # surface properties
# emissivity      = 0.05   # –
# h_convection    = 5.0    # W/(m²·K)

# # Geometry
# inner_diameter  = 2*0.2286   # central hole r (m)

# def radial_thermocline(
#     Ts, Tg, HTF_type, r,
#     numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
#     h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
#     M_dot, time, dt,
#     Tg0, T_ref,
#     NodesToTraceInTime, TraceTimeIncrement,
#     timeHistory,
#     MinUsefulTemp,            # mi temperature (C) for exergy 
#     optimized_insul_thickness,
#     rR
# ):
    
#     numHistNodes    = len(NodesToTraceInTime)
#     timeHistRow     = len(timeHistory)
#     numHistOutputs  = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
#     timeHistory     = np.vstack([
#         timeHistory,
#         np.zeros((numHistOutputs, 1 + 2 * numHistNodes))
#     ])
#     MinUsefulTemp += 273.15 # C to K
#     energy_in_left    = 0.0
#     energy_out_right  = 0.0
#     pumping_energy    = 0.0
#     charging_exergy   = 0.0    
#     discharging_exergy= 0.0    

#     NPoints = len(Ts)
#     dr      = r[1] - r[0]

#     Area   = 2.0 * np.pi * r * B
#     vol_s  = dr * B * 2.0 * np.pi * r
#     c_s, k_s, mu_s, rho_s = material_library_fn('rock', 20)

#     dp_time_integral = 0.0
#     sim_time_total   = 0.0
#     peak_pump_power  = 0.0
#     stored_energy = 0.0

#     lastHistOutput = time

#     for iStep in range(numTimeSteps):
#         Ts_old = Ts.copy()

#         time += dt

#         #Tg_r     = np.gradient(Tg, dr)
#         Ts_r     = np.gradient(Ts, dr)
#         rTs_r    = r * Ts_r
#         dr_rTs_r = np.gradient(rTs_r, dr)
#         dr_rTs_r_over_r = dr_rTs_r / r

#         c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
#         v_s     = M_dot / (rho_g * Area) # air velocity
#         v_local = v_s / eps # local air velocity in the void fraction
#         Re      = rho_g * np.abs(v_s) * Dp / mu_g # Re for the packed bed
#         Pr      = c_g * mu_g / k_g # prandtl number

#         Nu = ((7 - 10*eps + 5*eps**2)*(1 + 0.7*Re**0.2 * Pr**0.33)
#               + (1.33 - 2.4*eps + 1.2*eps**2)*Re**0.7 * Pr**0.33) # Nusselt number
#         h_sg = (k_g / Dp) * Nu # heat transfer coefficient between the gas and solid
#         H_sg = 6.0 * h_sg * (1 - eps) / (alpha * Dp) # total heat transfer coefficient
#         fv    = 610.0/Re + 13.0/(Re**0.11) if M_dot != 0 else 0.0 # friction factor for the packed bed
#         dp_dr = rho_g * v_s**2 * (1-eps)/(2*Dp*eps**3) * fv if M_dot != 0 else 0.0 # pressure drop per unit length 

#         if M_dot == 0:
#             v_local = 0.0
#             H_sg     = 0.0

#         k2 = H_sg / (c_g * rho_g * eps)
#         Tg_fwd = np.empty_like(Tg)
#         Tg_fwd[:-1] = (Tg[1:]   - Tg[:-1])/dr
#         Tg_fwd[-1]  = Tg_fwd[-2]

#         Tg_bwd = np.empty_like(Tg)
#         Tg_bwd[1:]  = (Tg[1:]   - Tg[:-1])/dr
#         Tg_bwd[0]   = Tg_bwd[1]

#         Tg_r_up = np.where(v_local > 0, Tg_bwd, Tg_fwd)

#         Tg += (k2 * (Ts - Tg) - v_local * Tg_r_up) * dt
#         #Tg += (k2 * (Ts - Tg) - v_local * Tg_r) * dt

#         if M_dot > 0:
#             Tg[0] = Tg0 # gas BC at left end of domain i
#         elif M_dot < 0:
#             Tg[-1] = Tg0 # gas BC at right end of domain o
        
#         Tg = np.maximum(Tg, T_ref)

#         k1 = H_sg / (rho_s * c_s * (1 - eps))
#         #Ts += k1 * (Tg - Ts) * dt # convective heat transfer from gas to solid
#         Ts = (Ts + k1*dt*Tg)/(1 + k1*dt)
#         Ts += alpha_gravel * dr_rTs_r_over_r * dt # radial thermal diffusion in packed bed
      
        
#         # remove axial heat loss from packed bed

#         # 1) your pure conduction resistances per unit area (scalars)
#         Rcond_roof_pa = ( optimized_insul_thickness/k_super_wool
#                         + thickness_sand      /k_sand
#                         + thickness_fiberglass/k_fiberglass
#                         + thickness_tank      /k_g_steel )
#         Rcond_bot_pa  = thickness_tank / k_g_steel
#         Rground_pa    = thickness_ground / k_concrete

#         # 2) the radiative + convective resistances at the roof surface (arrays)
#         Rrad_roof_pa = 1.0 / ( emissivity
#                             * stefan_boltzmann
#                             * (Ts**2 + T_ref**2)
#                             * (Ts + T_ref) )
#         Rconv_pa     = 1.0 / h_convection

#         # 3) assemble the total roof resistance per area (array) and bottom (scalar)
#         Rroof_tot_pa = Rcond_roof_pa + 1.0/(1.0/Rrad_roof_pa + 1.0/Rconv_pa)
#         Rbot_tot_pa  = Rcond_bot_pa   + Rground_pa

#         # 1) compute each ring’s area and mass
#         A_ring = 2*np.pi * r * dr                  # m²
#         m_i    = vol_s * rho_s * (1 - eps)         # kg
#         m_i[0] = m_i[1]

#         # 2) compute ring‐wise heat fluxes (W/m²)
#         qpp_roof = (Ts - T_ref)     / Rroof_tot_pa   # array, roof
#         qpp_bot  = (Ts - (T_ref+25)) / Rbot_tot_pa    # array, bottom

#         # 3) convert to ring‐wise total W, then ΔT
#         Q_roof_ring = qpp_roof * A_ring             # W per ring
#         Q_bot_ring  = qpp_bot  * A_ring             # W per ring

#         T_amb_roof = T_ref           # air/roof ambient
#         T_amb_bot  = T_ref + 25.0    # ground ambient (25 K above air, or whatever you choose)

#         β_roof = A_ring/(m_i*c_s) * (1.0/ Rroof_tot_pa)   # array, s⁻¹
#         β_bot  = A_ring/(m_i*c_s) * (1.0/ Rbot_tot_pa)    # array, s⁻¹
#         β_tot  = β_roof + β_bot                           # array, s⁻¹

#         Ts = ( Ts
#             + dt*( β_roof*T_amb_roof
#                 + β_bot *T_amb_bot )
#             ) / (1.0 + dt*β_tot)


#         # k3_equiv = (Q_roof + Q_bottom) / ((T_bed_avg - T_ref) * ms·cs)
#         # Ts -= k3_equiv*(Ts - T_ref)*dt

#         # if iStep % 100000 == 0:
#         #     print(f'deltaT_roof: {deltaT_roof},Q_roof: {Q_roof},deltaT_bot: {deltaT_bot},Q_bot: {Q_bottom}, dt: {dt}, mass_gravel: {mass_gravel}, c_s: {c_s}')

#         # heat loss from outer cylinder
#         beta_outer = h_surf_outer / (rho_s * c_s * (1 - eps) * dr)
#         Ts[-1] -= beta_outer * (Ts[-1] - T_ref) * dt
#         # heat loss from inner cylindrical boundary
#         beta_inner = h_surf_inner * r[0] / (2 * B * rho_s * c_s * (1 - eps) * dr)
#         Ts[0] -= beta_inner * (Ts[0] - T_ref) * dt

                  
#         if iStep % 100000 == 0:
#             cfl_max   = np.max(np.abs(v_local)) * dt / dr
#             k1dt_max  = np.max(k1) * dt
#             dTs    = Ts - Ts_old
#             idx       = np.argmax(dTs)
#             r_max     = r[idx]
#             dTs_max   = dTs[idx]

#             print(
#                 f"Step {iStep}, t={time:.1f}s → "
#                 f"CFL_max={cfl_max:.2f}, k1·dt_max={k1dt_max:.2f}, "
#                 f"ΔTs_max={dTs_max:.2f} K @ r={r_max:.3f} m"
#             )


#         if iStep % 100000 == 0:
#             print(f'convective HT g>s : {np.mean(k1 * (Tg - Ts))}, th diffusion: {np.mean(alpha_gravel * dr_rTs_r_over_r * dt)}, outer loss: {np.mean(beta_outer * (Ts[-1] - T_ref) * dt)}, inner loss: {np.mean(beta_outer * (Ts[0] - T_ref) * dt)}')

#         # energy input and output
#         power_in_left = M_dot * (Tg[0] - T_ref) * c_g[0]
#         energy_in_left += power_in_left * dt
#         power_out_right = M_dot * (Tg[-1] - T_ref) * c_g[-1]
#         energy_out_right += power_out_right * dt

#         # pumping power
#         pumping_power = np.sum(Area * np.abs(v_s) * dp_dr) * dr
#         rho_in = rho_g[0]
#         Q_inst = abs(M_dot) / rho_in         # m³/s
#         dp_inst = pumping_power / Q_inst
#         pumping_energy += pumping_power * dt
#         dp_time_integral += dp_inst * dt
#         peak_pump_power = max(peak_pump_power, pumping_power)

#         # exergy integration 
#         if M_dot > 0:
#             Psi = 1 if Tg[0] >= MinUsefulTemp else 0
#             charging_exergy += Psi * M_dot * c_g[0] * (Tg[0] - T_ref) * dt
#         elif M_dot < 0:
#             if iStep % 100000 == 0:
#               print(f"c_g[0]: {c_g[0]:.3f}, Tg[0]: {Tg[0]:.3f}, T_ref: {T_ref:.3f}")
#             Psi = 1 if Tg[0] >= MinUsefulTemp else 0
#             discharging_exergy += Psi* (-M_dot) * c_g[0] * (Tg[0] - T_ref) * dt # change to -1?

#         if iStep % 100000 == 0:
#               print(f"Step {iStep}/{numTimeSteps}: t = {time:.1f}s, M_dot = {M_dot:.4f}")
#               print(f"peak pump power: {peak_pump_power:.3f}, pumping_energy: {pumping_energy:.3f}, charging exergy: {charging_exergy:.2f}, discharging exergy: {discharging_exergy:.2f}")
              
 
#         energy_in_left   += M_dot     * (Tg[0]  - T_ref) * c_g[0]  * dt
#         energy_out_right += M_dot     * (Tg[-1] - T_ref) * c_g[-1] * dt

#         if (time - lastHistOutput) >= TraceTimeIncrement:
#             lastHistOutput = time
#             timeHistRow   += 1
#             if timeHistRow >= timeHistory.shape[0]:
#                 timeHistory = np.vstack([
#                     timeHistory,
#                     np.zeros((1, timeHistory.shape[1]))
#                 ])
#             timeHistory[timeHistRow,0] = time
#             for jNode in range(numHistNodes):
#                 timeHistory[timeHistRow,2*jNode+1:2*jNode+3] = [
#                     Tg[NodesToTraceInTime[jNode]],
#                     Ts[NodesToTraceInTime[jNode]]
#                 ]

        
#     energy_in       = energy_in_left - energy_out_right
#     stored_energy   = (
#         np.sum((Ts-T_ref)*c_s*rho_s*vol_s*(1-eps))
#       + np.sum((Tg-T_ref)*c_g*rho_g*vol_s*eps)
#     )
#     print(f'stored_energy: {stored_energy:.2f}')


#     dp_avg          = dp_time_integral / sim_time_total
#     peak_pump_kW    = peak_pump_power / 1e3
#     pumping_kWh     = pumping_energy   / 3.6e6
#     charge_ex_kWh   = charging_exergy   / 3.6e6
#     discharge_ex_kWh= discharging_exergy/ 3.6e6

#     return (
#         Ts, Tg, timeHistory, time,
#         energy_in, stored_energy,
#         pumping_kWh, dp_avg,
#         peak_pump_kW,
#         charge_ex_kWh, discharge_ex_kWh
#     )