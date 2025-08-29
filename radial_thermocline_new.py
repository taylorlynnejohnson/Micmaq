
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

# Geometry
inner_diameter  = 2*0.2286   # central hole r (m)

# Surface parameters
emissivity      = 0.05   # Surface emissivity
h_convection    = 5    # Convective coefficient top annulus (W/(m²·K))

k_concrete      = 1.2    # W/mK
thickness_ground = 0.5   # m



def radial_thermocline(
    Ts, Tg, HTF_type, r,
    numTimeSteps, alpha, Dp, eps, B, alpha_gravel,
    h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
    M_dot, time, dt,
    Tg0, T_ref,
    NodesToTraceInTime, TraceTimeIncrement,
    timeHistory,
    MinUsefulTemp,            # min temperature (°C) for exergy counting
    optimized_insul_thickness,
    rR
):
    
    Rcond_roof_pa = optimized_insul_thickness / (k_super_wool * (2*np.pi*(rR/2)**2 - np.pi*(inner_diameter/2)**2))
    Rcond_bot_pa  = thickness_tank       / (k_g_steel    * (2*np.pi*(rR/2)**2 - np.pi*(inner_diameter/2)**2))
    Rground_pa    = thickness_ground   / (k_concrete   * (2*np.pi*(rR/2)**2 - np.pi*(inner_diameter/2)**2))
    
    numHistNodes    = len(NodesToTraceInTime)
    timeHistRow     = len(timeHistory)
    numHistOutputs  = int(np.floor(numTimeSteps * dt / TraceTimeIncrement))
    timeHistory     = np.vstack([
        timeHistory,
        np.zeros((numHistOutputs, 1 + 2 * numHistNodes))
    ])
    MinUsefulTemp += 273.15 # C to K
    energy_in_left    = 0.0
    energy_out_right  = 0.0
    pumping_energy    = 0.0
    charging_exergy   = 0.0    
    discharging_exergy= 0.0    

    NPoints = len(Ts)
    dr      = r[1] - r[0]

    Area   = 2.0 * np.pi * r * B
    vol_s  = dr * B * 2.0 * np.pi * r
    c_s, k_s, mu_s, rho_s = material_library_fn('rock', Ts-273.15)

    dp_time_integral = 0.0
    sim_time_total   = 0.0
    peak_pump_power  = 0.0
    stored_energy = 0.0

    lastHistOutput = time


    for iStep in range(numTimeSteps):
        time += dt
        
        # radial temperature gradients
        Tg_r     = np.gradient(Tg, dr)
        Ts_r     = np.gradient(Ts, dr)
        rTs_r    = r * Ts_r
        dr_rTs_r = np.gradient(rTs_r, dr)
        dr_rTs_r_over_r = dr_rTs_r / r

        # fluid & convective material properties for the gas HTF

        c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
        v_s     = M_dot / (rho_g * Area) # air velocity
        v_local = v_s / eps # local air velocity in the void fraction
        Re      = rho_g * np.abs(v_s) * Dp / mu_g # Re for the packed bed
        Pr      = c_g * mu_g / k_g # prandtl number

        Nu = ((7 - 10*eps + 5*eps**2)*(1 + 0.7*Re**0.2 * Pr**0.33)
              + (1.33 - 2.4*eps + 1.2*eps**2)*Re**0.7 * Pr**0.33) # Nusselt number
        h_sg = (k_g / Dp) * Nu # heat transfer coefficient between the gas and solid
        H_sg = 6.0 * h_sg * (1 - eps) / (alpha * Dp) # total heat transfer coefficient
        fv    = 610.0/Re + 13.0/(Re**0.11) if M_dot != 0 else 0.0 # friction factor for the packed bed
        dp_dr = rho_g * v_s**2 * (1-eps)/(2*Dp*eps**3) * fv if M_dot != 0 else 0.0 # pressure drop per unit length 

        # update gas temp using convective heat transfer and advection
        k2 = H_sg / (c_g * rho_g * eps)
        Tg += (k2 * (Ts - Tg) - v_local * Tg_r) * dt
        if M_dot >= 0:
            Tg[0] = Tg0 # gas BC at left end of domain (inlet temperature)
        else:
            Tg[-1] = Tg0 # gas BC at right end of domain (inlet temperature)

        Tg = np.maximum(Tg, T_ref)
        # Update solid temperature (convection + diffusion)
        k1 = H_sg / (rho_s * c_s * (1 - eps))
        Ts += k1 * (Tg - Ts) * dt # convective heat transfer from gas to solid
        Ts += alpha_gravel * dr_rTs_r_over_r * dt # radial thermal diffusion in packed bed

        # Axial heat loss - conductive/convective heat transfer to environment (top and bottom surfaces)
        A_ring = 2 * np.pi * r * dr                             #  (NPoints,)
        m_i    = vol_s * rho_s * (1 - eps)                      #  (NPoints,)

        Rrad_roof_pa = 1.0 / (
            emissivity * stefan_boltzmann *
            (Ts**2 + T_ref**2) * (Ts + T_ref)
        )                                                    # (NPoints,)
        Rconv_pa = 1.0 / h_convection                           

        Rroof_tot_pa = Rcond_roof_pa + 1.0 / (1.0/Rrad_roof_pa + 1.0/Rconv_pa)

        Rbot_tot_pa = Rcond_bot_pa + Rground_pa                 

        qpp_roof = (Ts - T_ref) / Rroof_tot_pa                   # W/m2
        qpp_bot  = (Ts - (T_ref + 25)) / Rbot_tot_pa             # W/m2

        Q_roof_ring = qpp_roof * A_ring                         # W
        Q_bot_ring  = qpp_bot  * A_ring                         # W

        deltaT = (Q_roof_ring + Q_bot_ring) * dt / (m_i * c_s)   # K

        Ts -= deltaT

        # Radial loss through outer and inner cylinder 
        # heat loss from outer cylindrical boundary
        # beta_outer = h_surf_outer / (rho_s * c_s * (1 - eps) * dr)
        # Ts[-1] -= beta_outer * (Ts[-1] - T_ref) * dt
        beta_outer = h_surf_outer / (rho_s * c_s * (1 - eps) * dr)
        Ts_old = Ts[-1]
        β_out = beta_outer[-1]
        Ts[-1] = (Ts_old + β_out * dt * T_ref) / (1 + β_out * dt) # switched to explicit to avoid undershoot
        # heat loss from inner cylindrical boundary
        
        if M_dot != 0:
            beta_inner = h_surf_inner * r[0] / (2 * B * rho_s * c_s * (1 - eps) * dr)
            β_in = beta_inner[0]
            Ts[0] = (Ts_old + β_in*dt*T_ref) / (1 + β_in*dt) # samd ^ 
            if iStep % 100000 == 0:
                print(f'convective HT gas : {np.mean(k1 * (Tg - Ts))}, th diffusion: {np.mean(alpha_gravel * dr_rTs_r_over_r * dt)}, outer loss: {np.mean(beta_outer * (Ts[-1] - T_ref) * dt)}, inner loss: {np.mean(beta_outer * (Ts[0] - T_ref) * dt)}')


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

        # exergy integration 
        if M_dot > 0:
            
            Psi = 1 if Tg[0] >= MinUsefulTemp else 0
            charging_exergy += Psi * M_dot * c_g[0] * (Tg[0] - T_ref) * dt
        elif M_dot < 0:
            Psi = 1 if Tg[0] >= MinUsefulTemp else 0
            discharging_exergy += Psi* (-M_dot) * c_g[0] * (Tg[0] - T_ref) * dt # change to -1?

        if iStep % 100000 == 0:
              print(f"Step {iStep}/{numTimeSteps}: t = {time:.1f}s, M_dot = {M_dot:.4f}")
              print(f"peak pump power: {peak_pump_power:.3f}, pumping_energy: {pumping_energy:.3f}, charging exergy: {charging_exergy:.2f}, discharging exergy: {discharging_exergy:.2f}")
              
 
        energy_in_left   += M_dot     * (Tg[0]  - T_ref) * c_g[0]  * dt
        energy_out_right += M_dot     * (Tg[-1] - T_ref) * c_g[-1] * dt

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

        
    energy_in       = energy_in_left - energy_out_right
    stored_energy   = (
        np.sum((Ts-T_ref)*c_s*rho_s*vol_s*(1-eps))
      + np.sum((Tg-T_ref)*c_g*rho_g*vol_s*eps)
    )
    print(f'stored_energy: {stored_energy:.2f}')


    dp_avg          = dp_time_integral / sim_time_total
    peak_pump_kW    = peak_pump_power / 1e3
    pumping_kWh     = pumping_energy   / 3.6e6
    charge_ex_kWh   = charging_exergy   / 3.6e6
    discharge_ex_kWh= discharging_exergy/ 3.6e6

    return (
        Ts, Tg, timeHistory, time,
        energy_in, stored_energy,
        pumping_kWh, dp_avg,
        peak_pump_kW,
        charge_ex_kWh, discharge_ex_kWh
    )