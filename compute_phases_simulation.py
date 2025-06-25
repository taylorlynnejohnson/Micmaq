import numpy as np
import matplotlib.pyplot as plt

def compute_phases(GEOM_type, rL, rR, B, HTF_type, alpha, Dp, eps, alpha_gravel,
                   h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
                   simHist, time, runName, tracePositions, TraceTimeIncrement):
    # Simulate and plot temperature evolution in packed bed
    phaseStartTime = 0
    iPhase = 0
    A_left = 2 * np.pi * rL * B
    A_right = 2 * np.pi * rR * B
    A_avg = (A_left + A_right) / 2
    numPhases = simHist.shape[0]
    M_dot = simHist[0, 0]  # initial mass flow rate of gas, kg/s
    T_g0 = simHist[0, 2]  # initial gas entry temperature (C)
    T_s0 = simHist[0, 3]  # initial temperature of solid (C)
    T_ref = simHist[0, 3]  # initial reference temperature (C)
    Tavg = (T_g0 + T_s0) / 2

    estThermoclineLength, estSuperficialVelocity, estThermoclineVelocity = estimate_thermocline(
        M_dot, A_avg, eps, alpha, Dp, Tavg, HTF_type)
    dr = estThermoclineLength / 12  # segment length
    LengthDomain = rR - rL
    NSegs = int(np.floor(LengthDomain / dr))  # number of spatial segments
    NPoints = NSegs + 1  # number of spatial points
    dt = 0.025 * (dr / estSuperficialVelocity)

    # Initial Temperature Distributions (Celsius)
    Tg = np.full(NPoints, T_s0)
    Ts = np.full(NPoints, T_s0)
    Tg[0] = T_g0
    Ts[0] = T_g0
    r = np.linspace(rL, rR, NPoints)  # coordinates of nodes

    NodesToTraceInTime = []
    rNodesToTraceInTime = []
    timeHistory = [time]
    for Position in tracePositions:
        iNode = np.argmin(np.abs(r - Position))
        NodesToTraceInTime.append(iNode)
        rNodesToTraceInTime.append(r[iNode])
        timeHistory.extend([Tg[iNode], Ts[iNode]])

    # Plot settings
    plt.figure(figsize=(12, 8))
    plt.title(f"{GEOM_type} geometry with {HTF_type} heat transfer fluid")
    plt.xlabel("Radius (m)")
    plt.ylabel("Temperature (C)")
    plt.grid(True)

    # Simulate Phases
    while iPhase < numPhases:
        iPhase += 1
        M_dot = simHist[iPhase - 1, 0]
        duration = 3600 * simHist[iPhase - 1, 1]
        T_g0 = simHist[iPhase - 1, 2]  # gas temperature at entry point
        T_ref = simHist[iPhase - 1, 3]  # environmental temperature
        numSteps = int(np.floor(duration / dt))
        phaseStartTime += numSteps * dt

        if M_dot != 0:
            if M_dot > 0:
                inletTemp = T_g0
                lineColor = 'r'
            else:
                inletTemp = T_g0
                lineColor = 'b'

            if GEOM_type == 'Radial':
                Ts, Tg, timeHistory, time, energy_in, stored_energy, pumping_energy = radial_thermocline(
                    Ts, Tg, HTF_type, r, numSteps, alpha, Dp, eps, B, alpha_gravel,
                    h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
                    M_dot, time, dt, T_g0, T_ref, NodesToTraceInTime, TraceTimeIncrement, timeHistory)
            else:
                # Axial and spherical cases need to be implemented separately
                pass
        else:  # M_dot == 0
            inletTemp = T_s0
            lineColor = 'g'
            dt1 = 100 * dt
            Ts, Tg, timeHistory, time, energy_in, stored_energy, pumping_energy = transient_thermocline_degradation(
                Ts, Tg, T_ref, 'Radial', r, B, eps, alpha_gravel,
                h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
                duration, dt1, time, NodesToTraceInTime, TraceTimeIncrement, timeHistory)

        plt.plot(r, Ts, '-', linewidth=2, color=lineColor, linestyle='-')
        fileName = f"{runName}-{iPhase}"
        np.savez(fileName, Ts=Ts, Tg=Tg, timeHistory=timeHistory)
        plt.draw()
        plt.pause(0.01)

    plt.show()


def estimate_thermocline(M_dot, A_avg, eps, alpha, Dp, Tavg, HTF_type):
    # Placeholder function to estimate thermocline length, superficial velocity, and thermocline velocity
    estThermoclineLength = 10.0  # Example value
    estSuperficialVelocity = 0.1  # Example value
    estThermoclineVelocity = 0.05  # Example value
    return estThermoclineLength, estSuperficialVelocity, estThermoclineVelocity

def radial_thermocline(Ts, Tg, HTF_type, r, numSteps, alpha, Dp, eps, B, alpha_gravel,
                       h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
                       M_dot, time, dt, T_g0, T_ref, NodesToTraceInTime, TraceTimeIncrement, timeHistory):
    # Placeholder function for radial thermocline simulation
    # Implement the actual logic here or use previously converted code
    return Ts, Tg, timeHistory, time, 0, 0, 0

def transient_thermocline_degradation(Ts, Tg, T_ref, GEOM_type, r, B, eps, alpha_gravel,
                                      h_surf_upper, h_surf_lower, h_surf_outer, h_surf_inner,
                                      duration, dt1, time, NodesToTraceInTime, TraceTimeIncrement, timeHistory):
    # Placeholder function for thermocline degradation simulation
    # Implement the actual logic here
    return Ts, Tg, timeHistory, time, 0, 0, 0


