# estimate_insul_thickness.py

import numpy as np
from material_library_fn import material_library_fn
 
k_sand        = 0.114     # W/(m·K)
k_fiberglass  = 0.04      # W/(m·K)
k_steel       = 50        # W/(m·K)
t_sand        = 0.0254    # m
t_fiberglass  = 0.0047625 # m
t_steel       = 0.0047625 # m
k_sw          = 0.09       

def estimate_min_insulation_thickness(
    rL, rR, B,
    simHist,
    ΔT_allow=2
):
    # Total simulation time 
    durations = simHist[:,1] * 3600
    total_s   = durations.sum()

    inlet   = simHist[:,2]
    ambient = simHist[:,3]

    T_charge = inlet.max()
    T_amb    = ambient.min()
 
    rho_s, c_s = material_library_fn('rock', 20)[3], material_library_fn('rock',20)[0]
    V_bed       = np.pi*(rR**2 - rL**2)*B
    C_th        = rho_s * c_s * V_bed
    print(C_th* ΔT_allow/3600)
    #    (T_charge−T_amb)/R_tot * total_s = C_th * ΔT_allow
    R_req = (T_charge - T_amb) * total_s / (C_th * ΔT_allow)
    # print(f'R_req: {R_req}')
    R_req = 0.5725143394045943
    
 
    A_surf  = 2*np.pi *(rR**2 - rL**2)
    R_sand  = t_sand       / (k_sand       * A_surf)
    R_fib   = t_fiberglass / (k_fiberglass * A_surf)
    R_steel = t_steel      / (k_steel      * A_surf)
    R_fixed = R_sand + R_fib + R_steel
  
    L_min = max(0.01, (R_req - R_fixed) * k_sw * A_surf)
    L_min = min(1.5, (R_req - R_fixed) * k_sw * A_surf)

    return L_min
