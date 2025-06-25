# one_D_model.py

import numpy as np

# Physical constants
stefan_boltzmann = 5.67e-8  # W/(m²·K⁴)

# Material thermal conductivities (W/(m·K))
k_super_wool   = 0.09   # Super‑wool insulation
k_sand         = 0.114  # Sand
k_fiberglass   = 0.04   # Fiberglass blanket
k_g_steel      = 50     # Galvanized steel tank

# Layer thicknesses (m)
thickness_super_wool  = .15       # Super‑wool
thickness_sand        = 0.0254    # Sand
thickness_fiberglass  = 0.0047625 # Fiberglass
thickness_tank        = 0.0047625 # Steel tank

# Geometry
diameter     = 1.8288   # Packed‐bed outer diameter (m)
hole_radius  = 0.0254   # Central hole radius (m)
outer_area   = np.pi * (diameter/2)**2 - np.pi * hole_radius**2

# Surface parameters
emissivity      = 0.05   # Surface emissivity
h_convection    = 3.0    # Convective coefficient (W/(m²·K))


def thermal_resistance(thickness, k, area):
    # Series conduction resistance: R = L/(k·A)
    return thickness / (k * area)


def convective_resistance(h, area):
    # Convection resistance: R = 1/(h·A)
    return 1.0 / (h * area)


def radiative_resistance(T_surf, T_amb, emissivity, area):
    # Linearized radiative resistance at surface T_surf to ambient T_amb
    Ts_k = T_surf + 273.15
    Ta_k = T_amb  + 273.15
    # From dQ/dT ≈ ε·σ·A·(T⁴_surf − T⁴_amb)/(T_surf − T_amb)
    denom = emissivity * stefan_boltzmann * area * (
        (Ts_k**2 + Ta_k**2) * (Ts_k + Ta_k)
    )
    return 1 / denom


def heat_loss(temp_surface, temp_ambient, total_resistance, emissivity, area, h):
    # Radiative + convective in parallel
    R_rad = radiative_resistance(temp_surface, temp_ambient, emissivity, area)
    R_conv = convective_resistance(h, area)
    R_parallel = 1.0 / (1.0/R_rad + 1.0/R_conv)
    # Total network = conduction series with parallel branch
    R_total = total_resistance + R_parallel
    
    return (temp_surface - temp_ambient) / R_total  # W


def calculate_total_resistance():
     
    R1 = thermal_resistance(thickness_super_wool, k_super_wool, outer_area)
    R2 = thermal_resistance(thickness_sand,       k_sand,       outer_area)
    R3 = thermal_resistance(thickness_fiberglass, k_fiberglass, outer_area)
    R4 = thermal_resistance(thickness_tank,       k_g_steel,    outer_area)
    return R1 + R2 + R3 + R4


def calculate_heat_loss(temp_surface, temp_ambient):
   
    R_tot = calculate_total_resistance()
    # print('is this working lol')
    return heat_loss(temp_surface,
                     temp_ambient,
                     R_tot,
                     emissivity,
                     outer_area,
                     h_convection)

# ""
# # one_D_model.py
# import numpy as np
 
# stefan_boltzmann = 5.67e-8  # (W/m2K4)
 
# k_super_wool = 0.09  # super wool insulation
# k_sand = 0.114       # sand
# k_fiberglass = 0.04  # fiberglass blanket
# k_g_steel = 50       # galv steel tank
 
# thickness_super_wool = 1    # change
# thickness_sand = 0.0254        # Sand thickness
# thickness_fiberglass = 0.0047625  # fiberglass thickness
# thickness_tank = 0.0047625     # steel tank thickness
 
# diameter = 1.8288  # PB diameter (m)
# hole_radius = 0.0254  # center hole rad (m)
# emissivity = 0.05  # surface emissivity  
# h_convection = 3   # convective heat transfer coefficient (W/m2-K)
 
# outer_area = np.pi * (diameter / 2) ** 2 - np.pi * hole_radius ** 2  

# def thermal_resistance(thickness, k, area):
#     return thickness / (k * area)

# def convective_resistance(h, area):
#     return 1 / (h * area)

# def radiative_resistance(temp_surface, temp_ambient, emissivity, area):
#     temp_surface_k = temp_surface + 273.15
#     temp_ambient_k = temp_ambient + 273.15
#     return 1 / (emissivity * stefan_boltzmann * area * (
#         (temp_surface_k**2 + temp_ambient_k**2) * (temp_surface_k + temp_ambient_k)
#     ))

# def heat_loss(temp_initial, temp_final, total_resistance, emissivity, area, h):
#     r_radiation = radiative_resistance(temp_initial, temp_final, emissivity, area)
#     r_convection = convective_resistance(h, area)

#     r_parallel = 1 / (1 / r_radiation + 1 / r_convection)
    
#     r_total = 1 / (1 / total_resistance + 1 / r_parallel)
    
#     return (temp_initial - temp_final) / r_total 

# def calculate_total_resistance():
#     r_super_wool = thermal_resistance(thickness_super_wool, k_super_wool, outer_area)
#     r_sand = thermal_resistance(thickness_sand, k_sand, outer_area)
#     r_fiberglass = thermal_resistance(thickness_fiberglass, k_fiberglass, outer_area)
#     r_tank = thermal_resistance(thickness_tank, k_g_steel, outer_area)
     
#     return r_super_wool + r_sand + r_fiberglass + r_tank

# def calculate_heat_loss(temp_surface, temp_ambient):
 
#     total_resistance = calculate_total_resistance()
#     return heat_loss(temp_surface, temp_ambient, total_resistance, emissivity, outer_area, h_convection)""