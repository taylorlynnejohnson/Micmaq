# one_D_model.py

import numpy as np

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

def outer_area(diameter):
    return np.pi * (diameter/2)**2 - np.pi * (inner_diameter/2)**2

# Surface parameters
emissivity      = 0.05   # Surface emissivity
h_convection    = 5    # Convective coefficient top annulus (W/(m²·K))


def thermal_resistance(thickness, k, area):
    # Series conduction resistance: R = L/(k·A)
    return thickness / (k * area)


def convective_resistance(h, area):
    # Convection resistance: R = 1/(h·A)
    return 1.0 / (h * area)


def radiative_resistance(T_surf, T_amb, emissivity, area):
    # linearized radiative resistance at surface T_surf to ambient T_amb
    Ts_k = T_surf
    Ta_k = T_amb
    # From dQ/dT = εσA(T⁴_surf − T⁴_amb)/(T_surf − T_amb)
    denom = emissivity * stefan_boltzmann * area * (
        (Ts_k**2 + Ta_k**2) * (Ts_k + Ta_k)
    )
    return 1 / denom


def calculate_total_resistance(thickness_super_wool, diameter):

    A = outer_area(diameter)
    R1 = thermal_resistance(thickness_super_wool, k_super_wool, A)
    R2 = thermal_resistance(thickness_sand,       k_sand,       A)
    R3 = thermal_resistance(thickness_fiberglass, k_fiberglass, A)
    R4 = thermal_resistance(thickness_tank,       k_g_steel,    A)
    return R1 + R2 + R3 + R4


def calculate_axial_annular_roof_loss(temp_surface, temp_ambient,
                              thickness_super_wool, outer_diameter):
   
    A_roof = np.pi*(outer_diameter/2)**2 - np.pi*(inner_diameter/2)**2

    # conduction through all layers
    R_cond = ( thickness_super_wool  / (k_super_wool  * A_roof)
             + thickness_sand        / (k_sand        * A_roof)
             + thickness_fiberglass  / (k_fiberglass  * A_roof)
             + thickness_tank        / (k_g_steel     * A_roof)
             )

    # film + radiation at the outside of the roof
    R_rad  = radiative_resistance(temp_surface, temp_ambient, emissivity, A_roof)
    R_conv = convective_resistance(h_convection, A_roof)
    R_par  = 1.0 / (1.0/R_rad + 1.0/R_conv)

    R_tot  = R_cond + R_par
    return (temp_surface - temp_ambient) / R_tot  # W

def calculate_bottom_annular_loss(
    temp_surface,
    temp_ambient,
    outer_diameter
):
    # Annular area at bottom
    A_bot = np.pi*(outer_diameter/2)**2 - np.pi*(inner_diameter/2)**2
    R_cond = thickness_tank / (k_g_steel * A_bot)
    temp_ground = temp_ambient + 25
    k_concrete = 1.2  # W/mK
    thickness_concrete = 0.5 # m
    R_ground = thickness_concrete / (k_concrete * A_bot)
    R_total = R_cond + R_ground
    return (temp_surface - temp_ground) / R_total