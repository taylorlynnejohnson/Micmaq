import numpy as np

# Physical constants
stefan_boltzmann = 5.67e-8  # Stefan-Boltzmann constant (W/m^2-K^4)

# Thermal conductivities (W/m-K)
k_super_wool = 0.09  # Super wool insulation
k_sand = 0.114       # Sand
k_fiberglass = 0.04  # Fiberglass blanket
k_g_steel = 50       # Galvanized steel tank

# Layer thicknesses (meters)
thickness_super_wool = 1.016   # 40 inches in meters
thickness_sand = 0.0254        # Sand thickness
thickness_fiberglass = 0.0047625  # Fiberglass thickness
thickness_tank = 0.0047625     # Steel tank thickness

# Default geometry and emissivity
diameter = 1.8288  # Diameter of packed bed (m)
hole_radius = 0.0254  # Central hole radius (m)
emissivity = 0.05  # Surface emissivity (dimensionless)
h_convection = 3   # Convective heat transfer coefficient (W/m^2-K)

# Derived properties
outer_area = np.pi * (diameter / 2) ** 2 - np.pi * hole_radius ** 2  # Area (m^2)

def thermal_resistance(thickness, k, area):
    """Calculate the thermal resistance for a layer."""
    return thickness / (k * area)

def convective_resistance(h, area):
    """Calculate the convective resistance."""
    return 1 / (h * area)

def radiative_resistance(temp_surface, temp_ambient, emissivity, area):
    """Calculate the radiative resistance."""
    temp_surface_k = temp_surface + 273.15
    temp_ambient_k = temp_ambient + 273.15
    return 1 / (emissivity * stefan_boltzmann * area * (
        (temp_surface_k**2 + temp_ambient_k**2) * (temp_surface_k + temp_ambient_k)
    ))

def heat_loss(temp_initial, temp_final, total_resistance, emissivity, area, h):
    """
    Calculate the total heat loss rate (W) through conduction, convection, and radiation.
    """
    r_radiation = radiative_resistance(temp_initial, temp_final, emissivity, area)
    r_convection = convective_resistance(h, area)
    
    # Combine radiative and convective resistances in parallel
    r_parallel = 1 / (1 / r_radiation + 1 / r_convection)
    
    # Combine with conductive resistance in series
    r_total = 1 / (1 / total_resistance + 1 / r_parallel)
    
    return (temp_initial - temp_final) / r_total  # Heat loss rate (W)

def calculate_total_resistance():
    """
    Compute the total thermal resistance of the packed bed system
    based on its layer properties.
    """
    # Conductive resistances
    r_super_wool = thermal_resistance(thickness_super_wool, k_super_wool, outer_area)
    r_sand = thermal_resistance(thickness_sand, k_sand, outer_area)
    r_fiberglass = thermal_resistance(thickness_fiberglass, k_fiberglass, outer_area)
    r_tank = thermal_resistance(thickness_tank, k_g_steel, outer_area)
    
    # Total conductive resistance (series sum)
    return r_super_wool + r_sand + r_fiberglass + r_tank

def calculate_heat_loss(temp_surface, temp_ambient):
    """
    Calculate heat loss rate (W/m^2) for given surface and ambient temperatures.
    """
    total_resistance = calculate_total_resistance()
    return heat_loss(temp_surface, temp_ambient, total_resistance, emissivity, outer_area, h_convection)