
import numpy as np

def material_library_fn(material, TempC):

    norm = np.ones_like(TempC)

    if material[:3].lower() == 'air':
        c = 1.0e3 + (150 / 800) * TempC  # J/kg-C
        TgKelvin = TempC + 273.15
        k = 1.5207e-11 * TgKelvin**3 - 4.8574e-8 * TgKelvin**2 + 1.0184e-4 * TgKelvin - 3.9333e-4  # W/mk
        mu0 = 17.15e-6
        mu = mu0 * (TgKelvin / 273)**1.5 * ((273 + 113) / (TgKelvin + 113))
        rho_g0 = 1.0 * norm  # at Albuquerque elevation
        density = rho_g0 * (20 + 273.15) / (TempC + 273.15)  # kg/m^3
    elif material[:3].lower() == 'wat': # water 
 
        c = 4182 * norm  # J/kg-C
        k = 0.598 * norm  # W/m-C
        mu = 0.001002 * norm  # Pa-s
        density = 1000 * norm  # kg/m^3
    elif material[:3].lower() == 'roc': # rock 
 
        c = 840 * norm  # J/kg-C
        k = 1.5 * norm  # W/m-C
        mu = None  # Pa-s
        density = 2500 * norm  # kg/m^3
    elif material[:3].lower() == 'lim': # limestone

        c = (920 + 0.835 * TempC) * norm  # J/kg-C  
        k = 1.5 * norm  # W/m-C
        mu = None  # Pa-s
        density = 2760 * norm  # kg/m^3
    else:
        c = None  # J/kg-C
        k = None  # W/m-C
        mu = None  # Pa-s
        density = 2500  # kg/m^3 of rock

    return c, k, mu, density