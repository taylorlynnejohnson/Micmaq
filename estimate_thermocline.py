# estimate_thermocline.py
import numpy as np
from material_library_fn import material_library_fn

def estimate_thermocline(M_dot, Area, eps, alpha, Dp, Tg, HTF_type):

    # get material properties for the heat transfer fluid (HTF) and solid
    c_g, k_g, mu_g, rho_g = material_library_fn(HTF_type, Tg)
    c_s, k_s, mu_s, rho_s = material_library_fn('rock', Tg)

    # calc superficial velocity and local velocity
    print(f'Mdot: {M_dot}')
    v_s = abs(M_dot) / (Area * rho_g) 
    print(f'v_s: {v_s}')
    v_local = v_s / eps 

    # Re and Prandtl number for packed bed 
    Re = rho_g * abs(v_s) * Dp / mu_g 
    Pr = c_g * mu_g / k_g

    # Nusselt number using Gunn correlation
    Nu_Gunn = ((7 - (10 * eps) + 5 * (eps ** 2)) * (1 + 0.7 * (Re ** 0.2) * (Pr ** 0.33)) +
               (1.33 - (2.4 * eps) + 1.2 * (eps ** 2)) * (Re ** 0.7) * (Pr ** 0.33))
    h_sg = (k_g / Dp) * Nu_Gunn

    # overall heat transfer coefficient between solid and gas
    H_sg = 6 * h_sg * (1 - eps) / (alpha * Dp)

    # estimate thermocline length
    estThermoclineLength = 50 * v_s * c_g * rho_g / H_sg

    # calculate time constant and thermocline velocity
    k1 = (6 * (1 - eps) * h_sg / (alpha * Dp)) * (1 - eps) / (c_s * rho_s)
    t_star = (10 / k1) + estThermoclineLength / v_s
    estThermoclineVelocity = estThermoclineLength / t_star

    return estThermoclineLength, v_s, estThermoclineVelocity
