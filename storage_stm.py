# -*- coding: utf-8 -*-
"""
storage.py: Storage system (BES,TES) definition
"""
import math
import numpy as np
from matplotlib import pyplot as plt


#BES constants
#capex_BES_BOS_USDpkW=1785.08-1080.88 #https://atb.nrel.gov/electricity/2023/utility-scale_battery_storage
capex_BES_INV_USDpkW=68 #updated with 2023 cost of bi-inverters
#capex_BES_capacity_USDpkWh=270
capex_BES_capacity_USDpkWh=87.4 #updated with 2023 cost balance minus bi-inverters
capex_BES_BOS_USDpkWh=226 #balance of systems (costs less inverter and cells)
ratio_BES_acres2MWh = 4.5/730 #https://www.teslarati.com/tesla-megapack-pge-monterey-ca-moss-landing/
DOD_f = 0.6 # Augmentation batteries depth of discharge
N_80 = 6.5*365-1 # Cycles to 80% rating at 80% DOD ~3% degradation
N_60 = 22*365-1 # Cycles to 60% rating at 60% DOD ~1.3% degradation

#TES constants
SC_skip = 0#18000 # [$/MWht] Skip hoist cost 
SC_TES = 10000 # [$/MWht] TES tank + insulation + media cost
SC_pb = 0#1683000 # [$/MWe] Powerblock cost
SC_hx = 0#300000 # [$/MWe] Particle-to-SCO2 PHX cost per MWe
SC_insulation = 4  # [$/sqft per inch thickness] for superwool insulation layer cost
SC_media = 0.5  # [$/kg] for pea gravel packed bed material

# TES Heater Replacement Constants
C_controller = 0#138 # [$/kW] heater contoller cost
C_shell = 0#75 # [$/kW] heater shell cost
C_rod = 20 #17 # [$/kW] heater rod cost

# Blower cost
C_blower = 500  # research later

# O&M Defined Input
C_OM_TES = 5000 # [$/MW-year] TES O&M cost
C_OM_Battery = 0#9870 # [$/MW-year] Battery O&M cost


class BES_System:
    def __init__(self, name, site, capacity_MWh_e, power_rating_MW_e, percent_discharge_depth, charge_efficiency, discharge_efficiency, systems_charging, start_full, off_grid_operation, max_total_power_MW):
        self.name = name
        self.site = site
        self.capacity_MWh_e = capacity_MWh_e
        self.power_rating_MW_e = power_rating_MW_e
        self.percent_discharge_depth = percent_discharge_depth
        self.charge_efficiency = charge_efficiency
        self.discharge_efficiency = discharge_efficiency
        self.systems_charging = systems_charging
        self.start_full = start_full
        self.off_grid_operation = off_grid_operation
        self.max_total_power_MW = max_total_power_MW
        self.area_land_acres = self.capacity_MWh_e * ratio_BES_acres2MWh
        self.capex_BES_capacity_USDpkWh = capex_BES_capacity_USDpkWh
        
        if self.power_rating_MW_e is not None:
            self.capex_USD = self.capex_calc(self.power_rating_MW_e)
        elif self.power_rating_MW_e is None:
            self.capex_USD = None

        self.annual_OM_USD = self.capacity_MWh_e * C_OM_Battery

    def capex_calc(self, power_rating):
        return int(self.capacity_MWh_e * 1000 * (capex_BES_capacity_USDpkWh + capex_BES_BOS_USDpkWh) + power_rating * 1000 * capex_BES_INV_USDpkW)

class TES_System_STM:
    def __init__(self, name, site, capacity_MWh_t, power_rating_MW_e, power_minimum_MW_e, efficiency_rating_t2e, percent_discharge_depth, percent_heat_loss_daily, charge_rate_CSP_MW_t, charge_efficiency_t2TES, charge_rate_resistive_MW_e, charge_efficiency_e2TES, systems_charging, start_full, off_grid_operation=True, DOE_2030_targets=False):

        self.name = name
        self.site = site
        self.capacity_MWh_t = capacity_MWh_t
        self.power_rating_MW_e = power_rating_MW_e
        self.power_minimum_MW_e = power_minimum_MW_e
        self.efficiency_rating_t2e = efficiency_rating_t2e
        self.percent_discharge_depth = abs(percent_discharge_depth)
        self.percent_heat_loss_daily = abs(percent_heat_loss_daily)
        self.charge_rate_CSP_MW_t = charge_rate_CSP_MW_t
        self.charge_efficiency_t2TES = charge_efficiency_t2TES
        self.charge_rate_resistive_MW_e = charge_rate_resistive_MW_e
        self.charge_efficiency_e2TES = charge_efficiency_e2TES
        self.systems_charging = systems_charging
        self.start_full = start_full
        self.DOE_2030_targets = DOE_2030_targets
        self.off_grid_operation = off_grid_operation

        self.area_land_acres = 0
        self.capacity_MWh_e = self.capacity_MWh_t * self.efficiency_rating_t2e

        if self.charge_rate_resistive_MW_e is not None:
            self.capex_USD = self.calculate_capex(self.charge_rate_resistive_MW_e)
            self.cpx_htr_aug_annual = self.calculate_htr_aug(self.charge_rate_resistive_MW_e)
        elif self.charge_rate_resistive_MW_e is None:
            self.capex_USD = None
            self.cpx_htr_aug_annual = None
        self.annual_OM_USD = self.power_rating_MW_e * C_OM_TES

    def calculate_htr_aug(self, charge_rate_resistive):
        prh = charge_rate_resistive
        years = 30

        # Determine Year N heater costs
        cpx_heater_aug_annual = np.zeros(years)
        for year in range(years):
            cpx_heater_aug_annual[year] = prh * 1000 * C_rod if year in [10, 20, 30, 40] else 0
        return list(cpx_heater_aug_annual)

    def calculate_capex(self, charge_rate_resistive):

        if self.charge_rate_CSP_MW_t >= self.charge_rate_resistive_MW_e:
            C_skip = self.charge_rate_CSP_MW_t * SC_skip
        else:
            C_skip = self.charge_rate_resistive_MW_e * self.charge_efficiency_e2TES * SC_skip

        C_heater = (C_rod + C_shell + C_controller) * charge_rate_resistive * 1000
        
        # Detailed TES cost calculations
        insulation_area = math.pi * (6 * 0.3048) ** 2
        C_insulation = SC_insulation * insulation_area * 2  # Assuming 2 inches superwool
        
        media_mass = insulation_area * 0.4572 * 2500  # density of basalt * packed bed volume
        C_media = media_mass * SC_media

        blower_cost = C_blower  # research later

        if self.DOE_2030_targets == False:
            C_TES = SC_TES * self.capacity_MWh_t + C_skip + C_insulation + C_media + blower_cost
            C_hx = SC_hx * self.power_rating_MW_e
            #C_pb = SC_pb * self.power_rating_MW_e
            C_pb = 0 
        else:
            C_TES = 32000 * self.capacity_MWh_t + C_skip + C_insulation + C_media + blower_cost
            C_hx = SC_hx * self.power_rating_MW_e
            #C_pb = 1312000 * self.power_rating_MW_e
            C_pb = 0

        return C_TES + C_hx + C_heater + C_pb