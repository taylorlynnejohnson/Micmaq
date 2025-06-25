import os
import pvlib
import numpy as np
import sys
sys.path.insert(0, '/projects/wg-ASGARD/bin/spt_dev/SolarPILOT/deploy/api/')
from copylot import CoPylot
import pandas as pd
import glob
from natsort import natsorted, ns
import math
import multiprocess as mp
import functools
from matplotlib import pyplot as plt

#PV fixed values
ratio_MW2acres = {'SAT':0.24,'fixed':0.35} #https://eta-publications.lbl.gov/sites/default/files/land_requirements_for_utility-scale_pv.pdf
ratio_ground_coverage = {'SAT':0.3,'fixed':0.45}
#capex_ratio_PV_USD2MW_DC = {'SAT':860000,'fixed':830000} #https://www.ny-engineers.com/blog/breaking-down-the-price-of-solar-power-systems#:~:text=The%20November%202021%20technical%20report,scale%20system%20without%20solar%20tracking
capex_ratio_PV_USD2MW_DC = {'SAT':1161000,'fixed':1150000}
#OM_ratio_PV_USD2MW_AC = {'SAT':23000, 'fixed':23000} # https://atb.nrel.gov/electricity/2023/utility-scale_pv
#OM_ratio_PV_USD2MW_DC = {'SAT':16580,'fixed':15000} #https://www.nrel.gov/docs/fy23osti/87303.pdf
OM_ratio_PV_USD2MW_DC = {'SAT':18290,'fixed':17560} #updated 2023 cost benchmark
OM_ratio_CSP_USD2MW = 40000 #is this MWt or MWe
SAT_max_tilt = 55
gamma = -0.003
temp_ref = 25

#CSP fixed values
""" Constants for Falling Particle Receiver Style 1"""
A = 0.8481
B = 0.2498
C = -1.0116
D = -7.9429e-5
E = -1.4575e-7
F = 5.5
G = 7.5
H = 5000

class PV_System: 

    def __init__(self, name, site, capacity_MW_DC, PV_array_type, PV_tilt, PV_azimuth, ratio_DC2AC, off_grid_operation, power_priority_load_MW_AC=None):
        """initiate object"""
        #set by system definition
        self.name = name
        self.site=site
        self.capacity_MW_DC = capacity_MW_DC
        #print(capacity_MW_DC)
        self.PV_array_type = PV_array_type
        self.PV_tilt = PV_tilt
        self.PV_azimuth = PV_azimuth
        self.ratio_DC2AC = ratio_DC2AC
        self.off_grid_operation = off_grid_operation
        self.power_priority_load_MW_AC = power_priority_load_MW_AC
        
        #automatic calculations
        self.capacity_MW_AC = round(self.capacity_MW_DC/self.ratio_DC2AC,3)
        self.capex_USD = capex_ratio_PV_USD2MW_DC[self.PV_array_type]*self.capacity_MW_DC
        #self.annual_OM_USD = OM_ratio_PV_USD2MW_AC[self.PV_array_type]*self.capacity_MW_AC
        self.annual_OM_USD = OM_ratio_PV_USD2MW_DC[self.PV_array_type]*self.capacity_MW_DC
        self.area_land_acres = int(round(self.capacity_MW_DC/ratio_MW2acres[self.PV_array_type],1))
        self.power_timeseries = self.calculate_power()

    def __str__(self):
        return self.name
    
    def calculate_poa(self):
        """calculate POA irradiance for the PV system"""
        solar_position = self.site.location.get_solarposition(self.site.weather_data.dataframe.index)

        if self.PV_array_type=='fixed':
            #Calculating POA irradiance for fixed tilt system
            poa = pvlib.irradiance.get_total_irradiance(
                surface_tilt=self.PV_tilt,
                surface_azimuth=self.PV_azimuth,
                dni=self.site.weather_data.dataframe.dni,
                ghi=self.site.weather_data.dataframe.ghi,
                dhi=self.site.weather_data.dataframe.dhi,
                solar_zenith=solar_position.apparent_zenith,
                solar_azimuth=solar_position.azimuth,
                model='isotropic')
            
        elif self.PV_array_type=='SAT':
            #Calculating tracking angle for SAT system
            tracking_angle = pvlib.tracking.singleaxis(
                apparent_zenith=solar_position.apparent_zenith,
                apparent_azimuth=solar_position.azimuth,
                axis_tilt=self.PV_tilt,
                axis_azimuth=self.PV_azimuth,
                max_angle=SAT_max_tilt,
                backtrack=True,
                gcr=ratio_ground_coverage[self.PV_array_type])
            #Calculating POA irradiance for SAT system
            poa = pvlib.irradiance.get_total_irradiance(
                surface_tilt=tracking_angle.surface_tilt,
                surface_azimuth=tracking_angle.surface_azimuth,
                dni=self.site.weather_data.dataframe.dni,
                ghi=self.site.weather_data.dataframe.ghi,
                dhi=self.site.weather_data.dataframe.dhi,
                solar_zenith=solar_position.apparent_zenith,
                solar_azimuth=solar_position.azimuth,
                model='isotropic')
        return poa.poa_global

    def calculate_power(self):
        """Function to calculate PV power production"""
        dfpv = self.site.weather_data.dataframe.copy()
        #dc calculation
        dfpv['irradiance_poa_global_Wpm2'] = self.calculate_poa()
        dfpv['temperature_module_C'] = pvlib.temperature.faiman(dfpv['irradiance_poa_global_Wpm2'], dfpv.temp_air, dfpv.wind_speed)
        dfpv['power_MW_DC'] = pvlib.pvsystem.pvwatts_dc(dfpv['irradiance_poa_global_Wpm2'], dfpv['temperature_module_C'], self.capacity_MW_DC, gamma, temp_ref)
        
        #ac calculation for for all systems
        dfpv['power_MW_DC'] = dfpv['power_MW_DC'].replace(np.nan,0)
        dfpv['power_MW_AC'] = dfpv['power_MW_DC'].clip(upper=self.capacity_MW_AC,inplace=False) 
        return dfpv.copy()

class CSP_System:
    def __init__(self, name, site, config, off_grid_operation, DOE_2030_targets = False):
        loc = os.getcwd()
        os.chdir('/projects/wg-ASGARD/bin/spt_dev/SolarPILOT/deploy/api')
        
        self.name = name
        self.site = site
        self.thermal_power_MW_t = config['thermal_power_MW_t']
        self.power_rating_MW_e = config['power_rating_MW_e']
        self.dni_des = config['dni_des']
        self.tower_height_m = config['tower_height_m']
        self.receiver_width_m = config['receiver_width_m']
        self.receiver_height_m = config['receiver_height_m']
        self.accept_ang_y = config['accept_ang_y']
        self.accept_ang_x = config['accept_ang_x']
        self.heliostat_width_m = config['heliostat_width_m']
        self.heliostat_height_m = config['heliostat_height_m']
        self.field_max_scaled_rad = config['field_max_scaled_rad']
        self.layout_method = config['layout_method']
        self.field_shape = config['field_shape']
        self.interaction_limit = config['interaction_limit']
        self.row_spacing_x = config['row_spacing_x']
        self.row_spacing_y = config['row_spacing_y']
        self.off_grid_operation = off_grid_operation
        self.field, self.r, self.area_land_acres, self.heliostat_area_m2 = self.create_csp_field()
        self.plot_field = self.show_field()
        

        self.par, self.pah = self.csp_power()
        self.plot_pah = self.show_PAH()
        self.plot_par = self.show_PAR()

        os.chdir(loc)

        if DOE_2030_targets == False:
            # Tower Cost
            fixed_tower = 1194000 # [$] Fixed tower cost 
            exp_tower = 0.0124 # [n/a] exponent in tower cost calculation
            C_tower = fixed_tower * math.exp(exp_tower*self.tower_height_m)
            
            # Heliostat cost
            SC_helio = 75 # [$/m2] heliostat cost per heliostat area
            SP_helio = 10 # [$/m2] heliostat site prep cost per area
            C_helio = self.heliostat_area_m2 * (SC_helio+SP_helio) #[$] heliostat field cost
            
            # Receiver Cost
            R_area = 31 #[m2] Reference receiver area
            R_C = 1172864 #[$] Referecne receiver cost
            C_rec = (R_C * ((self.receiver_height_m*self.receiver_width_m)/(R_area))**0.7)

            # Duct cost
            SC_duct = 28000
            C_duct = SC_duct * self.tower_height_m
            
            C_tot = (C_tower + C_helio + C_rec + C_duct) # [$] Total CAPEX overnight 
        else:
            fixed_tower = 3000000 # [$] Fixed tower cost 
            exp_tower = 0.0113 # [n/a] exponent in tower cost calculation
            C_tower = fixed_tower * math.exp(exp_tower*(self.tower_height_m - self.receiver_height_m/2 + (self.heliostat_area_m2**0.5)/2))
            
            # Heliostat cost
            SC_helio = 75 # [$/m2] heliostat cost per heliostat area
            SP_helio = 10 # [$/m2] heliostat site prep cost per area
            C_helio = self.heliostat_area_m2 * (SC_helio+SP_helio) #[$] heliostat field cost
            
            # Receiver Cost
            R_C = 37400 #[$/kWt] Referecne receiver cost
            C_rec = (R_C * self.receiver_height_m*self.receiver_width_m)

            # Duct Cost
            SC_duct = 28000
            C_duct = SC_duct * self.tower_height_m
            
            C_tot = (C_tower + C_helio + C_rec + C_duct) # [$] Total CAPEX overnight
    
        
        self.capex_USD = C_tot
        self.annual_OM_USD = OM_ratio_CSP_USD2MW*self.power_rating_MW_e


    
        # self.pah = self.read_csv()

    # def read_csv(self):
    #     csp_mods = glob.glob('/projects/wg-ASGARD/hybrid_data/rec_eff_arrays_data*.csv')
    #     csp_mods = natsorted(csp_mods, alg=ns.IGNORECASE)
    #     csp_mod = [m for m in csp_mods if str(self.thermal_power_MW_t) in m][0]
    #     print(csp_mod)
    #     df = pd.read_csv(csp_mod,index_col=0)
    #     return df.PAH.values


    def create_csp_field(self):
         # Field Setup
        cp = CoPylot();
        r = cp.data_create();
        
        ## Ambient Setup
        """ Setup the model & optimize field"""
        cp.data_set_string(r,'ambient.0.weather_file','/projects/wg-ASGARD/weather_data/KAFB/443653_35.09_-106.66_2013.csv')
        cp.data_set_string(r, "ambient.0.sun_type", "Gaussian sun")
        cp.data_set_string(r, "ambient.0.atm_model", "DELSOL3 clear day")
        cp.data_set_number(r, "ambient.0.sun_rad_limit", 2.73)
        
        ## Flux Sim Setup
        cp.data_set_string(r, "fluxsim.0.flux_model", "Hermite (analytical)")
        cp.data_set_string(r, "fluxsim.0.aim_method", "Simple aim points")
        
        
        ## Heliostat setup
        cp.data_set_string(r,'heliostat.0.focus_method',"At slant")
        cp.data_set_number(r,'heliostat.0.width',self.heliostat_width_m)
        cp.data_set_number(r,'heliostat.0.height',self.heliostat_height_m)
        cp.data_set_number(r,'heliostat.0.reflectivity',0.94)
        cp.data_set_number(r,'heliostat.0.track_period',1)
        #cp.data_set_number(r,'heliostat.0.cant_radius',200)
        
        
        ## Land Setup
        cp.data_set_number(r,'land.0.max_scaled_rad',self.field_max_scaled_rad)
        cp.data_set_number(r,'land.0.min_scaled_rad',0.05)
        cp.data_set_number(r,'land.0.land_mult',1.0)
        cp.data_set_number(r,'land.0.land_const',7)
        
        ## Receiver setup
        cp.data_set_string(r, "receiver.0.rec_type", "Flat plate")
        cp.data_set_number(r, "receiver.0.therm_loss_base", 0)
        cp.data_set_number(r, "receiver.0.piping_loss_coef", 0)
        cp.data_set_number(r, "receiver.0.peak_flux", 1000)
        cp.data_set_number(r, "receiver.0.rec_width", self.receiver_width_m)
        cp.data_set_number(r, "receiver.0.rec_height", self.receiver_height_m)
        cp.data_set_number(r, "receiver.0.absorptance", 1)
        cp.data_set_string(r, "receiver.0.accept_ang_type", "Rectangular")
        cp.data_set_number(r, "receiver.0.accept_ang_x", self.accept_ang_x)
        cp.data_set_number(r, "receiver.0.accept_ang_y", self.accept_ang_y)
        cp.data_set_number(r, "receiver.0.rec_azimuth", 0)
        cp.data_set_number(r, "receiver.0.rec_offset_x", 0)
        
        ## Solar Field Setup
        cp.data_set_string(r, "solarfield.0.des_sim_detail", "Representative profiles")
        cp.data_set_string(r, "solarfield.0.hsort_method", "Power to receiver")
        cp.data_set_string(r, "solarfield.0.layout_method", self.layout_method)
        cp.data_set_string(r, "solarfield.0.sun_loc_des", "Summer solstice")
        cp.data_set_string(r, "solarfield.0.xy_field_shape", self.field_shape)
        cp.data_set_number(r, "solarfield.0.interaction_limit", self.interaction_limit)
        cp.data_set_number(r, "solarfield.0.des_sim_nhours", 5)
        cp.data_set_number(r,'solarfield.0.row_spacing_x',self.row_spacing_x)
        cp.data_set_number(r,'solarfield.0.row_spacing_y',self.row_spacing_y)
        cp.data_set_number(r,'solarfield.0.q_des',self.thermal_power_MW_t)
        cp.data_set_number(r,'solarfield.0.dni_des',self.dni_des)
        cp.data_set_number(r,'solarfield.0.tht',self.tower_height_m)
        cp.generate_layout(r)
        
        field = cp.get_layout_info(r)

        area_land_acres = cp.data_get_number(r,'land.0.land_area')
        area_heliostat_total_m2 = cp.data_get_number(r,'solarfield.0.sf_area')
        
        return field, r, area_land_acres, area_heliostat_total_m2

    def show_field(self):
        plt.scatter(self.field['x_location'], self.field['y_location'], s=1.5)
        plt.tight_layout()
        plt.show()

    def show_PAH(self):
        x = range(len(self.pah))
        plt.plot(x, self.pah)
        plt.tight_layout()
        plt.show()
    def show_PAR(self):
        x = range(len(self.par))
        plt.plot(x, self.par)
        plt.tight_layout()
        plt.show()

    def cp_eval(self,list):
        cp = CoPylot()
        dni=list[0]
        year=list[1]
        month=list[2]
        day=list[3]
        hour=list[4]
        wind_speed = list[5]
        wind_direction = list[6]
        cp.data_set_number(self.r,'fluxsim.0.flux_dni',dni)
        cp.data_set_number(self.r,'fluxsim.0.flux_year',year)
        cp.data_set_number(self.r,'fluxsim.0.flux_month',month)
        cp.data_set_number(self.r,'fluxsim.0.flux_day',day)
        cp.data_set_number(self.r,'fluxsim.0.flux_hour',hour)
        cp.simulate(self.r)  
        
        dict_ = cp.summary_results(self.r)
        
        try:
            val=np.array([dict_.get('Power absorbed by the receiver')])
            num_val = float(val)
            val=num_val
        except:
            val=0
        
        par = val/1000.

        rec_area = self.receiver_width_m*self.receiver_height_m
        
        eta = CSP_System.receiver_efficiency(par, rec_area, dni, wind_speed, wind_direction)
        pah = eta*par
        
        return par,pah, eta

    def receiver_efficiency(par, rec_area, dni, wind_speed, wind_direction):
        #single receiver
        omega = ((180-abs(180-wind_direction))**F)*math.exp(-((180-(abs(180-wind_direction)))/G))/H
        q = math.exp(-par/rec_area)
        eta = A + B*q + C*q**2 + D*q*wind_speed*omega + E*wind_speed**2*omega
        if eta < 0 :
            eta = 0
        
        return eta
        
    def csp_power(self):
        with mp.Pool(processes=mp.cpu_count()) as pool:
            PAR,PAH,ETA=zip(*pool.map(self.cp_eval, self.site.weather_data.dataframe[['dni','Year','Month','Day','Hour','wind_speed','wind_direction']].values.tolist()))

        return PAR,PAH
