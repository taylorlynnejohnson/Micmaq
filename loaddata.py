import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_greenhouse_load(weather_files, 
                              T_target=30, 
                              T_max=40,
                              footprint=200, 
                              height=4,
                              wall_thickness=0.0127, 
                              roof_thickness=0.0127, 
                              rock_thickness=0.0227, 
                              transmission_roof=0.8,
                              U_factor=0.62,          
                              initial_rock_temp=23, 
                              timestep_hours=1):

    Cp_air = 1005    
    rho_air = 1.225  
    Cp_rock = 800    
    rho_rock = 1600  
    k_wall = 0.0202  
    k_roof = 0.21    
  
    volume = footprint * height          
    wall_area = 4 * np.sqrt(footprint) * height  
    roof_area = footprint                
    ACH = 0  
    V_exchanged = volume * ACH  

    T_rock = initial_rock_temp
    T_greenhouse = T_target  

    log_data = []

    for weather_file in weather_files:
        df_weather = pd.read_csv(weather_file, skiprows=[0, 1], header=0)  
        temperature_array = df_weather['Temperature'].values
        wind_speed_array = df_weather['Wind Speed'].values
        solar_irradiance_array = df_weather['GHI'].values
        
        for i in range(len(temperature_array)):
            T_out = temperature_array[i]
            wind_speed = wind_speed_array[i]
            solar_irradiance = solar_irradiance_array[i]
            
            h_outside = 10.45 - wind_speed + 10 * np.sqrt(wind_speed)
            h_outside = max(h_outside, 0.0)
            
            R_wall = (wall_thickness / k_wall) + (1.0 / h_outside)
            R_roof = (roof_thickness / k_roof) + (1.0 / h_outside)
            
            Q_wall = (T_out - T_greenhouse) * (wall_area / R_wall)  
            Q_roof = (T_out - T_greenhouse) * (roof_area / R_roof)  
            Q_air = Q_wall + Q_roof                         

            Q_solar = solar_irradiance * roof_area * transmission_roof  
               
            Q_leakage = - (V_exchanged * rho_air * Cp_air * (T_greenhouse - T_out) * 0.000278)  
       
            Q_net = Q_air + Q_solar + Q_leakage  
            
            delta_T = (Q_net / (rho_air * volume * Cp_air)) * timestep_hours * 3600  
            T_greenhouse += delta_T

            if T_greenhouse > T_max:
                T_greenhouse = T_max

            if T_greenhouse < T_target:
                net_load = (T_target - T_greenhouse) * rho_air * volume * Cp_air / (3600 * timestep_hours)  
                T_greenhouse = T_target  
            else:
                net_load = 0
            
            log_data.append({
                'Timestep': i,
                'T_out': T_out,
                'WindSpeed': wind_speed,
                'SolarIrradiance': solar_irradiance,
                'Q_air': Q_air,
                'Q_solar': Q_solar,
                'Q_leakage': Q_leakage,
                'Q_net': Q_net,
                'T_greenhouse': T_greenhouse,
                'T_rock': T_rock,  
                'NetLoad': net_load,
                'T_amb': temperature_array
            })
     
    load_df = pd.DataFrame(log_data)
     
    load_df['NetLoad_MW'] = load_df['NetLoad'] / 1e6  

    total_Q_net_Wh = load_df['Q_net'].sum() / 1e6  
    print(f"Total Q_net for the year: {total_Q_net_Wh:.2f} MWh")
  
    num_hours = len(load_df)
    date_range = pd.date_range(start='2023-01-01', periods=num_hours, freq=f'{timestep_hours}H')
    load_df.index = date_range
    load_df.index.name = 'DateTime'
    
    return load_df