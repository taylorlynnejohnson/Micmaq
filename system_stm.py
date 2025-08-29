# to-do:
# check curtailment calc, rerun PV sweep
# load/generation/TES profiles
# 50 debt fraction 13 COE 30 ITC loan percentage of 0 base case for PV curtailment 
# CAD engineering drawings of the assembly, subassemblies, and main parts 
# final BOM (parts, estimate of price) in excel 
# final report 


from storage_stm import *
from generation_stm import *
from itertools import cycle, islice, chain
from random import randrange
import datetime as DT
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy_financial as npf
import multiprocess as mp
from random import randint

# system_stm.py
# vary financing assumptions, make heat map plot    
# LCOE (x), multiple y axes  

#L = 30 # [yrs] Operational Life
#n = 30 # [yrs] Analysis period
tax = 0.257 # [frac] State and federal tax rate  # 
inflation = 0.028 # [frac] Inflation rate
ITC = 0.5 # [frac] Ineternal tax credit ########### vary [0, 0.3, 0.5]
insurance = 0.004   #[frac] Insurance rate
property_tax = 0.0084 # [frac] Property Tax Rate
I = 0.08 # [frac] Nominal Interest rate 
COE = 0.13 # [frac] Cost of equity      ############ vary [.065, 0.13, 0.26]
DF=0.5 #  [frac] Debt fraction          ############ vary [0, 0.5, 1]
MACRS_yrs = 7 # [yrs] MACRS Depreciation Period 
PVD = 0.73281282777303 # [frac] Present Value of Depreciation based on MACRS_yrs
esc = 0.02 # [frac] Escalation rate 
grant_frac = 0                           ############ vary [0, .1, .25]

def set_finance_params(**kwargs):
    allowed = {
        'tax','inflation','insurance','property_tax',
        'I','COE','DF','MACRS_yrs','PVD','esc','grant_frac'
    }
    for k, v in kwargs.items():
        if k in allowed:
            globals()[k] = v
        else:
            raise KeyError(f"Unknown finance param: {k}")


class System: #need to set dict objects to exclude timeseries for metrics
  
    def __init__(self, load_MW, systems_load_order, analysis_period = 30, ITC = 0.5,fuel_type='LPG'):
        self.load_MW = load_MW
        self.systems = systems_load_order
        self.system_names = [sys.name for sys in self.systems]
        self.sites = list(set([sys.site for sys in self.systems]))
        self.analysis_period = analysis_period
        self.ITC = ITC
        self.fuel_type = fuel_type

        #separating systems by type
        self.pv_systems = [sys for sys in self.systems if isinstance(sys,PV_System)]
        self.csp_systems = [sys for sys in self.systems if isinstance(sys,CSP_System)]
        self.bes_systems = [sys for sys in self.systems if isinstance(sys,BES_System)]
        self.tes_systems = [sys for sys in self.systems if isinstance(sys,TES_System_STM)]

        #defining off-grid systems
        self.off_grid_systems = [sys for sys in self.systems if sys.off_grid_operation==True]

        self.pv_systems_off_grid = [sys for sys in self.off_grid_systems if isinstance(sys,PV_System)]
        self.csp_systems_off_grid = [sys for sys in self.off_grid_systems if isinstance(sys,CSP_System)]
        self.bes_systems_off_grid = [sys for sys in self.off_grid_systems if isinstance(sys,BES_System)]
        self.tes_systems_off_grid = [sys for sys in self.off_grid_systems if isinstance(sys,TES_System_STM)]
        self.e_sale = 0.055

        #saving pv power to system variable
        for sys in self.pv_systems:
            setattr(self, sys.name, sys.power_timeseries.copy().reset_index(inplace=True))

        #saving csp heat to system variables
        for csp in self.csp_systems:
            setattr(self, csp.name, csp.pah)
        
        self.timeseries = self.operation()

        #determination of TES and BES power ratings and capex if set to None

        for tes in self.tes_systems:
            if tes.charge_rate_resistive_MW_e is None:
                pv_cols=[]
                for sys in self.pv_systems:#tes.systems_charging:
                    if sys.name in tes.systems_charging:# sys.__class__.__name__==self.pv_systems[0].__class__.__name__:
                        pv_cols.append(sys.name)
                pv_cols = [pv_name+'_to_'+tes.name+'_MWh_e' for pv_name in pv_cols]
                    
                tmax = self.timeseries[pv_cols].sum(axis=1).max()
                tes.charge_rate_resistive_MW_e = tmax
                tes.capex_USD = tes.calculate_capex(tmax)
                tes.cpx_htr_aug_annual=tes.calculate_htr_aug(tmax)

        for bes in self.bes_systems:
            if bes.power_rating_MW_e is None:
                pv_cols=[]
                for sys in self.pv_systems:#bes.systems_charging:
                    if sys.name in bes.systems_charging:
                        pv_cols.append(sys.name)
                pv_cols = [pv_name+'_to_'+bes.name+'_MWh_e' for pv_name in pv_cols]
                cmax = self.timeseries[pv_cols].sum(axis=1).max()
                bmax = max(self.timeseries[bes.name+'_to_load_MWh_e'].max(),cmax)
                bes.power_rating_MW_e = bmax
                bes.capex_USD = bes.capex_calc(bmax)

        self.metrics = self.system_metrics()

    def operational_systems(self, operation):
        if operation=='normal':
            pv_systems = self.pv_systems
            csp_systems = self.csp_systems
            bes_systems = self.bes_systems
            tes_systems = self.tes_systems

        elif operation=='off_grid':
            pv_systems = self.pv_systems_off_grid
            csp_systems = self.csp_systems_off_grid
            bes_systems = self.bes_systems_off_grid
            tes_systems = self.tes_systems_off_grid
        
        return pv_systems, csp_systems, bes_systems, tes_systems
        

    def initialize_df(self, operation='normal', critical_load_MW = None, df_res = None, res_initial = None):

        #determining which systems are being used (dynamic for off-grid scenario)
        pv_systems, csp_systems, bes_systems, tes_systems = self.operational_systems(operation)

        #load handling
        if operation == 'normal':
            df = self.load_MW.copy().reset_index(inplace=False)
            df.load_MW = [abs(l) for l in df.load_MW.values] #load is positive
        elif operation == 'off_grid':
            df = df_res.load_MW.copy().reset_index(inplace=False)
            df.load_MW = critical_load_MW
            
        numrec = len(df)

        df['fuel_to_load_MWh_th'] = 0
        #df['RES_to_load_MWh_th'] = 0

        #initialize pv system power values and flows
        for pv in pv_systems:
            if operation == 'normal':
                ts = pv.power_timeseries.power_MW_AC.values
            elif operation == 'off_grid':
                ts = df_res[pv.name+'_power_MW_AC'].values
            #matching pv power series length to load length
            if len(ts)>=numrec:
                ts = ts[:numrec] #truncate if power timeseries is longer than load timeseries
            elif len(ts)<numrec:
                ts = list(islice(cycle(ts), numrec)) #repeat power timeseries if shorter than load timeseries
            df[pv.name+'_power_MW_AC'] = ts

            df[pv.name+'_to_load_MWh_e'] = 0
            df[pv.name+'_to_grid_MWh_e'] = 0
            df[pv.name+'_curtailed_MWh_e'] = 0

            #initialize charging of bes by pv
            for bes in bes_systems:
                if pv.name in bes.systems_charging: #check if charging allowed by pv system
                    df[pv.name+'_to_'+bes.name+'_MWh_e'] = 0

            #initialize charging of tes by pv
            for tes in tes_systems:
                df[pv.name+'_to_'+tes.name+'_MWh_e'] = 0

        #initializing bes variables based on system parameters
        for bes in bes_systems:
            if operation == 'normal':
                #initializing BES full or empty
                if bes.start_full:
                    start_bes = bes.capacity_MWh_e
                elif not bes.start_full:
                    start_bes = 0
            elif operation == 'off_grid':
                if res_initial == 'actual':
                    start_bes = df_res[bes.name+'_MWh_DC'].values[0]
                elif res_initial == 'full':
                    start_bes = bes.capacity_MWh_e
            df[bes.name+'_MWh_DC'] = start_bes
            df.loc[1:,bes.name+'_MWh_DC'] = 0

            df[bes.name+'_to_load_MWh_e'] = 0

        #initializing csp variables
        for csp in self.csp_systems:
            if operation == 'normal':
                ts = csp.pah
            elif operation == 'off_grid':
                ts = df_res[csp.name+'_heat_MW_t'].values
            #matching csp ts length to load length
            if len(ts)>=numrec:
                ts = ts[:numrec] #truncate if power timeseries is longer than load timeseries
            elif len(ts)<numrec:
                ts = list(islice(cycle(ts), numrec)) #repeat power timeseries if shorter than load timeseries
            df[csp.name+'_heat_MW_t'] = ts
            df[csp.name+'_heat_unused_MWh_t'] = ts

            #initializing power flow csp -> tes
            for tes in self.tes_systems:
                if csp.name in tes.systems_charging:
                    df[csp.name+'_to_'+tes.name+'_MWh_t'] = 0

        #initializing tes variables
        for tes in tes_systems:
            if operation == 'normal':
                #initializing TES full or empty
                if tes.start_full: #change to variable starting capacity?
                    start_tes = tes.capacity_MWh_t
                elif not tes.start_full:
                    start_tes = 0
            elif operation == 'off_grid':
                if res_initial == 'actual':
                    start_tes = df_res[tes.name+'_MWh_t'].values[0]
                elif res_initial == 'full':
                    start_tes = tes.capacity_MWh_e
            df[tes.name+'_MWh_t'] = start_tes
            df.loc[1:,tes.name+'_MWh_t'] = 0
            df[tes.name+'_thermal_loss_MWh_t'] = 0

            df[tes.name+'_to_load_MWh_t'] = 0
            df[tes.name+'_to_load_MWh_e'] = 0

        for site in self.sites:
            df[site.name+'_POI_MW'] = 0
        
        return pv_systems, csp_systems, bes_systems, tes_systems, df

    def operation(self, operation='normal', critical_load_MW = None, df_res = None, res_initial = None):
    
        # fuel properties
        fuels = {
            'LPG': {'heating_value': 46.1, 'efficiency': 0.83}, 
        }
        fuel_properties = fuels[self.fuel_type]

        if operation == 'normal':
            pv_systems, csp_systems, bes_systems, tes_systems, df = self.initialize_df(operation)
        elif operation == 'off_grid':
            pv_systems, csp_systems, bes_systems, tes_systems, df = self.initialize_df(operation, critical_load_MW = critical_load_MW, df_res = df_res, res_initial = res_initial)
            to_end = 0
            endseries = df.index[-1]
    
        for i in range(len(df)-1):

            #local variable to track load satisfaction
            unmet_load_MW = df.loc[i+1,'load_MW']
            
            #tes heat losses first
            for tes in tes_systems:
                #local variable to track charge from previous
                tes_MWh_t = df.loc[i, tes.name+'_MWh_t']
                #percent daily loss converted to fraction for hourly loss
                tloss = tes_MWh_t*tes.percent_heat_loss_daily/2400
                #update local variable with thermal loss, ensure nonnegative
                tes_MWh_t = max(tes_MWh_t-tloss,0)
                df.loc[i+1,tes.name+'_MWh_t'] = tes_MWh_t
                #thermal losses for this timestep after csp charging contribution
                df.loc[i+1,tes.name+'_thermal_loss_MWh_t'] = tloss

                tes_charge_csp_remaining = tes.charge_rate_CSP_MW_t
                tes_cols = []

                for csp in csp_systems:
                    #ensures each csp system is checked for charging each tes
                    if csp.name not in tes.systems_charging:
                        continue
                    #tracking in loop, updated end of loop
                    csp_remaining = df.loc[i+1,csp.name+'_heat_unused_MWh_t']
                    #amount of csp heat transferred. limited by availability, charge rate, and tes headroom

                    #calculating tes headroom
                    headroom_tes = tes.capacity_MWh_t-tes_MWh_t
                    
                    csp2tes = max(0,min([val for val in [tes_charge_csp_remaining, csp_remaining, headroom_tes/tes.charge_efficiency_t2TES] if val is not None])) #check charge rate efficiency

                    #csp heal left over
                    csp_remaining -= csp2tes
                    
                    #local tes variable updated with heat transfer scaled by charging efficiency
                    tes_MWh_t += csp2tes*tes.charge_efficiency_t2TES

                    #tracking charge_csp_remaining
                    if tes_charge_csp_remaining is not None:
                        tes_charge_csp_remaining -= csp2tes
                    
                    #csp heat at this timestep
                    df.loc[i+1,csp.name+'_to_'+tes.name+'_MWh_t'] = csp2tes
                    tes_cols.append(csp.name+'_to_'+tes.name+'_MWh_t')
                    
                    #tracking remaining csp heat after charging TES
                    df.loc[i+1,csp.name+'_heat_unused_MWh_t'] = csp_remaining
                    
                #updated tes charge at this timestep after charging by all csp
                df.loc[i+1,tes.name+'_MWh_t'] = tes_MWh_t
                if len(tes_cols)>0:
                    df.loc[i+1,csp.__class__.__name__+'_to_'+tes.name+'_MWh_t'] = sum([df.loc[i+1,c] for c in tes_cols])

                    
            #next allocating pv power
            for pv in pv_systems:
                #local variable to track pv power remaining
                pv_remaining = df.loc[i+1, pv.name+'_power_MW_AC']
                #local variable to track site POI
                if pv.site.POI_limit != None:
                    poi_remaining = pv.site.POI_limit - df.loc[i+1,pv.site.name+'_POI_MW']
                else:
                    poi_remaining = np.inf
                pv_poi=0
                #pv to load, limited by unmet load, power available, and power limits (if directing to BES systems)
                pv2l = min([val for val in [unmet_load_MW, pv_remaining, pv.power_priority_load_MW_AC, pv.site.POI_limit] if val is not None])
                #pv to load for this timestep
                df.loc[i+1,pv.name+'_to_load_MWh_e'] = pv2l
                #updating local variable for pv availability
                pv_remaining -= pv2l
                poi_remaining -= pv2l
                pv_poi+= pv2l
                #updating unmet load after pv contribution
                unmet_load_MW -= pv2l

                #battery charging with pv
                for bes in bes_systems:
                    #check if pv is allowed to charge bes
                    if pv.name not in bes.systems_charging: 
                        continue

                    #local variable to track bes charge, checks if charged by another pv system already
                    if df.loc[i+1,bes.name+'_MWh_DC']>df.loc[i,bes.name+'_MWh_DC']:
                        bes_MWh = df.loc[i+1,bes.name+'_MWh_DC']
                        delta_bes = df.loc[i+1,bes.name+'_MWh_DC']-df.loc[i,bes.name+'_MWh_DC']
                    else:
                        bes_MWh = df.loc[i,bes.name+'_MWh_DC']
                        delta_bes=0
                    #calculate unused bes capacity
                    headroom_bes = bes.capacity_MWh_e - bes_MWh

                    #calculate pv contribution to battery, limited by bes charge rate, pv remaining, and headroom
                    if pv.site.name != bes.site.name:
                        poi_lim = poi_remaining
                    else:
                        poi_lim = np.inf
                    pv2batt = max(0,min([val-delta_bes for val in [bes.power_rating_MW_e, pv_remaining+delta_bes, poi_lim+delta_bes, (headroom_bes/bes.charge_efficiency)+delta_bes] if val is not None])) # assumes bes power rating = bes charge rating
                    #update local variable for pv availability, tracking 
                    pv_remaining -= pv2batt
                    if pv.site.name != bes.site.name:
                        poi_remaining -= pv2batt
                        df.loc[i+1,bes.site.name+'_POI_MW']-= pv2batt
                        pv_poi += pv2batt

                    #tracking pv contribution to bes
                    df.loc[i+1,pv.name+'_to_'+bes.name+'_MWh_e'] = pv2batt
                    #bes at this time step after charging contributions
                    df.loc[i+1,bes.name+'_MWh_DC'] = bes_MWh+pv2batt*bes.charge_efficiency
                    
                #charging tes from pv after bes is satisfied
                for tes in tes_systems:
                    #check if charging by pv is allowed
                    if pv.name not in tes.systems_charging:
                        continue

                    #local variable to track tes charge, from updated df variable, checking if updated for multiple charging sources
                    if len(csp_systems)>0:
                        csp_cont = df.loc[i+1,csp_systems[0].__class__.__name__+'_to_'+tes.name+'_MWh_t']
                    else:
                        csp_cont = 0
                    tloss = df.loc[i+1,tes.name+'_thermal_loss_MWh_t']
                    tes_MWh_t = df.loc[i+1,tes.name+'_MWh_t']
                    delta_tes = max(0, tes_MWh_t - (df.loc[i,tes.name+'_MWh_t'] - tloss + csp_cont))

                    if tes.charge_rate_resistive_MW_e is not None:
                        tes_charge_res_remaining = tes.charge_rate_resistive_MW_e - delta_tes
                    elif tes.charge_rate_resistive_MW_e is None:
                        tes_charge_res_remaining = None
                    
                    #recalculating tes headroom
                    headroom_tes = tes.capacity_MWh_t-tes_MWh_t
                    
                    #calculating pv contribution to tes, limited by charge rate, headroom, pv availability
                    if pv.site.name != tes.site.name:
                        poi_lim = poi_remaining
                    else:
                        poi_lim = np.inf
                    pv2tes = min([val for val in [tes_charge_res_remaining, headroom_tes/tes.charge_efficiency_e2TES, pv_remaining, poi_lim] if val is not None])
                    
                    #tracking pv contribution to tes
                    df.loc[i+1,pv.name+'_to_'+tes.name+'_MWh_e'] = pv2tes

                    #update pv still available
                    pv_remaining -= pv2tes
                    #tracking poi contributions
                    if pv.site.name != tes.site.name:
                        poi_remaining -= pv2tes
                        df.loc[i+1,tes.site.name+'_POI_MW']-= pv2tes
                        pv_poi += pv2tes
                    
                    #update tes charge state with pv scaled by efficiency
                    tes_MWh_t += pv2tes*tes.charge_efficiency_e2TES
                    #tes charge state at this timestep after csp and pv contributions
                    df.loc[i+1,tes.name+'_MWh_t'] = tes_MWh_t

                #second chance pv to load if bes and/or TES is full
                if pv_remaining>0:
                    pv2l2 = max(0,min(unmet_load_MW, pv_remaining, poi_remaining))
                    unmet_load_MW -= pv2l2
                    pv_remaining -= pv2l2
                    poi_remaining -= pv2l2
                    pv_poi += pv2l2
                    df.loc[i+1,pv.name+'_to_load_MWh_e'] = df.loc[i+1,pv.name+'_to_load_MWh_e'] + pv2l2
                
                #pv after all has been directed to pv, bes, or tes
                df.loc[i+1,pv.name+'_to_grid_MWh_e'] = min(pv_remaining,max(0,poi_remaining))
                # print(f'curatiled PV: {max([0,pv_remaining-max([0,poi_remaining])])}')
                # df.loc[i+1,pv.name+'_curtailed_MWh_e'] = max([0,pv_remaining-max([0,poi_remaining])])
                df.loc[i+1,pv.name+'_curtailed_MWh_e'] = max([0,min(pv_remaining,max([0,poi_remaining]))]) # updated curtailment 
                # print(pv_remaining-max([0,poi_remaining]))
                df.loc[i+1,pv.site.name+'_POI_MW'] = pv_poi

            #load satisfaction by tes (2nd after pv)
            for tes in tes_systems:
                #local variable to track tes charge state
                tes_MWh_t = df.loc[i+1,tes.name+'_MWh_t']
                if tes.site.POI_limit is not None:
                    poi_remaining = tes.site.POI_limit - df.loc[i+1,tes.site.name+'_POI_MW']
                else:
                    poi_remaining = np.inf
                #how much tes MWh_t are available for use, limited by discharge depth in normal operation
                if operation=='normal':
                    tes_avail = max(0,tes_MWh_t-tes.capacity_MWh_t*(1.-tes.percent_discharge_depth/100.))
                elif operation=='off_grid':
                    tes_avail = max(0,tes_MWh_t)
                #how much tes MWh_e are available for use
                tes_avail_e = tes_avail*tes.efficiency_rating_t2e

                #tes_power_minimum ensures no power cycling for small loads in normal operation
                if operation=='normal':
                    tes_power_minimum = tes.power_minimum_MW_e
                elif operation=='off_grid':
                    tes_power_minimum = 0 #no minimum in off-grid operation

                #only turns on if load meets power minimum
                if unmet_load_MW >= tes.power_minimum_MW_e:
                    #calculating tes to load, limited by tes availability, unmet load thermal equivalent, and tes power rating thermal equivalent
                    tes2load = min([val/tes.efficiency_rating_t2e for val in [tes_avail*tes.efficiency_rating_t2e, unmet_load_MW,tes.power_rating_MW_e,poi_remaining] if val is not None])
                else:
                    tes2load = 0

                #updating local tes variable with energy delivered to load - thermal
                tes_MWh_t -= tes2load
                #df.loc[i+1,'RES_to_load_MWh_th'] = tes2load 
                #updating local load variable with energy delivered by tes - electrical
                unmet_load_MW -= tes2load*tes.efficiency_rating_t2e

                #updating tes at this timestamp after load contribution
                df.loc[i+1,tes.name+'_MWh_t'] = tes_MWh_t
                #tracking load contribution in thermal
                df.loc[i+1,tes.name+'_to_load_MWh_t'] = tes2load
                #tracking load contribution in electrical
                df.loc[i+1,tes.name+'_to_load_MWh_e'] = tes2load*tes.efficiency_rating_t2e
                #tracking poi_limit
                df.loc[i+1,tes.site.name+'_POI_MW'] += tes2load*tes.efficiency_rating_t2e

            # bes contributions to load are last
            for bes in bes_systems:
                 
                if (bes.max_total_power_MW is not None) and operation=='normal':
                    pvs, total = bes.max_total_power_MW
                    if not isinstance(pvs, list):
                        pvs = [pvs]
                    #calculate total power of all counted systems for that timestep
                    p = sum([df.loc[i+1,pv+'_to_load_MWh_e'] for pv in pvs])
                    #maximum bes contribution, ensure nonnegative
                    bes_max = max(0,total-p)
                else:
                    bes_max = bes.power_rating_MW_e
                #local variable to track bes charge state, updated with charging above
                bes_MWh = df.loc[i+1,bes.name+'_MWh_DC']
                #calculate bes available, limited by discharge depth, ensure nonnegative
                #if operation=='normal':
                #    bes_avail = max(0,bes_MWh-bes.capacity_MWh_e*(1.-bes.percent_discharge_depth/100.))
                #elif operation=='off_grid':
                bes_avail = max(0,bes_MWh) #comment out previous 3 lines to account for usable BES capacity changes
                if bes.site.POI_limit is not None:
                    poi_remaining = bes.site.POI_limit - df.loc[i+1,bes.site.name+'_POI_MW']
                else:
                    poi_remaining = np.inf
                #calculating bes to load, limited by availability, unmet load, and bes power rating
                bes2load = max([0,min([val/bes.discharge_efficiency for val in [bes_avail*bes.discharge_efficiency, unmet_load_MW, bes_max, poi_remaining] if val is not None])])
                #updating local variable for bes charge state
                bes_MWh -= bes2load
                #updating local variable to track unmet load
                unmet_load_MW -= bes2load*bes.discharge_efficiency

                #bes charge state at this timestep
                df.loc[i+1,bes.name+'_MWh_DC'] = bes_MWh
                #tracking bes to load at this timestep
                df.loc[i+1,bes.name+'_to_load_MWh_e'] = bes2load*bes.discharge_efficiency
                #tracking poi
                df.loc[i+1,bes.site.name+'_POI_MW'] += bes2load*bes.discharge_efficiency

            #after pv, tes, and csp, any unmet load satisfied by fuel
            df.loc[i+1,'fuel_to_load_MWh_th'] = unmet_load_MW

            if operation=='off_grid':
                if unmet_load_MW>0:
                    lasthour=i
                    break
                elif i+1==endseries:
                    lasthour=i+1
                    to_end=1
                    break
        
        # update fuel usage w unmet load
            df.loc[i + 1, 'fuel_to_load_MWh_th'] = unmet_load_MW 
            df.loc[i + 1, 'Energy provided by fuel [MWh_t]'] = unmet_load_MW

            if unmet_load_MW > 0: 
                unmet_load_MJ = unmet_load_MW * 3600  # MWh to MJ
                fuel_mass_kg = unmet_load_MJ/(fuel_properties['heating_value'] * fuel_properties['efficiency'])
            else:
                unmet_load_MJ = 0
                fuel_mass_kg = 0

            df.loc[i + 1, 'fuel_mass_kg'] = fuel_mass_kg
            df.loc[i + 1, 'fuel_energy_MJ'] = unmet_load_MJ

        df.set_index(self.load_MW.index.name,inplace=True)

        if operation=='normal':
            return df
        elif operation=='off_grid':
            return lasthour, to_end


########################################################################


    def group_sum_metrics(self, metrics, label):
        #summing metrics by system type and full system
        syscols = []
        
        for cat in [self.pv_systems, self.csp_systems, self.bes_systems, self.tes_systems]:
            if len(cat)==0:
                continue
            cols = []
            for sys in cat:
                cols.append(sys.name+'_'+label)
                metrics[sys.name+'_'+label] = getattr(sys,label)

            metrics[sys.__class__.__name__+'_'+label]=sum([metrics[n] for n in cols])

            syscols.append(sys.__class__.__name__+'_'+label)

        metrics['system_'+label] = sum([metrics[n] for n in syscols])

        return metrics

    def storage_power_metrics(self, metrics):
        for sys in self.tes_systems:
            metrics[sys.name+'_charge_rate_resistive_MW_e'] = sys.charge_rate_resistive_MW_e
            metrics[sys.name+'_power_rating_MW_e'] = sys.power_rating_MW_e
        for sys in self.bes_systems:
            metrics[sys.name+'_power_rating_MW_e'] = sys.power_rating_MW_e
        return metrics

    def capacity_metrics(self, metrics):
        #summing capacity metrics by system type
        for cat in [self.pv_systems, self.csp_systems, self.bes_systems, self.tes_systems]:
            if len(cat)==0:
                continue
            cols = []
            types = [i for i in cat[0].__dict__ if 'capacity' in i]
            for type in types:
                for sys in cat:
                    cols.append(sys.name+'_'+type)
                    metrics[sys.name+'_'+type] = getattr(sys, type)
                metrics[sys.__class__.__name__+'_'+type] = sum([metrics[n] for n in cols])
        return metrics
        
    def production_metrics(self, metrics):
        #calculating production metrics for the timeseries
        df = self.timeseries
        #generating list of production columns, assumes contains MWh
        prodcols = [c for c in df.columns if 'MWh' in c]
        #initializing local list to track metrics grouped by generation system but not storage system
        halfcols = []
        for cat in [self.pv_systems, self.csp_systems, self.bes_systems, self.tes_systems]:
            #skip if no systems in category
            if len(cat)==0:
                continue
            #list to track variables for power receipt
            prods = []
            #cycling through systems in the category
            for sys in cat:
                #find all variables that start with that system name
                prods.extend([c.replace(sys.name,'') for c in prodcols if c.startswith(sys.name)])
            #find the list of unique variables
            prods = set(prods)
            #cycle through the unique set
            for prod in prods:
                #tracking the contributing systems for each power receipt variable
                syscols = []   
                for sys in cat:
                    if sys.name+prod in df.columns:
                        syscols.append(sys.name+prod)
                        #cataloging each flow
                        metrics[sys.name+prod]=df[sys.name+prod].sum()
                #cataloging flows from system type to receipt variable
                halfcols.append(sys.__class__.__name__+prod)
                metrics[sys.__class__.__name__+prod] = sum([metrics[n] for n in syscols])

        #grouping receipt variables
        for cat in [self.bes_systems, self.tes_systems]:
            if len(cat)==0:
                continue
            modcols = [h for h in halfcols if any(sys.name in h for sys in cat)]
            for source in [self.pv_systems, self.csp_systems]:
                if len(source)==0:
                    continue
                mets = []
                for sys in cat:
                    mets.extend([m.replace(sys.name,'') for m in modcols if m.startswith(source[0].__class__.__name__+'_to_'+sys.name)])
                for m in mets:
                    syscols = []
                    for sys in cat:
                        syscols.append(source[0].__class__.__name__+'_to_'+sys.name+m.replace(source[0].__class__.__name__+'_to_',''))
                    #cataloging flows from system type to system type
                    metrics[source[0].__class__.__name__+'_to_'+sys.__class__.__name__+m.replace(source[0].__class__.__name__+'_to_','')] = sum([metrics[n] for n in syscols])

        return metrics

    def generation_metrics(self, metrics):
        #summing generation by each system and system type
        df = self.timeseries

        #for pv
        pvcols = []
        for sys in self.pv_systems:
            pvcols.append(sys.name+'_energy_MWh_AC')
            metrics[sys.name+'_energy_MWh_AC'] = df[sys.name+'_power_MW_AC'].sum()
        metrics[sys.__class__.__name__+'_energy_MWh_AC'] = sum([metrics[n] for n in pvcols])

        #for csp
        cspcols = []
        for sys in self.csp_systems:
            cspcols.append(sys.name+'_energy_MWh_t')
            metrics[sys.name+'_energy_MWh_t'] = df[sys.name+'_heat_MW_t'].sum()
        metrics[sys.__class__.__name__+'_energy_MWh_t'] = sum([metrics[n] for n in cspcols])

        return metrics

    def load_satisfaction_annual(self):
        df = self.timeseries.copy()
        conts = []
         
        for sys in [self.pv_systems, self.tes_systems, self.bes_systems]:
            if len(sys) == 0:
                continue
            cols = [c.name + '_to_load_MWh_e' for c in sys]
            df[sys[0].__class__.__name__ + '_to_load_MWh_e'] = df[cols].sum(axis=1)
            conts.append(sys[0].__class__.__name__ + '_to_load_MWh_e')
         
        df['system_to_load_MWh_e'] = df[conts].sum(axis=1)
         
        annual = df['system_to_load_MWh_e'].groupby(by=df.index.year).sum().astype(int).to_list()
         
        hourly = df['system_to_load_MWh_e']
        #print(hourly)
        return annual, hourly


    def LCOE_metrics(self, metrics):

        def MACRS_pvd(input_year, WACC_n):
             
            macrs_rates = {
            3: [0.3333, 0.4445, 0.1481, 0.0741],  # 3-year property
            5: [0.2000, 0.3200, 0.1920, 0.1152, 0.1152, 0.0576],  # 5-year property
            7: [0.1429, 0.2449, 0.1749, 0.1249, 0.0893, 0.0892, 0.0893, 0.0446]  # 7-year property
        }
            if input_year in macrs_rates:
                macrs_schedule = macrs_rates[input_year]
                return sum(rate * npf.pv(rate=WACC_n, nper=year, pmt=0, fv=-1) for year, rate in enumerate(macrs_schedule, start=1))

         
        metrics['ITC'] = self.ITC
        metrics['analysis_period'] = self.analysis_period
        metrics['WACC_n'] = DF*I*(1-tax) + (1-DF)*COE # nom LCOE
        metrics['WACC_r'] = ((1+metrics['WACC_n'])/(1+inflation))-1 # real
        metrics['PVD'] = MACRS_pvd(MACRS_yrs, metrics['WACC_n'])
        metrics['CRF'] = metrics['WACC_r']/(1-(1+metrics['WACC_r'])**(-metrics['analysis_period']))
        metrics['FCR'] =((metrics['CRF']* ((1 -(tax* metrics['PVD'])*(1-metrics['ITC']/2) - metrics['ITC'])) ) + property_tax + insurance )/ (1 - tax)
        # pv_cost = sum(pv.capex_USD for pv in self.pv_systems) #***********************************************
        # metrics['system_capex_USD'] -= pv_cost#***********************************************

         
        fuel_price_per_gal = 1.8 # $/gal
        rho_liq_propane = 493 # kg/m3
        annual_fuel_mass = np.array(metrics['annual_fuel_mass_kg'])
        annual_fuel_gals = (annual_fuel_mass/rho_liq_propane) * 264.172 #m3 to gallons
        annual_fuel_cost = annual_fuel_gals * fuel_price_per_gal
        metrics['annual_fuel_cost_USD'] = annual_fuel_cost


        def augment_array_np(arr, L):
            l_val = arr.size  
            if L > l_val:
                shortfall = L - l_val
                repeat_times = (shortfall + l_val - 1) // l_val
                repeated_section = np.tile(arr, (repeat_times,))[-shortfall:]
                augmented_array = np.concatenate((arr, repeated_section))
                return augmented_array
            else:
                return arr
            
        OM_NPV_arr=[]
        annual_renewables_NPV_arr=[]
        aug_batt_NPV_arr=[]
        aug_htr_NPV_arr=[]
    
        annual_renewables = np.array(metrics['system_to_load_annual_MWh_th'])
        annual_grid = np.array(metrics['fuel_to_load_annual_MWh_th'])

        add = np.zeros(metrics['analysis_period'])
        for bes in self.bes_systems:
            add = add + np.array(augment_array_np(np.array(bes.cpx_batt_aug_annual),metrics['analysis_period'])[:metrics['analysis_period']])
        aug_batt = add
        metrics['system_cpx_batt_aug_annual'] = list(aug_batt)
        
        add = np.zeros(metrics['analysis_period'])
        for tes in self.tes_systems:
            add = add + np.array(augment_array_np(np.array(tes.cpx_htr_aug_annual),metrics['analysis_period'])[:metrics['analysis_period']])
        aug_htr = add
        metrics['system_cpx_htr_aug_annual'] = list(aug_htr)
        
        annual_renewables = augment_array_np(annual_renewables, metrics['analysis_period'])*1000
        annual_fuel_cost = augment_array_np(annual_fuel_cost, metrics['analysis_period'])
        annual_grid = augment_array_np(annual_grid,metrics['analysis_period'])*1000*3.29*21.63/1000
        annual_electricity_sales = (metrics['export_energy_MWh_e']*self.e_sale * 1000) # get curtailed PV, set e_sale price 
        
        total_energy_supplied = np.array(metrics['load_annual_MWh']) 
        annual_load = augment_array_np(total_energy_supplied, metrics['analysis_period'])*1000

        # Net Present Value Over Analysis Period
        NPV_batt_ARMO_N = npf.npv(metrics['WACC_r'], aug_batt[:metrics['analysis_period']])
        NPV_htr_ARMO_N = npf.npv(metrics['WACC_r'], aug_htr[:metrics['analysis_period']])
        NPV_renewables_N = npf.npv(metrics['WACC_r'], annual_renewables[:metrics['analysis_period']])
        NPV_grid_N = npf.npv(metrics['WACC_r'], annual_grid[:metrics['analysis_period']])
        NPV_electricity_sales_N = npf.npv(metrics['WACC_r'], annual_electricity_sales)  #********* 
        NPV_fuel_cost_N = npf.npv(metrics['WACC_r'], annual_fuel_cost[:metrics['analysis_period']])
        NPV_load_N = npf.npv(metrics['WACC_r'], annual_load[:metrics['analysis_period']])
    
        
        OM_NPV_arr = []
        for v in range(metrics['analysis_period']):
            OM_esc = metrics['system_annual_OM_USD'] * ((1 + esc) ** v)
            OM_NPV_arr.append(OM_esc)
        #OM_NPV_arr[0] = 0
        NPV_OM_N = npf.npv(metrics['WACC_r'], np.array(OM_NPV_arr))
         
        NPV_batt_ARMO_L = npf.npv(metrics['WACC_r'], aug_batt)
        NPV_htr_ARMO_L = npf.npv(metrics['WACC_r'], aug_htr)
        NPV_renewables_L = npf.npv(metrics['WACC_r'], annual_renewables)
        NPV_grid_L = npf.npv(metrics['WACC_r'], annual_grid)
        NPV_electricity_sales_L = npf.npv(metrics['WACC_r'], annual_electricity_sales)
        NPV_fuel_cost_L = npf.npv(metrics['WACC_r'], annual_fuel_cost)
    

        OM_NPV_arr = []
        for v in range(len(aug_batt)): 
            OM_esc = metrics['system_annual_OM_USD'] * ((1 + esc) ** v)
            OM_NPV_arr.append(OM_esc)
        OM_NPV_arr[0] = 0
        NPV_OM_L = npf.npv(metrics['WACC_r'], np.array(OM_NPV_arr))
    
        annualized_CAPEX = metrics['system_capex_USD'] * metrics['FCR'] * (1-grant_frac)
        #print(f"system capex {metrics['system_capex_USD']}")
        annualized_OM = NPV_OM_N * metrics['CRF']
        annualized_ARMO = (NPV_batt_ARMO_N + NPV_htr_ARMO_N) * metrics['CRF']
        annualized_grid = NPV_grid_N * metrics['CRF']
        annualized_fuel = NPV_fuel_cost_N * metrics['CRF']
        annualized_sales = NPV_electricity_sales_N * metrics['CRF'] *-1

        annual_ARR = annualized_CAPEX + annualized_OM + annualized_ARMO + annualized_fuel + annualized_sales

        # Calculate residual value for longer project life compared to analysis period (e.g. = 0 for n=L) 
        Rv = ((((1 + metrics['WACC_r']) ** metrics['analysis_period']) * ((1 - NPV_renewables_N / NPV_renewables_L) * (metrics['system_capex_USD'] * (1 - (tax * metrics['PVD']) * (1 - metrics['ITC'] / 2) - metrics['ITC'])) + (NPV_OM_N + NPV_htr_ARMO_N + NPV_batt_ARMO_N) -
(NPV_renewables_N / NPV_renewables_L) * (NPV_OM_L + NPV_htr_ARMO_L + NPV_batt_ARMO_L)))) / ((1 + metrics['WACC_r']) ** metrics['analysis_period'])
         
        cash_flow = np.ones(metrics['analysis_period']) * annual_ARR
        NPV_ARR = npf.npv(metrics['WACC_r'], cash_flow)
         
          
        if NPV_renewables_N == 0:
            metrics['LCOE_real_USD_kWh'] = (NPV_ARR) / NPV_load_N
            #print(NPV_ARR,NPV_fuel_cost_N,total_energy_supplied)
            #print(f'NPV_ARR: {NPV_ARR}')
            #print(f'NPV_load_N: {NPV_load_N}')
        else:
            #metrics['LCOE_real_USD_kWh'] = (NPV_ARR - Rv) / NPV_renewables_N 
            metrics['LCOE_real_USD_kWh'] = (NPV_ARR) / NPV_load_N
         
        metrics['LCOE_residual_value'] = Rv
        metrics['NPV_fuel_cost_USD'] = NPV_fuel_cost_L

        return metrics
   
    def system_metrics(self):
        df = self.timeseries.copy()
        df.index = pd.to_datetime(df.index)
   
        metrics = {}

        metrics['hourly_PV_to_load_MWh'] = df[[pv.name+'_to_load_MWh_e' for pv in self.pv_systems]].sum(axis=1).tolist()
        metrics['hourly_TES_to_load_MWh'] = df[[tes.name+'_to_load_MWh_e' for tes in self.tes_systems]].sum(axis=1).tolist()
        metrics['hourly_Fuel_to_load_MWh'] = df['fuel_to_load_MWh_th'].tolist()
 
        metrics['PV_systems'] = [pv.name for pv in self.pv_systems]
        metrics['CSP_systems'] = [csp.name for csp in self.csp_systems]
        metrics['BES_systems'] = [bes.name for bes in self.bes_systems]
        metrics['TES_systems'] = [tes.name for tes in self.tes_systems]
 
        metrics['years'] = len(df) / (24 * 365)
 
        pv_capex = sum([getattr(pv, 'capex_USD', 0) for pv in self.pv_systems])
        csp_capex = sum([getattr(csp, 'capex_USD', 0) for csp in self.csp_systems])
        tes_capex = sum([getattr(tes, 'capex_USD', 0) for tes in self.tes_systems])
        bes_capex = sum([getattr(bes, 'capex_USD', 0) for bes in self.bes_systems])
 
        metrics['PV_capex_USD'] = pv_capex
        metrics['CSP_capex_USD'] = csp_capex
        metrics['TES_capex_USD'] = tes_capex
        metrics['BES_capex_USD'] = bes_capex
        metrics['system_capex_USD'] = pv_capex + csp_capex + tes_capex + bes_capex
 
        #print(f"CAPEX Details:")
        #print(f"  PV CAPEX: {pv_capex}")
        #print(f"  CSP CAPEX: {csp_capex}")
        #print(f"  TES CAPEX: {tes_capex}")

        #systems that work off-grid
        metrics['off_grid_systems'] = [sys.name for sys in list(chain.from_iterable([self.pv_systems_off_grid, self.csp_systems_off_grid,self.bes_systems_off_grid, self.tes_systems_off_grid]))]
 
        metrics['years'] = len(df)/(24*365)

        metrics['load_MWh'] = df.load_MW.sum()
        metrics['load_annual_MWh'] = df.load_MW.groupby(df.index.year).sum().astype(int).to_list()

        metrics = self.capacity_metrics(metrics)
        metrics = self.storage_power_metrics(metrics)
        
        for label in ['area_land_acres','capex_USD','annual_OM_USD']:
            metrics = self.group_sum_metrics(metrics, label)
            
        metrics = self.generation_metrics(metrics)
        
        metrics = self.production_metrics(metrics)

        metrics['fuel_to_load_MWh_th'] = df.fuel_to_load_MWh_th.sum()
        #############################################
        #metrics['RES_to_load_MWh_th'] = df.RES_to_load_MWh_th.sum() 

        metrics['fuel_to_load_annual_MWh_th'] = df.fuel_to_load_MWh_th.groupby(df.index.year).sum().astype(int).to_list()
        metrics['system_to_load_annual_MWh_th'], metrics['hourly_RES_to_load'] = self.load_satisfaction_annual()
   
        load_series = df.load_MW
        
        hourly_RES_LMP = 100 * np.array(metrics['hourly_RES_to_load']) / np.array(load_series)
        #print(hourly_RES_LMP)
        metrics['hourly_RES_LMP'] = hourly_RES_LMP
        metrics['annual_fuel_energy_MJ'] = df['fuel_energy_MJ'].groupby(df.index.year).sum().astype(int).to_list()
        metrics['annual_fuel_mass_kg'] = df['fuel_mass_kg'].groupby(df.index.year).sum().astype(int).to_list()

        # Calculate the start and end indices for March 1 and June 1
        start_index = 60 * 24  # March 1
        end_index = 151 * 24  # June 1
        load_MWh_Mar1_to_Jun1 = df.load_MW[start_index:end_index].sum()
        fuel_to_load_Mar1_to_Jun1 = df.fuel_to_load_MWh_th[start_index:end_index].sum()
        pct_sys_window = 100*(1-(fuel_to_load_Mar1_to_Jun1/load_MWh_Mar1_to_Jun1))
        metrics['pct_sys_window'] = pct_sys_window

        metrics['percent_load_by_grid'] = 100*metrics['fuel_to_load_MWh_th']/metrics['load_MWh']
    
        #print(f"load_MWh: {metrics['load_MWh']}")
        #print(f"fuel to load: {metrics['fuel_to_load_MWh_th']}")
        #print(f"system_to_load_annual_MWh_th: {metrics['system_to_load_annual_MWh_th']}")
        #print(metrics['load_MWh'],metrics['fuel_to_load_MWh_th'],metrics['system_to_load_annual_MWh_th'])
        #metrics['percent_load_by_system_net'] = 100*(metrics['load_MWh'] - metrics['fuel_to_load_MWh_th'] + metrics['system_to_load_annual_MWh_th'])/metrics['load_MWh']
        metrics['percent_load_by_system_net'] = 100*(metrics['system_to_load_annual_MWh_th'][0])/metrics['load_MWh']
        #print(f"percent_load_by_system_net: {metrics['percent_load_by_system_net']}")
        
        metrics['percent_time_50pct_load'] = 100*((df.fuel_to_load_MWh_th/df.load_MW)<.5).sum()/len(df)
        
        
        # for metric in metrics:
        #     if isinstance(metrics[metric], float):
        #         metrics[metric] = int(round(metrics[metric],0))
        
        # List TEA Metrics
        metrics = self.LCOE_metrics(metrics)
        
        return metrics


    def randomstart(self, critical_load_mw):

        start = randint(1, len(self.timeseries)) #start anywhere
    
        df_res = pd.concat([self.timeseries.copy().iloc[start::,:],self.timeseries.copy().iloc[::start,:]]) #loop through the full timeseries
    
        res_results = {}
    
        lasthour_actual, to_end_actual = self.operation(operation='off_grid', critical_load_MW = 25, df_res = df_res, res_initial = 'actual')
    
        lasthour_full, to_end_full = self.operation(operation='off_grid', critical_load_MW = 25, df_res = df_res, res_initial = 'full')
    
        res_results['start_hour']=start
        res_results['actual_duration_hours']=lasthour_actual
        res_results['full_duration_hours']=lasthour_full
        res_results['to_end_actual']=to_end_actual
        res_results['to_end_full']=to_end_full
            
        return res_results

    def resilience_cases(self, critical_load_MW, n_starts = 100):
        
        #random sampling and resiliency analysis
        df_summary = pd.DataFrame(columns=['start_hour','actual_duration_hours','full_duration_hours','to_end_actual','to_end_full'])
    
        async_results = []
        with mp.Pool(processes=mp.cpu_count()) as pool:
            for _ in range(n_starts):
                async_results.append(pool.apply_async(self.randomstart, args=[critical_load_MW]))
            pool.close()
            pool.join()
        
        for num, async_result in enumerate(async_results):
            df_summary=pd.concat([df_summary, pd.DataFrame(async_result.get(),index=[int(num),])])
        
        return df_summary


    def resilience_metrics(self, df ,target_hours):
        results = {}
        results['pct_meets_actual'] = (df['actual_duration_hours']>=target_hours).sum()/len(df)
        results['pct_meets_full']=(df['full_duration_hours']>=target_hours).sum()/len(df)
        results['tenth_pctile_actual_hrs'] = df['actual_duration_hours'].quantile(.1)
        results['fiftieth_pctile_actual_hrs'] = df['actual_duration_hours'].quantile(.5)
        results['tenth_pctile_full_hrs'] = df['full_duration_hours'].quantile(.1)
        results['fiftieth_pctile_full_hrs'] = df['full_duration_hours'].quantile(.5)
        results['actual_reaches_end_pct'] = df['to_end_actual'].sum()/len(df)
        results['full_reaches_end_pct'] = df['to_end_full'].sum()/len(df)
        
        return results

    def timeseries_plot_source(self, start_date='random', days=365*2, type='area', color_seq=None, system_seq=None):
        # Plotting contributions to load by individual systems
    
        df = self.timeseries.copy()
    
        # Starting at random midnight
        ndays = len(df) / 24
        if start_date == 'random':
            start_date = df.index[randrange(int(ndays - days)) * 24]
        else:
            start_date = pd.to_datetime(start_date)
        end_date = start_date + pd.Timedelta(days=days)
    
        df2 = df.loc[start_date:end_date].copy()
    
        if system_seq is None:
            conts = [c.name + '_to_load_MWh_e' for c in list(chain.from_iterable([self.pv_systems, self.tes_systems, self.bes_systems]))]
            conts.append('Energy provided by fuel [MWh_t]')
        else:
            conts = system_seq
     
        df2['Energy provided by fuel [MWh_t]'] = df2['Energy provided by fuel [MWh_t]'].abs()
        sum_grid = np.sum(df2['Energy provided by fuel [MWh_t]'])
        sum_pv = np.sum(df2['ASGARD_PV_to_load_MWh_e'])
        sum_TES = np.sum(df2['ASGARD_TES_to_load_MWh_e'])
        sum_total = sum_grid + sum_pv + sum_TES
    
        fig, ax = plt.subplots()
        ax.plot(df2[['load_MW']], linestyle='-', color='black')
    
        if type == 'bar':
            df2[conts].plot(kind='bar', stacked=True, width=1, ax=ax)
    
        elif type == 'area':
            if color_seq is None:
                df2[conts].plot.area(stacked=True, ax=ax, lw=0)
            else:
                df2[conts].plot.area(stacked=True, ax=ax, lw=0, color=color_seq)
    
        df2[['load_MW']].plot(linestyle='-', color='black', ax=ax, label=None)
    
        plt.ylabel('Hourly Energy (MWh)')
        plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
        plt.show()


    def timeseries_plot_group(self,start_date='random',days=7):
        #plotting contributions to load by system group (pv, tes, bes)
    
        df=self.timeseries.copy()
        
        ndays = len(df)/24
        if start_date == 'random':
            start_date = df.index[randrange(ndays-days)*24]
        else:
            start_date = pd.to_datetime(start_date)
        end_date = start_date + DT.timedelta(days=days)

        df2 = df.loc[start_date:end_date].copy()
        conts = []
        for sys in [self.pv_systems, self.tes_systems, self.bes_systems]:
            if len(sys)==0:
                continue
            cols = [c.name+'_to_load_MWh_e' for c in sys]
            df2[sys[0].__class__.__name__+'_to_load_MWh_e'] = df2[cols].sum(axis=1)
            conts.append(sys[0].__class__.__name__+'_to_load_MWh_e')
        
        #df.load_MW.plot()
        df2.plot.area(y=conts, stacked=True,use_index=True,lw=0)
        plt.plot(df2['load_MW'],color='black')
        plt.xlim(start_date,end_date)
        plt.ylabel('Hourly Energy (MWh)')
        plt.show()


    def plot_tes_capacity(self, start_date='random', days=365*2, type='area', color_seq=None, system_seq=None):
        # Plotting TES capacity over time
    
        df = self.timeseries.copy()
    
        # Starting at random midnight
        ndays = len(df) / 24
        if start_date == 'random':
            start_date = df.index[randrange(int(ndays - days)) * 24]
        else:
            start_date = pd.to_datetime(start_date)
        end_date = start_date + pd.Timedelta(days=days)
    
        df2 = df.loc[start_date:end_date].copy()
    
        if system_seq is None:
            conts = [tes.name + '_MWh_t' for tes in self.tes_systems]
        else:
            conts = system_seq
        
        fig, ax = plt.subplots()
        if type == 'bar':
            df2[conts].plot(kind='bar', stacked=False, width=1, ax=ax)
    
        elif type == 'area':
            if color_seq is None:
                df2[conts].plot.area(stacked=False, ax=ax, lw=0)
            else:
                df2[conts].plot.area(stacked=False, ax=ax, lw=0, color=color_seq)
    
        plt.ylabel('Stored Energy (MWh_thermal)')
        plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
        plt.show()

    import matplotlib.pyplot as plt

    def plot_fuel_mass(self, start_date='random', days=365 * 2):

        df = self.timeseries.copy()

        ndays = len(df) / 24
        if start_date == 'random':
            start_date = df.index[randrange(int(ndays - days)) * 24]
        else:
            start_date = pd.to_datetime(start_date)
        end_date = start_date + pd.Timedelta(days=days)

        df2 = df.loc[start_date:end_date].copy()

        if 'fuel_mass_kg' in df2:
            fig, ax = plt.subplots(figsize=(10, 5))

            ax.plot(df2.index, df2['fuel_mass_kg'], color='red', linewidth=2)

            plt.xlabel('Time', fontsize=12)
            plt.ylabel('Fuel Mass (kg)', fontsize=12)

            plt.show()
        else:
            print("check")
