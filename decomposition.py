class decomposition():
    def decomposition(df, timeslice, location = None):
        
        """
        This function breaks down the change in [H+] into various drivers:
            - Temperature
            - Salinity
            - TA mixing
            - TA bio
            - DIC mixing
            - DIC bio
            - Total mixing
            - Total bio
            - DIC air sea flux
        The function also provides the error for each component.
        
        Required Inputs:
            - DataFrame (df), which includes
                - "Temperature [degC]"
                - "Salinity [PSU]"
                - "pH final"
                - "DateTime"
            - 'timeslice' or the period over which you average
                - 'D' for daily averages
                - 'W' for weekly averages
                - 'M' for monthly averages
                - 'Y' for yearly averages
                
        Option Inputs:
            - Column in DataFrame for "avgWindSpeed"
                - If windspeed not provided, 3 m/s will be used.
        
        Outputs (in order):
            - DataFrame containing the individual components of the change in [H+]
                - 'DateTime' 
                - 'Temperature'
                - 'Salinity'
                - 'TA mixing'
                - 'TA bio',
                - 'DIC mixing'
                - 'DIC_bio'
                - 'DIC air sea flux'
                - 'Total bio',
                - 'Total mixing'
                - 'Total'
            - DataFrame containing the errors for the respective components
                - 'DateTime' (not an error, just the index)
                - 'Temperature'
                - 'Salinity'
                - 'TA mixing'
                - 'TA bio',
                - 'DIC mixing'
                - 'DIC bio'
                - 'DIC ASF'
                - 'Total bio',
                - 'Total mixing'
                - 'Total error'
        
        This function employs supplementary functions, modified from MATLAB functions written by Zelun Wu 
        (University of Delaware & Xiamen University; zelunwu@outlook.com; zelunwu.github.io).
            - kt = co_gas_transfer_velocity(sst,wspd,c,unit)
            - K0 = co_K0_Weiss(sst,sss)
            - co2_flux = co_co2flux(pco2_sea, pco2_air, sst, sss, wspd)
            
        References:
            - Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. 
            Marine Chemistry, 2(3), 203–215. https://doi.org/10.1016/0304-4203(74)90015-2.
            - Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. 
            Limnology and Oceanography: Methods, 12(6), 351–362. https://doi.org/10.4319/lom.2014.12.351.
            - Pimenta, A., Oczkowski, A., McKinney, R., & Grear, J., 2023. 
            Geographical and seasonal patterns in the carbonate chemistry of Narragansett Bay, RI, 
            Regional Studies in Marine Science, 62, 102903.
            
        Author:
        Abby Baskind (they/she)
        PhD Student
        URI Graduate School of Oceanography
        Wang Lab/Ocean Carbon Group
        abaskind@uri.edu
        GitHub: abby_baskind
        30 July 2023
        """
        # Basic packages
        import numpy as np
        from matplotlib import pyplot as plt
        import xarray as xr
        import pandas as pd
        import math

        # DateTime packages
        from matplotlib.dates import DateFormatter
        from datetime import datetime, timedelta
        import time
        import matplotlib.dates as mdates

        # Stats/science packages
        import scipy
        import PyCO2SYS as pyco2
        import gsw
        
    
        # Check if `timeslice` is an appropriate selection-------------------------------------------------
        if timeslice != 'D' and timeslice != 'W' and timeslice != 'Y' and timeslice != 'M':
            print('Please select "D," "W," "M," or "Y" for day, week, month, or year.')
            return
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # Check that the DataFrame has all the necessary inputs--------------------------------------------
        # # DateTime
        if not 'DateTime' in df.keys():
            print('Please ensure your dataframe has a column named "DateTime".')
            return
        # # Temperature
        if not 'Temperature [degC]' in df.keys():
            print('Please ensure your dataframe has a column for temperature named "Temperature [degC]".')
            return
        # # Salinity
        if not 'Salinity [PSU]' in df.keys():
            print('Please ensure your dataframe has a column for salinity named "Salinity [PSU]".')
            return
        # # pH
        if not 'pH final' in df.keys():
            print('Please ensure your dataframe has a column for pH named "pH final".')
            return
        # # Location
        if location == 'PLT':
            MLD = 7
        elif location == 'GB':
            MLD = 9
        else:
            MLD = 8
        # -------------------------------------------------------------------------------------------------
    
        # ************************************************************************************************* 
    
        # Take daily/weekly/monthly/yearly mean and standard error-----------------------------------------
        df['Alkalinity [umol/kg]'] = 477.62 + 51.99 * df['Salinity [PSU]']
        MN = df.resample(timeslice, on='DateTime').mean()            # MN contains daily/monthly means
        STD = df.resample(timeslice, on='DateTime').std(ddof = 1)    # STD contains stan devs of means
        count = df.resample(timeslice, on='DateTime').count()        # N for each mean
        COUNT = count['DateTime']
        # -------------------------------------------------------------------------------------------------
    
        # ************************************************************************************************* 
    
        # Optional input for wind speed--------------------------------------------------------------------
        # # If wind speed not provided, default U10 to t m/s
        if not 'avgWindSpeed' in MN.keys():
            print('Your dataframe does not include a column for wind speed labeled as "avgWindSpeed." Default U10 wind speed is set to 5 m/s.')
            U10 = 5
        # # If wind speed is provided, use that for U10 speed
        else:
            U10 = MN['avgWindSpeed']
        # -------------------------------------------------------------------------------------------------
    
        # ************************************************************************************************* 
    
        # pH ERROR-----------------------------------------------------------------------------------------
        # # Calculate error in pH arising from both instrumental error and error from averaging
        STD['pH instrumemnt error'] = 0.1
        ERR_PH = 0.1/(np.sqrt(COUNT))
        STD['pH total error'] = np.sqrt(STD['pH final']**2 + ERR_PH**2)
        # np.sqrt(STD['pH final']**2 + STD['pH instrumemnt error']**2)
        # -------------------------------------------------------------------------------------------------
        
        # Temperature + Salinity ERROR---------------------------------------------------------------------
        # # Calculate error in T and S arising from both instrumental error and error from averaging
        ERR_T = 0.01/(np.sqrt(COUNT))
        STD['Temperature [degC]'] = np.sqrt(STD['Temperature [degC]']**2 + ERR_T**2)
        ERR_S = 0.01/(np.sqrt(COUNT))
        STD['Salinity [PSU]'] =np.sqrt(STD['Salinity [PSU]']**2 + ERR_T**2)
        # -------------------------------------------------------------------------------------------------
        
        # TA ERROR-----------------------------------------------------------------------------------------
        # # Calculate error in TA arising from both instrumental error and error from averaging
        ERR_TA = 10/(np.sqrt(COUNT))
        STD['Alkalinity [umol/kg]'] = np.sqrt(STD['Alkalinity [umol/kg]']**2 + ERR_TA**2)
        # -------------------------------------------------------------------------------------------------
    
        # ************************************************************************************************* 
    
        # SOLVE CARBONATE SYSTEM---------------------------------------------------------------------------
        pH = MN['pH final']                             # pH
        dpH = STD['pH total error']                     # pH error
        dT = STD['Temperature [degC]']                  # Temperature error
        dS = STD['Salinity [PSU]']                      # Salinity error 
        S = MN['Salinity [PSU]']                        # Salinity
        T = MN['Temperature [degC]']                    # Temperature
    
        # Alkalinity calculation   
        TA = MN['Alkalinity [umol/kg]']                 # Alkalinity
        dTA = STD['Alkalinity [umol/kg]']               # Alkalinity error (arising from salinity average)

        results = pyco2.sys(par1 = TA, par2 = pH, par1_type = 1, par2_type = 3, temperature = T, salinity = S,
                            # Calculate uncertainities of...
                            uncertainty_into =["alkalinity", "dic","Hfree", 'pCO2'],
                            # Sources of uncertainty: pH, T, S, and TA
                            uncertainty_from ={"par2": dpH, 'temperature': dT, 'salinity': dS, 'par1': dTA})
    
        MN['H+ [umol/kg]'] = results['Hfree']           # H+ concentration
        MN['DIC [umol/kg]'] = results['dic']            # DIC
        MN['pCO2 [ppm]'] = results['pCO2']              # pCO2

        STD['H+ [umol/kg]'] = results['u_Hfree']        # H+ error (from TA, pH, S, and T)
        STD['DIC [umol/kg]'] = results['u_dic']         # DIC error (from TA, pH, S, and T)
        STD['pCO2'] = results['u_pCO2']                 # pCO2 error (from TA, pH, S, and T)
        # -------------------------------------------------------------------------------------------------
    
        # ************************************************************************************************* 

        # Re-solve carb system using TA and DIC (total)----------------------------------------------------
        # # Calculate gradients of H+
        DIC = MN['DIC [umol/kg]']
        dDIC = STD['DIC [umol/kg]']
        results = pyco2.sys(par1 = DIC, par2 = TA, par1_type = 2, par2_type = 1, temperature = T, salinity = S,
                            # Take grad of H+
                            grads_of=["Hfree"],
                            # In terms of par1 (which is DIC)
                            grads_wrt=["par1", 'temperature', 'salinity', 'par2'],
                            uncertainty_into =["pCO2", "Hfree"],
                            # Uncertainty from DIC already includes uncertainty from T, S, pH, and TA
                            uncertainty_from ={"par1": dDIC})

        MN['∂H/∂DIC'] = results['d_Hfree__d_par1']      # ∂H/∂DIC
        MN['∂H/∂TA'] = results['d_Hfree__d_par2']       # ∂H/∂TA
        MN['∂H/∂T'] = results['d_Hfree__d_temperature'] # ∂H/∂T
        MN['∂H/∂S'] = results['d_Hfree__d_salinity']    # ∂H/∂S
    
        # Error in H+ arising from uncertainty in TA, DIC, T, S, and pH
        STD['H+ [umol/kg]'] = results['u_Hfree']
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # Finite differences-------------------------------------------------------------------------------
        MN['∆DIC'] = MN['DIC [umol/kg]'].diff()
        MN['∆TA'] = MN['Alkalinity [umol/kg]'].diff()
        MN['∆T'] = MN['Temperature [degC]'].diff()
        MN['∆S'] = MN['Salinity [PSU]'].diff()
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # DIC mixing---------------------------------------------------------------------------------------
        # # Based on the relationship presented in Pimenta et al (2023)
        # # 397.65 + 50.59 * S = DIC
        MN['∆DIC_mix'] = (MN['Salinity [PSU]'] * 50.59).diff()
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # DIC air sea flux---------------------------------------------------------------------------------
        # # t for conversion based on timeslice
        if timeslice == 'D':
            t = 365.25
        elif timeslice == 'M':
            t = 12
        elif timeslice == 'W':
            t = 52
        
        # # Requires functions 
        # #        - co2flux(pco2_sea, pco2_air, sst, sss, wspd)
        # #        - K0_Weiss(sst,sss)
        # #        - gas_transfer_velocity(sst,wspd)
    
        # molC/m2/yr
        fgco2_ann = decomposition.co2flux(MN['pCO2 [ppm]'], 410, MN['Temperature [degC]'], MN['Salinity [PSU]'], U10)
        # molC/m2/month or molC/m2/day or...
        fgco2_month = fgco2_ann/t
        # Assumed mixed layer depth of 9m
        # MLD = 9
        # Density [kg/m3]
        rho = gsw.rho(MN['Salinity [PSU]'].to_numpy(), MN['Temperature [degC]'].to_numpy(), np.zeros(len(MN)))
        # umolC/kg/month
        MN['DIC_flux [umolC/L]'] = -(fgco2_month/MLD)*(1e6)/rho
        # finite difference
        MN['∆DIC_flux [umolC/L]'] = MN['DIC_flux [umolC/L]'].diff()
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************

        # DIC bio------------------------------------------------------------------------------------------
        MN['∆DIC_bio [umolC/L]'] = MN['∆DIC'] - MN['∆DIC_flux [umolC/L]'] - MN['∆DIC_mix']
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # TA bio-------------------------------------------------------------------------------------------
        # Converts DIC bio to TA bio using the Redfield ratio (i.e. -16/107)
        MN['∆TA_bio [umolC/L]']=MN['∆DIC_bio [umolC/L]']*(-16/107)
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # TA mixing----------------------------------------------------------------------------------------
        # Based on the relationship presented in Pimenta et al (2023)
        # 477.62 + 51.99 * S = TA
        MN['∆TA_mix [umolC/L]'] = (MN['Salinity [PSU]'] * 51.59).diff()
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # Calculate components-----------------------------------------------------------------------------
        temp = MN['∂H/∂T'] * MN['∆T']                          # Temperature: ∂H/∂T * ∆T
        sal = MN['∂H/∂S'] * MN['∆S']                           # salinity: ∂H/∂S * ∆S
        alk_mix = MN['∂H/∂TA'] * MN['∆TA_mix [umolC/L]']       # TA mixing: ∂H/∂TA * ∆TAmix
        alk_bio = MN['∂H/∂TA'] * MN['∆TA_bio [umolC/L]']       # TA bio: ∂H/∂TA * ∆TAbio
        dic_bio = MN['∂H/∂DIC'] * MN['∆DIC_bio [umolC/L]']     # DIC bio: ∂H/∂DIC * ∆DICbio
        ASF = MN['∂H/∂DIC'] * MN['∆DIC_flux [umolC/L]']        # DIC air sea flux: ∂H/∂DIC * ∆DICflux
        dic_mix = MN['∂H/∂DIC'] * MN['∆DIC_mix']               # DIC mixing: ∂H/∂DIC * ∆DICmix
        BIO = dic_bio + alk_bio                                # Total bio = DIC bio + TA bio
        MIX = dic_mix + alk_mix                                # Total mixing = DIC mixing + TA mixing
        TOT = temp + sal + ASF + BIO + MIX                     # Total
    
        # DataFrame containing all the components
        data = {'DateTime': MN.index,
                'Temperature': temp,
                'Salinity': sal,
                'TA mixing': alk_mix,
                'TA bio': alk_bio,
                'DIC mixing': dic_mix,
                'DIC_bio': dic_bio,
                'DIC air sea flux': ASF,
                'Total bio': BIO,
                'Total mixing': MIX,
                'Total': TOT}
        dH_component = pd.DataFrame(data)
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        # Calculate errors of components-------------------------------------------------------------------
        YERR = {'DateTime': [],
                'Temperature': [],
                'Salinity': [],
                'TA mixing': [],
                'TA bio': [],
                'DIC mixing': [],
                'DIC bio': [],
                'DIC ASF': [],
                'Total bio': [],
                'Total mixing': [],
                'Total error': []}
        # Error of each component
        # ∆z = z √[(∆x/x)^2 + (∆y/y)^2]
        YERR['Temperature'] = temp * np.sqrt((STD['Temperature [degC]']/MN['Temperature [degC]'])**2
                                             + (STD['H+ [umol/kg]']/MN['H+ [umol/kg]'])**2)
        YERR['Salinity'] = sal * np.sqrt((STD['Salinity [PSU]']/MN['Salinity [PSU]'])**2
                                         + (STD['H+ [umol/kg]']/MN['H+ [umol/kg]'])**2)
        YERR['TA mixing'] = alk_mix * np.sqrt((STD['Alkalinity [umol/kg]']/MN['Alkalinity [umol/kg]'])**2
                                              + (STD['H+ [umol/kg]']/MN['H+ [umol/kg]'])**2)
        YERR['TA bio'] = alk_bio * np.sqrt((STD['Alkalinity [umol/kg]']/MN['Alkalinity [umol/kg]'])**2
                                           + (STD['H+ [umol/kg]']/MN['H+ [umol/kg]'])**2)
        YERR['DIC bio'] = dic_bio * np.sqrt((STD['DIC [umol/kg]']/MN['DIC [umol/kg]'])**2
                                            + (STD['H+ [umol/kg]']/MN['H+ [umol/kg]'])**2)
        YERR['DIC mixing'] = dic_mix * np.sqrt((STD['DIC [umol/kg]']/MN['DIC [umol/kg]'])**2
                                               + (STD['H+ [umol/kg]']/MN['H+ [umol/kg]'])**2)
        YERR['DIC ASF'] = ASF * np.sqrt((STD['DIC [umol/kg]']/MN['DIC [umol/kg]'])**2
                                        + (STD['H+ [umol/kg]']/MN['H+ [umol/kg]'])**2)
        YERR['Total bio'] = np.sqrt(YERR['DIC bio']**2 + YERR['TA bio']**2)
        YERR['Total mixing'] = np.sqrt(YERR['DIC mixing']**2 + YERR['TA mixing']**2)
        YERR['DateTime'] = MN.index
    
        TOTERR_ = 0
        for k in YERR.keys():
            if k != 'Total bio' and k != 'Total mixing' and k!= 'Total error' and k != 'DateTime':
                TOTERR_ += YERR[k]**2
        TOTAL_ERROR = np.sqrt(TOTERR_)
        YERR['Total error'] = TOTAL_ERROR
    
        dH_errors = pd.DataFrame(YERR)
        # -------------------------------------------------------------------------------------------------
    
        # *************************************************************************************************
    
        return dH_component, dH_errors
    
    def co2flux(pco2_sea, pco2_air, sst, sss, wspd):
        """
        This function calculates the CO2 surface flux (positive out of the water column) from:
            - pCO2 ocean [uatm] or [ppm]
            - pCO2 air [uatm] or [ppm]
            - sea surface temperature [°C]
            - sea surface salinity [PSU]
            - wind speed [m/s]
        Inputs
            - pCO2_sea [uatm] or [ppm]
            - pCO2_air [uatm] or [ppm]
            - SST [°C]
            - SSS [PSU]
            - wind speed [m/s?]
        Output
            - CO2 upward flux [molC/m2/yr] 
            
        This function is modified from a MATLAB function called co_co2flux(pco2_sea, pco2_air, sst, sss, wspd).
        Original MATLAB function author:
            Zelun Wu
            University of Delaware & Xiamen University
            zelunwu@outlook.com
            zelunwu.github.io
            
        Author:
        Abby Baskind (they/she)
        PhD Student
        URI Graduate School of Oceanography
        Wang Lab/Ocean Carbon Group
        abaskind@uri.edu
        GitHub: abby_baskind
        30 July 2023
        """
        kt = decomposition.gas_transfer_velocity(sst,wspd)
        K0 = decomposition.K0_Weiss(sst,sss)
        dpco2 = pco2_sea - pco2_air
        co2_flux = kt * K0 * dpco2 * (24*365/100000) 
        return co2_flux
    
    def K0_Weiss(sst,sss):
        """
        A function calculated the CO2 solubility with SST [°C] and SSS [PSU].
        Uncertainty is 2% according to Weiss (1974).
        
        Output:
            - K0: CO2 solubility, unit: mol/l/atm
            
        References:
            - Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. 
            Marine Chemistry, 2(3), 203–215. https://doi.org/10.1016/0304-4203(74)90015-2.
            
        This function is modified from a MATLAB function called co_co2flux(pco2_sea, pco2_air, sst, sss, wspd).
        Original MATLAB function author:
            Zelun Wu
            University of Delaware & Xiamen University
            zelunwu@outlook.com
            zelunwu.github.io
            
        Author:
        Abby Baskind (they/she)
        PhD Student
        URI Graduate School of Oceanography
        Wang Lab/Ocean Carbon Group
        abaskind@uri.edu
        GitHub: abby_baskind
        30 July 2023
        """
        import numpy as np
        sst = sst + 273.15 # transfer to Kelvin degree
        A1 =-58.0931
        A2 = 90.5069
        A3 = 22.294
        B1 = 0.027766
        B2 = -0.025888
        B3 = 0.0050578
        ln_K0 = A1 + A2 * (100/sst) + A3 * np.log(sst/100) + sss * (B1 + B2 * (sst/100) + B3 * (sst/100)**2)
        K0 = np.exp(ln_K0)
        return K0
    
    def gas_transfer_velocity(sst,wspd):
        """
        Calculate the gas transfer velocity from sst, sss, and wind speed
        Input:
            - sst: sea surface temperature
            - wspd: wind speed
        Defaults:
            - c: coefficient, default is 0.251 (Wanninkhof, 2014)
        Output:
            - kt: gas transfer velocity [cm/hr]
            
        References:
            - Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. 
            Limnology and Oceanography: Methods, 12(6), 351–362. https://doi.org/10.4319/lom.2014.12.351.
        
        This function is modified from a MATLAB function called co_co2flux(pco2_sea, pco2_air, sst, sss, wspd).
        Original MATLAB function author:
            Zelun Wu
            University of Delaware & Xiamen University
            zelunwu@outlook.com
            zelunwu.github.io
            
        Author:
        Abby Baskind (they/she)
        PhD Student
        URI Graduate School of Oceanography
        Wang Lab/Ocean Carbon Group
        abaskind@uri.edu
        GitHub: abby_baskind
        30 July 2023
        """
        A = 2116.8 
        B = -136.25
        C = 4.7353
        D = -0.092307
        E = 0.000755
        c = 0.251
        Sc = A + B*(sst) + C*(sst**2) + D*(sst**3) + E*(sst**4) # Jähne et al. (1987), Wanninkhof 2014
        kt = c * wspd**2 *((Sc/660)**(-0.5)) # unit: cm/hour
        return kt