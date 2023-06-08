class PLT():
    
    def get_hydrocat(start_date, end_date, buoy):
        """
        Retrieve Hydrocat data from Andy Davies' sensors using the API
        Drop 0 or negative pH values
        Change time zone from UTC to EST (-5 hours)
        Convert pH from NBS to total scale using PyCO2SYS
        
        INPUTS:
            - start_date: date as a string, formatted as YYYY-MM-DD
            - end_date: date as a string, formatted as YYYY-MM-DD
            - buoy: Jamestown/620/PLT or 720/Greenwich Bay/BG
        
        RETURNS pandas DataFrame
        
        Add '%cd directory' and '%run PLT.py' to your script, making sure the path is correct
        """
        import re
        import time
        import PyCO2SYS as pyco2
        import requests
        import pandas as pd
        import numpy as np
        from datetime import datetime, timedelta
        import math
    
        # Check that start and end dates are formatted correctly
        # Return error message if format is wrong
        matched_start = re.match("\d{4}-\d{2}-\d{2}", start_date)
        matched_end = re.match("\d{4}-\d{2}-\d{2}", end_date)
        if bool(not matched_start) or bool(not matched_end):
            print('Error: Input date must be a string of the form YYYY-MM-DD')
            return
    
        # Check that start date is before end date
        # Switch dates if start date after end date
        a = time.strptime(start_date, '%Y-%m-%d')
        b = time.strptime(end_date, '%Y-%m-%d')
        if a > b:
            start = start_date
            start_date = end_date
            end_date = start
    
        # create URL string using start date, end date, and buoy
        if buoy == 'Jamestown' or buoy == '620' or buoy == 'PLT': 
            URL = 'https://api.riddc.brown.edu/telemetry/Buoy-620/CoreMetrics/range?start=' + start_date + '&end=' + end_date
        elif buoy == 'Greenwich Bay' or buoy == '720' or buoy == 'GB': 
            URL = 'https://api.riddc.brown.edu/telemetry/Buoy-720/CoreMetrics/range?start=' + start_date + '&end=' + end_date
        else:
            print('Buoy options are 620 for Jamestown/PLT or 720 for Greenwich Bay/GB.')
            return 
    
        # Ingest data from API
        r = requests.get(URL)
        JSON = r.json()
        df = pd.DataFrame(JSON['Hydrocat'])
    
        # Drop 0 and negative values and nan values
        for ind in df.index:
            if df['hydrocatPH'][ind] <= 0:
                df = df.drop(ind)
            elif math.isnan(df['hydrocatPH'][ind]):
                df = df.drop(ind)
            
        # Get DateTime object from TmStamp string and change time zones from UTC            
        df['DateTime'] = np.zeros(len(df['TmStamp']))
        for ind in df.index:
            df['DateTime'][ind] = datetime.strptime(df['TmStamp'][ind], '%Y-%m-%dT%H:%M:%S.%fZ') - timedelta(hours = 5)
    
        # Convert pH from NBS to total scale using PyCO2SYS
        results = pyco2.sys(par1=df['hydrocatPH'], par1_type=3, 
                            temperature = df['hydrocatTemperature'], salinity = df['hydrocatSalinity'],opt_pH_scale = 4)
        df['pH total'] = results['pH_total']
    
        return df
    
    def dic_to_uM(dic, S, T):
        """
        This function converts from units of umol/kg to uM [umol/L] by calculating density using an equation of state.
        
        INPUTS:
            - DIC or TA in umol/kg
            - Salinity in PSU
            - Temperature in degrees C
        
        OUTPUTS:
            - DIC or TA in uM
        
        """
        rho =(999.842594 + 0.06793952*T 
              - 0.00909529*T**2 
              + 0.0001001685*T**3  
              - 0.000001120083*T**4 
              + 0.000000006536332*T**5 
              + (0.824493 - 0.0040899*T + 0.000076438*T**2 - 0.00000082467*T**3 + 0.0000000053875*T**4) * S 
              + (-0.00572466 + 0.00010227*T - 0.0000016546*T**2) * S**(1.5) + 0.00048314*S**2)/1000
        dic_out = rho * dic
        
        return dic_out