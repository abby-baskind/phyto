class PLT():
    
    def get_hydrocat(start_date, end_date, buoy):
        """
        Retrieve Hydrocat data from Andy Davies' sensors using the API
        Drop 0 or negative pH values
        Time zone kept in UTC
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
            elif math.isnan(df['hydrocatDissOxygen'][ind]):
                df = df.drop(ind)
            elif math.isnan(df['hydrocatSalinity'][ind]):
                df = df.drop(ind)
            elif math.isnan(df['hydrocatTemperature'][ind]):
                df = df.drop(ind)
        df = df.reset_index(drop=True)
            
        # Get DateTime object from TmStamp string 
        # Keep timezone in UTC
        df['DateTime'] = np.zeros(len(df['TmStamp']))
        for ind in df.index:
            df['DateTime'][ind] = datetime.strptime(df['TmStamp'][ind], '%Y-%m-%dT%H:%M:%S.%fZ') 
    
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
    def get_lab_samples(googlesheet_url, wks_name):
        
        """
        Retrieve lab data from Wang Lab Google folder using the API
        Convert data to appropriate data types
        Drop rows where TA, DIC, and/or salinity is not recorded
        Time zone converted to UTC
        Convert lab measured TA [uM] and DIC [uM] to TA [umol/kg] and DIC [umol/kg]...
            ...using gsw.rho(S,T,P)
        Calculate total pH using PyCO2SYS
        
        INPUTS:
            - URL of the Google sheet as a string
            - worksheet name as a string
        
        RETURNS pandas DataFrame with
            - Sample Name
            - Sample Location
            - DateTime in UTC
            - Sample Location
            - Sample depth
            - DIC [umol/kg]
            - TA [umol/kg]
            - Salinity
            - In Situ Temperature [degC]
            - pH in total scale
        
        Add '%cd directory' and '%run PLT.py' to your script, making sure the path is correct
        """
    
        # # PACKAGES REQUIRED --------------------------------------------------------------------------------------
        import gspread
        from df2gspread import df2gspread as d2g
        from oauth2client.service_account import ServiceAccountCredentials
        import gsw
        import PyCO2SYS as pyco2
        import pandas as pd
        import math
        from datetime import datetime, timedelta
        # # --------------------------------------------------------------------------------------------------------
    
        # # API SETUP ----------------------------------------------------------------------------------------------
        scope = ['https://spreadsheets.google.com/feeds',
                 'https://www.googleapis.com/auth/drive']

        # Name of our Service Account Key
        google_key_file = 'servicecredentials.json'
        credentials = ServiceAccountCredentials.from_json_keyfile_name(google_key_file, scope)
        gc = gspread.authorize(credentials)
        # # --------------------------------------------------------------------------------------------------------
    
    
        # # PARSE URL STRING TO GET THE SPREADSHEET KEY-------------------------------------------------------------
        # strip url of google spreadsheet
        spreadsheet_key = googlesheet_url.strip('https://docs.google.com/spreadsheets/d/')
        # split remainder of URL at the /
        spreadsheet_key = spreadsheet_key.split('/')
        # take the split with the spreadsheet key
        spreadsheet_key = spreadsheet_key[0]
        # # --------------------------------------------------------------------------------------------------------
    
        # # PULL DATA-----------------------------------------------------------------------------------------------
        #Opening the worksheet by using Worksheet ID
        workbook = gc.open_by_key(spreadsheet_key)
        #Selecting which sheet to pulling the data
        sheet = workbook.worksheet(wks_name)

        #Pulling the data and transform it to the data frame
        values = sheet.get_all_values()
        labdata = pd.DataFrame(values[1:], columns = values[0])
        # # --------------------------------------------------------------------------------------------------------
    
        # # DATA CLEANUP -------------------------------------------------------------------------------------------
        # Since I have notes written in some of the columns, 
        # we will first select which columns we would like from the dataframe. 
        # I've also found that often when pulling from Google Sheets, all the data is written as strings, 
        # so here we also convert our data to appropriate data types.
        df = labdata[['Sample', 'Time', 'Location', 'depth', 'TA Temp (degC)', 'TA (uM)',
                      'DIC Temp (degC)', 'DIC (uM)', 'Salinity', 'In Situ Temperature']]
        df['TA (uM)'] = pd.to_numeric(df['TA (uM)'])
        df['TA Temp (degC)'] = pd.to_numeric(df['TA Temp (degC)'])
        df['DIC (uM)'] = pd.to_numeric(df['DIC (uM)'])
        df['DIC Temp (degC)'] = pd.to_numeric(df['DIC Temp (degC)'])
        df['Salinity'] = pd.to_numeric(df['Salinity'])
        df['In Situ Temperature'] = pd.to_numeric(df['In Situ Temperature'])
        df["DateTime"] = pd.to_datetime(df["Time"], errors = 'coerce')
    
        # # Drop missing data
        for ind in df.index:
            # if TA or DIC or salinity or in situ tempertuare at that index is Nan/empty...
            if math.isnan(df['TA (uM)'][ind]) or math.isnan(df['DIC (uM)'][ind]) or math.isnan(df['Salinity'][ind]) or math.isnan(df['In Situ Temperature'][ind]):                                                                      
                # Drop that index
                df = df.drop(ind)
        # # --------------------------------------------------------------------------------------------------------
            
        # # CONVERT TIME ZONE TO UTC--------------------------------------------------------------------------------
        # Eastern time daylights savings is 4 hours behind UTC
        # Eastern standard time is 5 hours behind UTC
        daylight_start_date23 = pd.Timestamp(datetime(2023, 3, 12, 0, 0,0))     # Begin 2023 Daylight Time
        daylight_end_date23 = pd.Timestamp(datetime(2023, 11, 5, 0, 0,0))       # End 2023 Daylight Time
        daylight_start_date22 = pd.Timestamp(datetime(2022, 3, 13, 0, 0,0))     # Begin 2022 Daylight Time
        daylight_end_date22 = pd.Timestamp(datetime(2022, 11, 6, 0, 0,0))       # End 2022 Daylight Time

        # For each index in dataframe df...
        for ind in df.index:
            if (df['DateTime'][ind] >= daylight_start_date23) & (df['DateTime'][ind] < daylight_end_date23):
                df['DateTime'][ind] = df['DateTime'][ind] + timedelta(hours = 4)
            elif (df['DateTime'][ind] >= daylight_start_date22) & (df['DateTime'][ind] < daylight_end_date22):
                df['DateTime'][ind] = df['DateTime'][ind] + timedelta(hours = 4)
            else:
                df['DateTime'][ind] = df['DateTime'][ind] + timedelta(hours = 5)
        # # --------------------------------------------------------------------------------------------------------
            
        # # CONVERT FROM uM TO umol/kg------------------------------------------------------------------------------
        # Pressure
        df['P'] = 0

        # Density
        df['rho_DIC'] = gsw.rho(df['Salinity'], df['DIC Temp (degC)'], df['P'])
        df['rho_TA'] = gsw.rho(df['Salinity'], df['TA Temp (degC)'], df['P'])

        # Converted DIC and TA
        df['DIC (umol/kg)'] = df['DIC (uM)']/df['rho_DIC']/0.001
        df['TA (umol/kg)'] = df['TA (uM)']/df['rho_TA']/0.001
    
        # Calculate pH with PyCO2SYS
        results = pyco2.sys(par1=df['TA (umol/kg)'],par2=df['DIC (umol/kg)'],par1_type=1,par2_type=2,
                            salinity = df['Salinity'], temperature = df['In Situ Temperature'])
        df['pH'] = results['pH']
        # # --------------------------------------------------------------------------------------------------------
        dfout = df[['Sample', 'DateTime', 'Location', 'depth', 'Salinity', 'In Situ Temperature', 'DIC (umol/kg)', 'TA (umol/kg)', 'pH']]
        
        return dfout
    
    def get_NBFSMN(googlesheet_url, wks_name):
        
        """
        Retrieve NBFSMN data from Wang Lab Google folder using the API
        Convert data to appropriate data types
        Drop rows where TA, DIC, and/or DO is not recorded
        Drop columns with unnecessary data
        Convert NBS pH to total pH
        
        INPUTS:
            - URL of the Google sheet as a string
            - worksheet name as a string
        
        RETURNS pandas DataFrame with
            - Sample Date
            - DateTime in UTC
            - Surface Temperature [degC]
            - Surface Salinity
            - Surface DO concentration [mg/L]
            - pH in total scale
        
        Add '%cd directory' and '%run PLT.py' to your script, making sure the path is correct
        """
    
        # # PACKAGES REQUIRED --------------------------------------------------------------------------------------
        import gspread
        from df2gspread import df2gspread as d2g
        from oauth2client.service_account import ServiceAccountCredentials
        import gsw
        import PyCO2SYS as pyco2
        import pandas as pd
        import math
        from datetime import datetime, timedelta
        # # --------------------------------------------------------------------------------------------------------
    
        # # API SETUP ----------------------------------------------------------------------------------------------
        scope = ['https://spreadsheets.google.com/feeds',
                 'https://www.googleapis.com/auth/drive']

        # Name of our Service Account Key
        google_key_file = 'servicecredentials.json'
        credentials = ServiceAccountCredentials.from_json_keyfile_name(google_key_file, scope)
        gc = gspread.authorize(credentials)
        # # --------------------------------------------------------------------------------------------------------
    
    
        # # PARSE URL STRING TO GET THE SPREADSHEET KEY-------------------------------------------------------------
        # strip url of google spreadsheet
        spreadsheet_key = googlesheet_url.strip('https://docs.google.com/spreadsheets/d/')
        # split remainder of URL at the /
        spreadsheet_key = spreadsheet_key.split('/')
        # take the split with the spreadsheet key
        spreadsheet_key = spreadsheet_key[0]
        # # --------------------------------------------------------------------------------------------------------
    
        # # PULL DATA-----------------------------------------------------------------------------------------------
        #Opening the worksheet by using Worksheet ID
        workbook = gc.open_by_key(spreadsheet_key)
        #Selecting which sheet to pulling the data
        sheet = workbook.worksheet(wks_name)

        #Pulling the data and transform it to the data frame
        values = sheet.get_all_values()
        MVdata = pd.DataFrame(values[3:], columns = values[1])
        # # --------------------------------------------------------------------------------------------------------
    
        # # DATA CLEANUP -------------------------------------------------------------------------------------------
        # NSFSMN sites record surface and bottom data. 
        # For now, I am just selecting the surface data in columns 0 through 15, 
        # rename some of the columns, and drop some of the unnecessary columns
        surfaceMV = MVdata.iloc[:, 0:15]
        surfaceMV = surfaceMV.drop(columns=['Site', "Agency ", 'surface FS'])
        surfaceMV = surfaceMV.rename(columns={"C": "DateTime", 'Date ': 'Date'})
        
        # Correct data types
        for key in surfaceMV.keys():
            if key == 'DateTime' or key == 'Date' or key == 'Time':
                surfaceMV[key] = pd.to_datetime(surfaceMV[key])
            else:
                surfaceMV[key] = pd.to_numeric(surfaceMV[key])
    
        # # Drop missing data
        for ind in surfaceMV.index:
            if math.isnan(surfaceMV['surface pH'][ind]) or math.isnan(surfaceMV['surface DO Conc'][ind]) or math.isnan(surfaceMV['surface Temp'][ind]) or math.isnan(surfaceMV['surface Salinity'][ind]):
                surfaceMV = surfaceMV.drop(ind)
        surfaceMV = surfaceMV.reset_index(drop=True)
        # # --------------------------------------------------------------------------------------------------------
            
        # # CONVERT FROM NBS PH TO TOTAL PH-------------------------------------------------------------------------
        resultsMV = pyco2.sys(par1=surfaceMV['surface pH'], par1_type=3, temperature = surfaceMV['surface Temp'], 
                              salinity = surfaceMV['surface Salinity'],opt_pH_scale = 4)
        surfaceMV['pH total'] = resultsMV['pH_total']
        # # --------------------------------------------------------------------------------------------------------
        df = surfaceMV[['Date', 'DateTime', 'surface Temp',
                        'surface Salinity', 'surface DO Conc',
                        'pH total']]
        
        return df
                                                                                                                    
    def get_buoy(start_date, end_date, buoy, sensor):
        """
        Retrieve data from Andy Davies' sensors using the API
        Choose from PAR, SUNA, or MetData sensors for, respectively,
            - light
            - nutrients
            - meteorological data
        Drop 0 or negative pH values
        Time zone kept in UTC
        
        INPUTS:
            - start_date: date as a string, formatted as YYYY-MM-DD
            - end_date: date as a string, formatted as YYYY-MM-DD
            - buoy: Jamestown/620/PLT or 720/Greenwich Bay/BG
            - PAR, SUNA, or MetData sensors
        
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
        df = pd.DataFrame(JSON[sensor])
    
        # Drop 0 and negative values and nan values
        for k in df.keys():
            for ind in df.index:
                if isinstance(df[k][ind], float) or isinstance(df[k][ind], int):
                    if math.isnan(df[k][ind]):
                        df = df.drop(ind)
                        df = df.reset_index(drop=True)
            
        # Get DateTime object from TmStamp string 
        # Keep timezone in UTC
        df['DateTime'] = np.zeros(len(df['TmStamp']))
        for ind in df.index:
            df['DateTime'][ind] = datetime.strptime(df['TmStamp'][ind], '%Y-%m-%dT%H:%M:%S.%fZ') 
    
        return df
    
    def o2sat(S,T):
        """
        CALCULATE OXYGEN CONCENTRATION AT SATURATION        
        adapted from o2sat.m by: Edward T Peltzer, MBARI (revised 2007 Apr 26)

        Source: The solubility of nitrogen, oxygen and argon in water and seawater - Weiss (1970) Deep Sea Research V17(4): 721-735.

        Molar volume of oxygen at STP obtained from NIST website on the thermophysical properties of fluid systems:
        http://webbook.nist.gov/chemistry/fluid/

        Inputs: 
        S = Salinity (PSU)
        T = Temperature (degC)

        Outputs:
            Oxygen saturation at one atmosphere (umol/kg)
        """
        import numpy as np
        T1 = (T + 273.15) / 100

        OSAT = -177.7888 + 255.5907 / T1 + 146.4813 * np.log(T1) - 22.2040 * T1
        OSAT = OSAT + S * (-0.037362 + T1 * (0.016504 - 0.0020564 * T1))
        OSAT = np.exp(OSAT)

        #     Convert from ml/kg to um/kg
        O2 = OSAT * 1000 / 22.392
        return O2