
# phyto
Narragansett Bay Long Term Phytoplankton Time Series pH analysis

![The bay!](bay.JPG)

Most recent analysis can be viewed in `QAQC_TimeSeries.ipynb`.
![Time Series](prelim_pH_literallyeverything.png)
![Time Series](QP_Buoy620_PLT.png)
![Time Series](MV_Buoy720_GB.png)
![Time Series](lab_minus_davies_anoms.png)


# Data

* Hydrocat sensors from Davies Lab
* Narragansett Bay Fixed Station Monitoring Network - MV
* Wang Lab SeaFET sensor
* NB Long Term Phytoplankton Time Series Bottle Samples, analyzed by Wang Lab

# Data Retrieval
I've opted to set up [API retrieval for Google sheets](https://towardsdatascience.com/from-google-sheet-to-your-jupyter-notebook-ccdbf28fbf1b). NBFSM data available [here](https://docs.google.com/spreadsheets/d/1m61xkdLMaMSvw533FmVDIBqJqF73QF6b/edit#gid=924440352). Bottle samples availabke [here](https://docs.google.com/spreadsheets/d/17FFbtUuhUS4UtxB-OjKIP2wCYJoEAmaW6VaHQPcup9U/edit#gid=0). Up-to-date SeaFET data is not yet available via Google Sheets.

## `PLT` class

The `PLT` class has a function that ingests Hydrocat sensor data from the Davies lab (`get_hydrocat(start_date, end_date, buoy)`).

# Required Dependencies

* `PyCO2SYS`
* `re`
* `time`
* `pandas`
* `numpy`
* `xarray`
* `requests`
* `datetime`
* `math`
* `gspread` 
* `oauth2client` 
* `df2gspread`

# Effect of rain on pH

I merged data for Hydrocat 620 (PLT) and Hydrocat 720 (GB), both of which had been subjected to the QA/QC process, with weather data from the same buoy. Weather data included maximum precipitation for every 15 minutes measurement interval. I differentiated pH over time for finite time differences, as the data is only every 15 minutes. I then filtered the data to only include days in the 99th percentile of precipitation; in other words, I selected only the rainiest days. I also selected data from the time period 24 hours after the rain event. The analysis for this can be seen in Precipitation.ipynb.

First, I plotted the discrete change in pH at the time of precipitation over time for both Hydrocat 620 and Hydrocat 720. The colormap corresponds to the amount of precipitation. The mean change in pH during the rain event, with standard error, is effectively 0 for both sensors, suggesting rain has no impact on pH.
![Hydrocat 620](H620_fromrain.png)
![Hydrocat 720](H720_fromrain.png)

Then, I plotted the change in pH over the 24 hours following the rain event over time for both Hydrocat 620 and Hydrocat 720. The colormap corresponds to the amount of precipitation. Similarly, the mean change in pH during the rain event, with standard error, is effectively 0 for both sensors, suggesting rain has no delayed impact on pH, at least within 24 hours. Indeed, I tried this with other time intervals, namely 12 and 48 hours, and similarly no significant change in pH manifested.
![Hydrocat 620](H620_24hr_afterrain.png)
![Hydrocat 720](H720_24hr_afterrain.png)

We can also get a perspective of pH following rain looking directly at the time series. The full time series of pH for each sensor is in red. The measurement intervals with the most precipitation are in blue. Intervals 24 hours after the rainiest intervals are in green. From this coarse perspective, we do not consistently see either increased on decreased pH measurements following rain. 
![Hydrocat 620](H620_rain.png)
![Hydrocat 720](H720_rain.png)

# Diurnal cycle

To detect the diurnal cycle of pH in the Narragansett Bay, we used Prophet, statistical software designed for time series analysis. Documentation is available [here](https://facebook.github.io/prophet/docs/quick_start.html#python-api). The associated publication for this software is [Taylor, S. J. & Letham, B., 2018. Forecasting at scale, The American Statistician, 72(1), 37–45.](https://www.tandfonline.com/doi/full/10.1080/00031305.2017.1380080?casa_token=tnhGKYTQIIYAAAAA%3AvaObcD-ZD_S8Ld-XhxryUePRNsDLlcaSCpIlhPycdRN_HLddkcSFGL00UJW-0SnpVlgJj9BToQxC) Using Prophet, we are able to tease out various timescales of variability.

We see a consistent diurnal cycle across all timeseries. pH decreases in the night, reaching its local minimum sometime after midnight local time prior to dawn. pH then increases to its maximum around midday. The diurnal cycle of pH closely follows that of oxygen, suggesting daily variability is closely tied to biological activity. Indeed, this makes sense, as during the night, respiration will outpace photosynthesis and in turns produces CO2 in the water column and lowers pH. Conversely, photosynthesis increases during daylight hours, reaching its maximum around the time of maximum sunlight. Photosynthesis consumes CO2 and consequently increases pH.
![Diurnal Cycle](diurnalcycle_ph.png)
![Diurnal Cycle](diurnalcycle_ph_o2.png)

# Seasonal cycle
With just over a year of data, we have *just* enough data to begin examining the seasonal cycle of pH in Narragansett Bay. This assumes the ~18 months of data are more or less typical in their variability. The seasonaly cycle that manifests agrees with that observed by Pimenta et al. (2023). pH is highest in the late winter and early spring, approximately the timing of the wintertime phytoplankton bloom. By summer, pH decreases to its annual minimum, after the collapse of the bloom and an increase in respiration relative to photosynthesis. The decrease in photosynthesis could be related to several things: an exhaustion of the bay's nutrient supply, too much sunlight, or simply a popultion boom of non-photosynthetic species (e.g. zooplankton), to name a few possibilities. The autumn brings slightly increased pH, as a (weaker) fall phytoplankton bloom commences. This seasonal cycle appears to be driven strongly by biological activity (among other things) and the bay's unique seasonality with its wintertime bloom. As the climate continues to change, this seasonality may change as the winter bloom weakens and temporally shifts. 
![Monthly Cycle](monthly_ph.png)
![Seasonal Cycle](seasonal_ph.png)

# Drivers of pH Variability

We have attempted to separate the various drivers of pH in Narragansett Bay into several components: temperature, salinity, mixing of alkalinity, biological production of alkalinity, biological production of DIC, mixing of DIC, and air-sea fluxes of DIC. We converted pH to H+ concentration to avoid issues with logarithms, and modelled the components as follows, based on Kwiatokowski and Orr:

The figure below shows the results of this separation of drivers.

![pH Drivers](ALL_components.png)