# phyto
Narragansett Bay Long Term Phytoplankton Time Series pH analysis

![A test image](IMG_6436.JPG)

 
# `PLT` class

The `PLT` class has a function that ingests Hydrocat sensor data from the Davies lab (`get_hydrocat(start_date, end_date, buoy)`).

# Data Retrieval
I've opted to set up [API retrieval for Google sheets](https://towardsdatascience.com/from-google-sheet-to-your-jupyter-notebook-ccdbf28fbf1b).

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
