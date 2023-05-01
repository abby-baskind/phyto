# phyto
Narragansett Bay Long Term Phytoplankton Time Series pH analysis

![A test image](C:/Users/akbaskind/Desktop/YoE/IMG_6436 2.JPG)

<p align="center">
  <img src="/Users/akbaskind/Desktop/YoE/IMG_6436 2.JPG" width="350" title="hover text">
</p>
 
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
