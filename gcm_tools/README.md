# common

Usage:
  Copy prep_gcm_wrevap_gcm, calc_gcm_wrevap_gcm, plot_gcm_wrevap_gcm and wrevap_gcm_gcm.py
    to common folder near study area's data files.

# prep_gcm_wrevap

December 2016 modifications:
+ Miscellaneous code cleanup.
+ Renamed some conversion variables.
+ Created variable 'output_monthly_prcp'
+ Added support of optional dynamic reservoir depths.
+ Extended solar constant from 4.92 to 4.9212 to be consistent with ET applications.

Usage:
  Modify code for input and output data folders (it does not use current working directory)
    with absolute paths or relative to code folder (assuming code is run there).
  Edit file "reservoir_parameters.csv" that lists following attributes:
    GCM_NAME,OUTPUT_NAME,DEPTH,SALINITY,ELEVATION,LATITUDE,LONGITUDE
    Check Thornton-Running for units.
    Lat/Lon are in decimal degrees with lon being negative for western hemisphere.
    Thornton-Running B calculation parameters named as TR_B0, TR_B1, TR_B2,
    Mean monthly dew point depression (ko) named as KO_JAN, KO_FEB, etc.
  Optional reservoir depth variables RDM (method) and RDU (units) can also exist in reservoir_paramters.csv.
  Reservoirs can be specified using a list in code or obtained automatically from data folders.
  Edit "default.ini" file for wrevap_gcm, parameters will be updated by script to 
    reservoir specific contents from parameters in "reservoir_parameters.csv".
  Run script which reads data for all reservoirs and gcm traces, estimates solar radiation
      and creates wrevap_gcm configuration (ini) and time series data file(s).

# calc_gcm_wrevap.py

December 2016 modifications:
+ Miscellaneous code cleanup.
+ Modified system calls and folder deletions.
+ Changed 'WREVAP' to 'wrevap' to support unix operating system.
+ Added support of optional dynamic reservoir depths.

Usage
  Modify code for input data folder (output of prep_gcm_wrevap is input of calc_gcm_wrevap).
  Run script which runs wrevap_gcm for each data file (gcm trace) in input_data folder.
  wrevap_gcm is run first with LK = 2, SOL file is copied to TGW, then wrevap_gcm is run with LK = 3.
  wrevap_gcm posts to an ensemble specific ".OUT" file which is a modified version of wrevap's RES output.
  ".OUT" file does not have mean monthly values or header row and is formatted as a CSV file.

# plot_gcm_wrevap.py

Usage:
  Modify code for input data folder (output of calc_gcm_wrevap_gcm is input of plot_gcm_wrevap_gcm).
  Run script.
