# WREVAP

Operational Estimates of Areal Evapotranspiration and Lake Evaporation - Program WREVAP  
Python implementation of original Fortran model.

# Original Documentation

It is highly recommended that user read original WREVAP model documentation.  Most questions about WREVAP model, parameters, or units, can be answered by referring to original WREVAP documentation (Morton et al. 1985).  This version of WREVAP model is functionally identical to original WREVAP, however input format and options have been significantly improved for easy data entry and formatting.  All changes to input data, parameters, and format are outlined below.  Where possible, section number of original documentation is listed in square brackets [].

# File Setup

## CSV data file

Only file that is required to run WREVAP model is a CSV (comma separated value) file of time period starting date, temperature (T), humidity (TD), and insolation (S) (i.e. solar radiation) data.  waterborne energy input to lake (HADD) could also be included in this file, but it is not a required field, and later computed if not input.  For each entry, starting date must be specified either by listing YEAR, MONTH, DAY (of month) or YEAR and DOY (day of year).  It is also necessary to include LENGTH of time period.  This format is similar to Record E [6.5] of File Tape 1 (the .DAT file) in original WREVAP model where T, TD, and S are listed along with start date.

Script will only read data below field names (YEAR, LENGTH, T, TD, S) in CSV data file, so additional metadata can be stored in header of file.  In following example, site name is included in file, but it will not be read.  In original model, C1, C2, & C3 fields (Record D [6.4] and E [6.5] of File Tape 1) could be used to specify whether a T, TD, or S value were observed or estimated.  These values were not used by model and were removed from this version.  Finally, fields can be arranged in any order, but field names must be capitalized and spelled exactly.

## Units

Units for data are specified through script interface or a separate parameter INI file.


## Example Input

Following is an example CSV data file for Lahontan reservoir, Nevada.  Details on source of input data for Lahontan reservoir can be found in [Huntington and McEvoy (2011)](docs/Huntington_and_McEvoy_2011.pdf).

```
LAHONTAN,,,,,,
YEAR,MONTH,STARTDAY,LENGTH,TD,T,S
2000,1,1,31,-3.3,3.7,8.52
2000,2,1,29,-1.4,6.2,10.51
2000,3,1,31,-5.9,8.1,17.29
2000,4,1,30,-2.5,13.8,23.14
2000,5,1,31,-0.2,17.3,26.55
2000,6,1,30,1.4,23.8,28.84
2000,7,1,31,3.2,25.9,28.5
2000,8,1,31,14.3,24.7,22.8
2000,9,1,30,1.4,18.9,20.04
2000,10,1,31,-1.6,12.7,13.75
2000,11,1,30,-5.2,3,8.9
2000,12,1,31,-6.9,2.2,8.49
2001,1,1,31,-6.8,0,8.31
2001,2,1,28,-7.1,2.2,12.4
2001,3,1,31,-4.7,9.7,17.09
2001,4,1,30,-3.4,10.8,21.32
2001,5,1,31,0.1,20.8,27.21
2001,6,1,30,0.4,22.9,26.95
2001,7,1,31,4.6,26.3,26.38
2001,8,1,31,4.6,26.8,26.03
2001,9,1,30,3.2,22.2,20.8
2001,10,1,31,-1.1,15.8,14.57
2001,11,1,30,-1.4,7.4,9.28
2001,12,1,31,-2.5,2,5.92
```

## INI Parameter File

The model parameter values can be set using a parameter INI file or through script command prompt.  This parameter file is functionally equivalent to Records A [6.1] and B [6.2] of File Tape 1 (the .PAR file).  Record C is not included in this file since date and time period information is included in data CSV file.  script will read model parameter values from a file if it has same name as data file but with an “.INI” extension (for example, sample1.csv -> sample1.ini).  If an INI parameter file is not set, script will prompt user to enter necessary parameter values and then offer to save values into a new INI file.  This is easiest way to get a properly formatted (and commented) parameter file.  To build a parameter file from scratch, following values need to be specified in an INI file with same name as data file.  first non-commented line of INI file must be: [INPUTS].  All other values can be listed in any order after this first line.  Refer to original documentation for description, units, and suitable values for each of parameters.  # indicate commented lines for information only; comments could be removed or excluded to make input file significantly shorter in length.

```
# WREVAP INPUTS FILE
# FIRST DATA LINE MUST BE "[INPUTS]"
[INPUTS]

# SITE NAME
SITE = LAHONTAN

# LATITUDE [DECIMAL DEGREES]
PHID = 39.46

# STATION ALTITUDE OR PRESSURE (SEE PARAMETER IP)
# AVERAGE ATMOSPHERIC PRESSURE AT STATION [MB]
# ALTITUDE OF STATION ABOVE MEAN SEA LEVEL [M]
P = 1264.00

# AVERAGE ANNUAL PRECIPITATION [MM/YEAR]
# USED IF LK MODEL = 0
PPN = 0.0

# AVERAGE DEPTH OF LAKE [M]
# USED IF LK MODEL > 0
DA = 7.0

# TOTAL DISSOLVED SOLIDS OR SALINITY [PPM]
# USED IF LK MODEL > 0
SALT = 300.0

# LK - MODEL OPTION
#  0 - CRAE (AREAL EVAPOTRANSPIRATION)
#  1 - CRWE (WET SURFACE EVAPORATION)
#  2 - CRLE (LAKE EVAPORATION WITHOUT ANTECEDENT INFORMATION
#            ON SOLAR AND WATER BORNE ENERGY INPUTS)
#  3 - CRLE (LAKE EVAPORATION WITH ANTECEDENT SOLAR AND WATER
#            BORNE ENERGY INPUTS
LK = 3

# ISUM - CONTROL PARAMETER FOR STATION SUMMARY
#  0 - TABULATION OF AVERAGED MONTHLY TOTALS IS NOT LISTED (DEFAULT)
#  1 - TABULATION OF AVERAGED MONTHLY TOTALS IS LISTED
ISUM = 1

# IT - CONTROL PARAMETER FOR TEMPERATURE DATA
#  0 - AIR TEMPERATURE [DEGREES CELSIUS] (DEFAULT)
#  1 - AIR TEMPERATURE [DEGREES FAHRENHEIT]
IT = 0

# IS - CONTROL PARAMETER FOR INSOLATION DATA
#  0 - SUNSHINE DURATION RATIO
#  1 - SUNSHINE DURATION [HOURS/DAY] (DEFAULT)
#  2 - INCIDENT GLOBAL RADIATION [LY/DAY]
#  3 - INCIDENT GLOBAL RADIATION [MJ/M^2/DAY]
IS = 3

# IV - CONTROL PARAMETER FOR HUMIDITY DATA
#  0 - TD IS DEW POINT IN DEGREES [SEE PARAMETER IT FOR UNITS] (DEFAULT)
#  1 - TD IS VAPOUR PRESSURE AT DEW POINT [MB]
#  2 - TD IS RELATIVE HUMIDITY
IV = 0

# IP - CONTROL PARAMETER FOR STATION ALTITUDE
#  0 - AVERAGE ATMOSPHERIC PRESSURE AT STATION [MB] (DEFAULT)
#  1 - ALTITUDE OF STATION ABOVE MEAN SEA LEVEL [M]
IP = 1
```

## TGW/SOL available water and solar energy files

Format for these files is unchanged.  If no antecedent solar and waterborne energy inputs are available, CRLE with LK=2 must be run first.  .SOL file must be copied and renamed .TGW, and CRLE must be run with LK=3.  Refer to documentation on Record F [6.6] of File Tape 2 and Record G [6.7] of File Tape 3 for more information.

# Running Model

WREVAP script can be run by calling it directly from command prompt.  No GUI (graphical user interface) currently exist and all interaction with script is through command prompt and parameter and data files.

The script can be run by either setting **--data** or **--ini** and **--data** command line arguments:
```
> python wrevap\wrevap.py --data example\lahontan.csv --ini example\lahontan.ini

WREVAP - Python
  Input file: D:\WREVAP\example\lahontan.ini
  Data file:  D:\WREVAP\example\lahontan.csv

The following WREVAP parameter INI file was found: D:\WREVAP\example\lahontan.ini

Do you want to enter WREVAP parameters from this INI file [Y/n]?
Please enter choice: Y
```

If script is be run with only **--data** parameter set, in which case user will be prompted to enter input parameters (see below):
```
> python wrevap\wrevap.py --data example\lahontan.csv

WREVAP - Python
  Input file: None
  Data file:  D:\WREVAP\example\lahontan.csv

The WREVAP parameter INI file was not found
Please enter WREVAP parameters manually

Enter a site name: Lahontan
```

The example below illustrates user prompts for Lahontan reservoir, where CRLE model is desired given inputs of TD, T, and S, and no heat storage estimates are available.

```
------------------------------------------------------------
  PARAMETER LK - MODEL OPTION
------------------------------------------------------------
  0 - CRAE (AREAL EVAPOTRANSPIRATION)
  1 - CRWE (WET SURFACE EVAPORATION)
  2 - CRLE (LAKE EVAPORATION WITHOUT ANTECEDENT INFORMATION
            ON SOLAR AND WATER BORNE ENERGY INPUTS)
  3 - CRLE (LAKE EVAPORATION WITH ANTECEDENT SOLAR AND WATER
            BORNE ENERGY INPUTS)
------------------------------------------------------------
Please enter choice: 2
```

The example below illustrates user prompts for specifying units of input data and water body latitude, altitude (m), average depth (m) and salinity (mg/l or ppm).

```
------------------------------------------------------------
  IT - CONTROL PARAMETER FOR TEMPERATURE DATA
------------------------------------------------------------
  0 - AIR TEMPERATURE IN CELSIUS (default)
  1 - AIR TEMPERATURE IN FAHRENHEIT
------------------------------------------------------------
Please enter choice: 0

------------------------------------------------------------
  IS - CONTROL PARAMETER FOR INSOLATION DATA
------------------------------------------------------------
  0 - SUNSHINE DURATION RATIO
  1 - SUNSHINE DURATION IN HOURS/DAY (default)
  2 - INCIDENT GLOBAL RADIATION IN LY/DAY
  3 - INCIDENT GLOBAL RADIATION IN MJ/M**2/DAY
------------------------------------------------------------
Please enter choice: 3

------------------------------------------------------------
  IV - CONTROL PARAMETER FOR HUMIDITY DATA
------------------------------------------------------------
  0 - TD IS DEW POINT IN DEG.C (default)
  1 - TD IS VAPOUR PRESSURE AT DEW POINT IN M
  2 - TD IS RELATIVE HUMIDITY
------------------------------------------------------------
Please enter choice: 0

------------------------------------------------------------
  IP - CONTROL PARAMETER FOR STATION ALTITUDE
------------------------------------------------------------
  0 - AVERAGE ATMOSPHERIC PRESSURE AT STATION IN (default)
  1 - ALTITUDE OF STATION ABOVE MEAN SEA LEVEL IN M

------------------------------------------------------------
Please enter choice: 1

ENTER LATITUDE IN DECIMAL DEGREES [dd.dddd]: 39.46
ENTER STATION ALTITUDE [m]: 1264
ENTER AVERAGE DEPTH OF LAKE [m]: 7
ENTER TOTAL DISSOLVED SOLIDS OR SALINITY [mg/L or PPM]: 300
```

Results data are saved to a file with same name as data CSV file but with an RES extension.  overall format of results file is very similar to original output file, but exact spacing of values was changed.

Some error checking of input data and parameters is done but it is possible to enter inappropriate values.  Please refer to original documentation for details about suitable inputs and limitations of model.

# Reservoir Depths

Original version of WREVAP used average (constant) reservoir depth in meters.  This version of WREVAP was extended to allow average depth in feet and/or use provide reservoir depth as a time series.  Two additional configuration (INI) specifications were created to inform model of user’s reservoir depth usage.  Parameter RDU tells model what units of reservoir depth are (0 (default) is meters, 1 is feet).  Parameter RDM tells model what reservoir depth method is used.  RDM options are:

0 - average reservoir depth (default)
1 - time series reservoir depth embedded in main data file using header "RD".
2 - time series reservoir depth in a separate file and specified using configuration file.

To support a RDM value of 2, configuration variable "site_fn" was created to allow more flexibility in matching reservoir data filenames to WREVAP run.  If site_fn is not specified, it is assigned value of SITE_ID.  

Four formats are available for non-embedded reservoir depths.  Reservoir depths that are provided in a separate file use standard monthly date specifications.  Values are interpolated to match main data file’s timestamps.
```

SITE = Lahontan
site_fn = lahontan
RDM = 2
RDU = 1

# time series reservoir depth specifications

rd_folder = C:\pythonApps\clre\WREVAP-denver\example

# data structure type, file type (xls, csv, cdf, rdb) and file name format

# station file format

# rd_file_type = csv
# rd_ds_type = SF P
# rd_name_format = %s_sfp_rd.csv

# Parameter formats

# rd_ds_type = PF S.P

# individual parameter files example - %p is wildcard for parameter name (rd_fnspec)

# additional file specs - used if %p in name_format for parameter formats

rd_fnspec = rd

# csf setup

# rd_file_type = csf
# rd_name_format = lahontan_csf_rd.csv

# rdb setup

# rd_file_type = rdb
# rd_name_format = lahontan_rdb_rd.txt

# workbook setup

# rd_file_type = xls
# rd_name_format = lahontan_clre_data.xlsx

# worksheet names

rd_wsspec = RD

# used setup

rd_file_type = csv
rd_ds_type = SF P
rd_name_format = %s_sfp_rd.csv

# 1's based indices

rd_header_lines = 1
rd_names_line = 1
# rd_delimiter = \t
rd_delimiter = ,

# Field names

rd_date = Date
# rd_year = YEAR
# rd_month = MONTH
rd_field = Reservoir Depth
```

# Requirements

Python 2.7 must be installed on system.  Python 2.7 is installed by default with ArcGIS 10.1+, but if it needs to be installed, best source is [Anaconda](https://www.continuum.io/downloads) distribution and package manager. Be sure to select installer that corresponds with your operating system.  Code was not developed or tested for Python 3.

Python numerical modeling module [NumPy](http://www.numpy.org/) must also be installed.  NumPy is installed by default with ArcGIS 10.X or Anaconda Python distributions.

Following assumes that WREVAP script is being run on a computer running Microsoft Windows 7 or newer operating system.  Script should be able to execute on any computer that has Python and NumPy installed, but it was not tested on any older Windows or non-Windows machines.

# Reproducibility

Latest Python version of model may not always generate outputs that are identical to original Fortran model, however differences are negligible.  original Fortran version of model used some 32-bit floating point values and math functions, while Python version uses entirely 64 bit floating point values and functions, which are more precise than original version.  Rounding errors in original model result in negligible differences between output data.

# Conversion from Fortran

Upper case text and variables in scripts are copied directly from original Fortran code.

# References

PDFs for selected references are available in docs folder.

## WREVAP

[Morton, F.I., Ricard, F., and Fogarasi, S. (1985). Operational estimates of areal evapotranspiration and lake evaporation – Program WREVAP. NHRI Paper No. 24. Inland Waters Directorate. Ottawa, Canada.](docs/WREVAP_documentation_and_appendix.pdf)

[Huntington, J.L. and McEvoy, D. (2011). Climatological Estimates of Open Water Evaporation from Selected Truckee and Carson River Basin Water Bodies, California and Nevada. Desert Research Publication 41254, 34pp.](docs/Huntington_and_McEvoy_2011.pdf)

## CRLE

[Morton, F.I. (1979). Climatological estimates of lake evaporation. Water Resources Research, 15:64-76.](docs/Morton_1979.pdf)

[Morton, F.I. (1983a). Operational estimates of lake evaporation. Journal of Hydrology, 66:77-100.](docs/Morton_1983a.pdf)

[Morton, F.I. (1983b). Operational estimates of areal evapotranspiration and their significance to science and practice of hydrology. Journal of Hydrology, 66:1–76.](docs/Morton_1983b.pdf)

[Morton, F.I. (1986). Practical Estimates of Lake Evaporation. Journal of Climate and Applied Meteorology, 25(3):371-387.](docs/Morton_1986.pdf)

Morton, F.I. (1994). Evaporation research – A critical review and its lessons for environmental sciences. Critical Reviews in Environmental Science Technology, 24(3):237-280.

