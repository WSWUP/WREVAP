# WREVAP

Operational Estimates of Areal Evapotranspiration and Lake Evaporation - Program WREVAP  
Python implementation of the original Fortran model.

# Original Documentation

It is highly recommended that the user read the original WREVAP model documentation.  Most questions about the WREVAP model, parameters, or units, can be answered by referring to the original WREVAP documentation (Morton et al. 1985).  This version of the WREVAP model is functionally identical to the original WREVAP, however input format and options have been significantly improved for easy data entry and formatting.  All changes to the input data, parameters, and format are outlined below.  Where possible, the section number of the original documentation is listed in square brackets [].

# File Setup

## CSV data file

The only file/data required to run the WREVAP model is a CSV (comma separated value) file of the time period starting date, temperature (T), humidity (TD), and insolation (S) (i.e. solar radiation) data.  The waterborne energy input to the lake (HADD) could also be included in this file, but it is not a required field, and later computed if not input.  For each entry, the starting date must be specified either by listing the YEAR, MONTH, DAY (of month) or the YEAR and DOY (day of year).  It is also necessary to include the LENGTH of the time period.  This format is similar to Record E [6.5] of File Tape 1 (the .DAT file) in the original WREVAP model where the T, TD, and S are listed along with the start date.

The script will only read data below the field names (YEAR, LENGTH, T, TD, S) in the CSV data file, so additional metadata can be stored in the header of the file.  In the following example, the site name is included in the file, but it will not be read.  In the original model, the C1, C2, & C3 fields (Record D [6.4] and E [6.5] of File Tape 1) could be used to specify whether a T, TD, or S value were observed or estimated.  These values were not used by the model and were removed from this version.  Finally, the fields can be arranged in any order, but the field names must be capitalized and spelled exactly.

## Units

The units for the data are specified through the script interface or a separate parameter INI file.


## Example Input

The following is a sample of the example CSV data file for Lahontan reservoir, Nevada.  Details on the source of the input data for Lahontan reservoir can be found in [Huntington and McEvoy (2011)](docs/Huntington_and_McEvoy_2011.pdf).

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

The model parameter values can be set using a parameter INI file or through the script command prompt.  This parameter file is functionally equivalent to Records A [6.1] and B [6.2] of File Tape 1 (the .PAR file).  Record C is not included in this file since the date and time period information is included in the data CSV file.  The script will read the model parameter values from a file if it has the same name as the data file but with an “.INI” extension (for example, sample1.csv -> sample1.ini).  If an INI parameter file is not set, the script will prompt the user to enter the necessary parameter values and then offer to save the values into a new INI file.  This is the easiest way to get a properly formatted (and commented) parameter file.  To build a parameter file from scratch, the following values need to be specified in an INI file with the same name as the data file.  The first non-commented line of the INI file must be: [INPUTS].  All other values can be listed in any order after this first line.  Refer to the original documentation for the description, units, and suitable values for each of the parameters.  # indicate commented lines for information only; comments could be removed or excluded to make the input file significantly shorter in length.

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

# AVERAGE DEPTH OF THE LAKE [M]
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

The format for these files is unchanged.  If no antecedent solar and waterborne energy inputs are available, the CRLE with LK=2 must be run first.  The .SOL file must be copied and renamed .TGW, and the CRLE must be run with LK=3.  Refer to the documentation on Record F [6.6] of File Tape 2 and Record G [6.7] of File Tape 3 for more information.

# Running The Model

The WREVAP script can be run by calling the script directly from the command prompt.  Currently there is no GUI (graphical user interface) and all interaction with the script is the through the command prompt and the parameter and data files.

The script can be run be either setting both the **--ini** and **--data** command line arguments:
```
> python wrevap\wrevap.py --data example\lahontan.csv --ini example\lahontan.ini

WREVAP - Python
  Input file: D:\WREVAP\example\lahontan.ini
  Data file:  D:\WREVAP\example\lahontan.csv

The following WREVAP parameter INI file was found: D:\WREVAP\example\lahontan.ini

Do you want to enter WREVAP parameters from this INI file [Y/n]?
Please enter choice: Y
```

The script can also be run with only the **--data** parameter set, in which case user will be prompted to enter the input parameters (see below):
```
> python wrevap\wrevap.py --data example\lahontan.csv

WREVAP - Python
  Input file: None
  Data file:  D:\WREVAP\example\lahontan.csv

The WREVAP parameter INI file was not found
Please enter WREVAP parameters manually

Enter a site name: Lahontan
```

The example below illustrates the user prompts for Lahontan reservoir, where the CRLE model is desired given inputs of TD, T, and S, and no heat storage estimates are available.

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

The example below illustrates the user prompts for specifying units of input data and water body latitude, altitude (m), average depth (m) and salinity (mg/l or ppm).

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

The results data is saved to a file with the same name as the data CSV file but with an RES extension.  The overall format of the results file is very similar to the original output file, but the exact spacing of the values was changed.

There is some error checking of the input data and parameters, but it is possible to enter inappropriate values.  Please refer to the original documentation for details about suitable inputs and the limitations of the model.

# Requirements

Python 2.7 must be installed on the system.  Python 2.7 is installed by default with ArcGIS 10.1+, but if it needs to be installed, the best source is the [Anaconda](https://www.continuum.io/downloads) distribution and package manager.  Make sure to select the installer that corresponds with your operating system.  The code was not developed or tested for Python 3.

The Python numerical modeling module [NumPy](http://www.numpy.org/) must also be installed.  NumPy is installed by default with ArcGIS 10.X or Anaconda Python distributions.

All of the following notes are assuming that the WREVAP script is being run on a computer running the Microsoft Windows 7 operating system or newer.  The script should be able to execute on any computer that has Python and NumPy installed, but it was not tested on any older Windows or non-Windows machines.

# Reproducibility

The new Python version of the model may not always generate outputs that are identical to the original Fortran model, however the differences are negligible.  The original Fortran version of the model used some 32-bit floating point values and math functions, while the Python version uses entirely 64 bit floating point values and functions, which are more precise than the original version.  Rounding errors in the original model result in negligible differences between the output data.

# Conversion from Fortran

Upper case text and variables in the scripts are copied directly from the original Fortran code.

# References

PDFs for selected references are available in the docs folder.

## WREVAP

[Morton, F.I., Ricard, F., and Fogarasi, S. (1985). Operational estimates of areal evapotranspiration and lake evaporation – Program WREVAP. NHRI Paper No. 24. Inland Waters Directorate. Ottawa, Canada.](docs/WREVAP_documentation_and_appendix.pdf)

[Huntington, J.L. and McEvoy, D. (2011). Climatological Estimates of Open Water Evaporation from Selected Truckee and Carson River Basin Water Bodies, California and Nevada. Desert Research Publication 41254, 34pp.](docs/Huntington_and_McEvoy_2011.pdf)

## CRLE

[Morton, F.I. (1979). Climatological estimates of lake evaporation. Water Resources Research, 15:64-76.](docs/Morton_1979.pdf)

[Morton, F.I. (1983a). Operational estimates of lake evaporation. Journal of Hydrology, 66:77-100.](docs/Morton_1983a.pdf)

[Morton, F.I. (1983b). Operational estimates of areal evapotranspiration and their significance to the science and practice of hydrology. Journal of Hydrology, 66:1–76.](docs/Morton_1983b.pdf)

[Morton, F.I. (1986). Practical Estimates of Lake Evaporation. Journal of Climate and Applied Meteorology, 25(3):371-387.](docs/Morton_1986.pdf)

Morton, F.I. (1994). Evaporation research – A critical review and its lessons for the environmental sciences. Critical Reviews in Environmental Science Technology, 24(3):237-280.

