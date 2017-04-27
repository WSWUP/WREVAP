# WREVAP

Output format in RES file was adjusted from original (Fortran)

December 2016 modifications to support dynamic reservoir depths:
+ Added RDU parameter to specify units (0 is Meters; 1 is Feet)
+ Added RDM parameter to specify reservoir depth type and input type.
  0 is average depth (default)
  1 is embedded in main time series file using attribute 'RD'
  2 is provided in a separate file specified using configuration (ini) file
+ Available formats for RDM of 2 are:
  SFP - Station file with station's parameters
  CSF - Column Slot File (text version of a worksheet)
  RDB - Text relational database
  XLS - Workbook

March 2017 modifications to support additional post processing:
+ Added ".OUT" csv output to wrevap as per wrevap_gcm.

April 2017 modifications to address issue number 1:
+ Added 'DECIMAL' before 'RELATIVE HUMIDITY' in code, example ini file and documentation.

# WREVAP - GCM

Original WREVAP code was modified to support being run using Global Climite Model (GCM) input data.

Differences between versions:
+ An input parameter file must be passed to the script (using the **--ini** command line argument)
+ A separate CSV output file (*.OUT) is created
+ Precipitation (PPT) values are carried through the model

December 2016 modifications to support dynamic reservoir depths:
+ Same as WREVAP

# run_wrevap

March 2017 Created run_wrevap to to support single trace multiple reservoirs runs:
+ Created run_wrevap for managing single trace runs for multiple stations.
  Uses a list of reservoirs and their respective ini and data files.
