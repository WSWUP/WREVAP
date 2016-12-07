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

# WREVAP - GCM

Original WREVAP code was modified to support being run using Global Climite Model (GCM) input data.

Differences between versions:
+ An input parameter file must be passed to the script (using the **--ini** command line argument)
+ A separate CSV output file (*.OUT) is created
+ Precipitation (PPT) values are carried through the model

December 2016 modifications to support dynamic reservoir depths:
+ Same as WREVAP
