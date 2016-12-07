#--------------------------------
# Name:         wrevap_gcm.py
# Authors:      Charles Morton and Justin Huntington
# Modified:     2016-11-02
# Python:       2.7
#--------------------------------

import argparse
import calendar
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import math
import os
import sys

import numpy as np
import pandas as pd

class paths():
    pass
class param():
    pass
class const():
    pass
class data():
    pass

def WREVAP(input_path, data_path):
    """Operational Estimates of Areal Evapotranspiration and Lake Evaporation

    Program WREVAP
    Python implementation of original Fortran model
    Modified to run on large numbers of GCM datasets

    Args:
        input_path (str): file path to input parameter file
        data_path (str): file pat to CSV data file
    """
    logging.info('WREVAP GCM - Python')
    logging.info('  Input file: {}'.format(input_path))
    logging.info('  Data file:  {}'.format(data_path))

    # Initialize file paths

    set_input_paths(input_path, data_path)

    # logging.info('  INI Filename:  {}'.format(paths.ini))
    # logging.info('  Data Filename: {}'.format(paths.csv))

    # Initialize Parameters

    get_parameters()

    # Read input time series data

    read_input_tsdata()

    # CHECK INPUT SPECIFICATIONS AND DO A GENERAL INITIALIZATION

    initialize()

    # Process data for each time period

    for i, dt_start in enumerate(data.DATE):
        dt_stop = dt_start + dt.timedelta(days=(data.LENGTH[i] - 1))
        doy_start = dt_start.timetuple().tm_yday
        doy_stop = dt_stop.timetuple().tm_yday
        doy_adjust = adjust_doy(dt_start)

        # COMPUTE RADIUS VECTOR OF EARTH'S ORBIT(neta) AND
        # DECLINATION OF SUN(theta)
        neta, theta = compute_orbit(doy_start, doy_stop, doy_adjust)

        # CHECK RAW DATA AND IF NECESSARY CONVERT UNITS AND THEN,
        # COMPUTE ABSORBED GLOBAL RADIATION OR
        # SOLAR & WATER-BORNE HEAT INPUT(GW).

        compute_budget_1(i, neta, theta)

        if param.LK in [0, 1]:
            # DO FINAL CALCULATION AND LIST INPUTS AND
            # RESULTS FOR EACH PERIOD.
            compute_budget_2(i)
        elif param.LK in [2, 3]:
            # IF LAKE EVAPORATION IS BEING ESTIMATED, RESULTS ARE
            # STORED FOR LATER USE.
            
            TRANLK(i)

    # FOLLOWING SECTION DEALS WITH LAKE EVAPORATION MODEL
    
    if param.LK in [2, 3]:
        # READ 12 SOLAR & WATER-BORNE HEAT INPUT VALUES(TGW) OF
        # PRECEDING TIME PERIODS, OR COPY MOST RECENT 12 VALUES OF
        # TGW INTO PRECEDING PERIODS.
        if param.LK == 3:
            read_tgw_file()
        if param.LK == 2:
            for i in xrange(12):
                data.TGW[i] = data.TGW[i + 12]

        # CALCULATE AVAILABLE SOLAR & WATER-BORNE HEAT
        # CALLING IT TABLE TGL.
        
        compute_available_heat()

        # DO FINAL CALCULATION AND LIST INPUTS AND RESULTS FOR EACH PERIOD.
        
        for i in xrange(param.NN):
            compute_budget_2(i)

        # WRITE 12 MOST RECENT VALUES OF TGW TO SOL FILE
        write_sol_file()
    print_output()

    # LIST MONTHLY TOTALS AVERAGED OVER 5 YEARS

    if param.ISUM == 1:
        print_monthly_averages()

def set_input_paths(input_path, data_path):
    """"""
    # For now, get workspace and name from data_path
    input_ws = os.path.dirname(data_path)
    input_name, input_ext = os.path.splitext(os.path.basename(data_path))
    paths.input_ws = input_ws

    paths.ini = input_path
    paths.csv = data_path
    paths.res = os.path.join(input_ws, input_name + '.RES')
    paths.tgw = os.path.join(input_ws, input_name + '.TGW')
    paths.sol = os.path.join(input_ws, input_name + '.SOL')
    
    # Write separate output file for GCM runs
    
    paths.out = os.path.join(input_ws, input_name + '.OUT')

def read_param(p_str, p_default, p_type, config):
    """"""
    try:
        if p_type is float:
            p_value = float(config.get('INPUTS', p_str))
        elif p_type is int:
            p_value = int(config.get('INPUTS', p_str))
        else:
            # String type or NoneType
            # p_value = config.get('INPUTS', p_str).upper()    # upper is problematic for file specs
            p_value = config.get('INPUTS', p_str)
    except:
        p_value = p_default
    return p_value

def check_param_in_list(p_value, p_list, p_name):
    """"""
    if p_value not in p_list:
        logging.error(
            '\nERROR: {} paramter must be set to {}'.format(
                p_name, ', '.join(['{}'.format(i) for i in p_list])))
        sys.exit()

def get_parameters():
    """"""
    # Parameter lists
    LK_list = [0, 1, 2, 3]
    ISUM_list = [0, 1]
    IP_list = [0, 1]
    IS_list = [0, 1, 2, 3]
    IT_list = [0, 1]
    IV_list = [0, 1, 2]
    RDM_list = [0, 1, 2]
    RDU_list = [0, 1]

    # Parameters must be read from file

    if not os.path.isfile(paths.ini):
        logging.error(
            '\nERROR: Input file {} does not exist\n'.format(paths.ini))
        sys.exit()

    # Read parameters from file and check for basic errors
    
    logging.info('\nReading CONFIGURATION file: {}'.format(
        os.path.basename(paths.ini)))

    # Check that INI file opens and has section entry ([INPUTS])
    # Get list of all parameter keys

    config = ConfigParser.SafeConfigParser()
    try:
        config.readfp(open(paths.ini))
        config.has_section('INPUTS')
        config_items = config.items('INPUTS')
    except ConfigParser.NoSectionError:
        logging.error((
            '\nERROR: {}\n' +
            '    Paramter INI file is missing a section line\n' +
            '    First data line in file needs to be: [INPUTS]\n' +
            '    Try removing INI file and rebuilding it\n').format(
            paths.ini))
        sys.exit()
    except ConfigParser.Error:
        logging.error((
            '\nERROR: {}\n' +
            '    There is an unknown problem with paramter INI file \n' +
            '    Try removing INI file and rebuilding it\n').format(
            paths.ini))
        sys.exit()

    # Check that all parameters are present in INI file

    param_list = ['SITE', 'PHID', 'P', 'LK', 'ISUM', 'IT', 'IS', 'IV', 'IP']
    config_keys = [i.upper() for i in dict(config_items).keys()]
    for param_str in param_list:
        if param_str not in config_keys:
            logging.error((
                '\nERROR: Parameter {} was not found in ' +
                'parameter INI file').format(param_str))
            sys.exit()
    if (('PPN' not in config_keys) or
        ('DA' not in config_keys and
         'SALT' not in config_keys)):
            logging.error(
                '\nERROR: Either PPN parameter or DA and SALT ' +
                'parameters must be specified')
            sys.exit()

    # Read parameters from file

    SITE = read_param('SITE', '', str, config)
    site_fn = read_param('site_fn', SITE, str, config)
    PHID = read_param('PHID', None, float, config)
    if not (-90 <= PHID <= 90):
        logging.error(
            '\nERROR: PHID paramter must be between -90 and +90')
        sys.exit()
    P = read_param('P', None, float, config)
    if not P or P < 0:
        logging.error("\nERROR: P paramter must be >= 0")
        sys.exit()
    LK = read_param('LK', None, float, config)
    check_param_in_list(LK, LK_list, 'LK')
    if LK == 0:
        PPN = read_param('PPN', None, float, config)
        if not PPN or PPN < 0:
            logging.error('\nERROR: PPN paramter must be >= 0 when LK==0')
            sys.exit()
        RDM = 0
        RDU = 0
        DA = 0.0
        SALT = 0.0
    elif LK in [1, 2, 3]:
        PPN = 0.0
        RDM = read_param('RDM', 0, int, config)
        check_param_in_list(RDM, RDM_list, 'RDM')
        RDU = read_param('RDU', 0, int, config)
        check_param_in_list(RDU, RDU_list, 'RDU')
        if RDM == 2:
            # read non embedded reservoir depth specifications
                
            param.rd_folder = read_param('rd_folder', paths.input_ws, str, config)
            param.rd_file_type = read_param('rd_file_type', 'csv', str, config).lower()
            param.rd_ds_type = read_param('rd_ds_type', 'SF P', str, config).upper()
            param.rd_name_format = read_param('rd_name_format', '%s_rd.csv', str, config)
            param.rd_header_lines = read_param('rd_header_lines', 2, int, config)
            param.rd_names_line = read_param('rd_names_line', 2, int, config)
            param.rd_delimiter = read_param('rd_delimiter', ',', str, config).lower()
            if param.rd_delimiter not in [' ', ',', '\\t']: param.rd_delimiter = ','
            if "\\" in param.rd_delimiter and "t" in param.rd_delimiter:
                param.rd_delimiter = param.rd_delimiter.replace('\\t', '\t')
            if RDU == 1:
                param.rd_units = 'Feet'
            else:
                param.rd_units = 'Meter'
                
            # Date can be read directly or computed from year, month, and day        

            param.rd_date = read_param('rd_date', None, str, config)
            param.rd_year = read_param('rd_year', None, str, config)
            param.rd_month = read_param('rd_month', None, str, config)
            param.rd_day = read_param('rd_day', None, str, config)
            if param.rd_date is not None:
                logging.info('  INMET: Reading date from date column')
            elif (param.rd_year is not None and
                  param.rd_month is not None and
                  param.rd_day is not None):
                logging.info('  INMET: Reading date from year, month, and day columns')
            else:
                logging.error('  ERROR: RESERVOIR DEPTH date_field (or year, month, and '+
                              'day fields) must be set in  INI')
                sys.exit()                  
            param.rd_field = read_param('rd_field', 'RD', str, config)
            param.rd_fnspec = read_param('rd_name', None, str, config)
            if param.rd_fnspec is None or param.rd_fnspec == 'None': 
                param.rd_fnspec = 'RD'
            if param.rd_file_type == 'xls' or param.rd_file_type == 'wb':
                param.rd_wsspec = read_param('rd_wsspec', None, str, config)
                if param.rd_wsspec is None or param.rd_wsspec == 'None':
                    logging.info('  INFO:  INMET: reservoir depth worksheet name set to RD')
                    param.rd_wsspec = 'RD'
            else:
                param.rd_wsspec = 'RD'
        else:
            param.rd_folder = None
            param.rd_file_type = None
            param.rd_ds_type = None
            param.rd_name_format = None
            param.rd_header_lines = None
            param.rd_names_line = None
            param.rd_delimiter = None
            param.rd_date = None
            param.rd_year = None
            param.rd_month = None
            param.rd_day = None
            param.rd_field = None
            param.rd_fnspec = None
            param.rd_wsspec = None
        DA = read_param('DA', None, float, config)
        if RDU == 1:
            # convert static avarge depth to meters
              
            DA /= 3.2808399
        SALT = read_param('SALT', None, float, config)
        if not DA or DA <= 0:
            logging.error(
                '\nERROR: DA paramter must be > 0 when LK > 0')
            sys.exit()
        elif not SALT or SALT < 0:
            logging.error(
                '\nERROR: SALT paramter must be >= 0 when LK > 0')
            sys.exit()
    ISUM = read_param('ISUM', 0, int, config)
    IP = read_param('IP', 0, int, config)
    IS = read_param('IS', 1, int, config)
    IT = read_param('IT', 0, int, config)
    IV = read_param('IV', 0, int, config)
    check_param_in_list(ISUM, ISUM_list, 'ISUM')
    check_param_in_list(IP, IP_list, 'IP')
    check_param_in_list(IS, IS_list, 'IS')
    check_param_in_list(IT, IT_list, 'IT')
    check_param_in_list(IV, IV_list, 'IV')

    # Save parameters

    param.SITE = SITE
    param.site_fn = site_fn
    param.PHID = PHID
    param.P = P
    param.PW = P       # Save original P value
    param.PPN = PPN
    param.RDM = RDM
    param.RDU = RDU
    param.DA = DA
    param.SALT = SALT
    param.LK = LK
    param.ISUM = ISUM
    param.IT = IT
    param.IS = IS
    param.IV = IV
    param.IP = IP

def read_input_tsdata():
    """PROCESS ALL INPUT DATA"""
    logging.info('\nReading DATA file: {}'.format(
        os.path.basename(paths.csv)))
    dat_f = open(paths.csv, 'rb')
    dat_lines = [line.strip() for line in dat_f.readlines()]
    dat_f.close()

    # First line may be SITE name or Fields

    test_line = dat_lines[0].split(',')
    for test_i, test_line in enumerate(dat_lines):
        if ('YEAR' in test_line and 'LENGTH' in test_line and
                'S' in test_line and 'T' in test_line and 'TD' in test_line):
            header_line = test_i
            data_line = test_i + 1
            logging.debug(
                '  Assuming line {} is field names\n  Fields: {}'.format(
                    (header_line + 1),
                    ' '.join(dat_lines[header_line].split(','))))
    dat_header = dict([
        (item.upper(), i) for i, item in enumerate(
            dat_lines[header_line].split(','))])
    dat_fields = dat_header.keys()

    # Build data arrays based on number of time steps

    param.NN = len(dat_lines) - data_line
    data.DATE = [''] * param.NN
    data.LENGTH = [0] * param.NN
    data.TD = [0.0] * param.NN
    data.T = [0.0] * param.NN
    data.S = [0.0] * param.NN
    data.TDW = [0.0] * param.NN
    data.TW = [0.0] * param.NN
    data.SW = [0.0] * param.NN
    data.HADD = [0.0] * param.NN
    data.PPT = [0.0] * param.NN
    data.reservoir_depth = [param.DA] * param.NN
    logging.info('  {} data points in file\n'.format(param.NN))

    # Get column numbers for each field

    def column_index(name_list):
        for name in name_list:
            if name in dat_fields:
                return dat_header[name]
        return None
    year_i = column_index(['YEAR'])
    month_i = column_index(['MONTH'])
    day_i = column_index(['DAY', 'STARTDAY', 'START_DAY'])
    doy_i = column_index(['DOY', 'STARTDOY', 'START_DOY'])
    length_i = column_index(['LENGTH'])
    td_i = column_index(['TD'])
    t_i = column_index(['T'])
    s_i = column_index(['S'])
    hadd_i = column_index(['HADD'])
    ppt_i = column_index(['PPT'])
    if param.RDM == 1:
        rd_i = column_index(['RD'])
    else:
        rd_i = -1

    # Check date fields

    if year_i >= 0 and doy_i >= 0 and length_i >= 0:
        dt_doy_flag = True
        dt_format = '%Y_%j'
    elif year_i >= 0 and day_i >= 0 and day_i >= 0 and length_i >= 0:
        dt_doy_flag = False
        dt_format = '%Y_%m_%d'
    else:
        logging.error('  Fields: {}'.format(dat_fields))
        logging.error(
            '\nERROR: Data file is missing start date or length fields' +
            '\nERROR: YEAR/STARTDOY/LENGTH or YEAR/MONTH/STARTDAY/LENGTH\n')
        sys.exit()

    # Check data fields

    if not (s_i and t_i and td_i):
        logging.error(
            '\nERROR: Data file is missing data fields' +
            '\nERROR: TD, T, or S\n')
        sys.exit()

    # Read data

    dat_list = [i.split(',') for i in dat_lines[data_line:]]

    # Process data for each time period

    for i, dat_line in enumerate(dat_list):
        # logging.debug('\n  {}'.format(dat_line))
        year = int(dat_line[year_i])
        if dt_doy_flag:
            dt_start = dt.datetime.strptime(
                '{}_{:03d}'.format(year, int(dat_line[doy_i])),
                dt_format)
        else:
            dt_start = dt.datetime.strptime(
                '{}_{:02d}_{:02d}'.format(
                    year, int(dat_line[month_i]), int(dat_line[day_i])),
                dt_format)
        data.DATE[i] = dt_start
        data.LENGTH[i] = int(dat_line[length_i])
        data.TD[i] = float(dat_line[td_i])
        data.T[i] = float(dat_line[t_i])
        data.S[i] = float(dat_line[s_i])
        # data.HADD[i] = float(dat_line[hadd_i]) if hadd_i >= 0 else 0.0
        if hadd_i >= 0:
            data.HADD[i] = float(dat_line[hadd_i])
        if param.RDM == 1:
            # reservoir depth is embedded in data file
                
            if rd_i >= 0:
                if param.RDU == 0:
                    data.reservoir_depth[i] = float(dat_line[rd_i])
                else:
                    data.reservoir_depth[i] = float(dat_line[rd_i]) / 3.2808399
            else:
                logging.error('\nERROR: Unable to find reservoir depth column')
                sys.exit()
        # data.PPT[i] = float(dat_line[ppt_i]) if ppt_i >= 0 else 0.0
        if ppt_i >= 0:
            data.PPT[i] = float(dat_line[ppt_i])
    param.start_dt = dt.datetime(data.DATE[0].year, data.DATE[0].month, 1, 0, 0, 0)
    param.end_dt = data.DATE[len(data.DATE) - 1]
    if param.RDM == 2:
        # reservoir depth is specified in configuration file

        if not read_reservoir_depths():
            logging.error('\nERROR: Unable to process reservoir depths')
            sys.exit()

def initialize():
    """
    CHECKS INPUT SPECIFICATIONS FOR, LOGICAL ERROR
    MODIFIES THESE SPECIFICATIONS IF NECESSARY,
    THEN INITIALIZE CONSTANTS AND TABLES
    """
    P = param.P
    if param.IP == 1:
        P = 1013 * (1 - 0.0065 * P / 288) ** 5.256
    if P < 0:
        logging.error(
            '\n*** INPUT ERROR NUMBER = 23\n' +
            'Value of P should be greater than or equal to 0\n')
        sys.exit()
    param.P = P

    # DIFFERENT CONSTANTS ARE USED DEPENDING ON WHETHER AREAL
    # EVAPOTRANSPIRATION MODEL IS BEING USED(LK=0) OR WET SURFACE
    # EVAPORATION OR LAKE EVAPORATION MODELS ARE BEING USED(LK>0)
    if param.LK == 0:
        # CONSTANTS FOR AREAL EVAPOTRANSPIRATION MODELS

        const.SB = 5.22 / 1E8
        const.AZZZ = (0.26 - 0.00012 * math.sqrt(P / 1013) * param.PPN *
                      (1 + abs(param.PHID) / 42 + (param.PHID / 42) ** 2))
        const.CONST2 = 1.20
        const.FZ = 28.0
        const.CONST1 = 14.0
    else:
        # CONSTANTS FOR WET SURFACE AND LAKE EVAPORATION MODELS

        const.SB = 5.5 / 1E8
        const.AZZZ = 0.05
        const.CONST2 = 1.12
        const.FZ = 25.0
        const.CONST1 = 13.0

    # VARIABLES FOR MONTHLY MEANS ARE INITIALIZED
    # IF ISUM=0 , NO MONTHLY SUMMARY IS PRINTED AT END
    # IF ISUM=1 , TABLE OF MONTHLY MEANS IS PRINTED

    if param.ISUM == 1:
        data.TETM = defaultdict(dict)
        data.TETPM = defaultdict(dict)
        data.TRTM = defaultdict(dict)
        for DATE in data.DATE:
            data.TETM[DATE.month][DATE.year] = []
            data.TETPM[DATE.month][DATE.year] = []
            data.TRTM[DATE.month][DATE.year] = []
    data.V = [0.0] * param.NN
    data.VD = [0.0] * param.NN
    data.GW = [0.0] * param.NN
    data.ETMM = [0.0] * param.NN
    data.RTMM = [0.0] * param.NN
    data.ETPMM = [0.0] * param.NN
    if param.LK >= 2:
        data.TDATE = [0] * param.NN
        data.TLENGTH = [0] * param.NN
        data.TJ = [0] * param.NN
        data.TT = [0.0] * param.NN
        data.TSW = [0.0] * param.NN
        data.TTDW = [0.0] * param.NN
        data.TTW = [0.0] * param.NN
        data.TS = [0.0] * param.NN
        data.THADD = [0.0] * param.NN
        data.TV = [0.0] * param.NN
        data.TVD = [0.0] * param.NN
        # For Lake Evap model TGW stores data for 12 previous time periods
        data.TGW = [0.0] * (param.NN + 12)
        data.TGL = [0.0] * param.NN
        data.GWT = [0.0] * param.NN

    # DEFINE TEMPERATURE DEPENDENT CONSTANTS
    # 0 - FOR TEMPERATURE GREATER THAN OR EQUAL TO 0 DEGREES CELSIUS.
    # 1 - FOR TEMPERATURE LESS THAN 0 DEGREES CELSIUS

    const.ALPHA = [17.27, 21.88]
    const.BETA = [237.3, 265.5]
    const.GAMMA = [(0.66 * P / 1013), ((0.66 * P / 1013) / 1.15)]
    const.FTZ = [(const.FZ * math.sqrt(1013 / P)),
                 ((const.FZ * math.sqrt(1013 / P)) * 1.15)]

    # DEFINE OTHER VALUES

    const.PI = math.pi
    const.CONV = math.pi / 180
    param.PHI = param.PHID * const.CONV

def read_tgw_file():
    """"""
    # if os.path.isfile(paths.sol) and not os.path.isfile(paths.tgw):
    #     logging.info('  NOTE: Copying data in SOL file to TGW file')
    #     shutil.copy2(paths.sol, paths.tgw)
    try:
        tgw_f = open(paths.tgw, 'r')
    except:
        logging.error(
            ('ERROR: TGW file does not exist or ' +
             'could not be opened\n    ({})').format(paths.tgw))
        logging.error(
            '  If running LK==3, try running LK==2 first\n' +
            '  Then copy data in SOL file to TGW file\n' +
            '  See section 6.6 and 6.7 in NHRI Paper No. 24')
        sys.exit()
    tgw_data = [float(i) for i in tgw_f.readlines()]
    tgw_f.closed
    if len(tgw_data) != 13:
        logging.error('ERROR: There should be 13 values in .TGW file')
        sys.exit()
    param.GLBGN = tgw_data.pop(0)
    # Values are reversed in file?
    # data.TGW[:12] = tgw_data[:]
    data.TGW[:12] = tgw_data[::-1]
    del tgw_data, tgw_f

def write_sol_file():
    """"""
    with open(paths.sol, 'w') as sol_f:
        sol_f.write('{0:9.4f}\n'.format(param.GLEND))
        # Values are reversed in file? (reverse sort then get last 12)
        sol_f.write('\n'.join(map('{0:9.4f}'.format, data.TGW[::-1][:12])))
        sol_f.closed

def adjust_doy(dt_start):
    """JULIAN DAY ADJUSTMENT FOR LEAP YEARS"""
    if calendar.isleap(dt_start.year):
        doy_adjust = -0.5
    else:
        doy_adjust = 0.5
    if dt_start.month <= 2:
        doy_adjust = 0.0
    # logging.debug('  DOY Adjust: {}'.format(doy_adjust))
    return doy_adjust

def compute_orbit(doy_start, doy_stop, doy_adjust):
    """
    CALCULATES ORBITAL PARAMETERS SUCH AS:
        DECLINATION OF SUN IN RADIAN(theta)
        NON-DIMENSIONNAL RADIUS VECTOR OF EARTH'S ORBIT
        AROUND SUN(NETA)
    """
    doy_array = (np.arange(doy_start, doy_stop + 1) + doy_adjust)
    orb1 = np.minimum((29.5 + doy_array / 270.0), 30.4)
    orb2 = ((doy_array + 0.5 * (orb1 - 1)) / orb1)
    neta = 1 + np.sin((29.5 * orb2 - 106.0) * const.CONV) / 60.0
    theta = (23.45 * const.CONV * np.sin((29.5 * orb2 - 94.0) * const.CONV))
    neta = np.mean(neta)
    theta = np.mean(theta)
    return neta, theta

def compute_budget_1(i, neta, theta):
    """
    CHECKS LOGICAL ERRORS IN RAW DATA INPUT
    CONVERTS RAW DATA  INTO STANDARD UNITS.
    CALCULATES V,VD,GW
    ASSIGNS VALUES TO  J,TDW,TW,SW
    """
    TD = data.TD[i]
    T = data.T[i]
    S = data.S[i]
    HADD = data.HADD[i]

    # STORE RAW DATA INPUT FOR PRINT OUT

    data.TDW[i] = TD
    data.TW[i] = T
    data.SW[i] = S

    # CHECK LOGICAL ERRORS OF RAW DATA INPUT

    if param.IV >= 1 and TD < 0:
        logging.error(
            '\n*** INPUT ERROR NUMBER = 5\n' +
            '  TD values should be > 0 when TD represents ' +
            'vapour pressure at dew point or relative humidity\n')
        sys.exit()
    if S < 0:
        logging.error(
            '\n*** INPUT ERROR NUMBER = 6\n  S values should be >= 0\n')
        sys.exit()

    # CONVERT FAHRENHEIT TEMPERATURE INTO CELSIUS
    
    if param.IT != 0:
        if param.IV == 0:
            TD = (TD - 32) * 5. / 9.
        T = (T - 32) * 5. / 9.
    J = 1 if T < 0 else 0

    # SATURATION VAPOUR PRESSURE(V) AT AIR TEMPERATURE
    
    V = 6.11 * math.exp(const.ALPHA[J] * T / (T + const.BETA[J]))

    # SATURATION VAPOUR PRESSURE(VD) AT DEW POINT TEMPERATURE
    # WHEN TD IS DEW POINT
    
    if param.IV == 0:
        VD = 6.11 * math.exp(17.27 * TD / (TD + 237.3))
    elif param.IV == 1:
        VD = TD
    else:
        VD = TD * V    # WHEN TD IS RELATIVE HUMIDITY

    # CALCULATE SUNSHINE DURATION RATIO, (S), WHEN
    # SUNSHINE DURATION IS GIVEN IN HOURS PER DAY
    
    if param.IS == 1:
        CZENAA = math.cos(param.PHI - theta) + 0.005
        if CZENAA < 0.001:
            CZENAA = 0.001
        ACOMA = 1 - CZENAA / math.cos(param.PHI) / math.cos(theta)
        if ACOMA < -1:
            ACOMA = -1.0
        OMEGAA = math.acos(ACOMA)
        MSD = 24 * OMEGAA / const.PI
        S = S / MSD

    # SOLAR ZENITH ANGLE (ZENA)
    
    CZENA = math.cos(param.PHI - theta)
    if CZENA < 0.001:
        CZENA = 0.001
    ZENA = math.acos(CZENA)

    # NO. OF DEGREES EARTH ROTATES BETWEEN SUNRISE AND NOON (OMEGA)
    
    ACOM = 1 - (CZENA / math.cos(param.PHI) / math.cos(theta))
    if ACOM < -1:
        ACOM = -1.0
    OMEGA = math.acos(ACOM)

    # COSINE OF AVERAGE ANGULAR ZENITH DISTANCE OF SUN (COSZ)
    
    COSZ = (CZENA + (math.sin(OMEGA) / OMEGA - 1) *
            math.cos(param.PHI) * math.cos(theta))

    # EXTRA-ATMOSPHERIC GLOBAL RADIATION(GE)
    
    GE = 1354.0 * COSZ * OMEGA / (const.PI * neta ** 2)

    # SNOW FREE, CLEAR SKY ALBEDOIF SUN WERE AT  ZENITH (AZZ)
    
    AZZ = const.AZZZ

    # CONSTRAINT FOR AZZ WHICH MAY APPLY DURING WET SEASON IN DRY AREA
    
    if param.LK == 0:
        if (AZZ > ((0.91 - VD / V) / 2)):
            AZZ = (0.91 - VD / V) / 2
            
        # TWO FURTHER CONSTRAINTS FOR AZZ
        
        if AZZ > 0.17:
            AZZ = 0.17
        elif AZZ < 0.11:
            AZZ = 0.11

    # WEIGHTING FACTOR FOR EFFECT OF SNOW ON ALBEDO (ARAT)
    
    VPDL = V - VD
    if VPDL < 0:
        VPDL = 0.0
    elif VPDL > 1:
        VPDL = 1.0
    ARAT = 1 - VPDL * VPDL

    # CLEAR SKY ALBEDO IF SUN WERE AT ZENITH(AZZ)

    AZ = AZZ + ARAT * (0.34 - AZZ)

    # CLEAR SKY ALBEDO (A0)

    A0 = (AZ * (math.exp(1.08) - math.exp(ZENA * 2.16 / const.PI) *
                (math.cos(ZENA) * 2.16 / const.PI + math.sin(ZENA))) /
          (1.473 * (1.0 - math.sin(ZENA))))

    # PRECIPITABLE WATER VAPOUR (W)

    W = VD / (0.49 + T / 129)

    # WEIGHTING FACTOR FOR EFFECTS OF TEMPERATURE ON TURBIDITY (TST)
    
    TST = 21 - T
    if TST < 0:
        TST = 0.0
    elif TST > 5:
        TST = 5.0

    # TURBIDITY COEFFICIENT (DUST)

    DUST = (0.5 + 2.5 * COSZ ** 2) * math.exp(TST * (param.P / 1013 - 1))

    # TRANSMITTANCY OF CLEAR SKIES TO DIRECT BEAM SOLAR RADIATION(TAUT)

    DUSTT = 0.083 * (DUST / COSZ) ** 0.9
    WVT = 0.029 * (W / COSZ) ** 0.6
    LNX = (-0.089 * (param.P / (COSZ * 1013)) ** 0.75 - DUSTT - WVT)
    if LNX < -675:
        LNX = -675.0
    TAUT = math.exp(LNX)

    # PART OF TAUT THAT IS RESULT OF ABSORBTION(TAUA)

    WVA = math.sqrt(WVT / 10.0)
    if WVA > WVT:
        WVA = WVT
    LNY = -DUSTT / 2 - WVA
    if LNY < -675:
        LNY = -675.0
    TAUA = math.exp(LNY)

    # CLEAR SKY GLOBAL RADIATION (G0)

    G0 = GE * (TAUT + TAUT * (1 - TAUT / TAUA) * (1 + A0 * TAUT))

    # ESTIMATE INCIDENT GLOBAL RADIATION FROM GIVEN VALUE OF
    # SUNSHINE DURATION RATIO (S)
    
    if param.IS in [0, 1]:
        G = G0 * S + (0.08 + 0.3 * S) * (1 - S) * GE
    # OR IF IS = 2,THEN S IS INCIDENT GLOBAL RADIATION(LANGLEY/DAY)
    # OR IF IS = 3,THEN S IS INCIDENT GLOBAL RADIATION(MJ/M**2/DAY)
    elif param.IS in [2, 3]:
        if param.IS == 2:
            G = S / 2.064
        elif param.IS == 3:
            G = S / 0.0864
        S = 0.53 * G / (G0 - 0.47 * G)
        if S > 1:
            S = 1.0
        elif S < 0:
            S = 0.0
    # ESTIMATE AVERAGE ALBEDO (A)
    
    A = A0 * (S + (1 - ZENA * 180 / const.PI / 330) * (1 - S))

    # COMPUTE GW WHICH IS ABSORBED GLOBAL RADIATION FOR AREAL
    # EVAPOTRANSPIRATION AND WET SURFACE EVAPORATION OPTIONS, AND
    # SOLAR & WATER-BORNE HEAT INPUT FOR LAKE EVAPORATION OPTION.
    
    GW = (1 - A) * G
    if param.LK > 1:
        GW += HADD

    # Save parameters

    data.TD[i] = TD
    data.T[i] = T
    data.S[i] = S
    data.V[i] = V
    data.VD[i] = VD
    data.GW[i] = GW

def TRANLK(i):
    """
    THIS SUBROUTINE IS USED ONLY FOR LAKE EVAPORATION MODEL.
    ALL THESE TRANLK ARGUMENTS ARE STORE IN TABLES FOR LATER
    APPLICATION IN BUDGT2.  TABULATED NAMES ARE PREFIXED WITH LETTER
    'T' AND FOLLOWED BY PREVIOUS NON-SUBSCRIPTED TRANLK ARGUMENTS.
    """
    data.TDATE[i] = data.DATE[i]
    data.TLENGTH[i] = data.LENGTH[i]
    data.TT[i] = data.T[i]
    data.TS[i] = data.S[i]
    data.TTDW[i] = data.TDW[i]
    data.TTW[i] = data.TW[i]
    data.TSW[i] = data.SW[i]
    data.TV[i] = data.V[i]
    data.TVD[i] = data.VD[i]
    data.THADD[i] = data.HADD[i]
    
    # Lake Evap model TGW stores data for 12 previous time periods
    
    data.TGW[i + 12] = data.GW[i]
    return

def compute_budget_2(i):
    """
    BUDGT2 CALCULATES THE:
        POWER EQUIVALENT OF AREAL EVAPOTRANSPIRATION, LAKE-SIZE
        WET SURFACE EVAPORATION OR LAKE EVAPORATION(ET),
        POWER EQUIVALENT OF POTENTIAL EVAPOTRANSPIRATION,
        PAN-SIZE WET SURFACE EVAPORATION OR POTENTIAL EVAPORATION(ETP),
        NET RADIATION OR NET AVALAIBLE ENERGY WITH SOIL-PLANT
        SURFACES, WET SURFACE OR LAKE SURFACE AT TEMPERATURE T(RT)
        AND CALLS print_output TO LIST INPUTS AND RESULTS FOR EACH TIME PERIOD
    """
    if param.LK <= 1:
        DATE = data.DATE[i]
        LENGTH = data.LENGTH[i]
        T = data.T[i]
        VD = data.VD[i]
        V = data.V[i]
        S = data.S[i]
        GW = data.GW[i]
        GL = data.GW[i]
    elif param.LK >= 2:
        DATE = data.TDATE[i]
        LENGTH = data.TLENGTH[i]
        T = data.TT[i]
        VD = data.TVD[i]
        V = data.TV[i]
        S = data.TS[i]
        GW = data.TGW[i + 12] + 1E-11
        GL = data.TGL[i]
    J = 1 if T < 0 else 0

    # WEIGHTING FACTOR FOR EFFECT OF CLOUDS ON ATMOSPHERIC RADIATION (ATM)

    ATM = 10 * (VD / V - S - 0.42)
    if ATM < 0:
        ATM = 0.0
    elif ATM > 1:
        ATM = 1.0

    # PROPORTIONAL INCREASE IN ATMOSPHERIC RADIATION DUE TO CLOUDS (RHO)

    RHO = (0.18 * 1013 / param.P *
           (ATM * math.sqrt(1 - S) + (1 - ATM) * (1 - S) ** 2))

    # NET LONG-WAVE RADIATION LOSS WITH SURFACE AT T (B)

    AK = T + 273
    B = (const.SB * AK ** 4 *
         (1 - (0.71 + 0.007 * VD * param.P / 1013) * (1 + RHO)))
    if B < (0.03 * const.SB * AK ** 4):
        B = 0.03 * const.SB * AK ** 4

    # LINE BEFORE CHANGE 0.03 WAS 0.O5
    # NET RADIATION OR NET AVALAIBLE ENERGY WITH SOIL-PLANT
    # SURFACES, WET SURFACE OR LAKE SURFACE AT TEMPERATURE T(RT)

    RT = GL - B

    # SLOPE OF SATURATION VAPOUR PRESSURE CURVE AT T(DELTA)

    DELTA = (const.ALPHA[J] * const.BETA[J] * V / (T + const.BETA[J]) ** 2)

    # STABILITY FACTOR (ZETA)

    EE = const.FTZ[J] * (V - VD)
    RTC = RT
    if RTC < 0:
        RTC = 0.0
    ZETA = (
        0.28 * (1 + VD / V) +
        const.FZ / 28. * DELTA * RTC / const.GAMMA[J] / EE) ** -1
    if ZETA < 1:
        ZETA = 1.0

    # VAPOUR TRANSFER COEFFICIENT (FT)

    FT = const.FTZ[J] / ZETA

    # HEAT TRANSFER COEFFICIENT (LAMDA)

    LAMDA = const.GAMMA[J] + (const.SB * 4) * (T + 273) ** 3 / FT

    # ITERATIVE PROCEDURE FOR COMPUTING POTENTIAL EVAPOTRANSPIRATION
    # BY COMBINING ENERGY BUDGET AND VAPOUR TRANSFER EQUATIONS.
    # ITERATION CONTINUES UNTIL INCREMENT (TDEL) OF POTENTIAL
    # EVAPOTRANSPIRATION EQUILIBRIUM TEMPERATURE (TP) IS < 0.01 DEGREES C.
    
    VP = V
    TP = T
    DELP = DELTA
    ALPHA = const.ALPHA[J]
    BETA = const.BETA[J]
    while True:
        TDEL = ((RT / FT) + VD + LAMDA * (T - TP) - VP) / (DELP + LAMDA)
        TP += TDEL
        VP = (6.11 * math.exp(ALPHA * TP / (TP + BETA)))
        DELP = (ALPHA * BETA * VP / (TP + BETA) ** 2)
        if abs(TDEL) < 0.01:
            break

    # POWER EQUIVALENT OF POTENTIAL EVAPOTRANSPIRATION,
    # PAN-SIZE WET SURFACE EVAPORATION OR POTENTIAL EVAPORATION(ETP)

    ETP = RT - FT * LAMDA * (TP - T)

    # NET RADIATION AT TP (RTP)

    RTP = ETP + FT * const.GAMMA[J] * (TP - T)

    # POWER EQUIVALENT OF WET ENVIRONMENT AREAL EVAPOTRANS-
    # PIRATION, LAKE-SIZE WET SURFACE EVAPORATION OR LAKE
    # EVAPORATION(ETW)

    ETW = (const.CONST1 + const.CONST2 * DELP *
           RTP / (DELP + const.GAMMA[J]))
    if param.LK == 0 and ETW < (0.5 * ETP):
        ETW = ETP / 2
    if not (GL <= GW or RT <= 0 or TP <= T):
        if ETP <= ETW:
            ETP = ETW
    if ETW > ETP:
        ETW = ETP

    # POWER EQUIVALENT OF AREAL EVAPOTRANSPIRATION, LAKE-SIZE
    # WET SURFACE EVAPORATION OR LAKE EVAPORATION(ET)

    if param.LK == 0:
        ET = ETW * 2 - ETP
    if param.LK > 0:
        ET = ETW

    # CONVERT UNIT FROM WATT PER SQUARE METRE INTO MILLIMETRE
    # DEADBEEF - This could probably be simplified

    LTHEAT = 28.5
    if T < 0:
        LTHEAT *= 1.15
    if not (param.LK == 0 or GL <= GW or RT <= 0 or TP <= T):
        LTHEAT = 28.5
    RTMM = LENGTH * RT / LTHEAT
    ETMM = (LENGTH * ET / LTHEAT) / (1 + param.SALT / 1E6)
    ETPMM = (LENGTH * ETP / LTHEAT) / (1 + param.SALT / 1E6)

    # ACCUMULATE ETMM,ETPMM AND RTMM IN TABLES FOR MONTHLY SUMMARY
    # (IF MEANS = 1)

    if param.ISUM == 1:
        data.TRTM[DATE.month][DATE.year].append(RTMM)
        data.TETPM[DATE.month][DATE.year].append(ETPMM)
        data.TETM[DATE.month][DATE.year].append(ETMM)

    data.ETMM[i] = ETMM
    data.RTMM[i] = RTMM
    data.ETPMM[i] = ETPMM
    return

def compute_available_heat():
    """
    CALCULATES DELAYED SOLAR & WATER-BORNE HEAT(GWT)
        BY USING SOLAR & WATER-BORNE HEAT INPUT(TGW)
    CALCULATES AVAILABLE SOLAR & WATER-BORNE HEAT(TGL)
        BY USING DELAYED SOLAR & WATER-BORNE HEAT (GWT)
    """

    # initialize month count

    J = 12

    # CALCULATE DELAYED SOLAR & WATER-BORNE HEAT(GWT) BY
    # USING SOLAR & WATER-BORNE HEAT INPUT(TGW)
    # This could probably be vectorized using numpy arrays

    for I in range(param.NN):
        SLT = 0.13 * data.reservoir_depth[I]
        if SLT > (0.96 + 0.013 * data.reservoir_depth[I]):
            SLT = 0.96 + 0.013 * data.reservoir_depth[I]
        if SLT < (0.039 * data.reservoir_depth[I]):
            SLT = 0.039 * data.reservoir_depth[I]
        LT = SLT / (1 + (param.SALT / 27000) ** 2)
        INTE = int(LT)
        INT1 = INTE + 1
        if LT == 0:
            INT1 = 0
        FRACT = LT - INTE
        FRACT = LT - INTE
        I1 = J - INT1
        II = J - INTE
        data.GWT[I] = data.TGW[II] + FRACT * (data.TGW[I1] - data.TGW[II])

        # Moved to end of loop to adjust indices from 1's based to 0's based

        J += 1

    # INITIALIZE VARIABLES K AND GLB AND SUBSCRIPTS NNFR, NNTO AND M
    
    NNTO = 12
    K = SLT / (1 + (param.DA / 93) ** 7)
    if param.LK == 2:
        MM = 0
        GLB = 50.0
    else:
        MM = 2
        GLB = param.GLBGN

    # CALCULATE AVAILABLE SOLAR & WATER-BORNE HEAT(TGL) BY USING
    # DELAYED SOLAR & WATER-BORNE HEAT(GWT)

    for I in range(MM, 3):
        if I == 2:
            NNTO = param.NN
        for J in range(NNTO):
            GLE = GLB + (data.GWT[J] - GLB) / (K + 0.5)
            data.TGL[J] = 0.5 * (GLB + GLE)
            GLB = GLE
        if I == 1:
            param.GLBGN = GLB
    param.GLEND = GLE

def print_output():
    """"""
    res_f = open(paths.res, 'w')

    # Write output to separate OUT file for CRLE runs

    if param.LK > 1:
        out_f = open(paths.out, 'w')

    SITE = '{0:20s}'.format(param.SITE)
    PHID = '{0:>6s}{1:<7.2f}'.format(' PHID= ', param.PHID)
    if param.IP == 0:
        LINE1 = '{0:>6s}{1:<7.2f}'.format(' P= ', param.PW)
    elif param.IP == 1:
        LINE1 = '{0:>6s}{1:<7.1f}'.format(' ALTI= ', param.PW)
    NET = '{0:>10s}'.format('NET ')
    if param.LK == 0:
        PPN = '{0:>6s}{1:<7.2f}'.format(' PPN= ', param.PPN)
        res_f.write(
            '{0}{1:15s}{2:15s}{3:15s}{6:15s}{4:10s}{5:^20s}\n'.format(
                SITE, PHID, LINE1, PPN, NET, 'EVAPOTRANSPIRATION', ' '))
    elif param.LK == 1:
        SALT = '{0:>6s}{1:<7.2f}'.format(' SALT= ', param.SALT)
        res_f.write(
            '{0}{1:15s}{2:15s}{3:15s}{6:15s}{4:10s}{5:^20s}\n'.format(
                SITE, PHID, LINE1, SALT, NET, 'WET SURFACE EVAP.', ' '))
    else:
        DA = '{0:>6s}{1:<7.1f}'.format(' DA= ', param.DA)
        SALT = '{0:>6s}{1:<7.2f}'.format(' SALT= ', param.SALT)
        res_f.write(
            '{0}{1:15s}{2:15s}{3:15s}{4:15s}{5:10s}{6:^20s}\n'.format(
                SITE, PHID, LINE1, DA, SALT, NET, '  EVAPORATION'))
    LINE2 = [''] * 4
    if (param.IV == 0 and param.IT == 0):
        LINE2[0] = 'TD'
    elif (param.IV == 0 and param.IT == 1):
        LINE2[0] = 'TDF'
    elif param.IV == 1:
        LINE2[0] = 'VD'
    elif param.IV == 2:
        LINE2[0] = 'RELH'
    if param.IT == 0:
        LINE2[1] = 'T'
    elif param.IT == 1:
        LINE2[1] = 'TF'
    if param.IS == 0:
        LINE2[2] = 'S'
    elif param.IS == 1:
        LINE2[2] = 'HS'
    elif param.IS == 2:
        LINE2[2] = 'GIL'
    elif param.IS == 3:
        LINE2[2] = 'GIJ'
    if param.LK > 1:
        LINE2[3] = 'HADD'
    else:
        LINE2[3] = ''
    LINE2_start = '{0:>10s}{1:>10s}{2:>10s}{3:>10s}'.format(
        'YEAR', 'MONTH', 'STARTDAY', 'LENGTH')
    LINE2_mid = ''.join(['{0:>10s}'.format(i) for i in LINE2])
    if param.LK == 0:
        LINE2_end = '{0:>10s}{1:>10s}{2:>10s}'.format(
            'RAD.', 'POTENT.', 'AREAL')
        res_f.write('{}{}{}\n'.format(LINE2_start, LINE2_mid, LINE2_end))
    elif param.LK == 1:
        LINE2_end = '{0:>10s}{1:>10s}{2:>10s}'.format(
            'RAD.', 'PAN-SIZE.', 'LAKE-SIZE')
        res_f.write('{}{}{}\n'.format(LINE2_start, LINE2_mid, LINE2_end))
    elif param.LK > 1:
        LINE2_end = '{0:>10s}{1:>10s}{2:>10s}{3:>10s}'.format(
            'RAD.', 'POTENT.', 'LAKE', 'GW(W/M*M)')
        res_f.write('{}{}{}\n'.format(LINE2_start, LINE2_mid, LINE2_end))
        
        # Write separate output file for GCM runs
        
        LINE2_end = [
            'TD_C', 'T_C', 'RS_MJ_M2_D', 'HADD', 'NET_RAD_MM',
            'ET_POT_MM', 'ET_LAKE_MM', 'GW_W_M2', 'PPT_MM']
        out_f.write(','.join(LINE2_start.split() + LINE2_end) + '\n')

    for i, dt_start in enumerate(data.DATE):
        MONTH = '{0:>10s}'.format(calendar.month_abbr[dt_start.month].upper())
        DAY = '{0:10d}'.format(dt_start.day)
        LENGTH = '{0:10d}'.format(data.LENGTH[i])

        if dt_start.year < 9900:
            YEAR = '{0:>10d}'.format(dt_start.year)
        else:
            YEAR = '{0:>10s}'.format('AVYR')
        if param.IS == 0:
            SW = '{0:10.1f}'.format(data.SW[i])
        else:
            SW = '{0:10.1f}'.format(data.SW[i])
        RTMM = '{0:10.1f}'.format(data.RTMM[i])
        ETPMM = '{0:10.1f}'.format(data.ETPMM[i])
        ETMM = '{0:10.1f}'.format(data.ETMM[i])
        if param.LK <= 1:
            TDW = '{0:10.1f}'.format(data.TDW[i])
            TW = '{0:10.1f}'.format(data.TW[i])
            res_f.write('{}\n'.format(''.join([
                YEAR, MONTH, DAY, LENGTH, TDW, TW, SW, ' ' * 10,
                RTMM, ETPMM, ETMM])))
        else:
            TDW = '{0:>10s}'.format('{0:9.1f}'.format(data.TTDW[i]))
            TW = '{0:>10s}'.format('{0:9.1f}'.format(data.TTW[i]))
            # DEADBEEF - What does this do?
            # if (i + 12) > param.NN:
            #     GW = '{0:10.2f}'.format(data.GW[i])
            # else:
            #     GW = ''
            HADD = '{0:10.2f}'.format(data.HADD[i])
            GW = '{0:10.2f}'.format(data.GW[i])
            res_f.write('{}\n'.format(''.join([
                YEAR, MONTH, DAY, LENGTH, TDW, TW, SW, HADD,
                RTMM, ETPMM, ETMM, GW])))
                
            # Write separate output file for GCM runs
            
            PPT = '{0:10.2f}'.format(data.PPT[i])
            out_f.write(','.join([
                YEAR, MONTH, DAY, LENGTH, TDW, TW, SW, HADD,
                RTMM, ETPMM, ETMM, GW, PPT]).replace(' ', '') + '\n')
    if param.LK >= 2:
        res_f.write(
            '\n**** GLBGN = {0:10.4f}{2:10s} ***** GLEND = {1:10.4f}\n'.format(
                param.GLBGN, param.GLEND, ' '))
    res_f.close()
    return True

def print_monthly_averages():
    """"""
    # Calculate monhtly averages

    TRTM = data.TRTM
    TETPM = data.TETPM
    TETM = data.TETM
    NYR = 0
    for month in range(1, 13):
        NYR = max(NYR, len(TRTM[month].keys()))
        # Sum values for each year/month
        for year in TRTM[month].keys():
            TRTM[month][year] = sum(TRTM[month][year])
            TETPM[month][year] = sum(TETPM[month][year])
            TETM[month][year] = sum(TETM[month][year])
        # For each month, average values for all years
        TRTM[month] = sum(TRTM[month].values()) / len(TRTM[month].values())
        TETPM[month] = sum(TETPM[month].values()) / len(TETPM[month].values())
        TETM[month] = sum(TETM[month].values()) / len(TETM[month].values())

    # Sum of monthly averages

    RTNYR = sum(TRTM.values())
    ETPNYR = sum(TETPM.values())
    ETNYR = sum(TETM.values())

    res_f = open(paths.res, 'a')
    res_f.write('\n\n\n{0:<20s}{1}\n'.format(
        param.SITE,
        '     MONTHLY TOTALS AVERAGED OVER {} YEARS\n'.format(NYR)))
    # NET = '{0:>10s}'.format('NET ')
    if param.LK == 0:
        res_f.write('{0:30s}{1:>10s}{2:^20s}\n'.format(
            ' ', 'NET ', 'EVAPOTRANSPIRATION'))
        res_f.write('{0:20s}{1:>10s}{2:>10s}{3:>10s}{4:>10s}\n\n'.format(
            ' ', 'MONTH', 'RAD.', 'POTENT.', 'AREAL'))
    elif param.LK == 1:
        res_f.write('{0:30s}{1:>10s}{2:^20s}\n'.format(
            ' ', 'NET ', 'WET SURFACE EVAP.'))
        res_f.write('{0:20s}{1:>10s}{2:>10s}{3:>10s}{4:>10s}\n\n'.format(
            ' ', 'MONTH', 'RAD.', 'PAN-SIZE', 'LAKE-SIZE'))
    else:
        res_f.write('{0:30s}{1:>10s}{2:^20s}\n'.format(
            ' ', 'NET ', '  EVAPORATION'))
        res_f.write('{0:20s}{1:>10s}{2:>10s}{3:>10s}{4:>10s}\n\n'.format(
            ' ', 'MONTH', 'RAD.', 'POTENT.', 'LAKE'))
    for i in range(1, 13):
        MONTH = '{}'.format(calendar.month_abbr[i].upper())
        res_f.write('{0:>30s}{1:10.1f}{2:10.1f}{3:10.1f}\n'.format(
            MONTH, TRTM[i], TETPM[i], TETM[i]))
    res_f.write('{0:30s}{1:>10s}{1:>10s}{1:>10s}\n\n'.format(' ', '_' * 8))
    res_f.write(
        '{0:>30s}{1:>10.1f}{2:>10.1f}{3:>10.1f}\n'.format(
            'TOTAL OF MONTHLY AVERAGES', RTNYR, ETPNYR, ETNYR))
    res_f.close()
    return

def read_reservoir_depths():
    """Read reservoir depths as monthly values and transform to input timestamps

    Args:

    Returns:
        success: True or False
    """
    logging.debug('Reading reservoir depths')
    if param.rd_ds_type == 'SF P':
        success, res_depth_df = SF_P_reservoir_depths()
    else:
        success, res_depth_df = DMI_reservoir_depths()
    if not success:
        logging.error('Unable to read reservoir depths.')
        return False
    # print "read res_depth_df\n", res_depth_df.head(2), "\n", res_depth_df.tail(2)
        
    # Check/modify units if necessary
        
    if param.RDU == 1:
        res_depth_df[param.rd_field] = res_depth_df[param.rd_field] / 3.2808399
    # print "converted res_depth_df\n", res_depth_df.head(2), "\n", res_depth_df.tail(2)

    # integrate with main data

    res_depth_df['year'] = res_depth_df.index.year
    res_depth_df['month'] = res_depth_df.index.month
    res_depth_df['day'] = 1
    res_depth_df['bom_date'] = res_depth_df[['year', 'month', 'day']].apply(
            lambda s : dt.datetime(*s), axis = 1)
    # print "final res_depth_df\n", res_depth_df.head(2), "\n", res_depth_df.tail(2)
    known_xs = pd.to_datetime(res_depth_df['bom_date'].values).to_julian_date()
    unknown_xs = pd.to_datetime(np.asarray(data.DATE, dtype = np.datetime64)).to_julian_date()
    data.reservoir_depth = np.interp(unknown_xs, known_xs, res_depth_df[param.rd_field].values)
    # print "data.reservoir_depth\n", data.reservoir_depth
    del res_depth_df, known_xs, unknown_xs
    return True

def SF_P_reservoir_depths():
    """Read reservoir depths in station files with all parameters

    Args:

    Returns:
        success: True or False
        res_depth_df: reservoir depths dataframe
    """
    res_depth_df = None
    res_depth_path = os.path.join(param.rd_folder, param.rd_name_format % param.site_fn)
    if not os.path.isfile(res_depth_path):
        logging.error('ERROR:  input met file {} does not exist'.format(res_depth_path))
        return False, res_depth_df
    logging.debug('  {0}'.format(res_depth_path))

    # Get list of 0 based line numbers to skip
    # Ignore header but assume header was set as 1's based index
    data_skip = [i for i in range(param.rd_header_lines) if i + 1 <> param.rd_names_line]
    res_depth_df = pd.read_table(res_depth_path, engine = 'python',
            header = param.rd_names_line - len(data_skip) - 1, 
            skiprows = data_skip, sep = param.rd_delimiter, 
            na_values = 'NaN')
    logging.debug('  Columns: {0}'.format(', '.join(list(res_depth_df.columns))))

    # verify that reservoir depth column exist
        
    if not param.rd_field in res_depth_df.columns:
        logging.error('Reservoir depth column ' + param.rd_field + ' was not found.')
        return False, res_depth_df

    # Convert date strings to datetimes and index on date
        
    if param.rd_date is not None:
        res_depth_df['date'] = pd.to_datetime(res_depth_df[param.rd_date])
    else:
        res_depth_df['date'] = res_depth_df[[param.rd_year, param.rd_month, param.rd_day]].apply(
            lambda s : dt.datetime(*s),axis = 1)
    res_depth_df.set_index('date', inplace = True)

    # verify period
        
    if param.start_dt is None: 
        pydt = res_depth_df.index[0]
        param.start_dt = pd.to_datetime(dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute))
    if param.end_dt is None: 
        pydt = res_depth_df.index[len(res_depth_df.index) - 1]
        param.end_dt = pd.to_datetime(dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute))
            
    # truncate period
        
    try:
        res_depth_df = res_depth_df.truncate(before = param.start_dt, after = param.end_dt)
    except:
        logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred truncating input reservoir depths')
        return False, res_depth_df
    if len(res_depth_df.index) < 1:
        logging.error('No values found reading reservoir depths')
        return False, res_depth_df
    return True, res_depth_df

def DMI_reservoir_depths():
    """Read reservoir depths for single station using specified DMI format

    Args:

    Returns:
        success: True or False and reservoir depths dataframe
    """

    # Read data from files by fields

    res_depth_df = None
    if '%p' in param.rd_name_format:    # one file per parameter
        res_depth_path = os.path.join(param.rd_folder, 
            param.rd_name_format.replace('%p', param.rd_fnspec))
    else:    # shared file
        res_depth_path = os.path.join(param.rd_folder, param.rd_name_format)
    if not os.path.isfile(res_depth_path):
        logging.error('ERROR:  Reservoir depth file is {} does not exist'.format(res_depth_path))
        return False, res_depth_df
    logging.debug('  Reservoir depth path is {}'.format(res_depth_path))
    if param.rd_ds_type == 'PF S.P':
        if param.rd_file_type == 'csf':
            res_depth_df = ReadOneColumnSlot(res_depth_path, param.rd_header_lines, 
                    param.rd_names_line, param.site_fn, param.rd_field, param.rd_units, 
                    1.0, 'month', 1, param.rd_delimiter, param.start_dt, param.end_dt)
        elif param.rd_file_type == 'rdb':
            res_depth_df = ReadOneTextRDB(res_depth_path, param.rd_header_lines, 
                    param.rd_names_line, param.site_fn, param.rd_field, param.rd_units, 
                    1.0, 'month', 1, param.rd_delimiter, param.start_dt, param.end_dt)
        elif param.rd_file_type == 'xls' or param.rd_file_type == 'wb':
            res_depth_df = ReadOneExcelColumn(res_depth_path, 
                    param.rd_wsspec, param.rd_header_lines, param.rd_names_line, 
                    param.site_fn, param.rd_field, param.rd_units, 
                    1.0, 'month', 1, param.start_dt, param.end_dt)
        else:
            logging.error('ERROR:  File type {} is not supported'.format(param.rd_file_type))
            return False, res_depth_df
        if res_depth_df is None:
            logging.error('ERROR:  unable to read reservoir depths from {}'.format(res_depth_path))
            return False, res_depth_df
        else:
            if param.start_dt is None:
                pydt = res_depth_df.index[0]
                param.start_dt = dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute)
            if param.end_dt is None: 
                pydt = res_depth_df.index[len(res_depth_df) - 1]
                param.end_dt = dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute)
    return True, res_depth_df

def water_year_agg_func(wyem):
    """Sets annual aggregation function for water year end month
    Args:
        wyem: Water Year End Month
    Return; Water year aggregation function
    """
    ann_freqs = ['A-JAN', 'A-FEB', 'A-MAR', 'A-APR', 'A-MAY', 'A-JUN'] + \
                ['A-JUL', 'A-AUG', 'A-SEP', 'A-OCT', 'A-NOV', 'A-DEC']
    return ann_freqs[wyem - 1]

def make_dt_index(time_step, ts_quantity, start_dt, end_dt, wyem = 12):
    """ Make a pandas DatetimeIndex from specified dates, time_step and ts_quantity
    
     Args:
        time_step: RiverWare style string timestep
        ts_quantity: Interger number of time_steps's in interval
        start_dt: starting date time
        end_dt: ending date time
        wyem: Water Year End Month

    Returns:
        dt_index:pandas DatetimeIndex
    """
    dt_index = None
    if time_step == 'day':
        dt_index =pd.date_range(start_dt, end_dt, freq = "D", name = "date")
    elif time_step == 'year':
        dt_index = pd.date_range(start_dt, end_dt, freq = water_year_agg_func(wyem), name = "date")
    elif time_step == 'month':
        dt_index = pd.date_range(start_dt, end_dt, freq = "M", name = "date")
    elif time_step == 'hour':
        if ts_quantity == 1:
            dt_index =pd.date_range(start_dt, end_dt, freq = "H", name = "date")
        else:
            dt_index = pd.date_range(start_dt, end_dt, freq = str(ts_quantity) + "H", name = "date")
    elif time_step == 'minute':
        if ts_quantity == 1:
            dt_index = pd.date_range(start_dt, end_dt, freq = "T", name = "date")
        else:
            dt_index = pd.date_range(start_dt, end_dt, freq = str(ts_quantity) + "T", name = "date")
    elif time_step == 'week':
        dt_index = pd.date_range(start_dt, end_dt, freq = "W", name = "date")
    else:
        logging.error('\nERROR: Timestep {} and ts quantity {} are an invalid combination', format(time_step, ts_quantity))
    return dt_index

def make_ts_dataframe(time_step, ts_quantity, start_dt, end_dt, wyem = 12):
    """ Make a pandas dataframe from specified dates, time_step and ts_quantity
    
     Args:
        time_step: RiverWare style string timestep
        ts_quantity: Interger number of time_steps's in interval
        start_dt: starting date time
        end_dt: ending date time
        wyem: Water Year End Month

    Returns:
        Empty pandas datafame wihh indexed dates
    """
    ts_dataframe = None
    if time_step == 'day':
        ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = "D", name = "date"))
    elif time_step == 'year':
        ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = water_year_agg_func(wyem), name = "date"))
    elif time_step == 'month':
        ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = "M", name = "date"))
    elif time_step == 'hour':
        if ts_quantity == 1:
            ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = "H", name = "date"))
        else:
            ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = str(ts_quantity) + "H", name = "date"))
    elif time_step == 'minute':
        if ts_quantity == 1:
            ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = "T", name = "date"))
        else:
            ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = str(ts_quantity) + "T", name = "date"))
    elif time_step == 'week':
        ts_dataframe = pd.DataFrame(index = pd.date_range(start_dt, end_dt, freq = "W", name = "date"))
    else:
        logging.error('\nERROR: Timestep {} and ts quantity {} are an invalid combination', format(time_step, ts_quantity))
    return ts_dataframe

def ReadOneColumnSlot(file_path, header_lines, names_line, stationToRead, 
        parameterToRead, units, scaleFactor, time_step, ts_quantity, 
        valuesSeparator, start_dt = None, end_dt = None, 
        mia_value = 'NaN', wyem = 12):
    """Reads column slot data for one station and parameter
    
     Args:
        file_path: fully specified file path
        header_lines: number of header lines
        names_line: line of header names
        stationToRead: station to read
        parameterToRead: parameter to read
        units: units of parameter
        scaleFactor: scale to apply to parameter
        time_step: RiverWare style string timestep
        ts_quantity: Interger number of time_steps's in interval
        valuesSeparator: separator of values
        start_dt: starting date time
        end_dt: ending date time
        mia_value: missing value
        wyem: Water Year End Month

    Returns:
        populated dataframe
    """
    return_df = None
    lc_station = stationToRead.lower()
    lc_param = parameterToRead.lower()
    try:
        # Get list of 0 based line numbers to skip - Ignore header but assume header was set as 1's based index
        data_skip = [i for i in range(header_lines) if i + 1 <> names_line]
        input_df = pd.read_table(file_path, engine = 'python',
                header = names_line - len(data_skip) - 1, skiprows = data_skip, 
                sep = valuesSeparator, na_values = mia_value)
        if input_df.empty:
            logging.error("No data read in file" + file_path)
            return return_df
        input_columns = list(input_df.columns)
        lc_columns = [x.lower() for x in input_columns]
        
        # determine date column

        try:
            date_column = lc_columns.index('date')
        except:
            date_column = 0
        date_column_name = input_columns[date_column]
        input_columns.remove(date_column_name)
        lc_columns.remove(date_column_name.lower())
            
        # set date column as index

        input_df = input_df.rename(columns = {date_column_name:'date'})
        
        # make sure that daily, monthly and annual data use end of period dates and do not include a time stamp
        
        input_df['date'] = pd.to_datetime(input_df['date'])
        input_df.set_index('date', inplace = True)
        if time_step == 'day' or time_step == 'month' or time_step == 'year':
            input_df['year'] = input_df.index.year
            input_df['month'] = input_df.index.month
            if time_step == 'day':
                input_df['day'] = input_df.index.day
            elif time_step == 'month':
                input_df['day'] = input_df.index.days_in_month
            else:
                pydt = input_df.index[len(input_df) - 1]
                pydt = pd.to_datetime(dt.datetime(2000, wyem, 1, pydt.hour, pydt.minute))
                pydt = pd.to_datetime(dt.datetime(pydt.year, pydt.month, pydt.days_in_month, pydt.hour, pydt.minute))
                input_df['day'] = pydt.days_in_month
            input_df['date'] = input_df[['year', 'month', 'day']].apply(
                lambda s : dt.datetime(*s),axis = 1)
            input_df['date'] = pd.to_datetime(input_df['date'])
            input_df.set_index('date', inplace = True)
        
        # verify period
        
        if start_dt is None:
            pydt = input_df.index[0]
            start_dt = pd.to_datetime(dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute))
        if end_dt is None: 
            pydt = input_df.index[len(input_df) - 1]
            end_dt = pd.to_datetime(dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute))
        try:
            input_df = input_df.truncate(before = start_dt, after = end_dt)
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred truncating input data')
            return return_df
        if len(input_df.index) < 1:
           input_df = make_ts_dataframe(time_step, ts_quantity, start_dt, 
                       end_dt, wyem)
        
        # adjust for missing rows
        
        full_index = make_dt_index(time_step, ts_quantity, start_dt, end_dt, wyem)
        full_index = full_index + pd.Timedelta(full_index[0] - input_df.index[0])
        input_df = input_df.reindex(index = full_index)
    
        # determine values column

        notFound = True
        for column, input_column in enumerate(lc_columns):
            if lc_station in input_column and lc_param in input_column:
                notFound = False
                break
        if notFound:
            logging.error("Unable to locate station " + stationToRead + " and parameter " + parameterToRead + " in file " + file_path + ".")
            return return_df
            
        # merge values
        
        column_name = input_columns[column]
        try:
            return_df = pd.merge(make_ts_dataframe(time_step, ts_quantity, start_dt, end_dt), 
                    input_df[[column_name]], left_index = True, right_index = True)
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred merging input data with return dataframe.\n')
            # raise
            return return_df
        del input_df, full_index
        return_df = return_df.rename(columns = {column_name:parameterToRead})
        
        # scale values

        return_df[parameterToRead] *= scaleFactor
        return return_df
    except:
        logging.error('\nERROR: ' + str(sys.exc_info()[0]) + ' Error occurred reading column slot data from file\n' + file_path)
        return return_df

def ReadOneTextRDB(file_path, header_lines, names_line, stationToRead, 
        parameterToRead, units, scaleFactor, time_step, ts_quantity, 
        valuesSeparator, start_dt = None, end_dt = None, 
        mia_value = 'NaN', wyem = 12):
    """Reads RDB Text Relational Database (Type C database design) for one station and parameter

     Args:
        file_path: fully specified file path
        header_lines: number of header lines
        names_line: line of header names
        stationToRead: station to read
        parameterToRead: parameter to read
        units: units of parameter
        scaleFactor: scale to apply to parameter
        time_step: RiverWare style string timestep
        ts_quantity: Interger number of time_steps's in interval
        valuesSeparator: separator of values
        start_dt: starting date time
        end_dt: ending date time
        mia_value: missing value
        wyem: Water Year End Month

    Returns:
        populated dataframe
    """
    return_df = None
    lc_station = stationToRead.lower()
    lc_param = parameterToRead.lower()
    try:
        if names_line == 0:
            # default column names and locations
            
            input_df = pd.read_table(file_path, engine = 'python', 
                    header = None, sep = valuesSeparator, na_values = mia_value)
            if input_df.empty:
                logging.error("No data read in file" + file_path)
                return return_df
            sta_column_name = "Station"
            param_column_name = "Parameter"
            date_column_name = "Date"
            values_column_name = "Value"
            input_columns = ['Station', 'Parameter', 'Date', 'Value']
            input_df.columns = input_columns
        else:
            # dynamic column names and location
            
            # Get list of 0 based line numbers to skip - Ignore header but assume header was set as 1's based index
            data_skip = [i for i in range(header_lines) if i + 1 <> names_line]
            input_df = pd.read_table(file_path, engine = 'python',
                    header = names_line - len(data_skip) - 1, skiprows = data_skip, 
                    na_values = mia_value, sep = valuesSeparator)
            if input_df.empty:
                logging.error("No data read in file" + file_path)
                return return_df
            input_columns = list(input_df.columns)
            lc_columns = [x.lower() for x in input_columns]
        
            # determine column types

            try:
                sta_column = lc_columns.index('station')
            except:
                try:
                    sta_column = lc_columns.index('object')
                except:
                    stacolumn = 0
            sta_column_name = input_columns[sta_column]
            try:
                param_column = lc_columns.index('parameter')
            except:
                try:
                    param_column = lc_columns.index('slot')
                except:
                    stacolumn = 1
            param_column_name = input_columns[param_column]
            try:
                date_column = lc_columns.index('date')
            except:
                date_column = 2
            date_column_name = input_columns[date_column]
            try:
                values_column = lc_columns.index('value')
            except:
                values_column = 3
            values_column_name = input_columns[values_column]
        
        # create new column of lower case station and parameter
        
        try:
            input_df["sta_param"] = (input_df[sta_column_name].map(str) + "." 
                                      + input_df[param_column_name].map(str)).str.lower()
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred creating sta_param column.\n')
            return return_df
    
        # locate requested station and parameter

        sta_params = list(pd.unique(input_df.sta_param.ravel()))
        notFound = True
        for sta_param in sta_params:
            if lc_station in sta_param and lc_param in sta_param:
                notFound = False
                break
        if notFound:
            logging.error("Unable to locate station " + stationToRead + " and parameter " + parameterToRead + " in file " + file_path + ".")
            return return_df
    
        # filter data to requested station and parameter values
        
        try:
            input_df = input_df[input_df.sta_param == sta_param]
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred filtering to ', sta_param, '.\n')
            return return_df

        # set date column as index

        input_df = input_df.rename(columns = {date_column_name:'date'})
        if time_step == 'year' and len(str(input_df['date'][0])) == 4:
            input_df['date'] = pd.to_datetime(input_df['date'], format = '%Y')
        else:
            input_df['date'] = pd.to_datetime(input_df['date'])
        input_df.set_index('date', inplace = True)
        
        # set starting and ending dates
        
        if start_dt is None:
            pydt = input_df.index[0]
            start_dt = pd.to_datetime(dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute))
        if end_dt is None: 
            pydt = input_df.index[len(input_df) - 1]
            end_dt = pd.to_datetime(dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute))
        try:
            input_df = input_df.truncate(before = start_dt, after = end_dt)
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred truncating input data')
            return return_df
        
        # adjust for missing rows
        
        full_index = make_dt_index(time_step, ts_quantity, start_dt, end_dt, wyem)
        full_index = full_index + pd.Timedelta(full_index[0] - input_df.index[0])
        input_df = input_df.reindex(index = full_index)

        # merge values
        
        try:
            return_df = pd.merge(make_ts_dataframe(time_step, ts_quantity, start_dt, end_dt), 
                    input_df[[values_column_name]], left_index = True, right_index = True)
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred merging input data with return dataframe.\n')
            # raise
            return return_df
        del input_df, full_index
        return_df = return_df.rename(columns = {values_column_name:parameterToRead})
        
        # scale values

        return_df[parameterToRead] *= scaleFactor
        return return_df
    except:
        logging.error('\nERROR: ' + str(sys.exc_info()[0]) + ' Error occurred reading RDB data from file\n' + file_path)
        return return_df

def ReadOneExcelColumn(file_path, ws_name, header_lines, names_line, stationToRead, 
        parameterToRead, units, scaleFactor, time_step, ts_quantity, 
        start_dt = None, end_dt = None, mia_value = 'NaN', wyem = 12):
    """Reads one Excel column station and parameter

     Args:
        file_path: fully specified file path
        ws_name: worksheet name
        header_lines: number of header lines
        names_line: line of header names
        stationToRead: station to read
        parameterToRead: parameter to read
        units: units of parameter
        scaleFactor: scale to apply to parameter
        time_step: RiverWare style string timestep
        ts_quantity: Interger number of time_steps's in interval
        valuesSeparator: separator of values
        start_dt: starting date time
        end_dt: ending date time
        mia_value: missing value
        wyem: Water Year End Month

    Returns:
        populated dataframe
    """
    return_df = None
    lc_station = stationToRead.lower()
    lc_param = parameterToRead.lower()
    try:
        # Get list of 0 based line numbers to skip - Ignore header but assume header was set as 1's based index
        data_skip = [i for i in range(header_lines) if i + 1 <> names_line]
        input_df = pd.read_excel(file_path, sheetname = ws_name, index_col = 0,
                    header = names_line - len(data_skip) - 1, 
                    skiprows = data_skip, na_values = mia_value)
        if input_df.empty:
            logging.error("No data read in file" + file_path)
            return return_df
                
        # deal with excess rows at bottom that show up as NaT dates

        input_df = input_df[pd.notnull(input_df.index)]
            
        # Deal with possible exitence of 23:59 hour:minute in workbook dates from RiverWare dates
            
        input_df.index.names = ['date']
        input_df['year'] = input_df.index.year
        input_df['month'] = input_df.index.month
        input_df['day'] = input_df.index.day
        if time_step == 'hour' or time_step == 'minute':
            input_df['hour'] = input_df.index.hour
            input_df['minute'] = input_df.index.minute
        else:
            input_df['hour'] = 0
            input_df['minute'] = 0
        input_df['Date'] = input_df[['year', 'month', 'day', 'hour', 'minute']].apply(lambda s : dt.datetime(*s),axis = 1)
        input_df.reset_index('date', inplace = True, drop = True)
        input_df.set_index('Date', inplace = True)
        input_df.drop(['year', 'month', 'day', 'hour', 'minute'], axis = 1, inplace = True)
        input_columns = list(input_df.columns)
        lc_columns = [x.lower() for x in input_columns]
        
        # set starting and ending dates

        if start_dt is None:
            pydt = input_df.index[0]
            if time_step == 'hour' or time_step == 'minute':
                sdt = dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute)
            else:
                sdt = dt.datetime(pydt.year, pydt.month, pydt.day, 0, 0)
            start_dt = pd.to_datetime(sdt)
        else:
            if time_step == 'day' or time_step == 'month' or time_step == 'year':
                start_dt = pd.to_datetime(dt.datetime(start_dt.year, start_dt.month, start_dt.day, 0, 0))
        if end_dt is None: 
            pydt = input_df.index[len(input_df) - 1]
            if time_step == 'hour' or time_step == 'minute':
                edt = dt.datetime(pydt.year, pydt.month, pydt.day, pydt.hour, pydt.minute)
            else:
                edt = dt.datetime(pydt.year, pydt.month, pydt.day, 23, 59)
            end_dt = pd.to_datetime(edt)
        else:
            if time_step == 'day' or time_step == 'month' or time_step == 'year':
                end_dt = pd.to_datetime(dt.datetime(end_dt.year, end_dt.month, end_dt.day, 23, 59))
        try:
            input_df = input_df.truncate(before = start_dt, after = end_dt)
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred truncating input data')
            return return_df
    
        # determine values column

        notFound = True
        for column, input_column in enumerate(lc_columns):
            if lc_station in input_column and lc_param in input_column:
                notFound = False
                break
        if notFound:
            logging.error("Unable to locate station " + stationToRead + " and parameter " + parameterToRead + " in file " + file_path + ".")
            return return_df
        
        # adjust for missing rows
        
        full_index = make_dt_index(time_step, ts_quantity, start_dt, end_dt, wyem)
        full_index = full_index + pd.Timedelta(full_index[0] - input_df.index[0])
        input_df = input_df.reindex(index = full_index)
        column_name = input_columns[column]
            
        # merge values
        
        try:
            return_df = pd.merge(make_ts_dataframe(time_step, ts_quantity, start_dt, end_dt), 
                    input_df[[column_name]], left_index=True, right_index=True)
        except:
            logging.error('\nERROR: ' + str(sys.exc_info()[0]) + 'occurred merging input data with return dataframe.\n')
            # raise
            return return_df
        del input_df, full_index
        return_df = return_df.rename(columns = {column_name:parameterToRead})
        
        # scale values

        return_df[parameterToRead] *= scaleFactor
        return return_df
    except:
        logging.error('\nERROR: ' + str(sys.exc_info()[0]) + ' Error occurred reading data from workbook')
        return return_df

def is_valid_file(parser, arg):
    """"""
    if os.path.isfile(arg):
        return arg
    elif os.path.isfile(os.path.abspath(arg)):
        return os.path.abspath(arg)
    else:
        parser.error('\nFile {} does not exist!'.format(arg))

def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description = 'WREVAP GCM - Python',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', metavar = 'INI', type = lambda x: is_valid_file(parser, x),
        required = True, help = 'Input ini file')
    parser.add_argument(
        '--data', metavar = 'CSV', type = lambda x: is_valid_file(parser, x),
        required = True, help = 'Data csv file')
    parser.add_argument(
        '-o', '--overwrite', default = False, action = "store_true",
        help = 'Force overwrite of existing files')
    parser.add_argument(
        '-d', '--debug', default = logging.INFO, const = logging.DEBUG,
        help = 'Debug level logging', action = "store_const", dest = "loglevel")
    args = parser.parse_args()

    # Convert relative paths to absolute paths

    if os.path.isfile(os.path.abspath(args.ini)):
        args.ini = os.path.abspath(args.ini)
    if os.path.isfile(os.path.abspath(args.data)):
        args.data = os.path.abspath(args.data)
    return args

if __name__ == '__main__':
    args = arg_parse()
    logging.basicConfig(level = args.loglevel, format = '%(message)s')
    WREVAP(args.ini, args.data)
