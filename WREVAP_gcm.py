#--------------------------------
# Name:         WREVAP_gcm.py
# Purpose:      WREVAP - Python version of Fortran model
#               Modified to run on large numbers of GCM datasets
# Author:       Charles Morton
# Co-Author     Justin Huntington
# Email:        charles.morton@dri.edu
# Email:        justin.huntington@dri.edu
# Created       2013/06/10
# Python:       2.6
#--------------------------------

import calendar
from collections import defaultdict
import ConfigParser
import csv
import datetime as dt
import logging
import math
import os
import re
import shutil
import struct
import sys

import numpy as np

class paths(): pass
class param(): pass
class const(): pass
class data(): pass

def isfloat(test_value):
    try:
        float(test_value)
        return True
    except ValueError:
        return False

def WREVAP(input_ws, input_name, data_name):  
    #### Initialize file paths
    set_input_paths(input_ws, input_name, data_name)
    #### Initialize Parameters
    get_parameters()
    #### Read in input data
    read_data_file()
    #### CHECK INPUT SPECIFICATIONS AND DO A GENERAL INITIALIZATION
    initialize()

    #### Process data for each time period
    for i, dt_start in enumerate(data.DATE):
        #### SET UP THE DOY INFORMATION
        doy_start = dt_start.timetuple().tm_yday
        dt_stop = dt_start + dt.timedelta(days=(data.LENGTH[i]-1))
        doy_stop = dt_stop.timetuple().tm_yday    
        doy_adjust = INIJUL(dt_start)
        #### COMPUTE THE RADIUS VECTOR OF THE EARTH'S ORBIT(neta) AND
        #### THE DECLINATION OF THE SUN(theta)
        neta, theta = compute_orbit(doy_start, doy_stop, doy_adjust)
        #### CHECK RAW DATA AND IF NECESSARY CONVERT UNITS AND THEN,
        #### COMPUTE THE ABSORBED GLOBAL RADIATION OR
        #### THE SOLAR & WATER-BORNE HEAT INPUT(GW).
        compute_budget_1(i, neta, theta)
        #### DO THE FINAL CALCULATION AND LIST THE INPUTS AND THE
        #### RESULTS FOR EACH PERIOD.
        if param.LK <= 1:
            compute_budget_2(i)
        #### IF LAKE EVAPORATION IS BEING ESTIMATED, THE RESULTS ARE
        #### STORED FOR LATER USE.
        elif param.LK >= 2:
            TRANLK(i)

    #### THE FOLLOWING SECTION DEALS WITH LAKE EVAPORATION MODEL
    #### WHEN LK = 2 OR 3
    if param.LK >= 2:
        #### READ THE 12 SOLAR & WATER-BORNE HEAT INPUT VALUES(TGW) OF THE
        #### PRECEDING TIME PERIODS, OR COPY THE MOST RECENT 12 VALUES OF
        #### TGW INTO THE PRECEDING PERIODS.
        if param.LK == 3:
            read_tgw_file()
        if param.LK == 2:
            for i in xrange(12):
                data.TGW[i] = data.TGW[i+12]            
        #### CALCULATE THE AVAILABLE SOLAR & WATER-BORNE HEAT CALLING IT
        #### TABLE TGL.
        compute_available_heat()
        #### DO THE FINAL CALCULATION AND LIST THE INPUTS AND THE
        #### RESULTS FOR EACH PERIOD.
        for i in xrange(param.NN):
            compute_budget_2(i)
        #### WRITE THE 12 MOST RECENT VALUES OF TGW TO SOL FILE
        write_sol_file()

    #### Print output to results file
    print_output()

    #### LIST THE MONTHLY TOTALS AVERAGED OVER 5 YEARS
    if param.ISUM == 1:
        print_monthly_averages()

def set_input_paths(input_ws, ini_name, data_name):
    input_name = os.path.splitext(data_name)[0]
    paths.csv = os.path.join(input_ws, data_name)
    paths.ini = os.path.join(input_ws, ini_name)
    paths.res = os.path.join(input_ws, input_name+'.RES')
    paths.tgw = os.path.join(input_ws, input_name+'.TGW')
    paths.sol = os.path.join(input_ws, input_name+'.SOL')
    #### DEADBEEF
    paths.out = os.path.join(input_ws, input_name+'.OUT')
     
def read_param(p_str, p_default, p_type, config):
    try:
        if p_type is float:
            p_value = float(config.get('INPUTS', p_str))
        elif p_type is int or p_type is int:
            p_value = int(config.get('INPUTS', p_str))
        else:  ## String type or NoneType
            p_value = config.get('INPUTS',p_str).upper()
    except:
        p_value = p_default
    return p_value

def check_param_in_list(p_value, p_list, p_name):
    if p_value not in p_list:
        error_str = (
            ("\n  ERROR: %s paramter must be set to %s") %
            (p_name, ', '.join(["{0}".format(i) for i in p_list])))
        raise SystemExit(error_str)                 

def get_parameters():
    #### Parameter lists
    LK_list = [0,1,2,3]
    ISUM_list = [0,1]
    IP_list = [0,1]
    IS_list = [0,1,2,3]
    IT_list = [0,1]
    IV_list = [0,1,2]
    
    #### Parameters must be read from file
    if not os.path.isfile(paths.ini):
        logging.error(
            '\nERROR: The input file {0} does not exist\n'.format(path.ini))

    #### Read parameters from file and check for basic errors
    logging.debug('  Reading INI file')
    #### Check that INI file opens and has section entry ([INPUTS])
    #### Get list of all parameter keys
    config = ConfigParser.SafeConfigParser()
    try:
        config.readfp(open(paths.ini))
        config.has_section('INPUTS')
        config_items = config.items('INPUTS')
    except ConfigParser.NoSectionError:
        error_str = (
            ("\n  ERROR: %s\n"+
             "    The paramter INI file is missing a section line\n"+
             "    The first data line in the file needs to be: [INPUTS]\n"+
             "    Try removing the INI file and rebuilding it\n") %
            paths.ini)
        raise SystemExit(error_str)
    except ConfigParser.Error:
        error_str = (
            ("\n  ERROR: %s\n"+
             "    There is an unknown problem with the paramter INI file \n"+
             "    Try removing the INI file and rebuilding it\n") %
            paths.ini)
        raise SystemExit(error_str)
    #### Check that all parameters are present in INI file
    param_list = ['SITE', 'PHID', 'P', 'LK','ISUM', 'IT', 'IS', 'IV', 'IP']
    config_dict = dict(config_items)
    config_keys = [i.upper() for i in dict(config_items).keys()]
    for param_str in param_list:
        if param_str not in config_keys:
            error_str = ((
                "\n  ERROR: The parameter %s was not found in the "+
                "parameter INI file") % param_str)
            raise SystemExit(error_str)
    if (('PPN' not in config_keys) or 
        ('DA' not in config_keys and
         'SALT' not in config_keys)):
            error_str = (
                "\n  ERROR: Either the PPN parameter or the DA and SALT "+
                "parameters must be specified")
            raise SystemExit(error_str)
    #### Read parameters from file
    SITE = read_param('SITE', '', str, config)
    PHID = read_param('PHID', None, float, config)
    if not (-90 <= PHID <= 90):
        error_str = (
            "\n  ERROR: %s paramter must be between -90 and +90" % 'PPN')
        raise SystemExit(error_str)                 
    P = read_param('P', None, float, config)
    if not P or P < 0:
        error_str = ("\n  ERROR: P paramter must be >= 0")
        raise SystemExit(error_str)                 
    LK = read_param('LK', None, float, config)
    check_param_in_list(LK, LK_list, 'LK')
    if LK == 0:
        SALT, DA = 0.0, 0.0 
        PPN  = read_param('PPN', None, float, config)
        if not PPN or PPN < 0:
            error_str = ("\n  ERROR: PPN paramter must be >= 0 when LK==0")
            raise SystemExit(error_str)                 
        DA   = 0.0
        SALT = 0.0
    elif LK in [1,2,3]:
        PPN  = 0.0
        DA   = read_param('DA', None, float, config)
        SALT  = read_param('SALT', None, float, config)
        if not DA or DA <= 0:
            error_str = (
                "\n  ERROR: DA paramter must be > 0 when LK > 0")
            raise SystemExit(error_str)                 
        elif not SALT or SALT < 0:
            error_str = (
                "\n  ERROR: SALT paramter must be >= 0 when LK > 0")
            raise SystemExit(error_str)
    ISUM = read_param('ISUM', 0, int, config)
    IP   = read_param('IP', 0, int, config)
    IS   = read_param('IS', 1, int, config)
    IT   = read_param('IT', 0, int, config)
    IV   = read_param('IV', 0, int, config)
    check_param_in_list(ISUM, ISUM_list, 'ISUM')
    check_param_in_list(IP, IP_list, 'IP')
    check_param_in_list(IS, IS_list, 'IS')
    check_param_in_list(IT, IT_list, 'IT')
    check_param_in_list(IV, IV_list, 'IV')    

    #### Save parameters 
    param.SITE = SITE
    param.PHID = PHID
    param.P = P
    param.PW = P       #### Save original P value
    param.PPN = PPN
    param.DA = DA
    param.SALT = SALT
    param.LK = LK
    param.ISUM = ISUM
    param.IT = IT
    param.IS = IS
    param.IV = IV
    param.IP = IP
            
def read_data_file():
    #### PROCESS ALL THE INPUT DATA
    logging.debug('  Reading DATA file')
    dat_f = open(paths.csv, 'rb')
    dat_lines = [line.strip() for line in dat_f.readlines()]
    dat_f.close()
    #### First line may be SITE name or Fields
    test_line = dat_lines[0].split(',')
    for test_i, test_line in enumerate(dat_lines):
        if ('YEAR' in test_line and 'LENGTH' in test_line and
            'S' in test_line and 'T' in test_line and 'TD' in test_line):
            header_line = test_i
            data_line = test_i+1
            logging.debug(
                '    Assuming line {0} is the field names'.format(header_line+1))
            logging.debug(
                '    Fields: {0}'.format(
                    ' '.join(dat_lines[header_line].split(','))))
    dat_header = dict([
        (item.upper(), i) for i, item in enumerate(
            dat_lines[header_line].split(','))])
    dat_fields = dat_header.keys()
    
    #### Build data arrays based on number of time steps
    param.NN = len(dat_lines) - data_line
    data.DATE   = [''] * param.NN
    data.LENGTH = [0] * param.NN
    data.TD     = [0.0] * param.NN
    data.T      = [0.0] * param.NN
    data.S      = [0.0] * param.NN
    data.TDW    = [0.0] * param.NN
    data.TW     = [0.0] * param.NN
    data.SW     = [0.0] * param.NN
    data.HADD   = [0] * param.NN
    #### DEADBEEF
    data.PPT      = [0.0] * param.NN
    logging.info('    {0} data points in file\n'.format(param.NN))
    
    #### Get column numbers for each field
    def column_index(name_list):
        for name in name_list:
            if name in dat_fields: return dat_header[name]
        return None
    YEAR_i = column_index(['YEAR'])
    MONTH_i = column_index(['MONTH'])
    DAY_i = column_index(['DAY', 'STARTDAY', 'START_DAY'])
    DOY_i = column_index(['DOY', 'STARTDOY', 'START_DOY'])
    LENGTH_i = column_index(['LENGTH'])
    TD_i = column_index(['TD'])
    T_i = column_index(['T'])
    S_i = column_index(['S'])
    HADD_i = column_index(['HADD'])
    PPT_i = column_index(['PPT'])
    
    #### Check date fields
    if YEAR_i >= 0 and DOY_i >= 0 and LENGTH_i >= 0:
        dt_doy_flag = True
        dt_format = '%Y_%j'
    elif YEAR_i >= 0 and DAY_i >= 0 and DAY_i >= 0 and LENGTH_i >= 0:
        dt_doy_flag = False
        dt_format = '%Y_%m_%d'
    else:
        logging.error('  Fields: {0}'.format(dat_fields))
        logging.error(
            '\n  ERROR: Data file is missing start date or length fields\n'+
            '  ERROR: YEAR/STARTDOY/LENGTH or YEAR/MONTH/STARTDAY/LENGTH\n')
        raise SystemExit()
    
    #### Check data fields
    if not (S_i and T_i and TD_i):
        logging.error(
            '\n  ERROR: Data file is missing data fields\n'+
            '  ERROR: TD, T, or S\n')
        raise SystemExit()
    
    #### Read in data
    dat_list = [i.split(',') for i in dat_lines[data_line:]]
    #### Process data for each time period
    for i, dat_line in enumerate(dat_list):
        ##logging.debug('\n  {0}'.format(dat_line))
        year = int(dat_line[YEAR_i])
        if dt_doy_flag:
            doy = int(dat_line[DOY_i])
            dt_start = dt.datetime.strptime(
                '%04d_%03d' % (year, doy), dt_format)
        else:
            month = int(dat_line[MONTH_i])
            day = int(dat_line[DAY_i])
            dt_start = dt.datetime.strptime(
                '%04d_%s_%s' % (year, month, day), dt_format)
        data.DATE[i] = dt_start
        data.LENGTH[i] = int(dat_line[LENGTH_i])
        data.TD[i] = float(dat_line[TD_i])
        data.T[i] = float(dat_line[T_i])
        data.S[i] = float(dat_line[S_i])
        data.HADD[i] = float(dat_line[HADD_i]) if HADD_i >= 0 else 0
        data.PPT[i] = float(dat_line[PPT_i]) if PPT_i >= 0 else 0

def initialize():
    #### INIT => CHECKS THE INPUT SPECIFICATIONS FOR, LOGICAL ERROR
    ####      => MODIFIES THESE SPECIFICATIONS IF NECESSARY,THEN
    ####      =>  INITIALIZE THE CONSTANTS AND TABLES
    P = param.P
    if param.IP == 1:
        P = 1013 * (1 - 0.0065 * P / 288) ** 5.256
    if P < 0:
        logging.error(
            '\n  *** INPUT ERROR NUMBER = 23\n'+
            'Value of P should be greater than or equal to 0\n')
        raise SystemExit()
    param.P = P
    #### DIFFERENT CONSTANTS ARE USED DEPENDING ON WHETHER THE AREAL
    #### EVAPOTRANSPIRATION MODEL IS BEING USED(LK=0) OR THE WET SURFACE
    #### EVAPORATION OR LAKE EVAPORATION MODELS ARE BEING USED(LK>0)
    #### CONSTANTS FOR THE AREAL EVAPOTRANSPIRATION MODELS
    if param.LK == 0:
        const.SB = 5.22 / 1E8
        const.AZZZ = (0.26 - 0.00012 * math.sqrt(P / 1013) * param.PPN *
                      (1 + abs(param.PHID) / 42 + (param.PHID / 42) ** 2))
        const.CONST2 = 1.20
        const.FZ = 28.0
        const.CONST1 = 14.0
    #### CONSTANTS FOR WET SURFACE AND LAKE EVAPORATION MODELS
    else:
        const.SB = 5.5 / 1E8
        const.AZZZ = 0.05
        const.CONST2 = 1.12
        const.FZ = 25.0
        const.CONST1 = 13.0
    #### VARIABLES FOR MONTHLY MEANS ARE INITIALIZED
    #### IF ISUM=0 , NO MONTHLY SUMMARY IS PRINTED AT THE END
    #### IF ISUM=1 , TABLE OF MONTHLY MEANS IS PRINTED
    if param.ISUM == 1:
        data.TETM  = defaultdict(dict)
        data.TETPM = defaultdict(dict)
        data.TRTM  = defaultdict(dict)
        for DATE in data.DATE:
            data.TETM[DATE.month][DATE.year] = []
            data.TETPM[DATE.month][DATE.year] = []
            data.TRTM[DATE.month][DATE.year] = []
    ####
    data.V       = [0.0] * param.NN
    data.VD      = [0.0] * param.NN
    data.GW      = [0.0] * param.NN
    data.ETMM    = [0.0] * param.NN
    data.RTMM    = [0.0] * param.NN
    data.ETPMM   = [0.0] * param.NN
    ####
    if param.LK >= 2:
        data.TDATE   = [0] * param.NN
        data.TLENGTH = [0] * param.NN
        data.TJ      = [0] * param.NN
        data.TT      = [0.0] * param.NN
        data.TSW     = [0.0] * param.NN
        data.TTDW    = [0.0] * param.NN
        data.TTW     = [0.0] * param.NN
        data.TS      = [0.0] * param.NN
        data.THADD   = [0.0] * param.NN
        data.TV      = [0.0] * param.NN
        data.TVD     = [0.0] * param.NN
        #### For Lake Evap model TGW stores data for 12 previous time periods
        data.TGW     = [0.0] * (param.NN + 12)
        data.TGL     = [0.0] * param.NN   
        data.GWT     = [0.0] * param.NN
    #### DEFINE TEMPERATURE DEPENDENT CONSTANTS
    #### 0 - FOR TEMPERATURE GREATER THAN OR EQUAL TO 0 DEGREES CELSIUS.
    #### 1 - FOR TEMPERATURE LESS THAN 0 DEGREES CELSIUS
    const.ALPHA = [17.27, 21.88]
    const.BETA  = [237.3, 265.5]
    const.GAMMA = [(0.66 * P / 1013), ((0.66 * P / 1013) / 1.15)]
    const.FTZ   = [(const.FZ * math.sqrt(1013 / P)),
                   ((const.FZ * math.sqrt(1013 / P)) * 1.15)]
    #### DEFINE OTHER VALUES
    const.PI = math.pi
    const.CONV = math.pi / 180
    param.PHI = param.PHID * const.CONV
    return
    
def read_tgw_file():
    ##if os.path.isfile(paths.sol) and not os.path.isfile(paths.tgw):
    ##    logging.info("  NOTE: Copying data in SOL file to TGW file")
    ##    shutil.copy2(paths.sol, paths.tgw)
    try: tgw_f = open(paths.tgw, 'r')
    except:
        logging.error(
            ("  ERROR: The TGW file does not exist or "+
             "could not be opened\n         ({0})").format(paths.tgw))
        logging.error(
            "         If running LK==3, try running LK==2 first\n"+
            "         Then copy the data in the SOL file to the TGW file")
        logging.error(
            "         See section 6.6 and 6.7 in NHRI Paper No. 24")
        raise SystemExit()
    tgw_data = [float(i) for i in tgw_f.readlines()]
    tgw_f.closed
    if len(tgw_data)<>13:
        logging.error("  ERROR: There should be 13 values in the .TGW file")
        raise SystemExit()
    param.GLBGN = tgw_data.pop(0)
    #### Values are reversed in file?
    ##data.TGW[:12] = tgw_data[:]
    data.TGW[:12] = tgw_data[::-1]
    del tgw_data, tgw_f

def write_sol_file():
    sol_f = open(paths.sol, 'w')
    sol_f.write("{0:9.4f}\n".format(param.GLEND))
    #### Values are reversed in file? (reverse sort then get last 12)
    sol_f.write('\n'.join(map('{0:9.4f}'.format, data.TGW[::-1][:12])))
    sol_f.closed
    del sol_f

def INIJUL(dt_start):
    #### JULIAN DAY ADJUSTMENT FOR LEAP YEARS
    if calendar.isleap(dt_start.year):
        doy_adjust = -0.5
    else:
        doy_adjust = 0.5
    if dt_start.month <= 2: doy_adjust = 0.0
    ##logging.debug('  DOY Adjust: {0}'.format(doy_adjust))
    return doy_adjust

def compute_orbit(doy_start, doy_stop, doy_adjust):
    #### CALCULATES THE ORBITAL PARAMETERS SUCH AS:
    ####   THE DECLINATION OF THE SUN IN RADIAN(theta)
    ####   THE NON-DIMENSIONNAL RADIUS VECTOR OF THE EARTH'S ORBIT
    ####   AROUND THE SUN(NETA)
    doy_array = (np.arange(doy_start, doy_stop+1) + doy_adjust)
    orb1 = np.minimum((29.5 + doy_array / 270.0), 30.4)
    orb2 = ((doy_array + 0.5 * (orb1 - 1)) / orb1)
    neta = 1 + np.sin((29.5 * orb2 - 106.0) * const.CONV) / 60.0
    theta = (23.45 * const.CONV * np.sin((29.5 * orb2 - 94.0) * const.CONV))
    neta = np.mean(neta)
    theta = np.mean(theta)
    return neta, theta
      
def compute_budget_1(i, neta, theta):
    #### => CHECKS LOGICAL ERRORS IN RAW DATA INPUT
    #### => CONVERTS RAW DATA  INTO STANDARD UNITS.
    #### => CALCULATES V,VD,GW
    #### => ASSIGNS VALUES TO  J,TDW,TW,SW
    TD = data.TD[i]
    T = data.T[i]
    S = data.S[i]
    HADD = data.HADD[i]

    #### STORE RAW DATA INPUT FOR PRINT OUT
    data.TDW[i] = TD
    data.TW[i] = T
    data.SW[i] = S

    #### CHECK LOGICAL ERRORS OF RAW DATA INPUT
    if param.IV >= 1 and TD < 0:
        logging.error(
            '\n  *** INPUT ERROR NUMBER = 5\n'+
            '  TD values should be > 0 when TD represents the '+
            'vapour pressure at dew point or relative humidity\n')
        raise SystemExit()
    if S < 0:
        logging.error(
            '\n  *** INPUT ERROR NUMBER = 6\n  S values should be >= 0\n')
        raise SystemExit()

    #### CONVERT FAHRENHEIT TEMPERATURE INTO CELSIUS
    if param.IT <> 0:
        if param.IV == 0: TD = (TD - 32) * 5. / 9.
        T = (T - 32) * 5. / 9.
    J = 1 if T < 0 else 0

    #### SATURATION VAPOUR PRESSURE(V) AT AIR TEMPERATURE
    V = 6.11 * math.exp(const.ALPHA[J] * T / (T + const.BETA[J]))

    #### SATURATION VAPOUR PRESSURE(VD) AT DEW POINT TEMPERATURE
    #### WHEN TD IS DEW POINT
    if param.IV == 0:
        VD = 6.11 * math.exp(17.27 * TD / (TD + 237.3))
    elif param.IV == 1:
        VD = TD
    else:
        VD = TD * V    #### WHEN TD IS THE RELATIVE HUMIDITY

    #### CALCULATE THE SUNSHINE DURATION RATIO, (S), WHEN
    #### THE SUNSHINE DURATION IS GIVEN IN HOURS PER DAY
    if param.IS == 1:
        CZENAA = math.cos(param.PHI - theta) + 0.005
	if CZENAA < 0.001: CZENAA = 0.001
	ACOMA = 1 - CZENAA / math.cos(param.PHI) / math.cos(theta)
	if ACOMA < -1: ACOMA = -1.0
	OMEGAA = math.acos(ACOMA)
	MSD = 24 * OMEGAA / const.PI
	S = S / MSD

    #### SOLAR ZENITH ANGLE (ZENA)
    CZENA = math.cos(param.PHI - theta)
    if CZENA < 0.001: CZENA = 0.001
    ZENA = math.acos(CZENA)

    #### NO. OF DEGREES THE EARTH ROTATES BETWEEN SUNRISE AND NOON (OMEGA)
    ACOM = 1 - (CZENA / math.cos(param.PHI) / math.cos(theta))
    if ACOM < -1: ACOM = -1.0
    OMEGA = math.acos(ACOM)

    #### COSINE OF THE AVERAGE ANGULAR ZENITH DISTANCE OF THE SUN (COSZ)
    COSZ = (CZENA + (math.sin(OMEGA) / OMEGA - 1) *
            math.cos(param.PHI) * math.cos(theta))

    #### THE EXTRA-ATMOSPHERIC GLOBAL RADIATION(GE)
    GE = 1354.0 * COSZ * OMEGA / (const.PI * neta ** 2)

    #### SNOW FREE, CLEAR SKY ALBEDOIF SUN WERE AT  ZENITH (AZZ)
    AZZ = const.AZZZ

    #### CONSTRAINT FOR AZZ WHICH MAY APPLY DURING WET SEASON IN DRY AREA
    if param.LK == 0:
	if (AZZ > ((0.91 - VD / V) / 2)):
            AZZ = (0.91 - VD / V) / 2
        #### TWO FURTHER CONSTRAINTS FOR AZZ
        if AZZ > 0.17: AZZ = 0.17
        if AZZ < 0.11: AZZ = 0.11

    #### WEIGHTING FACTOR FOR THE EFFECT OF SNOW ON ALBEDO (ARAT)
    VPDL = V - VD
    if VPDL < 0: VPDL = 0.0
    if VPDL > 1: VPDL = 1.0
    ARAT = 1 - VPDL * VPDL

    #### CLEAR SKY ALBEDO IF SUN WERE AT ZENITH(AZZ)
    AZ = AZZ + ARAT * (0.34 - AZZ)

    #### CLEAR SKY ALBEDO (A0)
    A0 = (AZ *(math.exp(1.08) - math.exp(ZENA * 2.16 / const.PI) *
               (math.cos(ZENA) * 2.16 / const.PI + math.sin(ZENA))) /
          (1.473 * (1.0 - math.sin(ZENA))))

    #### PRECIPITABLE WATER VAPOUR (W)
    W = VD / (0.49 + T / 129)

    #### WEIGHTING FACTOR FOR THE EFFECTS OF TEMPERATURE ON TURBIDITY (TST)
    TST = 21 - T
    if TST < 0: TST = 0.0
    if TST > 5: TST = 5.0

    #### TURBIDITY COEFFICIENT (DUST)
    DUST = (0.5 + 2.5 * COSZ ** 2) * math.exp(TST * (param.P / 1013 - 1))

    #### TRANSMITTANCY OF CLEAR SKIES TO DIRECT BEAM SOLAR
    #### RADIATION(TAUT)
    DUSTT = 0.083 * (DUST / COSZ) ** 0.9
    WVT = 0.029 * (W / COSZ) ** 0.6
    LNX = (-0.089 * (param.P / (COSZ * 1013)) ** 0.75 - DUSTT - WVT)
    if LNX < -675: LNX = -675.0
    TAUT = math.exp(LNX)

    #### PART OF TAUT THAT IS RESULT OF ABSORBTION(TAUA)
    WVA = math.sqrt(WVT / 10.0)
    if WVA > WVT: WVA = WVT
    LNY = -DUSTT / 2 - WVA
    if LNY < -675: LNY = -675.0
    TAUA = math.exp(LNY)

    #### CLEAR SKY GLOBAL RADIATION (G0)
    G0 = GE * (TAUT + TAUT * (1 - TAUT / TAUA) * (1 + A0 * TAUT))

    #### ESTIMATE THE INCIDENT GLOBAL RADIATION FROM THE GIVEN VALUE OF
    #### SUNSHINE DURATION RATIO (S)
    if param.IS in [0, 1]:
        G = G0 * S + (0.08 + 0.3 * S) * (1 - S) * GE
    #### OR IF IS = 2,THEN S IS THE INCIDENT GLOBAL RADIATION(LANGLEY/DAY)
    #### OR IF IS = 3,THEN S IS THE INCIDENT GLOBAL RADIATION(MJ/M**2/DAY)
    elif param.IS in [2, 3]:
        if param.IS == 2: G = S / 2.064
        if param.IS == 3: G = S / 0.0864
        S = 0.53 * G / (G0 - 0.47 * G)
	if S > 1: S = 1.0
	if S < 0: S = 0.0
    #### ESTIMATE THE AVERAGE ALBEDO (A)
    A = A0 * (S + (1 - ZENA * 180 / const.PI / 330) * (1 - S))

    #### COMPUTE GW WHICH IS THE ABSORBED GLOBAL RADIATION FOR AREAL
    #### EVAPOTRANSPIRATION AND WET SURFACE EVAPORATION OPTIONS, AND
    #### THE SOLAR & WATER-BORNE HEAT INPUT FOR THE LAKE EVAPORATION
    #### OPTION.
    GW = (1 - A) * G
    if param.LK > 1: GW += HADD

    #### Save parameters
    data.TD[i] = TD
    data.T[i] = T
    data.S[i] = S
    data.V[i] = V
    data.VD[i] = VD
    data.GW[i] = GW
      
def TRANLK(i):
    #### THIS SUBROUTINE IS USED ONLY FOR LAKE EVAPORATION MODEL.
    #### ALL THESE TRANLK ARGUMENTS ARE STORE IN TABLES FOR LATER
    #### APPLICATION IN BUDGT2.  THE TABULATED NAMES ARE PREFIXED WITH LETTER
    #### "T" AND FOLLOWED BY THE PREVIOUS NON-SUBSCRIPTED TRANLK ARGUMENTS.
    data.TDATE[i]   = data.DATE[i]
    data.TLENGTH[i] = data.LENGTH[i]
    data.TT[i]      = data.T[i]
    data.TS[i]      = data.S[i]
    data.TTDW[i]    = data.TDW[i]
    data.TTW[i]     = data.TW[i]
    data.TSW[i]     = data.SW[i]
    data.TV[i]      = data.V[i]
    data.TVD[i]     = data.VD[i]
    data.THADD[i]   = data.HADD[i]
    #### For Lake Evap model TGW stores data for 12 previous time periods
    data.TGW[i+12] = data.GW[i]
    return
      
def compute_budget_2(i):
    #### BUDGT2 CALCULATES THE:
    ####    POWER EQUIVALENT OF AREAL EVAPOTRANSPIRATION, LAKE-SIZE
    ####    WET SURFACE EVAPORATION OR LAKE EVAPORATION(ET),
    ####    POWER EQUIVALENT OF POTENTIAL EVAPOTRANSPIRATION,
    ####    PAN-SIZE WET SURFACE EVAPORATION OR POTENTIAL EVAPORATION(ETP),
    ####    NET RADIATION OR NET AVALAIBLE ENERGY WITH SOIL-PLANT
    ####    SURFACES, WET SURFACE OR LAKE SURFACE AT TEMPERATURE T(RT)
    ####    AND CALLS print_output TO LIST INPUTS AND RESULTS FOR EACH TIME PERIOD
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
        GW = data.TGW[i+12] + 1E-11
        GL = data.TGL[i]
    J = 1 if T < 0 else 0
    
    #### WEIGHTING FACTOR FOR THE EFFECT OF CLOUDS ON ATMOSPHERIC RADIATION (ATM)
    ATM = 10 * (VD / V - S - 0.42)
    if ATM < 0: ATM = 0.0
    if ATM > 1: ATM = 1.0

    #### PROPORTIONAL INCREASE IN ATMOSPHERIC RADIATION DUE TO CLOUDS (RHO)
    RHO = (0.18 * 1013 / param.P *
           (ATM * math.sqrt((1 - S)) + (1 - ATM) * (1 - S) ** 2))

    #### NET LONG-WAVE RADIATION LOSS WITH SURFACE AT T (B)
    AK = T + 273
    B = (const.SB * AK ** 4 *
         (1 - (0.71 + 0.007 * VD * param.P / 1013) * (1 + RHO)))
    if B < (0.03 * const.SB * AK ** 4):
        B = 0.03 * const.SB * AK ** 4
    #### LINE BEFOR CHAGE 0.03 WAS 0.O5
    #### NET RADIATION OR NET AVALAIBLE ENERGY WITH SOIL-PLANT
    #### SURFACES, WET SURFACE OR LAKE SURFACE AT TEMPERATURE T(RT)
    RT = GL - B

    #### SLOPE OF SATURATION VAPOUR PRESSURE CURVE AT T(DELTA)
    DELTA = (const.ALPHA[J] * const.BETA[J] * V / (T + const.BETA[J]) ** 2)

    #### STABILITY FACTOR (ZETA)
    EE = const.FTZ[J] * (V - VD)
    RTC = RT
    if RTC < 0: RTC = 0.0
    ZETA = (1. / (0.28 * (1 + VD / V) +
                  const.FZ / 28. * DELTA * RTC / const.GAMMA[J] / EE))
    if ZETA < 1: ZETA = 1.0

    #### VAPOUR TRANSFER COEFFICIENT (FT)
    FT = const.FTZ[J] / ZETA

    #### HEAT TRANSFER COEFFICIENT (LAMDA)
    LAMDA = const.GAMMA[J] + (const.SB * 4) * (T + 273) ** 3 / FT

    #### ITERATIVE PROCEDURE FOR COMPUTING POTENTIAL EVAPOTRANSPIRATION
    #### BY COMBINING THE ENERGY BUDGET AND VAPOUR TRANSFER EQUATIONS.
    #### THE ITERATION CONTINUES UNTIL THE INCREMENT (TDEL) OF THE POTENTIAL
    #### EVAPOTRANSPIRATION EQUILIBRIUM TEMPERATURE (TP) IS < 0.01 DEGREES C.
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
        if abs(TDEL) < 0.01: break
    
    #### POWER EQUIVALENT OF POTENTIAL EVAPOTRANSPIRATION,
    #### PAN-SIZE WET SURFACE EVAPORATION OR POTENTIAL EVAPORATION(ETP)
    ETP = RT - FT * LAMDA * (TP - T)

    #### NET RADIATION AT TP (RTP)
    RTP = ETP + FT * const.GAMMA[J] * (TP - T)

    #### POWER EQUIVALENT OF WET ENVIRONMENT AREAL EVAPOTRANS-
    #### PIRATION, LAKE-SIZE WET SURFACE EVAPORATION OR LAKE
    #### EVAPORATION(ETW)
    ETW = (const.CONST1 + const.CONST2 * DELP *
           RTP / (DELP + const.GAMMA[J]))
    if param.LK == 0 and ETW < (0.5 * ETP):
        ETW = ETP / 2
    if not (GL <= GW or RT <= 0 or TP <= T):
        if ETP <= ETW:  ETP = ETW
    if ETW > ETP: ETW = ETP

    #### POWER EQUIVALENT OF AREAL EVAPOTRANSPIRATION, LAKE-SIZE
    #### WET SURFACE EVAPORATION OR LAKE EVAPORATION(ET)
    if param.LK == 0: ET = ETW * 2 - ETP
    if param.LK > 0: ET = ETW

    #### CONVERT UNIT FROM WATT PER SQUARE METRE INTO MILLIMETRE
    #### DEADBEEF - This could probably be simplified
    LTHEAT = 28.5
    if T < 0: LTHEAT *= 1.15
    if not (param.LK == 0 or GL <= GW or RT <= 0 or TP <= T):
        LTHEAT = 28.5
    RTMM = LENGTH * RT / LTHEAT
    ETMM = (LENGTH * ET / LTHEAT) / (1 + param.SALT / 1E6)
    ETPMM = (LENGTH * ETP / LTHEAT) / (1 + param.SALT / 1E6)

    #### ACCUMULATE ETMM,ETPMM AND RTMM IN TABLES FOR MONTHLY SUMMARY
    #### (IF MEANS = 1)
    if param.ISUM == 1:
        data.TRTM[DATE.month][DATE.year].append(RTMM)
        data.TETPM[DATE.month][DATE.year].append(ETPMM)
        data.TETM[DATE.month][DATE.year].append(ETMM)

    data.ETMM[i] = ETMM
    data.RTMM[i] = RTMM
    data.ETPMM[i] = ETPMM
    return
      
def compute_available_heat():
    ####  => CALCULATES THE DELAYED SOLAR & WATER-BORNE HEAT(GWT)
    ####     BY USING THE SOLAR & WATER-BORNE HEAT INPUT(TGW)
    ####  => CALCULATES THE AVAILABLE SOLAR & WATER-BORNE HEAT(TGL)
    ####     BY USING THE DELAYED SOLAR & WATER-BORNE HEAT (GWT)

    #### INITIALIZE THE VARIABLE FRACT AND THE SUBSCRIPTS INTE AND INT1
    J = 12
    SLT = 0.13 * param.DA
    if SLT > (0.96 + 0.013 * param.DA):
        SLT = 0.96 + 0.013 * param.DA
    if SLT < (0.039 * param.DA):
        SLT = 0.039 * param.DA
    LT = SLT / (1 + (param.SALT / 27000) ** 2)
    INTE = int(LT)
    INT1 = INTE+1
    if LT == 0: INT1 = 0 
    FRACT = LT - INTE

    #### CALCULATE THE DELAYED SOLAR & WATER-BORNE HEAT(GWT) BY
    #### USING THE SOLAR & WATER-BORNE HEAT INPUT(TGW)
    #### This could probably be vectorized using numpy arrays
    for I in range(param.NN):
	 I1 = J - INT1
	 II = J - INTE
	 data.GWT[I] = data.TGW[II] + FRACT * (data.TGW[I1] - data.TGW[II])
         #### Moved to end of loop to adjust indices from 1's based to 0's based
	 J = J + 1

    #### INITIALIZE THE VARIABLES K AND GLB AND THE SUBSCRIPTS NNFR, NNTO AND M
    NNTO = 12
    K = SLT / (1 + (param.DA / 93) ** 7)
    if param.LK == 2:
        MM = 0
        GLB = 50.0
    else:
        MM = 2
        GLB = param.GLBGN

    #### CALCULATE THE AVAILABLE SOLAR & WATER-BORNE HEAT(TGL) BY  USING
    #### THE DELAYED SOLAR & WATER-BORNE HEAT(GWT)
    for I in range(MM, 3):
        if I == 2: NNTO = param.NN
        for J in range(NNTO):
            GLE = GLB + (data.GWT[J] - GLB) / (K + 0.5)
            data.TGL[J] = 0.5 * (GLB + GLE)
            GLB = GLE
        if I == 1: param.GLBGN = GLB
    param.GLEND = GLE
      
def print_output():
    res_f = open(paths.res, 'w')
    #### DEADBEEF
    if param.LK > 1:
        out_f = open(paths.out, 'w')

    SITE = '{0:20s}'.format(param.SITE)
    PHID = '{0:>6s}{1:<7.2f}'.format(' PHID= ', param.PHID)
    if param.IP == 0: LINE1 = '{0:>6s}{1:<7.2f}'.format(' P= ', param.PW)
    if param.IP == 1: LINE1 = '{0:>6s}{1:<7.1f}'.format(' ALTI= ', param.PW)
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

    LINE2 = ['']*4
    if (param.IV == 0 and param.IT == 0): LINE2[0] = 'TD'
    if (param.IV == 0 and param.IT == 1): LINE2[0] = 'TDF'
    if param.IV == 1: LINE2[0] = 'VD'
    if param.IV == 2: LINE2[0] = 'RELH'

    if param.IT == 0: LINE2[1] = 'T'
    if param.IT == 1: LINE2[1] = 'TF'
    if param.IS == 0: LINE2[2] = 'S'
    if param.IS == 1: LINE2[2] = 'HS'
    if param.IS == 2: LINE2[2] = 'GIL'
    if param.IS == 3: LINE2[2] = 'GIJ'

    if param.LK > 1: LINE2[3] = 'HADD'
    else: LINE2[3] = ''

    LINE2_start = '{0:>10s}{1:>10s}{2:>10s}{3:>10s}'.format(
        'YEAR', 'MONTH', 'STARTDAY', 'LENGTH')
    LINE2_mid = ''.join(['{0:>10s}'.format(i) for i in LINE2])
    if param.LK == 0:
        LINE2_end = '{0:>10s}{1:>10s}{2:>10s}'.format(
            'RAD.', 'POTENT.', 'AREAL')
        res_f.write('%s%s%s\n' % (LINE2_start, LINE2_mid, LINE2_end))
    if param.LK == 1:
        LINE2_end = '{0:>10s}{1:>10s}{2:>10s}'.format(
            'RAD.', 'PAN-SIZE.', 'LAKE-SIZE')
        res_f.write('%s%s%s\n' % (LINE2_start, LINE2_mid, LINE2_end))
    if param.LK > 1:
        LINE2_end = '{0:>10s}{1:>10s}{2:>10s}{3:>10s}'.format(
            'RAD.', 'POTENT.', 'LAKE', 'GW(W/M*M)')
        res_f.write('%s%s%s\n' % (LINE2_start, LINE2_mid, LINE2_end))
        #### DEADBEEF
        LINE2_end = [
            'TD_C', 'T_C', 'RS_MJ_M2_D', 'HADD', 'NET_RAD_MM', 
            'ET_POT_MM', 'ET_LAKE_MM', 'GW_W_M2', 'PPT_MM']
        out_f.write(','.join(LINE2_start.split() + LINE2_end)+'\n')
    del LINE2_start, LINE2_mid, LINE2_end, LINE2


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
            res_f.write('%s\n'% ''.join([
                YEAR, MONTH, DAY, LENGTH, TDW, TW, SW, ' '*10, RTMM, ETPMM, ETMM]))
        else:
            TDW = '{0:>10s}'.format('{0:9.1f}'.format(data.TTDW[i]))
            TW = '{0:>10s}'.format('{0:9.1f}'.format(data.TTW[i]))
            #### DEADBEEF - What does this do?
            ##if (i + 12) > param.NN:
            ##    GW = '{0:10.2f}'.format(data.GW[i])
            ##else:
            ##    GW = ''
            HADD = '{0:10.2f}'.format(data.HADD[i])
            GW = '{0:10.2f}'.format(data.GW[i])
            res_f.write('%s\n'% ''.join([
                YEAR, MONTH, DAY, LENGTH, TDW, TW, SW, HADD, RTMM, ETPMM, ETMM, GW]))
            #### DEADBEEF
            PPT = '{0:10.2f}'.format(data.PPT[i])
            out_f.write(','.join([
                YEAR, MONTH, DAY, LENGTH, TDW, TW, SW, HADD, 
                RTMM, ETPMM, ETMM, GW, PPT]).replace(' ','')+'\n')
            del HADD, GW
            del PPT
        del YEAR, MONTH, DAY, LENGTH, SW
        del RTMM, ETPMM, ETMM
        del TDW, TW

    if param.LK >= 2:
        res_f.write(
            '\n**** GLBGN = {0:10.4f}{2:10s} ***** GLEND = {1:10.4f}\n'.format(
                param.GLBGN, param.GLEND, ' '))

    res_f.close()
    return

def print_monthly_averages():
    #### Calculate monhtly averages
    TRTM = data.TRTM
    TETPM = data.TETPM
    TETM = data.TETM
    NYR = 0
    for month in range(1,13):
        NYR = max(NYR, len(TRTM[month].keys()))
        #### Sum values for each year/month
        for year in TRTM[month].keys():
            TRTM[month][year] = sum(TRTM[month][year])
            TETPM[month][year] = sum(TETPM[month][year])
            TETM[month][year] = sum(TETM[month][year])
        #### For each month, average values for all years   
        TRTM[month] = sum(TRTM[month].values())/len(TRTM[month].values())
        TETPM[month] = sum(TETPM[month].values())/len(TETPM[month].values())
        TETM[month] = sum(TETM[month].values())/len(TETM[month].values())
    #### Sum of monthly averages
    RTNYR  = sum(TRTM.values())
    ETPNYR = sum(TETPM.values())
    ETNYR  = sum(TETM.values())

    res_f = open(paths.res, 'a')
    res_f.write('\n\n\n{0:<20s}{1}\n'.format(
        param.SITE, '     MONTHLY TOTALS AVERAGED OVER {0} YEARS\n'.format(NYR)))
    NET = '{0:>10s}'.format('NET ')
    if param.LK == 0:
        res_f.write('{0:30s}{1:>10s}{2:^20s}\n'.format(
            ' ', 'NET ', 'EVAPOTRANSPIRATION'))
        res_f.write('{0:20s}{1:>10s}{2:>10s}{3:>10s}{4:>10s}\n\n'.format(
            ' ','MONTH','RAD.','POTENT.','AREAL'))
    elif param.LK == 1:
        res_f.write('{0:30s}{1:>10s}{2:^20s}\n'.format(
            ' ', 'NET ', 'WET SURFACE EVAP.'))
        res_f.write('{0:20s}{1:>10s}{2:>10s}{3:>10s}{4:>10s}\n\n'.format(
            ' ','MONTH','RAD.','PAN-SIZE','LAKE-SIZE'))
    else:
        res_f.write('{0:30s}{1:>10s}{2:^20s}\n'.format(
            ' ', 'NET ', '  EVAPORATION'))
        res_f.write('{0:20s}{1:>10s}{2:>10s}{3:>10s}{4:>10s}\n\n'.format(
            ' ','MONTH','RAD.','POTENT.','LAKE'))

    for i in range(1,13):
        MONTH = '{0}'.format(calendar.month_abbr[i].upper())
        res_f.write('{0:>30s}{1:10.1f}{2:10.1f}{3:10.1f}\n'.format(
            MONTH, TRTM[i], TETPM[i], TETM[i]))
    res_f.write('{0:30s}{1:>10s}{1:>10s}{1:>10s}\n\n'.format(' ', '_'*8))
    res_f.write(
        '{0:>30s}{1:>10.1f}{2:>10.1f}{3:>10.1f}\n'.format(
            'TOTAL OF MONTHLY AVERAGES', RTNYR, ETPNYR, ETNYR))
    res_f.close()
    return
      
################################################################################
if __name__ == '__main__':
    workspace = os.getcwd()

    #### Create Display Logger ####
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.info('  WREVAP - Python')

    #### Get user input
    try:
        ini_name = sys.argv[1]
        data_name = sys.argv[2]
    except IndexError:
        logging.error('ERROR: WREVAP_v4 must be passed an input and data file')
        raise SystemExit()

    #### Check user_input
    #### INI and dataname can be just name (use CWD) or full path
    if not os.path.isfile(os.path.join(workspace, ini_name)):
        logging.error('\nERROR: The INI name is invalid\n')
        raise SystemExit()
    if not os.path.isfile(os.path.join(workspace, data_name)):
        logging.error('\nERROR: The data name is invalid\n')
        raise SystemExit()

    #### Run WREVAP
    logging.info('  Workspace:     {0}'.format(workspace))
    logging.info('  INI Filename:  {0}'.format(ini_name))
    logging.info('  Data Filename: {0}'.format(data_name))
    WREVAP(workspace, ini_name, data_name)
