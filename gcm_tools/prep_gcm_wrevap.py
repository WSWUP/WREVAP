#--------------------------------
# Name:         prep_gcm_wrevap.py
# Authors:      Charles Morton and Justin Huntington
# Modified      2016-11-02
# Python:       2.7
#--------------------------------

import calendar
import datetime as dt
import ConfigParser
import logging
import math
import os
import re
import sys

import numpy as np
from scipy import ndimage

def prep_gcm_wrevap(workspace):
    """Prepare GCM WREVAP monthly data files

    Calculate and save daily Thornton-Running solar
    Calculate and save daily Tdew
    Calculate Monthly Tmean, Tdew, and Rs
    Save monthly data to WREVAP data files
    Build WREVAP input files
    """

    logging.info('\nPreparing GCM WREVAP data')

    # Folder with csv data files

    # input_ws = r'Z:\USBR_Ag_Demands_Project\reservoir_evap\GCM_data\BC\CSV'
    # input_ws = 'C:\fieldProjects\KlamathBasin\data\openwaterevap\data\HDe'
    input_ws = r'C:\fieldProjects\KlamathBasin\data\openwaterevap\data\HDe'
    # input_ws = '../data\HDe'
    input_daily_temp = 'F'
    input_daily_prcp = 'IN'
    output_daily_temp = 'F'

    # Output folder

    # output_ws = r'Z:\USBR_Ag_Demands_Project\reservoir_evap\WREVAP_GCM'
    # output_ws = r'C:\fieldProjects\KlamathBasin\data\openwaterevap\data'
    output_ws = 'C:\\fieldProjects\KlamathBasin\data\openwaterevap\data'
    output_ws = '../data'
    output_monthly_temp = 'C'
    output_monthly_prcp = 'MM'
    res_list = []
    res_list = ['Clear_Lake']
    # res_list = ['Lake_Tahoe', 'Lake_Mead', 'Lake_Powell']
    #     'American_Falls', 'Boysen', 'Canyon_Ferry', 'Elephant_Butte',
    #     'Grand_Coulee', 'Lahontan', 'Millerton_Lake', 'Shasta_Lake',
    #     'Upper_Klamath']

    # Adjust precipitation
    # PPT will be multiplied by these ratios
    #  (i.e. decrease estimated values that are read)
    ppt_month_factor_dict = dict()
    """
    # New factors from Justin (2014-09-23)
    ppt_month_factor_dict['Lake_Tahoe'] = [
        0.683793809, 0.696432660, 0.691828784, 0.712677417,
        0.743181881, 0.668170098, 0.543929977, 0.604882144,
        0.607345296, 0.675140118, 0.678983427, 0.779749794]
    """

    # Adjust incoming solar
    # From: Lake_Tahoe_Solar_Calibration_Check.xlxs
    # Ratios are estimated divided by measured
    #  (i.e. decrease estimated values that are read in)
    rs_month_factor_dict = dict()
    """
    rs_month_factor_dict['Lake_Tahoe'] = [
        1.233, 1.144, 1.084, 1.025, 1.015, 1.041,
        1.065, 1.031, 1.029, 1.013, 1.132, 1.396]
    rs_month_factor_dict['Lahontan'] = [
        1.108210828, 1.086, 1.066844029, 1.05, 1.030296257, 1.011742403,
        1.013442609, 1.013532861, 1.01, 1.017420093, 1.072664479, 1.133997425]
    rs_month_factor_dict['Lake_Mead'] = [
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    rs_month_factor_dict['Lake_Powell'] = [
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    """

    # Output solar is in W/m^2 for dailies, and MJ/m^2/d for WREVAP inputs

    # Ancillary files

    res_param_name = 'reservoir_parameters.csv'
    Ko_monthly_temp = 'C'
    default_ini_name = 'default.ini'

    # header_name = 'HEADER.txt'
    # bad_data_name = 'bad_data.csv'

    # Reservoir parameter fields

    # gcm_name_field = 'GCM_NAME'
    input_name_field = 'GCM_NAME'
    output_name_field = 'OUTPUT_NAME'
    rdm_field = 'RDM'
    rdu_field = 'RDU'
    depth_field = 'DEPTH'
    salin_field = 'SALINITY'
    elev_field = 'ELEVATION'
    lat_field = 'LATITUDE'
    # lon_field = 'LONGITUDE'
    tr_b0_field = 'TR_B0'
    tr_b1_field = 'TR_B1'
    tr_b2_field = 'TR_B2'
    ko_field = 'KO_'

    # Overwrite existing files

    overwrite_flag = True

    # DeltaT year range

    delta_t_year_range = (1950, 1999)

    # Check inputs

    if not os.path.isdir(input_ws):
        logging.error(
            '\nERROR: Input workspace {} does not exist'.format(input_ws))
        raise SystemExit
    if not os.path.isdir(output_ws):
        os.mkdir(output_ws)
    res_param_path = os.path.join(workspace, res_param_name)
    if not os.path.isfile(res_param_path):
        logging.error((
            '\nERROR: Reservoir parameter file (\'{}\')' +
            ' does not exist').format(res_param_path))
        raise SystemExit()
    default_ini_path = os.path.join(workspace, default_ini_name)
    if not os.path.isfile(default_ini_path):
        logging.error((
            '\nERROR: Default config file (\'{}\')' +
            ' does not exist').format(default_ini_path))
        raise SystemExit()
        
    # find some values that can come from default.ini or reservoir_parameters.csv
        
    config = ConfigParser.ConfigParser()
    try:
        config.readfp(open(default_ini_path))
        """
        # configuration debug
        cfg_path = os.path.join(os.getcwd(), "default.cfg")
        with open(cfg_path, 'wb') as cf: config.write(cf)
        """
        config.has_section('INPUTS')
        config_items = config.items('INPUTS')
    except ConfigParser.NoSectionError:
        logging.error((
            '\nERROR: {}\n' +
            '  Paramter INI file is missing a section line\n' +
            '  First data line in file needs to be: [INPUTS]\n' +
            '  Try removing INI file and rebuilding it\n').format(
                default_ini_path))
        sys.exit()
    except ConfigParser.Error:
        logging.error((
            '\nERROR: \n' + str(sys.exc_info()[0]) +
            '\noccurred in configuration (INI) file \n{}' +
            'Try removing INI file and rebuilding it\n').format(
            default_ini_path))
        sys.exit()
    config_keys = [i.upper() for i in dict(config_items).keys()]
    if "RDM" in config_keys:
        crdm = read_param('RDM', 0, int, config)
    else:
        crdm = 0
    if "RDU" in config_keys:
        crdu = read_param('RDU', 0, int, config)
    else:
        crdu = 0
    input_daily_temp = input_daily_temp.upper()
    input_daily_prcp = input_daily_prcp.upper()
    output_daily_temp = output_daily_temp.upper()
    output_monthly_temp = output_monthly_temp.upper()
    Ko_monthly_temp = Ko_monthly_temp.upper()
    if input_daily_temp not in ['F', 'C']:
        logging.error(
            '\nERROR: Input daily temperature must be F or C\n')
        raise SystemExit()
    if input_daily_prcp not in ['MM', 'IN']:
        logging.error(
            '\nERROR: Input daily prcipitation must be MM or IN\n')
        raise SystemExit()
    if output_daily_temp not in ['F', 'C']:
        logging.error(
            '\nERROR: Output daily temperature must be F or C\n')
        raise SystemExit()
    if output_monthly_temp not in ['F', 'C']:
        logging.error(
            '\nERROR: Output monthly temperature must be F or C\n')
        raise SystemExit()
    if output_monthly_prcp not in ['MM', 'IN']:
        logging.error(
            '\nERROR: Output monthly prcipitation must be MM or IN\n')
        raise SystemExit()
    if Ko_monthly_temp not in ['F', 'C']:
        logging.error(
            '\nERROR: Ko temperature must be F or C\n')
        raise SystemExit()

    # Input file regular expressions

    csv_re = re.compile('\w+_\w+.csv$')
    tmax_csv_re = re.compile('\w+_TMAX.csv$')
    tmin_csv_re = re.compile('\w+_TMIN.csv$')
    # dtrg_csv_re = re.compile('\w+_DTRG.csv$')
    prcp_csv_re = re.compile('\w+_PRCP.csv$')
    # tdew_csv_re = re.compile('\w+_TDEW.csv$')
    # rso_csv_re = re.compile('\w+_RSO.csv$')
    # rs_csv_re = re.compile('\w+_RS.csv$')

    # Build list of reservoirs with TMAX and TMIN data

    csv_list = [item for item in os.listdir(input_ws) if csv_re.match(item)]
    tmax_csv_dict = dict([
        (item.split('_')[0], item) for item in csv_list
        if tmax_csv_re.match(item)])
    tmin_csv_dict = dict([
        (item.split('_')[0], item) for item in csv_list
        if tmin_csv_re.match(item)])
    prcp_csv_dict = dict([
        (item.split('_')[0], item) for item in csv_list
        if prcp_csv_re.match(item)])
    # tdew_csv_dict = dict([
    #     (item.split('_')[0], item) for item in csv_list
    #     if tdew_csv_re.match(item)])
    # output_res_list = sorted(list(
    #     set(tmax_csv_dict.keys() + tmin_csv_dict.keys())))

    # Read reservoir parameters as structured array
    
    res_param_array = np.atleast_2d(np.genfromtxt(
        res_param_path, delimiter=',', names = True, dtype = None))
    columns = res_param_array.dtype.names
    # res_param_fields = res_param_array.dtype.names
    # DEADBEEF
    # res_param_names = res_param_array[gcm_name_field]
    input_res_array = res_param_array[input_name_field]
    output_res_array = res_param_array[output_name_field]

    # Use reservoir parameters file to set input and output name lists

    input_res_list = input_res_array[0]
    output_res_list = output_res_array[0]

    # Write bad data to a txt file
    # bad_data_f = open(bad_data_path, 'w')

    # Check keys on factors

    for output_res in ppt_month_factor_dict.keys():
        if output_res not in output_res_list:
            logging.error(
                'Reservoir name {} is invalid'.format(output_res))
            raise SystemExit()
    for output_res in rs_month_factor_dict.keys():
        if output_res not in output_res_list:
            logging.error(
                'Reservoir name {} is invalid'.format(output_res))
            raise SystemExit()

    # Process each reservoir

    for input_res, output_res in sorted(zip(input_res_list, output_res_list)):
        logging.info('  Evaporation Name: {}'.format(output_res))

        # Process target reservoirs

        if res_list and output_res not in res_list:
            continue

        # A TMAX and TMIN csv must be present to run

        try:
            tmax_csv_path = os.path.join(input_ws, tmax_csv_dict[input_res])
            tmin_csv_path = os.path.join(input_ws, tmin_csv_dict[input_res])
            prcp_csv_path = os.path.join(input_ws, prcp_csv_dict[input_res])
        except KeyError:
            continue
        logging.info('  GCM Name: {}'.format(input_res))

        # Daily output file names

        tdew_csv_name = '{}_TDEW.csv'.format(input_res)
        rso_csv_name = '{}_RSO.csv'.format(input_res)
        rs_csv_name = '{}_RS.csv'.format(input_res)

        # Daily output file paths

        tdew_csv_path = os.path.join(input_ws, tdew_csv_name)
        rso_csv_path = os.path.join(input_ws, rso_csv_name)
        rs_csv_path = os.path.join(input_ws, rs_csv_name)

        # Remove existing daily files if necessary

        if overwrite_flag and os.path.isfile(tdew_csv_path):
            os.remove(tdew_csv_path)
        if overwrite_flag and os.path.isfile(rso_csv_path):
            os.remove(rso_csv_path)
        if overwrite_flag and os.path.isfile(rs_csv_path):
            os.remove(rs_csv_path)

        # Get reservoir parameter data from file using output_res as key

        if output_res not in output_res_array:
            logging.info('    Reservoir does not have data in {}'.format(
                output_res))
            continue
        res_param = res_param_array[output_res_array == output_res]
        # output_name = res_param[output_name_field][0]
        depth_flt = float(res_param[depth_field][0])
        salin_flt = float(res_param[salin_field][0])
        elev_flt = float(res_param[elev_field][0])
        lat_flt = float(res_param[lat_field][0])
        # lon_flt = float(res_param[lon_field][0])
        tr_b0_flt = float(res_param[tr_b0_field][0])
        tr_b1_flt = float(res_param[tr_b1_field][0])
        tr_b2_flt = float(res_param[tr_b2_field][0])
        if rdm_field in columns:
            rdm = int(res_param[rdm_field][0])
        else:
            rdm = crdm
        if rdu_field in columns:
            rdu = int(res_param[rdu_field][0])
        else:
            rdu = crdu
        ko_month_array = np.array([
            float(res_param[ko_field + month_abbr.upper()][0])
            for month_i, month_abbr in enumerate(calendar.month_abbr)
            if month_abbr])
        del res_param

        # Convert ko temperature to celsius

        if Ko_monthly_temp == 'F':
            ko_month_array = fahrenheit_to_celsius(ko_month_array)
        elif Ko_monthly_temp == 'C':
            pass
        else:
            continue

        # Monthly output workspace

        res_output_ws = os.path.join(output_ws, output_res)
        if not os.path.isdir(res_output_ws):
            os.mkdir(res_output_ws)
        res_data_ws = os.path.join(res_output_ws, 'input_data')
        if not os.path.isdir(res_data_ws):
            os.mkdir(res_data_ws)

        # Read header row for each file to get GCM list
        # First field is date field, skip

        with open(tmax_csv_path, 'r') as f:
            header_str = f.readline().strip()
            # Header str has a comma at end
            if header_str[-1] == ',':
                header_str = header_str[:-1]
        f.closed
        gcm_name_list = [i.strip() for i in header_str.split(',') if i][1:]
        header_fmt_str = '%d,' + ','.join(['%7.4f'] * len(gcm_name_list))

        # Read daily Tmin and Tmax data, skip header row

        tmax_csv_array = np.loadtxt(
            tmax_csv_path, delimiter=',', skiprows=1)
        tmin_csv_array = np.loadtxt(
            tmin_csv_path, delimiter=',', skiprows=1)
        prcp_csv_array = np.loadtxt(
            prcp_csv_path, delimiter=',', skiprows=1)
        del tmax_csv_path, tmin_csv_path, prcp_csv_path

        # Get data subsets from full csv array

        date_array = tmax_csv_array[:, 0].astype(np.int)
        tmax_array = tmax_csv_array[:, 1:]
        tmin_array = tmin_csv_array[:, 1:]
        prcp_array = prcp_csv_array[:, 1:]
        del tmax_csv_array, tmin_csv_array, prcp_csv_array

        # Convert input temperature to celsius (for Thornton-Running)
        
        if input_daily_temp == 'F':
            tmax_array = fahrenheit_to_celsius(tmax_array)
            tmin_array = fahrenheit_to_celsius(tmin_array)
        elif input_daily_temp == 'C':
            pass
        else:
            continue

        # Convert input precipitation to mm (from inches)

        if input_daily_prcp == 'IN':
            prcp_array *= 25.4
        elif input_daily_prcp == 'MM':
            pass
        else:
            continue

        # Calculate Tmax - Tmin

        tmax_tmin_array = tmax_array - tmin_array

        # Where Tmax > Tmin, set to nan and save

        tmax_tmin_mask = (tmax_tmin_array >= 0)
        if not np.all(tmax_tmin_mask):
            logging.warning(
                '    WARNING: TMIN > TMAX for {} cells'.format(
                    np.sum(~tmax_tmin_mask)))
            tmax_tmin_mask_i = np.nonzero(~tmax_tmin_mask)
            for gcm, date in zip(
                    np.array(gcm_name_list)[tmax_tmin_mask_i[1]],
                    date_array[tmax_tmin_mask_i[0]]):
                # bad_data_f.write(
                #     ','.join([output_res, gcm, str(date)]) + '\n')
                logging.debug('    {} {}'.format(gcm, date))
                del date, gcm
            del tmax_tmin_mask_i
        tmax_tmin_array[~tmax_tmin_mask] = np.nan
        tmax_array[~tmax_tmin_mask] = np.nan
        tmin_array[~tmax_tmin_mask] = np.nan

        # Calculate mean temperature

        tmean_array = 0.5 * (tmax_array + tmin_array)
        del tmax_array

        # Get Year/Month for each DOY

        ym_daily_array = np.array([
            dt.datetime.strptime(str(d), '%Y%m%d').strftime('%Y_%m') for
            d in date_array])
        year_array = np.array([
            dt.datetime.strptime(str(d), '%Y%m%d').strftime('%Y') for
            d in date_array]).astype(np.int)
        month_array = np.array([
            dt.datetime.strptime(str(d), '%Y%m%d').strftime('%m') for
            d in date_array]).astype(np.int)

        # Get DOY from date

        doy_array = np.array([
            dt.datetime.strptime(str(d), '%Y%m%d').strftime('%j') for
            d in date_array]).astype(np.int)

        # Get Year/Month for each month

        ym_month_array = np.unique(ym_daily_array)

        # Calculate mean monthly temperature difference for each day

        delta_t_array = np.zeros(tmean_array.shape)

        # Only use "historical" temperatures

        year_mask = ((year_array >= delta_t_year_range[0]) &
                     (year_array <= delta_t_year_range[1]))
        for gcm_i, gcm_name in enumerate(gcm_name_list):

            # Only calculate mean for days with Tmax > Tmin

            gcm_mask = (year_mask & tmax_tmin_mask[:, gcm_i])
            mean_array = ndimage.mean(
                tmax_tmin_array[:, gcm_i][gcm_mask],
                labels=month_array[gcm_mask], index=range(1, 13))
            for month_i, month in enumerate(range(1, 13)):
                month_mask = (month_array == month)
                delta_t_array[:, gcm_i][month_mask] = mean_array[month_i]
                del month_mask
            del mean_array, gcm_mask

        # Map KO from mean monthlies to individual dailies

        ko_array = np.zeros(tmin_array.shape)
        for month_i, month in enumerate(range(1, 13)):
            ko_array[month_array == month] = ko_month_array[month_i]

        # Calculate Tdew

        tdew_array = (tmin_array - ko_array)
        del ko_array, tmin_array

        # Calculate vapor pressure

        ea_array = calc_vapor_pressure(tdew_array)

        # Calculate pressure

        pres_flt = float(calc_pressure(elev_flt))

        # Calculate Thornton-Running Solar

        rso_array, rs_array = calc_tr_rso(
            lat_flt, pres_flt, ea_array, doy_array[:, np.newaxis],
            delta_t_array, tmax_tmin_array,
            tr_b0_flt, tr_b1_flt, tr_b2_flt)
        del ea_array, pres_flt, delta_t_array, tmax_tmin_array
        del doy_array

        # Scale Rso values by monthly factor

        if output_res in rs_month_factor_dict.keys():
            for month_i, month in enumerate(range(1, 13)):
                logging.info('    Scaling Monthly Rs: {} {}'.format(
                    month, rs_month_factor_dict[output_res][month_i]))
                rs_array[month_array == month] /= rs_month_factor_dict[
                    output_res][month_i]

        # Make copies of output arrays for unit conversions

        tdew_daily_array = np.copy(tdew_array)
        rso_daily_array = np.copy(rso_array)
        rs_daily_array = np.copy(rs_array)

        # Convert output daily dew point temperature from celsius

        if output_daily_temp == 'F':
            tdew_daily_array = celsius_to_fahrenheit(tdew_daily_array)
        elif output_daily_temp == 'C':
            pass
        else:
            continue

        # Convert daily solar from MJ/m2/d to w/m2

        rso_daily_array *= 11.574
        rs_daily_array *= 11.574

        # Save daily values
        # Append date_array at beginning

        if not os.path.isfile(tdew_csv_path):
            np.savetxt(
                tdew_csv_path,
                np.hstack((date_array[:, np.newaxis], tdew_daily_array)),
                delimiter = ',', header = header_str,
                fmt = header_fmt_str, comments = '')
        if not os.path.isfile(rso_csv_path):
            np.savetxt(
                rso_csv_path,
                np.hstack((date_array[:, np.newaxis], rso_daily_array)),
                delimiter=',', header=header_str,
                fmt=header_fmt_str, comments='')
        if not os.path.isfile(rs_csv_path):
            np.savetxt(
                rs_csv_path,
                np.hstack((date_array[:, np.newaxis], rs_daily_array)),
                delimiter=',', header=header_str,
                fmt=header_fmt_str, comments='')
        del tdew_daily_array, rso_daily_array, rs_daily_array,
        del date_array, header_fmt_str

        # Calculate fixed monthly values

        year_array = np.array(
            [int(ym.split('_')[0]) for ym in ym_month_array]).astype(np.int)
        month_array = np.array(
            [int(ym.split('_')[1]) for ym in ym_month_array]).astype(np.int)
        start_array = np.ones(ym_month_array.size).astype(np.int)
        length_array = np.array([
            calendar.monthrange(*map(int, ym.split('_')))[1]
            # calendar.monthrange(int(ym.split('_')[0]), int(ym.split('_')[1]))[1]
            for ym in ym_month_array]).astype(np.int)

        # Calculate variable monthly values (Rso, Rs, Tdew, Tmean)

        tmean_month_array = np.zeros((ym_month_array.size, tmean_array.shape[1]))
        tdew_month_array = np.zeros((ym_month_array.size, tmean_array.shape[1]))
        rso_month_array = np.zeros((ym_month_array.size, tmean_array.shape[1]))
        rs_month_array = np.zeros((ym_month_array.size, tmean_array.shape[1]))
        prcp_month_array = np.zeros((ym_month_array.size, prcp_array.shape[1]))
        for gcm_i, gcm_name in enumerate(gcm_name_list):
            tmean_month_array[:, gcm_i] = ndimage.mean(
                tmean_array[:, gcm_i][tmax_tmin_mask[:, gcm_i]],
                labels=ym_daily_array[tmax_tmin_mask[:, gcm_i]],
                index=ym_month_array)
            tdew_month_array[:, gcm_i] = ndimage.mean(
                tdew_array[:, gcm_i][tmax_tmin_mask[:, gcm_i]],
                labels=ym_daily_array[tmax_tmin_mask[:, gcm_i]],
                index=ym_month_array)
            rso_month_array[:, gcm_i] = ndimage.mean(
                rso_array[:, gcm_i][tmax_tmin_mask[:, gcm_i]],
                labels=ym_daily_array[tmax_tmin_mask[:, gcm_i]],
                index=ym_month_array)
            rs_month_array[:, gcm_i] = ndimage.mean(
                rs_array[:, gcm_i][tmax_tmin_mask[:, gcm_i]],
                labels=ym_daily_array[tmax_tmin_mask[:, gcm_i]],
                index=ym_month_array)
            prcp_month_array[:, gcm_i] = ndimage.sum(
                prcp_array[:, gcm_i][tmax_tmin_mask[:, gcm_i]],
                labels=ym_daily_array[tmax_tmin_mask[:, gcm_i]],
                index=ym_month_array)
        del tmean_array, tdew_array, rso_array, rs_array, prcp_array
        del ym_month_array, ym_daily_array

        # Convert output monthly dew point temperature from celsius

        if output_monthly_temp == 'F':
            tdew_month_array = celsius_to_fahrenheit(tdew_month_array)
            tmean_month_array = celsius_to_fahrenheit(tmean_month_array)
        elif output_monthly_temp == 'C':
            pass
        else:
            continue

        # Scale PPT values by monthly factor

        if output_res in ppt_month_factor_dict.keys():
            for month_i, month in enumerate(range(1, 13)):
                logging.info('    Scaling Monthly PPT: {} {}'.format(
                    month, ppt_month_factor_dict[output_res][month_i]))
                prcp_month_array[month_array == month] *= ppt_month_factor_dict[
                    output_res][month_i]

        # Convert output precipitation to inches if requested

        if output_monthly_prcp == 'IN':
            prcp_month_array /= 25.4
        elif output_monthly_prcp == 'MM':
            pass
        else:
            continue

        # Save monthly data for each GCM to a separate data file

        wrevap_header_str = 'YEAR,MONTH,STARTDAY,LENGTH,TD,T,S,PPT'
        wrevap_fmt_str = '%d,%d,%d,%d,%f,%f,%f,%f'
        for gcm_i, gcm_name in enumerate(gcm_name_list):
            output_csv_name = '{}_{}.csv'.format(
                output_res, gcm_name.replace('.', '_'))
            output_csv_path = os.path.join(res_data_ws, output_csv_name)
            if overwrite_flag and os.path.isfile(output_csv_path):
                os.remove(output_csv_path)
            if not os.path.isfile(output_csv_path):
                np.savetxt(
                    output_csv_path,
                    np.vstack((
                        year_array, month_array, start_array, length_array,
                        tdew_month_array[:, gcm_i],
                        tmean_month_array[:, gcm_i],
                        rs_month_array[:, gcm_i],
                        prcp_month_array[:, gcm_i])).T,
                    delimiter=',', header=wrevap_header_str,
                    comments='', fmt=wrevap_fmt_str)
        del year_array, month_array, start_array, length_array
        del tdew_month_array, tmean_month_array
        del rs_month_array, rso_month_array, prcp_month_array

        # Build a Python WREVAP_gcm input file

        output_ini_name = '{}.ini'.format(output_res.lower())
        output_ini_path = os.path.join(res_output_ws, output_ini_name)
        if os.path.isfile(output_ini_path):
            os.remove(output_ini_path)
        default_ini_f = open(default_ini_path, 'r')
        output_ini_f = open(output_ini_path, 'w')
        for default_line in default_ini_f:
            if default_line.strip().startswith('SITE = '):
                default_line = 'SITE = {}\n'.format(
                    output_res.upper().replace(' ', '_'))
            elif default_line.strip().startswith('site_fn = '): continue
            elif default_line.strip().startswith('PHID = '):
                default_line = 'PHID = {}\n'.format(lat_flt)
            elif default_line.strip().startswith('P = '):
                default_line = 'P = {}\n'.format(elev_flt)
            elif default_line.strip().startswith('RDM = '): continue
            elif default_line.strip().startswith('RDU = '): continue
            elif default_line.strip().startswith('DA = '):
                default_line = 'RDM = {}\n'.format(rdm)
                output_ini_f.write(default_line)
                default_line = 'RDU = {}\n'.format(rdu)
                output_ini_f.write(default_line)
                default_line = 'DA = {}\n'.format(depth_flt)
            elif default_line.strip().startswith('SALT = '):
                default_line = 'SALT = {}\n'.format(salin_flt)
            elif default_line.strip().startswith('LK = '):
                default_line = 'LK = 2\n'
            elif default_line.strip().startswith('IT = '):
                if output_monthly_temp == 'C':
                    default_line = 'IT = 0\n'
                elif output_monthly_temp == 'F':
                    default_line = 'IT = 1\n'
            elif default_line.strip().startswith('IS = '):
                default_line = 'IS = 3\n'
            elif default_line.strip().startswith('IV = '):
                default_line = 'IV = 0\n'
            elif default_line.strip().startswith('IP = '):
                default_line = 'IP = 1\n'
            output_ini_f.write(default_line)
        default_ini_f.close()
        output_ini_f.close()

        # Cleanup

        del header_str, gcm_name_list
        del depth_flt, salin_flt, elev_flt, lat_flt
    del input_res_list, output_res_list

def fahrenheit_to_celsius(temp_array):
    temp_array -= 32
    temp_array *= (5. / 9)
    return temp_array

def celsius_to_fahrenheit(temp_array):
    temp_array *= (9. / 5)
    temp_array += 32
    return temp_array

def calc_pressure(elev):
    """From ASCE-EWRI 2005"""
    return 101.3 * np.power(((293. - 0.0065 * elev) / 293.), 5.26)
    # Wikipedia International Standard Atmosphere
    # return 101.325 * np.power(((288.15 - 0.0065 * elev) / 288.15), 5.2558)

def calc_vapor_pressure(tdew):
    """From METRIC code"""
    return 0.6108 * np.exp(17.27 * tdew / (tdew + 237.3))

def calc_tr_rso(lat, P, ea, J, deltaT, tmax_tmin,
                b0=0.031, b1=0.201, b2=-0.185):
    """Empirical Thornton-Running Rso (modified from Matlab code)

    keyword arguments:
    lat -- given lat at station (decimal deg)
    P -- pressure at station (in kPa)
    ea -- actual vapor pressure (in kPa)
    J -- day of year
    deltaT --
    tmax_tmin --

    """
    # LATITUDE IN RADIANS

    lat_rad = lat * math.pi / 180.

    # Solar constant [MJ m-2 h-1]
    # Changed to 4.9212 to be same as ref et modules's 0.08202 * 60
    # Gsc = 4.92
    Gsc = 4.9212

    # epn 24

    sinb24 = np.sin(
        0.85 + 0.3 * lat_rad * np.sin((2 * math.pi / 365) * J - 1.39) -
        0.42 * lat_rad ** 2)
    # eqn 22 PRECIPITABLE WATER
    W = 0.14 * ea * P + 2.1

    # atmospheric clearness coeff. ranges from 0-1. 1 for clean air

    Ktb = 0.9

    # eqn 21 index of atmospheric clearness for direct beam radiation

    kb = 0.98 * np.exp(
        ((-0.00146 * P) / (Ktb * sinb24)) - 0.075 * (W / sinb24) ** 0.4)

    # eqn 51 SOLAR DECLINATION ANGLE (radians)

    sda = 0.409 * np.sin(2 * math.pi * J / 365 - 1.39)

    # eqn 23 INDEX OF TRANSMISSIVITY FOR DIFFUSE RADIATION

    kd = np.copy(kb)
    kb_mask = kb > 0.15
    kd[kb_mask] = (0.35 - 0.36 * kd[kb_mask])
    kd[~kb_mask] = (0.18 + 0.82 * kd[~kb_mask])

    # eqn 20 ATMOSPHERIC TRANSMISSIVITY

    kt = kb + kd
    ws = np.arccos(-1 * math.tan(lat_rad) * np.tan(sda))

    # eqn 18 -- squared inverse relative distance factor

    dr = 1 + 0.033 * np.cos((2 * math.pi / 365) * J)

    # eqn 17 EXOATMOSPHERIC RADIATION

    Ra = 24. / math.pi * Gsc * dr * (
        ws * math.sin(lat_rad) * np.sin(sda) +
        math.cos(lat_rad) * np.cos(sda) * np.sin(ws))

    # eqn 16  THEORETICAL SOLAR RADIATION

    Rso = kt * Ra

    # eqn 15  EMPIRICAL FITTING COEFFICIENT %From ALLEN an Robison

    # B = 0.023 + 0.1 * math.exp(-0.2 * deltaT);
    # eqn 15  EMPIRICAL FITTING COEFFICIENT %TR
    # B = 0.031 + 0.201 * np.exp(-0.185 * deltaT)
    B = b0 + b1 * np.exp(b2 * deltaT)

    # eqn 14 EMPIRICAL SOLAR RADIATION

    Rs = Rso * (1 - 0.9 * np.exp(-1 * B * np.power(tmax_tmin, 1.5)))
    return Rso, Rs

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

if __name__ == '__main__':
    workspace = os.getcwd()
    logging.basicConfig(level = logging.INFO, format = '%(message)s')
    logging.info('\n{}'.format('#' * 80))
    log_f = '{:<20s} {}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    prep_gcm_wrevap(workspace)
