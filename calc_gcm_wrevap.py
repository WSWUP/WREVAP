#--------------------------------
# Name:         calc_gcm_wrevap.py
# Authors:      Charles Morton and Justin Huntington
# Modified:     2016-11-02
# Python:       2.7
#--------------------------------

import calendar
import datetime as dt
import logging
import os
import shutil
import sys

# import numpy as np


def calc_gcm_wrevap(workspace):
    """Calculate WREVAP for each reservoir and GCM input data set"""
    logging.info('Calculating GCM WREVAP')

    # Output folder
    wrevap_ws = r'Z:\USBR_Ag_Demands_Project\reservoir_evap\WREVAP_GCM'

    input_folder_name = 'input_data'
    output_folder_name = 'output_data'

    # Input parameters
    lk_model_str = 'LK = '
    wrevap_name = 'WREVAP_gcm.py'

    res_list = []
    # res_list = ['Lake_Mead']
    res_list = ['Lake_Tahoe', 'Lake_Mead', 'Lake_Powell']
    # res_list = [
    #     'American_Falls', 'Boysen', 'Canyon_Ferry', 'Elephant_Butte',
    #     'Grand_Coulee', 'Lahontan', 'Millerton_Lake', 'Shasta_Lake',
    #     'Upper_Klamath']

    # Ancillary files
    # res_param_name = 'reservoir_parameters.csv'

    # Reservoir parameter fields
    # gcm_name_field = 'GCM_NAME'
    # output_name_field = 'OUTPUT_NAME'

    # Overwrite existing files
    overwrite_flag = True

    # Adjust ET rates to account for advective heat
    # From: advection_addition_for_lake_tahoe.xlsx
    # This amount of ET (mm) will be added to ET_LAKE_MM values
    et_lake_month_adjust_dict = dict()
    # Tahoe rates are intentionally being subtracted (typically negative advected heat flux)
    et_lake_month_adjust_dict['Lake_Tahoe'] = [
        -16.581008, -14.740917, -10.949723, -4.945036, 1.720671, -0.201838,
        -3.128492, -5.005587, -4.137683, -2.972068, -1.463327, -5.109871]
    # Mead and Powell (on average net advected heat flux is positive)
    # These are the approximated factors that are much smoother from month to month
    et_lake_month_adjust_dict['Lake_Mead'] = [
        17.74612501, 0.974963521, -15.52785938, -7.513480347, 13.30959376, 31.12727573,
        57.67490629, 75.4210313, 51.5210081, 33.2739844, 23.61379538, 22.18265626]
    et_lake_month_adjust_dict['Lake_Powell'] = [
        17.74612501, 0.974963521, -15.52785938, -7.513480347, 13.30959376, 31.12727573,
        57.67490629, 75.4210313, 51.5210081, 33.2739844, 23.61379538, 22.18265626]

    # # These were the original factors from Moreo (average of 2 years)
    # # There is a lot of variation from month to month
    # et_lake_month_adjust_dict['Lake_Mead'] = [
    #     14.3548819, 0.788650266, -12.56052166, -6.077671771, 10.76616142, 25.17892591,
    #     56.52234747, 86.12929139, 22.57420944, 35.88720474, 19.10125414, 51.13926676]
    # et_lake_month_adjust_dict['Lake_Powell'] = [
    #     14.3548819, 0.788650266, -12.56052166, -6.077671771, 10.76616142, 25.17892591,
    #     56.52234747, 86.12929139, 22.57420944, 35.88720474, 19.10125414, 51.13926676]

    # Adjust evaporation values for times when reservoir is frozen
    evap_month_factor_dict = dict()
    # evap_month_factor_dict['Boca'] = [
    #     0.82, 0.90, 0.39, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25]


    # Check inputs
    if not os.path.isdir(wrevap_ws):
        logging.error(
            ('\nERROR: The WREVAP workspace {} does ' +
             'not exist\n'.format(wrevap_ws)))
        raise SystemExit()
    # res_param_path = os.path.join(workspace, res_param_name)
    # if not os.path.isfile(res_param_path):
    #     logging.error(
    #         ('\nERROR: The reservoir parameter file (\'{}\')' +
    #          ' does not exist\n').format(res_param_path))
    #     raise SystemExit()

    # Read in reservoir parameters as structured array
    # res_param_array = np.genfromtxt(
    #     res_param_path, delimiter=',', names=True, dtype=None)
    # res_param_fields = res_param_array.dtype.names
    # res_name_list = res_param_array[output_name_field]

    # Try to process each sub folder
    res_name_list = os.listdir(wrevap_ws)

    # Process each reservoir
    for res_name in res_name_list:
        # Reservoir folder must already exist
        if not os.path.isdir(res_name):
            continue
        # Process target reservoirs
        if res_list and res_name not in res_list:
            continue
        logging.info('{}'.format(res_name))
        res_ws = os.path.join(wrevap_ws, res_name)

        # Check if there is an input file
        ini_name = res_name.lower() + '.ini'
        ini_path = os.path.join(res_ws, ini_name)
        if not os.path.isfile(ini_path):
            logging.info(
                '  The input file {} is not present'.format(ini_name))
            continue

        # Make copies of input file with model option 2 and 3 set
        ini_option2_name = ini_name.replace('.ini', '_option2.ini')
        ini_option3_name = ini_name.replace('.ini', '_option3.ini')
        ini_option2_path = os.path.join(res_ws, ini_option2_name)
        ini_option3_path = os.path.join(res_ws, ini_option3_name)
        ini_option2_f = open(ini_option2_path, 'w')
        ini_option3_f = open(ini_option3_path, 'w')
        with open(ini_path, 'r') as input_f:
            for line in input_f:
                if line.strip().startswith(lk_model_str):
                    ini_option2_f.write(lk_model_str + '2\n')
                    ini_option3_f.write(lk_model_str + '3\n')
                else:
                    ini_option2_f.write(line)
                    ini_option3_f.write(line)
                del line
        ini_option2_f.close()
        ini_option3_f.close()

        # Check if there are input data files
        input_folder_ws = os.path.join(res_ws, input_folder_name)
        if not os.path.isdir(input_folder_ws):
            logging.info('  The input data folder {} is not present'.format(
                input_folder_name))
            continue

        # Build output data folder
        output_folder_ws = os.path.join(res_ws, output_folder_name)
        if os.path.isdir(output_folder_ws) and overwrite_flag:
            logging.info(
                '  Removing the existing output data folder {}'.format(
                    output_folder_name))
            shutil.rmtree(output_folder_ws)
            # DEADBEEF - I needed the code to pause here for a minute
            #   so our network drive could catch up
            raw_input('Press ENTER to continue')
        if not os.path.isdir(output_folder_ws):
            os.mkdir(output_folder_ws)

        # Copy WREVAP code into folder if necessary
        wrevap_path = os.path.join(wrevap_ws, wrevap_name)
        wrevap_res_path = os.path.join(res_ws, wrevap_name)
        if (os.path.isfile(wrevap_path) and
                (not os.path.isfile(wrevap_res_path) or overwrite_flag)):
            shutil.copy(wrevap_path, res_ws)
        if not os.path.isfile(wrevap_res_path):
            logging.info('  {} was not present'.format(wrevap_name))
            continue

        # Set working directory to reservoir folder for os.system calls
        os.chdir(res_ws)

        # Iteratively get data files
        for data_name in os.listdir(input_folder_ws):
            data_path = os.path.join(input_folder_ws, data_name)
            if not os.path.isfile(data_path):
                continue
            if not data_name.endswith('.csv'):
                continue
            sol_path = os.path.join(res_ws, data_name.replace('.csv', '.SOL'))
            tgw_path = os.path.join(res_ws, data_name.replace('.csv', '.TGW'))
            res_path = os.path.join(res_ws, data_name.replace('.csv', '.RES'))
            out_path = os.path.join(res_ws, data_name.replace('.csv', '.OUT'))
            csv_path = os.path.join(output_folder_ws, data_name)

            # Check if output file already eixsts
            # if os.path.isfile(csv_path) and not overwrite_flag:
            #     continue
            # elif os.path.isfile(csv_path) and overwrite_flag:
            #     try:
            #         os.remove(csv_path)
            #     except:
            #         continue

            logging.info('  {}'.format(
                data_name.split('.')[0].split('_', 1)[1]))

            # Copy data file to reservoir folder
            shutil.copy(data_path, res_ws)

            # Set model option to 2
            # for line in fileinput.input(ini_path, inplace=True):
            #     if line.strip().startswith(lk_model_str):
            #         logging.info(lk_model_str+'2\n')
            #     else:
            #         logging.info('{}'.format(line))
            #     del line

            # Run WREVAP
            os.system('{0} {1} {2}'.format(
                wrevap_name, ini_option2_name, data_name))
            # os.system('{0} {1} {2}'.format(wrevap_name, ini_name, data_name))
            # subprocess.call([wrevap_name, ini_name, data_name])

            # Set model option to 3
            # for line in fileinput.input(ini_path, inplace=True):
            #     if line.strip().startswith(lk_model_str):
            #         logging.info(lk_model_str+'3\n')
            #     else:
            #         logging.info('{}'.format(line))
            #     del line

            # Copy SOL to TGW
            shutil.copy(sol_path, tgw_path)

            # Run WREVAP
            os.system('{0} {1} {2}'.format(
                wrevap_name, ini_option3_name, data_name))
            # os.system('{0} {1} {2}'.format(wrevap_name, ini_name, data_name))
            # subprocess.call([wrevap_name, ini_name, data_name])

            # Copy OUT (modified RES as csv output) file to output_data folder
            # shutil.copy(out_path, output_folder_ws)
            shutil.copy(out_path, csv_path)

            # Cleanup
            os.remove(os.path.join(res_ws, data_name))
            os.remove(sol_path)
            os.remove(tgw_path)
            os.remove(res_path)
            os.remove(out_path)
            del data_path, sol_path, tgw_path, res_path, out_path

        # Remove ancillary files
        os.remove(ini_option2_path)
        os.remove(ini_option3_path)

    # Reset working directory
    os.chdir(wrevap_ws)

    # Scale evaporation values if necessary
    logging.info(
        '\nAdjusting evaporation values for advection' +
        '\nAdjusting evaporation values for frozen lakes')
    month_num_dict = dict([
        (month_abbr.upper(), (month_i)) for
        month_i, month_abbr in enumerate(calendar.month_abbr)])

    # Try to process each sub folder
    res_name_list = os.listdir(wrevap_ws)

    # Check keys on factors
    for res in evap_month_factor_dict.keys():
        if res not in res_name_list:
            logging.error('The reservoir name {} is invalid'.format(res_name))
            raise SystemExit()

    # Process each reservoir
    for res_name in res_name_list:
        # Reservoir folder must already exist
        if not os.path.isdir(res_name):
            continue
        # Process target reservoirs
        if res_list and res_name not in res_list:
            continue
        # Only process reservoirs needing adjustment
        if (res_name not in evap_month_factor_dict.keys() and
                res_name not in et_lake_month_adjust_dict.keys()):
            continue
        logging.info(res_name)
        res_ws = os.path.join(wrevap_ws, res_name)

        # Check if there are output data files
        output_folder_ws = os.path.join(res_ws, output_folder_name)
        if not os.path.isdir(output_folder_ws):
            logging.info('  The output data folder {} is not present'.format(
                output_folder_name))
            continue

        # Iteratively get data files
        for data_name in os.listdir(output_folder_ws):
            csv_path = os.path.join(output_folder_ws, data_name)
            if not os.path.isfile(csv_path):
                continue
            elif not csv_path.endswith('.csv'):
                continue
            logging.info('  {}'.format(
                data_name.split('.')[0].split('_', 1)[1]))

            # Read output CSV file
            csv_f = open(csv_path, 'rb')
            csv_lines = csv_f.readlines()
            csv_lines = [item.strip().split(',') for item in csv_lines]
            csv_f.close()

            # Scale evaporation values by monthly factor
            # This seems like a really messy way to do this
            month_col_i = csv_lines[0].index('MONTH')
            et_pot_col_i = csv_lines[0].index('ET_POT_MM')
            et_lake_col_i = csv_lines[0].index('ET_LAKE_MM')
            for csv_i, csv_line in enumerate(csv_lines[1:]):
                month = month_num_dict[csv_line[month_col_i]]
                et_pot = float(csv_line[et_pot_col_i])
                et_lake = float(csv_line[et_lake_col_i])
                # Reduce for advection then scale for frozen
                # Currently there are no reservoirs being adjusted for both
                if res_name in et_lake_month_adjust_dict.keys():
                    logging.debug('    Adding Monthly Evap: {} {}'.format(
                        month, et_lake_month_adjust_dict[res_name][month - 1]))
                    et_lake += et_lake_month_adjust_dict[res_name][month - 1]
                if res_name in evap_month_factor_dict.keys():
                    logging.debug('    Scaling Monthly Evap: {} {}'.format(
                        month, evap_month_factor_dict[res_name][month - 1]))
                    et_pot *= (1 - evap_month_factor_dict[res_name][month - 1])
                    et_lake *= (1 - evap_month_factor_dict[res_name][month - 1])
                csv_line[et_pot_col_i] = str(et_pot)
                csv_line[et_lake_col_i] = str(et_lake)
                csv_lines[csv_i + 1] = csv_line
                del month, et_pot, et_lake

            # Write data back to file
            csv_f = open(csv_path, 'wb')
            for line in csv_lines:
                csv_f.write(','.join(line) + '\n')
            csv_f.close()

    # os.chdir(wrevap_ws)


if __name__ == '__main__':
    workspace = os.getcwd()

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.info('\n{}'.format('#' * 80))
    log_f = '{:<20s} {}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    calc_gcm_wrevap(workspace)
