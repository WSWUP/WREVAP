#--------------------------------
# Name:         plot_gcm_wrevap.py
# Authors:      Charles Morton and Justin Huntington
# Modified:     2016-11-02
# Python:       2.7
#--------------------------------

import datetime as dt
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


def plot_gcm_wrevap(workspace):
    """Plot WREVAP data for each reservoir and GCM input data set"""
    logging.info('\nPlotting GCM WREVAP')

    # Output folder
    wrevap_ws = r'Z:\USBR_Ag_Demands_Project\reservoir_evap\WREVAP_GCM'

    year_range = (1950, 2099)

    res_list = []
    # res_list = ['Lake_Mead']
    res_list = ['Lake_Tahoe', 'Lake_Mead', 'Lake_Powell']
    # res_list = [
    #     'American_Falls', 'Boysen', 'Canyon_Ferry', 'Elephant_Butte',
    #     'Grand_Coulee', 'Lahontan', 'Millerton_Lake', 'Shasta_Lake',
    #     'Upper_Klamath']

    # WREVAP output CSV fields
    year_field = 'YEAR'
    month_field = 'MONTH'
    length_field = 'LENGTH'
    tmean_field = 'T_C'
    tdew_field = 'TD_C'
    rs_field = 'RS_MJ_M2_D'
    et_lake_field = 'ET_LAKE_MM'
    et_pot_field = 'ET_POT_MM'
    ppt_field = 'PPT_MM'

    # Subset ranges
    subset_dict = dict()
    subset_dict['American_Falls'] = (1992, 2011)
    subset_dict['Boysen'] = (1999, 2010)
    subset_dict['Canyon_Ferry'] = (1997, 2011)
    subset_dict['Elephant_Butte'] = (2002, 2009)
    subset_dict['Grand_Coulee'] = (2003, 2012)
    subset_dict['Lahontan'] = (2000, 2009)
    subset_dict['Lake_Mead'] = (2001, 2010)
    subset_dict['Lake_Powell'] = (1999, 2008)
    subset_dict['Lake_Tahoe'] = (2000, 2009)
    subset_dict['Millerton_Lake'] = (1989, 2010)
    subset_dict['Shasta_Lake'] = (2002, 2009)
    subset_dict['Upper_Klamath'] = (2000, 2011)

    # subfolder where WREVAP GCM output is stored
    output_folder_name = 'output_data'

    # Subfolder where WREVAP GCM plots are stored
    plots_folder_name = 'plots'

    # GIJ (calculated from Thornton Running is in MJ/m^2/day)
    # Correction factor of 11.574 will convert output W/m^2
    rs_correction_factor = 11.574

    # Overwrite existing files
    overwrite_flag = True

    #
    plot_si_units = False


    # Build plots folder
    # plots_ws = os.path.join(wrevap_ws, plots_folder_name)
    # if os.path.isdir(plots_ws) and overwrite_flag:
    #     for item in os.listdir(plots_ws):
    #         try: os.remove(os.path.join(plots_ws, item))
    #         except: pass
    # if not os.path.isdir(plots_ws): os.mkdir(plots_ws)

    year_index = np.arange(year_range[0], year_range[1] + 1)

    # Process each sub folder in wrevap_ws
    res_name_list = os.listdir(wrevap_ws)
    for res_name in res_name_list:
        # Reservoir folder must already exist
        if not os.path.isdir(res_name):
            continue
        # Process target reservoirs
        if res_list and res_name not in res_list:
            continue
        logging.info('{}'.format(res_name))
        res_ws = os.path.join(wrevap_ws, res_name)

        # Check output_data folder
        csv_ws = os.path.join(res_ws, output_folder_name)
        if not os.path.isdir(csv_ws):
            continue

        # Build plots folder
        plots_ws = os.path.join(res_ws, plots_folder_name)
        # if os.path.isdir(plots_ws) and overwrite_flag:
        #     for item in os.listdir(plots_ws):
        #         try: os.remove(os.path.join(plots_ws, item))
        #         except: pass
        if not os.path.isdir(plots_ws):
            os.mkdir(plots_ws)

        # Pre build list of output files
        csv_list = []
        for csv_name in os.listdir(csv_ws):
            csv_path = os.path.join(csv_ws, csv_name)
            if not os.path.isfile(csv_path):
                continue
            if not csv_path.endswith('.csv'):
                continue
            csv_list.append([csv_name, csv_path])

        # Prebuild arrays
        year_shape = (len(csv_list), len(year_index))
        yearly_length_array = np.zeros(year_shape)
        yearly_tmean_array = np.zeros(year_shape)
        yearly_tdew_array = np.zeros(year_shape)
        yearly_rs_array = np.zeros(year_shape)
        yearly_et_lake_array = np.zeros(year_shape)
        yearly_et_pot_array = np.zeros(year_shape)
        yearly_ppt_array = np.zeros(year_shape)
        yearly_et_net_array = np.zeros(year_shape)
        del year_shape

        # Process each file in output_data folder
        for csv_i, (csv_name, csv_path) in enumerate(csv_list):
            logging.info('  {}'.format(csv_name))

            # Read in data from each file
            # wrevap_array = np.loadtxt(
            #     csv_path, delimiter=',', skiprows=1)
            wrevap_array = np.genfromtxt(
                csv_path, delimiter=',', names=True, dtype=None)
            # wrevap_fields = wrevap_array.dtype.names
            year_array = wrevap_array[year_field]
            if not np.array_equal(np.unique(year_array), year_index):
                logging.warning('    Year ranges do not match, skipping')
            month_array = wrevap_array[month_field]
            monthly_length_array = wrevap_array[length_field]
            monthly_tmean_array = wrevap_array[tmean_field]
            monthly_tdew_array = wrevap_array[tdew_field]
            monthly_rs_array = wrevap_array[rs_field]
            monthly_et_lake_array = wrevap_array[et_lake_field]
            monthly_et_pot_array = wrevap_array[et_pot_field]
            monthly_ppt_array = wrevap_array[ppt_field]

            # Convert RS from MJ/m^2/day to W/m^2
            monthly_rs_array *= rs_correction_factor
            # np.times(
            #     monthly_rs_array, rs_correction_factor,
            #     out=monthly_rs_array)

            # Yearly means need to be weighted by the number of days in the month
            monthly_tmean_array *= monthly_length_array
            monthly_tdew_array *= monthly_length_array
            monthly_rs_array *= monthly_length_array

            # Calculate mean annual ET
            yearly_length_array = ndimage.sum(
                monthly_length_array, labels=year_array, index=year_index)
            yearly_tmean_array[csv_i, :] = ndimage.sum(
                monthly_tmean_array, labels=year_array, index=year_index)
            yearly_tdew_array[csv_i, :] = ndimage.sum(
                monthly_tdew_array, labels=year_array, index=year_index)
            yearly_rs_array[csv_i, :] = ndimage.sum(
                monthly_rs_array, labels=year_array, index=year_index)
            yearly_et_lake_array[csv_i, :] = ndimage.sum(
                monthly_et_lake_array, labels=year_array, index=year_index)
            yearly_et_pot_array[csv_i, :] = ndimage.sum(
                monthly_et_pot_array, labels=year_array, index=year_index)
            yearly_ppt_array[csv_i, :] = ndimage.sum(
                monthly_ppt_array, labels=year_array, index=year_index)

            # DEADBEEF Calculate for all GCMs
            # Net ET (ET - PPT)
            # yearly_et_net_array[csv_i, :] = (
            #     yearly_et_lake_array[csv_i, :] - yearly_ppt_array[csv_i, :])

            # Yearly means need to be weighted by the number of days in the month
            yearly_tmean_array[csv_i, :] /= yearly_length_array
            yearly_tdew_array[csv_i, :] /= yearly_length_array
            yearly_rs_array[csv_i, :] /= yearly_length_array

            # Cleanup
            del wrevap_array, month_array,
            del yearly_length_array, monthly_length_array
            del monthly_tmean_array
            del monthly_tdew_array
            del monthly_rs_array
            del monthly_et_lake_array
            del monthly_et_pot_array
            del monthly_ppt_array

        # Yearly means need to be weighted by the number of days in the month
        # yearly_tmean_array /= yearly_length_array
        # yearly_tdew_array /= yearly_length_array
        # yearly_rs_array /= yearly_length_array

        # Net ET (ET - PPT)
        yearly_et_net_array = yearly_et_lake_array - yearly_ppt_array

        # Calculate Medians
        median_tmean_array = np.median(yearly_tmean_array, axis=0)
        median_tdew_array = np.median(yearly_tdew_array, axis=0)
        median_rs_array = np.median(yearly_rs_array, axis=0)
        median_et_lake_array = np.median(yearly_et_lake_array, axis=0)
        median_et_pot_array = np.median(yearly_et_pot_array, axis=0)
        median_ppt_array = np.median(yearly_ppt_array, axis=0)
        median_et_net_array = np.median(yearly_et_net_array, axis=0)

        # Save yearly medians
        median_header_str = 'YEAR,TMEAN,TDEW,RS,ET_LAKE,ET_POT,PPT,ET_NET'
        median_fmt_str = '%d,%f,%f,%f,%f,%f,%f,%f'
        median_path = os.path.join(
            res_ws, res_name.lower() + '_medians.csv')
        if os.path.isfile(median_path) and overwrite_flag:
            os.remove(median_path)
        if not os.path.isfile(median_path):
            np.savetxt(
                median_path,
                np.vstack(
                    (year_index, median_tmean_array, median_tdew_array,
                     median_rs_array, median_et_lake_array,
                     median_et_pot_array, median_ppt_array,
                     median_et_net_array)).T,
                delimiter=',', header=median_header_str,
                comments='', fmt=median_fmt_str)

        # Save yearly medians subsets
        if subset_dict:
            # subset_range = subset_dict[res_name]
            subset_path = os.path.join(
                res_ws, res_name.lower() + '_medians_sub.csv')
            if os.path.isfile(subset_path) and overwrite_flag:
                os.remove(subset_path)
            if not os.path.isfile(subset_path):
                subset_mask = (
                    (year_index >= subset_dict[res_name][0]) &
                    (year_index <= subset_dict[res_name][1]))
                np.savetxt(
                    subset_path,
                    np.vstack(
                        (year_index[subset_mask],
                         median_tmean_array[subset_mask],
                         median_tdew_array[subset_mask],
                         median_rs_array[subset_mask],
                         median_et_lake_array[subset_mask],
                         median_et_pot_array[subset_mask],
                         median_ppt_array[subset_mask],
                         median_et_net_array[subset_mask])).T,
                    delimiter=',', header=median_header_str,
                    comments='', fmt=median_fmt_str)
                del subset_mask

        # Input data is in SI units
        if not plot_si_units:
            yearly_tmean_array = 1.8 * yearly_tmean_array + 32
            median_tmean_array = 1.8 * median_tmean_array + 32
            yearly_tdew_array = 1.8 * yearly_tdew_array + 32
            median_tdew_array = 1.8 * median_tdew_array + 32
            # yearly_rs_array *= 1
            # median_rs_array *= 1
            yearly_et_lake_array /= 25.4
            median_et_lake_array /= 25.4
            yearly_et_pot_array /= 25.4
            median_et_pot_array /= 25.4
            yearly_ppt_array /= 25.4
            median_ppt_array /= 25.4
            yearly_et_net_array /= 25.4
            median_et_net_array /= 25.4

        # Plot Mean Temperature
        plot_path = os.path.join(
            plots_ws, res_name.lower() + '_tmean.png')
        if os.path.isfile(plot_path) and overwrite_flag:
            os.remove(plot_path)
        if not os.path.isfile(plot_path):
            plt.figure()
            for tmean_array in yearly_tmean_array:
                plt.plot(
                    year_index, tmean_array,
                    '-', markersize=0, alpha=0.1, color='0.5')
            plt.plot(
                year_index, median_tmean_array,
                '-k', markersize=0, linewidth=2)
            plt.title('{}'.format(res_name.replace('_', ' ')))
            plt.xlabel('Year')
            if plot_si_units:
                plt.ylabel('Mean Annual Temperature [C]')
            else:
                plt.ylabel('Mean Annual Temperature [F]')
            plt.xlim([year_index[0], year_index[-1]])
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()

        # Plot Dew Point
        plot_path = os.path.join(
            plots_ws, res_name.lower() + '_tdew.png')
        if os.path.isfile(plot_path) and overwrite_flag:
            os.remove(plot_path)
        if not os.path.isfile(plot_path):
            plt.figure()
            for tdew_array in yearly_tdew_array:
                plt.plot(
                    year_index, tdew_array,
                    '-', markersize=0, alpha=0.1, color='0.5')
            plt.plot(
                year_index, median_tdew_array,
                '-k', markersize=0, linewidth=2)
            plt.title('{}'.format(res_name.replace('_', ' ')))
            plt.xlabel('Year')
            if plot_si_units:
                plt.ylabel('Mean Annual Dewpoint [C]')
            else:
                plt.ylabel('Mean Annual Dewpoint [F]')
            plt.xlim([year_index[0], year_index[-1]])
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()

        # Plot Thornton-Running Solar Radiation
        plot_path = os.path.join(
            plots_ws, res_name.lower() + '_rs.png')
        if os.path.isfile(plot_path) and overwrite_flag:
            os.remove(plot_path)
        if not os.path.isfile(plot_path):
            plt.figure()
            for rs_array in yearly_rs_array:
                plt.plot(
                    year_index, rs_array,
                    '-', markersize=0, alpha=0.1, color='0.5')
            plt.plot(
                year_index, median_rs_array,
                '-k', markersize=0, linewidth=2)
            plt.title('{}'.format(res_name.replace('_', ' ')))
            plt.xlabel('Year')
            if plot_si_units:
                plt.ylabel('Mean Annual Thornton-Running Solar [W/m^2]')
            else:
                plt.ylabel('Mean Annual Thornton-Running Solar [W/m^2]')
            plt.xlim([year_index[0], year_index[-1]])
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()

        # Plot Lake ET
        plot_path = os.path.join(
            plots_ws, res_name.lower() + '_et_lake.png')
        if os.path.isfile(plot_path) and overwrite_flag:
            os.remove(plot_path)
        if not os.path.isfile(plot_path):
            plt.figure()
            for et_lake_array in yearly_et_lake_array:
                plt.plot(
                    year_index, et_lake_array,
                    '-', markersize=0, alpha=0.1, color='0.5')
            plt.plot(
                year_index, median_et_lake_array,
                '-k', markersize=0, linewidth=2)
            plt.title('{}'.format(res_name.replace('_', ' ')))
            plt.xlabel('Year')
            if plot_si_units:
                plt.ylabel('Total Annual Lake ET [mm]')
            else:
                plt.ylabel('Total Annual Lake ET [in]')
            plt.xlim([year_index[0], year_index[-1]])
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()

        # Plot Potential ET
        plot_path = os.path.join(
            plots_ws, res_name.lower() + '_et_pot.png')
        if os.path.isfile(plot_path) and overwrite_flag:
            os.remove(plot_path)
        if not os.path.isfile(plot_path):
            plt.figure()
            for et_pot_array in yearly_et_pot_array:
                plt.plot(
                    year_index, et_pot_array,
                    '-', markersize=0, alpha=0.1, color='0.5')
            plt.plot(
                year_index, median_et_pot_array,
                '-k', markersize=0, linewidth=2)
            plt.title('{}'.format(res_name.replace('_', ' ')))
            plt.xlabel('Year')
            if plot_si_units:
                plt.ylabel('Total Annual Potential ET [mm]')
            else:
                plt.ylabel('Total Annual Potential ET [in]')
            plt.xlim([year_index[0], year_index[-1]])
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()

        # Plot Precipitation
        plot_path = os.path.join(
            plots_ws, res_name.lower() + '_ppt.png')
        if os.path.isfile(plot_path) and overwrite_flag:
            os.remove(plot_path)
        if not os.path.isfile(plot_path):
            plt.figure()
            for ppt_array in yearly_ppt_array:
                plt.plot(
                    year_index, ppt_array,
                    '-', markersize=0, alpha=0.1, color='0.5')
            plt.plot(
                year_index, median_ppt_array,
                '-k', markersize=0, linewidth=2)
            plt.title('{}'.format(res_name.replace('_', ' ')))
            plt.xlabel('Year')
            if plot_si_units:
                plt.ylabel('Total Annual Precipitation [mm]')
            else:
                plt.ylabel('Total Annual Precipitation [in]')
            plt.xlim([year_index[0], year_index[-1]])
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()

        # Plot Net ET (ET-PPT)
        plot_path = os.path.join(
            plots_ws, res_name.lower() + '_et_net.png')
        if os.path.isfile(plot_path) and overwrite_flag:
            os.remove(plot_path)
        if not os.path.isfile(plot_path):
            plt.figure()
            for et_net_array in yearly_et_net_array:
                plt.plot(
                    year_index, et_net_array,
                    '-', markersize=0, alpha=0.1, color='0.5')
            plt.plot(
                year_index, median_et_net_array,
                '-k', markersize=0, linewidth=2)
            plt.title('{}'.format(res_name.replace('_', ' ')))
            plt.xlabel('Year')
            if plot_si_units:
                plt.ylabel('Total Annual Net ET (ET - PPT) [mm]')
            else:
                plt.ylabel('Total Annual Net ET (ET - PPT) [in]')
            plt.xlim([year_index[0], year_index[-1]])
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()

        # Cleanup
        del plot_path
        del yearly_tmean_array
        del yearly_tdew_array
        del yearly_rs_array
        del yearly_et_lake_array
        del yearly_et_pot_array
        del yearly_ppt_array
        del yearly_et_net_array
        del res_ws, csv_ws, csv_list


if __name__ == '__main__':
    workspace = os.getcwd()

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.info('\n{}'.format('#' * 80))
    log_f = '{:<20s} {}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    plot_gcm_wrevap(workspace)
