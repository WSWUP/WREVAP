#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import shutil

import numpy as np

def main(res_list_file, data_ws = None, passes = 0):
    """Wrapper for wrevap

    Args:
        res_list_file: file with list of reservoirs

    Returns:
        None
    """
    # print 'res_list_file\n", res_list_file

    bin_ws = r'C:\pythonApps\clre\WREVAP-denver\wrevap'

    # scenarios DMI manager

    wrevap_path = os.path.join(bin_ws, 'wrevap_gcm.py')
    # print "wrevap_path\n", wrevap_path
    
    # set data workspace to cwd if not specified
    
    if not data_ws:
        data_ws = os.getcwd()
    # print "Data workspace\n", data_ws
        
    # verify data workspace existence
    
    if not os.path.isdir(data_ws):
        print "Data workspace", data_ws, "does not exist."
        sys.exit()

    # Check reservoir list file path
    
    if not os.path.isfile(res_list_file):
        print 'Reservoir list file', res_list_file, 'does not exist.'
        sys.exit()
    elif not os.path.isdir(bin_ws):
        print 'Code workspace', bin_ws, 'does not exist.'
        sys.exit()
    elif not os.path.isfile(wrevap_path):
        print 'WRECAP path', wrevap_path, 'does not exist.'
        sys.exit()

    # read list of reservoirs
    
    res_ini_data_array = np.atleast_2d(np.genfromtxt(
        res_list_file, delimiter=',', names = True, dtype = None))
    columns = res_ini_data_array.dtype.names
    if not "NAME" in columns:
        print "'NAME' column is missing from reservoir list file."
        sys.exit()
    elif not "INI" in columns:
        print "'INI' column is missing from reservoir list file."
        sys.exit()
    elif not "DATA" in columns:
        print "'DATA' column is missing from reservoir list file."
        sys.exit()
    reservoirs = res_ini_data_array['NAME'].tolist()[0]
    ini_files = res_ini_data_array['INI'].tolist()[0]
    data_files = res_ini_data_array['DATA'].tolist()[0]
        
    # loop thru reservoirs
    
    for rcount, reservoir in enumerate(reservoirs):
        # specify and verify ini and input data file existence

        print "Processing reservoir", reservoir
        data_basename, data_ext = os.path.splitext(os.path.basename(data_files[rcount]))
        sol_path = os.path.join(data_ws, data_basename + ".SOL")
        tgw_path = os.path.join(data_ws, data_basename + ".TGW")
        base_ini_path = os.path.join(data_ws, ini_files[rcount])
        if not os.path.isfile(base_ini_path):
            print "reservoir", reservoir, "ini file\n", base_ini_path, "\ndoes not exist."
            sys.exit()
        temp_ini_path = os.path.join(data_ws, 'temp.ini')
        data_path = os.path.join(data_ws, data_files[rcount])
        if not os.path.isfile(data_path):
            print "reservoir", reservoir, "data file\n", data_path, "\ndoes not exist."
            sys.exit()

        # set up wrevap batch running
    
        args_list = ['python', wrevap_path]
        args_list.append('-i')
        args_list.append(temp_ini_path)
        args_list.append('--data')
        args_list.append(data_path)
        print "WREVAP command line is\n", args_list
        
        if passes == 0 or passes == 1:
            # first pass

            temp_f = open(temp_ini_path, 'w')
            with open(base_ini_path, 'r') as input_f:
                for line in input_f:
                    if line.strip().startswith('LK = '):
                        temp_f.write('LK = 2\n')
                    else:
                        temp_f.write(line)
            temp_f.close()
            print "Running first pass"
            proc = subprocess.Popen(args_list)
            proc.wait()
            if proc.returncode <> 0:
                print 'First pass of WREVAP run exited with error code', proc.returncode
                sys.exit()
            shutil.copy(sol_path, tgw_path)
        if passes == 0 or passes == 2:
            # second pass

            temp_f = open(temp_ini_path, 'w')
            with open(base_ini_path, 'r') as input_f:
                for line in input_f:
                    if line.strip().startswith('LK = '):
                        temp_f.write('LK = 3\n')
                    else:
                        temp_f.write(line)
            temp_f.close()
            print "Running second pass"
            proc = subprocess.Popen(args_list)
            proc.wait()
            if proc.returncode <> 0:
                print 'Second pass of WREVAP run exited with error code', proc.returncode
                sys.exit()

def parse_args():  
    parser = argparse.ArgumentParser(
        description = 'Scenarios DMI',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-r', '--rlf', metavar = 'PATH',
        type = lambda x: is_valid_file(parser, x), help = 'Reservoir list file')
    parser.add_argument(
        '-w', '--dws', metavar = 'data_ws', default = None,
        help = "Data workspace")
    parser.add_argument(
        '-p', '--passes', metavar = 'passes', default = 0,
        help = "Passes (0 (both), 1, or 2")
    args = parser.parse_args()

    # Convert reservoir list path to an absolute path if necessary
    
    if args.rlf and os.path.isfile(os.path.abspath(args.rlf)):
        args.rlf = os.path.abspath(args.rlf)
    return args

def get_res_list_file(workspace):
    import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    res_list_file = tkFileDialog.askopenfilename(
        initialdir = workspace, parent = root, filetypes = [('INI files', '.ini')],
        title = 'Select INI file')
    root.destroy()
    return res_list_file

def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        local_path = os.path.join(os.getcwd(), arg)
        if os.path.isfile(local_path):
           return local_path
        else:
           parser.error('Neither {} or {} exists!'.format(arg, local_path))
    else:
        return arg

def is_valid_directory(parser, arg):
    if not os.path.isdir(arg):
        parser.error('The directory {} does not exist!'.format(arg))
    else:
        return arg

if __name__ == '__main__':
    args = parse_args()
    # print "args", args
    if args.rlf:
        if args.dws:
            main(res_list_file = args.rlf, data_ws = args.dws, passes = int(args.passes))
        else:
            main(res_list_file = args.rlf, data_ws = os.getcwd(), passes = int(args.passes))
    else:
        print "Need to specify reservoir list file."
