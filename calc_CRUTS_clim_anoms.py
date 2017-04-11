#!/usr/bin/env python
"""
Program     : Calculate mean climatology anomaly fields of CRUTS data for upload to the IPCC Data Distribution Centre
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 05/04/2017
Requires    : CDO, Python
"""

import os, sys
import argparse
import subprocess

from CRUTS_climatological_periods import CRUTS_climatological_periods
from CRUTS_var_list import CRUTS_var_list


def get_input_directory():
    """You have to copy the CRUTS files to this directory - they are gzipped in the archive"""
    return "/group_workspaces/jasmin/cedaproc/nrmassey/CRU_TS3.24.01/data"


def get_output_directory():
    return "/group_workspaces/jasmin/cedaproc/nrmassey/CRU_TS3.24.01/derived"

DEBUG = True

def cdo(cmd_params):
    """Run the cdo command using subprocess.check_output"""
    cmd = ["/usr/bin/cdo", "-s", "--no_warnings"]
    cmd.extend(cmd_params)
    if DEBUG:
        print " ".join(cmd)
    return subprocess.check_output(cmd)


def calculate_climatological_mean_anomaly(var, years):
    # get input and output file path
    input_filepath = get_input_directory() + "/" + var + "/cru_ts3.24.01.1901.2015." + var + ".dat.nc"
    output_path = get_output_directory() + "/clim-anom-ref5/netcdf/" + var
    output_filepath = output_path + "/cru_ts3.24.01." + str(years[0]) + "." + str(years[1]) + "." + var + ".clim-anom-ref5.nc"
    # check if path exists
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    # build the cdo command - subtract the 1986->2005 mean from the mean over the period years[0]->years[1]
    cmd = ["sub", "-ymonmean", "-selyear,"+str(years[0])+"/"+str(years[1]), input_filepath, "-ymonmean", "-selyear,1986/2005", input_filepath, output_filepath]
    cdo(cmd)


if __name__ == "__main__":
    # create the argument parser and the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", nargs=1, help="Variable:" + ("|").join(CRUTS_var_list.values()), default=None)
    parser.add_argument("-c", nargs=1, help="Climatological period: " + ("|").join(CRUTS_climatological_periods), default=None)
    args = parser.parse_args()

    # ensure that the arguments are valid
    assert(args.v == None or args.v[0] in CRUTS_var_list.values())
    assert(args.c == None or (args.c[0] in CRUTS_climatological_periods))

    # get the variable or use all
    if args.v == None:
        var_list = CRUTS_var_list.values()
    else:
        var_list = [args.v[0]]

    if args.c == None:
        clim_periods = CRUTS_climatological_periods.keys()
    else:
        clim_periods = [args.c[0]]

    # loop over all variables and years in the climatological period(s)
    for v in var_list:
        for cp in clim_periods:
            for years in CRUTS_climatological_periods[cp]:
                calculate_climatological_mean_anomaly(v, years)
