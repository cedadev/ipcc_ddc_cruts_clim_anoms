#!/usr/bin/env python
"""
Program     : Calculate the "sensitivity" to the differences in how the AR4 and AR5 data is calculated
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 05/04/2017
Requires    : CDO, Python
"""

import os, sys
import argparse
import subprocess

from CRUTS_var_list import CRUTS_var_list
from CRUTS_fix_metadata import fix_metadata
from calc_CRUTS_clim_anoms import get_input_directory, get_output_directory, cdo


def calculate_ref_period_sensitivity(var):
    # get input and output file path
    input_filepath = get_input_directory() + "/" + var + "/cru_ts3.24.01.1901.2015." + var + ".dat.nc"
    output_path = get_output_directory() + "/sens/" + var + "/netcdf"
    output_filepath = output_path + "/cru_ts3.24.01." + var + ".sens5-4.nc"
    # check if path exists
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    # build the cdo command - subtract the 1986->2005 mean from the 1961->1990 mean
    cmd = ["sub", "-ymonmean", "-selyear,1986/2005", input_filepath, "-ymonmean", "-selyear,1961/1990", input_filepath, output_filepath]
    cdo(cmd)
    # fix the metadata as we go
    fix_metadata(output_filepath, var, (1,20), "difference of monthly multi-year averages")


if __name__ == "__main__":
    # create the argument parser and the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", nargs=1, help="Variable:" + ("|").join(CRUTS_var_list.values()), default=None)
    args = parser.parse_args()

    # ensure that the arguments are valid
    assert(args.v == None or args.v[0] in CRUTS_var_list.values())

    # get the variable or use all
    if args.v == None:
        var_list = CRUTS_var_list.values()
    else:
        var_list = [args.v[0]]

    # loop over all variables and years in the climatological period(s)
    for v in var_list:
        calculate_ref_period_sensitivity(v)
