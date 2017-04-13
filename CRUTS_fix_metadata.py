#!/usr/bin/env python
"""
Program     : Fix the metadata in the files created by calc_CMIP5_clim_anoms
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 14/03/2017
Modified    : 23/03/2017
Requires    : CDO, NCO, Python, drslib
"""
import os, sys

from drslib.drs import CmipDRS
from drslib import cmip5

from CRUTS_climatological_periods import CRUTS_climatological_periods
from CRUTS_var_list import CRUTS_var_list
from calc_CRUTS_clim_anoms import get_output_directory

from netCDF4 import Dataset


def fix_metadata(fname, vname, years, proc_desc):
    """Fix the metadata for the CMIP5 anomaly file
       :param string fname: name of file to fix metadata for
    """
    # open the dataset
    nc_fh = Dataset(fname, 'a')
    # get the attributes
    nc_attrs = nc_fh.ncattrs()
    # fix the history attribute first
    hist_attr = nc_fh.getncattr("history").split("\n")
    for h in range(0, len(hist_attr)):
        # make the history string shorter
        if "cdo -s --no_warnings" in hist_attr[h]:
            s = hist_attr[h].split(" ")
            new_s = []
            for t in s:
                if "datacentre" in t or "group_workspaces" in t:
                    t = t.split("/")[-1]
                new_s.append(t)
                
            s = new_s[:]
            # fix an introduced mistake
            if "2005" in new_s:
                idx_2005 = new_s.index("2005")
                new_s[idx_2005] = "-selyear,1986/2005"
            if str(1900+years[1]) in new_s:
                idx_y1 = new_s.index(str(1900+years[1]))
                new_s[idx_y1] = "-selyear,"+str(1900+years[0])+"/"+str(1900+years[1])
            if str(2000+years[1]) in new_s:
                idx_y1 = new_s.index(str(2000+years[1]))
                new_s[idx_y1] = "-selyear,"+str(2000+years[0])+"/"+str(2000+years[1])

            cdo_pos = s.index("cdo")
            cdo_string = " ".join(new_s[cdo_pos:])
            hist_attr[h] = cdo_string
    new_hist_attr = "\n".join(hist_attr)

    # fix the cdo history
    nc_fh.setncattr("history", new_hist_attr)
    # add some other descriptive data about the content / processing
    nc_fh.setncattr("processing_description", proc_desc)
    nc_fh.setncattr("processing_type", "mean average")
    nc_fh.setncattr("average_dimension", "time [nested]")
    n_yrs = years[1] - years[0] + 1
    if (n_yrs == 10):
        avg_type = "decadal"
    elif (n_yrs == 20):
        avg_type = "20 years"
    elif (n_yrs == 30):
        avg_type = "30 years"
    else:
        avg_type = "yearly"
    nc_fh.setncattr("average_type_outer", avg_type)
    nc_fh.setncattr("average_type_inner", "monthly")
    nc_fh.setncattr("anomaly_reference_period", "1986-2005")
    nc_fh.setncattr("contact", "support@ceda.ac.uk")
    nc_fh.setncattr("reference", "IPCC DDC AR5 climatologies")
    nc_fh.close()


def fix_metadata_email(fname, vname, years):
    nc_fh = Dataset(fname, 'a')
    nc_attrs = nc_fh.ncattrs()
    nc_fh.setncattr("contact", "support@ceda.ac.uk")    
    nc_fh.close()


if __name__ == "__main__":

    # loop over all variables and years in the climatological period(s)
    for var in CRUTS_var_list.values():
        for cp in CRUTS_climatological_periods:
            for years in CRUTS_climatological_periods[cp]:
                clim_anom_path = get_output_directory() + "/clim-anom-ref5/" + var + "/netcdf"
                clim_anom_filepath = clim_anom_path + "/cru_ts3.24.01." + str(years[0]) + "." + str(years[1]) + "." + var + ".clim-anom-ref5.nc"
                print clim_anom_filepath
                fix_metadata(clim_anom_filepath, var, years,"anomalies of monthly multi-year averages")

                clim_path = get_output_directory() + "/clim/" + var + "/netcdf"
                clim_filepath = clim_path + "/cru_ts3.24.01." + str(years[0]) + "." + str(years[1]) + "." + var + ".clim.nc"
                print clim_filepath
                fix_metadata(clim_filepath, var, years, "monthly multi-year averages")
