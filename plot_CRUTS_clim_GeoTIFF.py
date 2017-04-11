
#!/usr/bin/env python
"""
Program     : Plot a the CRUTS derived climatologies as GeoTIFFs.
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 07/04/2017
Requires    : Python, GDAL
"""

from netCDF_to_GeoTIFF import netcdf_to_geotiff
import argparse
from netCDF4 import Dataset
import numpy
import os, sys
from copy import copy

def GetMinMaxValues(fname, vname):
    nc_fh = Dataset(fname, 'r')
    var = nc_fh.variables[vname]
    min = numpy.min(var[:])
    max = numpy.max(var[:])

    nc_fh.close()
    return min, max


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", nargs=1, help="netCDF variable name")
    args = parser.parse_args()

    error = False
    if not args.v:
        print "ERROR: input variable name not supplied"
        error = True

    if error:
        parser.print_help()
        sys.exit()

    v = args.v[0]
    extra_metadata={"CREATOR":"Centre for Environmental Data Analysis (CEDA), STFC",
                    "CONTACT":"support@ceda.ac.uk",
                    "REFERENCE":"IPCC DDC climatologies, CRU TS3.24.01"}
    if v == "tmp":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Near surface air temperature"
        extra_metadata["TITLE"] = "Near surface air temperature [celsius];"
    elif v == "tmx":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Near surface air temperature maximum"
        extra_metadata["TITLE"] = "Near surface air temperature maximum [celsius];"
    elif v == "tmn":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Near surface air temperature minimum"
        extra_metadata["TITLE"] = "Near surface air temperature minimum [celsius];"
    elif v == "cld":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Cloud cover"
        extra_metadata["TITLE"] = "Cloud cover [percentage];"
    elif v == "pet":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Potential Evapotranspiration"
        extra_metadata["TITLE"] = "Potential evapotranspiration [mm/day];"
    elif v == "dtr":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Near surface air temperature: diurnal range"
        extra_metadata["TITLE"] = "Near surface air temperature: diurnal range [celsius];"
    elif v == "frs":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Ground frost frequency"
        extra_metadata["TITLE"] = "Ground frost frequency [days];"
    elif v == "pre":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Precipitation"
        extra_metadata["TITLE"] = "Precipitation [kg m-2];"
    elif v == "vap":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Water vapor pressure"
        extra_metadata["TITLE"] = "Water vapor pressure [hPa];"
    elif v == "wet":
        extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] = "Wet day frequency"
        extra_metadata["TITLE"] = "Wet day frequency [days];"

    extra_metadata["TIFFTAG_IMAGEDESCRIPTION"] += " climatology"
    extra_metadata["TITLE"] += " CRU TS3.24.01"

    in_dir = "/group_workspaces/jasmin/cedaproc/nrmassey/CRU_TS3.24.01/derived/clim/"+v+"/netcdf/"
    out_dir = "/group_workspaces/jasmin/cedaproc/nrmassey/CRU_TS3.24.01/derived/clim/"+v+"/geotiff/"
    files = os.listdir(in_dir)
    # get the min and max values so that all the plots have the same scale
    min = 2e20
    max = -2e20
    for f in files:
        c_min, c_max = GetMinMaxValues(in_dir+f, v)
        if c_min < min:
            min = c_min
        if c_max > max:
            max = c_max

    # timestep reference list
    tstep_list = ["January", "February", "March", "April", "May", "June",
                  "July", "August", "September", "October", "November", "December"]

    for i in range(0, len(tstep_list)):
        tstep_list[i] = ", monthly mean for " + tstep_list[i]

    for f in files:
        output_fname = f[:-3] + ".tif"
        # get the year and the month and add to the metadata - we need a copy of the metadata
        copy_metadata = copy(extra_metadata)
        # get the start and end year
        fsplit = f.split(".")
        sy = fsplit[3]
        ey = fsplit[4]
        copy_metadata["TIFFTAG_IMAGEDESCRIPTION"] += ", "+sy + "-" + ey
        netcdf_to_geotiff(in_dir+f, v, out_dir+output_fname, True, min, max, extra_metadata=copy_metadata, tstep_list=tstep_list)
