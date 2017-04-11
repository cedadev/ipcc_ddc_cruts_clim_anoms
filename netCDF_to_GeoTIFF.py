#!/usr/bin/env python
"""
Program     : Plot a netCDF file as a GeoTIFF.  Very simplistic at the moment and will probably only work for global files
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 07/04/2017
Requires    : Python, GDAL
"""

import gdal
from gdalconst import *
from osgeo import osr
import argparse
import numpy
import sys
import os

def create_color_table(R,G,B):
    """Create a color table for use with a palettized GeoTIFF
       :param array(byte) R: array / list of reds, length n
       :param array(byte) G: array / list of greens, length n
       :param array(byte) B: array / list of blues, length n
    """
    assert(len(R) == len(G) and len(G) == len(B))

    # create a GDAL color table
    ct = gdal.ColorTable()
    for p in range(0, len(R)):
       ct.SetColorEntry(p, (R[p], G[p], B[p], 255)) 
    return ct
       

def load_color_table(nc_vname):
    """Load a color table for use with a palettized GeoTIFF
       :param string nc_vname: name of the variable which maps to a color table.  The color tables are stored in the directory color_tables with a filename extension of .ct
    """

    if nc_vname in ["tmp", "tmx", "tmn"]:
        ct_name = "tmp"
    elif nc_vname in ["cld"]:
        ct_name = "cld"
    elif nc_vname in ["wet", "vap", "pre", "pet"]:
        ct_name = "pre"
    elif nc_vname in ["frs"]:
        ct_name = "frs"
    elif nc_vname in ["dtr"]:
        ct_name = "dtr"
    else:
        raise Exception("No color table defined for variable: "+nc_vname)

    fname = "color_tables/"+ct_name+".ct"
    # file is just a text file, use readlines
    fh = open(fname, 'r')
    lines = fh.readlines()

    # create GDAL color table
    ct = gdal.ColorTable()

    # first line is length of color table
    ct_length = int(lines[0])
    for p in range(0, ct_length):
        c = lines[p+1].split(",")
        ct.SetColorEntry(p, (int(c[0]), int(c[1]), int(c[2]), int(c[3])))
    return ct


def get_scaled_data(src_ds, nc_vname, mini=None, maxi=None):
    """Get the data from the source dataset and scale it to a byte
      :param GDALDataset src_ds: source dataset
      :param string nc_vname: name of variable
      :mini: minimum value to use in scaling (or None to calculate from data)
      :maxi: maximum value to use in scaling (or None to calculate from data)
    """
    # get the data, as a masked array
    data = src_ds.ReadAsArray()
    mv = numpy.float64(src_ds.GetMetadataItem(nc_vname+"#_FillValue"))
    src_data_ma = numpy.ma.masked_equal(data, mv)

    # scale the data - reserve 255 as nodata
    # get the min and max values first - either from the parameters (if not None) or calculate from the array
    if mini == None:
        min = numpy.min(src_data_ma)
    else:
        min = float(mini)
    if maxi == None:
        max = numpy.max(src_data_ma)
    else:
        max = float(maxi)

    offset = min
    scale = (max - min) / 254
    dst_data_ma = (src_data_ma - offset) / scale
    # fill in the mask as 255
    dst_data = numpy.ma.filled(dst_data_ma, 255)
    dst_data_byte = dst_data.astype("u8")
    return dst_data_byte, offset, scale


def write_geotiff(filename, nc_vname, data, geotransform, projection, metadata):
    # get the size of the rasters and number of rasters to write out
    if data.ndim == 2:
        Xsize = data.shape[1]
        Ysize = data.shape[0]
        nR = 1
    elif data.ndim == 3:
        Xsize = data.shape[2]
        Ysize = data.shape[1]
        nR = data.shape[0]
    else:
        raise Exception("Data does not have correct number of dimensions (2 or 3)")

    # write out the GeoTIFF
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(filename, Xsize, Ysize, nR, gdal.GDT_Byte,
                           ["COMPRESS=LZW", "INTERLEAVE=BAND"])
    
    # copy most of the information from the netCDF file
    dst_ds.SetGeoTransform(geotransform)

    # set a projection - WGS 84 if a projection isn't given
    if projection == "":
        proj = osr.GetWellKnownGeogCSAsWKT( "WGS84" )
        dst_ds.SetProjection(proj)
    else:
        dst_ds.SetProjection(projection)

    # set the color table
    rb = dst_ds.GetRasterBand(1)
    rb.WriteArray(data)

    # set the nodata value
    rb.SetNoDataValue(255)

    # only 2D data can have a color table
    if data.ndim == 2:
        # create the color table for the variable
        ct = load_color_table(nc_vname)
        rb.SetColorTable(ct)

    # set the metadata
    dst_ds.SetMetadata(metadata)
    dst_ds = None
    

def get_metadata(metadata, nc_vname):
    """Manipulate or select the input metadata to form a subset."""
    output_metadata = {}
    output_metadata["UNITS"] = metadata[nc_vname+"#units"]
    return output_metadata


def netcdf_to_geotiff(nc_fname, nc_vname, geotiff_name, one_file_per_tstep=True, min=None, max=None, extra_metadata={}, tstep_list=[]):
    """Convert a netcdf file to a GeoTIFF file using GDAL.
       :param string nc_fname: filename of the netCDF file
       :param string nc_vname: variable to extract from the netCDF file
       :param boolean one_file_per_tstep: output each timestep in the netCDF file as a separate GeoTIFF file (default=True)
       :param string geotiff_name: filename to save GeoTIFF to (optional)
       :param float min: minimum value to use in the scaling - if None then calculate from netCDF file
       :param float max: maximum value to use in the scaling - if None then calculate from netCDF file
    """
    # open the netCDF file
    src_ds = gdal.Open('NETCDF:"'+nc_fname+'":'+nc_vname, GA_ReadOnly)
    
    # get the scaled data - need the offset and scale to write the metadata
    dst_data, offset, scale = get_scaled_data(src_ds, nc_vname, min, max)

    # create the metadata
    metadata = get_metadata(src_ds.GetMetadata(), nc_vname)
    metadata["SCALE"] = str(scale)
    metadata["OFFSET"] = str(offset)

    # add any extra metadata
    for e in extra_metadata.keys():
        metadata[e] = extra_metadata[e]  # this will overwrite any existing metadata, which is desirable

    # write out the GeoTIFF - have to choose here between packing multiple bands into a strictly grayscale image
    # or writing one image per timestep with an associated colormap
    if one_file_per_tstep and dst_data.ndim > 2:
        # calculate the number of zeros needed to do the padding
        n_zeros = len(str(dst_data.shape[0]))
        # create a directory to put the individual files in
        path = geotiff_name[:-4]
        if not os.path.exists(path):
            os.makedirs(path)
        # create the path to the file
        fname_stub = path.split("/")[-1]
        # loop over and create each filename and then each TIFF
        for i in range(0, dst_data.shape[0]):
            # create the filename with _timestep suffix
            fmt_str = "{:0"+str(n_zeros)+"d}"
            t_step_number = fmt_str.format(i+1)
            fname = path + "/" + fname_stub + "_" + t_step_number + ".tif"
            # get the month / day / etc. out of the timestep list
            if tstep_list != [] and i < len(tstep_list) and "TIFFTAG_IMAGEDESCRIPTION" in metadata.keys():
                tstep_ref = tstep_list[i]
                metadata["TIFFTAG_IMAGEDESCRIPTION"] += tstep_ref
            # write the geotiff
            write_geotiff(fname, nc_vname, dst_data[i], src_ds.GetGeoTransform(), src_ds.GetProjection(), metadata)
                
    else:
        write_geotiff(geotiff_name, nc_vname, dst_data, src_ds.GetGeoTransform(), src_ds.GetProjection(), metadata)

if __name__ == "__main__":
    # create the argument parser and the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", nargs=1, help="netCDF input file name")
    parser.add_argument("-v", nargs=1, help="netCDF variable name")
    parser.add_argument("-o", nargs=1, help="GeoTIFF output file name (optional)", default=None)
    parser.add_argument("-m", nargs=1, help="Minimum value to use in scaling netCDF values to bytes (optional - if not given then calculate from the netCDF file)", default=None)
    parser.add_argument("-x", nargs=1, help="Maximum value to use in scaling netCDF values to bytes (optional - if not given then calculate from the netCDF file)", default=None)
    parser.add_argument("-t", action="store_true", help="Write each timestep in the netCDF file as a separate GeoTIFF file (optional)", default=True)
    args = parser.parse_args()

    # ensure that the arguments are valid
    error = False
    if not args.f:
        print "ERROR: input filename not supplied"
        error = True

    if not args.v:
        print "ERROR: input variable name not supplied"
        error = True
        
    if error:
        parser.print_help()
        sys.exit()

    # get the output filename
    if args.o == None:
        output_fname = args.f[0][:-3]+".tif"
    else:
        output_fname = args.o[0]

    # run the conversion
    if args.m:
        min = args.m[0]
    else:
        min = None

    if args.x:
        max = args.x[0]
    else:
        max = None

    netcdf_to_geotiff(args.f[0], args.v[0], output_fname, args.t, min, max, extra_metadata={"CREATOR":"Neil Massey, CEDA"})
