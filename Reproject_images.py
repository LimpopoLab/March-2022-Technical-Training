#!/bin/bash

import glob
from osgeo import gdal

ls = (glob.glob("/Volumes/hydro3-raid/Technical_Training/hyp3/*/*_clip.tif"))

for fn in ls:
    split = fn.split(".")
    out = split[0]+"_rp.tif"
    gdal.Warp(out,fn,srcSRS='EPSG:32736',dstSRS='EPSG:4326')
