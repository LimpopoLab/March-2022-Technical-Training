#!/bin/bash

##Script for clipping to a smaller sub area Xai-Xai _clip_rp.tif files
import glob
from osgeo import gdal
import numpy as np

ls = (glob.glob("/Volumes/hydro3-raid/Technical_Training/hyp3/*/*_clip_rp.tif"))

for fn in ls:
    split = fn.split(".")
    out = split[0]+"_AOI.tif"
    shpin = "/Volumes/hydro3-raid/GIS/XaiXai_Area/Xai_Xai_Munic_Shape.shp"
    gdal.Warp(out,fn, cutlineDSName = shpin, cropToCutline = True)


