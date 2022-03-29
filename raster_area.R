# Limpopo Resilience Lab
# Technical Training - March 2022
# Remote sensing image analysis (water quantification) with GDAL

# This example code will measure the area of water based on NDWI values within 
# a pre-programmed areas.

# This work was supported by the United States Agency for International 
# Development, Southern Africa Regional Mission, Fixed Amount Award 
# 72067419FA00001. This work reflects the work of the authors and does not 
# necessarily reflect the views of USAID or the United States Government.  
# Further information is available at: 
# www.duq.edu/limpopo 
# https://github.com/LimpopoLab 

library(XML) # For reading xml metadata file
library(raster) # to import image
library(sp)
library(rgeos) # gCovers - test extents, not needed to crop
library(dplyr) # used first in rename - general data frame management 
library(ggplot2) # plotting tools (histogram)

# Set your working directory.  Avoid spaces in your path and OneDrive
setwd("/Volumes/T7/planet/weir/")

# ------------------------------------------------------------------------------
## IMPORT IMAGE, CROP, ADJUST TO TOA(if needed)
# Import raster image
pic <- stack("20210802_072057_91_2430_3B_AnalyticMS_SR_8b_harmonized.tif") # Aug
# pic <- stack("20210903_071821_66_2421_3B_AnalyticMS_SR_8b_harmonized.tif") # Sep
# pic <- stack("20211009_071655_67_2459_3B_AnalyticMS_SR_8b_harmonized.tif") # Oct
# pic <- stack("20211202_071519_29_2420_3B_AnalyticMS_SR_8b_harmonized.tif") # Dec
# pic <- stack("20220102_071539_49_242d_3B_AnalyticMS_SR_8b_harmonized.tif") # Jan

# To crop the image, first, check the extent, then crop
## Mutale Weir
crp_ext <- as(extent(246825, 247075, 7479485, 7479710), 'SpatialPolygons') # Desired extent and crop limits
crs(crp_ext) <- "+proj=utm +zone=36 +datum=WGS84"
pic_ext <- as(extent(pic), 'SpatialPolygons') # Extent of image
crs(pic_ext) <- "+proj=utm +zone=36 +datum=WGS84"

if (gCovers(pic_ext,crp_ext)) { # returns TRUE if no point in spgeom2 (e, needed) is outside spgeom1 (test, image extent) # used to be (gWithin(e, test, byid = FALSE))
     pic <- crop(pic, crp_ext)
}

# If you need the reflectance coefficients from Planet metadata xml file
# md <- xmlParse("20210219_074647_1009_3B_AnalyticMS_metadata.xml")
# rc <- setNames(xmlToDataFrame(node=getNodeSet(md, "//ps:EarthObservation/gml:resultOf/ps:EarthObservationResult/ps:bandSpecificMetadata/ps:reflectanceCoefficient")),"reflectanceCoefficient")
# rc <- as.matrix(rc)
# # 1 Red
# # 2 Green
# # 3 Blue
# # 4 Near infrared
# rc2 <- as.numeric(rc[2]) # Green
# rc4 <- as.numeric(rc[4]) # NIR
# rm(rc, md)
# surface reflectance (SR) images do not need adjustment

# Calculate NDWI using the green and NIR bands: (G-NIR)/(G+NIR)
# 4-band: green (band 2) and NIR (band 4)
#ndwi <- ((pic[[2]]) - (pic[[4]])) / ((pic[[2]]) + (pic[[4]]))
# 8-band: green (band 4) and NIR (band 8)
ndwi <- ((pic[[4]]) - (pic[[8]])) / ((pic[[4]]) + (pic[[8]]))

# To see new NDWI image
rm(pic) # to free up RAM
ndwi_df <- as.data.frame(ndwi, xy = TRUE)
ggplot() + 
     geom_raster(data=ndwi_df, aes(x=x,y=y,fill=layer)) +
     scale_fill_gradientn(colours = topo.colors(5), limits = c(-1,0.1)) +
     theme(axis.text = element_blank(), axis.title = element_blank())

writeRaster(ndwi, 
            filename= "20220102_ndwi.tif",
            format = "GTiff", # save as a tif, save as a FLOAT if not default, not integer
            overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.

# ndwi_dry <- stack("20210802_ndwi.tif")
# ndwi_wet <- stack("20220102_ndwi.tif")

# ------------------------------------------------------------------------------
## TO DETERMINE WATER THRESHOLD
# Build histogram 
h = hist(ndwi, # built-in histogram function
         breaks=seq(-1,1,by=0.01), # span of possible NDWI values
         plot=FALSE) 

bins <- h$mids # positions    number
v <- h$counts # counts        integer

# Allocate arrays used in analysis
windows <- 10 # number of averaging windows
avg <- array(0, dim = c(length(bins),windows)) # moving average filter
peaks <- array(0, dim = c(length(bins),windows)) # peak locations
nop <- array(0, dim = c(windows)) # number of peaks found

# Loop around number averaging window values
for (w in 1:10){
     # filter values (v=h$counts) with the averaging window size 2*w+1; therefore, w is the radius of the window
     for (k in (w+1):(200-w)){ # assume tails are unimportant, default to zero based on preallocation
          avg[k,w] <- ((sum(v[(k-w):(k+w)]))/((2*w)+1))
     }
     # identify and number peaks
     cnt <- 0
     for (j in (w+1):(200-w)){
          if ( ((avg[j-1,w])<(avg[j,w])) & ((avg[j+1,w])<(avg[j,w])) ) {
               cnt <- (cnt+1) # if a peak, add one to count
               peaks[j,w] <- cnt # add peak count to location
               nop[w] <- cnt # count peaks
               
          }
     }
}

# Set error values for the result vectors in case neither two nor three peaks are found:
threepeak <- -1
twopeak <- -1

# Detect peaks
for (w in 1:10){
     # testing in three peaks
     # due to the order of the w variable, only the 'smoothest' result will be kept
     if ((nop[w])==3){
          # finds the second and third peak
          for (j in 1:200){
               if ((peaks[j,w])==2){
                    sec <- j # stores the index of the second peak
               }
               if ((peaks[j,w])==3){
                    thr <- j # stores the index of the third peak
               }
          }
          # finds minimum between second and third peak
          m <- max(v) # create variable for minimum, initially set higher than any value
          for (j in (sec):(thr)){
               if ((avg[j,w])<m){
                    goal <- j
                    m <- avg[j,w]
               }
          }
          threepeak <- (bins[(goal)])
     }
     # test in case exactly three peaks were not found
     if ((nop[w])==2){
          # find the position of the first and second (the only) peaks
          for (j in 1:200){
               if ((peaks[j,w])==1){
                    fst <- j # stores the index of the second peak
               }
               if ((peaks[j,w])==2){
                    sec <- j # stores the index of the third peak
               }
          }
          # finds minimum between first and second peak
          m <- max(v) # create variable for minimum, initially set higher than any value
          for (j in (fst):(sec)){
               if ((avg[j,w])<m){
                    goal <- j
                    m <- avg[j,w]
               }
          }
          twopeak <- (bins[(goal)])
     }
}

# Store values
ndwi_values <- data.frame(ndwi@data@values)
ndwi_values <- rename(ndwi_values, data=ndwi.data.values)

# histogram visualization
h <- ggplot(ndwi_values, aes(x=data)) +
     geom_histogram(breaks = (c(0:200)/100-1), color = "black", fill = "gray", na.rm = TRUE) +
     geom_vline(aes(xintercept = twopeak), color = "green") +
     geom_vline(aes(xintercept = threepeak), color = "blue") +
     xlab("NDWI") +
     ylab("Count") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 24)) +
     theme(axis.title = element_text(face = "plain", size = 24))
ggsave("20210102_hist.eps", h, device = "eps", dpi = 72)

# Analysis
# The two-peak analysis provides a more stable threshold value by visual 
# inspection; therefore, we will set the threshold to twopeak:
ndwi_threshold <- twopeak

# ------------------------------------------------------------------------------
## CALCULATE AREA

# If it is sufficient to compare number of pixels:
no_of_pixels <- sum(ndwi[] >= ndwi_threshold)

# With raster:area
# Need to reproject:
ndwi_WGS84 <- projectRaster(ndwi, crs = "+proj=longlat +datum=WGS84 +no_defs")
sz <- area(ndwi_WGS84) # finds areas
sm <- 0
for (i in 1:length(sz@data@values)) {
     if (is.na(ndwi_WGS84@data@values[i]) == FALSE) {
          if (ndwi_WGS84@data@values[i] >= ndwi_threshold) {
               sm <- sm + sz@data@values[i]
          }
     }
}

