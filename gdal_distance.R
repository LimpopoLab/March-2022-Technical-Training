# Limpopo Resilience Lab
# Technical Training - March 2022
# Remote sensing image analysis (water quantification) with GDAL

# This example code will measure the width of a river at a pre-programmed 
# location based on NDWI values.

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


library(rgdal) # measurement tools.  Must change to GDAL and PROJ: sf/stars/terra, by 2023; although, these are not available for all R versions at this time.

library(methods)

library(MASS)

library(stringr)

library(lubridate)
library(readr)


# Set your working directory.  Avoid spaces in your path and OneDrive
setwd("/Volumes/T7/planet/mutale/2021/20210219_074647_1009/analytic_udm2/")

# ------------------------------------------------------------------------------
## IMPORT IMAGE, CROP, ADJUST TO TOA(if needed)
# Import raster image
pic <- stack("20210219_074647_1009_3B_AnalyticMS.tif")

# To crop the image, first, check the extent, then crop
## Mutale River downstream
crp_ext <- as(extent(245850, 246350, 7478700, 7479200), 'SpatialPolygons') # Desired extent and crop limits
crs(crp_ext) <- "+proj=utm +zone=36 +datum=WGS84"
pic_ext <- as(extent(pic), 'SpatialPolygons') # Extent of image
crs(pic_ext) <- "+proj=utm +zone=36 +datum=WGS84"

if (gCovers(pic_ext,crp_ext)) { # returns TRUE if no point in spgeom2 (e, needed) is outside spgeom1 (test, image extent) # used to be (gWithin(e, test, byid = FALSE))
     pic <- crop(pic, crp_ext)
}

# If you need the reflectance coefficients from Planet metadata xml file
md <- xmlParse("20210219_074647_1009_3B_AnalyticMS_metadata.xml")
rc <- setNames(xmlToDataFrame(node=getNodeSet(md, "//ps:EarthObservation/gml:resultOf/ps:EarthObservationResult/ps:bandSpecificMetadata/ps:reflectanceCoefficient")),"reflectanceCoefficient")
rc <- as.matrix(rc)
# 1 Red
# 2 Green
# 3 Blue
# 4 Near infrared
rc2 <- as.numeric(rc[2]) # Green
rc4 <- as.numeric(rc[4]) # NIR
rm(rc, md)

# Calculate NDWI using the green (band 2) and nir (band 4) bands
ndwi <- ((rc2*pic[[2]]) - (rc4*pic[[4]])) / ((rc2*pic[[2]]) + (rc4*pic[[4]]))

# To see new NDWI image
plot(ndwi)
rm(pic) # to free up RAM

writeRaster(x = ndwi, 
            filename= "20210219_ndwi.tif",
            format = "GTiff", # save as a tif, save as a FLOAT if not default, not integer
            overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.

# ndwi <- stack("20210219_ndwi.tif")

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
ggsave("20210219_hist.eps", h, device = "eps", dpi = 72)

# Analysis
# The two-peak analysis provides a more stable threshold value by visual 
# inspection; therefore, we will set the threshold to twopeak:
ndwi_threshold <- twopeak

# ------------------------------------------------------------------------------
## SEARCH ALONG THE PREDETERMINED TRANSECT

# Transect for Mutale River near Tshandama
# 246130, 7478894
# 246066, 7479006, UTM Zone 36S
x1 <- (246066)
x2 <- (246130)
y1 <- (7479006)
y2 <- (7478894)

# Determine coordinates for the transect line
ma <- (y2-y1)/(x2-x1) # this relies on UTM (the coordinates are in meters)
ra <- (0.1) # r is the step size of each point along our width
t <- sqrt(((x2-x1)^2) + ((y2-y1)^2)) # length along search transect
f <- ceiling(t/ra) # number of nodes along the transect
pointers <- array(999.999, dim = c(f,2)) # preallocation for coordinates
pointers[1,1] <- x1
pointers[1,2] <- y1
for (i in 2:f){
     a <- 1
     b <- (-2)*pointers[i-1,1]
     c <- (pointers[i-1,1]^2) - (ra^2)/((ma^2)+1)
     pointers[i,1] <- ((-b)+(sqrt((b^2)-4*a*c)))/(2*a)
     pointers[i,2] <- ((pointers[i-1,2])+(ma*((pointers[i,1])-(pointers[i-1,1]))))
}
spat <- SpatialPoints(pointers)

# Determine the NDWI value along the transect line
alng <- extract(ndwi, spat, method='simple')
vals <- data.frame(pointers,alng)
n <- ggplot(vals) +
     geom_point(aes(x=X1, y=alng)) +
     geom_hline(aes(yintercept=ndwi_threshold), color = "blue") +
     xlab("Easting") +
     ylab("NDWI") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 24)) +
     theme(axis.title = element_text(face = "plain", size = 24))
ggsave("20210219_trans.eps", n, device = "eps", dpi = 72)

# Determine pixel midpoints
RDB <- -9999 # preallocate 
LDB <- -9999
alng_per <- array(-9, dim=c(f,2)) #allocation for the midpoints
# when you reach -9 in that array, you've reached the end of the midpoints/values found
restart <- 2 #initial start for i search
for (j in 1:f){
     cnt <- 1
     for (i in restart:f) {
          if (is.na(alng[i])==FALSE) {
               if (alng[i]==alng[i-1]) { # determines if the next value is equal
                    cnt <- cnt + 1 # counts how many values there are
               } else {
                    restart <- i + 1 #to keep moving forward from the last section without causing a loop at the end of it
                    break # breaks from current for loop.
               }
          }
     }
     mp <- ((cnt*ra)/2) #ra is the spacing, and mp gives the midpoint of the current distance section
     if (is.na(alng[i-1])==FALSE) {
          if (i<(f)) {
               alng_per[j,1] <- (((i-1)*ra)-mp) 
               alng_per[j,2] <- alng[i-1] 
          } else {
               alng_per[j,1] <- ((f*ra)-mp) 
               alng_per[j,2] <- alng[i-1] 
          }
     }
     if (i>=(f)) {
          break
     } 
}

# Determine where the NDWI values cross the threshold
for (i in (2:f)){
     if (alng_per[i,2]>ndwi_threshold){
          if (alng_per[i-1,2]<ndwi_threshold){
               i1 <- alng_per[i-1,1]
               i2 <- alng_per[i,1]
               j1 <- alng_per[i-1,2]
               j2 <- alng_per[i,2]
               n <- ndwi_threshold    
               RDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
               break
          }
     }
}
for (i in 1:(f-1)){
     if (alng_per[f-i,2]>ndwi_threshold){        #expressing the index such that when i = 1, f, and when i = 2, f-1.
          if (alng_per[f-i+1,2]<ndwi_threshold){
               i1 <- alng_per[f-i+1,1]
               i2 <- alng_per[f-1,1]
               j1 <- alng_per[f-i+1,2]
               j2 <- alng_per[f-1,2]
               n <- ndwi_threshold    
               LDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
               break
          }
     }
}

# ------------------------------------------------------------------------------
# Cross-river distance
LDB-RDB # meters
