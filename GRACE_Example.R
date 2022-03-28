#Sophia Bakar
# Limpopo Resilience Lab Technical Training
# March 27th, 2022
#GRACE NetCDF Data Extraction
# Necessary Packages
library(tidyr)
library(ggplot2)
library(dplyr)
library(ncdf4.helpers)
library(lubridate)
library(ncdf4) #https://www.rdocumentation.org/packages/ncdf4/versions/1.10/topics/ncdf4

example_data <- nc_open("Grace_ExampleDataset.nc") # entire global dataset


lon <-ncvar_get(example_data, "lon") # longitude values
lat <-ncvar_get(example_data, "lat") # latitude values
Days_Since_2002_01_01<- ncvar_get(example_data,"time") # time
WET <- ncvar_get(example_data, "lwe_thickness") # water equivalent thickness
##Generate dates from days since 2002_01_01
D <- as.Date(Days_Since_2002_01_01, origin = '2002-01-01') # create date column from time
Date <- as_date(parse_date_time(D,"ymd")) 

#Create new data frame with only WET values from lon(E): 29-30 (cell indices: 59-61), lat(S) (-22) - (-23) (cell index _____)
WET1 <- WET[59,134,]
WET2 <- WET[60,135,]
WET3 <- WET[61, 136,]
WET <- data.frame(Date,WET1, WET2,WET3) # 3 mascons that lie within the limpopo river basin in South Africa
WET$mean <- rowMeans(WET_SM[,c(2,3,4)], na.rm = TRUE) # average 3 mascon values

ggplot(WET, aes(x = Date, y = mean)) + # plot data
  geom_point()+
  geom_smooth(method=lm) + # add linear trendline
  ggtitle("GRACE Time Series Limpopo River Basin") +
  xlab ("Time(monthly)") + ylab("Water Equivalent Thickness (cm)")



write.csv(WET_SM,file="SMGRACE.csv",row.names = TRUE) # create csv with water equivalent thickness data

