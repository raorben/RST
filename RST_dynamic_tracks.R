#########################################################
##																															##
## DYNAMIC SCALING CODE FOR USING THE RESIDUAL METHOD OVER A TRACKING DATASET    ##
## Modified from Irina Tolkova, November 2015																				##
## Rachael Orben, November 2015																															##
#########################################################

# --------------------- LOADING FUNCTIONS AND LIBRARIES --------------------- #

# clear workspace and set working directory
rm(list = ls())
#setwd("C:\\Users\\leigh.torres\\Dropbox (HMSC - OSU)\\Data\\Campbell\\Grey-Heads\\Irina_Tolkova\\GeneralResidualMethod\\final_code")
setwd("/Users/rachaelorben/Dropbox/Research/RST")

# load C function, secondary functions, plotting functions
dyn.load("RST_residenceRadii.so") # if using a mac
#dyn.load("RST_residenceRadii.dll") # if using Windows
source("RST_functions_all.R")

# load libraries
library(ggplot2) # for plotting with ggplot
library(mapproj) # for calculating map projections

# --------------------- PREPARING DATASET --------------------- #
# read in csv data file
origData <- read.csv("GHAL_23059.csv")

# construct a dataset for analysis (take the band identifiers, longitude, latitude, and time columns from the original data)
# the unique identifier "band" needs to be numeric
dataset <- data.frame("band" = origData$band,
                       "lat" = origData$lat,
                       "lon" = origData$lon,
                       "datetime" = strptime(origData$Datetime_GMT, format = "%d/%m/%Y %H:%M:%S"))

# remove NA lon/lat rows
dataset <- dataset[!is.na(dataset$lon) & !is.na(dataset$lat), ]


# create grid x- and y- coordinates from the longitude/latitude coordinates
library(mapproj)
lambert <- mapproject(dataset$lon, dataset$lat, projection = "lambert", parameters = c(mean(dataset$lon), mean(dataset$lat)))
scale <- haversineDist(min(dataset$lon), min(dataset$lat), max(dataset$lon), max(dataset$lat)) / projectDist(min(lambert$x), min(lambert$y), max(lambert$x), max(lambert$y), 1)
dataset$x <- lambert$x * scale
dataset$y <- lambert$y * scale
plot(dataset$x,dataset$y,'l')

# --------------------- CALCULATING RESIDENCE VALUES --------------------- #
# create a time array: duration of trip in given time units
time_units = "mins" # one of "secs", "mins", "hours", "days", "weeks"

# set the desired radius (km) and threshold values
# make sure this is in sequential order with no duplicate values
# radius should be chosen around the expected movement scale. Its is also helpful to increase the resolution of radii increments
# near where the proportion of transit points falls below 0.05. You will probably need to run this a couple of times to refine.
radius = c(0.001, 0.01, 0.05, seq(.1,.9,by=.1), seq(1,5,by=.1)) 
threshold <- rep(0, length(radius))

# determine the different individuals
bandIDs <- unique(dataset$band)

# for each track, calculate residence values
# If a track never leaves the chosen radius this will result in NAs for residence time, distance and residuals and an associated warning.

# For residence and residual calculations, NA are assigned to 
# (1) locations at the beginning of tracks until the animal moves beyond R from the initial point, 
# and (2) to those locations at the end that are all within R of the last radius built.

all_tracks = data.frame()

for (i in 1:length(bandIDs)) {
  subdata = dataset[dataset$band == bandIDs[i], ]
  
  subdata$time_diff = as.numeric(subdata$datetime - subdata$datetime[1], units = time_units)
  
  result <- residenceCalc(subdata$x, subdata$y, subdata$time_diff, radius, threshold)
  subdata = cbind(subdata, result)
  all_tracks = rbind(all_tracks, subdata)
}


# Save results of all scales
MoveResid<-list(all_tracks, radius, threshold, time_units, bandIDs)
saveRDS(MoveResid,file = ("GHAL_23059_allscales.rda"))
#MoveResid<-readRDS("GHAL_23059_allscales.rda")
#all_tracks<-MoveResid[[1]]; radius<-MoveResid[[2]]; threshold<-MoveResid[[3]]
#time_units<-MoveResid[[4]]; bandIDs<-MoveResid[[5]]
#rm(MoveResid)


# --------------------- PLOTTING SCALE OF RESIDUALS & CHOOSING BEST ONE--------------------- #
#Resti= estimated prefered radius value 
library(dplyr)
pdf("GHAL_23059_MultiScale.pdf")
scales<-plotMultiScale1(all_tracks, radius, xmax = max(radius), Resti = 1.9)
dev.off()

#for multiple tracks
#pdf("MultiScale.pdf")
#scalesum<-plotMultiScale(all_tracks, radius) 
#dev.off()

#if you don't choose radi large enough for the % transit points to be <0.05 then chooseDynScale will return an error
dynscale<-chooseDynScale(scales, radius) #animal/track ids must be numeric
write.csv(x = dynscale,file = "DynScaleResults.csv")
detach("package:dplyr", unload=TRUE)

# --------------------- CALCULATING DYNAMICALLY SCALED RESIDENCE VALUES  --------------------- #

all_tracksD = data.frame()

# for each track, calculate residence values
for (i in 1:length(bandIDs)) {
  subdata = dataset[dataset$band == bandIDs[i], ]
  print(bandIDs[i])
  
  #dynamically choose radius scale using ouput from plotMultiScale & DynScale
  radius = as.numeric(dynscale[dynscale[,1]== bandIDs[i],3])
  threshold = 0 ##MODIFY (optional)
  
  # manipulate time to measure time passed from beginning of trip, in minutes (for residence computation)
  ref_time <- strptime(subdata$datetime[1], "%Y-%m-%d %H:%M:%S")
  subdata_time <- strptime(subdata$datetime, "%Y-%m-%d %H:%M:%S")
  subdata$time_diff <- as.numeric(subdata_time - ref_time, units = time_units)
  
  subdata$dist_diff = numeric(nrow(subdata))
  for (j in 2:nrow(subdata)) {
    subdata$dist_diff[j] = subdata$dist_diff[j - 1] + projectDist(subdata$x[j], subdata$y[j], subdata$x[j - 1], subdata$y[j - 1], 1)
  }
  
  # call residence metric function and obtain residence distance and time	
  result <- residenceCalc(subdata$x, subdata$y, subdata$time_diff, radius, threshold)
  colnames(result)<-cbind("rows","RT","RD","nRT","nRD","res") #Use for choosing radius dynamically -rename column names in result so scale independent
  subdata = cbind(subdata, result)
  
  all_tracksD = rbind(all_tracksD, subdata)
}

# Dynamic Scaling for Radius, save results
MoveResid_DynScale<-list(all_tracksD, dynscale, threshold, time_units, bandIDs)
saveRDS(MoveResid_DynScale,file = paste0("GHAL_DynScale.rda"))
#MoveResid_DynScale<-readRDS(paste0(dir,"GHAL_DynScale.rda"))
#all_tracksD<-MoveResid_DynScale[[1]]; dynscale<-MoveResid_DynScale[[2]]; threshold<-MoveResid_DynScale[[3]]
#time_units<-MoveResid_DynScale[[4]]; bandIDs<-MoveResid_DynScale[[5]]
#rm(MoveResid_DynScale)


# --------------------- PLOTTING DYMAICALLY SCALED RESIDUALS --------------------- #
# make sure to detach dplyr before running (that might help if you get an error). 
# If not the other fix is to restart R at this point and then reload the plotting functions and run the codes.  
# This function will return a warning from geom_point when the residual values are equal to NA.
# This occurs at the beginning and end of a track when points stay within the chosen radius (usually 3 rows).
# So expect "Removed 3 rows containing missing values (geom_point)" repeated 6 times.  
# If the function removes 6+ rows, check the beginng and end of tracks to see what data is 
# not included in the RST calculations.

for (i in 1:length(bandIDs)) {
  id<-bandIDs[i]
  quartz(width=16, height=8) #use for saving plot - Mac only
  #windows()#for Windows users
  plotTrackResDyn(id, all_tracksD, dynscale, time_units, ps=2)
  quartz.save(file=(paste0('GHAL_',radius,'.png')), type = "png", device = dev.cur(), dpi = 100)
  #savePlot(file=(paste0('GHAL_',radius,'.png')), type = "png", device = dev.cur()) #for Windows users
  dev.off()
}
