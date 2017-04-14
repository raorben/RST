#########################################################
##																															##
## GENERAL CODE FOR USING THE RESIDUAL METHOD OVER A TRACKING DATASET    ##
## Irina Tolkova, November 2015																				##
##																															##
#########################################################

# --------------------- LOADING FUNCTIONS AND LIBRARIES --------------------- #

# clear workspace and set working directory
rm(list = ls())
setwd("/Users/ira/Dropbox/Irina_Tolkova/GeneralResidualMethod/final_code")

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
dataset <- data.frame(	"band" = origData$band,
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



# --------------------- CALCULATING RESIDENCE VALUES --------------------- #
# create a time array: duration of trip in given time units
time_units = "mins" # one of "secs", "mins", "hours", "days", "weeks"

# set the desired radius and threshold values
radius <- c(0.5, 1.935, 3, 5)
threshold <- rep(0, length(radius))

# determine the different individuals
bandIDs <- unique(dataset$band)

# for each track, calculate residence values
all_tracks = data.frame()

for (i in 1:length(bandIDs)) {
	subdata = dataset[dataset$band == bandIDs[i], ]
	
	subdata$time_diff = as.numeric(subdata$datetime - subdata$datetime[1], units = time_units)
	
	result <- residenceCalc(subdata$x, subdata$y, subdata$time_diff, radius, threshold)
	subdata = cbind(subdata, result)
	all_tracks = rbind(all_tracks, subdata)
}



# --------------------- PLOTTING TRACKS AND RESIDUALS --------------------- #
# The plotting functions below take the dataset with residence values and and a string naming the residual. To look at the full dataset, type "track = all_tracks". To "zoom in" on a section of the track, pass in the specific indeces.

# This function plots the track contours, color-coded by both residence time and distance:
for (i in 1:length(radius)) {
	quartz()
	#windows()#for Windows users
	plotTrackContour(all_tracks, radius[i], threshold[i])
}

# This function plots the residual values and the track contour, color-coded by sign of residual:
for (i in 1:length(radius)) {
	quartz()
	#windows()#for Windows users
	plotTrackRes(all_tracks, radius[i], threshold[i], time_units)
}

