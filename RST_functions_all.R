## RESIDENCE CALC METHOD
# This method calculates the residence time and distance of a GPS track dataset. 

residenceCalc <- function(x, y, time, radius, threshold) {
	if (!is.numeric(x) || !is.numeric(y) || !is.numeric(time)) {
		stop("The x, y, and time array must be numeric.")
	}
	
	if (length(x) != length(y) | length(time) != length(y)) {
		stop("x, y, and time arrays must be of equal length")
	}
	
	if (length(threshold) != length(radius) & length(threshold) != 1) {
		stop("The threshold array should either be of length 1, or have equal length to the radius array.")
	}
	
	if (length(threshold) == 1) {
		threshold = rep(threshold, length(radius))
	}
	
	n <- length(x)
	r <- length(radius)
	
	resid_time <- numeric(n * r)
	resid_dist <- numeric(n * r)
	
	i1 = -1;
	i2 = -1;

	result = .C("residenceTD", as.integer(n), as.double(x), as.double(y), as.double(time), as.integer(r), as.double(radius), as.double(threshold), as.integer(i1), as.integer(i2), as.double(resid_time), as.double(resid_dist))
	res_time = result[[10]]
	res_time[res_time == -1] = NA
	res_dist = result[[11]]
	res_dist[res_dist == -1] = NA
	
	res_values = data.frame("Row" = numeric(n))
	
	for (i in 0:(r - 1)) {
		times <- res_time[(i * n + 1):((i + 1) * n)]
		dists <- res_dist[(i * n + 1):((i + 1) * n)]
		res_values[paste0("scale_RT_rad_", radius[i + 1], "_th_", threshold[i + 1])] = times
		res_values[paste0("scale_RD_rad_", radius[i + 1], "_th_", threshold[i + 1])] = dists
		res_values[paste0("norm_RT_rad_", radius[i + 1], "_th_", threshold[i + 1])] = times / max(times, na.rm = T)
		res_values[paste0("norm_RD_rad_", radius[i + 1], "_th_", threshold[i + 1])] = dists / max(dists, na.rm = T)
		res_values[paste0("res_rad_", radius[i + 1], "_th_", threshold[i + 1])] = dists / max(dists, na.rm = T) - times / max(times, na.rm = T)
	}
	
	return(res_values)
}


# HAVERSINE FORMULA DISTANCE METHOD
# This method accepts the longitude and latitute of two coordinates and calculates the distance between them using the haversine formula. It returns this distance in km.
haversineDist <- function(lon1, lat1, lon2, lat2, rad = 6371) {
	lon1 <- lon1 * pi / 180
	lon2 <- lon2 * pi / 180
	lat1 <- lat1 * pi / 180
	lat2 <- lat2 * pi / 180
	
	t1 <- sin((lat2 - lat1)/2)^2
	t2 <- cos(lat1) * cos(lat2)
	t3 <- sin((lon2 - lon1)/2)^2
	
	t <- t1 + t2 * t3   # by the haversine formula, this is equal to haversine(d/r)
	frac <- 2 * asin(sqrt(t))   # compute d/r
	
	dist <- frac * rad    # final distance between the given points
	dist
}


# DISTANCE FORMULA
# Accepts grid coordinates (after a transformation) and calculates diagonal length.
projectDist <- function(x1, y1, x2, y2, scale) {
	dist = sqrt((x1 - x2)^2 + (y1 - y2)^2)
	dist * scale
}


# PLOT TRACK CONTOUR FUNCTION
# This function accepts a dataset with residence values (time, distance, and residuals) and the radius value. It plots the spatial tracks, color-coded by value of residence time and distance.

plotTrackContour <- function(tracks, radius, thresh, lwd = 1) {
	rt = paste0("norm_RT_rad_", radius, "_th_", thresh)
	rd = paste0("norm_RD_rad_", radius, "_th_", thresh)
	res = paste0("res_rad_", radius, "_th_", thresh)
	title = paste0("(radius of ", radius, ")")
	
	RT1 <- ggplot() +
		ylab("Latitude") +
		xlab("Longitude") +
		ggtitle(paste0("Track Colored by Residence Time ", title)) +
		geom_path(data = tracks, aes_string(x = "lon", y = "lat", colour = rt), lwd = lwd) +
		scale_colour_gradientn(colours = c("cyan", "blue", "navy"), name = "Time") +
		theme(axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 1.5), plot.title = element_text(vjust = 1.5, face = "bold", size = 13))
	
	RD1 <- ggplot() +
		ylab("Latitude") +
		xlab("Longitude") +
		ggtitle(paste0("Track Colored by Residence Distance ", title)) +
		geom_path(data = tracks, aes_string(x = "lon", y = "lat", colour = rd), lwd = lwd) +
		scale_colour_gradientn(colours = c("cyan", "blue", "navy"), name = "Distance") +
		theme(axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 1.5), plot.title = element_text(vjust = 1.5, face = "bold", size = 13))
		
	multiplot(RT1, RD1, cols = 2)
}


# PLOT TRACK RES FUNCTION
# This function accepts a dataset with residence values (time, distance, and residuals) and a radius value. It plots the color-coded residual values, as well as the track, color-coded by residuals.

plotTrackRes <- function(tracks, radius, thresh, time_units, ps = 1.5) {
	res = paste0("res_rad_", radius, "_th_", thresh)
	title = paste0("(radius of ", radius, ")")
	res_ind = which(colnames(tracks) == res)
	tracks = tracks[!is.na(tracks[ , res_ind]), ]
	
	## residuals colored by sign of residuals	
	res_plot <- ggplot() +
	geom_point(data = tracks[tracks[ ,res_ind]  == 0,], aes_string(x = "time_diff", y = res), size = ps, colour = "black") +
	geom_point(data = tracks[tracks[ ,res_ind] > 0,], aes_string(x = "time_diff", y = res), size = ps, colour = "blue") +
	geom_point(data = tracks[tracks[ ,res_ind]  < 0,], aes_string(x = "time_diff", y = res), size = ps, colour = "red") +
	ggtitle(paste("Residuals", title)) +
	xlab(paste0("Trip duration (", time_units, ")")) + ylab("Residuals") + 
	theme(axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 1.5), plot.title = element_text(vjust = 1.5, size = 13, face = "bold"))
	
	## spatial track colored by sign of residuals
	track_plot <- ggplot() +
	geom_path(data = tracks, aes(x = lon, y = lat), lwd = 0.5, colour = "grey") +
	geom_point(data = tracks[tracks[,res_ind] == 0,], aes(x = lon, y = lat), size = ps, colour = "black") +
	geom_point(data = tracks[tracks[,res_ind] > 0,], aes(x = lon, y = lat), size = ps, colour = "blue") +
	geom_point(data = tracks[tracks[,res_ind] < 0,], aes(x = lon, y = lat), size = ps, colour = "red") +
	ggtitle("Track, colored by sign of residuals") +
	xlab("Longitude") + ylab("Latitude")+ 
	theme(axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 1.5), plot.title = element_text(vjust = 1.5, size = 13, face = "bold"))
	
	multiplot(res_plot, track_plot, cols = 2)
}


# MULTIPLOT FUNCTION (from "Cookbook for R")
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# PLOT MULTISCALE FUNCTION
# This function accepts a dataset with residence values (time, distance, and residuals) and a vector of radius values. It plots the the percentage of Residuals for each type of point (Positive, Negative, Zeros) at each radius scale.  Radius vector must be the one that was used to create dataset. 
# This function makes multiple plots in a grid with one plot per track

plotMultiScale <- function(tracks,radius){ 
	
	require(dplyr)
	
#proportion of positive values
	Ppos<- function(x) {a<-length(subset(x,x>0)); b<-length(na.omit(x)); c<-(a/b); return(c)}
#proportion of negative values
	Pneg<- function(x) {a<-length(subset(x,x<0)); b<-length(na.omit(x)); c<-(a/b); return(c)}
#proportion of zeros
	Ptran <- function(x) {a<-length(subset(x,x==0)); b<-length(na.omit(x)); c<-(a/b); return(c)}
	
#makes vector of just res values
	subsetData<- tracks[,grep("res", colnames(tracks))] 
	res<-cbind(tracks[,1:8],subsetData) 
	a<-res%>%group_by(band)%>%summarise_each(funs(Ppos = Ppos(.),Pneg = Pneg(.),Ptran = Ptran(.)),c(9:length(res)))
	
	par(mfrow = c(3, 3))
	par(cex = 0.6)
	par(mar = c(2, 2, 0, 0), oma = c(4, 4, 0.5, 0.5))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))
	ix<-length(radius)
	for (i in 1:length(bandIDs)) {
		
		plot(radius,a[i,2:(ix+1)], col="blue", ylim=c(0, 1), pch=20); lines(radius,a[i,2:(ix+1)], col="blue")
		points(radius,a[i,(ix+2):((ix*2)+1)], col="red",pch=20); lines(radius,a[i,(ix+2):((ix*2)+1)], col="red")
		points(radius,a[i,(length(a)-ix+1):length(a)], col="black", pch=20); lines(radius,a[i,(length(a)-ix+1):length(a)], col="black")
		text((max(radius) - (max)(radius)*.15),.9,labels=bandIDs[i])
	}
	#a is the percentage of each type of point at each scale 
	return(a) 
}


# PLOT MULTISCALE FUNCTION
# This function accepts a dataset with residence values (time, distance, and residuals) and a vector of radius values. It plots the the percentage of Residuals for each type of point (Positive, Negative, Zeros) at each radius scale.  Radius vector must be the one that was used to create dataset. 
# This function make one plot and input should be values from one track

plotMultiScale1<-function(track, radius, xmax, Resti=-10){ 
	
	require(dplyr)
	
#proportion of positive values
	Ppos<- function(x) {a<-length(subset(x,x>0)); b<-length(na.omit(x)); c<-(a/b); return(c)}
#proportion of negative values
	Pneg<- function(x) {a<-length(subset(x,x<0)); b<-length(na.omit(x)); c<-(a/b); return(c)}
#proportion of zeros
	Ptran <- function(x) {a<-length(subset(x,x==0)); b<-length(na.omit(x)); c<-(a/b); return(c)}
	
#makes vector of just res values
	subsetData<- track[,grep("res", colnames(track))] 
	res<-cbind(track[,1:8],subsetData) 
	a<-res%>%group_by(band)%>%summarise_each(funs(Ppos = Ppos(.),Pneg = Pneg(.),Ptran = Ptran(.)),c(9:length(res)))
	
    ix<-length(radius)  
	
    plot(radius,a[i,2:(ix+1)], col="blue", xlim=c(0, xmax),ylim=c(0, 1), pch=20, ylab="Percentage of residuals",xlab="Radius (km)"); lines(radius,a[i,2:(ix+1)], col="blue",lwd=2)
    points(radius,a[i,(ix+2):((ix*2)+1)], col="red",pch=20); lines(radius,a[i,(ix+2):((ix*2)+1)], col="red",lwd=2)
    points(radius,a[i,(length(a)-ix+1):length(a)], col="black", pch=20); lines(radius,a[i,(length(a)-ix+1):length(a)], col="black",lwd=2)
    points(rep(Resti,3),seq(-.1,4,by=2),type = "l",col=rgb(0,0,0, alpha=.2),lwd=6)
    points(seq(-1,xmax+5,by=xmax+5),rep(0.05,2),type = "c", lty="dotted",col="black",lwd=1)
    text((xmax - xmax*.15),.9,labels=bandIDs[i])
	
	return(a) #a is the percentage of each type of point at each scale 
}

# PLOT MULTISCSCALE FUNCTION
#input from the plotMultiScale.  This function determines the scale that first falls below 0.05 to dynamically select a radius.  
chooseDynScale<-function(scales,radius,bandIDs){
	S1<-NULL
	bandIDs=unique(scales$band)
	for (i in 1:length(bandIDs)) {
		Row<-subset(scales, band==bandIDs[i])
		TransRes<-c(Row[,1],Row[,(length(scales)-length(radius)+1):length(scales)])
		index<-min(which(TransRes < 0.05))
		N<-names(TransRes)
		ll<-N[index]
		s<-as.numeric(unique(na.omit(as.numeric(unlist(strsplit(unlist(ll), "[^0-9 & .]+"))))))
		info<-cbind(bandIDs[i],s[1],s[1]/2)
		S1<-rbind(S1,info)
	}
	colnames(S1)<-cbind("id","diameter(km)","radius(km)")
	dynrad<-S1
	return(dynrad)
}

# PLOT TRACKRESDYN FUNCTION
# This function accepts a dataset with residence values (time, distance, and residuals) and a dynamically chosen radius values. It plots the color-coded residual values, as well as the track, color-coded by residuals.

plotTrackResDyn <- function(id, all_tracksD, dynscale, time_units, ps = 1) {

  tracks<-all_tracksD[all_tracksD$band == id, ]
  radius = as.numeric(dynscale[dynscale[,1]== id,3])
	res_ind = 14
	#tracks = tracks[!is.na(tracks[ , res_ind]), ] #gets rid of NA rows and warning, but good to have in to double check
	
	res_plot <- ggplot() +
    geom_point(data = tracks[tracks[ ,res_ind]  == 0,], aes_string(x = "time_diff", y = "res"), size = ps, colour = "black") +
    geom_point(data = tracks[tracks[ ,res_ind] > 0,], aes_string(x = "time_diff", y = "res"), size = ps, colour = "blue") +
    geom_point(data = tracks[tracks[ ,res_ind]  < 0,], aes_string(x = "time_diff", y = "res"), size = ps, colour = "red") +
    xlab(paste0("Trip duration (", time_units, ")")) + ylab("Residuals") + ggtitle(paste("Dynamic Scaled Radius (km)", radius)) +
    theme(axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 1.5), plot.title = element_text(vjust = 1.5, size = 13, face = "bold"))
	
## spatial track colored by sign of residuals
	track_plot <- ggplot() +
    geom_path(data = tracks, aes(x = lon, y = lat), lwd = 0.5, colour = "grey") +
    geom_point(data = tracks[tracks[,res_ind] == 0,], aes(x = lon, y = lat), size = ps, colour = "black") +
    geom_point(data = tracks[tracks[,res_ind] > 0,], aes(x = lon, y = lat), size = ps, colour = "blue") +
    geom_point(data = tracks[tracks[,res_ind] < 0,], aes(x = lon, y = lat), size = ps, colour = "red") +
    ggtitle("Track, colored by sign of residuals") +
    xlab("Longitude") + ylab("Latitude")+ 
    theme(axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 1.5), plot.title = element_text(vjust = 1.5, size = 13, face = "bold"))
	#require(gridExtra)
	#grid.arrange(res_plot, track_plot, ncol = 2)
	multiplot(res_plot, track_plot, cols = 2)
}



