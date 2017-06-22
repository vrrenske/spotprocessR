###################using Utrack###########################
##28-10-2015 Renske van Raaphorst

###processing output into useful files

######################################################################

#might be good to clean up first:
rm(list=ls(all=TRUE))

#libraries you need:
library(tidyr)
library(sp)
library(ggplot2)
library(ggthemes)
library(Cairo)
#########HERE COMES A SHITLOAD OF FUNCTIONS###########################

##############take 1 track as a starter:

onetrack <- function(Info,t,Startrow,num, ind){
  n <- Startrow[t,]
  y <- n + num[t,]-1
  track1 <- Info[n:y,]
  track1 <- t(track1)
  track1 <- data.frame(track1)
  if(y-n == 4){
    colnames(track1) <- c("track", "split1", "split2", "split3", "split4")
  }
  if(y-n == 3){
    colnames(track1) <- c("track", "split1", "split2", "split3")
  }
  if(y-n == 2){
    colnames(track1) <- c("track", "split1", "split2")
  }
  if(y-n == 1){
    colnames(track1) <- c("track", "split1")
  }
  if(y-n == 0){
    colnames(track1) <- "track"
  }
  track1$track <- as.numeric(track1$track)
  track1 <- assigntype(track1, ncol(ind), "tr", "t")
  track1$time <- c(rep(1:ncol(ind), times=1, each=8))
  if("split1"%in% colnames(track1)){
    s1 <- data.frame(track1$split1, track1$tr, track1$time)
    colnames(s1) <- c("split1", "tr", "time")
    s1$type <- 2
    track1$split1 <- NULL
    s1 <- spread(s1, tr, split1)
  }
  if("split2"%in% colnames(track1)){
    s2 <- data.frame(track1$split2, track1$tr, track1$time)
    colnames(s2) <- c("split2", "tr", "time")
    s2$type <- 3
    s2 <- spread(s2, tr, split2)
    track1$split2 <- NULL
  }
  if("split3"%in% colnames(track1)){
    s3 <- data.frame(track1$split3, track1$tr, track1$time)
    colnames(s3) <- c("split3", "tr", "time")
    s3$type <- 4
    s3 <- spread(s3, tr, split3)
    track1$split3 <- NULL
  }
  if("split4"%in% colnames(track1)){
    s4 <- data.frame(track1$split4, track1$tr, track1$time)
    colnames(s4) <- c("split4", "tr", "time")
    s4$type <- 5
    s4 <- spread(s4, tr, split4)
    track1$split4 <- NULL
  }
  track1$type <- 1
  track1 <- spread(track1, tr, track)
  if(exists("s1"))
    track1 <- merge(track1, s1, all=T)
  if(exists("s2"))
    track1 <- merge(track1, s2, all=T)
  if(exists("s3"))
    track1 <- merge(track1, s3, all=T)
  if(exists("s4"))
    track1 <- merge(track1, s4, all=T)
  track1 <- track1[!is.na(track1$t_A),]
  track1 <- track1[order(track1$time),]
  return(track1)
}

assigntype <- function(track1, z, w, u){
  track1$exp <- paste(u,"coord_x", sep="_")
  track1$exp[c(((1:z)*8)-6)] <- paste(u, "coord_y", sep="_")
  track1$exp[c(((1:z)*8)-5)] <- paste(u, "coord_z", sep="_")
  track1$exp[c(((1:z)*8)-4)] <- paste(u, "A", sep="_")
  track1$exp[c(((1:z)*8)-3)] <- paste(u, "std_x", sep= "_")
  track1$exp[c(((1:z)*8)-2)] <- paste(u, "std_y", sep="_")
  track1$exp[c(((1:z)*8)-1)] <- paste(u, "std_z", sep="_")
  track1$exp[c(((1:z)*8))] <- paste(u, "std_A", sep="_")
  colnames(track1)[colnames(track1)=="exp"] <- w
  return(track1)
}



###############mark splitting/merging events################

##here: 
#- matrix in Matlab tracksFinal -> seqOfEvents has the splitting and merging events

#import the Matlab file using >vertcat(tracksFinal.seqOfEvents); >save('seqofevents.txt', 'tmp', '-ASCII'); and load it into the R workspace
#as "seqofevents"

#bit of manipulation:
proseq <- function(seqofevents){  
  colnames(seqofevents) <- c("time", "event", "type", "merge")
  x <- 0
  seqofevents$trackn <- 0
  for (n in 1:nrow(seqofevents)){
        if(seqofevents$type[n]==1&seqofevents$event[n]==1){x <- x+1}
        seqofevents$trackn[n] <- x   
  }
  for(n in 1:max(seqofevents$trackn)){
    if(max(seqofevents$time[seqofevents$trackn==n])-min(seqofevents$time[seqofevents$trackn==n])<max(seqofevents$time)/3){
      seqofevents <- seqofevents[seqofevents$trackn!=n,]
    }
  }
  return(seqofevents)
}

##now you have a file telling you whether a track splits/merges etc.
#- manual addition until now + use insertrow below, z = which track, y is which event - to insert rows to visualize merging/splitting
#insertRow <- function(existingDF, newrow, r,z,y) {
  #existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  #existingDF[r,] <- newrow
  #existingDF[r,]$type <- z
  #existingDF[r,]$event <- y
  #existingDF
#}

##so combining those things:
combotrack <- function(track1, s, z){
  track1$event <- NA
  if (8 %in% s$trackn){
    s1 <- s[s$trackn==z,]

    for(n in 1:nrow(s1)){
      if(s1$event[n]==1){
        if(is.na(s1$merge[n])){
          track1$event[track1$time==s1$time[n]&track1$type==s1$type[n]] <- "start"
        }
      }
      if(s1$event[n]==2){
        if(is.na(s1$merge[n])){
          track1$event[track1$time==s1$time[n]&track1$type==s1$track[n]] <- "stop"
        }
      }  
    }
  }
  return(track1)
}

##plot one track
onetrackplot <- function(track2, M, cnum){ 
  trackplot <- ggplot(M, aes(x=x0, y=y0)) + geom_polygon(color="black") + geom_path(data= track2, aes(x=t_coord_x, y=t_coord_y, colour=time, line=type, group=track, label=track)) + 
    # + geom_point(data=track2[!is.na(track2$event),], 
                           #aes(x=t_coord_x, y= t_coord_y, shape=event),  size=3, colour="#0072B2") + 
    scale_colour_gradient2(low = "#56B4E9", mid= "#F0E442", high = "#CC79A7", 
                         midpoint = 30, space = "Lab", guide = "colourbar", name="Time(s)") + 
    ggtitle(paste("Cell", cnum, sep=" ")) +
    theme_minimal() + theme(plot.title=element_text(size=21, face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=16)) +
    xlab("") + ylab("") + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y=element_blank()) +
    geom_rect(xmin=min(M$x0), xmax=(min(M$x0)+15.625), ymin=min(M$y0)-2, ymax=min(M$y0)-1.75,fill="#999999") +
    annotate("text", label="1 ?m", x=min(M$x0)+7.8125, y=min(M$y0)-2.5, size=6)
  return(trackplot)
}

##plot one track
onetrackonlyplot <- function(track2, tnum){ 
  trackplot <- ggplot(track2, aes(x=t_coord_x, y=t_coord_y)) + geom_path(aes(colour=time, line=type, group=track, label=track)) + 
    # + geom_point(data=track2[!is.na(track2$event),], 
    #aes(x=t_coord_x, y= t_coord_y, shape=event),  size=3, colour="#0072B2") + 
    scale_colour_gradient2(low = "#56B4E9", mid= "#F0E442", high = "#CC79A7", 
                           midpoint = 30, space = "Lab", guide = "colourbar", name="Time(s)") + 
    ggtitle(paste("track", tnum, sep=" ")) +
    theme_minimal() + theme(plot.title=element_text(size=21, face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=16)) +
    xlab("") + ylab("") + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y=element_blank()) +
    geom_rect(xmin=min(track2$t_coord_x), xmax=(min(track2$t_coord_x)+3.90625), ymin=min(track2$t_coord_y)-1.1, ymax=min(track2$t_coord_y)-1,fill="#999999") +
    annotate("text", label="0.25 ?m", x=min(track2$t_coord_x)+2, y=min(track2$t_coord_y)-1.4, size=6)
  return(trackplot)
}

##plot one track in complete meshfile to find the right cell:
##simplify meshfile so it can be used to draw cells in polygon
simplemesh <- function(MESH, o){
  if(o=="n")
  {MESH <- MESH[MESH$slice==1,]}
  MESH2 <- MESH[c("slice","cell", "x0","y0", "length", "Xmid", "Ymid", "angle")]
  MESH3 <- MESH[c("slice","cell", "x1","y1", "length", "Xmid", "Ymid", "max_length", "angle")]
  colnames(MESH3) <- c("slice", "cell", "x0", "y0", "length", "Xmid", "Ymid", "max_length", "angle")
  MESH3$n <- MESH3$max_length + MESH3$max_length-MESH3$length
  MESH2$n <- MESH2$length
  MESH4 <- merge(MESH2, MESH3, all=T)
  MESH4 <- MESH4[order(MESH4$slice,MESH4$cell, MESH4$n),]
  return(MESH4)
}

#function to plot either a single track or all tracks (alltrack) on all cells.
checkplot <- function(track2, MESH){
  p <- ggplot(MESH, aes(x=x0, y=(y0))) + geom_polygon(aes(group=as.factor(cell)), alpha=0.1, colour="black") +  
    geom_path(data=track2, aes(x=t_coord_x, y=t_coord_y, group=as.factor(track), colour=time), alpha=0.5) +
    scale_colour_gradient2(low = "#0072B2", mid= "#F0E442", high = "#CC79A7", 
                           midpoint = 30, space = "Lab", guide = "colourbar", name="Time(s)") +
    geom_text(data=MESH[MESH$length==0,], aes(x=Xmid, y=Ymid, label=cell),hjust=0, vjust=0)  +theme_minimal() +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y=element_blank()) + xlab("") + ylab("") +
    geom_rect(xmin=(max(MESH$x0) - 30), xmax=(max(MESH$x0) -14.375), ymin= (min(MESH$y0)+ 15), ymax=(min(MESH$y0)+ 18), fill="#999999") +
    annotate("text", label="1 ?m", x =(max(MESH$x0)-22.1875), y=(min(MESH$y0) + 10))
  return(p)
}

#calculate MSD for a single track
#MSD: average(r(t)-r(0))^2 where r(t) is position at time t and r(0) is position at time 0.

MSD <- function(combtrack){
  ndata <- nrow(combtrack)
  combtrack$t_coord_x <- combtrack$t_coord_x*0.064
  combtrack$t_coord_y <- combtrack$t_coord_y*0.064
  numberOfdeltaT <- floor(ndata/4)
  dtX <- NULL
  msd <- NULL
  distlist <- NULL
  for(dt in 1:numberOfdeltaT){
    for(n in 1:(nrow(combtrack)-dt)){
      distlist[n] <- (combtrack$t_coord_x[n]-combtrack$t_coord_x[n+dt])^2 + (combtrack$t_coord_y[n]-combtrack$t_coord_y[n+1])^2
    }
    msd$meansd[dt] <- mean(distlist)
    msd$stdv[dt] <- sd(distlist)
    msd$se[dt] <- msd$stdv[dt]/sqrt(length(distlist))
    msd$length[dt] <- length(distlist)
  }
  msd <- data.frame(msd)
  msd$deltat <- 1:nrow(msd)
  return(msd)
}

vanHove <- function(combtrack){
  ndata <- nrow(combtrack)
  combtrack$t_coord_x_u <- combtrack$t_coord_x*0.064
  combtrack$t_coord_y_u <- combtrack$t_coord_y*0.064
  numberOfdeltaT <- floor(ndata/4)
  dtX <- NULL
  msd <- NULL
  combtrack$distlist <- NA
  if(nrow(combtrack)>1){
  for(dt in 1:numberOfdeltaT){
    for(n in 1:(nrow(combtrack)-dt)){
      combtrack$distlist[n] <- sqrt((combtrack$t_coord_x_u[n]-combtrack$t_coord_x_u[n+dt])^2 + (combtrack$t_coord_y_u[n]-combtrack$t_coord_y_u[n+1])^2)
        if((combtrack$t_coord_x_u[n]-combtrack$t_coord_x_u[n+dt])<0){combtrack$distlist[n] <- combtrack$distlist[n]*-1}
    }
  }
  }
  return(combtrack)
}
  
##plot the MSD over delta Time
MSDplot<- function(tt){
  return(ggplot(tt, aes(x=deltat, y=meansd, group=as.factor(track))) + geom_point() + geom_errorbar(aes(ymin=meansd-se, ymax=meansd+se), width=0.1) + theme_minimal() + geom_smooth(se=FALSE) + geom_text(data=tt[tt$deltat==1,], aes(label=track)))
}

##############for time-spans where the cell is growing: turning co?rdinates to relative ones.
reltrack <- function(at, M){
  M$time <- M$slice
  at <- merge(at, M[c("time", "cell", "Xmid", "Ymid", "angle")], all=T)
  at <- unique(at)
  at <- at[!is.na(at$t_coord_x),]
  at <- at[order(at$cell, at$time),]
  at$Xcor <- at$t_coord_x - at$Xmid
  at$Ycor <- at$t_coord_y - at$Ymid
  at$xrot <- at$Xcor * cos(at$angle) - at$Ycor * sin(at$angle)
  at$yrot <- at$Xcor * sin(at$angle) + at$Ycor * cos(at$angle)
  at$t_coord_x <- at$xrot
  at$t_coord_y <- at$yrot
  return(at)
}

##and also the M then:
relmesh <- function(M){
  M$x0 <- M$x0_rot
  M$y0 <- M$y0_rot
  M$y1 <- M$y1_rot
  M$x1 <- M$x1_rot
  return(M)
}
##########################################HERE IS THE CODE!#########################################################################################

##FIRST GET YOUR DATA:

##set your working directory
setwd(choose.dir(default = "F:/microscope files 2015/", caption = "Choose working directory"))

##get your files:
##saved as .txt files:
#- trackedFeatureInfo as "Info"
#- trackedFeatureIndx as "ind"
#- trackStartRow as "Startrow"
#- numSegments as "num"
#- seqofevents as seqofevents
#- MESH
Infofile = file.choose()
Info <-read.table(Infofile,header=F,sep="", dec = ".")
indfile= file.choose()
ind <-read.table(indfile,header=F,sep="", dec = ".")
Startrowfile = file.choose()
Startrow <-read.table(Startrowfile,header=F,sep="", dec = ".")
numfile = file.choose()
num <-read.table(numfile,header=F,sep="", dec = ".")
seqfile = file.choose()
seqofevents <-read.table(numfile,header=F,sep="", dec = ".")
U <- readline("Do you want to add meshes? y/n")
Mfile = file.choose()
Z <- readline("Need to correct for growing cells? y/n")
interval <- as.numeric(readline("What is the time interval (sec)?"))
if(U=="y"){
  load(Mfile)
  if(Z=="y"){
    MESH <- relmesh(MESH)
  }
  MESH <- simplemesh(MESH,Z)
}


####function concatenating all tracks & all MSD's (with row telling which track it is)
for(n in 1:nrow(Startrow)){
  track <- onetrack(Info, n, Startrow, num, ind)  
  track <- combotrack(track, seqofevents, n)
  track$track <- n
  if(Z=="y"){
    track <- reltrack(track,MESH)
  }
  track$time <- track$time*interval
  track <- vanHove(track)
  MSDt <- MSD(track)
  MSDt$track <- n
  if(n==1){
    alltracks <- track
    allMSDs <- MSDt
  }
  if(n>1){
    alltracks <- rbind(track,alltracks)
    allMSDs <- rbind(MSDt, allMSDs)
  }
}

#get rid of tracks shorter than 10 frames
counts <- data.frame(table(alltracks$track))
colnames(counts) <- c("track", "Freq")
alltracks <- merge(alltracks, counts, all=T)
alltracks <- alltracks[alltracks$Freq>10,]
allMSDs <- merge(allMSDs, counts, all=T)
allMSDs <- allMSDs[allMSDs$Freq>10,]

if(U=="y"){
  alltracks$cell <- NA
  for(i in unique(MESH$cell)){
   t <- point.in.polygon(alltracks$t_coord_x, alltracks$t_coord_y, MESH$x0[MESH$cell==i], (MESH$y0[MESH$cell==i]))
    alltracks$cell[t==1] <- i
  }
  alltracks <- alltracks[!is.na(alltracks$cell),]
}

#plot the histogram of all Ddistances with the variance (=MSD of the total) displayed.
#make normal distribution prediction to superimpose
grid <- with(alltracks, seq(min(distlist, na.rm=T), max(distlist,na.rm=T), length = 100))
normpred <- data.frame(predicted = grid, density = dnorm(grid, mean(alltracks$distlist, na.rm=T), sd(alltracks$distlist, na.rm=T)))
#and plot
totalMSD <- ggplot(alltracks, aes(x=distlist)) + geom_histogram(aes(y=..density..), fill="#0072B2", alpha=0.7, binwidth =0.12) + theme_minimal() + xlab("\u0394d (um)") + annotate("text", label = paste("MSD(\u0394t=", alltracks$time[2]-alltracks$time[1],"s) = \n", round(sd(alltracks$distlist, na.rm=T)^2, digits=4), "\u03BCm/s\u00B2", sep=" "), x=0.7*min(alltracks$distlist,na.rm=T), y=1.4, size=5) + geom_line(data = normpred, aes(x = predicted, y = density), colour = "#D55E00", size=1)

#######################################################OUTPUT##############################################################################################

#Now the output. We'll make a new folder in the folder you're working on because this will be quite some files.
foldername <- readline("Name your new directory - no spaces please:")
dir.create(foldername)
curdir <- getwd()
setwd(paste(curdir,foldername,sep="/"))

#and here we go:
if(U=="y"){
ggsave(checkplot(alltracks, MESH)+coord_fixed(), file="Cells_tracks.pdf")
}
#all MSD's seperate
for(i in unique(allMSDs$track)){
  ggsave(MSDplot(allMSDs[allMSDs$track==i,]), file=paste("MSD_", i, ".pdf", sep=""), width=5, height = 5)
}

#all cells seperate
#to get rid of some oufti-processing bug I didnt solve yet:
if(U=="y"){
  MESH$max_length <- NULL
  MESH <- na.omit(MESH)

for(i in unique(MESH$cell)){
  if(i %in% alltracks$cell){
  ggsave(onetrackplot(alltracks[alltracks$cell==i,], MESH[MESH$cell==i,], i) + coord_fixed(), file=paste("Cell_", i, ".pdf", sep=""))
  name <- paste("track_cell_", i, ".pdf", sep="")
  cairo_pdf(file=name , height=6, width=3)
  print(ggplot(alltracks[alltracks$cell==i,], aes(x=distlist, y=time, color=as.factor(track))) + geom_path() + theme_minimal() +
                  xlab("\u0394d (\u03BCm)") + ylab("time (sec)") + scale_colour_colorblind(name="track") +
                   theme(legend.title=element_text(size=14)))
  dev.off()
  }
}
}

if(U=="n"){
  for(i in unique(alltracks$track)){
    ggsave(onetrackonlyplot(alltracks[alltracks$track==i,], i), file=paste("Track_", i, ".pdf", sep=""), width=(max(alltracks[alltracks$track==i,]$x0)-min(alltracks[alltracks$track==i,]$x0)), height=(max(alltracks[alltracks$track==i,]$y0)-min(alltracks[alltracks$track==i,]$y0)))
    name <- paste("track_", i, ".pdf", sep="")
    cairo_pdf(file=name , height=6, width=3)
    print(ggplot(alltracks[alltracks$track==i,], aes(x=distlist, y=time, color=as.factor(track))) + geom_path() + theme_minimal() +
            xlab("\u0394d (\u03BCm)") + ylab("time (sec)") + scale_colour_colorblind(name="track") +
            theme(legend.title=element_text(size=14)) )
    dev.off()
  }
}

#plot the MSD distribution + total final MSD
cairo_pdf(file="MSD_histogram.pdf", height=4, width=5)
totalMSD
dev.off()

save(alltracks, file="alltracks.Rda")
save(allMSDs, file="allMSDs.Rda")

##plot MSD track as calculated in Jacobs-Wagner, Cell 2014:
for(i in unique(allMSDs$deltat)){
  if(i==1){
    MSDplotframe <- data.frame(MSD = mean(allMSDs$meansd[allMSDs$deltat==1]), se = mean(allMSDs$se[allMSDs$deltat==1]), time = 1)
  }
  if(i>1){
    if(sum(table(unique(allMSDs$track[allMSDs$deltat==i])))>1)
    MSDplotframe <- rbind(MSDplotframe,c(mean(allMSDs$meansd[allMSDs$deltat==i]),mean(allMSDs$se[allMSDs$deltat==i]), i))
  }
}

MSDtrackplot <- ggplot(MSDplotframe, aes(x=time, y=MSD)) + geom_point() + geom_errorbar(aes(ymin=MSD-se, ymax=MSD+se), width=0.1) + theme_minimal() + xlab("time (s)") + ylab("MSD (\u03BCm/s\u00B2)")
cairo_pdf(file="MSD_alltrackscombined.pdf", height=5,width=6)
MSDtrackplot
dev.off()