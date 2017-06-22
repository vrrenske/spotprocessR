####################COMBINING DATA UTRACK###############################################

##Renske van Raaphorst 25-04-2016

##To combine outputs from different movies; same experiment; into 1 plot and/or combine different conditions to produce comparison plots
rm(list=ls(all=TRUE))
library(ggplot2)
library(ggthemes)
##FUNCTIONS##################################################################################

#####load object in seperate environment so you can name it yourself#########################
load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

##load files from Utrack:
pickmovies <- function(C, Z){
  for(X in 1:C){
    MSDfile <- choose.files(caption=paste("Condition ", Z, "; file ", X, sep=""))
    allMSDs <- load_obj(MSDfile)
    allMSDs$mov <- X
    allMSDs$track <- as.numeric(allMSDs$track)
    if(X==1){
      totMSD <- allMSDs
    }
    if(X>1){
      allMSDs$track <- allMSDs$track + max(totMSD$track)
      totMSD <- rbind(totMSD, allMSDs)
    }
  }
  return(totMSD)
}


########LM display function from: http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph
##by Jayden

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

###########estimate Diffusion coefficient from histogram
D_K <- function(allcontracks, par){
  ##make histogram
  for(n in unique(allcontracks$time[!is.na(allcontracks$time)&allcontracks$con==par])){
    histdat <- hist(allcontracks$distlist[allcontracks$con==par&allcontracks$time==n], breaks=50, plot=FALSE)
    histdat <- data.frame(x=histdat$mids, y=histdat$density)
    ##make formula and fit the formula to the histogram
    f <- function(x, k=kvalue, D=Dvalue){k*exp(-x^2/(4*D))/((pi*4*D)^(-1/2))}
    fit <- nls(y~f(x, k, D), data=histdat, start=list(k=30, D=0.003))
    ##get D and K out of the fit
    kvalue <- as.numeric(summary(fit)$parameters[,"Estimate"][1])
    Dvalue <- as.numeric(summary(fit)$parameters[,"Estimate"][2])/n
    kse <- as.numeric(summary(fit)$parameters[,"Std. Error"][1])
    Dse <- as.numeric(summary(fit)$parameters[,"Std. Error"][2])
    hisframe <- data.frame(k=kvalue, D=Dvalue, kse=kse, Dse=Dse, con=par, time=n)
    if(n==(min(unique(allcontracks$time),na.rm=T))){allhis <- hisframe}
    if(n>(min(unique(allcontracks$time),na.rm=T))){allhis <- rbind(allhis, hisframe)}
  }
  return(data.frame(k = mean(kvalue), D=mean(Dvalue), kse=mean(kse), Dse=mean(Dse), con=par))
}

#plot the histogram of all Ddistances with the variance (=MSD of the total) displayed.
#make normal distribution prediction to superimpose and plot
totalMSDplot<-function(alltracks, sumframe){
  sumframe$con <- as.character(sumframe$con)
  plist <- list()
  Z <- 0
  for(X in 1:nrow(sumframe)){
    #make normal distribution prediction to superimpose
    grid <- with(allcontracks[allcontracks$con==sumframe$con[X],], seq(min(distlist, na.rm=T), max(distlist,na.rm=T), length = 100))
    normpred <- data.frame(predicted = grid, density = dnorm(grid, mean(allcontracks$distlist[allcontracks$con==sumframe$con[X]], na.rm=T), sd(allcontracks$distlist[allcontracks$con==sumframe$con[X]], na.rm=T)))
    if(max(normpred$density)>Z){Z<-max(normpred$density)}
    plot <- ggplot(alltracks[alltracks$con==sumframe$con[X],], aes(x=distlist)) + geom_histogram(aes(y=..density..), fill="#0072B2", alpha=0.7, binwidth =0.05) + theme_minimal() + xlab("\u0394d (um)") + 
      ggtitle(paste("D = ", round(sumframe$D[X], digits=4), "+/-", round(sumframe$Dse[X],digits=4), "\u03BCm/s\u00B2\nMSD(\u0394t=", alltracks$time[2]-alltracks$time[1], "s) = ", round(sumframe$MSD_vanHove[X], digits=4), "\u03BCm\u00B2", sep="")) + 
      geom_line(data = normpred, aes(x = predicted, y = density), colour = "#D55E00", size=1) + 
      #coord_fixed() + 
      xlim(c(-1,1)) + ylim(c(0,(Z+1.5)))
    plist[[X]] <- plot
  }
  return(plist)
}

##plot MSD track as calculated in Jacobs-Wagner, Cell 2014:
MSDtrackplot <- function(allMSDs){
  Z <- table(allMSDs$track)
  Z <- data.frame(Z)
  colnames(Z) <- c("track", "tracklength")
  allMSDs <- merge(Z, allMSDs)
  allMSDs <- allMSDs[allMSDs$tracklength>9,]
  for(i in unique(allMSDs$deltat[!is.na(allMSDs$deltat)])){
    if(i==1){
      MSDplotframe <- data.frame(MSD = mean(allMSDs$meansd[allMSDs$deltat==1]), se = mean(allMSDs$se[allMSDs$deltat==1]), time = 1)
    }
    if(i>1){
      if(sum(table(unique(allMSDs$track[allMSDs$deltat==i])))>1)
        MSDplotframe <- rbind(MSDplotframe,c(mean(allMSDs$meansd[allMSDs$deltat==i]),mean(allMSDs$se[allMSDs$deltat==i]), i))
    }
  }
  return(MSDplotframe)
}

##Plot all MSD combotracks in one plot (if you have more than one ;) )
MSDfinplotframe <- function(allconMSD, U=U){
  listcon <- unique(allconMSD$con)
  for(X in 1:U){
    MSDplotframe <- MSDtrackplot(allconMSD[allconMSD$con==listcon[X],])
    MSDplotframe$con <- listcon[X]
    if(X==1){finalplotframe <- MSDplotframe}
    if(X>1){finalplotframe <- merge(MSDplotframe, finalplotframe,all=TRUE)}
  }
  return(finalplotframe)
}
#################Summary of MSDs ########################################################################################
MSDsummary <- function(allcontracks, allconMSD, U=U){
  listM <- c()
  listW <- c()
  listS <- c()
  listcon <- unique(allcontracks$con)
  for(X in 1:U){
    listM[X] <- sd(allcontracks$distlist[allcontracks$con==listcon[X]],na.rm=T)^2
    listW[X] <- median(allconMSD$meansd[allconMSD$con==listcon[X]])
    listS[X] <- median(allconMSD$se[allconMSD$con==listcon[X]])
    if(X==1){DKs <- D_K(allcontracks,listcon[X])}
    if(X>1){DKs <- rbind(DKs, D_K(allcontracks, listcon[X]))}
  }
  MSDsummary <- data.frame(con=unique(allcontracks$con), MSD_vanHove = listM, MSD_Wagner = listW, SE_Wagner=listS)
  MSDsummary <- merge(MSDsummary, DKs, all=T)
  return(MSDsummary)
}  

############plot together:#######################
# Multiple plot function
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
###################### code ##################################################################

##set your working directory
setwd(choose.dir(default = "F:/microscope files 2015/", caption = "Choose working directory"))

##First question: how many conditions; then how many files for a single condition
U <- readline("How many conditions?: ")


##Per condition: load the files into one data frame:
for(X in 1:U){
  C <- readline(paste("How many files for condition ", X, "?: ", sep=""))
  totMSD <- pickmovies(C, X)
  tottracks <- pickmovies(C, X)
  condition <- readline(paste("Name condition ", X, ": "))
  totMSD$con <- condition
  tottracks$con <- condition
  if(X==1){
    allconMSD <-totMSD
    allcontracks <- tottracks
  }
  if(X>1){
    allconMSD <- rbind(allconMSD, totMSD)
    allcontracks <- rbind(allcontracks, tottracks)
  }
} 

##Make summary table
MSDsum <- MSDsummary(allcontracks, allconMSD, U)

##Make plot of all histograms
cairo_pdf("MSDhists.pdf")
multiplot(plotlist=totalMSDplot(allcontracks, MSDsum), cols=3)
dev.off()

##and of all tracks
MSDplotframe <- MSDfinplotframe(allconMSD, U)
MSDalltracksoneplot <- ggplot(MSDplotframe, aes(x=time, y=MSD, color=con)) + geom_point() + geom_errorbar(aes(ymin=MSD-se, ymax=MSD+se), width=0.1) + theme_minimal() + xlab("time (s)") + ylab("MSD (\u03BCm\u00B2)") + geom_smooth(data=MSDplotframe[MSDplotframe$time<5,],method="lm") + 
  scale_colour_brewer(palette="Set2") + scale_x_continuous(limits=(c(0, max(MSDplotframe$time)*.67))) + scale_y_continuous(limits=c(0, (max(MSDplotframe$MSD[MSDplotframe$time<max(MSDplotframe$time*0.67)])+0.01)))

cairo_pdf("alltracks.pdf")
print(MSDalltracksoneplot)
dev.off()



###################################################################################

##plot fit with the original histogram
#ggplot(allcontracks[allcontracks$con=="O",], aes(x=distlist)) + geom_histogram(aes(y=..density..), bins=50, colour="#56B4E9", fill="#56B4E9", alpha=0.8) + stat_function(fun=f, colour="#E69F00", size=1) + theme_minimal() + xlab("displacements (um)") + ggtitle(paste("D = ", round(Dvalue, digits=4), "+/-", round(Dse,digits=4), "\u03BCm/s\u00B2\nMSD = ", round(MSDsum$MSD_vanHove[MSDsum$con=="O"], digits=4), "\u03BCm\u00B2", sep=""))

              
