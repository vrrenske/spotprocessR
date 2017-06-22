###Renske van Raaphorst##
##27-1-2015

#NOW REAAALLYY for timelapses ;) ##

#cleanup
rm(list=ls(all=TRUE))

#some useful packages
library(ggplot2)
library(ggthemes)
library(tidyr)

#load the famous MR.Rda dataset. either the one saved after prep/plotting or after the spots are counted. 
#also load the dataset MF.Rda which contains the locations of the FtsZ-ring "spots"

######################################quick glance################################################

#order per cell/frame. just for easy looking into the dataset.
MR <- MR[order(MR$cell, MR$frame, MR$length),]
#MF <- MF[order(MF$cell, MF$frame, MF$length),]

#function to change the frame in minutes where R = the "frame" column, timeinterval=number of minutes interval per frame:
minute <- function(R,timeinterval){
  R <- (R-1)*timeinterval
  return(R)
}

#so in case of a 2 min interval:
MR$time <- minute(MR$frame,4)

#MF$time <- minute(MF$frame,2)



##################################averaging time######################################################

#mark the nth division. 
#division <- function(MR){
  MR$division <- 1
  for(n in 1:(nrow(MR)-1)){
    while(MR$length[n] == MR$length[n+1]&is.na(MR$Lcor[n+1]) &is.na(MR$Lcor[n])){
      MR <- MR[-(n+1),]
     print(n)
    }
  }

  for(n in 1:nrow(MR)){
    if(MR$cell[n] == MR$cell[n+1]){
      if(MR$length[n+1]<3/4*MR$length[n]){
       MR$division[n+1] <- MR$division[n]+1
     } 
     else{
        MR$division[n+1]<- MR$division[n]
      }
   }
   else{ 
      MR$division[n+1] <- 1
   }
  }
  MR$division <- as.numeric(MR$division)
 # return(MR)
#}

MR$num2 <- c(1:nrow(MR))

#per nth division: put the time from "0" to "100%" divided. 
#percentages <- function(MR){
  listcel <- as.integer(levels(as.factor(MR$cell)))
  MR$percent <- 0
  MR$percentL <- 0
  for (n in listcel){
    d <- c(1:max(MR[MR$cell==n,]$division))
    for(x in d){
      minL <- min(MR[MR$cell==n&MR$division==x,]$time)
      maxL <- max(MR[MR$cell==n&MR$division==x,]$time)
      minLL <- min(MR[MR$cell==n&MR$division==x,]$length)
      maxLL <- max(MR[MR$cell==n&MR$division==x,]$length)
      sub <- MR[MR$cell==n&MR$division==x,]$num2
      print(maxL)
      print(sub[1])
      for(q in sub){
        MR$percent[sub] <- (MR$time[sub]-minL)/(maxL-minL)*100
        MR$percentL[sub] <- (MR$length[sub]-minLL)/(maxLL-minLL)*100
     }
    }
  }
 # return(MR)
#}

#MR <- division(MR)
MR <- percentages(MR)

#quick plot function for any cell
cellprofile <- function(cellnum, MR){
  intdat <- MR[c("time","division","color", "pole1","pole2", "Lmid")][MR$cell==cellnum,]
  intdat <- gather(intdat, spot, position, pole1:Lmid)
  intdat$color[intdat$spot=="pole1"] <- "pole1"
  intdat$color[intdat$spot=="pole2"] <- "pole2"
  cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
  return(ggplot(intdat[!is.na(intdat$position),], aes(x=time,y=position,colour=color))  + geom_point() + theme_solarized()  + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white")) + xlab("time (min)") + ylab("length(um)") + ggtitle(paste0("cell ", cellnum)) + theme(text=element_text(size=20)) +  scale_colour_manual(values= c("GFP" = "#009E73","RFP" = "red", "pole1"="#000000", "pole2"="#000000"))) 
}

#+ geom_path(aes(linetype=division))
#example:
cellprofile(2,MR)


#so now it's possible to plot everything according to their percentage of division
cellprofile2 <- function(cellnum, MR){
  return(ggplot(MR[MR$cell==cellnum,], aes(x=percent,y=Lmid, color=division))+geom_point() + geom_point(data=MR[MR$cell==cellnum,], aes(x=percent,y=pole1, color=division)) + geom_point(data=MR[MR$cell==cellnum,], aes(x=percent,y=pole2,color=division))+ theme_solarized() + xlab("time (min)") + ylab("length(um)") + ggtitle(paste0("cell ", cellnum)) + theme(text=element_text(size=20)))
}


ggplot(MR, aes(x=percent,y=Lmid))+geom_point(colour="green") + geom_point(data=MR, aes(x=percent,y=pole1)) + geom_point(data=MR, aes(x=percent,y=pole2))+ theme_solarized() + xlab("percent of division") + ylab("length(um)") + ggtitle("plot!") + theme(text=element_text(size=20)) + geom_smooth()

#bit of a mess, so binning the data. for instance in 100 time-groups:
#functions adapted from meshtransform code
mfunTL <- function(points, b, binlist){
  means <- c()
  std <- c()
  bins <- c(1:b)*100/b
  for(q in 1:b){
    mq <- mean(points[binlist==q], na.rm=T)
    sq <- sd(points[binlist==q], na.rm=T)
    means[q] <- mq
    std[q] <- sq
  }
  return(data.frame(means,std, bins))
}

########adapted "superfun" for timelapses
superfunTL <- function(dat, bins, varname, varname2){
  dat$binned <- 0
  dat$binned <-cut(dat$percent, breaks=bins, include.lowest=TRUE, labels = 1:bins)
  BL<- summary(as.factor(dat$binned))
  print(BL[bins])
  if(missing(varname2)){
    mL <- mfunTL(dat$Lcor, bins, dat$binned)
    mL$name <- varname
  }
  else{ 
    mL <- mfunTL(dat$Lcor[dat$color=="GFP"], bins, dat$binned)
    mL$name <- varname
    mL2 <- mfunTL(dat$Lcor[dat$color=="RFP"], bins, dat$binned)
    mL2$name <- varname2
  }
  mPole <- mfunTL(dat$pole1, bins, dat$binned)
  mPole$name <- "Pole"
  mPole2 <- mfunTL(dat$pole2, bins, dat$binned)
  mPole2$name <- "Pole"
  if(missing(varname2)){
    meanframe <- rbind(mL, mPole, mPole2)
  }
  else{
    meanframe <- rbind(mL, mL2, mPole, mPole2)
  }
  #add standard error & 95% confidence interval
  meanframe$n <- BL[1:bins]
  #estimated standard error
  meanframe$SE <- meanframe$std/sqrt(meanframe$n)
  #mean + E and mean - E are the 95 interval margins
  meanframe$E <- qt(.975,df=meanframe$n-1)*meanframe$SE
  return(meanframe)
}

#so now the finale:
finalframe <- superfunTL(MR, 50, "REP")

#plotting: choose betw std and E as error bars.
plotmeans <- function(dat, errors, title, var1, var2){
  dat$errors <- dat[, errors]
return(ggplot(dat, aes(x=bins, y=means, color=name)) + geom_point() + ylim(-1.2, 1.2) + geom_errorbar(aes(ymin=means-errors, ymax=means+errors), width=0.1) + theme_minimal()  + xlab("percentage of division") + ylab("location on length-axis") + ggtitle(title))
}

#or plot it again as a densityplot

#######FUNCTIONS for joining spots and making a proper tracking plot#########################
takecell <- function(MR, cellnum, threshold1, threshold2){
  #make seperate data frame & reshuffle
  test <- MR[c("time","division", "color", "pole1","pole2", "Lmid")][MR$cell==cellnum,]
  test$Lmid[(test$pole2-(abs(test$Lmid)))<(test$pole2/10)]<- NA
  if(missing(threshold1)){
    threshold1 <- max(test$Lmid[test$color=="GFP"],na.rm=T)/1.5
  }
  if(missing(threshold2)){
    threshold2 <- max(test$Lmid[test$color=="RFP"],na.rm=T)/1.5
  }
  test <- gather(test, spot, position, pole1:Lmid)
  test$color[test$spot=="pole1"] <- "pole1"
  test$color[test$spot=="pole2"] <- "pole2"
  test <- test[order(test$color,test$time, test$position),]
  #just in case: mark spotnrs
  test$spot<-1
  test <- test[!duplicated(test),]
  test <- test[!is.na(test$position),]
  test$spot[is.na(test$position)] <- 0
  for(n in 1:(nrow(test)-1)){
    if(test$time[n+1]==test$time[n]){
        test$spot[n+1] <- test$spot[n] + 1
      } }
  #make column for notation of "GFP up" and "GFP down"
  test$color3 <- test$color
  test$num <- 1:nrow(test)
  #function testtest does notate GFP/RFP up and down, threshold is the distance between previous spot which
  #defines 2 spots as "splitted"
  test$color3 <- vers2(test, "GFP", threshold1)
  if("RFP" %in% test$color){
    test$color3 <- vers2(test, "RFP", threshold2)
  }
  #splitfun is notating the spot befor splitting and copying this one. in this way, both traces will originate from this
  #spot in a plot
  test <- splitfun(test)
  return(test)
}


########HERE's the testtest from above: spreads the two spots##################
testtest <- function(test, kleur, threshold){
  for(n in 1:(max(as.numeric(test$division)))){
    frst <- test$position[test$time==min(test$time[test$color==kleur&test$division==n])&test$color==kleur&test$spot==1]
    minfirst <- min(test$num[test$division==n&test$color==kleur])
    maxfirst <- max(test$num[test$division==n&test$color==kleur])
    print(n)
    for(x in (minfirst+1):maxfirst){
      print(abs(frst-test$position[x]))
      if(test$color3[x-1]==kleur){
         if(abs(test$position[x]-test$position[x-1])>threshold)
           test$color3[x] <- paste(test$color[x],"up")
      }
      else{
        if(abs(test$position[x]-test$position[x-1])<=threshold)
          test$color3[x] <- paste(test$color[x],"up")
      }
     }
   }
  return(test$color3)
}

#######INSERTROW function (from stackoverflow, Ali B. Friedman, http://stackoverflow.com/questions/11561856/add-new-row-to-dataframe)
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

##########FINAL FUNCTION: for the splitting moment #####################################################
splitfun <- function(test){
  numframe <- data.frame(division=0, color=0, num=0)
  for(n in 2:(nrow(test)-1)){
    print(n)
    con1 <- as.numeric(test$division[test$num==n])
    con2 <- as.numeric(test$division[test$num==(n+1)])
    if(con1==con2){
      col1 <- test$color3[test$num==n]
      col2 <- test$color3[test$num==(n+1)]
      if(col1!=col2){
        nr <- test[c("division", "color", "num")][test$num==n-1,]
        if(numframe$division[end(numframe$division)] != nr$division|numframe$color[end(numframe$color)] !=nr$color)
        numframe <- rbind(numframe, nr)
      }
    }
  }

  for(n in 2:(nrow(numframe))){
    numss <- numframe$num[n]
    line <- test[test$num==numss,]
    if(line$time!=0){
     line$color3 <- test$color3[test$num==(numss+2)]
      test <- insertRow(test, line, numss)
   }
  }
  return(test)
}

#################################so if you want to plot a cell:###################################3
celltoplot <- takecell(MR,12)
finalplot(celltoplot,12)
plotlines(finalplot(celltoplot,12))
####################Quick plot function "Finalplot"####################################


finalplot <- function(celltoplot,cellnum){
  plotpalette <- c("#009E73", "#009E73" , "#009E73","#000000", "#000000", "#D55E00", "#D55E00", "#D55E00")
  return(ggplot(celltoplot, aes(x=time, y=position, fill=color3, shape=color)) + geom_point(size=3) + scale_fill_manual(values=c("GFP"="#009E73", "GFP up" = "#009E73", "GFP down" = "#009E73", "pole1" = "#000000", "pole2"="#000000", "RFP" = "#D55E00", "RFP up" = "#D55E00", "RFP down"="#D55E00"))+ scale_shape_manual(values=c("pole1"=20, "pole2"=20, "GFP"=21, "RFP"=24)) + theme_minimal()  + theme(text=element_text(size=16,family="Arial")) + xlab("time (min)") + ylab("length(um)") + ggtitle(paste("cell ",cellnum)))
}

plotlines <- function(plot){
  return(plot + geom_path(aes(linetype=as.factor(division))))
}


vers2 <- function(test, kleur, threshold){

  for(n in 1:(max(as.numeric(test$division)))){
    #frst <- test$position[test$time==min(test$time[test$color==kleur&test$division==n])&test$color==kleur&test$spot==1]
    minfirst <- min(test$num[test$division==n&test$color==kleur])
    if(abs(threshold-abs(test$position[minfirst])) < abs(test$position[minfirst])){
      if(test$position[minfirst]>0)
        test$color3[minfirst] <- paste(test$color3[minfirst],"up")
      else
        test$color3[minfirst] <- paste(test$color3[minfirst], "down")
    }
    maxfirst <- max(test$num[test$division==n&test$color==kleur])
    print(n)
    for(x in (minfirst+1):maxfirst){
      if(abs(test$position[x]-test$position[x-1])>threshold){
        if(test$color3[x-1]==kleur){
          if(abs(threshold-abs(test$position[x])) < abs(test$position[x])){
            if(test$position[x]>0)
              test$color3[x] <- paste(kleur,"up")
            if(test$position[x]<0)
              test$color3[x] <- paste(kleur,"down")
          }
        }
        if(test$color3[x-1]==paste(kleur,"up")){
          if(abs(threshold-abs(test$position[x])) < abs(test$position[x]))
            test$color3[x] <- paste(test$color[x], "down")
        }
        if(test$color3[x-1]==paste(kleur,"down")){
          if(abs(threshold-abs(test$position[x])) < abs(test$position[x]))
            test$color3[x] <- paste(test$color[x],"up")
        }
      }
      else{
        if(test$color3[x-1]==paste(kleur,"up")){
          if(abs(threshold-abs(test$position[x])) < abs(test$position[x]))
            test$color3[x] <- paste(test$color[x],"up")
        }
        if(test$color3[x-1]==paste(kleur,"down")){
          if(abs(threshold-abs(test$position[x])) < abs(test$position[x]))
            test$color3[x] <- paste(test$color[x], "down")
        }
      }
    }
  }
  return(test$color3)
}

######################################plot location heatmap finally##############################################

densityplot <- function(plot){
  return(plot + stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE)) 
}

heatmap <- function(pdens, mp){
  return(pdens + scale_fill_gradient2(low = "#000000", mid= "#FF0000", high = "#FFFF00", midpoint = mp, space = "Lab", guide = "colourbar"))
}

heatplot <- ggplot(MR[MR$division<4&MR$color=="GFP",], aes(x=percent, y=Lmid))
heatplot <- densityplot(heatplot)
heatplot <- heatmap(heatplot, 0.010)

heatplot <- heatplot + geom_point(data=finalframe[finalframe$name=="Pole",], aes(x=bins, y=means), color="white")
