##15-05-13## 
##Renske van Raaphorst##

##to merge GFP/RFP data, so we can see both localizations in one plot.

#cleanup
rm(list=ls(all=TRUE))

#packages
library(ggplot2)
 
library(ggthemes)
###############################################merging GFP and RFP data############################################

##load both dataframes. note: both are saved as "merged.Rda" from the dataframe "MR". so we have to load them 
#seperately and rename them.

#set your working directory as you wish using setwd()

#first load the GFP one manually:

GFPfile <- file.choose()
load(GFPfile)
GFP<-MR
rm(MR)
GFPname = basename(GFPfile)

#then RFP:
RFPfile <- file.choose()
load(RFPfile)
RFP<-MR
rm(MR)
RFPname = basename(RFPfile)

##we need to get both datasets in one. only, we need to keep them apart. that's why we need to add an extra column
#indicating which spot is gfp and which is RFP.
GFP$color <- "GFP"
RFP$color <- "RFP"
GR <- merge(GFP, RFP, all=T)

###################################################################################################################

#example simple densityplot:
allcombined <- ggplot(GR, aes(x=Lcor, colour=color)) + geom_density()

##now you might want to go back to the prep and plotting script to make the quartile partitions.
##however, we have more cells now because of the GFP & RFP. we can circumvent that by cutting 
##the four partitions by length:

#measure one, by length: 
#GR$q1 <- cut(GR$length, breaks=4, labels = 1:4)

GR <- GR[order(GR$length),]
GR$num2 <- c(1:nrow(GR))
#or two, by quartiles of the number of cells:
GR$q1 <- cut(GR$num2, breaks=4, labels = 1:4)

#or if you just want to cut it into four even parts ordered by length:
#GR <- GR[order(GR$length),]
#GR$numall <- c(1:nrow(GR))
#GR$q1 <- cut(GR$numall,breaks=4, labels = 1:4)

#prep for the code: define maxima.
xmax <- 0.5*max(GR$length, na.rm=TRUE)
ymax <- 0.5*max(GR$max.width, na.rm=TRUE)

####################################################################################################################
#now this is modification of the code of prep/plotting. it works for "q1" so choose which from above will be q1.

#Length and width corrected for average length per quartile for plotting the coordinate plots
meansL <- c(mean(GR$length[GR$q1==1], na.rm=TRUE),mean(GR$length[GR$q1==2], na.rm=TRUE), mean(GR$length[GR$q1==3], na.rm=TRUE), mean(GR$length[GR$q1==4], na.rm=TRUE))
meansW <- c(mean(GR$max.width[GR$q1==1], na.rm=TRUE), mean(GR$max.width[GR$q1==2], na.rm=TRUE), mean(GR$max.width[GR$q1==3], na.rm=TRUE), mean(GR$max.width[GR$q1==4], na.rm=TRUE))

#length
GR$Lcor <- (GR$Lmid/GR$length*meansL[1])

#width
GR$Dcor <- (GR$Dum/GR$max.width*meansW[1])

#seperate frames: 
Q1 <- GR[GR$q1==1,] 
Q2 <- GR[GR$q1==2,] 
Q3 <- GR[GR$q1==3,] 
Q4 <- GR[GR$q1==4,]

#length/width per group correction.
Q2$Lcor <- (Q2$Lmid/Q2$length*meansL[2])
Q3$Lcor <- (Q3$Lmid/Q3$length*meansL[3])
Q4$Lcor <- (Q4$Lmid/Q4$length*meansL[4])

Q2$Dcor <- (Q2$Dum/Q2$max.width*meansW[2])
Q3$Dcor <- (Q3$Dum/Q3$max.width*meansW[3])
Q4$Dcor <- (Q4$Dum/Q4$max.width*meansW[4])

#########################################################################################################################

#now plotting. 
#identify basic plots:
p1 <- ggplot(Q1, aes(x=Lcor, y=Dcor))
p2 <- ggplot(Q2, aes(x=Lcor, y=Dcor))
p3 <- ggplot(Q3, aes(x=Lcor, y=Dcor))
p4 <- ggplot(Q4, aes(x=Lcor, y=Dcor))

#optional: plot as 2dimensional densityplot.
#p1 + stat_density2d(geom="density2d", aes(colour=color, alpha=..level..), contour = TRUE) + theme_minimal()

#plotfunction for a 2-color dotplot
cdot2 <- function(plot){
  return(plot + geom_point(aes(colour=color), size=2, alpha=0.6) + theme_minimal() + xlab("Length(?m)") + ylab("Width(?m)") + scale_color_manual(values=c("GFP"="#009E73", "RFP"="#D55E00")))
}

#plotting the 4 groups
p1 <- cdot2(p1)
p2 <- cdot2(p2)
p3 <- cdot2(p3)
p4 <- cdot2(p4)

###############################################################################################################################
#adding. here's a copy of the allplot function. I modified the histograms to become 2 color density plots instead, having the same color
#codes.

allplot <- function(plot, data, xmax, ymax, empty, xqmax){
  
  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + theme(legend.position = "none") + coord_cartesian(xlim = c(-xmax, xmax), ylim = c(-ymax,ymax)) + geom_vline(xintercept=xqmax, alpha=0.4) + geom_vline(xintercept=-xqmax, alpha=0.4)
  p1hL <- ggplot(data, aes(x=Lcor, colour=color), alpha=0.6) + geom_density(aes(fill=color), alpha = 0.3) + coord_cartesian(xlim = c(-xmax, xmax)) + theme_minimal() +theme(axis.title.x = element_blank(), legend.position="none") + scale_color_manual(values=c("GFP"="#009E73", "RFP"="#D55E00")) + scale_fill_manual(values=c("GFP"="#009E73", "RFP"="#D55E00")) 
  p1hD <- ggplot(data, aes(x=Dcor, colour=color), alpha=0.6) + geom_density(aes(fill=color), alpha = 0.3) + coord_flip(xlim = c(-ymax, ymax)) + theme_minimal() + theme(axis.title.y = element_blank(), legend.position="none") + scale_color_manual(values=c("GFP"="#009E73", "RFP"="#D55E00")) + scale_fill_manual(values=c("GFP"="#009E73", "RFP"="#D55E00"))
  
  #align the plots properly before putting them together
  p1Dg <- ggplotGrob(p1D)
  p1hLg <- ggplotGrob(p1hL)
  p1hDg <- ggplotGrob(p1hD)
  
  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)
  
  #put the grids together using gridarrange
  return(arrangeGrob(p1hLg, empty, p1Dg, p1hD, ncol=2, nrow=2, widths=c(3, 1), heights=c(1, 2))) 
}

#again you need the mockup plot in the corner:
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

#plotting
p1G <- allplot(p1, Q1, xmax, ymax, empty, max(Q1$length)*0.5)
p2G <- allplot(p2, Q2, xmax, ymax, empty, max(Q2$length)*0.5)
p3G <- allplot(p3, Q3, xmax, ymax, empty, max(Q3$length)*0.5)
p4G <- allplot(p4, Q4, xmax, ymax, empty, max(Q4$length)*0.5)

#saving
#setwd("Y:/staff/fwn/MOLGEN/USERS/Morten_Renske/MK351/combined")

ggsave(arrangeGrob(p1G, p2G, p3G, p4G, ncol=1), filename=paste(GFPname, RFPname,"combined_allplots.pdf", sep="_"), width = 10, height = 30)

##################################################################################################################
#only histograms
hisfun <- function(data){
  return(ggplot(data, aes(x=Lcor, colour=color), alpha=0.6) + geom_density(aes(fill=color), alpha = 0.3) + coord_cartesian(xlim = c(-xmax, xmax)) + theme_minimal() +theme(legend.position="none") + scale_color_manual(values=c("GFP"="#009E73", "RFP"="#D55E00")) + scale_fill_manual(values=c("GFP"="#009E73", "RFP"="#D55E00"))  + xlab("Location length-axis (um from mid-point)"))
}
#make them
p1H <- hisfun(Q1)
p2H <- hisfun(Q2)
p3H <- hisfun(Q3)
p4H <- hisfun(Q4)

#saving
ggsave(arrangeGrob(p1H, p2H, p3H, p4H, ncol=1), filename=paste(GFPname, RFPname,"combinedhis_quart.pdf", sep="_"), width = 10, height = 20)

##############################################################################################################

#save the data file "GR" for further use
save(GR, file = paste(GFPname, RFPname,"dataframe.Rda"))

#when using 2 different mutants: compare length, mutant first
Lhist2 <- ggplot(GR, aes(x=length, fill=color, color = color), alpha=0.3) + geom_density(alpha=0.3) + scale_fill_manual(values=c("GFP"="#009E73", "RFP"="#D55E00"), labels=c(GFPfile, RFPfile)) + scale_colour_manual(values=c("GFP"="#009E73", "RFP"="#D55E00"), labels = c(GFPfile, RFPfile)) + theme_solarized()

ggsave(Lhist2, filename=paste("lengths",GFPname,RFPname, ".pdf"))


################################################STATISTICSSSSSS#######################################################################

##compare using a two-sample Kolmogorov-Smirnov test

##output: a table.

ksfun <- function(frame, cond, col1="GFP", col2="RFP"){
  frame$cond <- frame[,cond]
  ksp <- ks.test(frame$cond[frame$color==col1], frame$cond[frame$color==col2])[2]
  tp <- t.test(frame$cond[frame$color==col1], frame$cond[frame$color==col2])[3]
  wp <- wilcox.test(frame$cond[frame$color==col1], frame$cond[frame$color==col2], alternative="two.sided")[3]
  
  return(data.frame(c("Quartile"= unique(frame$q1), "Kolmogorov-Smirnov"= ksp, "Welch t-test"= tp,"Wilcoxin-rank sum test"=wp)))
}

##here: KS test is more powerful to detect shape differences, while Wilcox look at differences in mean & median.
##NB t-test only works when the distribution is nomewhat normal. so not with FtsZ-3 peaks stuffs.
##also it's not the best 
u <- rbind(ksfun(Q1, "Lcor"), ksfun(Q2, "Lcor"), ksfun(Q3, "Lcor"), ksfun(Q4, "Lcor"))
y <- rbind(ksfun(Q1, "Dcor"), ksfun(Q2, "Dcor"), ksfun(Q3, "Dcor"), ksfun(Q4, "Lcor"))
save(u, file=paste("statistics_", GFPname,"_", RFPname, ".Rda", sep=""))

