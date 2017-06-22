##9-1-2015
##Renske van Raaphorst

#preparation of dataset containing spots for plotting.

#might be good to clean up first:
rm(list=ls(all=TRUE))

#before running the script, open datasets:
#"REP": output of peakfitter. 
#"M": the excel output file of the meshes used for peakfitter 
#aand start
library(ggplot2)
library(gridExtra)
library(scales)
library(ggthemes)
library(MASS)
library(rJava)
library(xlsx)

##set your working directory
setwd(choose.dir(default = "Y:/staff/fwn/MOLGEN/USERS/Morten_Renske", caption = "Choose working directory"))

##### GFP or RFP file ##########
fpfile = file.choose()
#if(substr(fpfile, start=(nchar(fpfile)-3), stop=(nchar(fpfile))) == ".xls")
  #REP <- read.xlsx(fpfile, sheetName="Sheet1")
#if(substr(fpfile, start=(nchar(fpfile)-3), stop=(nchar(fpfile))) == ".txt")
  REP <-read.table(fpfile,header=T,sep="\t", dec = ".")
name=basename(fpfile)

########### MESH file #############
cellfile = file.choose()
if(substr(cellfile, start=(nchar(cellfile)-3), stop=(nchar(cellfile))) == ".xls")
  M <- read.xlsx(cellfile, sheetName="Sheet1")
if(substr(cellfile, start=(nchar(cellfile)-3), stop=(nchar(cellfile))) == ".txt")
  M <- read.table(cellfile,header=T,sep="\t", dec=",")
cellname=basename(cellfile)

#!! --> if using .txt extension: make sure the "length" cells are notated as numeric, with the same amount of decimal points
#       for both data frames before importing.
#!! --> also check your decimal seperator for the tab delimited .txt files.

###############################################################################################
#wat basisplotfuncties 
densityplot <- function(plot){
  return(plot + stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE)) 
}

heatmap <- function(pdens, mp){
  return(pdens + scale_fill_gradient2(low = "#000000", mid= "#FF0000", high = "#FFFF00", midpoint = mp, space = "Lab", guide = "colourbar"))
}

#function for goodlooking x/y coordinate plot:
#makes a plot sized as the max cell (width/length: xmax) of the quartile inside a plot 
#which has a grey background as large as the largest cell in the dataset
#so all quartile plots will have the same dimensions. 
#the title, y axis and x axis will also be drawn.
coplot <- function(pheat, xmax, ymax, xqmax){
 return(pheat + xlab("Length (?m)") + ylab("Width (?m)") + coord_cartesian(xlim = c(-xmax,xmax), ylim=c(-ymax,ymax)) + geom_vline(xintercept = xqmax) + geom_vline(xintercept=-xqmax) + geom_hline(yintercept = ymax) + geom_hline(yintercept = -ymax) + theme(panel.background = element_rect(fill = "dark grey"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}
  
#############################OPTIONAL: SPOTFINDER TRANSLATION######################################
#merge dataset M with dataset REP
#might give problems because of differences in digits, so let's change that
specify_decimal <- function(x, k) as.numeric(format(round(x, k), nsmall=k))

oijen <- function(dat){
  dat$d <- dat$D
  dat$l <- dat$L
  dat$rel.l <- dat$L_normalized
  dat$frame <- dat$slice
  dat <-dat[order(dat$slice, dat$cell,dat$length),]
  dat$spot <- 1
  for(n in 1:(nrow(dat)-1)){
    if(dat$length[n+1]==dat$length[n]){
      dat$spot[n+1] <- dat$spot[n] + 1
    } }
  
  return(dat)
}

decimals <- function(dat, num){ 
  dat$length <- specify_decimal(dat$length, num)
  dat$area<- specify_decimal(dat$area, num)
  dat$volume<- specify_decimal(dat$volume, num)
  return(dat)
}


REP <- oijen(REP)
REP <- decimals(REP,1)
M <- decimals(M,1)


#######################START CODE##################################################################

#remove all spots where L (y-axis location of the spot) > length (total cell length) and all cells
#where L < 0. these are spots outside of the cell.

REP<- REP[(0<REP$l),]
REP <-REP[(REP$rel.l<1),]

#summary(REP)

M <- M[order(M$length),]
M$cellnum <- c(1:nrow(M))

#merging
MR <- merge(M,REP, all=T)

#NA in spots replaced by "0"
MR$spot[is.na(MR$spot)] <- 0

#remove MR's cells which have NA's in the cell area.
MR <- MR[!is.na(MR$length),]

#if needed: remove smallest and largest ones (cutoff: smaller than 1/2 mean and larger than 2x mean)
MR <- MR[MR$length<(2*mean(MR$length)),]
MR <- MR[MR$length>(0.5*mean(MR$length)),]

MR <- MR[order(MR$length),]

#make column with row numbers per cell length. 
MR$num <- c(1:nrow(MR))

#summary(MR)


################################################################################################
#plot preparation
#quartiles, maxima, etc.
MR$Lmid<-(MR$l-0.5*MR$length)*0.064
MR$pole1<- -MR$length*0.032
MR$pole2<- -MR$pole1
MR$Dum <- MR$d*0.064

##make quartile partitions
#measure one, by length: 
#MR$q1 <- cut(MR$length, breaks=4, labels = 1:4)

#or two, by quartiles of the number of cells:
MR$q1 <- cut(MR$cellnum, breaks=4, labels = 1:4)

MR$length <- MR$length*0.064
MR$max.width <- MR$max.width*0.064

xmax <- 0.5*max(MR$length, na.rm=TRUE)
ymax <- 0.5*max(MR$max.width, na.rm=TRUE)


#Length and width corrected for average length per quartile for plotting the coordinate plots
meansL <- c(mean(MR$length[MR$q1==1], na.rm=TRUE),mean(MR$length[MR$q1==2], na.rm=TRUE), mean(MR$length[MR$q1==3], na.rm=TRUE), mean(MR$length[MR$q1==4], na.rm=TRUE))
meansW <- c(mean(MR$max.width[MR$q1==1], na.rm=TRUE), mean(MR$max.width[MR$q1==2], na.rm=TRUE), mean(MR$max.width[MR$q1==3], na.rm=TRUE), mean(MR$max.width[MR$q1==4], na.rm=TRUE))

#length
MR$Lcor <- (MR$Lmid/MR$length*meansL[1])
#width
MR$Dcor <- (MR$Dum/MR$max.width*meansW[1])

#seperate frames: 
Q1 <- MR[MR$q1==1,] 
Q2 <- MR[MR$q1==2,] 
Q3 <- MR[MR$q1==3,] 
Q4 <- MR[MR$q1==4,]

#length
Q2$Lcor <- (Q2$Lmid/Q2$length*meansL[2])
Q3$Lcor <- (Q3$Lmid/Q3$length*meansL[3])
Q4$Lcor <- (Q4$Lmid/Q4$length*meansL[4])

#width
Q2$Dcor <- (Q2$Dum/Q2$max.width*meansW[2])
Q3$Dcor <- (Q3$Dum/Q3$max.width*meansW[3])
Q4$Dcor <- (Q4$Dum/Q4$max.width*meansW[4])

###############################################################################################################################################
##plotting! -> coordinate plots
p1 <- ggplot(Q1, aes(x=Lcor, y=Dcor))
p1 <- densityplot(p1)
p1 <- coplot(p1, xmax, ymax, max(Q1$length,na.rm=T)*0.5)


p2 <- ggplot(Q2, aes(x=Lcor, y=Dcor))
p2 <- densityplot(p2)
p2 <- coplot(p2, xmax, ymax, max(Q2$length,na.rm=T)*0.5)

p3 <- ggplot(Q3, aes(x=Lcor, y=Dcor))
p3 <- densityplot(p3)
p3 <- coplot(p3, xmax, ymax, max(Q3$length,na.rm=T)*0.5)

p4 <- ggplot(Q4, aes(x=Lcor, y=Dcor))
p4 <- densityplot(p4)
p4 <- coplot(p4, xmax, ymax, max(Q4$length, na.rm=TRUE)*0.5)

#plotting! -> L and D ordered by cell length
pL <- ggplot(MR, aes(x=num, y=Lmid))
pLpoint <- pL + geom_point() + ggtitle("Spot location on length axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("Y-position (?m)") + theme_bw()
pLD <- densityplot(pL) + ggtitle("Spot location on length axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("Y-position (?m)") + geom_line(data=MR, aes(x=num,y=pole1),colour="white") + geom_line(data=MR, aes(x=num,y=pole2),colour="white")

pW <- ggplot(MR, aes(x=num, y=Dum))
pWpoint <- pW + geom_point() + ggtitle("Spot location on width axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("X-position (?m)") + theme_bw()
pWD <- densityplot(pW) + ggtitle("Spot location on width axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("X-position (?m)") + geom_hline(yintercept=ymax) + geom_hline(yintercept=-ymax) + coord_cartesian(ylim=c(-ymax,ymax))

#make heatmap using the half max densities:

mp1 <- kde2d(Q1$Lmid[!is.na(Q1$Lmid)], Q1$Dum[!is.na(Q1$Dum)])
mp2 <- kde2d(Q2$Lmid[!is.na(Q2$Lmid)], Q2$Dum[!is.na(Q2$Dum)])
mp3 <- kde2d(Q3$Lmid[!is.na(Q3$Lmid)], Q3$Dum[!is.na(Q3$Dum)])
mp4 <- kde2d(Q4$Lmid[!is.na(Q4$Lmid)], Q4$Dum[!is.na(Q4$Dum)])
mplist <- c(median(range(mp1$z)), median(range(mp2$z)), median(range(mp2$z)),median(range(mp2$z)))

mp <- max(mplist)
p1 <- heatmap(p1, mp)
p2 <- heatmap(p2, mp)
p3 <- heatmap(p3, mp)
p4 <- heatmap(p4, mp)

mpL <- kde2d(MR$num[!is.na(MR$Lmid)], MR$Lmid[!is.na(MR$Lmid)])
mpL1 <- median(range(mpL$z))

pLD <- heatmap(pLD, mpL1)

mpW <- kde2d(MR$num[!is.na(MR$Dum)], MR$Dum[!is.na(MR$Dum)])
mpW1 <- median(range(mpW$z))

pWD <- heatmap(pWD, mpW1)

###########################MESH INCORPORATION!####################################################
#Take output of meshtransform file. this way you only have to transform it once per GFP/RFP combo

#tmp_env <- new.env()
#load(file.choose(), tmp_env)
#MESH<- get(ls(tmp_env), envir=tmp_env)

######################### mid-points ################################################################################
#mfun <- function(points, b, binlist){
  #means <- c()
  #for(q in 1:b){
   # mq <- mean(points[binlist==q])
    #means[q] <- mq
 # }
 # return(means)
#}

#superfun <- function(dat, bins){
  #dat$av <- 0
  #dat$av <- dat$y0_rot/dat$max_length*100
  
  #cutpoints<-quantile(dat$av,(0:bins)/bins)
  #dat$binned <-cut(dat$av,cutpoints, include.lowest=TRUE, labels = 1:bins)
  
  #x0means <- mfun(dat$x0_rot, bins, dat$binned)
  #x1means <- mfun(dat$x1_rot, bins, dat$binned)
  #y0means <- mfun(dat$y0_rot, bins, dat$binned)
  #y1means <- mfun(dat$y1_rot, bins, dat$binned)
  
  #meanframe <- data.frame(x0means, x1means, y0means, y1means)
  #meanframe <- meanframe * 0.064
  #colnames(meanframe) <- c("x0", "x1", "y0", "y1")
  #return(meanframe)
#}


#MESH$q <- "0"
#MESH$q[MESH$max_um<=max(MR$length[MR$q1=="1"], na.rm=T)] <- "Q1"
#MESH$q[MESH$max_um<=max(MR$length[MR$q1==2], na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==1],na.rm=T)] <- "Q2"
#MESH$q[MESH$max_um<=max(MR$length[(MR$q1)==3],na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==2],na.rm=T)] <- "Q3"
#MESH$q[MESH$max_um<=max(MR$length[(MR$q1)==4],na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==3],na.rm=T)] <- "Q4"

#meanq1 <- superfun(MESH[MESH$q=="Q1",], 30)
#meanq2 <- superfun(MESH[MESH$q=="Q2",], 30)
#meanq3 <- superfun(MESH[MESH$q=="Q3",], 30)
#meanq4 <- superfun(MESH[MESH$q=="Q4",], 30)

#made meanq's by running script above where MESH is replaced by MESH[MESH$q=="Q1",] etc and meanq <- meanframe
#p1 <- p1 + geom_point(data=meanq1, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq1, aes(x=y1,y=x1), colour="white")
#p2 <- p2 + geom_point(data=meanq2, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq2, aes(x=y1,y=x1), colour="white")
#p3 <- p3 + geom_point(data=meanq3, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq3, aes(x=y1,y=x1), colour="white")
#p4 <- p4 + geom_point(data=meanq4, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq4, aes(x=y1,y=x1), colour="white")

##############################save baseplots###########################################################
#save all plots. here make sure you made a seperate folder or change the filenames.
ggsave(p1 + ggtitle("first quartile"), file=paste(name, "first.pdf", sep="_"), width = 10*xmax, height=10*ymax)
ggsave(p2 + ggtitle("second quartile"), file=paste(name, "second.pdf", sep="_"), width = 10*xmax, height=10*ymax)
ggsave(p3 + ggtitle("third quartile"), file=paste(name, "third.pdf", sep="_"),width = 10*xmax, height=10*ymax)
ggsave(p4 + ggtitle("fourth quartile"), file=paste(name, "fourth.pdf", sep="_"),width = 10*xmax, height=10*ymax)
ggsave(pLD, file=paste(name, "Lheat.pdf", sep="_"))
ggsave(pWD, file=paste(name, "Wheat.pdf", sep="_"))
ggsave(pLpoint, file=paste(name, "Lpoint.pdf", sep="_"))
ggsave(pWpoint, file=paste(name, "Wpoint.pdf", sep="_"))

#save the data frame "MR" for further use.
save(MR, file=paste(name, "merged.Rda", sep="_"))

###########################################double histograms!! yay!!!######################################3333
##allplot function combines histograms and density plot.

allplot <- function(plot, data, xmax, ymax, empty){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + theme_bw() + theme(legend.position = "none")
  p1hL <- ggplot(data, aes(x=Lcor)) + geom_histogram() + coord_cartesian(xlim = c(-xmax, xmax)) + theme_bw() +theme(axis.title.x = element_blank())
  p1hD <- ggplot(data, aes(x=Dcor)) + geom_histogram() + coord_flip(xlim = c(-ymax, ymax)) + theme_bw() + theme(axis.title.y = element_blank())

  #align the plots properly before putting them together
  p1Dg <- ggplotGrob(p1D)
  p1hLg <- ggplotGrob(p1hL)
  p1hDg <- ggplotGrob(p1hD)

  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)

  #put the grids together using gridarrange
  return(arrangeGrob(p1hLg, empty, p1Dg, p1hD, ncol=2, nrow=2, widths=c(10*xmax, 2.5), heights=c(2, 10*ymax))) 
}

##before using the function:
#create mockup plot to make space
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



##############################################now plotting#############################################################

#and create the plots:
p1_all <- allplot(p1, Q1, xmax, ymax, empty)
p2_all <- allplot(p2, Q2, xmax, ymax, empty)
p3_all <- allplot(p3, Q3, xmax, ymax, empty)
p4_all <- allplot(p4, Q4, xmax, ymax, empty)

#in case you want it for the whole thing instead of only quartiles:
pall <- ggplot(MR, aes(x=Lcor, y=Dcor))
pall <- densityplot(pall)
pall <- coplot(pall,xmax, ymax, max(MR$length)*0.5)
mppall <- kde2d(MR$Lcor[!is.na(MR$Lcor)&!is.na(MR$Dcor)],MR$Dcor[!is.na(MR$Dcor)&!is.na(MR$Lcor)])
mpp <- mean(range(mppall$z))
pall <- heatmap(pall, mpp)
pall_all <- allplot(pall, MR, xmax, ymax, empty)

#and finally putting the four quartiles below each other: 
ggsave(arrangeGrob(p1_all, p2_all, p3_all, p4_all, ncol=1), filename=paste(name,"allplots_quartiles.pdf", sep="_"), width = 9*xmax, height = 42*ymax)
ggsave(pall_all, filename = paste(name,"allcellsallplot.pdf",sep="_"), width = 11*xmax, height = 11*ymax)

#save all histograms (L coordinates) of the quartiles too:
p1his <- ggplot(Q1, aes(x=Lcor)) + geom_histogram() + theme_bw() + labs(x="Length(?m)")
p2his <- ggplot(Q2, aes(x=Lcor)) + geom_histogram() + theme_bw() + labs(x="Length(?m)")
p4his <- ggplot(Q4, aes(x=Lcor)) + geom_histogram() + theme_bw() + labs(x="Length(?m)")
p3his <- ggplot(Q3, aes(x=Lcor)) + geom_histogram() + theme_bw() + labs(x="Length(?m)")
ggsave(arrangeGrob(p1his, p2his, p3his, p4his, ncol=1), filename=paste(name,"allhis_quartiles.pdf",sep="_"), width=10, height=30)

#and the total amount of cells:
ggsave(ggplot(MR, aes(x=Lcor)) + geom_histogram() + theme_bw() + labs(x="Length(?m)") + ggtitle("all cells"), filename=paste(name,"totalhist.pdf", sep="_"))


###############################################by size instead of quartiles###########################################
#of course you could do other things than evenly sized quartiles. cut by length for instance. if you keep 4 catagories
#you only would have to change the list q into your new cutoff values and run that part again. for instance:
#evenly distributed by length:






