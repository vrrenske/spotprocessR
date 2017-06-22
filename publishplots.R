#########03-04-12###########################################

#Renske van Raaphorst

#Re-vamping previous analyses to plot for publishing. only length plot, dotplot & overlay of 2 groups/quartile
##NECESARRY PACKAGES####################
library(ggplot2)
library(ggthemes)
library(extrafont)
library(gridExtra)
library(MASS)
library(cowplot)
#library(rJava)
#####LOADING YOUR FILE################################################################################################
OneOrTwo <- readline("Does your file contain information on 1 spot analysis or 2 (result of data merging script)? type 1/2 + ENTER: ")
setwd(choose.dir(default = "F:/microscope files 2015/", caption = "Choose working directory"))

fpfile <- file.choose()
load(fpfile)


##############FUNCTIONS TO PLOT########################################################################################
#general density
densityplot <- function(plot){
  return(plot + stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE)) 
}

#for length/width plot: making y-axis ALWAYS -1.5 to 1.5
LWplot <- function(plot, u="black", MR){
  return(plot + geom_rect(xmin=0, xmax=max(MR$num, na.rm=T), ymin=-1.5, ymax=1.5, fill=u) + theme_minimal()  + scale_x_continuous(limits=c(0,max(MR$num,na.rm=T))) + scale_y_continuous(limits=c(-1.5,1.5)))
}

#heatmapping
heatmap <- function(pdens, mp, colscheme){
  return(pdens + scale_fill_gradient2(low = colscheme[1], mid= colscheme[2], high = colscheme[3], midpoint = mp, space = "Lab") + theme_minimal() + theme(legend.position="none"))
}

#adding the right axis labels + dimensions. Default labels specified. Change by defining labelx and labely manually.
labelmap <- function(pcolor, MR, labelx = "Cell (small to large)", labely = "Location from midcell (\U03BCm)"){
  return(pcolor + coord_fixed(ratio=max(MR$num,na.rm=T)/5) + xlab(labelx) + ylab(labely))
}

#adding the font you want. Default is Arial. You want something else? Specify fontfam manually.
fontchange <- function(pheat, fontfam = "Arial"){
  return(pheat + theme(text = element_text(family=fontfam, size=12)))
}


#function plotting 6 different color schemes. just for fun! so you can choose which one you want.
colfun <- function(plot, mp, MR){
  testplot <- densityplot(plot)
  R <- c("#000000", "#FF0000", "#FFFF00", "#FFFFFF")
  Rp <- heatmap(testplot, mp, R)
  Y <- c("#000000", "#E69F00", "#FFFFFF", "#FFFFFF")
  Yp <- heatmap(testplot, mp,Y)
  O <- c("#000000", "#D55E00", "#F0E442", "#FFFFFF")
  Op <- heatmap(testplot, mp, O)
  G <- c("#000000", "#009E73", "#F0E442", "#FFFFFF")
  Gp <- heatmap(testplot, mp,G)
  B <- c("#000000", "#0072B2", "#FFFFFF", "#FFFFFF")
  Bp <- heatmap(testplot, mp, B)
  W <- c("#FFFFFF", "#D55E00", "#F0E442", "#000000")
  Wp <- heatmap(testplot, mp, W)
  print(plot_grid(Rp,Yp, Op, Gp, Bp, Wp, labels=c("R", "Y", "O", "G", "B", "W")))
  Choice <- readline("Pick color (R/Y/O/G/B/W): ")
  pick <- get(Choice)
  return(heatmap(densityplot(LWplot(plot, u=pick[1], MR)), mp, pick) + geom_line(data=MR, aes(x=num,y=pole1),colour=pick[4]) + geom_line(data=MR, aes(x=num,y=pole2),colour=pick[4]))
}

plotfin <- function(MR){
  p <- ggplot(MR, aes(num, Lmid))
  mpL <- kde2d(MR$num[!is.na(MR$Lmid)], MR$Lmid[!is.na(MR$Lmid)])
  mpL1 <- median(range(mpL$z))  
  p <- colfun(p, mpL1, MR)
  p <- labelmap(p,MR)
  p <- fontchange(p)
  return(p)
}
###################Actual Code#####################################################################  
if(OneOrTwo=="1"){
  p <- plotfin(MR)
  unumin <- max(MR$num[MR$length<(mean(MR$length,na.rm=T)*1.15)], na.rm=T)
  unumax <- min(MR$num[MR$length>(mean(MR$length,na.rm=T)*1.33)], na.rm=T)

####SAVING######
  outfilename <- readline("Name your heatmap: ")
  ggsave(p, file=paste(outfilename, ".pdf", sep=""), width=4, height=3)
  ggsave(p + xlab("") + ylab(""), file=paste(outfilename, "_no_axis_labels.pdf", sep=""), width=8, height=6)
  ggsave(p + xlab("") + ylab("") + geom_rect(xmin=unumin, xmax=unumax, ymin=-1.5, ymax=1.5, alpha=0.003, fill="grey", color=NA), file=paste(outfilename, "_no_axis_labels_forkpoint.pdf", sep=""), width=8, height=6)
}

#########If you have 2 different colors (or WT vs Mutant, etc) ######################################
if(OneOrTwo=="2"){
  ##seperate heatmaps
  unumin <- max(GR$num[GR$length<(mean(GR$length,na.rm=T)*1.115)], na.rm=T)
  unumax <- min(GR$num[GR$length>(mean(GR$length,na.rm=T)*1.33)], na.rm=T)
  c1 <- levels(factor(GR$color))[1]
  c2 <- levels(factor(GR$color))[2]
  p1 <- plotfin(GR[GR$color==c1,])
  print(p1+ ggtitle(c1))
  outfilename1 <- readline("Name heatmap 1: ")
  ggsave(p1, file= paste(outfilename1, ".pdf", sep=""), width=4, height=3)
  ggsave(p1 + xlab("") + ylab(""), file=paste(outfilename1, "_no_axis_labels.pdf", sep=""), width=4, height=3)
  ggsave(p1 + xlab("") + ylab("") + geom_rect(xmin=unumin, xmax=unumax, ymin=-1.5, ymax=1.5, alpha=0.003, fill="grey", color=NA), file=paste(outfilename1, "_no_axis_labels_forkpoint.pdf", sep=""), width=8, height=6)
  p2 <- plotfin(GR[GR$color==c2,])
  print(p2+ ggtitle(c2))
  outfilename2 <- readline("Name heatmap 2: ")
  ggsave(p2, file=paste(outfilename2, ".pdf", sep=""), width=4, height=3)
  ggsave(p2 + xlab("") + ylab(""), file=paste(outfilename2, "_no_axis_labels.pdf", sep=""), width=4, height=3)
  ggsave(p2 + xlab("") + ylab("") + geom_rect(xmin=unumin, xmax=unumax, ymin=-1.5, ymax=1.5, alpha=0.003, fill="grey", color=NA), file=paste(outfilename2, "_no_axis_labels_forkpoint.pdf", sep=""), width=8, height=6)
  bothontop <- plot_grid(p1 + xlab(""), p2, ncol=1)
  save_plot(bothontop, file=paste(outfilename1, outfilename2, "together.pdf", sep="_"), ncol=1, nrow=2, base_height=4, base_width=7)
  
  ##density overlays of quartiles
  #prep
  xmax <- 0.5*max(GR$length, na.rm=TRUE)
  ymax <- 0.5*max(GR$max.width, na.rm=TRUE)
  meansL <- c(mean(GR$length[GR$q1==1], na.rm=TRUE),mean(GR$length[GR$q1==2], na.rm=TRUE), mean(GR$length[GR$q1==3], na.rm=TRUE), mean(GR$length[GR$q1==4], na.rm=TRUE))
  meansW <- c(mean(GR$max.width[GR$q1==1], na.rm=TRUE), mean(GR$max.width[GR$q1==2], na.rm=TRUE), mean(GR$max.width[GR$q1==3], na.rm=TRUE), mean(GR$max.width[GR$q1==4], na.rm=TRUE))
  Q1 <- GR[GR$q1==1,] 
  Q2 <- GR[GR$q1==2,] 
  Q3 <- GR[GR$q1==3,] 
  Q4 <- GR[GR$q1==4,]
  Q2$Lcor <- (Q2$Lmid/Q2$length*meansL[2])
  Q3$Lcor <- (Q3$Lmid/Q3$length*meansL[3])
  Q4$Lcor <- (Q4$Lmid/Q4$length*meansL[4])
  Q2$Dcor <- (Q2$Dum/Q2$max.width*meansW[2])
  Q3$Dcor <- (Q3$Dum/Q3$max.width*meansW[3])
  Q4$Dcor <- (Q4$Dum/Q4$max.width*meansW[4])
  
  #histogram function
  hisfun <- function(data, cscheme, xlabel1="Length axis location from midcell (\U03BCm)"){
    return(ggplot(data, aes(x=Lcor, colour=color), alpha=0.8) + geom_density(aes(fill=color), alpha = 0.5) + coord_cartesian(xlim = c(-xmax, xmax)) + theme_minimal() +theme(legend.position="none", panel.grid.minor=element_blank()) + scale_color_manual(values=cscheme) + scale_fill_manual(values=cscheme)  + xlab(xlabel1))
  }
  
  #make them
  colorpick2 <- function(){
    dat <- data.frame(x=c(1:5), y=1, color = c("R", "Y", "O", "G", "B"))
    R <- "#D55E00"
    Y <- "#F0E442"
    O <- "#E69F00"
    G <- "#009E73"
    B <- "#56B4E9"
    print(ggplot(dat, aes(xmin=x, xmax=x+1, ymin=y-1, ymax=y, fill=color)) + geom_rect() + geom_text(aes(x=x+0.5, y=y-0.5, label=color)) + coord_fixed() + theme(legend.position="none",axis.text = element_blank(), axis.ticks = element_blank()) + xlab("") + ylab("") + scale_fill_manual(values =c(B,G,O, R, Y)))
    col1 <- get(readline("Pick 1th color (R/Y/O/G/B): "))
    col2 <- get(readline("Pick 2nd color (R/Y/O/G/B): "))
    return(cscheme <- c(col1, col2))
  }
  scheme <- colorpick2()
  p1H <- hisfun(Q1, scheme) + xlab("")
  p2H <- hisfun(Q2, scheme) + xlab("")
  p3H <- hisfun(Q3, scheme) + xlab("")
  p4H <- hisfun(Q4, scheme)
  p1H <- fontchange(p1H)
  p2H <- fontchange(p2H)
  p3H <- fontchange(p3H)
  p4H <- fontchange(p4H)
  totalH <- plot_grid(p1H, p2H, p3H, p4H, ncol=1)
  save_plot(totalH, file=paste(outfilename1, outfilename2, "hist.pdf", sep="_"), ncol=1, base_height=12, base_width=5)
  save_plot(plot_grid(p1H + ylab("") + theme(axis.text=element_blank(), axis.ticks=element_blank()), p2H + ylab("") + theme(axis.text=element_blank(), axis.ticks=element_blank()), p3H + ylab("") + theme(axis.text=element_blank(), axis.ticks=element_blank()), p4H+ylab("")+xlab("") + theme(axis.text.y=element_blank(), axis.ticks=element_blank()), ncol=1), file=paste(outfilename1, outfilename2, "hist_no_axis_labels.pdf", sep="_"), base_height=12, base_width=5)

  
}

