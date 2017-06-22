#30-05-2016
#Renske van Raaphorst

#might be good to clean up first:
rm(list=ls(all=TRUE))


##packages needed
library(ggplot2)
library(shotGroups)
library(ggthemes)

##########################Little extra code to put the middle of the cell into the MESH 
##functions

#get centre of object
centrefun <- function(dat, xie, yie){
  dat$centre_x <- NA
  dat$centre_y <- NA
  for(n in unique(dat$frame)){
    for(u in unique(dat$cell)){
      if(nrow(dat[dat$cell==u&dat$frame==n,])>0){
        for(i in unique(dat$obnum[dat$cell==u&dat$frame==n])){
          woei <- as.matrix(dat[c(xie, yie)][dat$cell==u&dat$frame==n&dat$obnum==i,])
          cen <- getMinCircle(woei)$ctr
          dat$centre_x[dat$cell==u&dat$frame==n&dat$obnum==i] <- cen[1]
          dat$centre_y[dat$cell==u&dat$frame==n&dat$obnum==i] <- cen[2]
        }
      }
    }
  }
  return(dat)
}

##turning the mesh
meshturn <- function(MESH){
  MESH <- MESH[!is.na(MESH$x0),]
  if("frame" %in% colnames(MESH)){ MESH$slice <- MESH$frame}
  MESH <- MESH[order(MESH$slice, MESH$cell, MESH$length),]
  MESH$Xmid <- 0
  MESH$Ymid <- 0
  MESH <- MESH[MESH$x0 >1,]
  
  for(n in 1:max(MESH$slice)){
    for(z in as.integer(levels(as.factor(MESH$cell[MESH$slice==n])))){
      X_1 <- MESH$x0[MESH$slice==n&MESH$cell==z][1]
      Y_1 <- MESH$y0[MESH$slice==n&MESH$cell==z][1]
      X <- MESH$x0[MESH$slice==n&MESH$cell==z]
      X_2 <- tail(X[!is.na(X)], n=1)
      Y <- MESH$y0[MESH$slice==n&MESH$cell==z]
      Y_2 <- tail(Y[!is.na(Y)], n=1)
      Xmid <- X_1 + 0.5*(X_2-X_1)
      Ymid <- Y_1 + 0.5*(Y_2-Y_1)
      MESH$Xmid[MESH$slice==n&MESH$cell==z] <- Xmid
      MESH$Ymid[MESH$slice==n&MESH$cell==z] <- Ymid
      
    }
  }
  
  ######################### mid-point correction ######################################################################
  
  
  MESH$x0_cor <- MESH$x0 - MESH$Xmid
  MESH$y0_cor <- MESH$y0 - MESH$Ymid
  MESH$x1_cor <- MESH$x1 - MESH$Xmid
  MESH$y1_cor <- MESH$y1 - MESH$Ymid
  
  ########################## rotation #################################################################################
  
  #adding a column with the rotation angle for each cell
  MESH$angle <- 0
  for(n in 1:max(MESH$slice)){
    for(z in as.integer(levels(as.factor(MESH$cell[MESH$slice==n])))){
      Xmin_c <- MESH$x0_cor[MESH$slice==n&MESH$cell==z][1]
      Ymin_c <- MESH$y0_cor[MESH$slice==n&MESH$cell==z][1]
      MESH$angle[MESH$cell==z&MESH$slice==n] <- atan(Xmin_c/Ymin_c)
    }
  }
  
  #rotating the coordinates
  MESH$x0_rot <- MESH$x0_cor * cos(MESH$angle) - MESH$y0_cor * sin(MESH$angle)
  MESH$x1_rot <- MESH$x1_cor * cos(MESH$angle) - MESH$y1_cor * sin(MESH$angle)
  MESH$y0_rot <- MESH$x0_cor * sin(MESH$angle) + MESH$y0_cor * cos(MESH$angle)
  MESH$y1_rot <- MESH$x1_cor * sin(MESH$angle) + MESH$y1_cor * cos(MESH$angle)
  MESH <- MESH[!is.na(MESH$y1_rot),]
  #average of all cells
  MESH$x0_rot <- -1*abs(MESH$x0_rot)
  MESH$x1_rot <- abs(MESH$x1_rot)
  
  MESH$max_um <- MESH$max_length*magnification
  
  return(MESH)
}

#add object centre to mesh file and turn accordingly
midobject <- function(MESH, OBJ){
  MESH <- merge(MESH, OBJ, all=T)
  MESH$xccor <- MESH$centre_x - MESH$Xmid
  MESH$yccor <- MESH$centre_y - MESH$Ymid
  MESH$Dum <- MESH$xccor * cos(MESH$angle) - MESH$yccor * sin(MESH$angle)
  MESH$Lmid <- MESH$xccor * sin(MESH$angle) + MESH$yccor * cos(MESH$angle)
  MO <- MESH[,c("frame", "cell", "obnum", "max_um", "max.width", "Dum", "Lmid")]
  MO$max.width <- MO$max.width*0.064
  MO$Dum <- MO$Dum * 0.064
  MO$Lmid <- MO$Lmid * 0.064
  MO <- unique(MO)
  MO <- MO[order(MO$max_um),]
  MO$num <- 1:nrow(MO)
  MO$length <- MO$max_um
  MO$pole1 <- 1/2*MO$length
  MO$pole2 <- -1*MO$pole1
  return(MO)
}

##code################################################################
##set your working directory
setwd(choose.dir(default = "Y:/staff/fwn/MOLGEN/USERS/Morten_Renske", caption = "Choose working directory"))
fpfile <- file.choose()
load(fpfile)
Mfile <- file.choose()
load(Mfile)

magnification <- as.numeric(readline("Magnification: (. delimeter)"))

OBJ <- OBJ[!is.na(OBJ$ob_x),]
OBJs <- centrefun(OBJ, "ob_x", "ob_y")
OBJs <- OBJs[,c("frame", "cell", "obnum", "centre_x", "centre_y")]
OBJs <- unique(OBJs)
MESH <- meshturn(MESH)
MR <- midobject(MESH, OBJs)

print(ggplot(MR, aes(x=num, y=Lmid)) + geom_point())
filename <- readline("Name file: ")
save(MR, file=paste("Objectmiddle_", filename, ".Rda", sep=""))

#now you can use publishplots or other things to plot it nicely
