###########FTSZ ANGLE MEASUREMENT#########################################################################

#3-2-2016 
#Renske van Raaphorst

##Goal: to obtain the angle of the longest axis of a polygon compared to the polygon it is sitting in. 
##Practical goal: to measure the angle of FtsZ perpundicular to the long cell axis to see if the localization
##is skewed.

#usual cleanup
rm(list=ls(all=TRUE))

#librarys
library(ggplot2)
library(shotGroups)
library(ggthemes)

#take the MESH file from the ouftimeshtransform. Also the object file can come from there
#MESH processed a little bit
simplemesh <- function(MESH){
  MESH2 <- MESH[c("frame","cell", "x0","y0", "length")]
  MESH3 <- MESH[c("frame","cell", "x1","y1", "length", "max_length")]
  colnames(MESH3) <- c("frame", "cell", "x0", "y0", "length", "max_length")
  MESH3$n <- MESH3$max_length + MESH3$max_length-MESH3$length
  MESH2$n <- MESH2$length
  MESH4 <- merge(MESH2, MESH3, all=T)
  MESH4 <- MESH4[order(MESH4$frame,MESH4$cell, MESH4$n),]
  return(MESH4)
}


anglefun <- function(dat, xie, yie){
  dat$angle <- NA
  dat$centre_x <- NA
  dat$centre_y <- NA
  for(n in unique(dat$frame)){
    for(u in unique(dat$cell)){
      if(nrow(dat[dat$cell==u&dat$frame==n,])>0){
        woei <- as.matrix(dat[c(xie, yie)][dat$cell==u&dat$frame==n,])
        an <- getMinBBox(woei)$angle
        dat$angle[dat$cell==u&dat$frame==n] <- an
        cen <- getMinCircle(woei)$ctr
        dat$centre_x[dat$cell==u&dat$frame==n] <- cen[1]
        dat$centre_y[dat$cell==u&dat$frame==n] <- cen[2]
        }
    }
  }
  return(dat)
}

compareangle <- function(M, O){
  M <- unique(M[c("cell", "frame", "angle")])
  O <- unique(O[c("cell", "frame", "angle")])
  colnames(O)[3] <- "obj_an"
  MO <- merge(M, O, all=T)
  MO$dif <- MO$angle - MO$obj_an
  return(MO)
}

#########CODE###########################################################################
finalfun <- function (MESH, OBJ){
  MESH <- simplemesh(MESH)
  MESH <- MESH[!is.na(MESH$x0),]
  MESH <- anglefun(MESH, "x0", "y0")
  OBJ <- OBJ[!is.na(OBJ$ob_x),]
  OBJ <- anglefun(OBJ, "ob_x", "ob_y")
  fincomp <- compareangle(MESH, OBJ)
  fincomp$dif <- abs(fincomp$dif)
  fincomp <- fincomp[fincomp$dif<165&fincomp$dif>15,]
  fincomp$dif90 <- abs(90-fincomp$dif)
  return(fincomp)
}

setwd(choose.dir(caption = "Choose working directory"))
U <- readline("put output file name: ")
load(file.choose())
load(file.choose())
fincomp <- finalfun(MESH, OBJ)
save(fincomp, file = paste(U, ".Rda", sep=""))

