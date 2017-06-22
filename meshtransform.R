##20-05-2015##
##Renske van Raaphorst##

#might be good to clean up first:
rm(list=ls(all=TRUE))

######### Transformation of meshes to get the "mean shape" of the cells. ############################################

#The file used for this transformation is the meshes file you can extract from the microbeTracker GUI.
#This contains the coordinates of all meshes, the maximum length & width of all meshes. The coordinates are
#ordered in a matrix, measuring from one pole to the other in both ways (see microbetracker website). This
#code calculates the mid point of each cell and translates the coordinates in such a way that the mid point is zero. 
#Then the mesh is rotated so the poles are both located on the X axis. 

library(xlsx)

##Since this code uses for-loops, which is not reccomended for R, it's slow. I'm sorry.
##set your working directory
setwd(choose.dir(default = "Y:/staff/fwn/MOLGEN/USERS/Morten_Renske", caption = "Choose working directory"))

#pick a file
Mfile = file.choose()
if(substr(Mfile, start=(nchar(Mfile)-3), stop=(nchar(Mfile))) == ".xls"){
  MESH <-read.table(Mfile,header=T,sep="\t", dec = ".")
  meshname=basename(Mfile)} 
if(substr(Mfile, start=(nchar(Mfile)-3), stop=(nchar(Mfile))) == ".txt"){
  MESH <-read.table(Mfile,header=T,sep="\t", dec = ".")
  meshname=basename(Mfile)} 
#else{load(Mfile)} #here you can load the MESH From ouftitransform

######################### mid-points ################################################################################

meshturn <- function(MESH){
if("frame" %in% colnames(MESH)){ MESH$slice <- MESH$frame}
MESH <- MESH[order(MESH$slice, MESH$cell, MESH$length),]
MESH$Xmid <- 0
MESH$Ymid <- 0
MESH <- MESH[MESH$x0 >1,]

for(n in 1:max(MESH$slice)){
  for(z in as.integer(levels(as.factor(MESH$cell[MESH$slice==n])))){
    X_1 <- MESH$x0[MESH$slice==n&MESH$cell==z][1]
    Y_1 <- MESH$y0[MESH$slice==n&MESH$cell==z][1]
    X_2 <- tail(MESH$x0[MESH$slice==n&MESH$cell==z], n=1)
    Y_2 <- tail(MESH$y0[MESH$slice==n&MESH$cell==z], n=1)
    Xmid <- X_1 + 0.5*(X_2-X_1)
    Ymid <- Y_1 + 0.5*(Y_2-Y_1)
    MESH$Xmid[MESH$slice==n&MESH$cell==z] <- Xmid
    MESH$Ymid[MESH$slice==n&MESH$cell==z] <- Ymid
    print(Xmid)
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

#average of all cells
MESH$x0_rot <- -1*abs(MESH$x0_rot)
MESH$x1_rot <- abs(MESH$x1_rot)

MESH$max_um <- MESH$max_length*0.064

return(MESH)
}

MESH2 <- meshturn(MESH)
save(MESH, file= paste(meshname, ".Rda"))