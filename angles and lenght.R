###########################plotting and comparing multiple angle analyses#################################
##26-02-2016
##Renske van Raaphorst

library(ggplot2)
library(ggthemes)

##If you source, you can compare the variable you indicate in a plot quickly 

########################################################## 
plotone <- function(fincomp, var){
  plotlist <- fincomp[var]
  colnames(plotlist) <- "z"
  ggplot(plotlist, aes(x=z)) + geom_histogram(aes(y=..density..), fill="#0072B2", alpha=0.7) + theme_minimal() + annotate(geom="text", label=paste("Median =", round(median(plotlist$z,na.rm=T), digits=2), "degrees\nVariance =", round(sqrt(sd(plotlist$z,na.rm=T)),digits=2), sep=" "), x=0, y=0.02)
}

#####load object in seperate environment so you can name it yourself#########################
load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

#####################################################################
U <- readline("How many data sets do you want to compare? ")
var <- readline("What is the name of the variable you want to compare? ")
for(n in 1:U){
  firstfile <- file.choose()
  ffile <- load_obj(firstfile)
  FN <- readline("Name this dataset: ")
  ffile$ex <- FN
  ffile <- unique(ffile)
  ffile <- ffile[order(ffile$cell),]
  ffile$n <- 1:nrow(ffile)
  plot1 <- plotone(ffile, var)
  ggsave(plot1, filename = paste(FN, ".pdf", sep=""))
  if(n == 1){
    fincomp <- ffile
    list1 <- c(nrow(ffile))
  }
  if(n > 1){
    fincomp <- merge(fincomp, ffile, all=T)
    list1 <- append(list1, nrow(ffile))
  }
  
}

fincomp <- fincomp[fincomp$n<=min(list1),]

##########COMPARE DISTRIBUTION OF ANGLES
#Load different outcomes, name the type with ex
plotall <- function(fincomp, var, B=0.02){
  finfile <- fincomp[c("ex", var)]
  colnames(finfile) <- c("ex", "V")
  return( ggplot(finfile, aes(x=ex, y=V)) + geom_dotplot(aes(fill=ex, color=ex), binaxis="y", stackdir="center", binwidth=B, alpha=0.5)  + theme_minimal() + scale_colour_tableau() + scale_fill_tableau() + ylab(var) + geom_boxplot(alpha=0.5, outlier.colour="NA", fill="NA"))
}

print(plotall(fincomp, var))
q <- readline("Happy with the size of the dots? y/n: ")
newB <- 0.02
while(q=="n"){ 
  newB <- as.numeric(readline("Give new bin size (default=0.02): "))
  print(plotall(fincomp,var,newB))
  q <- readline("Happy now? y/n: ")
}

if(q=="y"){
  ggsave(plotall(fincomp,var,newB), filename= paste("Boxall_", var,".pdf", sep=""))
}