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
 ffile$ex <- factor(FN)
 ffile$order <- n
 ffile <- unique(ffile)
 ffile <- ffile[order(ffile$cell),]
 ffile$n <- 1:nrow(ffile)
 plot1 <- plotone(ffile, var)
 #ggsave(plot1, filename = paste(FN, ".pdf", sep=""))
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
 plotall <- function(fincomp, var, B=0.02,palettechoice=c("#999999", "#999999", "#999999", "#999999", "#999999")){
 finfile <- fincomp[c("ex", var)]
 colnames(finfile) <- c("ex", "V")
 return( ggplot(finfile, aes(x=ex, y=V)) + geom_dotplot(aes(fill=ex, color=ex), binaxis="y", stackdir="center", binwidth=B, alpha=0.8)  + theme_minimal() + scale_fill_manual(values=palettechoice) + scale_color_manual(values=palettechoice) + ylab(var) + geom_boxplot(fill=NA,alpha=0.8, size=0.8, outlier.colour="NA", width=0.3))
 }
 
 plotallviolin <- function(fincomp, var){
   finfile <- fincomp[c("ex", var)]
   colnames(finfile) <- c("ex", "V")
   return(ggplot(finfile, aes(x=ex, y=V)) + geom_violin(fill="#999999", color=NA, adjust=0.8) + geom_boxplot(width=0.14, size=1) + theme_minimal())

 }
 
 #color picking function
 colorpick2 <- function(ex)
   {
   dat <- data.frame(x=c(1:5), y=1, color = c("R", "Y", "O", "G", "B"))
   R <- "#D55E00"
   Y <- "#F0E442"
   O <- "#E69F00"
   G <- "#009E73"
   B <- "#56B4E9"
   print(ggplot(dat, aes(xmin=x, xmax=x+1, ymin=y-1, ymax=y, fill=color)) + geom_rect() + geom_text(aes(x=x+0.5, y=y-0.5, label=color)) + coord_fixed() + theme(legend.position="none",axis.text = element_blank(), axis.ticks = element_blank()) + xlab("") + ggtitle(paste("Pick color for ", ex, sep="")) + ylab("") + scale_fill_manual(values =c(B,G,O, R, Y)))
   color <- get(readline(paste("Pick the color for )", ex, " (R/Y/O/G/B): ", sep="")))
   return(color)
   }
 
 
 print(plotall(fincomp, var))
 q <- readline("Happy with the size of the dots? y/n: ")
 newB <- 0.02
 while(q=="n"){ 
 newB <- as.numeric(readline("Give new bin size (default=0.02): "))
 print(plotall(fincomp,var,newB))
 q <- readline("Happy now? y/n: ")
 }
 
 col <- readline("Want to change colors? y/n: ")
 if(col=="y"){
   pal <- c()
   for(n in 1:U){
     val <- colorpick2(levels(fincomp$ex)[n])
     pal[n] <- val
   }
 }
 
 if(q=="y"){
  print(plotall(fincomp,var,newB, pal))
  name <- readline("Name your final output: ")
  if(col=="y"){
  ggsave(plotall(fincomp,var,newB, pal), filename= paste(name, var,".pdf", sep=""))
  }
  ggsave(plotall(fincomp,var,newB), filename= paste(name, var,"_BW.pdf", sep=""))
  ggsave(plotallviolin(fincomp, var), filename=paste(name, var, "_Violin.pdf", sep=""))
  save(fincomp, file="testje.Rda")
  }
 
 