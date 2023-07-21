###Code Newly Rewritten by Wim, received 3/25/21

## Allison's trial with her data
## After some errors and notes from Wim, Success!! on 5/17/21
## Worked on it again 5/27/21

#=============================================================
#  This code brings in raw data from image spreadsheets and 
#  calculates carbon per bug, then calculates slopes over time
#  and picks the best slope, i.e., if the slope does not change over 
#  the experiment, that is the growth rate, but if the slope changes
#  the growth rate is the slope between the first 2 time points.
#  NOTE that this will have to be rejiggered for more than 2 time points.
#=============================================================
rm(list=ls())

library(readxl)
library(segmented)
install.packages("ggplot2")
library("ggplot2")
install.packages("writexl")
library(writexl)
## Also need to load the package for reading xlsx

#------------------------
# Replace Rdata file, Directory, and fn with new information as needed
# Then replace sample names and spreadsheet tabs
# The sample names are in the file


#Load Volume Calibration
load("/Users/allisonadams/Documents/Thesis/Growth Rates/Master Data/VolumeCalibration 2021.RData")
dir <- "/Users/allisonadams/Documents/Thesis/Growth Rates/Master Data/"

##Don't need to load this since it is already loaded above and it worked.
#load("VolumeCalibration 2021.RData")

## I got an error code 'path' does not exist. So trying setwd instead
## But later it worked, so I don't use this anymore
setwd("/Users/allisonadams/Documents/Thesis/Growth Rates/Master Data")
getwd()

fn <- "HG 2019 Growth Rate Image Analysis All Stations.xlsx"

samples <- paste0 (rep(c("LSZ","SJR","WLD","YBP"), each=2), rep(1:2, 4) )
tabs    <- paste0 (samples, " All days" )
## ?? The above are station names, correct? And does the excel data sheet consist of one spreadsheet
### tab for each station? And this command will tell R to use each sheet to calculate the growth rates
### of each station separately?


x       <- NULL

for (i in 1:8)
{
  xi <- read_xlsx(paste0(dir,fn))
  xi <- data.frame(xi)
  xi <- xi[,c("Station","days","rep","Volume.mm3", "X..of.Bugs")]
  names(xi) <- c("sample","days","rep","volmm3", "nBugs")
  x  <- rbind(x,xi)
}
summary(x)   # Check no missing values

#-------------------
# Calculate log of volume, then use the calibration to get log carbon

x$LnV <- log(x$volmm3)
x$LnC <- predict(VolumeCalibration$brokenModel,newdata=x ) # convert volume to carbon content
x$ugC <- exp(x$LnC) 

#-------------------
# Add Bquick function since it is used below
Bquick <- function(x, ind, FUN,  ...)
{
  # Uses aggregate but converts index columns back to same mode as source
  # Revised 6/15/02 to fix output from aggregate, which converts indices to factors or character
  # x is a data frame
  # ind is the number of columns (starting from the left) over which aggregation occurs
  # In either case factors in the input are returned as factors
  #   and ordered factors are returned as ordered factors
  # The ellipses ... can mean any argument passed to fun
  if(ind > 1.) w <- as.list(x[, 1.:ind])
  else 
  {
    w <- vector("list", length=1)
    w[[1]] <- x[, 1.]
  }
  xd <- data.frame(x[, (ind + 1.):dim(x)[2.]], stringsAsFactors = F)
  y <- aggregate(xd, w, FUN, ...)
  for(i in 1.:ind) {
    if(is.ordered(x[, i]))
      y[, i] <- ordered(y[, i], levels = levels(x[, i]))
    else if(is.factor(x[, i]))
      y[, i] <- factor(y[, i], levels = levels(x[, i]))
    else {
      y[, i] <- as.character(y[, i])  	
      if(is.numeric(x[,i])) y[, i] <- as.numeric(y[,i])
    }}
  
  for(i in seq(ind, 1., -1.))
    y <- y[order(y[, i]),  ]
  names(y) <- names(x)
  row.names(y) <- 1.:dim(y)[1.]
  y
}

# Calculate medians of volume and carbon

xmn <- Bquick(x[,c("sample","days","rep","volmm3", "LnC","ugC")], 3, median)

# Add rowindex function since it is used below
rowindex <- function(x, sort = F, global = T, countGroup = !global )
{
  # Returns an index of rows for each unique combination of columns in input.
  # Modified 3/16/03 to change the algorithm completely, 
  # and to allow for the data  not to be sorted first, 
  # MODIFIED 7/19/03 to make it more efficient, and to add as.numeric to get correct sorting
  # in which case the index will change each time the inputs change whether sorted or not.
  # MODIFIED 9/1/03 to convert character data to factor, then numeric
  # MODIFIED 7/15/08 to fix a problem with the indexing.  Should work OK now
  # Modified 8/6/09 to trap NA's in the indices
  # If sort is T the data frame is sorted from left to right.  The indices are then
  #  returned in the correct order for x (this is not too meaningful for global=F)
  # If global is T the indices are the same for each unique combination of columns, and run from
  #  1 to length(unique combination).
  # If global is F the indices run within each combination from 1 to the length of the combination
  #  		(this is also specified as countGroup=T, which overrides global)
  
  if (!is.data.frame(x)) x <- as.data.frame(x)
  row.names(x) <- seq(x[, 1])
  nx <- dim(x)[1]
  dx <- dim(x)[2]
  if (sort)
    x <- wsort(x, 1:dx)
  for(i in 1:dx)
    if(!is.numeric(x[, i])) x[, i] <- as.integer(as.factor(x[, i]))	
  if (nrow(na.omit(x)) != nx) stop ("Idiot!  No missing values in indices")
  xm <- as.matrix(x)
  if (dx == 1)
    v <- x[, 1]
  else
    v <- cumsum(c(T, apply(diff(xm) !=0,1,any)))
  names(v) <- row.names(x)
  rL <- rle(v)$lengths
  if (!countGroup)
    vv <- rep(seq(length(rL)), rL)
  else vv <- unlist(lapply(rL, seq))
  if(sort)
    vv <- vv[order(as.numeric(names(v)))]
  vv
}

xmn$seq <- rowindex(xmn$sample)
xmn$nBugs <- Bquick(x[,c("sample","days","rep","nBugs")], 3, length)$CopepodN

#-------------------
# Set up output files for analysis of growth rate

out <- unique(xmn[,c("seq","sample")])

# Add makedf funtion sincd it is used below

makedf <- function (names, first.col=NULL, nrow=1, fill=NA )
{
  # Function to set up a data frame to be filled.
  # Input:
  #	names		Names of the columns
  #	first.col	If supplied, the values for the first column of the data frame
  #	nrow		Number of rows (not needed if first.col is supplied)
  #	fill		Value to fill the df with
  
  if (!is.null (first.col))
    nrow <- length(first.col)
  x <- data.frame(matrix(fill, nrow=nrow, ncol=length(names)))
  names(x) <- names
  if (!is.null (first.col)) 		x[,1] <- first.col
  x
}

out <- data.frame(out, makedf(c("dfL","intL","slopeL","ciL",
                                "dfB","intB","slopeB","ciB", 
                                "dfL2","intL2","slopeL2","ciL2","p_Bent","Model"), 
                                 nrow=nrow(out)))
zi <- vector("list",3)
names(zi) <- c("Linear","Bent", "Linear02")
models <- vector("list", nrow(out) )
names(models) <- out$sample

# dayz <- data.frame(days=seq(0, max(xmn$days), by=0.1))
# dayz$d1  <- pmin(2, dayz$days)   # Dummy variable for the first time interval
# dayz$d2  <- pmax(0, dayz$days-2) # Second time interval

fitAll <-Bquick(xmn[,c("sample","days","rep")], 2, first)[,1:2]
## The above code gave me an error because "first" is a function Wim made.
## He said it could be replaced with the following code, see my notes doc for 
## detailed explanation

fitAll <-Bquick(xmn[,c("sample","days","rep")], 2, head, n=1)[,1:2]

fitAll$fit.ci <- fitAll$fit <- fitAll$d2 <- fitAll$d1 <- rep(NA, nrow(fitAll))

#-------------------
# Loop through the samples and calculate parameters for slope in each.
# Three models are fitted: 
#    Model 1 is linear in all data
#    Model 2 is a broken line
#    Model 3 is linear only over the first time interval
#    The code picks the best model based on the p value for the bent model
#       (p value < 0.05 implies the bent model is a better fit)

for (i in 1:nrow (out))
{
  samplei <- samples[i]
  xi <- xmn[xmn$seq==i,]
  Day2 <- median(unique(xi$days))
  xi$d1  <- pmin(Day2, xi$days)   # Dummy variable for the first time interval
  xi$d2  <- pmax(0, xi$days-Day2) # Second time interval
  zi[[1]] <- lm (LnC ~ days, data= xi) # Linear across all data
  zi[[2]] <- lm (LnC ~ d1 + d2, data= xi) # Broken line
  zi[[3]] <- lm (LnC ~ days, data= xi[xi$days <= Day2,])  # Linear in first segment
  models[[i]] <- zi
  pValue <- anova(zi[[1]],zi[[2]])$"Pr(>F)"[2]  # p value for Bent model 
  whichZ <- 1 + (pValue  < 0.05)# Pick "best"
  for (j in 1:3)
  {
    z <- zi[[j]]   
    df <- z$df.residual
    cf <- coef(z)[1:2]  # Coefficients 1 and 2 only (in Bent model, first time interval)
    se <- sqrt(diag(vcov(z)))[2]  # Standard errors of coefficient
    ci <- se * -qt(.025, df)  
    out[i,3:6  + (j-1)*4] <- c(df, cf, ci)
  }
  out$p_Bent[i] <- pValue 
  out$Model[i] <- names(zi)[whichZ]
  cat("------------------\n",i, xi$sample[1], " Model:",names(zi)[whichZ], 
      "\nLinear Model:\n")
  print(summary(zi[[1]]))
  cat ("-------\nBent Model:\n")
  print(summary(zi[[2]]))
  cat ("-------\nLinear Model 2 Days:\n")
  print(summary(zi[[3]]))
  
  
  fitRows <- which(fitAll$sample==samples[i])
  if (whichZ==2)
    {
      fitAll$d1[fitRows]  <- pmin(Day2, fitAll$days[fitRows]) 
      fitAll$d2[fitRows]  <- pmax(0, fitAll$days[fitRows]-Day2) 
    } 
  modelFit <- predict(zi[[whichZ]] , newdata=fitAll[fitRows,],
                      se.fit=T)
  modelFit <- data.frame(modelFit$fit, modelFit$se.fit)
   
  fitAll[fitRows, c("fit","fit.ci")] <- modelFit  
} # End Out loop

out$Model <- factor(out$Model)

#----------------------------
# Pick best model and put it in the last 4 columns.
 
addcol <- ifelse(out$Model=="Bent",8,0)
outx <- makedf(c("df","int","slope","ci"), nrow(out))
for (i in 1:nrow(out)) outx[i, ] <- out[i, (3:6) + addcol[i]]
outAll <- data.frame(out, outx)
outAll$label <- paste(round(outAll$slope, 2), "+-", round(outAll$ci,2))

#---------------------------
# Save the data
# summaryStats is for reporting growth rates
# fitAll gives the values in log-Carbon (ug) for the lines on plots of median data
# models contains all three models for each unique sample.s

growthCalculations <- list(
      rawData       = x, 
      sampleMedians = xmn,
      summaryStats  = outAll,
      fitAll        = fitAll,
      models        = models)
save(growthCalculations, 
     file=paste0(dir, "growthCalculations_HungerGames.Rdata"))

#### Make sure you've installed ggplot2, wimPalettes and wimGraph before running this code.

wimPalettes <- function (pal="point", ncolors=12, 
                         red=seq(0,1, length=ncolors), 
                         green=0, blue=seq(1,0,length=ncolors))
{
  # Selections are "gradient","point","line", "slide" as of 10/17/2016
  # Added "month" for monthly colors from linear interpolation 
  #     of blue,green,red,orange, and back to blue (January April July)
  # Gradient uses rgb for a selected number of colors ncolors.
  # Gradient also uses the red, green, and blue values above
  # Default gradient is blue to red through purple. 
  require(grDevices)
  nc <- ceiling(ncolors/7)
  palG <- rgb(red=red, green=green, blue=blue)
  palP <- rep(c("black","red","green4","blue","magenta","cyan","orange","gray"),nc)
  palL <- rep(c("red","black","green4","magenta","blue","gray","cyan","orange"),nc)
  palM <- c("#0000ff", "#0055aa", "#00aa55", "#00ff00", "#55aa00", "#aa5500", 
            "#ff0000", "#ff4400", "#ff8800", "#ffcc00", "#aa8855", "#5544aa")  # Blue..Green..Red..Orange..(blue)
  palO <- c("seagreen4", "red","darkgoldenrod4", "royalblue4","coral","gray50", "lightblue1","turquoise4","darkblue")
  palS <- rep(c("yellow","green","cyan","magenta","gray90","pink","orange"),nc)
  switch(tolower(substring(pal,1,1)), "g"=palG, "p"=palP, "l"=palL,"m"=palM, "o"=palO, "s"=palS, palP)
}


wimGraph <- function(background="white", panelBackground=NA,
                     panelBorder = "black",
                     textCol="black", titleCol="black", 
                     axisCol="black",gridCol="gray90",
                     stripCol="black", stripFill="white",
                     title.rel=1.2, axis.title.rel=1.6, axis.text.rel=1.0,
                     strip.text.rel=0.7, textSize=10)
{ 
  theme(panel.border = element_rect(colour = panelBorder, fill=NA, size=0.5),
        plot.background = element_rect(fill = background, linetype=0),
        panel.background = element_rect(fill = panelBackground, linetype=1, color=axisCol),
        legend.key = element_rect(fill = NA),
        legend.background=element_rect(fill=NA),
        text= element_text(color=textCol, size=textSize),
        panel.grid.major = element_line(color=gridCol, size=0.2),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(color=stripCol, fill=stripFill),
        strip.text = element_text(color=textCol, size=rel(strip.text.rel)),
        plot.title = element_text(size=rel(title.rel), color=titleCol),
        axis.title  = element_text(size=rel(axis.title.rel), color=titleCol, vjust=0),
        axis.text = element_text(size=rel(axis.text.rel),color=axisCol),
        axis.line.x =  element_line(color=axisCol),			# Changed these two
        axis.line.y =  element_line(color=axisCol),
        axis.ticks=element_line(color=axisCol))
}


ggplot(growthCalculations$sampleMedians, aes((days), ugC)) +
  geom_point(color="blue", alpha=0.5, size=1) +
  geom_smooth(method = "lm", se = TRUE, color="red", lwd=0.2) +
  scale_y_log10(limits=c(0.05, 1), breaks=c(0.1, 0.2,.5,1)) +
  scale_color_manual(values=wimPalettes()) +
  geom_text(data=growthCalculations$sampleMedians, 
            aes(x=0, y=1, label=sample, group=NULL), size=2, col=2, adj=0) + 
  xlab("Time, d") +
  #facet_wrap(ncol = 2) +
  facet_wrap(~ sample, ncol= 2) +
  ylab(expression(paste("Carbon, ", mu,"g"))) +	
  wimGraph(textSize = 10)



install.packages("writexl")
library(writexl)

write_xlsx(outAll, "/Users/allisonadams/OneDrive/Documents/Thesis Work/Growth Rates/Master Data/outAll")


## Now try to plot the growth rates. Look at the cache slough examples Wim gave me