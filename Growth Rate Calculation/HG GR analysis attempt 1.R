#01/29/19 Script for Allison. Vol-Carbon Conversions for growth rate experiments. Bquick() and makedf() functions pulled from wimFunctions.txt. 
rm(list=ls())

install.packages("segmented")
library(segmented)
#setwd("/Users/allisonadams/OneDrive/Documents/Thesis Work/Growth Rates/Allison")
getwd()

#Load Volume Calibration
load("VolumeCalibration.RData")

HGGR_Test <-read.csv("HGGR_Test.csv")
HGGR_Test

### Example of how to save, if the file is named "x"  Save(x, file="Nameofyourfile.Rdata")
# save(HGGR_Test, file = "HGGR_Test.RData")


#Load 2018 Hunger Games Growth Rate Analysis, YBP2, WLD2, SJR1
load("HGGR_Test.RData")
str(HGGR_Test)

grow <- HGGR_Test
#Create df with date, station, approx.days, rep, label, & volume


x <- grow[,c("Date","Station","exp","days","rep","Volume","X..of.Bugs")]
colnames(x) <- c("date","sta","exp","days","rep","volmm3","CopepodN")
str(x)
x


x$LnV <- log(x$volmm3) #take natural log of volume measurement
x$LnC <- predict(VolumeCalibration$brokenModel,newdata=x ) #load segmented, convert volume to carbon content
x$ugC <- exp(x$LnC) 
str(x)

#load Bquick fxn 
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

xmn <- Bquick(x[,c("sta","exp","date","days","rep","volmm3", "ugC")], 5, median)
xmn$nBugs <- Bquick(x[,c("exp","sta","date","days","rep","CopepodN")], 5, length)$CopepodN
xmn$LnC <- log(xmn$ugC)
str(xmn) 

out <- Bquick(xmn[,c("exp","date","sta")], 1, head)
str(out)

#makedf
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
                                "dfB","intB","slopeB","ciB", "dfL2","intL2","slopeL2","ciL2","p_Bent","Model"), nrow=nrow(out)))
zi <- vector("list",3)
names(zi) <- c("Linear","Bent", "Linear02")
models <- vector("list", 38 )
names(models) <- 1:38

fitAll <- xmn[xmn$days==0 & xmn$rep==1, c("exp","sta","date")]
dayz <- data.frame(days=seq(0, max(xmn$days), by=0.1))
dayz$d1  <- pmin(2, dayz$days)   # Dummy variable for the first time interval
dayz$d2  <- pmax(0, dayz$days-2) # Second time interval
fitAll <- merge(fitAll, data.frame(Exp=rep(1:38,each=nrow(dayz)), days=rep(dayz$days, 38)),1)
fitAll$fit.ci <- fitAll$fit <- rep(NA, nrow(fitAll))


for (i in out$exp)
{
  xi <- xmn[xmn$exp==i,]
  lab <- paste(unlist(xi[1,c("date","sta")]), collapse=" ")
  Day2 <- median(unique(xi$days))
  xi$d1  <- pmin(Day2, xi$days)   # Dummy variable for the first time interval
  xi$d2  <- pmax(0, xi$days-Day2) # Second time interval
  zi[[1]] <- lm (LnC ~ days, data= xi)  # Linear
  zi[[2]] <- lm (LnC ~ d1 + d2, data= xi) # Broken line
  zi[[3]] <- lm (LnC ~ days, data= xi[xi$days <= Day2,])  # Linear
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
    out[i,4:7 + (j-1)*4] <- c(df, cf, ci)
  }
  out$p_Bent[i] <- pValue 
  out$Model[i] <- names(zi)[whichZ]
  cat("------------------\n",i, lab, " Model:",names(zi)[whichZ], 
      "\nLinear Model:\n")
  print(summary(zi[[1]]))
  cat ("-------\nBent Model:\n")
  print(summary(zi[[2]]))
  cat ("-------\nLinear Model 2 Days:\n")
  print(summary(zi[[3]]))
  
  modelFit <- predict(zi[[whichZ]] , newdata=data.frame(dayz), se.fit=T)
  modelFit<- data.frame(modelFit$fit, modelFit$se.fit)
  fitAll[fitAll$Exp==i, 5:6] <- modelFit  
}
out$Model <- factor(out$Model)
ymax <- max(out$slopeL2 + out$ciL2)

#pick best model 
addcol <- ifelse(out$Model=="Bent",8,0)
out <- data.frame(out, makedf(c("df","int","slope","ci"), nrow(out)))
for (i in 1:nrow(out)) out[i, 18:21] <- out[i, (4:7) + addcol[i]]

out$label <- paste(round(out$slope, 2), "+-", round(out$ci,2))
print(out)
View(out)
