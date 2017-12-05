
rm(list=ls())

setwd("C:/Users/Gebruiker/Documents/ANLIM/Assignments/A2/")


packages <- c("plyr","data.table","mgcv", "evtree", "classInt", "rgdal", "ggplot2", "ggmap","grid","gridExtra","xtable","sas7bdat")
suppressMessages(packages <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}))

################################
# 1) Load and inspect the data #
################################

# Load and rearrange data
source("arrange_data_table.R")

# Structure and summary statistics
dim(DT)
names(DT)
str(DT)
summary(DT)
head(DT)

# Visualization parameters and functions
col <- "#003366"
fill <- "#99CCFF"
ylab <- "Relative exposure"

# This function creates the barplots
ggplot.bar <- function(DT,variable,xlab){
  ggplot(data=DT, aes(as.factor(variable))) + theme_bw() + 
    geom_bar(aes(y = (..count..)/sum(..count..)),col=col,fill=fill) + labs(x=xlab,y=ylab)
}

# This function creates the histograms
ggplot.hist <- function(DT,variable,xlab,binwidth){
  ggplot(data=DT, aes(variable)) + theme_bw() + 
    geom_histogram(aes(y = (..count..)/sum(..count..)),binwidth=binwidth,col=col,fill=fill) + 
    labs(x=xlab,y=ylab)
}
ggplot.hist.group <- function(DT,variable,factor,xlab,binwidth){
  ggplot(data=DT, aes(x=variable, color=factor, fill=factor)) + theme_bw() + 
    geom_histogram(aes(y = (..count..)/sum(..count..)),binwidth=binwidth, alpha=0.5, position="dodge") + 
    labs(x=xlab,y=ylab) + scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") 
}

# Getting insight in the claim frequency and severity components
plot.eda.nclaims <- ggplot.bar(DT,DT$NCLAIMS,"nclaims")
plot.eda.exp <- ggplot.hist(DT,DT$EXP,"exp",0.05)
plot.eda.amount <- ggplot(data=DT, aes(AMOUNT)) + geom_density(adjust=3,col=col,fill=fill) + xlim(0,1e4) + ylab(ylab) + xlab("amount")  + theme_bw()

# Bar plot of factor variables
plot.eda.coverage <- ggplot.bar(DT,DT$COVER,"coverage")
plot.eda.fuel <- ggplot.bar(DT,DT$FUEL,"fuel")
plot.eda.sex <- ggplot.bar(DT,DT$SEX,"sex")
plot.eda.use <- ggplot.bar(DT,DT$USE,"use")
plot.eda.fleet <- ggplot.bar(DT,DT$FLEET,"fleet")
plot.eda.power <- ggplot.bar(DT,DT$POWER,"power")
plot.eda.agec <- ggplot.bar(DT,DT$AGEc,"agec")
plot.eda.sport <- ggplot.bar(DT,DT$SPORT,"sport")
plot.eda.split <- ggplot.bar(DT,DT$SPLIT,"split")

# Histograms of continuous variables
plot.eda.ageph <- ggplot.hist(DT,DT$AGEph,"ageph",2)

# Mutiple layer histograms for interaction effects
plot.eda.agephsex <- ggplot.hist.group(DT,DT$AGEph,DT$SEX,"ageph",2)
plot.eda.agephsport <- ggplot.hist.group(DT,DT$AGEph,DT$SPORT,"ageph",2)
plot.eda.agephsplit <- ggplot.hist.group(DT,DT$AGEph,DT$SPLIT,"ageph",2)
plot.eda.agephcover <- ggplot.hist.group(DT,DT$AGEph,DT$COVER,"ageph",2)
plot.eda.agephpower <- ggplot.hist.group(DT,DT$AGEph,DT$POWER,"ageph",2)


# Putting these together
grid.arrange(plot.eda.nclaims,plot.eda.exp,plot.eda.amount,plot.eda.coverage,plot.eda.fuel,plot.eda.sex,plot.eda.use,plot.eda.fleet,plot.eda.ageph,plot.eda.power,plot.eda.agec,ncol=4)

# Spatial visualizations
DT_PC = DT[,.(sum(EXP)), by=PC][,freq:=V1/sum(V1)][,V1:=NULL]

# This function reads in the shape file of Belgium
readShapefile = function(){
  belgium_shape <- readOGR(dsn = path.expand(paste(dirname(getwd()),"/A2/Shape file Belgie postcodes",sep="")), layer = "npc96_region_Project1")
  belgium_shape <- spTransform(belgium_shape, CRS("+proj=longlat +datum=WGS84"))
  belgium_shape$id <- row.names(belgium_shape)
  return(belgium_shape)
}
belgium_shape = readShapefile()
belgium_shape@data <- merge(belgium_shape@data, DT_PC, by.x = "POSTCODE", by.y = "PC", all.x = TRUE)
belgium_shape@data$freq_class <- cut(as.numeric(belgium_shape@data$freq), breaks = quantile(belgium_shape@data$freq,c(0,0.2,0.8,1),na.rm=TRUE), right = FALSE,include.lowest = TRUE) 
belgium_shape_f <- fortify(belgium_shape)
belgium_shape_f <- merge(belgium_shape_f, belgium_shape@data, all.x = TRUE)

plot.eda.map <- ggplot(belgium_shape_f, aes(long,lat, group = group)) + geom_polygon(aes(fill = belgium_shape_f$freq_class),colour="grey",lwd=0.0000)
plot.eda.map <- plot.eda.map + theme_nothing(legend = TRUE) + labs(fill = "Relative exposure") + scale_fill_brewer(palette="Blues", na.value = "white")

layout <- rbind(c(1,1,2,2,3,3,4,4),
                c(5,5,6,6,7,7,8,8),
                c(9,9,10,10,11,11,NA,NA),c(12,12,12,12,12,12,NA,NA),c(12,12,12,12,12,12,NA,NA))
grid.arrange(plot.eda.nclaims,plot.eda.exp,plot.eda.amount,plot.eda.coverage,plot.eda.fuel,plot.eda.sex,plot.eda.use,plot.eda.fleet,plot.eda.ageph,plot.eda.power,plot.eda.agec,plot.eda.map,layout_matrix=layout)

layout <- rbind(c(1,1,2,2,3,3),c(4,4,5,5,NA,NA))
grid.arrange(plot.eda.agephpower,plot.eda.agephsplit,plot.eda.agephcover,plot.eda.agephsex,plot.eda.agephsport,layout_matrix=layout)

