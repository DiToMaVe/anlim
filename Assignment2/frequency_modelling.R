
###################################################
# 2) Fit a flexible GAM for frequency to the data #
###################################################

# Function to plot the GAM effects
ggplot.gam <- function(model,variable,gam_term,xlabel,ylabel){
  pred <- predict(model, type = "terms", se = TRUE)
  col_index <- which(colnames(pred$fit)==gam_term)
  x <- variable
  b <- pred$fit[, col_index]
  l <- pred$fit[, col_index] - qnorm(0.975) * pred$se.fit[, col_index]
  u <- pred$fit[, col_index] + qnorm(0.975) * pred$se.fit[, col_index]
  df <- unique(data.frame(x, b, l, u))
  p <- ggplot(df, aes(x = x))
  p <- p + geom_line(aes(y = b), size = 1,col="#003366")
  p <- p + geom_line(aes(y = l), size = 0.5, linetype = 2,col="#99CCFF")
  p <- p + geom_line(aes(y = u), size = 0.5, linetype = 2,col="#99CCFF")
  p <- p + xlab(xlabel) + ylab(ylabel) + theme_bw()
  p
}

# Choose reference levels for the GAM (most exposure)
DT$COVER = relevel(DT$COVER,ref="MTPL") 
DT$FUEL = relevel(DT$FUEL,ref="Petrol")
DT$POWER = relevel(DT$POWER,ref="<66")
DT$SPLIT = relevel(DT$SPLIT,ref="Once")
# 
# # Base model 
# gam.freq1 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + s(AGEph)+ s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq1)
# 
# # Base model + SPORT
# gam.freq2 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + SPORT + s(AGEph)+ s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq2)
# 
# # Base model + SPLIT
# gam.freq3 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + SPLIT + s(AGEph)+ s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq3)
# 
# # Base model + SPORT + SPLIT
# gam.freq4 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + SPORT + SPLIT + s(AGEph)+ s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq4)
# 
# # Base model + AGEph_POWER
# gam.freq5 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + s(AGEph) + s(AGEph,by=POWER, bs="fs") + s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq5)
# 
# # Base model + SPLIT + AGEph_POWER
# gam.freq6 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + SPLIT + s(AGEph) + s(AGEph,by=POWER, bs="fs") + s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq6)
# 
# # Base model + SPLIT + AGEph_SPLIT
# gam.freq7 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + SPLIT + s(AGEph) + s(AGEph,by=SPLIT, bs="fs") + s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq7)
# 
# # Base model + SPLIT + AGEph_SPORT
# gam.freq8 <- gam(NCLAIMS  ~  COVER + FUEL + POWER + SPLIT + s(AGEph) + s(AGEph,by=SPORT, bs="fs") + s(LONG,LAT), offset = log(EXP), data = DT, family = poisson(link = "log"))
# summary(gam.freq8)
# 
# save(gam.freq1, file = "base.rda")
# save(gam.freq2, file = "base_sport.rda")
# save(gam.freq3, file = "base_split.rda")
# save(gam.freq4, file = "base_split_sport.rda")
# save(gam.freq5, file = "base_agepower.rda")
# save(gam.freq6, file = "base_split_agepower.rda")
# save(gam.freq7, file = "base_split_agesplit.rda")
# save(gam.freq8, file = "base_split_agesport.rda")
# 
# BIC(gam.freq1)
# BIC(gam.freq2)
# BIC(gam.freq3)
# BIC(gam.freq4)
# BIC(gam.freq5)
# BIC(gam.freq6)
# BIC(gam.freq7)
# BIC(gam.freq8)

load("base_split.rda")

gam.freq <- gam.freq3

plot.gam.freq.ageph <- ggplot.gam(gam.freq,DT[, AGEph],"s(AGEph)","ageph",expression(hat(f)[1](AGEph)))

belgium_shape = readShapefile()
DT_maps <- data.table(PC=belgium_shape@data$POSTCODE,LONG=coordinates(belgium_shape)[,1],LAT=coordinates(belgium_shape)[,2])
DT_maps[,c("COVER","FUEL","AGEph","POWER","SPLIT","EXP"):=DT[1,.(COVER,FUEL,AGEph,POWER,SPLIT,EXP)]]
pred = predict(gam.freq,newdata = DT_maps,type = "terms",terms = "s(LONG,LAT)")
DT_pred = data.table("PC"=DT_maps$PC,"LONG"=DT_maps$LONG,"LAT"=DT_maps$LAT,pred)
belgium_shape@data <- merge(belgium_shape@data,DT_pred,by.x="POSTCODE",by.y="PC",all.x=TRUE)
belgium_shape_f <- fortify(belgium_shape)
belgium_shape_f <- merge(belgium_shape_f, belgium_shape@data,by="id",all.x=TRUE)

plot.gam.freq.map <- ggplot(belgium_shape_f, aes(long, lat, group = group)) + geom_polygon(aes(fill = belgium_shape_f$`s(LONG,LAT)`))
plot.gam.freq.map <- plot.gam.freq.map + theme_nothing(legend = TRUE) + scale_fill_gradient(low="#99CCFF",high="#003366") + labs(fill = expression(hat(f)[2](long,lat)))

layout <- rbind(c(1,1),
                c(2,2))
grid.arrange(plot.gam.freq.ageph,plot.gam.freq.map,layout_matrix=layout)

#############################
# 3) Bin the spatial effect #
#############################

belgium_shape = readShapefile()
DT_maps <- data.table(PC=belgium_shape@data$POSTCODE,LONG=coordinates(belgium_shape)[,1],LAT=coordinates(belgium_shape)[,2])
DT_maps[,c("COVER","FUEL","AGEph","POWER","SPLIT","EXP"):=DT[1,.(COVER,FUEL,AGEph,POWER,SPLIT,EXP)]]
pred = predict(gam.freq,newdata = DT_maps,type = "terms",terms = "s(LONG,LAT)")
GAM.freq.LONGLAT = data.table("pc"=factor(DT_maps$PC),"long"=DT_maps$LONG,"lat"=DT_maps$LAT,pred)
names(GAM.freq.LONGLAT) = c("pc","long","lat","s")
GAM.freq.LONGLAT <- GAM.freq.LONGLAT[order(pc)]

num_bins = 8

classint.fisher = classIntervals(GAM.freq.LONGLAT$s,num_bins,style="fisher")

# Plot empirical cdf with class intervals
crp=colorRampPalette(c("#99CCFF","#003366"))  
plot_ecdf <- plot(classint.fisher,crp(num_bins),xlab=expression(hat(f)[2](long,lat)),main="Fisher")

# Plot maps
belgium_shape <- readShapefile()
belgium_shape@data <- merge(belgium_shape@data,GAM.freq.LONGLAT[,.(pc,s)],by.x="POSTCODE",by.y="pc",all.x=TRUE)
belgium_shape@data$class_fisher <- cut(as.numeric(belgium_shape@data$s), breaks = classint.fisher$brks, right = FALSE,include.lowest=TRUE,dig.lab = 2) 
belgium_shape_f <- fortify(belgium_shape)
belgium_shape_f <- merge(belgium_shape_f, belgium_shape@data,by="id",all.x=TRUE)

Blues_plus <- colorRampPalette(brewer.pal(9,"Blues"))(11)
plot.bin.map.fisher <- ggplot(belgium_shape_f, aes(long,lat, group = group)) + geom_polygon(aes(fill = belgium_shape_f$class_fisher)) + theme_nothing(legend = TRUE) + labs(fill = "Fisher") + scale_fill_manual(values=Blues_plus, na.value = "white") 

#################################################
# 4) Refit the GAM with a binned spatial effect #
#################################################

DT.geo <- DT[,.(ID,NCLAIMS,EXP,COVER,FUEL,AGEph,POWER,SPLIT,PC)]
DT.geo <- merge(DT.geo,GAM.freq.LONGLAT,by.x="PC",by.y="pc",all.x = TRUE)
DT.geo[,GEO:=as.factor(cut(DT.geo$s, breaks = classint.fisher$brks, right = FALSE,include.lowest=TRUE,dig.lab=2))]

# Refit GAM
gam.freq.geo <- gam(NCLAIMS  ~  COVER + FUEL + POWER + SPLIT + s(AGEph) +
                      GEO, offset=log(EXP) , data = DT.geo, family = poisson(link = "log"))

plot.gam.freq.geo.ageph <- ggplot.gam(gam.freq.geo,DT.geo[, AGEph],"s(AGEph)","ageph",expression(hat(f)[1](ageph)))
plot.gam.freq.geo.ageph

# save(gam.freq.geo, file = "base_split_kbin=10.rda")
# load("base_split_kbin=10.rda")

AIC(gam.freq.geo)
BIC(gam.freq.geo)


#################################
# 5) Bin the continuous effects #
#################################

# Get GAM data which will serve as input for the trees
getGAMdata_single = function(GAMmodel,term,var,varname){
  pred = predict(GAMmodel,type = "terms",terms = term)
  DT_pred = data.table("x"=var,pred)[order(x)]
  names(DT_pred) = c("x","s")
  DT_unique = unique(DT_pred[,.(x,s)])
  DT_exp = DT_pred[,.N, by=x]
  GAM_data = merge(DT_unique,DT_exp,by="x")
  names(GAM_data) = c(varname,"s","exp")
  GAM_data = GAM_data[which(GAM_data$exp!=0)]
  return(GAM_data)
}

GAM.freq.AGEph = getGAMdata_single(gam.freq.geo,"s(AGEph)",DT.geo$AGEph,"ageph")

source("evtree.R")

# Set the control parameters
ctrl.freq = evtree.control(minbucket=0.05*nrow(DT),alpha=100,maxdepth=5)

# Fit the trees
evtree.freq.AGEph <- evtree(s ~ ageph,data = GAM.freq.AGEph,weights=exp,control=ctrl.freq)

# EVTREE splits
splits_evtree = function(evtreemodel,GAMvar,DTvar){
  preds=predict(evtreemodel,type="node")
  nodes=data.table("x"=GAMvar,"nodes"=preds)
  nodes[,change:=c(0,pmin(1,diff(nodes)))]
  splits_evtree=unique(c(min(DTvar),nodes$x[which(nodes$change==1)],max(DTvar)))
  return(splits_evtree)
}
splits2D_evtree = function(evtreemodel,GAMdata,GAMdata_X,GAMdata_Y){
  pred = predict(evtreemodel,GAMdata,type="response")
  values = data.table("X"=GAMdata_X,"Y"=GAMdata_Y,"pred"=pred)
  min = values[,lapply(.SD, min),by=pred]
  max = values[,lapply(.SD, max),by=pred]
  splits_2D_evtree = data.table("xmin"=min$X,"xmax"=max$X,"ymin"=min$Y,"ymax"=max$Y)
  return(splits_2D_evtree)
}
splits.freq.AGEph = splits_evtree(evtree.freq.AGEph,GAM.freq.AGEph$ageph,DT$AGEph)
splits.freq.AGEph

# Plot the bins
plot.bin.freq.ageph = ggplot.gam(gam.freq,DT[, AGEph],"s(AGEph)","ageph",expression(hat(f)[1](ageph))) + geom_vline(xintercept = splits.freq.AGEph)

######################################
# 6) Fit a GLM with binned variables #
######################################

# Create binned dataset
DT.freq.bin <- DT.geo[,.(ID,NCLAIMS,EXP,COVER,FUEL,POWER,SPLIT,GEO)]
DT.freq.bin[,AGEph:=cut(DT.geo$AGEph,splits.freq.AGEph,right=FALSE,include.lowest=TRUE)]

# Choose reference levels for the GLM (most exposure)
DT.freq.bin$COVER = relevel(DT.freq.bin$COVER,ref="MTPL") 
DT.freq.bin$FUEL = relevel(DT.freq.bin$FUEL,ref="Petrol")
DT.freq.bin$POWER = relevel(DT.freq.bin$POWER,ref="<66")
DT.freq.bin$SPLIT = relevel(DT.freq.bin$SPLIT,ref="Once")

# Bar plot of factor variables to check which age category has the most exposure
plot.eda.ageph <- ggplot.bar(DT.freq.bin,DT.freq.bin$AGEph,"ageph")
# Bar plot of factor variables to check which district has the most exposure
plot.eda.geo <- ggplot.bar(DT.geo,DT.geo$GEO,"spatial effect")

DT.freq.bin$AGEph = relevel(DT.freq.bin$AGEph,ref="[38,51)") 
DT.freq.bin$GEO = relevel(DT.freq.bin$GEO,ref="[-0.06,0.018)")

# Fit GLM
glm.freq <- gam(NCLAIMS ~ COVER + FUEL + POWER + SPLIT + AGEph + GEO, offset=log(EXP) , data = DT.freq.bin, family = poisson(link = "log"))

AIC(glm.freq)
BIC(glm.freq)

summary(glm.freq)
anova(glm.freq)

# Get GLM coefficients of single effects
GLM.freq.AGEph = getGAMdata_single(glm.freq,"AGEph",DT.geo$AGEph,"ageph")

# Shift GLM to get on scale of GAM
shift_GLM = function(GLMdata){
  shift = weighted.mean(GLMdata$s,GLMdata$exp)
  GLMdata = GLMdata[,s_shifted:=s-shift]
  return(GLMdata)
}

GLM.freq.AGEph <- shift_GLM(GLM.freq.AGEph)

# Get GLM coefficients of spatial effect
DT_maps <- data.table(PC=GAM.freq.LONGLAT$pc,GEO=cut(GAM.freq.LONGLAT$s, breaks = classint.fisher$brks, right = FALSE,include.lowest=TRUE,dig.lab = 2))
DT_maps[,c("COVER","FUEL","AGEph","POWER","SPLIT"):=DT.freq.bin[1,.(COVER,FUEL,AGEph,POWER,SPLIT)]]
GEO <- as.factor(round(predict(glm.freq,newdata = DT_maps,type = "terms",terms = "GEO"),digits=3))
DT_pred = data.table("PC"=DT_maps$PC,GEO)
belgium_shape <- readShapefile()
belgium_shape@data <- merge(belgium_shape@data,DT_pred,by.x="POSTCODE",by.y="PC",all.x=TRUE)
belgium_shape_f2 <- fortify(belgium_shape)
belgium_shape_f2 <- merge(belgium_shape_f2, belgium_shape@data,by="id",all.x=TRUE)

# Plot GLM effects on top of GAM effects
ggplot.glm <- function(model,variable,gam_term,xlabel,ylabel,glm_term,glmvar,splits){
  pred <- predict(model, type = "terms", se = TRUE)
  col_index <- which(colnames(pred$fit)==gam_term)
  x <- variable
  b <- pred$fit[, col_index]
  l <- pred$fit[, col_index] - qnorm(0.975) * pred$se.fit[, col_index]
  u <- pred$fit[, col_index] + qnorm(0.975) * pred$se.fit[, col_index]
  df <- unique(data.frame(x, b, l, u))
  p <- ggplot(df, aes(x = x))
  p <- p + geom_line(aes(y = b), size = 1, col="#003366")
  p <- p + geom_line(aes(y = l), size = 0.5, linetype = 2, col="#99CCFF")
  p <- p + geom_line(aes(y = u), size = 0.5, linetype = 2,col="#99CCFF")
  p <- p + geom_vline(xintercept = splits)
  p <- p + geom_point(data=glm_term,aes(x=glmvar,y=s_shifted),colour="#6666FF",size=3)
  p <- p + xlab(xlabel) + ylab(ylabel) + theme_bw()
  p
}

plot.glm.freq.ageph <- ggplot.glm(gam.freq,DT[, AGEph],"s(AGEph)","ageph",expression(hat(f)[1](ageph)),GLM.freq.AGEph,GLM.freq.AGEph$ageph,splits.freq.AGEph)  
plot.glm.freq.map <- ggplot(belgium_shape_f2, aes(long,lat, group = group)) + geom_polygon(aes(fill = belgium_shape_f2$GEO)) + theme_nothing(legend = TRUE) + labs(fill = "GLM") + scale_fill_manual(values=Blues_plus, na.value = "white") 



