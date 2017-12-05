########################
# 7) Severity modeling #
########################

# Choose reference levels for the GAM (most exposure)
DT.sev$COVER = relevel(DT.sev$COVER,ref="MTPL") 
DT.sev$FUEL = relevel(DT.sev$FUEL,ref="Petrol")
DT.sev$POWER = relevel(DT.sev$POWER,ref="<66")
DT.sev$SPLIT = relevel(DT.sev$SPLIT,ref="Once")
# 
# # Base model Gamma GAM
# gam.sev <- gam(AVG ~ COVER + s(AGEph),weights = NCLAIMS, data = DT.sev, family = Gamma(link = "identity"))
# summary(gam.sev)
# 
# # Base model LogNormal GAM
# gam.lognormal <- gam(log(AVG) ~ COVER + s(AGEph), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.lognormal)
# 
# # Base model LogNormal
# gam.sev1 <- gam(log(AVG) ~ COVER + FUEL + POWER + s(AGEph)+ s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev1)
# 
# # Base model + SPORT 
# gam.sev2 <- gam(log(AVG) ~ COVER + FUEL + POWER + SPORT + s(AGEph)+ s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev2)
# 
# # Base model + SPLIT 
# gam.sev3 <- gam(log(AVG) ~ COVER + FUEL + POWER + SPLIT + s(AGEph)+ s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev3)
# 
# # Base model + SPORT 
# gam.sev4 <- gam(log(AVG) ~ COVER + FUEL + POWER + SPORT + SPLIT + s(AGEph)+ s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev4)
# 
# # Base model + SPORT + SPLIT 
# gam.sev5 <- gam(log(AVG) ~ COVER + FUEL + POWER + SPORT + SPLIT + s(AGEph)+ s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev5)
# 
# # Base model COVER + FUEL + POWER + SPLIT + s(AGEph) + s(AGEph,by=POWER, bs="fs") + s(LONG,LAT) + SPLIT + AGEph_POWER
# gam.sev6 <- gam(log(AVG) ~  COVER + FUEL + POWER + SPLIT + s(AGEph) + s(AGEph,by=POWER, bs="fs") + s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev6)
# 
# # Base model + SPLIT + AGEph_SPLIT
# gam.sev7 <- gam(log(AVG) ~  COVER + FUEL + POWER + SPLIT + s(AGEph) + s(AGEph,by=SPLIT, bs="fs") + s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev7)
# 
# # Base model + SPLIT + AGEph_SPORT
# gam.sev8 <- gam(log(AVG) ~  COVER + FUEL + POWER + SPLIT + s(AGEph) + s(AGEph,by=SPORT, bs="fs") + s(LONG,LAT), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev8)
# 
# # Base model + SPLIT - LONG_LAT 
# gam.sev9 <- gam(log(AVG) ~ COVER + FUEL + POWER + SPLIT + s(AGEph), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev9)
# 
# # Base model + SPLIT - LONG_LAT - AGEph
# gam.sev10 <- gam(log(AVG) ~ COVER + FUEL + POWER + SPLIT , weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev10)
# 
# # Base model + SPLIT - LONG_LAT - POWER
# gam.sev11 <- gam(log(AVG) ~ COVER + FUEL + SPLIT + s(AGEph), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev11)
# 
# # Base model + SPLIT - LONG_LAT - POWER - FUEL
# gam.sev12 <- gam(log(AVG) ~ COVER + SPLIT + s(AGEph), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev12)
# 
# # Base model - LONG_LAT - POWER - FUEL
# gam.sev13 <- gam(log(AVG) ~ COVER + s(AGEph), weights = NCLAIMS, data = DT.sev, family = gaussian(link = "identity"))
# summary(gam.sev13)
# 
# save(gam.sev1, file = "base_sev.rda")
# save(gam.sev2, file = "base_sport_sev.rda")
# save(gam.sev3, file = "base_split_sev.rda")
# save(gam.sev4, file = "base_split_sport_sev.rda")
# save(gam.sev5, file = "base_agepower_sev.rda")
# save(gam.sev6, file = "base_split_agepower_sev.rda")
# save(gam.sev7, file = "base_split_agesplit_sev.rda")
# save(gam.sev8, file = "base_split_agesport_sev.rda")
# save(gam.sev9, file = "base_split__longlat_sev.rda")
# save(gam.sev10, file = "base_split__longlat__ageph_sev.rda")
# save(gam.sev11, file = "base_split__longlat__power_sev.rda")
# save(gam.sev12, file = "base_split__longlat__power__fuel_sev.rda")
# save(gam.sev13, file = "base__longlat__power__fuel_sev.rda")
# 
# load("base_sev.rda")
# load("base_sport_sev.rda")
# load("base_split_sev.rda")
# load("base_split_sport_sev.rda")
# load("base_agepower_sev.rda")
# load("base_split_agepower_sev.rda")
# load("base_split_agesplit_sev.rda")
# load("base_split_agesport_sev.rda")
# load("base_split__longlat_sev.rda")
# load("base_split__longlat_sev__ageph.rda")
# load("base_split__longlat_sev__power.rda")
# load("base_split__longlat_sev__power__fuel.rda")
# load("base__longlat_sev__power__fuel.rda")
# 
# AIC(gam.sev1)
# AIC(gam.sev2)
# AIC(gam.sev3)
# AIC(gam.sev4)
# AIC(gam.sev5)
# AIC(gam.sev6)
# AIC(gam.sev7)
# AIC(gam.sev8)
# AIC(gam.sev9)
# AIC(gam.sev10)
# AIC(gam.sev11)
# AIC(gam.sev12)
# AIC(gam.sev13)
# 
# BIC(gam.sev1)
# BIC(gam.sev2)
# BIC(gam.sev3)
# BIC(gam.sev4)
# BIC(gam.sev5)
# BIC(gam.sev6)
# BIC(gam.sev7)
# BIC(gam.sev8)
# BIC(gam.sev9)
# BIC(gam.sev10)
# BIC(gam.sev11)
# BIC(gam.sev12)
# BIC(gam.sev13)

# Optimal model: COVER + SPLIT + s(AGEph)
load("base_split__longlat_sev__power__fuel.rda")
gam.sev <- gam.sev12
plot.gam.sev.ageph <- ggplot.gam(gam.sev,DT.sev[, AGEph],"s(AGEph)","ageph",expression(hat(g)[1](ageph)))

# Bin the continuous risk factors
GAM.sev.AGEph = getGAMdata_single(gam.sev,"s(AGEph)",DT.sev$AGEph,"ageph")

ctrl.sev <- evtree.control(minbucket=0.06*nrow(DT.sev),alpha=90,maxdepth=5)

evtree.sev.AGEph <- evtree(s ~ ageph,data = GAM.sev.AGEph,weights=exp,control=ctrl.sev)

splits.sev.AGEph <- splits_evtree(evtree.sev.AGEph,GAM.sev.AGEph$ageph,DT.sev$AGEph)

plot.bin.sev.ageph <- ggplot.gam(gam.sev,DT.sev[, AGEph],"s(AGEph)","ageph",expression(hat(g)[1](ageph))) + geom_vline(xintercept = splits.sev.AGEph)


# Fit a GLM

# Choose reference levels for the GLM (most exposure)
DT.sev.bin = DT.sev[,.(ID,AVG,NCLAIMS,COVER,SPLIT)]
DT.sev.bin[,AGEph:=cut(DT.sev$AGEph,splits.sev.AGEph,right=FALSE,include.lowest=TRUE)]
DT.sev.bin$COVER = relevel(DT.sev.bin$COVER,ref="MTPL") 
DT.sev.bin$SPLIT = relevel(DT.sev.bin$SPLIT,ref="Once")

# Bar plot of factor variables to check which age category has the most exposure
plot.eda.sev.ageph <- ggplot.bar(DT.sev.bin,DT.sev.bin$AGEph,"AGEph")

DT.sev.bin$AGEph = relevel(DT.sev.bin$AGEph,ref="[42,64)")

glm.sev <- gam(log(AVG) ~ COVER + SPLIT + AGEph, weights = NCLAIMS , data = DT.sev.bin, family = gaussian(link = "identity"))

GLM.sev.AGEph <- getGAMdata_single(glm.sev,"AGEph",DT.sev$AGEph,"ageph")
GLM.sev.AGEph <- shift_GLM(GLM.sev.AGEph)

plot.glm.sev.ageph <- ggplot.glm(gam.sev,DT.sev[, AGEph],"s(AGEph)","ageph",expression(hat(g)[1](ageph)),GLM.sev.AGEph,GLM.sev.AGEph$ageph,splits.sev.AGEph)  

AIC(glm.sev)
BIC(glm.sev)


