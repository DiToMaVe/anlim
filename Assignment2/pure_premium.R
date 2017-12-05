
###################
# 8) Pure premium #
###################

# GAM components
DT.gam <- DT[,.(ID,NCLAIMS,COVER,FUEL,AGEph,POWER,SPLIT,LONG,LAT)]

gam.freq.pred <- predict(gam.freq, newdata = DT.gam, type = "response")
gam.sev.pred <- predict(gam.sev, newdata = DT.gam, type = "response",se.fit = TRUE)

# GLM components
DT.glm.freq <- DT.freq.bin[,EXP:=NULL]
DT.glm.sev <- DT[,.(ID,NCLAIMS,COVER,SPLIT)]
DT.glm.sev[,AGEph:=cut(DT$AGEph,splits.sev.AGEph,right=FALSE,include.lowest=TRUE)]
DT.glm.sev$AGEph[is.na(DT.glm.sev$AGEph)] <- levels(DT.glm.sev$AGEph)[length(levels(DT.glm.sev$AGEph))]

glm.freq.pred <- predict(glm.freq, newdata = DT.glm.freq, type = "response")
glm.sev.pred <- predict(glm.sev, newdata = DT.glm.sev, type = "response",se.fit = TRUE)

# Compare the obtained premiums
DT.Premium <- DT[,.(ID,NCLAIMS,COVER,FUEL,AGEph,POWER,SPLIT,LONG,LAT)]
DT.Premium[,"PPgam":=gam.freq.pred*exp(gam.sev.pred[[1]]+gam.sev.pred[[2]]^2/2)]
DT.Premium[,"PPglm":=glm.freq.pred*exp(glm.sev.pred[[1]]+glm.sev.pred[[2]]^2/2)]

sum(DT.Premium$PPgam)
sum(DT.Premium$PPglm)

fivenum(DT.Premium$PPgam)
fivenum(DT.Premium$PPglm)

DT.Premium.plot <- DT.Premium[,.(PPgam,PPglm)]
DT.Premium.plot <- melt(DT.Premium.plot)
ggplot(DT.Premium.plot,aes(x=variable, y=value, fill=variable)) + geom_boxplot() + theme_bw()


