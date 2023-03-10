##################################################################################
##################################################################################
###############                 GA ponds HABS Analysis              ##############
##################################################################################
##################################################################################


###################################################################
####### 1) read master file and make any necessary edits  #########
###################################################################


library(lme4)
library(lmerTest)
library(nlme)
library(gplots)
library(lubridate)
library(devtools)
library(pracma)

setwd("C:/Users/straussa/Documents/Research/Strauss Lab - UGA/GA zoop surveys/Whitehall data")
alldata <- read.csv("GAponds.master.summary.csv", na.strings=c("", "NA"))


###################################################################
###################### 2) prepare egg ratio summary  ##############
###################################################################

egg.summary <- data.frame(taxa = c("dam", "dla", "dpa", "sim"), healthy.eggs=NA, 
                          healthy.eggs.se=NA, l.eggs=NA, l.eggs.se=NA)
# fill in for D ambigua
egg.summary[egg.summary$taxa=="dam",]$healthy.eggs <- mean(alldata$dam.healthy.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="dam",]$healthy.eggs.se <- sd(alldata$dam.healthy.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$dam.healthy.eggs),]))
egg.summary[egg.summary$taxa=="dam",]$l.eggs <- mean(alldata$dam.l.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="dam",]$l.eggs.se <- sd(alldata$dam.l.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$dam.l.eggs),]))

# fill in for D laevis
egg.summary[egg.summary$taxa=="dla",]$healthy.eggs <- mean(alldata$dla.healthy.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="dla",]$healthy.eggs.se <- sd(alldata$dla.healthy.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$dla.healthy.eggs),]))
egg.summary[egg.summary$taxa=="dla",]$l.eggs <- mean(alldata$dla.l.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="dla",]$l.eggs.se <- sd(alldata$dla.l.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$dla.l.eggs),]))

# fill in for D parvula
egg.summary[egg.summary$taxa=="dpa",]$healthy.eggs <- mean(alldata$dpa.healthy.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="dpa",]$healthy.eggs.se <- sd(alldata$dpa.healthy.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$dpa.healthy.eggs),]))
egg.summary[egg.summary$taxa=="dpa",]$l.eggs <- mean(alldata$dpa.l.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="dpa",]$l.eggs.se <- sd(alldata$dpa.l.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$dpa.l.eggs),]))

# fill in for Simocephalus
egg.summary[egg.summary$taxa=="sim",]$healthy.eggs <- mean(alldata$sim.healthy.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="sim",]$healthy.eggs.se <- sd(alldata$sim.healthy.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$sim.healthy.eggs),]))
egg.summary[egg.summary$taxa=="sim",]$l.eggs <- mean(alldata$sim.l.eggs, na.rm=T)
egg.summary[egg.summary$taxa=="sim",]$l.eggs.se <- sd(alldata$sim.l.eggs, na.rm=T)/
  sqrt(nrow(alldata[!is.na(alldata$sim.l.eggs),]))

egg.summary$perc.red <- 1-(egg.summary$l.eggs/egg.summary$healthy.eggs)
egg.summary



###################################################################
###################### 3) plot egg ratio results  #################
###################################################################


setwd("C:/Users/straussa/Documents/Research/Strauss Lab - UGA/GA zoop surveys/Whitehall data/figures")
png("virulence.png", width = 3, height = 3, res = 600, units='in')

par(mfrow=c(1,1),mar=c(2.5,2.5,2,1), oma=c(0,0,0,0))

barplot(egg.summary$healthy.eggs, space=c(1,2,2,2), col="light grey", border=NA,
        ylim=c(0,3), xlim=c(0.5,9)) # to include simocephalus: xlim=c(0.5,12))
plotCI(c(1.5,4.5,7.5,10.5), egg.summary$healthy.eggs, uiw=egg.summary$healthy.eggs.se,
       add=T, sfrac=0, gap=0, pch=NA)
barplot(egg.summary$l.eggs, space=c(2,2,2,2),col="red", border=NA, add=T)
plotCI(c(2.5,5.5,8.5,11.5), egg.summary$l.eggs, uiw=egg.summary$l.eggs.se,
       add=T, sfrac=0, gap=0, pch=NA)

dev.off()




###################################################################
######################## 4) Egg ratio statistics  #################
###################################################################

##############################
# prepping data for D ambigua:

dam.healthy.eggs <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                       "year", "dam.healthy.eggs", "dam.nhealthymoms")]
dam.healthy.eggs$infect <- "no"
dam.l.eggs <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                               "year", "dam.l.eggs", "dam.nlmoms")]
dam.l.eggs$infect <- "yes"
colnames(dam.healthy.eggs) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                               "year", "eggs", "moms", "infect")
colnames(dam.l.eggs) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                                "year", "eggs", "moms", "infect")
dam.eggs <- rbind(dam.l.eggs, dam.healthy.eggs)
dam.eggs <- dam.eggs[!is.na(dam.eggs$eggs),]
# seasonal effects are not strong here, so I'm leaning towards removing them from the models
# dam.eggs$season <- "winter"
# dam.eggs$season <- ifelse(dam.eggs$jday > 79 & dam.eggs$jday < 172, "spring", dam.eggs$season)
# dam.eggs$season <- ifelse(dam.eggs$jday > 171 & dam.eggs$jday < 264, "summer", dam.eggs$season)
# dam.eggs$season <- ifelse(dam.eggs$jday > 263 & dam.eggs$jday < 355, "fall", dam.eggs$season)
# for inclusion as a random effect; don't want to force linear relationship:
dam.eggs$jday_cum <- as.factor(dam.eggs$jday_cum)

##############################
# prepping data for D parvula:

dpa.healthy.eggs <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                               "year", "dpa.healthy.eggs", "dpa.nhealthymoms")]
dpa.healthy.eggs$infect <- "no"
dpa.l.eggs <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                         "year", "dpa.l.eggs", "dpa.nlmoms")]
dpa.l.eggs$infect <- "yes"
colnames(dpa.healthy.eggs) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                                "year", "eggs", "moms", "infect")
colnames(dpa.l.eggs) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                          "year", "eggs", "moms", "infect")
dpa.eggs <- rbind(dpa.l.eggs, dpa.healthy.eggs)
dpa.eggs <- dpa.eggs[!is.na(dpa.eggs$eggs),]
dpa.eggs$jday_cum <- as.factor(dpa.eggs$jday_cum)


##############################
# prepping data for D laevis:

dla.healthy.eggs <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                               "year", "dla.healthy.eggs", "dla.nhealthymoms")]
dla.healthy.eggs$infect <- "no"
dla.l.eggs <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                         "year", "dla.l.eggs", "dla.nlmoms")]
dla.l.eggs$infect <- "yes"
colnames(dla.healthy.eggs) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                                "year", "eggs", "moms", "infect")
colnames(dla.l.eggs) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                          "year", "eggs", "moms", "infect")
dla.eggs <- rbind(dla.l.eggs, dla.healthy.eggs)
dla.eggs <- dla.eggs[!is.na(dla.eggs$eggs),]
dla.eggs$jday_cum <- as.factor(dla.eggs$jday_cum)


###############################
# check response distributions:

hist(dam.eggs$eggs)
hist(sqrt(dam.eggs$eggs))
hist(dam.eggs$eggs^(1/2))
hist(dam.eggs$eggs^(1/3))
hist(log(dam.eggs$eggs+1))

hist(dpa.eggs$eggs)
hist(sqrt(dpa.eggs$eggs))
hist(dpa.eggs$eggs^(1/2))
hist(dpa.eggs$eggs^(1/3))
hist(log(dpa.eggs$eggs+1))

hist(dla.eggs$eggs)
hist(sqrt(dla.eggs$eggs))
hist(dla.eggs$eggs^(1/2))
hist(dla.eggs$eggs^(1/3))
hist(log(dla.eggs$eggs+1))

##################
# fit the models:

dam.eggs.mod1 <- lme(eggs^(1/3) ~ infect + year, random= list(~1|lake_id, ~1|jday_cum), 
                     data=dam.eggs, na.action=na.omit, weights = ~ I(1/moms))
summary(dam.eggs.mod1)
hist(resid(dam.eggs.mod1, type="p"))
plot(dam.eggs.mod1, form = resid(., type = "p") ~ fitted(.) | infect, abline=0)
shapiro.test(resid(dam.eggs.mod1, type="p")) # hmm, not great. 

dpa.eggs.mod1 <- lme(eggs^(1/3) ~ infect + year, random= list(~1|lake_id, ~1|jday_cum), 
                     data=dpa.eggs, na.action=na.omit, weights = ~ I(1/moms))
summary(dpa.eggs.mod1)
hist(resid(dpa.eggs.mod1, type="p"))
plot(dpa.eggs.mod1, form = resid(., type = "p") ~ fitted(.) | infect, abline=0)
shapiro.test(resid(dpa.eggs.mod1, type="p")) # hmm, not great. 

dla.eggs.mod1 <- lme(eggs^(1/3) ~ infect + year, random= list(~1|lake_id, ~1|jday_cum), 
                     data=dla.eggs, na.action=na.omit, weights = ~ I(1/moms))
summary(dla.eggs.mod1)
hist(resid(dla.eggs.mod1, type="p"))
plot(dla.eggs.mod1, form = resid(., type = "p") ~ fitted(.) | infect, abline=0)
shapiro.test(resid(dla.eggs.mod1, type="p")) # hmm, not great. 

# # negative binomial model seems good in theory, but would require troubleshooting 
# library(glmmTMB)
# dam.eggs.mod2 <- glmmTMB(eggs ~ season + infect + year + (1|lake_id) + (1|jday_cum), 
#                      data=dam.eggs, na.action=na.omit, weights = moms,
#                      ziformula = ~0, family = nbinom2)
# warnings()
# summary(dam.eggs.mod2)




###################################################################
################### 5) plot infection prevalence  #################
###################################################################

# subset by pond and remove observations with fewer than x individuals counted 
mincount <- 50 # (value for threshold over which we feel good about interpreting prevalence)

# ambigua ponds:
deans.trim.dam <- alldata[alldata$lake_id=="Deans" & alldata$dam.ncount > mincount,]
vip.trim.dam <- alldata[alldata$lake_id=="VIP" & alldata$dam.ncount > mincount,]
bs.trim.dam <- alldata[alldata$lake_id=="Big Sister" & alldata$dam.ncount > mincount,]
deer.trim.dam <- alldata[alldata$lake_id=="Deer" & alldata$dam.ncount > mincount,]
#laevis ponds:
deer.trim.dla <- alldata[alldata$lake_id=="Deer" & alldata$dla.ncount > mincount,]
cat.trim.dla <- alldata[alldata$lake_id=="Catfish" & alldata$dla.ncount > mincount,]
NN3.trim.dla <- alldata[alldata$lake_id=="NN3" & alldata$dla.ncount > mincount,]
sis1.trim.dla <- alldata[alldata$lake_id=="Sister 1" & alldata$dla.ncount > mincount,]
# parvula:
deans.trim.dpa <- alldata[alldata$lake_id=="Deans" & alldata$dpa.ncount > mincount,]
nrow(deans.trim.dpa)



setwd("C:/Users/straussa/Documents/Research/Strauss Lab - UGA/GA zoop surveys/Whitehall data/figures")
png("microsporidian.time.series.png", width = 6, height = 4, res = 600, units='in')

par(mfrow=c(3,1),mar=c(0.5,0.5,0,0), oma=c(2,3,0.5,0.5))

plot(cat.trim.dla$jday_cum, cat.trim.dla$dla.lprev, type="o", 
     ylim=c(0,.8), xlim=c(50,550), xaxt="n", yaxt="n", ann=F, col="black", lwd=2)
plotCI(cat.trim.dla$jday_cum, cat.trim.dla$dla.lprev, uiw=cat.trim.dla$dla.lprev.se,
       add=T, col="black", lty=2, gap=0, sfrac=0)
points(NN3.trim.dla$jday_cum, NN3.trim.dla$dla.lprev, col="blue", type="o", lwd=2)
plotCI(NN3.trim.dla$jday_cum, NN3.trim.dla$dla.lprev, uiw=NN3.trim.dla$dla.lprev.se,
       add=T, col="blue", lty=2, gap=0, sfrac=0)
points(sis1.trim.dla$jday_cum, sis1.trim.dla$dla.lprev, col="forestgreen", type="o", lwd=2)
plotCI(sis1.trim.dla$jday_cum, sis1.trim.dla$dla.lprev, uiw=sis1.trim.dla$dla.lprev.se,
       add=T, col="forestgreen", lty=2, gap=0, sfrac=0)
points(deer.trim.dla$jday_cum, deer.trim.dla$dla.lprev, col="orange", type="o", lwd=2)
plotCI(deer.trim.dla$jday_cum, deer.trim.dla$dla.lprev, uiw=deer.trim.dla$dla.lprev.se,
       add=T, col="orange", lty=2, gap=0, sfrac=0)
#axis(side=1, at=c(90,180,270,360,450,540))
axis(side=2, at=c(0,.35,.7), las=1)
abline(v=32, lty=2); abline(v=151, lty=2);  abline(v=32+365, lty=2); abline(v=151+365, lty=2)
# showing Feb 1 (32) - May 31 (151)

plot(deans.trim.dam$jday_cum, deans.trim.dam$dam.lprev, type="o", 
     ylim=c(0,.2), xlim=c(50,550), xaxt="n", yaxt="n", ann=F, col="red", lwd=2)
plotCI(deans.trim.dam$jday_cum, deans.trim.dam$dam.lprev, uiw=deans.trim.dam$dam.lprev.se,
       add=T, col="red", lty=2, gap=0, sfrac=0)
points(vip.trim.dam$jday_cum, vip.trim.dam$dam.lprev, col="purple", type="o", lwd=2)
plotCI(vip.trim.dam$jday_cum, vip.trim.dam$dam.lprev, uiw=vip.trim.dam$dam.lprev.se,
       add=T, col="purple", lty=2, gap=0, sfrac=0)
points(bs.trim.dam$jday_cum, bs.trim.dam$dam.lprev, col="deepskyblue3", type="o", lwd=2)
plotCI(bs.trim.dam$jday_cum, bs.trim.dam$dam.lprev, uiw=bs.trim.dam$dam.lprev.se,
       add=T, col="deepskyblue3", lty=2, gap=0, sfrac=0)
points(deer.trim.dam$jday_cum, deer.trim.dam$dam.lprev, col="orange", type="o", lwd=2)
plotCI(deer.trim.dam$jday_cum, deer.trim.dam$dam.lprev, uiw=deer.trim.dam$dam.lprev.se,
       add=T, col="orange", lty=2, gap=0, sfrac=0)
abline(v=32, lty=2); abline(v=151, lty=2);  abline(v=32+365, lty=2); abline(v=151+365, lty=2)
# showing Feb 1 (32) - May 31 (151)
#axis(side=1, at=c(90,180,270,360,450,540))
axis(side=2, at=c(0,.1,.2), las=1)


plot(deans.trim.dpa$jday_cum, deans.trim.dpa$dpa.lprev, type="o", 
     ylim=c(0,0.08), xlim=c(50,550), xaxt="n", yaxt="n", ann=F, col="red", lwd=2)
plotCI(deans.trim.dpa$jday_cum, deans.trim.dpa$dpa.lprev, uiw=deans.trim.dpa$dpa.lprev.se,
       add=T, col="red", lty=2, gap=0, sfrac=0)
abline(v=32, lty=2); abline(v=151, lty=2);  abline(v=32+365, lty=2); abline(v=151+365, lty=2)
# showing Feb 1 (32) - May 31 (151)
axis(side=1, at=c(90,180,270,360,450,540))
axis(side=2, at=c(0,.03,.06), las=1)



dev.off()


###################################################################
######### 6) statistics for seasonality of outbreaks  #############
###################################################################

# prep data frame to ask how risk depends on taxa and season:
# D ambigua:
dam.lprev <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                               "year", "dam.ncount", "dam.lprev")]
dam.lprev$tax <- "dam"
colnames(dam.lprev) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                         "year", "ncount", "lprev", "tax")

# D. parvula
dpa.lprev <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                        "year", "dpa.ncount", "dpa.lprev")]
dpa.lprev$tax <- "dpa"
colnames(dpa.lprev) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                         "year", "ncount", "lprev", "tax")

# and D. laevis:
dla.lprev <- alldata[,c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                        "year", "dla.ncount", "dla.lprev")]
dla.lprev$tax <- "dla"
colnames(dla.lprev) <- c("lakeday", "jday_cum", "date", "lake_id", "jday", 
                         "year", "ncount", "lprev", "tax")

lprev.summary <- rbind(dam.lprev, dpa.lprev, dla.lprev)
lprev.summary <- lprev.summary[!is.na(lprev.summary$lprev),]
# # add season
lprev.summary$month <- month(as.POSIXlt(lprev.summary$date, format="%m/%d/%Y"))
# official definition of "season" based on julian day: (but not as useful since lots of infection
# happen right at the 'official' boundary of winter and spring) 
# lprev.summary$season <- "winter"
# lprev.summary$season <- ifelse(lprev.summary$jday > 79 & lprev.summary$jday < 172, "spring", lprev.summary$season)
# lprev.summary$season <- ifelse(lprev.summary$jday > 171 & lprev.summary$jday < 264, "summer", lprev.summary$season)
# lprev.summary$season <- ifelse(lprev.summary$jday > 263 & lprev.summary$jday < 355, "fall", lprev.summary$season)
lprev.summary$season <- ifelse(lprev.summary$month %in% c(2,3,4,5)==T, "spring",
                           ifelse(lprev.summary$month %in% c(6,7,8)==T, "summer",
                                  ifelse(lprev.summary$month %in% c(9,10,11)==T, "fall", "winter")))


lprev.mod1 <- glmer(lprev ~ season + tax + 
                (1|lake_id) + (1|jday_cum),
              family=binomial(link = "logit"), weights = ncount, 
              data=lprev.summary, na.action=na.omit,
              glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(lprev.mod1)


lprev.mod2 <- glmer(lprev ~ season * tax + 
                      (1|lake_id) + (1|jday_cum),
                    family=binomial(link = "logit"), weights = ncount, 
                    data=lprev.summary, na.action=na.omit,
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(lprev.mod2)
# interaction leads to no new significant terms, so stick with simpler model



# https://stats.oarc.ucla.edu/r/dae/logit-regression/
devtools::install_github("remkoduursma/bootpredictlme4")
library(bootpredictlme4)

predicted <- data.frame(season = rep(unique(lprev.summary$season),3))
predicted$tax <- c(rep("dam",4), rep("dla",4), rep("dpa",4))
mod1predict <- predict(lprev.mod1, newdata = predicted, type = "response", 
                       re.form=NA, se.fit=T, nsim=10)
predicted$mod1.est <- mod1predict$fit
predicted$mod1.se <- mod1predict$se.boot

# trying to bootstrap the model with interaction between tax and season
# leads lots of convergence warnings... and results don't look that interesting anyway

par(mfrow=c(1,1),mar=c(0.5,0.5,0,0), oma=c(2,3,0.5,0.5))

barplot(predicted$mod1.est, ylim=c(0,0.06), space= rep(c(1,0.2,0.2,0.2),3),
        col=c("light grey"), border=NA)
plotCI(c(1.5,2.7,3.9,5.1, 7.1,8.3,9.5,10.7, 12.7,13.9,15.1,16.3),
       predicted$mod1.est, predicted$mod1.se, add=T,
       sfrac=0, gap=0, pch=NA)

# barplot(predicted$mod3.est, ylim=c(0,0.12), space= rep(c(1,0.2,0.2,0.2),3),
#         col=c("light grey"), border=NA)
# plotCI(c(1.5,2.7,3.9,5.1, 7.1,8.3,9.5,10.7, 12.7,13.9,15.1,16.3),
#        predicted$mod3.est, predicted$mod3.se, add=T,
#        sfrac=0, gap=0, pch=NA)



###################################################################
### 7) time series for disease outbreaks vs algae in each pond  ###
###################################################################

cat <- alldata[alldata$lake_id=="Catfish",]
vip <- alldata[alldata$lake_id=="VIP",]
deer <- alldata[alldata$lake_id=="Deer",]
deans <- alldata[alldata$lake_id=="Deans",]
bs <- alldata[alldata$lake_id=="Big Sister",]
sis1 <- alldata[alldata$lake_id=="Sister 1",]
NN3 <- alldata[alldata$lake_id=="NN3",]

par(mfrow=c(2,4),mar=c(0.5,0.5,0,0), oma=c(2,3,0.5,0.5))
par(mfrow=c(1,1),mar=c(0.5,0.5,0,0), oma=c(2,3,0.5,0.5))

# catfish, dla:
plot(cat[!is.na(cat$t.min),]$jday_cum, cat[!is.na(cat$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(cat[!is.na(cat$dla.lprev),]$jday_cum, cat[!is.na(cat$dla.lprev),]$dla.lprev*30, col="red", type="o")
points(cat[!is.na(cat$chl.tot),]$jday_cum, cat[!is.na(cat$chl.tot),]$chl.tot*5, col="green", type="o")
points(cat[!is.na(cat$pc.tot),]$jday_cum, cat[!is.na(cat$pc.tot),]$pc.tot*2, col="orange", type="o")

# deer, dla:
plot(deer[!is.na(deer$t.min),]$jday_cum, deer[!is.na(deer$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(deer[!is.na(deer$dla.lprev),]$jday_cum, deer[!is.na(deer$dla.lprev),]$dla.lprev*30, col="red", type="o")
points(deer[!is.na(deer$chl.tot),]$jday_cum, deer[!is.na(deer$chl.tot),]$chl.tot*5, col="green", type="o")
points(deer[!is.na(deer$pc.tot),]$jday_cum, deer[!is.na(deer$pc.tot),]$pc.tot*2, col="orange", type="o")

# NN3, dla:
plot(NN3[!is.na(NN3$t.min),]$jday_cum, NN3[!is.na(NN3$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(NN3[!is.na(NN3$dla.lprev),]$jday_cum, NN3[!is.na(NN3$dla.lprev),]$dla.lprev*30, col="red", type="o")
points(NN3[!is.na(NN3$chl.tot),]$jday_cum, NN3[!is.na(NN3$chl.tot),]$chl.tot*5, col="green", type="o")
points(NN3[!is.na(NN3$pc.tot),]$jday_cum, NN3[!is.na(NN3$pc.tot),]$pc.tot*2, col="orange", type="o")

# SIS1, dla:
plot(sis1[!is.na(sis1$t.min),]$jday_cum, sis1[!is.na(sis1$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(sis1[!is.na(sis1$dla.lprev),]$jday_cum, sis1[!is.na(sis1$dla.lprev),]$dla.lprev*30, col="red", type="o")
points(sis1[!is.na(sis1$chl.tot),]$jday_cum, sis1[!is.na(sis1$chl.tot),]$chl.tot*5, col="green", type="o")
points(sis1[!is.na(sis1$pc.tot),]$jday_cum, sis1[!is.na(sis1$pc.tot),]$pc.tot*2, col="orange", type="o")

# deer, dam:
plot(deer[!is.na(deer$t.min),]$jday_cum, deer[!is.na(deer$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(deer[!is.na(deer$dam.lprev),]$jday_cum, deer[!is.na(deer$dam.lprev),]$dam.lprev*30, col="red", type="o")
points(deer[!is.na(deer$chl.tot),]$jday_cum, deer[!is.na(deer$chl.tot),]$chl.tot*5, col="green", type="o")
points(deer[!is.na(deer$pc.tot),]$jday_cum, deer[!is.na(deer$pc.tot),]$pc.tot*2, col="orange", type="o")

# vip, dam:
plot(vip[!is.na(vip$t.min),]$jday_cum, vip[!is.na(vip$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(vip[!is.na(vip$dam.lprev),]$jday_cum, vip[!is.na(vip$dam.lprev),]$dam.lprev*30, col="red", type="o")
points(vip[!is.na(vip$chl.tot),]$jday_cum, vip[!is.na(vip$chl.tot),]$chl.tot*5, col="green", type="o")
points(vip[!is.na(vip$pc.tot),]$jday_cum, vip[!is.na(vip$pc.tot),]$pc.tot*2, col="orange", type="o")

# bs, dam:
plot(bs[!is.na(bs$t.min),]$jday_cum, bs[!is.na(bs$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(bs[!is.na(bs$dam.lprev),]$jday_cum, bs[!is.na(bs$dam.lprev),]$dam.lprev*30, col="red", type="o")
points(bs[!is.na(bs$chl.tot),]$jday_cum, bs[!is.na(bs$chl.tot),]$chl.tot*5, col="green", type="o")
points(bs[!is.na(bs$pc.tot),]$jday_cum, bs[!is.na(bs$pc.tot),]$pc.tot*2, col="orange", type="o")

# deans, dam:
plot(deans[!is.na(deans$t.min),]$jday_cum, deans[!is.na(deans$t.min),]$t.min, col="black", type="o", xlim=c(80,630), ylim=c(0,30))
points(deans[!is.na(deans$dam.lprev),]$jday_cum, deans[!is.na(deans$dam.lprev),]$dam.lprev*30, col="red", type="o")
points(deans[!is.na(deans$chl.tot),]$jday_cum, deans[!is.na(deans$chl.tot),]$chl.tot*5, col="green", type="o")
points(deans[!is.na(deans$pc.tot),]$jday_cum, deans[!is.na(deans$pc.tot),]$pc.tot*2, col="orange", type="o")



###################################################################
####### 8) specific timing of outbreaks within the spring  ########
###################################################################

# compile 'main' lprev based on dominant taxa in the pond
core.lakes <- c("VIP", "Catfish", "Deer", "Deans", "Big Sister", "Sister 1", "NN3")
dam.lakes <- c("VIP", "Deans", "Big Sister")
dla.lakes <- c("Catfish", "Deer", "NN3", "Sis1")

outbreak.data <- alldata[alldata$lake_id %in% core.lakes==T,]
outbreak.data$main.lprev <- ifelse(outbreak.data$lake_id %in% dam.lakes ==T, outbreak.data$dam.lprev, 
                                   outbreak.data$dla.lprev)
outbreak.data$main.ncount <- ifelse(outbreak.data$lake_id %in% dam.lakes ==T, outbreak.data$dam.ncount, 
                                   outbreak.data$dla.ncount)
outbreak.data$lake_yr <- paste(outbreak.data$lake_id, outbreak.data$year, sep="_")

par(mfrow=c(1,1),mar=c(2,2,3,1), oma=c(0,0,0,0))

lakeyears <- unique(outbreak.data$lake_yr)

i=1
outbreak.summary <- data.frame(lake_yr = lakeyears)
for(i in 1:length(lakeyears)){
    temp <- outbreak.data[outbreak.data$lake_yr==lakeyears[i],]
    outbreak.summary$lake_id[i] <- temp$lake_id[1]
    outbreak.summary$year[i] <- temp$year[1]
    
    temp.inf <- temp[is.na(temp$main.ncount)==F & temp$main.ncount>50,]
    max.prev <- max(temp.inf$main.lprev)
    outbreak.summary$max.prev[i] <- max.prev
    
    outbreak.summary$day.max.prev[i] <- min(temp[temp$main.lprev==max.prev,]$jday, na.rm=T)
    outbreak.summary$tmin.max.prev[i] <- mean(temp[temp$main.lprev==max.prev,]$t.min, na.rm=T)
    outbreak.summary$tsurf.max.prev[i] <- mean(temp[temp$main.lprev==max.prev,]$t.surf, na.rm=T)
  
    outbreak.summary$first.day.over.15[i] <- min(temp[temp$t.min>15,]$jday, na.rm=T)
    min.temp <- min(temp$t.min, na.rm=T)
    outbreak.summary$coldest.temp[i] <- min.temp
    outbreak.summary$coldest.day[i] <- min(temp[temp$t.min == min.temp,]$jday, na.rm=T)
    
    # integrate degree days until jday 151
    warmperiod <- temp[temp$jday <152 & is.na(temp$t.min)==F,]
    warmperiod <- warmperiod[order(warmperiod$jday),]
    plot(warmperiod$jday, warmperiod$t.min, type="o", main=paste(lakeyears[i]), ylim=c(4,25))
    outbreak.summary$degree.days[i] <- trapz(warmperiod$jday, warmperiod$t.min)
    
    # and just days 100 - 151, since didn't start sampling in 2021 until ~day 100
    warmperiod2 <- temp[temp$jday < 152 & temp$jday  > 100 & is.na(temp$t.min)==F,]
    warmperiod2 <- warmperiod2[order(warmperiod2$jday),]
    #plot(warmperiod2$jday, warmperiod2$t.min, type="o", main=paste(lakeyears[i]), ylim=c(4,25))
    outbreak.summary$degree.days2[i] <- trapz(warmperiod2$jday, warmperiod2$t.min)
    
    
  }

# lakes without outbreak get NA for day max prev
outbreak.summary[outbreak.summary$max.prev==0,]$day.max.prev <- NA
outbreak.summary[outbreak.summary$max.prev==0,]$tmin.max.prev <- NA
outbreak.summary[outbreak.summary$max.prev==0,]$tsurf.max.prev <- NA

# catfish peak in 2022 happened on week we didn't do YSI, so fill in as average of YSI before and after
outbreak.summary[outbreak.summary$lake_yr=="Catfish_2022",]$tmin.max.prev <- 
  mean(alldata[alldata$lake_id == "Catfish" & alldata$jday_cum > 429 & alldata$jday_cum < 454 & is.na(alldata$jday)==F,]$t.min, na.rm=T)
outbreak.summary[outbreak.summary$lake_yr=="Catfish_2022",]$tsurf.max.prev <- 
  mean(alldata[alldata$lake_id == "Catfish" & alldata$jday_cum > 429 & alldata$jday_cum < 454 & is.na(alldata$jday)==F,]$t.surf, na.rm=T)

# same for big sister 2021:
outbreak.summary[outbreak.summary$lake_yr=="Big Sister_2021",]$tmin.max.prev <- 
  mean(alldata[alldata$lake_id == "Big Sister" & alldata$jday_cum < 110 & is.na(alldata$jday)==F,]$t.min, na.rm=T)
outbreak.summary[outbreak.summary$lake_yr=="Big Sister_2021",]$tsurf.max.prev <- 
  mean(alldata[alldata$lake_id == "Big Sister" & alldata$jday_cum < 110 & is.na(alldata$jday)==F,]$t.surf, na.rm=T)

outbreak.summary



##############################
# things i've learned so far:

#outbreaks happened earlier in 2022?
hist(outbreak.summary$day.max.prev)
hist(log(outbreak.summary$day.max.prev))
summary(lme(log(day.max.prev)~year,
            random=~1|lake_id, data=outbreak.summary,
            na.action=na.omit))
# yes, marginally


# was it warmer in spring in 2022? 
hist(outbreak.summary$degree.days2)
hist(log(outbreak.summary$degree.days2))
summary(lme(log(degree.days2)~year,
            random=~1|lake_id, data=outbreak.summary,
            na.action=na.omit))
# yes. so, it was warmer in 2022, and outbreaks happened sooner... but hard to feel 
# confident when i didn't start sampling until later in 2021 anyway. 


# interesting. outbreaks that start earlier tend to be bigger
# statistics are iffy with low sample size, but consistent between years
plot(outbreak.summary$day.max.prev, outbreak.summary$max.prev, 
     bg=ifelse(outbreak.summary$year==2021, "black", "red"),
     pch=ifelse(outbreak.summary$year==2021, 21,22)) 
plot(log(outbreak.summary$day.max.prev), log(outbreak.summary$max.prev), 
     bg=ifelse(outbreak.summary$year==2021, "black", "red"),
     pch=ifelse(outbreak.summary$year==2021, 21,22)) 
summary(lm(log(max.prev)~log(day.max.prev)+year, data=outbreak.summary))
summary(lm(log(max.prev)~log(day.max.prev), data=outbreak.summary[outbreak.summary$year==2021,]))
summary(lm(log(max.prev)~log(day.max.prev), data=outbreak.summary[outbreak.summary$year==2022,]))
summary(lme(log(max.prev)~log(day.max.prev) + year, random=~1|lake_id, data=outbreak.summary,
            na.action=na.omit))
# the dreaded p = 0.058! maybe wait for another year of data and then revisit this pattern??


# are outbreaks bigger when water is colder???
plot(outbreak.summary$tmin.max.prev, outbreak.summary$max.prev, 
     bg=ifelse(outbreak.summary$year==2021, "black", "red"),
     pch=ifelse(outbreak.summary$year==2021, 21,22)) 
plot(outbreak.summary$tmin.max.prev, outbreak.summary$max.prev, 
     bg=ifelse(outbreak.summary$year==2021, "black", "red"),
     pch=ifelse(outbreak.summary$year==2021, 21,22),
     ylim=c(0,.2)) 
# hmmm no clear pattern. looks unimodal if catfish excluded though. 


# outbreak size does not seem related to integrated degree days in spring.
plot(outbreak.summary$degree.days2, outbreak.summary$max.prev,
     bg=ifelse(outbreak.summary$year==2021, "black", "red"),
     pch=ifelse(outbreak.summary$year==2021, 21,22)) 


colnames(alldata)

