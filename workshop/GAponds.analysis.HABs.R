##################################################################################
##################################################################################
###############                 GA ponds HABS Analysis              ##############
##################################################################################
##################################################################################


###################################################################
####### 1) read master file and make any necessary edits  #########
###################################################################


library(plyr)
#library(devtools)

setwd("C:/Users/straussa/Documents/Research/Strauss Lab - UGA/GA zoop surveys/Whitehall data")
alldata <- read.csv("GAponds.master.summary.csv", na.strings=c("", "NA"))


###################################################################
####### 2) how correlated are different metrics of algae?  ########
###################################################################

# will need to deal with negative values of PC before attempting log transform. 
# add value to each that is difference between most negative and smallest positive:
pc.int.unf.add <- min(alldata[alldata$pc.int.unf>0,]$pc.int.unf, na.rm=T) 
pc.tot.add <- min(alldata[alldata$pc.tot>0,]$pc.tot, na.rm=T) 
pc.surf.add <- min(alldata[alldata$pc.surf>0,]$pc.surf, na.rm=T) 

par(mfrow=c(1,1),mar=c(2,2,2,0.5), oma=c(0,0,0,0))

# total from integrated sample is correlated pretty well with area-under-curve of profile for chlorophyll:
plot(log(alldata$chl.int.unf), log(alldata$chl.tot))
# and phycocyanins:
plot(log(alldata$pc.int.unf+pc.int.unf.add), log(alldata$pc.tot+pc.tot.add))

# integrated unfiltered sample is even better correlated with chl on surface..
plot(log(alldata$chl.int.unf), log(alldata$chl.surf))
# okay for pc
plot(log(alldata$pc.int.unf+pc.int.unf.add), log(alldata$pc.surf+pc.surf.add))

# and surface vs. integrated trapz:
plot(log(alldata$chl.tot), log(alldata$chl.surf))
plot(log(alldata$pc.tot+pc.tot.add), log(alldata$pc.surf+pc.surf.add))

# okay, so, in general, these metrics are at least reasonably well correlated within chl or pc.
# But is chl generally correlated with pc?
plot(log(alldata$pc.int.unf+pc.int.unf.add), log(alldata$chl.int.unf))
plot(log(alldata$pc.surf+pc.surf.add), log(alldata$chl.surf))
plot(log(alldata$pc.tot+pc.tot.add),log(alldata$chl.tot))
# hmm all are pretty good relationships. good to know. so, more cyanobacteria usually means more
# total alage too


###################################################################
############## 3) is inedible algae assay working?  ###############
###################################################################

# i.e., is fluorescence always lower after pouring through 30um mesh??
edined <- alldata[!is.na(alldata$chl.int.30),]
edined$chl.remain <- 1-(edined$chl.int.30 / edined$chl.int.unf)
hist(edined$chl.remain)
plot(edined$chl.int.unf, edined$chl.remain); abline(h=0, lty=2)
edined[edined$chl.remain < -0.1,1:5]
# so for chlorophyll, it usually works, with one big outlier. okay. 

edined$pc.remain <- 1-(edined$pc.int.30 / edined$pc.int.unf)
hist(edined$pc.remain)
plot(edined$pc.int.unf, edined$pc.remain); abline(h=0, lty=2)
edined[edined$pc.remain < -0.1,1:5]
# and similar for pc. okay, seems promising. note that the outlier is a different
# sample though. need to get more data before assessing. 



###################################################################
############## 4) plot seasonal changes in algae  #################
###################################################################

# for plotting
colnames(alldata)
min(alldata$jday_cum)
max(alldata$jday_cum)

# re-order so lakes with more data plotted first
sumlakes <- unique(alldata$lake_id)
sumlakes.ord <- sumlakes[rev(order(table(alldata$lake_id)))]
colnames(alldata)

setwd("C:/Users/straussa/Documents/Research/Strauss Lab - UGA/GA zoop surveys/Whitehall data")
png("whitehall.HABs.png", width = 12, height = 5, res = 600, units='in')

par(mfrow=c(2,6),mar=c(0.5,0.5,0,0), oma=c(2.5,2.5,0.5,2.5))

for(i in 1:length(sumlakes.ord)){
  print(sumlakes.ord[i])
  temp <- alldata[alldata$lake_id==sumlakes.ord[i],]
  plot(temp$jday_cum, temp$t.surf, cex=0,ylim=c(0,60), xlim=c(70,620),
       xaxt="n", yaxt="n", ann=F)
  points(temp$jday_cum, temp$chl.surf, type="o", col="green", pch=0)
  points(temp$jday_cum, temp$t.surf, type="o", col="black", pch=1)
  main=sumlakes.ord[i]
  if(i>6){axis(side=1, at=c(100,300,500))}
  if(i==1 | i==7){axis(side=2, at=c(0,20,40,60))}
  
  # separate scale for phycocyanins 
  par(new=T)
  plot(temp$jday_cum, temp$pc.surf, ylim=c(0,22), xlim=c(70,620),
       xaxt="n", yaxt="n", ann=F)
  points(temp$jday_cum, temp$pc.surf, type="o", col="red", pch=2)
  if(i==6 | i==12){axis(side=4, at=c(0,10,20), col="red", col.axis="red")}
}

dev.off()


# # just catfish:
# setwd("C:/Users/straussa/Documents/Research/Strauss Lab - UGA/GA zoop surveys/Whitehall data")
# png("catfish.png", width = 6, height = 4, res = 600, units='in')
# par(mfrow=c(1,1),mar=c(0.5,0.5,0,0), oma=c(2.5,2.5,0.5,2.5))
# temp <- alldata[alldata$lake_id=="Catfish",]
# plot(temp$jday_cum, temp$t.surf, cex=0,ylim=c(0,60), xlim=c(70,500),
#      xaxt="n", yaxt="n", ann=F)
# points(temp$jday_cum, temp$chl.surf, type="o", col="green", pch=0)
# points(temp$jday_cum, temp$t.surf, type="o", col="black", pch=1)
# axis(side=1, at=c(100,300,500))
# axis(side=2, at=c(0,20,40,60))
# par(new=T)
# plot(temp$jday_cum, temp$pc.surf, ylim=c(0,22), xlim=c(70,500),
#      xaxt="n", yaxt="n", ann=F)
# points(temp$jday_cum, temp$pc.surf, type="o", col="red", pch=2)
# axis(side=4, at=c(0,10,20), col="red", col.axis="red")
# dev.off()


# # just vip:
# setwd("C:/Users/straussa/Documents/Research/Strauss Lab - UGA/GA zoop surveys/Whitehall data")
# png("vip.png", width = 6, height = 4, res = 600, units='in')
# par(mfrow=c(1,1),mar=c(0.5,0.5,0,0), oma=c(2.5,2.5,0.5,2.5))
# temp <- alldata[alldata$lake_id=="VIP",]
# plot(temp$jday_cum, temp$t.surf, cex=0,ylim=c(0,60), xlim=c(70,500),
#      xaxt="n", yaxt="n", ann=F)
# points(temp$jday_cum, temp$chl.surf, type="o", col="green", pch=0)
# points(temp$jday_cum, temp$t.surf, type="o", col="black", pch=1)
# axis(side=1, at=c(100,300,500))
# axis(side=2, at=c(0,20,40,60))
# par(new=T)
# plot(temp$jday_cum, temp$pc.surf, ylim=c(0,22), xlim=c(70,500),
#      xaxt="n", yaxt="n", ann=F)
# points(temp$jday_cum, temp$pc.surf, type="o", col="red", pch=2)
# axis(side=4, at=c(0,10,20), col="red", col.axis="red")
# dev.off()


###################################################################
########### 5) plot HABs project reservoirs on one figure #########
###################################################################

bear <- alldata[alldata$lake_id=="Bear Creek",]
chap <- alldata[alldata$lake_id=="Chapman",]
ogle <- alldata[alldata$lake_id=="Oglethorpe",]
her <- alldata[alldata$lake_id=="Herrick",]
mem <- alldata[alldata$lake_id=="Memorial",]

par(mfrow=c(1,1),mar=c(3,3,3,1), oma=c(0,0,0,0))

plot(bear$jday, bear$pc.int.unf, typ="o", col="blue",
     xlim=c(90,260), ylim=c(0,7),
     xlab="Julian day", ylab="Phycocyanin integrated depth",
     main="Algae in reservoirs summer 2022")
points(chap$jday, chap$pc.int.unf, typ="o", col="green")
points(ogle$jday, ogle$pc.int.unf, typ="o", col="red")
points(mem$jday, mem$pc.int.unf, typ="o", col="brown")
points(her$jday, her$pc.int.unf, typ="o", col="orange")



###################################################################
########### 6) extract some general summary statistics  ###########
###################################################################


max(alldata[alldata$lake_id=="VIP",]$max.depth, na.rm=T)
max(alldata[alldata$lake_id=="Deans",]$max.depth, na.rm=T)
max(alldata[alldata$lake_id=="Big Sister",]$max.depth, na.rm=T)

min(alldata[alldata$lake_id=="Catfish",]$max.depth, na.rm=T)
min(alldata[alldata$lake_id=="NN3",]$max.depth, na.rm=T)
min(alldata[alldata$lake_id=="Deer",]$max.depth, na.rm=T)
min(alldata[alldata$lake_id=="Sister 1",]$max.depth, na.rm=T)



###################################################################
####### 7) plot some dominant zooplankton taxa over time  #########
###################################################################

cat <- alldata[alldata$lake_id=="Catfish",]
sis1 <- alldata[alldata$lake_id=="Sister 1",]
bs <- alldata[alldata$lake_id=="Big Sister",]
vip <- alldata[alldata$lake_id=="VIP",]

colnames(cat)

plot(cat$jday_cum, cat$laevis, type="o", xlim=c(80,500))
plot(sis1$jday_cum, sis1$laevis, type="o", xlim=c(80,500))
plot(bs$jday_cum, bs$ambig, type="o", xlim=c(80,500))


