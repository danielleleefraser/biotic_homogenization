
#### Steps to re-producing Figure 1 and Figure 4######################################################################################

#ALL SPECIES WITH ROCKY MOUNTAINS########################################################################################################
###########################################################################################################################################

PAtable<-read.csv("Complete PAtable_new.csv",header=T,row.names=1)
sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
matches<-intersect(rownames(PAtable),rownames(sites_ages))
sites_ages<-sites_ages[matches,]
sites_ages<-sites_ages[sites_ages$timeybp<50000,]
matches<-intersect(rownames(PAtable),rownames(sites_ages))
PAtable<-PAtable[matches,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
matches<-intersect(rownames(PAtable),rownames(latlongs))
latlongs<-latlongs[matches,]
latlongs$Machine.Number<-NULL

PAtable<-t(PAtable)

maxi<-max(sites_ages$timeybp)
mani<-min(sites_ages$timeybp)

# These are placeholders used for plotting but the date ranges are changed appropriately when plotting
age_bins<-seq(from=0,to=maxi,by=5000)
age_bins<-c(100,age_bins)

# Jaccard dissimilarity (converted downstread to similarity)
age_bins<-age_bins[1:8]
Jaccard_all<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Jaccard",calcmean="TRUE")
Jaccard_all[[1]][,3]<-as.numeric(Jaccard_all[[1]][,3])/sqrt(as.numeric(Jaccard_all[[2]][,2]))

# Geographic range size (these take some time to run)
Range_all<-Binned_rangeMCP(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")

# Occupancy
All_occ<-Occupancy_grid(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")


######Null model###################################################################################################################
#################################################################################################################################

# Null jaccard similarity
test<-Simulatechange(sites_ages,age_bins,PAtable,niter=1000,raw="FALSE",calcmeans="TRUE")
test_1<-summarize.sim(test)
test_1SD<-summarize.sim.sd(test)
Null_SD<-cbind(test[[1]][,3],test_1SD)
Null_final<-cbind(Jaccard_all[[1]][,3],test2)

# Null range size
simulated<-Null_rangeMCP(PAtable,age_bins,sites_ages,latlongs,niter=100,calcmeans="TRUE")
simfinal<-summarize.sim.range(simulated,calcmeans="TRUE")

# Null occupancy
null_occ<-Null_rangeGrid(PAtable,age_bins,sites_ages,latlongs,100,calcmeans="TRUE",method="Sites")
sim_final<-summarize.sim.rangeGrid(null_occ,raw="FALSE",calcmeans="TRUE")


# All species without rocky mountains Lyons PAtable no rockies#########################################################################################
##############################################################################################################################################

latlongs<-latlongs[latlongs$longitude>-109,]
matches<-intersect(rownames(latlongs),colnames(PAtable))
PAtable<-PAtable[,matches]
sites_ages<-sites_ages[matches,]

# Jaccard dissimilarity
age_bins<-age_bins[1:8]
Jaccard_nomt<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=Lyons_latlong,method="Jaccard",calcmean="TRUE")
Jaccard_nomt[[1]][,3]<-as.numeric(Jaccard_nomt[[1]][,3])/sqrt(as.numeric(Jaccard_nomt[[2]][,2]))

# Geographic range size
Range_nomt<-Binned_rangeMCP(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")

# Occupancy
Nomt_occ<-Occupancy_grid(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")

# Without extinct species########################################################################################################################
################################################################################################################################################

PAtable<-read.csv("Complete PAtable_new.csv",header=T,row.names=1)
extinct_mammals<-read.csv("Extinct_mammals.csv",header=T)
matches<-intersect(colnames(PAtable),extinct_mammals$Species)
diffs<-setdiff(colnames(PAtable),matches)
PAtable<-PAtable[,diffs]
PAtable<-PAtable[rowSums(PAtable)>0,]
PAtable<-PAtable[,colSums(PAtable)>0]

sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
matches<-intersect(rownames(PAtable),rownames(sites_ages))
sites_ages<-sites_ages[matches,]
sites_ages<-sites_ages[sites_ages$timeybp<50000,]
matches<-intersect(rownames(PAtable),rownames(sites_ages))
PAtable<-PAtable[matches,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
matches<-intersect(rownames(PAtable),rownames(latlongs))
latlongs<-latlongs[matches,]
latlongs$Machine.Number<-NULL

PAtable<-t(PAtable)

# Jaccard dissimilarity
Jaccard_noext<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=Lyons_latlong,method="Jaccard",calcmean=TRUE)
Jaccard_noext[[1]][,3]<-as.numeric(Jaccard_noext[[1]][,3])/sqrt(as.numeric(Jaccard_noext[[2]][,2]))

# Geographic range size
Range_noext<-Binned_rangeMCP(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")

# Occupancy
Noext_occ<-Occupancy_grid(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")

#### No occurrences North of the Canadian Border###########################################################################
#########################################################################################################################

PAtable<-read.csv("Complete PAtable_new.csv",header=T,row.names=1)
sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
matches<-intersect(rownames(PAtable),rownames(sites_ages))
sites_ages<-sites_ages[matches,]
sites_ages<-sites_ages[sites_ages$timeybp<50000,]
matches<-intersect(rownames(PAtable),rownames(sites_ages))
PAtable<-PAtable[matches,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
matches<-intersect(rownames(PAtable),rownames(latlongs))
latlongs<-latlongs[matches,]
latlongs$Machine.Number<-NULL

latlongs<-latlongs[latlongs$latitude<=50,]
matches<-intersect(rownames(latlongs),rownames(PAtable))
PAtable<-t(PAtable)
PAtable<-PAtable[,matches]
sites_ages<-sites_ages[matches,]

# Jaccard dissimilarity
Jaccard_noalaska<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Jaccard",calcmean="TRUE")
Jaccard_noalaska[[1]][,3]<-as.numeric(Jaccard_noalaska[[1]][,3])/sqrt(as.numeric(Jaccard_noalaska[[2]][,2]))

# Geographic range size
Range_noalaska<-Binned_rangeMCP(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")

# Occupancy
Noalaska_occ<-Occupancy_grid(PAtable,age_bins,sites_ages,latlongs,calcmeans="TRUE")

##### Removing species with small body masses#############################################################################
#########################################################################################################################

# Only species larger than 1 kg

PAtable_large<-read.csv("PAtable larger 1kg.csv",header=T,row.names=1)

sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
matches<-intersect(rownames(PAtable_large),rownames(sites_ages))
sites_ages<-sites_ages[matches,]
sites_ages<-sites_ages[sites_ages$timeybp<50000,]
matches<-intersect(rownames(PAtable_large),rownames(sites_ages))
PAtable_large<-PAtable_large[matches,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
matches<-intersect(rownames(PAtable_large),rownames(latlongs))
latlongs<-latlongs[matches,]
latlongs$Machine.Number<-NULL

PAtable_large<-t(PAtable_large)

maxi<-max(sites_ages$timeybp)
mani<-min(sites_ages$timeybp)

age_bins<-seq(from=0,to=maxi,by=5000)
age_bins<-c(100,age_bins)
age_bins<-age_bins[1:8]

# Jaccard dissimilarity
Jaccard_1<-Mega.beta.diversity(PAtable=PAtable_large,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Jaccard",calcmean="TRUE")
Jaccard_1[[1]][,3]<-as.numeric(Jaccard_1[[1]][,3])/sqrt(as.numeric(Jaccard_1[[2]][,2]))

# Geographic range size
Range_1<-Binned_rangeMCP(PAtable_large,age_bins,sites_ages,latlongs,calcmeans="TRUE")

# Occupancy
One_occ<-Occupancy_grid(PAtable_large,age_bins,sites_ages,latlongs,calcmeans="TRUE")


# Only species larger than 5 kg

PAtable_large<-read.csv("PAtable larger 5kg.csv",header=T,row.names=1)

sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
matches<-intersect(rownames(PAtable_large),rownames(sites_ages))
sites_ages<-sites_ages[matches,]
sites_ages<-sites_ages[sites_ages$timeybp<50000,]
matches<-intersect(rownames(PAtable_large),rownames(sites_ages))
PAtable_large<-PAtable_large[matches,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
matches<-intersect(rownames(PAtable_large),rownames(latlongs))
latlongs<-latlongs[matches,]
latlongs$Machine.Number<-NULL

PAtable_large<-t(PAtable_large)

maxi<-max(sites_ages$timeybp)
mani<-min(sites_ages$timeybp)

age_bins<-seq(from=0,to=maxi,by=5000)
age_bins<-c(100,age_bins)
age_bins<-age_bins[1:8]

# Jaccard dissimilarity
Jaccard_5<-Mega.beta.diversity(PAtable=PAtable_large,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Jaccard",calcmean="TRUE")
Jaccard_5[[1]][,3]<-as.numeric(Jaccard_5[[1]][,3])/sqrt(as.numeric(Jaccard_5[[2]][,2]))

# Geographic range size
Range_5<-Binned_rangeMCP(PAtable_large,age_bins,sites_ages,latlongs,calcmeans="TRUE")

#Occupancy
five_occ<-Occupancy_grid(PAtable_large,age_bins,sites_ages,latlongs,calcmeans="TRUE")


####### Create Figure 1#####################################################################################################
##########################################################################################################################


forPlot <- data.frame(Jaccard_all[[1]], Jaccard_noext[[1]], Jaccard_nomt[[1]],Jaccard_1[[1]],Jaccard_5[[1]],
  Jaccard_noalaska[[1]],Null_final)
df2 <- mutate_all(forPlot, function(x) as.numeric(as.character(x)))
colnames(df2)<-c("Age","Jaccardall","JaccardallSE",
  "Age2","Jaccardnoext","Jaccardnoext_SE",
  "Age3","Jaccardnomt","Jaccardnomt_SE",
  "Age4","Jaccard_large","Jaccard_large_SE",
  "Age6","Jaccard_5","Jaccard_5_SE",
  "Age7","Jaccard_alaska","Jaccard_a_SE",
  "Age5","Null","Lower","Upper")
df2$Jaccardall<-1-df2$Jaccardall
df2$Jaccardnoext<-1-df2$Jaccardnoext
df2$Jaccardnomt<-1-df2$Jaccardnomt
df2$Jaccard_large<-1-df2$Jaccard_large
df2$Jaccard_5<-1-df2$Jaccard_5
df2$Jaccard_alaska<-1-df2$Jaccard_alaska


library(ggplot2)
ggplot(df2, aes(x=Age7,y=Jaccard_alaska))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  scale_x_reverse(breaks=c(30000,25000,20000,15000,10000,5000,0),labels=c("30 ka-25 ka","25 ka-20 ka","20 ka-15 ka","15 ka-10 ka","10 ka - 5 ka","5 ka - modern","Modern"))+  
  ylim(0.0,0.4)+
  xlab("Time bins")+
  ylab("Mean proportion of shared taxa")+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="grey70")+
  geom_point(aes(x=Age,y=Null))+
  geom_line(aes(x=Age,y=Null))+
  geom_point(aes(x=Age,y=Jaccardnoext),colour="tan4")+
  geom_line(aes(x=Age,y=Jaccardnoext),colour="tan4",size=2)+
  geom_errorbar(aes(ymin=Jaccardnoext-Jaccardnoext_SE,ymax=Jaccardnoext+Jaccardnoext_SE),colour="tan4",width=0.2)+
  geom_point(aes(x=Age,y=Jaccardnomt),colour="tan3")+
  geom_line(aes(x=Age,y=Jaccardnomt),colour="tan3",size=2)+
  geom_errorbar(aes(ymin=Jaccardnomt-Jaccardnomt_SE,ymax=Jaccardnomt+Jaccardnomt_SE),colour="tan3",width=0.2)+
  geom_point(aes(x=Age,y=Jaccard_large),colour="tan2")+
  geom_line(aes(x=Age,y=Jaccard_large),colour="tan2",size=2)+
  geom_errorbar(aes(ymin=Jaccard_large-Jaccard_large_SE,ymax=Jaccard_large+Jaccard_large_SE),colour="tan2",width=0.2)+
  geom_point(aes(x=Age,y=Jaccard_5),colour="tan")+
  geom_line(aes(x=Age,y=Jaccard_5),colour="tan",size=2)+
  geom_errorbar(aes(ymin=Jaccard_5-Jaccard_5_SE,ymax=Jaccard_5+Jaccard_5_SE),colour="tan",width=0.2)+
  geom_point(aes(x=Age,y=Jaccardall),colour="black")+
  geom_line(aes(x=Age,y=Jaccardall),size=2,colour="black")+
  geom_errorbar(aes(ymin=Jaccardall-JaccardallSE,ymax=Jaccardall+JaccardallSE),colour="black",width=0.2)+
  geom_point(aes(x=Age7,y=Jaccard_alaska),colour="black")+
  geom_line(aes(x=Age7,y=Jaccard_alaska),size=2,colour="black")+
  geom_errorbar(aes(ymin=Jaccard_alaska-Jaccard_a_SE,ymax=Jaccard_alaska+Jaccard_a_SE),colour="black",width=0.2)




##### Reproduce Figure 2###################################################################################################
###########################################################################################################################
# Warning: takes a long time to run!

library(vegan)
library(sp)

Speciesmatrix<-read.csv("C:\\Users\\Dfraser\\Dropbox\\Peter Buck Postdoc\\NRI latitudinal gradient\\New matrix Behrmann projection.csv",header=TRUE,row.names=1)

locs<-read.csv("C:\\Users\\Dfraser\\Dropbox\\Peter Buck Postdoc\\NRI latitudinal gradient\\Coordinates behrmann projection.csv",header=T,row.names=1)

Speciesmatrix<-cbind(locs,Speciesmatrix)
Speciesmatrix<-Speciesmatrix[rowSums(Speciesmatrix[,3:ncol(Speciesmatrix)])>0,1:ncol(Speciesmatrix)]
Speciesmatrix<-Speciesmatrix[rowSums(Speciesmatrix[,3:ncol(Speciesmatrix)])>0,1:ncol(Speciesmatrix)]
Speciesmatrix_lats<-as.matrix(Speciesmatrix[,1:2])

nsamp<-100
results<-list()


for(i in 1:nsamp){
  results<-list()
  # do this multiple times and make a list
  Speciesmatrix_list<-matrix(nrow=nrow(Speciesmatrix),ncol=1)
  for(b in 1:nrow(Speciesmatrix)){
    print(b)
    Speciesmatrix_lats1<-Speciesmatrix_lats[b,]
    dists_Speciesmatrix<-spDistsN1(Speciesmatrix_lats,Speciesmatrix_lats1)
    names(dists_Speciesmatrix)<-rownames(Speciesmatrix_lats)
    probs<-1/dists_Speciesmatrix
    probs[probs==Inf]<-0
    samples_Speciesmatrix<-sample(dists_Speciesmatrix,99,prob=probs)
    Speciesmatrix_temp2<-Speciesmatrix[names(samples_Speciesmatrix),]
    Speciesmatrix_temp2<-rbind(Speciesmatrix[b,],Speciesmatrix_temp2)
    Speciesmatrix_temp2<-Speciesmatrix_temp2[,3:ncol(Speciesmatrix_temp2)]
    #Speciesmatrix_temp<-rbind(Speciesmatrix_temp,Speciesmatrix_temp2)
    factor_grp<-rep(1,100)
    groups<-factor(factor_grp,labels=1)    
    turnover<-vegdist(Speciesmatrix_temp2,method="jaccard",binary="TRUE")
    mod<-betadisper(turnover,groups)
    Speciesmatrix_list[b,1]<-mod$distances[1]
    
  }
  results[[i]]<-Speciesmatrix_list
}

Speciesmatrix_list2<-cbind(Speciesmatrix[,1:2],Speciesmatrix_list)
coordinates(Speciesmatrix_list2) = c("V1", "V2")
gridded(Speciesmatrix_list2) <- TRUE
x <- as(Speciesmatrix_list2, "SpatialGridDataFrame")
cols<-rev(heat.colors(n=200))
plot(x,col=cols)

plot(Speciesmatrix_list2[,3]~Speciesmatrix_list2[,2])

########################################################################################################################

###Reproduce Figure 3###################################################################################################
########################################################################################################################
#data from: http://www.datadryad.org/resource/doi:10.5061/dryad.1597g.2

require(rgdal)
library(tiff)
library(raster)
library(sp)
library(geoR) 
library(spsurvey)
library(shapefiles)

setwd("C:\\ccsm3_22-0k_all_tifs\\CCSM\\")
list.climate <- list.files(pattern =NULL, full.names=TRUE)
all.files <- list.files(path=list.climate, pattern =".tif$", full.names=TRUE)

sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
sites_ages<-sites_ages[sites_ages$timeybp<50000,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
latlongs$Machine.Number<-NULL

matches<-intersect(rownames(sites_ages),rownames(latlongs))
sites_ages<-sites_ages[matches,]
latlongs<-latlongs[matches,]
PAtable.Kate<-cbind(sites_ages,latlongs)
PAtable.Kate$sitekey<-NULL
PAtable.Kate$sitekey<-NULL
head(PAtable.Kate)
PAtable.Kate$time.bins <- cut(PAtable.Kate$timeybp, breaks=c(0,100,5000,10000,15000,20000,25000,30000), labels=c(0,5000,10000,15000,20000,25000,30000) )
head(PAtable.Kate)

ages<-seq(0,30000, by=5000)
ages<-c(100,ages)# I added this because you were combining the last time bin with 0-5000?
ages<-ages[order(ages)]
results<-list()
setwd("C:\\ccsm3_22-0k_all_tifs\\CCSM\\")
list.climate <- list.files(pattern =NULL, full.names=TRUE)
names<-substr(list.climate, 3, nchar(list.climate)-2)
name<-gsub("[^[:digit:]]", "", names)
name<-as.numeric(name)
climvars<-c("an_avg_TMAX.tif","an_avg_TMIN.tif","an_sum_AET.tif","an_sum_PRCP.tif")

overall_results<-list()
for(p in 1:length(climvars)){
  clim.temp<-all.files[grep(climvars[p], all.files)]
  climate.raster<-lapply(clim.temp, raster)
  names(climate.raster)<-paste(name)
  results<-list()
  for(i in 1:(length(ages)-1)){
    print(i)
    upper<-ages[i+1] 
    lower<-ages[i]
    sites_age_new<-PAtable.Kate[PAtable.Kate$timeybp<=upper,]
    sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
    coords_temp<-sites_age_new[,2:3]
    coordinates<-data.frame(x=coords_temp$longitude,y=coords_temp$latitude)
    coord<-SpatialPoints(coordinates, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    name_temp<-name[name<=upper]
    name_temp<-name_temp[name_temp>=lower]
    if(length(name_temp)<1){
      next
    }
    if(length(name_temp)==1){
      raster_temp2<-climate.raster[[as.character(name_temp)]]
      results[[i]]<-extract(raster_temp2, coord)
      results[[i]]<-cbind(coords_temp,results[[i]])
      colnames(results[[i]])<-c("latitude","longitude","results_temp3")
    }else{
      results_temp<-list()
      for(k in 1:length(name_temp)){
        raster_temp2<-climate.raster[[as.character(name_temp[k])]]
        results_temp[[k]]<-extract(raster_temp2, coord)
      }
      results_temp2<-results_temp[[1]]
      for(j in 2:length(results_temp)){
        results_temp2<-cbind(results_temp2,results_temp[[j]])
      }
      results_temp3<-rowMeans(results_temp2)
      results[[i]]<-cbind(coords_temp,results_temp3)
    }
  }
  names(results)<-ages[2:((length(ages))-1)]
  overall_results[[p]]<-results
}

names(overall_results)<-climvars

all_value_columns<-overall_results[[1]]
all_values<-do.call("rbind", all_value_columns)
rownames(all_values)<-c(rownames(all_value_columns[[1]]),rownames(all_value_columns[[2]]),rownames(all_value_columns[[3]]),
                        rownames(all_value_columns[[4]]),rownames(all_value_columns[[5]]),rownames(all_value_columns[[6]]))
colnames(all_values)<-c("latitude","longitude","an_avg_TMAX.tif")

names_all<-names(overall_results)

for(i in 2:length(overall_results)){
  all_value_columns<-overall_results[[i]]
  all_values_temp<-do.call("rbind", all_value_columns)
  rownames(all_values_temp)<-c(rownames(all_value_columns[[1]]),rownames(all_value_columns[[2]]),rownames(all_value_columns[[3]]),
                               rownames(all_value_columns[[4]]),rownames(all_value_columns[[5]]),rownames(all_value_columns[[6]]))
  colnames(all_values_temp)<-c("latitude","longitude",names_all[[i]])
  all_values<-cbind(all_values_temp,all_values)
}

PCA_results<-prcomp(~an_avg_TMAX.tif+an_avg_TMIN.tif+an_sum_AET.tif+an_sum_PRCP.tif,data=all_values,scale.=TRUE)
loadings<-as.data.frame(PCA_results$rotation)
biplot(PCA_results)
PCA_results<-PCA_results$x

all_results<-Climate_PCA_turnover(sites_ages,ages,PCA_results)
colnames(all_results)<-c("Mean_PC1","SE_PC1","Mean_PC2","SE_PC2")

library(dplyr)
forPlot <- data.frame(all_results)
df2<-cbind(c(0,5000,10000,15000,20000,25000,30000),forPlot)
colnames(df2)[1]<-c("Age")
df2 <- mutate_all(df2, function(x) as.numeric(as.character(x)))

library(ggplot2)
ggplot(df2, aes(x=Age,y=Mean_PC1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  scale_x_reverse(breaks=c(30000,25000,20000,15000,10000,5000,0),labels=c("30 ka-25 ka","25 ka-20 ka","20 ka-15 ka","15 ka-10 ka","10 ka - 5 ka","5 ka - modern","Modern"))+  
  ylim(0.0,3.0)+
  xlab("Time bins")+
  ylab("Mean climate difference")+
  geom_point(aes(x=Age,y=Mean_PC1),colour="black")+
  geom_line(aes(x=Age,y=Mean_PC1),size=2,colour="black")+
  geom_errorbar(aes(ymin=Mean_PC1-SE_PC1,ymax=Mean_PC1+SE_PC1),colour="black",width=0.2)+
  geom_point(aes(x=Age,y=Mean_PC2),colour="brown")+
  geom_line(aes(x=Age,y=Mean_PC2),size=2,colour="brown")+
  geom_errorbar(aes(ymin=Mean_PC2-SE_PC2,ymax=Mean_PC2+SE_PC2),colour="brown",width=0.2)

####### Create Figure 4a#####################################################################################################
##########################################################################################################################

library(ggplot2)
library(dplyr)

forPlot <- data.frame(Range_all[[2]][,1],Range_all[[2]][,3],Range_all[[1]],Range_noext[[1]],Range_nomt[[1]],Range_noalaska[[1]],Range_1[[1]],
                      Range_5[[1]],simfinal)
df2 <- mutate_all(forPlot, function(x) as.numeric(as.character(x)))
colnames(df2)<-c("Age","Nspeciesmeasured","Range","Range_SD","Range_noext","Range_SD_noext","Range_nomt","Range_SD_nomt",
                 "Range_noalaska","Range_noalaska_SD","Range_1","Range_1_SD","Range_5","Range_5_SD","Null_range","Lower","Upper")

ggplot(df2, aes(x=Age,y=Range))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  scale_x_reverse(breaks=c(30000,25000,20000,15000,10000,5000,0),labels=c("30 ka-25 ka","25 ka-20 ka","20 ka-15 ka","15 ka-10 ka","10 ka - 5 ka","5 ka - modern","Modern"))+  
  xlab("Time bins")+
  ylim(c(0,4e06))+
  ylab("Mean range size (km2)")+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="grey70")+
  geom_point(aes(x=Age,y=Null_range))+
  geom_line(aes(x=Age,y=Null_range))+
  geom_point(aes(x=Age,y=Range))+
  geom_line(aes(x=Age,y=Range,colour="black"))+
  geom_errorbar(aes(ymin=Range-(Range_SD/sqrt(Nspeciesmeasured)),ymax=Range+(Range_SD/sqrt(Nspeciesmeasured)),colour="black"),width=0.2)+
  geom_point(aes(x=Age,y=Range_nomt))+
  geom_line(aes(x=Age,y=Range_nomt,colour="green"))+
  geom_errorbar(aes(ymin=Range_nomt-(Range_SD_nomt/sqrt(Nspeciesmeasured)),ymax=Range_nomt+(Range_SD_nomt/sqrt(Nspeciesmeasured)),colour="green"),width=0.2)+
  geom_point(aes(x=Age,y=Range_noext))+
  geom_line(aes(x=Age,y=Range_noext,colour="red"))+
  geom_errorbar(aes(ymin=Range_noext-(Range_SD_noext/sqrt(Nspeciesmeasured)),ymax=Range_noext+(Range_SD_noext/sqrt(Nspeciesmeasured)),colour="red"),width=0.2)+
  geom_point(aes(x=Age,y=Range_noalaska))+
  geom_line(aes(x=Age,y=Range_noalaska,colour="red"))+
  geom_errorbar(aes(ymin=Range_noalaska-(Range_noalaska_SD/sqrt(Nspeciesmeasured)),ymax=Range_noalaska+(Range_noalaska/sqrt(Nspeciesmeasured)),colour="red"),width=0.2)+
  geom_point(aes(x=Age,y=Range_1))+
  geom_line(aes(x=Age,y=Range_1,colour="red"))+
  geom_errorbar(aes(ymin=Range_1-(Range_1_SD/sqrt(Nspeciesmeasured)),ymax=Range_1+(Range_1_SD/sqrt(Nspeciesmeasured)),colour="red"),width=0.2)+
  geom_point(aes(x=Age,y=Range_5))+
  geom_line(aes(x=Age,y=Range_5,colour="red"))+
  geom_errorbar(aes(ymin=Range_5-(Range_5_SD/sqrt(Nspeciesmeasured)),ymax=Range_5+(Range_5_SD/sqrt(Nspeciesmeasured)),colour="red"),width=0.2)


####### Create Figure 4b#####################################################################################################
#########################################################################################################################

All_occ[[1]][1:7,2]<-as.numeric(All_occ[[1]][1:7,2])/sqrt(as.numeric(All_occ[[2]][1:7,3]))
Nomt_occ[[1]][,2]<-as.numeric(Nomt_occ[[1]][,2])/sqrt(as.numeric(Nomt_occ[[2]][,3]))
Noext_occ[[1]][1:7,2]<-as.numeric(Noext_occ[[1]][1:7,2])/sqrt(as.numeric(Noext_occ[[2]][1:7,3]))
Noalaska_occ[[1]][1:7,2]<-as.numeric(Noalaska_occ[[1]][1:7,2])/sqrt(as.numeric(Noalaska_occ[[2]][1:7,3]))
One_occ[[1]][1:6,2]<-as.numeric(One_occ[[1]][1:6,2])/sqrt(as.numeric(One_occ[[2]][1:6,3]))
five_occ[[1]][1:6,2]<-as.numeric(five_occ[[1]][1:6,2])/sqrt(as.numeric(five_occ[[2]][1:6,3]))

library(dplyr)

forPlot <- data.frame(All_occ[[2]][,1],All_occ[[1]],Nomt_occ[[1]],Noext_occ[[1]],Noalaska_occ[[1]],One_occ[[1]],five_occ[[1]],sim_final)
df2<-cbind(forPlot,rownames(All_occ[[1]]))
df2 <- mutate_all(df2, function(x) as.numeric(as.character(x)))
colnames(df2)<-c("Age","All","Allse",
                 "Nomt","Nomtse",
                 "Noext","Noextse",
                 "NoAlaska","NoAlaskase",
                 "One","Onese",
                 "Five","Fivese",
                 "Null","Lower","Upper")


library(ggplot2)
ggplot(df2, aes(x=Age,y=All))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  scale_x_reverse(breaks=c(30000,25000,20000,15000,10000,5000,0),labels=c("30 ka-25 ka","25 ka-20 ka","20 ka-15 ka","15 ka-10 ka","10 ka - 5 ka","5 ka - modern","Modern"))+  
  ylim(0.0,0.008)+
  xlab("Time bins")+
  ylab("Mean geographic range size (number of grid cells)")+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="grey70")+
  geom_point(aes(x=Age,y=Null))+
  geom_line(aes(x=Age,y=Null))+
  geom_point(aes(x=Age,y=All),colour="black")+
  geom_line(aes(x=Age,y=All),colour="black",size=2)+
  geom_errorbar(aes(ymin=All-Allse,ymax=All+Allse),colour="black",width=0.2)+
  geom_point(aes(x=Age,y=Noext),colour="tan4")+
  geom_line(aes(x=Age,y=Noext),colour="tan4",size=2)+
  geom_errorbar(aes(ymin=Noext-Noextse,ymax=Noext+Noextse),colour="tan4",width=0.2)+
  geom_point(aes(x=Age,y=Nomt),colour="tan3")+
  geom_line(aes(x=Age,y=Nomt),colour="tan3",size=2)+
  geom_errorbar(aes(ymin=Nomt-Nomtse,ymax=Nomt+Nomtse),colour="tan3",width=0.2)+
  geom_point(aes(x=Age,y=One),colour="tan2")+
  geom_line(aes(x=Age,y=One),colour="tan2",size=2)+
  geom_errorbar(aes(ymin=One-Onese,ymax=One+Onese),colour="tan2",width=0.2)+
  geom_point(aes(x=Age,y=Five),colour="tan")+
  geom_line(aes(x=Age,y=Five),colour="tan",size=2)+
  geom_errorbar(aes(ymin=Five-Fivese,ymax=Five+Fivese),colour="tan",width=0.2)


#######Calculating values for Figure S2#####################################################################################################
#########################################################################################################################

PAtable<-read.csv("Complete PAtable_new.csv",header=T,row.names=1)
sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
matches<-intersect(rownames(PAtable),rownames(sites_ages))
sites_ages<-sites_ages[matches,]
sites_ages<-sites_ages[sites_ages$timeybp<50000,]
matches<-intersect(rownames(PAtable),rownames(sites_ages))
PAtable<-PAtable[matches,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
matches<-intersect(rownames(PAtable),rownames(latlongs))
latlongs<-latlongs[matches,]
latlongs$Machine.Number<-NULL

PAtable<-t(PAtable)

maxi<-max(sites_ages$timeybp)
mani<-min(sites_ages$timeybp)

# These are placeholders used for plotting but the date ranges are changed appropriately when plotting
age_bins<-seq(from=0,to=maxi,by=5000)
age_bins<-c(100,age_bins)

#### All values must be subtracted from 1 except for the distance decay of similarity metric

Betapar<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Betapair",calcmean = "TRUE")

Sor<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Sorensen",calcmean = "TRUE")

Forb<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Forbes",calcmean = "TRUE")

Dispersion<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="Betadisper",calcmean = "TRUE")

DDecay<-Mega.beta.diversity(PAtable=PAtable,sites_ages=sites_ages,age_bins=age_bins,latlongs=latlongs,method="DistDecay",calcmean = "TRUE")

#######Calculating values for Figure S3#####################################################################################################
#########################################################################################################################

subsamp_results<-SubsampleBeta(PAtable,sites_ages,age_bins,15,100,latlongs,method="Jaccard")
subsamp_final<-summarize.resamp(subsamp_results)

colnames(subsamp_final)<-c("Age","Mean","Upper","Lower")
subsamp_final<-data.frame(subsamp_final)
subsamp_final$Age<-as.numeric(as.character(subsamp_final$Age))
subsamp_final$Mean<-as.numeric(as.character(subsamp_final$Mean))
subsamp_final$Upper<-as.numeric(as.character(subsamp_final$Upper))
subsamp_final$Lower<-as.numeric(as.character(subsamp_final$Lower))

library(ggplot2)

ggplot(subsamp_final, aes(x=Age,y=Mean))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  scale_x_reverse(breaks=c(30000,25000,20000,15000,10000,5000,0),labels=c("30 ka-25 ka","25 ka-20 ka","20 ka-15 ka","15 ka-10 ka","10 ka - 5 ka","5 ka - modern","Modern"))+  
  ylim(0.0,0.35)+
  xlab("Time bins")+
  ylab("Mean proportion of shared taxa")+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="grey70")+
  geom_point(aes(x=Age,y=Mean))+
  geom_line(aes(x=Age,y=Mean))


#######Calculating values for Figure S4#####################################################################################################
#########################################################################################################################

rand_dates<-Random_ages(sites_ages,100)

Jaccard_results<-list()

for(i in 1:length(rand_dates)){
  print(i)
  Jaccard_results[[i]]<-Mega.beta.diversity(PAtable=PAtable,sites_ages=rand_dates[[i]],age_bins=age_bins,latlongs=latlongs,method="Jaccard",calcmean = "TRUE")
}

Dating_final<-summarize.daterands(Jaccard_results)

colnames(Dating_final)<-c("Age","Mean","Upper","Lower")
Dating_final<-data.frame(Dating_final)
Dating_final$Age<-as.numeric(as.character(Dating_final$Age))
Dating_final$Mean<-as.numeric(as.character(Dating_final$Mean))
Dating_final$Upper<-as.numeric(as.character(Dating_final$Upper))
Dating_final$Lower<-as.numeric(as.character(Dating_final$Lower))

library(ggplot2)

ggplot(Dating_final, aes(x=Age,y=Mean))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  scale_x_reverse(breaks=c(30000,25000,20000,15000,10000,5000,0),labels=c("30 ka-25 ka","25 ka-20 ka","20 ka-15 ka","15 ka-10 ka","10 ka - 5 ka","5 ka - modern","Modern"))+  
  ylim(0.0,0.35)+
  xlab("Time bins")+
  ylab("Mean proportion of shared taxa")+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="grey70")+
  geom_point(aes(x=Age,y=Mean))+
  geom_line(aes(x=Age,y=Mean))


#######Figure S5#####################################################################################################
#########################################################################################################################

# Panel A

plot(Jaccard_all[[2]][,3],(1-as.numeric(Jaccard_all[[1]][,2])),ylim=c(0,0.3))
model<-lm((1-as.numeric(Jaccard_all[[1]][,2]))~Jaccard_all[[2]][,3])
summary(model)

# Panel B

plot(Jaccard_all[[2]][,2],(1-as.numeric(Jaccard_all[[1]][,2])),ylim=c(0,0.3))
model2<-lm((1-as.numeric(Jaccard_all[[1]][,2]))~Jaccard_all[[2]][,2])
summary(model2)

#######Figure S6#####################################################################################################
#########################################################################################################################


mam<-read.csv("Complete PAtable_new.csv", header = TRUE, row.names = 1)
mam2<-t(mam)
mam.mat<-as.matrix(mam2)

age<-read.csv("absolutef_all.csv", header = TRUE)

latlong<-read.csv("Complete latlongs.csv", header = TRUE)

age$time.bins <- cut(age$timeybp, breaks=c(0,500,11600,20000,300000), labels=c("Modern", "Holocene", "Late_Glacial","Full_Glacial"))

time.bins<-unique(age$time.bins)
my.order<-c("Modern", "Holocene", "Full_Glacial", "Late_Glacial")

library(sp)
library(raster)

spec_results<-matrix(nrow=nrow(mam.mat),ncol=4) #to put results in
rownames(spec_results)<-rownames(mam.mat) 
colnames(spec_results)<-my.order

polygon.list<-list()
coord.list<-list()

num<-nrow(mam.mat)
site<-colnames(mam2)

for(i in 1:num){
  spec_poly<-list()
  coord.result<-list()
  print(i)
  df<-data.frame(x=mam.mat[i,], sitekey=site)
  df.1<-subset(df, x==1)
  new<-merge(df.1, latlong, by="sitekey")
  new.1<-merge(new, age, by="sitekey")
  for(k in 1:length(time.bins)){
    print(k)
    period<-subset(new.1, time.bins==levels(time.bins)[k])
    if(nrow(period)==0){
      spec_results[i,k]<-"NA"
      next 
    }else{
      coordinates<-data.frame(x=period$longitude,y=period$latitude)
      ch <- chull(x=coordinates$x, y=coordinates$y)
      coords <- coordinates[c(ch, ch[1]), ]
      coord<-SpatialPoints(coords, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
      sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coord)), ID=1)))
      spec_results[i,k]<-area(sp_poly) / 1000000
      spec_poly[[k]]<-sp_poly
      coord.result[[k]]<-coordinates
    }
  }
  polygon.list[[i]]<-spec_poly
  coord.list[[i]]<-coord.result
}

edited.results<-as.data.frame(spec_results)

edited.results$Modern[edited.results$Modern == "NA"]<-0
edited.results$Modern<-as.numeric(as.character(edited.results$Modern))
edited.results$Holocene[edited.results$Holocene == "NA"]<-0
edited.results$Holocene<-as.numeric(as.character(edited.results$Holocene))

edited.results$Late_Glacial[edited.results$Late_Glacial == "NA"]<-0
edited.results$Late_Glacial<-as.numeric(as.character(edited.results$Late_Glacial))

edited.results$Full_Glacial[edited.results$Full_Glacial == "NA"]<-0
edited.results$Full_Glacial<-as.numeric(as.character(edited.results$Full_Glacial))

edited.results$survivor<-ifelse(as.numeric(as.character(edited.results$Modern))+as.numeric(as.character(edited.results$Holocene))>0, paste("survivor"), paste("extinct"))
edited.results$survivor.avg<-ifelse(edited.results$survivor=="survivor", rowSums(edited.results[,1:4])/4, NA)
edited.results$extinct.avg<-ifelse(edited.results$survivor=="extinct", rowSums(edited.results[,3:4])/2, NA)

library(ggplot2)
ggplot(data=edited.results, aes(survivor, log10(as.numeric(as.character(Full_Glacial)))))+geom_boxplot()+theme_bw()+ylab("log range size")+xlab("")+ylim(c(2.5,7))
ggplot(data=edited.results, aes(survivor, log10(as.numeric(as.character(Late_Glacial)))))+geom_boxplot()+theme_bw()+ylab("log range size")+xlab("")+ylim(c(2.5,7))

fit <- aov(Full_Glacial~survivor, data=edited.results)
summary(fit)
fit2 <- aov(Late_Glacial~survivor, data=edited.results)
summary(fit2) # Check what this says in the text


#######Figure S7#####################################################################################################
#########################################################################################################################

PAtable<-read.csv("Complete PAtable_new.csv",header=T,row.names=1)

sites_ages<-read.csv("absolutef_all.csv",header=T)
rownames(sites_ages)<-sites_ages$sitekey
matches<-intersect(rownames(PAtable),rownames(sites_ages))
sites_ages<-sites_ages[matches,]
sites_ages<-sites_ages[sites_ages$timeybp<50000,]
matches<-intersect(rownames(PAtable),rownames(sites_ages))
PAtable<-PAtable[matches,]

latlongs<-read.csv("Complete latlongs.csv",header=T,row.names=1)
latlongs<-unique(latlongs)
rownames(latlongs)<-latlongs$sitekey
matches<-intersect(rownames(PAtable),rownames(latlongs))
latlongs<-latlongs[matches,]
latlongs$Machine.Number<-NULL

PAtable<-t(PAtable)

maxi<-max(sites_ages$timeybp)
mani<-min(sites_ages$timeybp)

age_bins<-seq(from=0,to=maxi,by=5000)
age_bins<-c(100,age_bins)

sites<-merge(latlongs, sites_ages, by="sitekey")
sites$time.bins <- cut(sites$timeybp, breaks=c(0,25,5000,10000,15000,20000,25000,32500), labels=c("0-25","500-5000","5000-10000","10000-15000","15000-20000","20000-25000","25000-30000"))

age.label<-unique(sites$time.bins)
num.age<-length(unique(sites$time.bins))

library(sp) 
library(raster)
library("mapdata")
library("maps")
library(ggplot2)

usa <- map_data("usa")
canada <- map_data("world", "Canada")
mexico <- map_data("world", "Mexico")

for(i in 1:num.age){
  print(i)
  new<-subset(sites, time.bins==age.label[i])
  if(nrow(new)==0){
    next
  }
  coordinates<-data.frame(x=new$longitude,y=new$latitude)
  coord<-SpatialPoints(coordinates, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  pdf(paste("Mammal sites ", unique(new$time.bins), " ybp.pdf", sep = ""))
  p<-ggplot() + 
    geom_polygon(data = usa, alpha=0.5, aes(x=long, y = lat, group = group))+
    geom_polygon(data = canada, alpha=0.5,aes(x=long, y = lat, group = group))+
    geom_polygon(data = mexico, alpha=0.5,aes(x=long, y = lat, group = group))+
    geom_point(data = new, aes(longitude, latitude), color="red", position = position_jitter(w = 0.3, h = 0.5)) +
    ggtitle(label=paste(unique(new$time.bins),"YBP"), subtitle = NULL)+
    coord_map(proj="gilbert")+theme_bw()
  print(p)
  dev.off()
}


##########################################################################################################################



