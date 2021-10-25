# Accessory functions


Summarize.mega.beta<-function(results_object){
  summary_object<-matrix(ncol=3,nrow=length(results_object))
  summary_object[,1]<-names(results_object)
  i<-1
  for(i in 1:length(results_object)){
    temp<-results_object[[i]]
    temp<-na.omit(temp)
    summary_object[i,2]<-mean(temp)
    summary_object[i,3]<-sd(temp)
  }
  return(summary_object)
}




nsite_list<-function(sites_ages,age_bins,PAtable){
  if(missing(age_bins)){
    ages<-unique(sites_ages$timeybp)
    ages<-ages[order(ages)]
    nsites<-matrix(nrow=(length(ages)),ncol=3)
    nsites[,1]<-ages[1:(length(ages))]
    colnames(nsites)<-c("Age","Nsites","Nspecies")
    # ages<-unique(sites_ages$timeybp)
    # ages<-ages[order(ages)]
    for(i in 1:length(ages)){
      age<-ages[i] 
      sites_age_new<-sites_ages[sites_ages$timeybp==age,]
      #PAtable_temp<-PAtable[,colSums(PAtable)!=0]
      PAtable_temp2<-PAtable[,rownames(sites_age_new)]
      PAtable_temp2<-PAtable_temp2[rowSums(PAtable_temp2)>0,]
      PAtable_temp2<-t(PAtable_temp2)
      nsites[i,2]<-nrow(PAtable_temp2)
      nsites[i,3]<-ncol(PAtable_temp2)
    }
  }else{
    ages<-age_bins[order(age_bins)]
    nsites<-matrix(nrow=(length(ages)-1),ncol=3)
    nsites[,1]<-ages[2:(length(ages))]
    colnames(nsites)<-c("Age","Nsites","Nspecies")
    for(i in 1:(length(ages)-1)){
      upper<-ages[i+1] 
      lower<-ages[i]
      sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
      sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
      if(nrow(sites_age_new)>1){
        PAtable_temp<-PAtable[,rownames(sites_age_new)]
        PAtable_temp<-PAtable_temp[rowSums(PAtable_temp)>0,]
        PAtable_temp<-t(PAtable_temp)
        nsites[i,2]<-nrow(PAtable_temp)
        nsites[i,3]<-ncol(PAtable_temp)
      }
    }
  }
  return(nsites)
}



non_binned<-function(sites_ages,PAtable,method=method){
  results<-list()
  ages<-unique(sites_ages$timeybp)
  ages<-ages[order(ages)]
  for(i in 1:length(ages)){
    print(i)
    age<-ages[i] 
    sites_age_new<-sites_ages[sites_ages$timeybp==age,]
    if(nrow(sites_age_new) < 1){
      results[[i]]<-"NA"
      next
    }
    PAtable_temp<-PAtable[,rownames(sites_age_new)]
    PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
    PAtable_temp<-t(PAtable_temp)
    if(method=="NA"){
      results[[i]]<-t(PAtable_temp)
    }
    if(method=="Sorensen"){
      results[[i]]<-vegdist(PAtable_temp,"bray")
    }
    if(method=="Jaccard"){
      results[[i]]<-vegdist(PAtable_temp,"jaccard",binary=TRUE)
    }
    if(method=="Forbes"){
      source("http://bio.mq.edu.au/~jalroy/Forbes.R")
      results[[i]]<-as.dist(1-forbesMatrix(PAtable_temp,corrected=TRUE))
    }
    if(method=="Betapair"){
      library(betapart)
      results[[i]]<-beta.pair(PAtable_temp)
    }
  }
  if(method=="NA"){
    results_new<-results
  }
  if(method=="Betapair"){
    names(results)<-ages
    results_new<-Summarize.beta.part(results,method="Betapair")
  }
  if(method!="Betapair"& method!="NA"){
    names(results)<-ages[2:(length(results)+1)]
    results_new<-Summarize.mega.beta(results)
  }
  return(results_new)
}






Binned<-function(age_bins,PAtable,sites_ages,method=method,calcmean=calcmean){
  ages<-age_bins[order(age_bins)]
  results<-list()
  for(i in 1:(length(ages)-1)){
    print(i)
    upper<-ages[i+1] 
    lower<-ages[i]
    sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
    sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
    if(nrow(sites_age_new)<2){
      results[[i]]<-"NA" #next
      next
    }
    if(nrow(sites_age_new)>2){
      PAtable_temp<-PAtable[,rownames(sites_age_new)]
      PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
      PAtable_temp<-t(PAtable_temp)
      if(method=="NA"){
        results[[i]]<-t(PAtable_temp)
      }
      if(method=="Sorensen"){
        results[[i]]<-vegdist(PAtable_temp,"bray")
      }
      if(method=="Jaccard"){
        results[[i]]<-vegdist(PAtable_temp,"jaccard",binary=TRUE)
      }
      if(method=="Forbes"){
        source("http://bio.mq.edu.au/~jalroy/Forbes.R")
        results[[i]]<-as.dist(1-forbesMatrix(PAtable_temp,corrected=TRUE))# I do not know if this is working correctly
      }
      if(method=="Betapair"){
        library(betapart)
        results[[i]]<-beta.pair(PAtable_temp)
      }
    }
  }
  if(method=="NA"){
    results_new<-results
    }
  if(method=="Betapair"){
    names(results)<-ages[2:(length(results)+1)]
    results_new<-Summarize.beta.part(results,method="Betapair")
  }
  if(method!="Betapair"& method!="NA" & calcmean=="TRUE"){
    names(results)<-ages[2:(length(results)+1)]
    results_new<-Summarize.mega.beta(results)
  }
  if(method!="Betapair"& method!="NA" & calcmean=="FALSE"){
    names(results)<-ages[2:(length(results)+1)]
    results_new<-results
  }
  return(results_new)
}



Summarize.betadisper<-function(betadisper_object){
  Beta_temp<-betadisper_object[[1]]
  level_beta<-unique(Beta_temp[,2])
  result<-matrix(nrow=length(level_beta),ncol=2)
  for(k in  1:length(level_beta)){
    Beta_temp2<-Beta_temp[Beta_temp[,2]==level_beta[k],]
    result[k,1]<-mean(Beta_temp2[,1])
    result[k,2]<-sd(Beta_temp2[,1])
  }
  return(result)
}




Calculate.ranges<-function(PAtable,sites_ages,age_bins){
  results_new<-Binned(age_bins,PAtable,method="NA")
  results<-list()
  for(i in 1:length(results_new)){
    PAtable_temp<-results_new[[i]]
    PAtable_temp<-PAtable_temp[rowSums(PAtable_temp)!=0,]
    spec<-rownames(PAtable_temp)
    results_temp<-matrix(nrow=length(spec),ncol=1)
    for(j in 1:length(spec)){
      PAspec<-results_new[[i]][spec[j],]
      results_temp[j,1]<-length(PAspec[PAspec>=1])
    }
    results[[i]]<-results_temp
  }
  return(results)
}



Calculate.extents<-function(PAtable,age_bins,sites_ages,latlongs,binned=c("TRUE","FALSE")){
  if(binned=="TRUE"){
    results_new<-Binned(age_bins,PAtable,sites_ages,method="NA")
  }
  if(binned=="FALSE"){
    results_new<-non_binned(sites_ages=sites_ages,PAtable=PAtable,method="NA")
  }
  results_temp<-matrix(nrow=length(results_new),ncol=10)
  for(j in 1:length(results_new)){
    PAtable_temp<-results_new[[j]]
    PAtable_temp<-PAtable_temp[rowSums(PAtable_temp)!=0,]
    PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)!=0]
    lat_longs_new<-latlongs[colnames(PAtable_temp),]
    results_temp[j,1]<-mean(lat_longs_new$latitude)
    results_temp[j,2]<-max(lat_longs_new$latitude)
    results_temp[j,3]<-min(lat_longs_new$latitude)
    results_temp[j,4]<-var(lat_longs_new$latitude)
    results_temp[j,5]<-mean(lat_longs_new$longitude)
    results_temp[j,6]<-max(lat_longs_new$longitude)
    results_temp[j,7]<-min(lat_longs_new$longitude)
    results_temp[j,8]<-var(lat_longs_new$longitude)
    results_temp[j,9]<-max(lat_longs_new$latitude)-min(lat_longs_new$latitude)
    results_temp[j,10]<-max(lat_longs_new$longitude)-min(lat_longs_new$longitude)
  }
  colnames(results_temp)<-c("Lat_mean","Lat_max","Lat_min","Lat_var",
                            "Long_mean","Long_max","Long_min","Long_var","Lat_extent","Long_extent")
  return(results_temp)
}





Summarize.beta.part<-function(results_object,method=method){
  summary_object<-matrix(nrow=length(results_object),ncol=7)
  summary_object[,1]<-names(results_object)
  if(method=="Betapair"){
    for(i in 1:length(results_object)){
      summary_object[i,2]<-mean(results_object[[i]][[1]])
      summary_object[i,3]<-sd(results_object[[i]][[1]])
      summary_object[i,4]<-mean(results_object[[i]][[2]])
      summary_object[i,5]<-sd(results_object[[i]][[2]])
      summary_object[i,6]<-mean(results_object[[i]][[3]])
      summary_object[i,7]<-sd(results_object[[i]][[3]])
    }
  }
  colnames(summary_object)<-c("Age","Beta.SIM","Beta.SIM.SD","Beta.SNE","Beta.SNE.SD","Beta.SOR","Beta.SOR.SD")
  return(summary_object)
}



Bin.climate.data<-function(age_bins,climate_data){
  ages<-age_bins[order(age_bins)]
  results<-matrix(nrow=(length(ages)-1), ncol=2)
  for(i in 1:(length(ages)-1)){
    upper<-ages[i+1] 
    lower<-ages[i]
    results[i,1]<-lower
    climate_new<-climate_data[climate_data$Time<=upper,]# find the sites with those ages
    climate_new<-climate_new[climate_new$Time>=lower,]
    if(nrow(climate_new)>1){
      results[i,2]<-mean(climate_new$d18O)
    }else{
      results[i,2]<-climate_new$d18O
    }
  }
  colnames(results)<-c("Age","d18O")
  return(results)
}


Climate.data<-function(sites_ages,climate_data){
  ages<-unique(sites_ages$timeybp)
  ages<-ages[order(ages)]
  results<-matrix(nrow=length(ages), ncol=2)
  for(i in 1:length(ages)){
    age<-ages[i]
    results[i,1]<-age
    if(age=="60"){
      climate_new<-climate_data[climate_data$Time==0,]# find the sites with those ages
      results[i,2]<-climate_new$d18O
    }else{
      climate_new<-climate_data[climate_data$Time==age,]# find the sites with those ages
      results[i,2]<-climate_new$d18O
    }
  }
  colnames(results)<-c("Age","d18O")
  return(results)
}




SubsampleBeta<-function(PAtable,sites_ages,age_bins,minsites,niters,latlongs,method=c("Sorensen","Jaccard","Forbes","Betadisper","DistDecay","Betapart")){
  library(vegan)
  library(fossil)
  library(glm2)
  min_sites<-minsites
  if(missing(age_bins)){
    ages<-unique(sites_ages$timeybp)
    ages<-ages[order(ages)]
    nsites<-nsite_list(sites_ages=sites_ages,PAtable=PAtable)
    results<-list()
    for(j in 1:niters){
      print(j)
      if(method=="Betadisper"){
        results[[j]]<-Betadisper.subsamp(PAtable=PAtable,sites_ages=sites_ages,latlongs=latlongs)
        next
      }
      if(method=="DistDecay"){
        results[[j]]<-Distdecay.subsamp(PAtable=PAtable,sites_ages=sites_ages,latlongs=latlongs)
        next
      }
      results_age<-list()
      for(i in 1:length(ages)){
        age<-ages[i] 
        sites_age_new<-sites_ages[sites_ages$timeybp==age,]
        PAtable_temp<-PAtable[,rownames(sites_age_new)]
        PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
        PAtable_temp<-t(PAtable_temp)
        PAtable_temp<-PAtable_temp[rowSums(PAtable_temp)>0,]
        if(nrow(PAtable_temp)<min_sites){
          next
        }
        sites_sampled<-sample(row.names(PAtable_temp),min_sites,replace=FALSE)
        PAtable_temp2<-PAtable_temp[sites_sampled,]
        PAtable_temp2<-PAtable_temp2[,colSums(PAtable_temp2)>0]
        if(method=="Sorensen"){
          results_age[[i]]<-vegdist(PAtable_temp2,"bray")
        }
        if(method=="Jaccard"){
          results_age[[i]]<-vegdist(PAtable_temp2,"jaccard",binary=TRUE)
        }
        if(method=="Forbes"){
          source("http://bio.mq.edu.au/~jalroy/Forbes.R")
          results_age[[i]]<-as.dist(1-forbesMatrix(PAtable_temp2,corrected=TRUE))
        }
        if(method=="Betapart"){
          library(betapart)
          results_age[[i]]<-beta.multi(PAtable_temp2)
        }
        if(method=="Betapart"){
          names(results_age)<-ages[1:length(results_age)]
          results[[j]]<-Summarize.beta.part(results_age)
        }
        if(method=="Betadisper"){
          results<-results
        }
        if(method!="Betapart" & method!="Betadisper"){
          names(results_age)<-ages[1:length(results_age)]
          results[[j]]<-Summarize.mega.beta(results_age)
        }
      }
    }
  }else{
    nsites_object<-nsite_list(sites_ages=sites_ages,age_bins=age_bins,PAtable=PAtable)
    ages<-age_bins[order(age_bins)]
    results<-list()
    for(j in 1:niters){
      print(j)
      if(method=="Betadisper"){
        results[[j]]<-Betadisper.subsamp(PAtable,sites_ages,age_bins,latlongs)
        next
      }
      if(method=="DistDecay"){
        results[[j]]<-Distdecay.subsamp(PAtable,sites_ages,age_bins,latlongs)
        next
      }
      results_age<-list()
      for(i in 1:(length(ages)-1)){
        upper<-ages[i+1] 
        lower<-ages[i]
        sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
        sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
        PAtable_temp<-PAtable[,rownames(sites_age_new)]
        PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
        PAtable_temp<-t(PAtable_temp)
        PAtable_temp<-PAtable_temp[rowSums(PAtable_temp)>0,]
        if(nrow(PAtable_temp)<min_sites){
          next
        }
        sites_sampled<-sample(row.names(PAtable_temp),min_sites,replace=FALSE)
        PAtable_temp2<-PAtable_temp[sites_sampled,]
        PAtable_temp2<-PAtable_temp2[,colSums(PAtable_temp2)>0]
        if(method=="NA"){
          results_age[[i]]<-t(PAtable_temp2)
        }
        if(method=="Sorensen"){
          results_age[[i]]<-vegdist(PAtable_temp2,"bray")
        }
        if(method=="Jaccard"){
          results_age[[i]]<-vegdist(PAtable_temp2,"jaccard",binary=TRUE)
        }
        if(method=="Forbes"){
          source("http://bio.mq.edu.au/~jalroy/Forbes.R")
          results_age[[i]]<-as.dist(1-forbesMatrix(PAtable_temp2,corrected=TRUE))
        }
        if(method=="Betapart"){
          library(betapart)
          results_age[[i]]<-beta.multi(PAtable_temp2)
        }
      }
      if(method=="Betapart"){
        names(results_age)<-ages[2:(length(results_age)+1)]
        results[[j]]<-Summarize.beta.part(results_age)
      }
      if(method=="Betadisper"){
        results<-results
      }
      if(method!="Betapart" & method!="Betadisper"){
        names(results_age)<-ages[2:(length(results_age)+1)]
        results[[j]]<-Summarize.mega.beta(results_age)
      }
    }
  }
  return(results)
}




Betadisper.subsamp<-function(PAtable,sites_ages,age_bins,latlongs){
  if(missing(age_bins)){
    ages<-unique(sites_ages$timeybp)
    ages<-ages[order(ages)]
    age<-ages[1]
    nsites<-nsite_list(sites_ages=sites_ages,PAtable=PAtable)
    sites_age_new<-sites_ages[sites_ages$timeybp==age,]
    PAtable_temp<-PAtable[,rownames(sites_age_new)]
    PAtable_temp<-t(PAtable_temp)
    PAtable_temp<-PAtable_temp[rowSums(PAtable_temp)>0,]
    sites_sampled<-sample(row.names(PAtable_temp),min_sites,replace=FALSE)
    PAtable_temp<-PAtable_temp[sites_sampled,]
    for(i in 2:length(ages)){
      age<-ages[i] 
      sites_age_new<-sites_ages[sites_ages$timeybp==age,]
      PAtable_temp2<-PAtable[,rownames(sites_age_new)]
      PAtable_temp2<-t(PAtable_temp2)
      if(nrow(PAtable_temp2)<min_sites){
        nsites<-nsites[-c(i), ]
        next
      }
      PAtable_temp2<-PAtable_temp2[rowSums(PAtable_temp2)>0,]
      sites_sampled<-sample(row.names(PAtable_temp2),min_sites,replace=FALSE)
      PAtable_temp2<-PAtable_temp2[sites_sampled,]
      PAtable_temp<-rbind(PAtable_temp,PAtable_temp2)
    }
  }else{
    nsites<-nsite_list(sites_ages=sites_ages,age_bins=age_bins,PAtable=PAtable)
    nsites<-na.omit(nsites)
    ages<-age_bins[order(age_bins)]
    upper<-ages[2]
    lower<-ages[1]
    sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
    sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
    PAtable_temp<-PAtable[,rownames(sites_age_new)]
    PAtable_temp<-t(PAtable_temp)
    sites_sampled<-sample(row.names(PAtable_temp),min_sites,replace=FALSE)
    PAtable_temp<-PAtable_temp[sites_sampled,]
    for(i in 2:(length(ages)-1)){
      upper<-ages[i+1]
      lower<-ages[i]
      sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
      sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
      if(nrow(sites_age_new)>1){
        PAtable_temp2<-PAtable[,rownames(sites_age_new)]
        PAtable_temp2<-t(PAtable_temp2)
        if(nrow(PAtable_temp2)<min_sites){
          nsites<-nsites[-c(i), ]
          next
        }
        sites_sampled<-sample(row.names(PAtable_temp2),min_sites,replace=FALSE)
        PAtable_temp2<-PAtable_temp2[sites_sampled,]
        PAtable_temp<-rbind(PAtable_temp,PAtable_temp2)
      }
    }
  }
  dist<-vegdist(PAtable_temp,method="jaccard",binary=TRUE)
  #nsite_temp<-nsites[1,2]# check if this works
  factors<-rep(1,15)
  for(i in 2:nrow(nsites)){
    #nsite_temp<-nsites[i,2]
    factor_grp<-rep(i,15)
    factors<-c(factors,factor_grp)
  }
  labels_temp<-as.factor(nsites[,1])
  groups<-factor(factors,labels=labels_temp)
  mod<-betadisper(dist,groups)
  results_new<-cbind(mod$distances,groups)
  results_new<-list(results_new,nsites)
  results_final<-Summarize.betadisper(results_new)
  results_final<-cbind(nsites[,1],results_final)
  return(results_final)
}




Distdecay.subsamp<-function(PAtable,sites_ages,age_bins,latlongs){
  library(fossil)
  library(glm2)
  if(missing(age_bins)){
    nsites_object<-nsite_list(sites_ages=sites_ages,PAtable=PAtable)
    results_new<-non_binned(sites_ages,PAtable,method="NA")
    nsites_object<-nsites_object[1:length(results_new),]
    for(i in 1:length(results_new)){
      if(ncol(results_new[[i]])<min_sites){
        results_new[[i]]<-results_new[[i]]
        next
      }
      sites_sampled<-sample(colnames(results_new[[i]]),min_sites,replace=FALSE)
      results_new[[i]]<-results_new[[i]][,sites_sampled]
    }
  }else{
    nsites_object<-nsite_list(sites_ages=sites_ages,age_bins=age_bins,PAtable=PAtable)
    nsites_object<-na.omit(nsites_object)
    results_new<-Binned(age_bins=age_bins,sites_ages=sites_ages,PAtable=PAtable,method="NA")
    nsites_object<-nsites_object[1:length(results_new),]
    for(i in 1:length(results_new)){
      if(ncol(results_new[[i]])<min_sites){
        next
      }
      sites_sampled<-sample(colnames(results_new[[i]]),min_sites,replace=FALSE)
      results_new[[i]]<-results_new[[i]][,sites_sampled]
    }
  }
  beta_matrix<-matrix(nrow=(length(results_new)),ncol=7)
  beta_matrix[,1]<-nsites_object[,1]
  colnames(beta_matrix)<-c("Ages","K","Median lat","Median long","n sites","n species","SE")
  for(j in 1:length(results_new)){
    PAtable_temp<-results_new[[j]]
    PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
    PAtable_temp<-PAtable_temp[rowSums(PAtable_temp)>0,]
    turnover<-1-vegdist(t(PAtable_temp),"jaccard",binary=TRUE)
    lat_longs_new<-latlongs[colnames(PAtable_temp),]
    lat<-lat_longs_new[,3]
    long<-lat_longs_new[,4]
    lat_longs_new<-cbind(long,lat)
    geo_dist<-earth.dist(lat_longs_new)
    coef.start <- coef(glm(turnover ~ geo_dist, family=binomial(link=log)))
    m <- glm2(turnover ~ geo_dist, family=binomial(link=log),start = coef.start)
    beta_matrix[j,2] <- -coef(m)[2]
    beta_matrix[j,3]<-median(lat_longs_new[,2])
    beta_matrix[j,4]<-median(lat_longs_new[,1])
    beta_matrix[j,5]<-ncol(PAtable_temp)
    beta_matrix[j,6]<-nrow(PAtable_temp[rowSums(PAtable_temp)>0,])
    beta_matrix[j,7]<-summary(m)[12]$coefficients[2,2]
  }
  results_new<-beta_matrix
  return(results_new)
}





metrics_only<-function(PAtable,method=method){
  library(vegan)
  if(method=="Sorensen"){
    results_dist<-vegdist(PAtable,"bray")
  }
  if(method=="Jaccard"){
    results_dist<-vegdist(PAtable,"jaccard",binary=TRUE)
  }
  if(method=="Forbes"){
    source("http://bio.mq.edu.au/~jalroy/Forbes.R")
    results_dist<-as.dist(1-forbesMatrix(PAtable,corrected=TRUE))
  }
  if(method=="Betapair"){
    library(betapart)
    results_dist<-beta.pair(PAtable)
  }
  return(results_dist)
}



Binned_range<-function(age_bins,PAtable,sites_ages){
  ages<-age_bins[order(age_bins)]
  specs<-nsite_list(sites_ages,age_bins,PAtable)
  results<-matrix(nrow=length(ages-1),ncol=3)
  num<-nrow(PAtable)
  nspecs<-c()
  for(i in 1:(length(ages)-1)){
    print(i)
    upper<-ages[i+1] 
    lower<-ages[i]
    sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]# this is where it's fucking up
    sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
    if(nrow(sites_age_new)<2){
      results[[i]]<-"NA"
      next
    }
    if(nrow(sites_age_new)>2){
      PAtable_temp<-PAtable[,rownames(sites_age_new)]
      PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
      PAtable_temp<-t(PAtable_temp)
      occups<-colSums(PAtable_temp)# gotta get rid of the zeroes
      reloccs<-occups[occups>0]
      nspecs[i]<-length(occups)
      #reloccs<-occups/nrow(PAtable_temp)
      results[i,1]<-mean(reloccs)
      results[i,2]<-sd(reloccs)
    }
  }
  results[,3]<-ages[2:(nrow(results)+1)]
  results<-na.omit(results)
  specs<-na.omit(specs)
  results_new<-list(results,specs,nspecs)
  return(results_new)
}




Non_binned_range<-function(PAtable,sites_ages){
  ages<-unique(sites_ages$timeybp)
  ages<-ages[order(ages)]
  specs<-nsite_list(sites_ages=sites_ages,PAtable=PAtable)
  results<-matrix(nrow=length(ages),ncol=3)
  for(i in 1:(length(ages))){
    print(i)
    age<-ages[i] 
    sites_age_new<-sites_ages[sites_ages$timeybp==age,]
    if(nrow(sites_age_new)<2){
      results[[i]]<-"NA"
      next
    }
    if(nrow(sites_age_new)>2){
      PAtable_temp<-PAtable[,rownames(sites_age_new)]
      PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
      PAtable_temp<-t(PAtable_temp)
      occups<-colSums(PAtable_temp)
      reloccs<-occups/nrow(PAtable_temp)
      results[i,1]<-mean(reloccs)
      results[i,2]<-sd(reloccs)
    }
  }
  results[,3]<-ages
  results<-na.omit(results)
  specs<-na.omit(specs)
  results_new<-list(results,specs)
  return(results_new)
}



Summarize.randomizations<-function(results_object){
  summary_object<-matrix(ncol=2,nrow=length(results_object))
  i<-1
  for(i in 1:length(results_object)){
    temp<-results_object[[i]]
    temp<-na.omit(temp)
    summary_object[i,1]<-mean(temp)
    summary_object[i,2]<-sd(temp)
  }
  return(summary_object)
}


Simulatechange<-function(sites_ages,age_bins,PAtable,niter,raw=c("TRUE","FALSE"),calcmeans=c("TRUE","FALSE")){
  iterations<-list()
  for(j in 1:niter){
    matrices<-list()
    ages<-age_bins[order(age_bins)]
    nages<-length(ages)-1
    specs<-nsite_list(sites_ages,age_bins,PAtable)
    sites<-colnames(PAtable)
    nsites_time<-specs[,2][1]
    if(is.na(nsites_time)=="TRUE"){
      next
    }
    selected_sites<-sample(sites,nsites_time,replace=FALSE)
    PAtable_temp<-PAtable[,selected_sites]
    #PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
    matrices[[1]]<-PAtable_temp
    #PAtable_temp[rowSums(PAtable_temp)>0,]
    diffs<-setdiff(colnames(PAtable),colnames(PAtable_temp))
    PAtable_new<-PAtable[,diffs]
    for(i in 2:nages){
      nsites_time<-specs[,2][i]
      if(is.na(nsites_time)=="TRUE"){
        next
      }
      selected_sites<-sample(diffs,nsites_time,replace=FALSE)
      PAtable_temp<-PAtable[,selected_sites]
      # PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
      matrices[[i]]<-PAtable_temp
      #PAtable_temp[rowSums(PAtable_temp)>0,]
      diffs<-setdiff(colnames(PAtable_new),colnames(PAtable_temp))
      PAtable_new<-PAtable_new[,diffs]
    }
    if(raw=="TRUE"){
      iterations[[j]]<-matrices
    }else{
      results<-list()
      for(k in 1:length(matrices)){
        results[[k]]<-metrics_only(t(matrices[[k]]),method="Jaccard")
      }
      if(calcmeans=="FALSE"){
        iterations[[j]]<-results
      }
      if(calcmeans=="TRUE"){
        iterations[[j]]<-cbind(Summarize.randomizations(results),na.omit(specs))
      }
    }
  }
  rm(j)
  #rm(k)
  rm(i)
  return(iterations)
}


Simulatechange_matrixswap<-function(sites_ages,age_bins,PAtable,niter,raw=c("TRUE","FALSE"),calcmeans=c("TRUE","FALSE")){
  iterations<-list()
  ages<-age_bins[order(age_bins)]
  nages<-length(ages)-1
  specs<-nsite_list(sites_ages,age_bins,PAtable)
  real_matrices<-list()
  for(i in 1:nages){
    upper<-ages[i+1] 
    lower<-ages[i]
    sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
    sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
    if(nrow(sites_age_new)>1){
      PAtable_temp<-PAtable[,rownames(sites_age_new)]
      real_matrices[[i]]<-PAtable_temp[rowSums(PAtable_temp)>0,]
    }
  }
  for(j in 1:niter){
    print(j)
    matrices<-sample(real_matrices,nages,replace=FALSE)
    if(raw=="TRUE"){
      iterations[[j]]<-matrices
    }else{
      results<-list()
      for(k in 1:length(matrices)){
        results[[k]]<-metrics_only(t(matrices[[k]]),method="Jaccard")
      }
      if(calcmeans=="FALSE"){
        iterations[[j]]<-results
      }
      if(calcmeans=="TRUE"){
        iterations[[j]]<-cbind(Summarize.randomizations(results),specs)
      }
    }
  }
  return(iterations)
}




summarize.sim<-function(results_object){
  ages<-results_object[[1]][,3]
  all_results<-matrix(nrow=length(ages),ncol=3)
  for(i in 1:length(ages)){
    results_temp<-matrix(nrow=length(results_object),ncol=1)
    for(j in 1:length(results_object)){
      results_temp[j,1]<-results_object[[j]][i,1]
    }
    temp_mean<-1-mean(results_temp[,1])
    all_results[i,1]<-temp_mean
    #all_results[i,2]<-sd(results_temp[,2])# do quantiles here instead.
    all_results[i,2]<-(1-quantile(results_temp[,1], c(0.025), type = 1))
    all_results[i,3]<-(1-quantile(results_temp[,1], c(0.975), type = 1))
  }
  return(all_results)
}


summarize.sim.sd<-function(results_object){
  ages<-results_object[[1]][,3]
  all_results<-matrix(nrow=length(ages),ncol=2)
  for(i in 1:length(ages)){
    results_temp<-matrix(nrow=length(results_object),ncol=1)
    for(j in 1:length(results_object)){
      results_temp[j,1]<-results_object[[j]][i,1]
    }
    temp_mean<-1-mean(results_temp[,1])
    all_results[i,1]<-temp_mean
    all_results[i,2]<-sd(results_temp[,1])# do quantiles here instead.
  }
  return(all_results)
}



summarize.resamp<-function(results_object){
  ages<-results_object[[1]][,1]
  all_results<-matrix(nrow=length(ages),ncol=3)
  for(i in 1:length(ages)){
    results_temp<-matrix(nrow=length(results_object),ncol=1)
    for(j in 1:length(results_object)){
      results_temp[j,1]<-results_object[[j]][i,2]
    }
    results_temp<-as.numeric(results_temp)
    temp_mean<-1-mean(results_temp)
    all_results[i,1]<-temp_mean
    #all_results[i,2]<-sd(results_temp[,2])# do quantiles here instead.
    all_results[i,2]<-(1-quantile(results_temp, c(0.025), type = 1))
    all_results[i,3]<-(1-quantile(results_temp, c(0.975), type = 1))
  }
  all_results<-cbind(ages,all_results)
  return(all_results)
}


summarize.daterands<-function(results_object){
  ages<-results_object[[1]][[2]][,1]
  all_results<-matrix(nrow=length(ages),ncol=3)
  for(i in 1:length(ages)){
    results_temp<-matrix(nrow=length(results_object),ncol=1)
    for(j in 1:length(results_object)){
      results_temp[j,1]<-results_object[[j]][[1]][i,2]
    }
    results_temp<-as.numeric(results_temp)
    temp_mean<-1-mean(results_temp)
    all_results[i,1]<-temp_mean
    #all_results[i,2]<-sd(results_temp[,2])# do quantiles here instead.
    all_results[i,2]<-(1-quantile(results_temp, c(0.025), type = 1))
    all_results[i,3]<-(1-quantile(results_temp, c(0.975), type = 1))
  }
  all_results<-cbind(ages,all_results)
  return(all_results)
}



Binned_rangeMCP<-function(PAtable,age_bins,sites_ages,latlongs,calcmeans=c("TRUE","FALSE")){
  library(sp)
  library(raster)
  library(adehabitatHR)
  #library(matrixStats)
  PAtable<-as.matrix(PAtable)
  site<-colnames(PAtable)
  num<-nrow(PAtable)
  ages<-age_bins[order(age_bins)]
  spec_results<-matrix(nrow=nrow(PAtable),ncol=length(ages)-1)
  rownames(spec_results)<-rownames(PAtable)
  colnames(spec_results)<-ages[2:length(ages)]
  specs<-nsite_list(sites_ages,age_bins=ages,PAtable)
  results<-matrix(nrow=length(ages)-1,ncol=2)
  rownames(results)<-ages[2:length(ages)]
  for(i in 1:num){# for each species
    print(i)
    df<-data.frame(x=PAtable[i,], sitekey=site)
    df.1<-subset(df, x==1)
    new<-merge(df.1, latlongs, by="sitekey")
    new.1<-merge(new, sites_ages, by="sitekey")
    for(k in 1:(length(ages)-1)){#for each age
      print(k)
      upper<-ages[k+1]
      lower<-ages[k]
      period<-subset(new.1, timeybp>lower)
      period<-subset(period, timeybp<=upper)
      if(nrow(period)<5){
        next
      }else{
        coordinates<-data.frame(x=period$longitude,y=period$latitude)
        coord<-SpatialPoints(coordinates, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
        coordinates <- spTransform(coord, CRS("+init=epsg:27573"))
        ch<-clusthr(coordinates,unin="m",unout="km2")
        ch2<-MCHu2hrsize(ch)
        spec_results[i,k]<-ch2$df[8]
      }
      
    }
  }
  if(calcmeans=="FALSE"){
    results_new<-spec_results
  }
  if(calcmeans=="TRUE"){
    nspecs<-apply(spec_results,2,function(x) num-sum(is.na(x)))
    results[,1]<-colMeans(spec_results,na.rm=TRUE)
    results[,2]<-apply(spec_results, 2, sd,na.rm=TRUE)
    results_new<-list(results,specs,nspecs)
  }
  return(results_new)
}



Null_rangeMCP<-function(PAtable,age_bins,sites_ages,latlongs,niter,calcmeans=c("TRUE","FALSE")){
  library(sp)
  library(raster)
  library(adehabitatHR)
  PAtable<-as.matrix(PAtable)
  num<-nrow(PAtable)
  ages<-age_bins[order(age_bins)]
  specs<-nsite_list(sites_ages,age_bins=ages,PAtable)
  results<-matrix(nrow=length(ages)-1,ncol=2)
  rownames(results)<-ages[2:length(ages)]
  all_results<-list()
  for(g in 1:niter){
    matrixdata<-Simulatechange(sites_ages,age_bins,PAtable,1,raw="TRUE")
    spec_results<-matrix(nrow=num,ncol=length(ages)-1)
    rownames(spec_results)<-rownames(PAtable)
    colnames(spec_results)<-ages[2:length(ages)]
    for(k in 1:(length(ages)-1)){
      print(k)
      matrixtemp<-matrixdata[[1]][[k]]
      site<-colnames(matrixtemp)
      specnames<-rownames(matrixtemp)
      for(i in 1:nrow(matrixtemp)){
        print(i)
        focalspec<-specnames[i]
        df<-data.frame(x=matrixtemp[i,], sitekey=site)
        df.1<-subset(df, x==1)
        new<-merge(df.1, latlongs, by="sitekey")
        period<-merge(new, sites_ages, by="sitekey")
        if(nrow(period)<3){
          next
        }
        else{
          coordinates<-data.frame(x=period$longitude,y=period$latitude)
          coord<-SpatialPoints(coordinates, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
          coordinates <- spTransform(coord, CRS("+init=epsg:27573"))
          #ch<-mcp(coordinates,unin="m",unout="km2")
          ch<-clusthr(coordinates,unin="m",unout="km2")
          ch2<-MCHu2hrsize(ch)
          spec_results[rownames(spec_results)==focalspec,k]<-ch2$df[8]
        }
      }
    }
    if(calcmeans=="TRUE"){
      nspecs<-apply(spec_results,2,function(x) num-sum(is.na(x)))
      results[,1]<-colMeans(spec_results,na.rm=TRUE)
      results[,2]<-apply(spec_results, 2, sd,na.rm=TRUE)
      all_results[[g]]<-list(results,specs,nspecs)
    }
    if(calcmeans=="FALSE"){
      all_results[[g]]<-spec_results
    }
  } 
  return(all_results)
}




Null_rangeMCPswap<-function(PAtable,age_bins,sites_ages,latlongs,niter,calcmeans=c("TRUE","FALSE")){
  library(sp)
  library(raster)
  library(adehabitatHR)
  PAtable<-as.matrix(PAtable)
  num<-nrow(PAtable)
  ages<-age_bins[order(age_bins)]
  specs<-nsite_list(sites_ages,age_bins=ages,PAtable)
  results<-matrix(nrow=length(ages)-1,ncol=2)
  rownames(results)<-ages[2:length(ages)]
  all_results<-list()
  for(g in 1:niter){
    matrixdata<-Simulatechange_matrixswap(sites_ages,age_bins,PAtable,1,raw=TRUE)
    spec_results<-matrix(nrow=num,ncol=length(ages)-1)
    rownames(spec_results)<-rownames(PAtable)
    colnames(spec_results)<-ages[2:length(ages)]
    for(k in 1:(length(ages)-1)){
      print(k)
      matrixtemp<-matrixdata[[1]][[k]]
      site<-colnames(matrixtemp)
      specnames<-rownames(matrixtemp)
      for(i in 1:nrow(matrixtemp)){
        print(i)
        focalspec<-specnames[i]
        df<-data.frame(x=matrixtemp[i,], sitekey=site)
        df.1<-subset(df, x==1)
        new<-merge(df.1, latlongs, by="sitekey")
        period<-merge(new, sites_ages, by="sitekey")
        if(nrow(period)<3){
          next
        }
        else{
          coordinates<-data.frame(x=period$longitude,y=period$latitude)
          coord<-SpatialPoints(coordinates, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
          coordinates <- spTransform(coord, CRS("+init=epsg:27573"))
          ch<-clusthr(coordinates,unin="m",unout="km2")
          ch2<-MCHu2hrsize(ch)
          spec_results[rownames(spec_results)==focalspec,k]<-ch2$df[8]
        }
      }
    }
    if(calcmeans=="FALSE"){
      all_results[[g]]<-spec_results
    }
    if(calcmeans=="TRUE"){
      nspecs<-apply(spec_results,2,function(x) num-sum(is.na(x)))
      results[,1]<-colMeans(spec_results,na.rm=TRUE)
      results[,2]<-apply(spec_results, 2, sd,na.rm=TRUE)
      all_results[[g]]<-list(results,specs,nspecs)
    }
  } 
  return(all_results)
}




summarize.sim.range<-function(results_object,calcmeans=c("TRUE","FALSE")){
  if(calcmeans=="TRUE"){
    ages<-results_object[[1]][[2]][,1]
    all_results<-matrix(nrow=length(ages),ncol=3)
    for(i in 1:length(ages)){
      results_temp<-matrix(nrow=length(results_object),ncol=1)
      for(j in 1:length(results_object)){
        results_temp[j,1]<-results_object[[j]][[1]][i,1]
      }
      all_results[i,1]<-mean(results_temp[,1],na.rm=TRUE)
      #all_results[i,2]<-sd(results_temp[,2])# do quantiles here instead.
      all_results[i,2]<-(quantile(results_temp[,1], c(0.025), type = 1,na.rm=TRUE))
      all_results[i,3]<-(quantile(results_temp[,1], c(0.975), type = 1,na.rm=TRUE))
    }
  }
  if(calcmeans=="FALSE"){
    first<-results_object[[1]]
    for(i in 2:length(results_object)){
      second<-results_object[[2]]
      first<-rbind(first,second)
    }
    all_results<-first
  }
  return(all_results)
}




Occupancy_grid<-function(PAtable,age_bins,sites_ages,latlongs,calcmeans=c("TRUE","FALSE")){
  library(sp)
  library(raster)
  library(adehabitatHR)
  library(maptools)
  library(sp)
  library(raster)
  #library(SDMTools)
  data(wrld_simpl)
  world <- crop(wrld_simpl, extent(-170, -15, -60, 90))
  behrmannCRS <- CRS("+proj=cea +lat_ts=30")
  world<-spTransform(world, behrmannCRS)
  r <- raster(world)
  res(r) <- 50000
  r[] <- rnorm(ncell(r))
  r[] <- 0
  world.raster <- rasterize(world, r)
  world.raster[world.raster>=1] <- 1
  world.raster[world.raster==0] <- NA
  world.raster<-aggregate(world.raster, fact=2, fun=max)
  Americas_grid<-rasterToPolygons(world.raster)
  Americas_grid<-as(Americas_grid, "SpatialPolygons")
  coords<-coordinates(Americas_grid) 
  PAtable<-as.matrix(PAtable)
  site<-colnames(PAtable)
  num<-nrow(PAtable)
  ages<-age_bins[order(age_bins)]
  spec_results<-matrix(nrow=nrow(PAtable),ncol=length(ages)-1)
  rownames(spec_results)<-rownames(PAtable)
  colnames(spec_results)<-ages[2:length(ages)]
  specs<-nsite_list(sites_ages,age_bins=ages,PAtable)
  results<-matrix(nrow=length(ages)-1,ncol=2)
  rownames(results)<-ages[2:length(ages)]
  for(i in 1:num){
    print(i)
    df<-data.frame(x=PAtable[i,], sitekey=site)
    df.1<-subset(df, x==1)
    new<-merge(df.1, latlongs, by="sitekey")
    new.1<-merge(new, sites_ages, by="sitekey")
    for(k in 1:(length(ages)-1)){
      PA_table1<-matrix(0,ncol=1,nrow=length(Americas_grid))
      rownames(PA_table1)<-1:nrow(PA_table1)
      print(k)
      upper<-ages[k+1]
      lower<-ages[k]
      period<-subset(new.1, timeybp>lower)
      period<-subset(period, timeybp<=upper)
      if(nrow(period)==0){
        print("No data")
        next
      }else{
        coordinates<-data.frame(x=period$longitude,y=period$latitude)
        coord<-SpatialPoints(coordinates, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
        coordinates_fossil <- spTransform(coord, behrmannCRS)
        for(j in 1:nrow(coordinates)){
          distVec <- spDistsN1(as.matrix(coords),coordinates_fossil[j],longlat = FALSE)
          minDistVec <- min(distVec)
          closestSiteVec <- which.min(distVec)
          PA_table1[closestSiteVec,1]<-1
        }
        spec_results[i,k]<-colSums(PA_table1)/length(coords)
      }
      
    }
  }
  if(calcmeans=="FALSE"){
    results_new<-spec_results
  }
  if(calcmeans=="TRUE"){ # edit this
    nspecs<-apply(spec_results,2,function(x) num-sum(is.na(x)))
    results[,1]<-colMeans(spec_results,na.rm=TRUE)
    results[,2]<-apply(spec_results, 2, sd,na.rm=TRUE)
    results_new<-list(results,specs,nspecs)
  }
  return(results_new)
}



Null_rangeGrid<-function(PAtable,age_bins,sites_ages,latlongs,niter,calcmeans=c("TRUE","FALSE"),method=c("Sites","Matrices")){
  library(sp)
  library(raster)
  library(adehabitatHR)
  library(maptools)
  library(sp)
  library(raster)
  #library(SDMTools)
  data(wrld_simpl)
  world <- crop(wrld_simpl, extent(-170, -15, -60, 90))
  behrmannCRS <- CRS("+proj=cea +lat_ts=30")
  world<-spTransform(world, behrmannCRS)
  r <- raster(world)
  res(r) <- 50000
  r[] <- rnorm(ncell(r))
  r[] <- 0
  world.raster <- rasterize(world, r)
  world.raster[world.raster>=1] <- 1
  world.raster[world.raster==0] <- NA
  world.raster<-aggregate(world.raster, fact=2, fun=max)
  Americas_grid<-rasterToPolygons(world.raster)
  Americas_grid<-as(Americas_grid, "SpatialPolygons")
  coords<-coordinates(Americas_grid) 
  PAtable<-as.matrix(PAtable)
  num<-nrow(PAtable)
  ages<-age_bins[order(age_bins)]
  specs<-nsite_list(sites_ages,age_bins=ages,PAtable)
  results<-matrix(nrow=length(ages)-1,ncol=3)
  rownames(results)<-ages[2:length(ages)]
  all_results<-list()
  for(g in 1:niter){
    print(g)
    if(method=="Sites"){
      matrixdata<-Simulatechange(sites_ages,age_bins,PAtable,1,raw="TRUE")
    }
    if(method=="Matrices"){
      matrixdata<-Simulatechange_matrixswap(sites_ages,age_bins,PAtable,1,raw=TRUE)
    }
    spec_results<-matrix(nrow=num,ncol=length(ages)-1)
    rownames(spec_results)<-rownames(PAtable)
    colnames(spec_results)<-ages[2:length(ages)]
    for(k in 1:length(matrixdata[[1]])){
      matrixtemp<-matrixdata[[1]][[k]]
      for(i in 1:nrow(matrixtemp)){
        PA_table1<-matrix(0,ncol=1,nrow=length(Americas_grid))
        rownames(PA_table1)<-1:nrow(PA_table1)
        site<-colnames(matrixtemp)
        df<-data.frame(x=matrixtemp[i,], sitekey=site)
        df.1<-subset(df, x==1)
        new<-merge(df.1, latlongs, by="sitekey")
        period<-merge(new, sites_ages, by="sitekey")
        if(nrow(period)==0){
          next
        }
        else{
          coordinates<-data.frame(x=period$longitude,y=period$latitude)
          coord<-SpatialPoints(coordinates, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
          coordinates_fossil <- spTransform(coord, behrmannCRS)
          for(j in 1:nrow(coordinates)){
            distVec <- spDistsN1(as.matrix(coords),coordinates_fossil[j],longlat = FALSE)
            minDistVec <- min(distVec)
            closestSiteVec <- which.min(distVec)
            PA_table1[closestSiteVec,1]<-1
          }
          spec_results[i,k]<-colSums(PA_table1)/length(coords)
        }
      }
    }
    if(calcmeans=="TRUE"){
      nspecs<-apply(spec_results,2,function(x) num-sum(is.na(x)))
      results[,1]<-colMeans(spec_results,na.rm=TRUE)
      for(i in 2:ncol(spec_results)){
        results[i-1,2]<-(quantile(spec_results[,i], c(0.025), type = 1,na.rm=TRUE))
        results[i-1,3]<-(quantile(spec_results[,i], c(0.975), type = 1,na.rm=TRUE))
      }
      all_results[[g]]<-list(results,specs,nspecs)
    }
    if(calcmeans=="FALSE"){
      all_results[[g]]<-spec_results
    }
  } 
  return(all_results)
}



summarize.sim.rangeGrid<-function(results_object,raw=c("TRUE","FALSE"),calcmeans=c("TRUE","FALSE")){
  if(calcmeans=="TRUE"& raw=="FALSE"){
    ages<-nrow(results_object[[1]][[1]])
    first<-results_object[[1]][[1]]
    for(i in 2:length(results_object)){
      second<-results_object[[2]][[1]]
      first<-cbind(first,second)
    }
    all_results<-matrix(nrow=ages,ncol=3)
    all_results[,1]<-rowMeans(first,na.rm=TRUE)
    for(i in 1:ages){
      all_results[i,2]<-(quantile(first[i,], c(0.025), type = 1,na.rm=TRUE))
      all_results[i,3]<-(quantile(first[i,], c(0.975), type = 1,na.rm=TRUE))
    }
  }
  if(calcmeans=="TRUE"& raw=="TRUE"){
    ages<-ncol(results_object)
    first<-na.omit(results_object[[1]])
    for(i in 2:length(results_object)){
      second<-results_object[[2]]
      first<-rbind(first,second)
    }
    all_results<-matrix(nrow=ages,ncol=3)
    all_results[,1]<-colMeans(first,na.rm=TRUE)
    for(i in 1:ages){
      all_results[i,2]<-(quantile(first[,i], c(0.025), type = 1,na.rm=TRUE))
      all_results[i,3]<-(quantile(first[,i], c(0.975), type = 1,na.rm=TRUE))
    }
  }
  if(calcmeans=="FALSE"){
    first<-na.omit(results_object[[1]])
    for(i in 2:length(results_object)){
      second<-results_object[[2]]
      first<-rbind(first,second)
    }
    all_results<-first
  }
  return(all_results)
}



linear_random_ages<-function(sites_ages,nsims){
  tradvsams<-read.csv("TradvsAMS.csv")
  model<-glm(AMS~Trad-1,data=tradvsams)
  results<-list()
  for(i in 1:nsims){
    sites_ages2<-sites_ages[order(sites_ages$timeybp),]
    colnames(sites_ages2)<-c("Site","Trad")
    nvar<-nrow(sites_ages2)
    new_ages<-predict.glm(model,sites_ages2)
    # add some random error here that is normally distributed
    for(b in 1:length(new_ages)){
      new_ages[b]<-new_ages[b]+rnorm(1,mean=0,sd=(summary(model)$coefficients[, 2])*sqrt(length(new_ages)))
    }
    sites_ages2$Trad<-new_ages
    colnames(sites_ages2)<-c("sitekey","timeybp")
    results[[i]]<-sites_ages2
  }
  return(results)
}



Random_ages<-function(sites_ages,nsims){
  results<-list()
  for(i in 1:nsims){
    sites_ages2<-sites_ages[order(sites_ages$timeybp),]
    modern<-sites_ages2[sites_ages2$timeybp<50,]
    diffs<-setdiff(rownames(sites_ages2),rownames(modern))
    sites_ages3<-sites_ages2[diffs,]
    nvar<-nrow(sites_ages3)
    for(b in 1:nvar){
      sites_ages3[b,2]<-sites_ages3[b,2]+abs(rnorm(1,mean=0,sd=2000))# makes some negative values
    }
    sites_ages3<-rbind(modern,sites_ages3)
    results[[i]]<-sites_ages3
  }
  return(results)
}



Climate_PCA_turnover<-function(sites_ages,age_bins,PCA_results){
  ages<-age_bins[order(age_bins)]
  results<-matrix(nrow=length(ages)-1,ncol=4)
  rownames(results)<-ages[1:(length(ages)-1)]
  for(i in 1:(length(ages)-2)){
    upper<-ages[i+1] 
    lower<-ages[i]
    sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
    sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
    if(nrow(sites_age_new)>1){
      PCA_results_temp<-PCA_results[rownames(sites_age_new),]
      temp12<-dist(PCA_results_temp[,1],method="euclidean")
      results[i,1]<-mean(temp12)
      results[i,2]<-sd(temp12)/sqrt(length(temp12))
      temp13<-dist(PCA_results_temp[,2],method="euclidean")
      results[i,3]<-mean(temp13)
      results[i,4]<-sd(temp13)/sqrt(length(temp13))
      
    }
  }
  return(results)
}

