Mega.beta.diversity<-function(PAtable,sites_ages,age_bins,latlongs,method=c("Sorensen","Jaccard","Forbes","Betadisper","DistDecay","Betapair","ABC"),calcmean=c("TRUE","FALSE")){
  library(vegan)
  if(method!="Betadisper" & method!="DistDecay"){
    if(missing(age_bins)){
      nsites_object<-nsite_list(sites_ages=sites_ages,PAtable=PAtable)
      results_new<-non_binned(sites_ages,PAtable,method)
    }else{
      nsites_object<-nsite_list(sites_ages=sites_ages,age_bins=age_bins,PAtable=PAtable)
      results_new<-Binned(age_bins,PAtable,sites_ages,method,calcmean=calcmean)
    }
  }
  if(method=="Betadisper"){
    if(missing(age_bins)){
      ages<-unique(sites_ages$timeybp)
      ages<-ages[order(ages)]
      age<-ages[1] 
      sites_age_new<-sites_ages[sites_ages$timeybp==age,]
      PAtable_temp<-PAtable[,rownames(sites_age_new)]
      PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)>0]
      PAtable_temp<-t(PAtable_temp)
      for(i in 2:length(ages)){
        age<-ages[i] 
        sites_age_new<-sites_ages[sites_ages$timeybp==age,]
        PAtable_temp2<-PAtable[,rownames(sites_age_new)]
        PAtable_temp2<-PAtable_temp2[,colSums(PAtable_temp2)>0]
        PAtable_temp2<-t(PAtable_temp2)
        PAtable_temp<-rbind(PAtable_temp,PAtable_temp2)
      }
      sites_ages_new<-sites_ages[rownames(PAtable_temp),]
      nsites<-nsite_list(sites_ages=sites_ages_new,PAtable=PAtable)
    }else{
      nsites<-nsite_list(sites_ages=sites_ages,age_bins=age_bins,PAtable=PAtable)
      nsites<-na.omit(nsites)
      ages<-age_bins[order(age_bins)]
      upper<-ages[2]
      lower<-ages[1]
      sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
      sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
      if(nrow(sites_age_new)>1){
        PAtable_temp<-PAtable[,rownames(sites_age_new)]
        PAtable_temp<-PAtable_temp[,colSums(PAtable_temp)!=0]
        PAtable_temp<-t(PAtable_temp)
      }
      for(i in 2:(length(ages)-1)){
        upper<-ages[i+1]
        lower<-ages[i]
        sites_age_new<-sites_ages[sites_ages$timeybp<=upper,]
        sites_age_new<-sites_age_new[sites_age_new$timeybp>lower,]
        if(nrow(sites_age_new)>1){
          PAtable_temp2<-PAtable[,rownames(sites_age_new)]
          PAtable_temp2<-PAtable_temp2[,colSums(PAtable_temp2)!=0]
          PAtable_temp2<-t(PAtable_temp2)
          PAtable_temp<-rbind(PAtable_temp,PAtable_temp2)
        }
      }
    }
    dist<-vegdist(PAtable_temp,method="jaccard",binary=TRUE)
    nsite_temp<-nsites[1,2]# check if this works
    factors<-rep(1,nsite_temp)
    for(i in 2:nrow(nsites)){
      nsite_temp<-nsites[i,2]
      factor_grp<-rep(i,nsite_temp)
      factors<-c(factors,factor_grp)
    }
    labels_temp<-as.factor(nsites[,1])
    groups<-factor(factors,labels=labels_temp)
    mod<-betadisper(dist,groups)
    results_new<-cbind(mod$distances,groups)
    nsites_object<-nsites
  }
  if(method=="DistDecay"){
    library(fossil)
    library(glm2)
    if(missing(age_bins)){
      nsites_object<-nsite_list(sites_ages=sites_ages,PAtable=PAtable)
      results_new<-non_binned(sites_ages,PAtable,method="NA")
    }else{
      nsites_object<-nsite_list(sites_ages=sites_ages,age_bins=age_bins,PAtable=PAtable)
      nsites_object<-na.omit(nsites_object)
      results_new<-Binned(age_bins=age_bins,sites_ages=sites_ages,PAtable=PAtable,method="NA",calcmean=calcmean)
      nsites_object<-nsites_object[1:length(results_new),]
    }
    beta_matrix<-matrix(nrow=(length(results_new)),ncol=7)
    beta_matrix[,1]<-nsites_object[,1]
    colnames(beta_matrix)<-c("Ages","K","Median lat","Median long","n sites","n species","SE")
    for(j in 1:length(results_new)){
      PAtable_temp<-results_new[[j]]
      turnover<-1-vegdist(t(PAtable_temp),"jaccard",binary=TRUE)
      lat_longs_new<-latlongs[colnames(PAtable_temp),]
      lat<-lat_longs_new[,2]
      long<-lat_longs_new[,3]
      lat_longs_new<-cbind(long,lat)
      geo_dist<-earth.dist(lat_longs_new)
      coef.start <- coef(glm(turnover ~ geo_dist, family=binomial(link=log)))
      m <- glm2(turnover ~ geo_dist, family=binomial(link=log),start = coef.start)
      if(inherits(glm2, "try-error")){ 
        stop("glm2 will often spit out an error related to the starting coefficient. Try new starting coefficients.")
      }
      beta_matrix[j,2] <- coef(m)[2]
      beta_matrix[j,3]<-median(lat_longs_new[,2])
      beta_matrix[j,4]<-median(lat_longs_new[,1])
      beta_matrix[j,5]<-ncol(PAtable_temp)
      beta_matrix[j,6]<-nrow(PAtable_temp[rowSums(PAtable_temp)>0,])
      beta_matrix[j,7]<-summary(m)[12]$coefficients[2,2]
    }
    results_new<-beta_matrix
  }
  return(list(results_new,nsites_object))
}
