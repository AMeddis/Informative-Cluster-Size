##---------------function for testing Informative Cluster Size------------------
#
#test for ICS with clustered survival data allowing for right cenoring
#employing the Nelson-Aalen estimator 

##source the functions in Rcpp to calculate
#N_tom, Y_tom, N_aom and N_tom
Rcpp::sourceCpp('Zstat_fun_cpp.cpp')

#@importFrom dplyr filter,group_by,mutate
#@importFrom data.table rleid

#@param family: vector indicating the cluster id
#@param status: vector indicating the status of event (1:event, 0: censored)
#@param time: vector of times of event

#@return object is an object that contains:
#Zstat: the test statistic standardized (follows Normal(0,1))
#test_stat: the test statistic with variance sd^2
#sd: standard deviation of the test statistic
#p.value: pvalue of the test (bilateral test)


Zstat_opt<-function(family,status,time)
{
  library(dplyr)
  library(data.table)
  
  #create the data.frame 
  sim<-data.frame(family=family,status=status,time=time)
  
  #create the index for the cluster
  sim<-sim[order(sim$family),]
  sim$ind<-rleid(sim$family)
  
  #number of clusters
  sim<-sim[order(sim$time),]
  n_cluster<-length(unique(sim$family))
  
  #sample size
  ni<-(group_by(sim,family) %>% dplyr::count())$n
  
  
  #vector of times
  nf<-which.max(sim$time)-1
  tempi<-unique(sim$time[1:nf])
  

  sumT<-sumall(tempi, sim$time,sim$ind, sim$status, ni)
  #AOM
  Na<-sumT[,1]
  Ya<-sumT[,2]
  dNa<-c(Na[1],diff(Na))
  #TOM
  Nt<-sumT[,3]
  Yt<-sumT[,4]
  dNt<-c(Nt[1],diff(Nt))
  
  #test statistic
  Z<-(1/sqrt(n_cluster))*sum((dNt/Yt - dNa/Ya)*Yt*Ya/n_cluster) 
  
  ##---------Variance estimation
  #YTOM and AOM for the observed failure times
  Ytobs<-rep(Yt,sumT[,5])
  Yaobs=rep(Ya,sumT[,5])
  
  sim2<-sim[1:nf,]
  ni2<-table(sim2$family)
  n_cluster2<-length(ni2)
  cid<-names(ni2)
  
  #sum_eps= \sum_{ijj'}\epsilon_ij \epsilon_ij'
  sum_eps=0  
  for(i in 1:n_cluster2)
  { 
    #for each cluster i
    ind<-which(sim2$family==cid[i])
    #omega_i=(Y_aom/(Ni*n_cluster) - Y_tom/N_i)
    Wi<-(Yaobs/(n_cluster*ni2[i])) - (Ytobs/n_cluster) 
    #num=Delta(time)\omega_i*Y(time_i)/Y(time)
    p2<-sapply(sim2$time[ind],FUN=function(t){sum(sim2$status *Wi*(sim2$time<=t)/Yaobs)})
    #\epsilon_ij
    epsilon<-(sim$status[ind]*Wi[ind])- p2 
    
  
    #\sum_{jj'}\epsilon_ij \epsilon_ij'
    prod<-sum(epsilon %*% t(epsilon))
    #\sum_i
    sum_eps<-sum_eps + prod
    
  }
  
  #variance
  V<-sum_eps/n_cluster
  
  #normalized test statistic
  Zstat=Z/sqrt(V)
  
 
  
  data.frame(Zstat=Zstat, test_stat=Z, sd=sqrt(V), p.value=2*(1-pnorm(abs(Zstat))) )
}
