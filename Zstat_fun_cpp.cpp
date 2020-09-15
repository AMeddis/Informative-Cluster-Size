#include <Rcpp.h>
using namespace Rcpp;
#define SHOW(x) Rcpp::Rcout << #x <<"="<<x<< std::endl; 

/*
 * this function calculates N_AOM(t*) and Y_AOM(t*)
 * imput parameters
 * time:the vector of times of event 
 * tstar: time at which we want to calculate the sum (defines at-risk subjects)
 * status: the vector of status (going with time)
 */

// [[Rcpp::export]]
NumericVector sumA(NumericVector time, double tstar, NumericVector status) {
  double sumN=0;
  double sumY=0;
  int n=time.size();
  int t=0;
  
  while(time[t]<=tstar && t<n)
  { sumN+=status[t];
    t++;
  }
  for(;t<n;t++)
  { 
    sumY+=1;
  }
  NumericVector out = NumericVector::create(sumN,sumY);
  return out;
}

/*
 * this function calculates N_TOM(t*) and Y_TOM(t*)
 * imput parameters
 * time:the vector of times of event 
 * family: vector of cluster id
 * tstar: time at which we want to calculate the sum (defines at-risk subjects)
 * status: the vector of status (going with time)
 * ss: sample size for each cluster
 * 
 * We first sum in each cluster (sumN,sumY)
 * and after we sum over the cluster (sum2N,sum2Y) weighting with
 * the inverse of the cluster sample size
 */


// [[Rcpp::export]]
NumericVector sumT(NumericVector time, NumericVector family, double tstar, NumericVector status, NumericVector ss){
  
  int n_cluster=ss.size();
  NumericVector sumY(n_cluster);
  NumericVector sumN(n_cluster);

  int t=0;
  int n=time.size();
  while(time[t]<= tstar && t<n)
  {
    sumN[(family[t]-1)]+=status[t] ;
    t++;
  }
  for(int l=t; l<n;l++)
  { 
    sumY[(family[l]-1)]+=1;
  }
  

  double sum2N=0;
  double sum2Y=0;
  for(int k=0; k< n_cluster; k++)
  {
    sum2N+= sumN[k]/ss[k];
    sum2Y+= sumY[k]/ss[k];
  }
  
  
  NumericVector out = NumericVector::create(sum2N,sum2Y);
  //SHOW(t);
  return out;
}



/*
 * In this function we calculate sumA and sumT for t* in tempi
 * imput parameters
 * tempi :the vector of unique times of event 
 * time: the vector of observed times of events
 * family: vector of cluster id (according with time)
 * status: the vector of status (according with time)
 * ss: sample size
 * 
 * As output: N_AOM(t), Y_AOM(t), N_TOM(t), Y_TOM(t) with t in tempi
 * keep tracks on the repeated time of event with rep (for the variance calculation )
 */


// [[Rcpp::export]]
NumericMatrix sumall(NumericVector tempi, NumericVector time, NumericVector family, NumericVector status, NumericVector ss)
{
  int n=tempi.size();
  int n2=time.size();
  int l=0;
  NumericVector rep(n);
  NumericMatrix outA (n,2);
  NumericMatrix outT (n,2);
  NumericMatrix output(n,5);
  for(int t=0;t<n;t++)
  {
    outA(t,_)=sumA(time,tempi[t],status);
    outT(t,_)=sumT(time,family,tempi[t],status,ss);
    while(time[l]==tempi[t]  && l< n2)
    {
      rep[t]+=1;
      l++;
    }
  }
  output(_,0)=outA(_,0);
  output(_,1)=outA(_,1);
  output(_,2)=outT(_,0);
  output(_,3)=outT(_,1);
  output(_,4)=rep;
  return output;
}