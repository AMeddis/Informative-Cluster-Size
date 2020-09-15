#----Example for application on ICS
source("Zstat_fun.R")

#library(dplyr)

#---Example 1
library(pec)
data(Pbc3)
# PBC3 was a multi-centre randomized clinical trial conducted in six European hospitals.
# Between 1 Jan. 1983 and 1 Jan. 1987, 349 patients with the liver disease primary biliary cirrhosis 
# (PBC) were randomized to either treatment with Cyclosporin A (CyA, 176 patients) or placebo (173 patients). 
# The purpose of the trial was to study the effect of treatment on 
# "failure of medical treatment" defined as either death or liver transplantation

#event: failure of medical treatment
status<-ifelse(Pbc3$status==0,0,1)

#Test of ICS
Zst<-Zstat_opt(Pbc3$unit,status,Pbc3$days)
Zst

#---Example 2
library(MST)
data("Teeth")
##We consider data of patients treated at the Creighton University School of Dentistry from August 2007 to March 2013.
#A total of 5336 patients with periodontal disease were collected with 65228 teeth. 
#We excluded from the analysis individuals with only one tooth resulting in a sample size of 65034 teeth
#outcome of interest: time to tooth loss

#variable we need
dent<-select(Teeth,id,tooth,event,time)
dent<-group_by(dent,id) %>% mutate(ss=n())

dent<-filter(dent, ss>1)

#test of ICS
resSt<-Zstat_opt(dent$id,dent$event,dent$time)
resSt

