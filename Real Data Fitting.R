    rm(list=ls(all=TRUE))
    
    #We will obtain it in two steps
    #1. Obtain the observed chi-square value using ML estimates and different categories of data
    #2. For the same set of ML estimates, we calculate chi-square in each simulation (say; 1000) 
    #3. The p-value is the proportion of simulated chi-square values exceeding the observed value of chi-square obtained in step 1.   
    
    #===================================================================#
    #                                                                   #
    #             Importing and Filtering Data                          #
    #                                                                   #
    #===================================================================#  
    
    #loading libraries
    #library("VGAM");library("zipfR");library("readxl");library("coda");
    #library("MCMCpack");library("stats");library("stats4")
    
    #Importing and Filtering Data
    #Data_2<-read.csv("location",header=T)
    #View(Data_2)
    names(Data_2)
    # S1: monitoring time for recall category
    # z1: time for recall category      
    
    # S2: monitoring time for non-recall
    # S3: monitoring time for censored 
    
    # n1: sample size for recall
    # n2: sample size for non-recall   
    # n3: sample size for censored
    # n: sample size
    
    #summary(S1);summary(S2);summary(S3)     #summary of monitoring times
    lower1=0; 
    
    #===================================================================#
    #                                                                   #
    #             Calculate Observed Chi-Square                         #
    #                                                                   #
    #===================================================================#  
    
    # Obtain ML Estimates from real data using EM steps
    # a=MLE of alpha; b=MLE of beta; lmd=MLE of lambda 
    
    #Define Required Quantities for each category of recall, non-recall and right-censored
    
    #1. If observation is exact recall
    zz1=a*b*(z1^(a-1))*exp(-lmd*(S1-z1));zz1
    
    #2. If observation is non-recall
    dd1=exp(-b*(lower1^a))-exp(-b*(S2^a));dd1 # Denominator
    
    int1=function(u){a*b*(u^a)*exp(-b*(u^a))}
    intt1=array(0,length(S2))
    for (i in 1:length(S2)) {
      l=lower1
      u=S2[i]
      intt1[i]=integrate(int1, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt1
    xi1=intt1/dd1;xi1
    zz2=a*b*(xi1^(a-1))*exp(-b*(xi1^a))*(1-exp(-lmd*(S2-xi1)));zz2
    
    #3. If observation is right censored 
    dd2=exp(-b*(S3^a));dd2 
    int2=function(u){a*b*(u^a)*exp(-b*(u^a))}
    intt2=array(0,length(S3))
    for (i in 1:length(S3)) {
      l=S3[i]
      u=Inf
      intt2[i]=integrate(int2, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt2
    xi2=intt2/dd2;xi2
    zz3<-a*b*(xi2^(a-1))*exp(-b*(xi2^a));zz3
    
    # Find Modified Chi Square based on above quanities
    ch1<-(1-zz1)^2/(zz1*(1-zz1));ch1 #for recall category
    ch2<-(1-zz2)^2/(zz2*(1-zz2));ch2 #for non-recall category
    ch3<-(0-zz3)^2/(zz3*(1-zz3));ch3 #for right-censored category
    
    obs_chisq<-sum(ch1,ch2,ch3);obs_chisq  #observed chi-square 
    
    #===================================================================#
    #                                                                   #
    #             Calculate Simulated Chi-Square                        #
    #                                                                   #
    #===================================================================#  
    
    # n: Sample Size same as data
    # a=MLE of alpha; b=MLE of beta; lmd=MLE of lambda 
    
    ff=function(){
      #generate simulated data using ML parameters earlier explained in simulation 
      TT=rweibull(n,shape=b,scale=a^(-1/b));TT                         
      S=runif(n,min(TT),max(TT));S
      cc=ifelse(TT<S,1,0);cc            # Event Status
      
      tt=TT[which(cc!=0)];tt
      ss=S[which(cc!=0)];ss
      
      p1=exp(-(ss-tt)*lmd);mean(p1)      #recall-prob
      zz=rbinom(length(tt),1,p1);zz
      
      z1=tt[which(zz==1)];z1
      S1=ss[which(zz==1)];S1
      
      S2=ss[which(zz==0)];S2
      S3=S[which(cc==0)];S3
      
      n1=length(z1);n2=length(S2);n3=length(S3)
      n1;n2;n3;n=n1+n2+n3;n
      
      lower1=min(S2)/50         #set lower limit value as per data

      #---------------- Define Required Quantities-----------------------#
      # If observation is exact recall
      zz1=a*b*(z1^(a-1))*exp(-lmd*(S1-z1));zz1
      
      # If observation is non-recall
      int1=function(S,u){a*b*(u^(a-1))*(1-exp(-lmd*(S-u)))}
      zz2=sapply(1:length(S2),function(ii){integrate(function(u){int1(S2[ii],u)},lower1,S2[ii],subdivisions=100L,stop.on.error=FALSE)$value})  ;zz2
      
      # If observation is right censored 
      int2=function(u){exp(-b*(u^a))}
      intt2=array(0,length(S3))
      for (i in 1:length(S3)) {
        l=S3[i]
        u=Inf
        intt2[i]=integrate(int2, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      zz3<-intt2;zz3
      
      # Find Modified Chi Square
      ch1<-(1-zz1)^2/(zz1*(1-zz1));ch1
      ch2<-(1-zz2)^2/(zz2*(1-zz2));ch2
      ch3<-(0-zz3)^2/(zz3*(1-zz3));ch3
      
      ch_simu<-sum(ch1,ch2,ch3);ch_simu  #simulated chi-square 
      
      return(ch_simu)
    }
    ff()
    rr1=replicate(1000,ff());rr1
    pvalue<-length(which(rr1>obs_chisq))/1000;pvalue 
    
    