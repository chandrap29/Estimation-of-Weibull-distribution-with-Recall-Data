    
  rm(list=ls(all=TRUE))
  
    gg=function(l){
    #===================================================================#
    #                                                                   #
    #                        Sample Generation                          #
    #                                                                   #
    #===================================================================#  
      
    library("VGAM");library("coda") # loading necessary packages
    try({ #We applied try as the calculation for a few samples are absurd due to divergence in integrals   
      
      n=50;alp=alp1=1.35;bet=bet1=1.85;lmd=lmd1=0.20;alpha=0.05 # designing simulation parameters
      #Note: The code can be run for other set of parameters used in the article
      
      # n2=0
      # while(n2<n*0.10){                             #condition for preserving at least a minimum number in non-recall category
      TT=rweibull(n,shape=alp,scale=bet^(-1/alp));TT  #generate Time to event                         
      S=rexp(n,0.25);S                                #generate monitoring time from Exponential(rate=0.25)
      #S=runif(n,0,4);S                               #generate monitoring time from Uniform(a=0,b=4)
      cc=ifelse(TT<S,1,0);cc                          # Event Status; 1 = occurred, 0 = right censored
      
      tt=TT[which(cc!=0)];tt                #time to event for non-censored observations
      ss=S[which(cc!=0)];ss                 #monitoring time for non-censored observations
      
      p1=exp(-(ss-tt)*lmd)                  #defining recall probability 
      
      zz=rbinom(length(tt),1,p1);zz         #generate a binomial random variate to categorize the uncensored 
                                            #observations into recall and non-recall
      
      t1=tt[which(zz==1)];t1                #time-to-event for recalled subjects
      s1=ss[which(zz==1)];s1                #monitoring time for recalled subjects
      
      s2=ss[which(zz==0)];s2                #monitoring time for non-recalled subjects
      s3=S[which(cc==0)];s3                 #monitoring time for censored subjects
      
      z1=t1;S1=s1;S2=s2;S3=s3                   #re-define the variables 
      n1=length(z1);n2=length(S2);n3=length(S3) #check the length of each category
      #}
      z1;S1;S2;S3         
      #z1: time-to-event for recall,       #S1: monitoring time for recall
      #S2: monitoring time for non-recall, #S3: monitoring time for censored data
      n1;n2;n3  #check the length of each category
      
      
      #===================================================================#
      #                                                                   #
      #                        ML Estimate (Using EM)                     #
      #                                                                   #
      #===================================================================#  
      
      m.thh1=m.thh2=m.thh3=rep()          #defining vectors to store estimates of parameters at each E-M iteration  
      m.th1=1.18;m.th2=0.82;m.th3=0.12    #arbitrary initial values to start E-M algorithm 
      it=50;                              #set no. of iteration 
      
      for(j in 1:it){
        
        aa=m.th2*(S2^m.th1);aa
        bb=m.th2*(S3^m.th1);bb
        
        D1=1-exp(-aa);D1  #denominator used in expectation for non-recall 
        D2=exp(-bb);D2    ##denominator used in expectation for right-censored
        
        #--------------------------------------------------#
        #Define the Integrals Used in the Expectation Terms  
        #first integral 
        int1=function(x){log(x)*exp(-x)}
        integd1=array(0,length(aa))
        for (i in 1:length(aa)) {
          l=0
          u=aa[i]
          integd1[i]=integrate(int1, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-5)$value
        }
        i11=integd1;i11   
        
        #3rd integral
        int3=function(x){x*log(x)*exp(-x)}
        integd3=array(0,length(aa))
        for (i in 1:length(aa)) {
          l=0
          u=aa[i]
          integd3[i]=integrate(int3, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-5)$value
        }
        i33=integd3;i33    
        
        #4th integral
        int4=function(x){x*(log(x)^2)*exp(-x)}
        integd4=array(0,length(aa))
        for (i in 1:length(aa)) {
          l=0
          u=aa[i]
          integd4[i]=integrate(int4, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-5)$value
        }
        i44=integd4;i44    
        
        #second integral
        int2=function(x){(1/x)*exp(-x)}
        integd2=array(0,length(bb))
        for (i in 1:length(bb)) {
          l=bb[i]
          u=Inf
          integd2[i]=integrate(int2, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-5)$value
        }
        i22=log(bb)*exp(-bb)+integd2;i22    
        
        i55=bb*log(bb)*exp(-bb)+exp(-bb)+i22;i55   #5th integral 
        
        #6th integral
        int6=function(x){(1/x)*log(x)*exp(-x)}
        integd6=array(0,length(bb))
        for (i in 1:length(bb)) {
          l=bb[i]
          u=Inf
          integd6[i]=integrate(int6, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-5)$value
        }
        i66=bb*((log(bb))^2)*exp(-bb)+2*i22+((log(bb))^2)*exp(-bb)+2*integd6;i66   
        
        #-----------Define the Expectation terms---------------------#
        xi0=((m.th2^(-1/m.th1))*(gamma((1/m.th1)+1)*pgamma((1/m.th1)+1,aa,lower.tail = TRUE)))/D1;xi0
        xi1=(i11-(log(m.th2)*D1))/(m.th1*D1);xi1
        xi2=(1-(1+aa)*(1-D1))/(m.th2*D1);xi2
        xi3=(i22-(log(m.th2)*D2))/(m.th1*D2);xi3
        xi4=(1+bb)/m.th2;xi4
        cc1=m.th3*(S2-xi0);cc1
        xi5=(1-((1+cc1)*exp(-cc1)))/(m.th3*(1-exp(-cc1)));xi5
        xi6=(i33-log(m.th2)*(1-(1+m.th2*(S2^m.th1))*(1-D1)))/(m.th1*m.th2*D1)
        xi8=(i55-log(m.th2)*((1+m.th2*(S3^m.th1))*D2))/(m.th1*m.th2*D2)
        
        #M step of E-M 
        m.thh1[j]=m.th1;m.thh2[j]=m.th2;m.thh3[j]=m.th3 
        m.thh1[j+1]=n/((m.thh2[j]*sum((z1^m.thh1[j])*log(z1)))-sum(log(z1))-sum(xi1)-sum(xi3)+m.thh2[j]*(sum(xi6)+sum(xi8)))
        m.thh2[j+1]=n/(sum(z1^m.thh1[j+1])+sum(xi2)+sum(xi4))
        m.thh3[j+1]=n2/(sum(S1-z1)+sum(xi5))
        m.th1=m.thh1[j+1];m.th2=m.thh2[j+1];m.th3=m.thh3[j+1] #Storing the final value of estimates
      }
      m.thh1;m.thh2;m.thh3
      #par(mfrow=c(2,2))   #checking the convergence through plots
      #plot(m.thh1,type='l');plot(m.thh2,type='l');plot(m.thh3,type='l')
      ml=c(m.thh1[it],m.thh2[it],m.thh3[it]);ml  #store the final estimates
      mse=c((ml[1]-alp1)^2,(ml[2]-bet1)^2,(ml[3]-lmd1)^2);mse #calculate the Mean Square Error
      ab=c(abs(ml[1]-alp1),abs(ml[2]-bet1),abs(ml[3]-lmd1));ab #calculate the Absolute Bias
      
      
      #===================================================================#
      #                                                                   #
      #                 Asymptotic Confidence Intervals (ACIs)            #
      #                                                                   #
      #===================================================================# 
      #Calculating the Information Matrix using Missing Information Principle
      
      #Calculate the expressions as in EM steps on Maximum likelihood estimates (MLEs)   
      aa1=ml[2]*(S2^ml[1]);aa1    
      bb1=ml[2]*(S3^ml[1]);bb1
      
      D11=1-exp(-aa1);D11  #based on S2
      D22=exp(-bb1);D22    #based on S3
      
      #-----------------------------------------#
      intt1=function(x){log(x)*exp(-x)}
      #first integral
      integd11=array(0,length(aa1))
      for (i in 1:length(aa1)) {
        l=0
        u=aa1[i]
        integd11[i]=integrate(intt1, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      ii11=integd11;ii11    
      
      intt3=function(x){x*log(x)*exp(-x)}
      #3rd integral
      integd33=array(0,length(aa1))
      for (i in 1:length(aa1)) {
        l=0
        u=aa1[i]
        integd33[i]=integrate(intt3, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      ii33=integd33;ii33    
      
      intt4=function(x){x*(log(x)^2)*exp(-x)}
      #4th integral 
      integd44=array(0,length(aa1))
      for (i in 1:length(aa1)) {
        l=0
        u=aa1[i]
        integd44[i]=integrate(intt4, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      ii44=integd44;ii44   
      
      #second integral
      intt2=function(x){(1/x)*exp(-x)}
      integd22=array(0,length(bb1))
      for (i in 1:length(bb1)) {
        l=bb1[i]
        u=Inf
        integd22[i]=integrate(intt2, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-5)$value
      }
      ii22=log(bb1)*exp(-bb1)+integd22;ii22    
      
      ii55=bb1*log(bb1)*exp(-bb1)+exp(-bb1)+ii22;ii55   #5th integral 
      
      #6th integral
      intt6=function(x){(1/x)*log(x)*exp(-x)}
      integd66=array(0,length(bb1))
      for (i in 1:length(bb1)) {
        l=bb1[i]
        u=Inf
        integd66[i]=integrate(intt6, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-5)$value
      }
      ii66=bb1*((log(bb1))^2)*exp(-bb1)+2*ii22+((log(bb1))^2)*exp(-bb1)+2*integd66;ii66
      
      
      #===================================================================#
      #                                                                   #
      #                        Expectation terms                          #
      #                                                                   #
      #===================================================================# 
      xi00=((ml[2]^(-1/ml[1]))*(gamma((1/ml[1])+1)*pgamma((1/ml[1])+1,aa1,lower.tail = TRUE)))/D11;xi00
      xi11=(ii11-log(ml[2])*D11)/(ml[1]*D11);xi11
      xi22=(1-(1+aa1)*(1-D11))/(ml[2]*D11);xi22
      xi33=(ii22-log(ml[2])*D22)/(ml[1]*D22);xi33
      xi44=(1+bb1)/ml[2];xi44
      cc1=ml[3]*(S2-xi00);cc1
      xi55=(1-((1+cc1)*exp(-cc1)))/(ml[3]*(1-exp(-cc1)));xi55
      xi66=(ii33-log(ml[2])*(1-(1+ml[2]*(S2^ml[1]))*(1-D11)))/(ml[1]*ml[2]*D11);xi66
      xi77=(ii44+(log(ml[2])^2)*(1-(1+ml[2]*(S2^ml[1]))*(1-D11))-2*log(ml[2])*ii33)/(ml[1]*(ml[2]^2)*D11);xi77
      xi88=(ii55-log(ml[2])*((1+ml[2]*(S3^ml[1]))*D22))/(ml[1]*ml[2]*D22);xi88
      xi99=(ii66+(log(ml[2]^2))*((1+ml[2]*(S3^ml[1]))*D22)-2*log(ml[2])*ii55)/(ml[1]*(ml[2]^2)*D22);xi99
      
      #Elements of Information Matrix
      
      #Elements of I1 is obtained by taking the expectation of the negative 
      #second derivative of the log-likelihood
      a11= -n/(ml[1]^2)-sum((z1^ml[1])*(log(z1))^2)-ml[2]*(sum(xi77)-sum(xi99));a11
      a12= -sum((z1^ml[1])*log(z1))-sum(xi66)-sum(xi88);a12
      a13=0;a13
      a21=a12;a21
      a22= -n/(ml[2]^2);a22
      a23=0;a23  
      a31=0;a31
      a32=0;a32
      a33=-n2/(ml[3]^2);a33
      
      mat1=-matrix(c(a11,a12,a13,a21,a22,a23,a31,a32,a33),nrow=3,byrow = TRUE);mat1
      
      #Elements of I2 is obtained by taking the expectation of the 
      #first derivative of the log-likelihood
      gd1= n/ml[1]+sum(log(z1))-ml[2]*sum((z1^ml[1])*log(z1))+sum(xi11)-ml[2]*sum(xi66)+sum(xi33)-ml[2]*sum(xi88);gd1
      gd2= n/ml[2]-sum(z1^ml[1])-sum(xi22)-sum(xi44);gd2
      gd3= n2/ml[3]-sum(S1-z1)-sum(xi55);gd3
      
      gd=c(gd1,gd2,gd3);gd
      mat2=gd%*%t(gd);mat2
      
      var_cov=solve(mat1-mat2);var_cov #variance-covariance matrix
      
      #calculate the Confidence Intervals
      ci_l=ci_u=rep()
      for (i in 1:3) {
        ci_l[i]=ml[i]-qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
        ci_u[i]=ml[i]+qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
      }
      
      lng=c(ci_u[1]-ci_l[1],ci_u[2]-ci_l[2],ci_u[3]-ci_l[3]);lng #length of Intervals
      
      #calculating the coverage
      cp1=ifelse(ci_l[1]<alp1 && alp1<ci_u[1],1,0)
      cp2=ifelse(ci_l[2]<bet1 && bet1<ci_u[2],1,0)
      cp3=ifelse(ci_l[3]<lmd1 && lmd1<ci_u[3],1,0)
      cp=c(cp1,cp2,cp3);cp   
      
      ML=c(ml,mse,ab,lng,cp);round(ML,4) #vector of estimates, mse, ab, length and coverage probability
      
      
      #===================================================================#
      #                                                                   #
      #                        Bayesian Estimation                        #
      #                                                                   #
      #===================================================================#  
      aaa=ccc=eee=1.0;bbb=ddd=fff=0.001 #Setting the Hyper-Parameters based on Moment Matching Criteria
      b.thh1=b.thh2=b.thh3=rep()
      b.th1=alp1;b.th2=bet1;b.th3=lmd1
      it1=100000    #number of iteration 
      for(j in 1:it1){
        ww=runif(n2);ww1=runif(n3)
        tt1=(-(1/b.th2)*log(1-ww*(1-exp(-b.th2*(S2^b.th1)))))^(1/b.th1);tt1  #generating observations based on non-recall latent variable
        tt2=(S3^b.th1-(1/b.th2)*log(1-ww1))^(1/b.th1);tt2                    #generating observations based on right-censored latent variable
        
        w_i=-(1/b.th3)*log(1-ww*(1-exp(-b.th3*(S2-tt1))));w_i  #generating observations from truncated 
        
        #posterior sample from alpha
        post=function(tth1,tth2){
          yy1=(tth1^(n+aaa-1))*prod(z1^(tth1-1))*prod(tt1^(tth1-1))*prod(tt2^(tth1-1))*exp(-tth2*(sum(z1^tth1)+sum(tt1^tth1)+sum(tt2^tth1)+bbb))
          return(yy1)}
        
        last=b.th1
        cand=abs(rnorm(1,alp1,sqrt(var_cov[1,1])));cand
        r1=ifelse(is.finite(post(cand,b.th2)/post(last,b.th2))==TRUE,post(cand,b.th2)/post(last,b.th2),0.01)
        if (runif(1)<min(r1,1)) last<-cand
        b.thh1[j+1]=last
        
        #posterior samples from beta
        b.thh2[j+1]<-rgamma(1,shape=n+ccc,scale=1/(sum(z1^b.thh1[j+1])+sum(tt1^b.thh1[j+1])+sum(tt2^b.thh1[j+1])+ddd))
        #posterior samples from lambda
        b.thh3[j+1]<-rgamma(1,shape=n2+eee,scale=1/(sum(S1-z1)+sum(w_i)+fff))
        b.th1=b.thh1[j+1];b.th2=b.thh2[j+1];b.th3=b.thh3[j+1]    
      }
      b.thh1;b.thh2;b.thh3
      length(b.thh1);length(b.thh2);length(b.thh3)
      
      #discard first 10000 values to attain the stationary chain  
      ch1=b.thh1[10001:it1];ch2=b.thh2[10001:it1];ch3=b.thh3[10001:it1]
      #par(mfrow=c(2,2))
      #acf(ch1);acf(ch2);acf(ch3)          #Checking the Lag 
      zz1=seq(15,length(ch1),15)           #Taking the values at proper thin 
      ch11=ch1[zz1];ch22=ch2[zz1];ch33=ch3[zz1]  #final Markov Chain
      length(ch11);length(ch22);length(ch33)
      
      #Bayes Estimates under SELF 
      b.tth1=b.tth2=b.tth3=rep()
      b.tth1=mean(ch11);b.tth2=mean(ch22,na.rm = TRUE);b.tth3=mean(ch33);
      b.th=c(b.tth1,b.tth2,b.tth3);b.th  
      b.mse=c((b.th[1]-alp1)^2,(b.th[2]-bet1)^2,(b.th[3]-lmd1)^2);b.mse #mse
      b.bias=c(abs(b.th[1]-alp1),abs(b.th[2]-bet1),abs(b.th[3]-lmd1));b.bias #bias
      
      #Highest Posterior Density (HPD) Intervals using coda package
      h.th1=HPDinterval(mcmc(ch11));h.th2=HPDinterval(mcmc(ch22));h.th3=HPDinterval(mcmc(ch33))
      #HPD_lim=c(h.th1[,1],h.th1[,2],h.th2[,1],h.th2[,2],h.th3[,1],h.th3[,2]);HPD_lim
      
      #Length of HPD Intervals
      h.th1_l=h.th1[,2]-h.th1[,1];h.th2_l=h.th2[,2]-h.th2[,1];h.th3_l=h.th3[,2]-h.th3[,1]
      HPD_L=c(h.th1_l,h.th2_l,h.th3_l);HPD_L                               
      
      #Coverage of HPD Intervals
      coverage_HPD1=ifelse(alp1>h.th1[,1] && alp1<h.th1[,2],1,0)
      coverage_HPD2=ifelse(bet1>h.th2[,1] && bet1<h.th2[,2],1,0)                      
      coverage_HPD3=ifelse(lmd1>h.th3[,1] && lmd1<h.th3[,2],1,0)
      coverage_HPD=c(coverage_HPD1,coverage_HPD2,coverage_HPD3);coverage_HPD              
      
      #Shape of HPD Intervals
      HPD_sp1=(h.th1[,2]-b.tth1)/(b.tth1-h.th1[,1]);HPD_sp1
      HPD_sp2=(h.th2[,2]-b.tth2)/(b.tth2-h.th2[,1]);HPD_sp2
      HPD_sp3=(h.th3[,2]-b.tth3)/(b.tth3-h.th3[,1]);HPD_sp3
      HPD_sp=c(HPD_sp1,HPD_sp2,HPD_sp3);HPD_sp  
      
      Bayes=c(b.th,b.mse,b.bias,HPD_L,HPD_sp,coverage_HPD);round(Bayes,4)
    }, silent = TRUE)
    return(c(ML,Bayes))
  }
  gg()
  
  #-----------------------------------------------------------#
  #Parallel Computation of Simulation 
  
  library(parallel)
  #Determine number of cores to use
  num_cores <- detectCores() - 1
  
  #Create a cluster
  cl <- makeCluster(num_cores)
  
  start=Sys.time()
  #Run the simulation in parallel
  l=50# number simulations to run (say; 500, 1000, 2000)
  results <- parLapply(cl, 1:l, gg)
  
  #Stop the cluster
  stopCluster(cl)
  end=Sys.time()
  end-start
  mm=unlist(results)
  rr1=matrix(mm,nrow = l,ncol = 33,byrow = T) # here ncol will be equal to number of values in return 
  rr2=colMeans(rr1);rr2 #estimates based on desired number of iterations
  
  #===================================================================#
  #                                                                   #
  #                        End of R code                              #
  #                                                                   #
  #===================================================================# 
