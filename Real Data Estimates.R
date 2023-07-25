    
    rm(list=ls(all=TRUE))
    #===================================================================#
    #                                                                   #
    #             Importing and Filtering Data                          #
    #                                                                   #
    #===================================================================#  
    
    
    library("VGAM");library("zipfR");library("coda") # loading necessary packages
    #dat<-read.csv("location",header=T) #importing data
    #View(dat)
    
    # S1: monitoring time for recall category
    # z1: time for recall category      
    
    # S2: monitoring time for non-recall
    # S3: monitoring time for censored 
    
    # n1: sample size for recall
    # n2: sample size for non-recall   
    # n3: sample size for censored
    # n: sample size
    
    #summary(S1);summary(S2);summary(S3)     #summary of monitoring times
    aa=15 # constant for scaling data
    z1=z1/aa;S1=S1/aa;S2=S2/aa;S3=S3/aa  # Scaling of data
    
    #===================================================================#
    #                                                                   #
    #                        ML Estimate (Using EM)                     #
    #                                                                   #
    #===================================================================#  
    
    #Fit Weibull distribution on time to event to initialize parameters
    #install.packages("fitdistrplus")
    library(fitdistrplus)
    est=fitdist(z1,"weibull")
    sc=est$estimate;sc
    a1=sc[1];b1=sc[2] # initial value of parameters 
    
    lower1=min(S2)/10;lower1 #set the lower limit for non-recall
    
    m.thh1=m.thh2=m.thh3=rep()
    m.th1=a1;m.th2=b1;m.th3=0.10
    it=150
    for(j in 1:it){
      #-----------Expectation terms---------------------#
      dd1=exp(-m.th2*(lower1^m.th1))-exp(-m.th2*(S2^m.th1));dd1  # Denominator 1
      dd2=exp(-m.th2*(S3^m.th1));dd2                             # Denominator 2
      
      int0=function(u){m.th1*m.th2*(u^m.th1)*exp(-m.th2*(u^m.th1))}
      intt0=array(0,length(S2))
      for (i in 1:length(S2)) {
        l=lower1
        u=S2[i]
        intt0[i]=integrate(int0, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt0
      xi0=intt0/dd1;xi0
      
      int1=function(u){log(u)*m.th1*m.th2*(u^(m.th1-1))*exp(-m.th2*(u^m.th1))}
      intt1=array(0,length(S2))
      for (i in 1:length(S2)) {
        l=lower1
        u=S2[i]
        intt1[i]=integrate(int1,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt1
      xi1=intt1/dd1;xi1   
      
      int2=function(u){m.th1*m.th2*(u^(2*m.th1-1))*exp(-m.th2*(u^m.th1))}
      intt2=array(0,length(S2))
      for (i in 1:length(S2)) {
        l=lower1
        u=S2[i]
        intt2[i]=integrate(int2,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt2
      xi2=intt2/dd1;xi2   
      
      intt3=array(0,length(S3))
      for (i in 1:length(S3)) {
        l=S3[i]
        u=Inf
        intt3[i]=integrate(int1, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt3
      xi3=intt3/dd2;xi3   
      
      intt4=array(0,length(S3))
      for (i in 1:length(S3)) {
        l=S3[i]
        u=Inf
        intt4[i]=integrate(int2, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt4
      xi4=intt4/dd2;xi4   
      
      S4=(S2-xi0);S4
      dd3=exp(-m.th3*lower1)-exp(-m.th3*S4);dd3
      
      int5=function(u){u*m.th3*exp(-m.th3*u)}
      intt5=array(0,length(S4))
      for (i in 1:length(S4)) {
        l=lower1
        u=S4[i]
        intt5[i]=integrate(int5, l, u, subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt5
      xi5=intt5/dd3;xi5   
      
      int6=function(u){log(u)*m.th1*m.th2*(u^(2*m.th1-1))*exp(-m.th2*(u^m.th1))}
      intt6=array(0,length(S2))
      for (i in 1:length(S2)) {
        l=lower1
        u=S2[i]
        intt6[i]=integrate(int6, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt6
      xi6=intt6/dd1;xi6
      
      int8=function(u){log(u)*m.th1*m.th2*(u^(2*m.th1-1))*exp(-m.th2*(u^m.th1))}
      intt8=array(0,length(S3))
      for (i in 1:length(S3)) {
        l=S3[i]
        u=Inf
        intt8[i]=integrate(int8, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
      }
      intt8
      xi8=intt8/dd2;xi8
      
      m.thh1[j]=m.th1;m.thh2[j]=m.th2;m.thh3[j]=m.th3
      m.thh1[j+1]=n/((m.thh2[j]*sum((z1^m.thh1[j])*log(z1)))-sum(log(z1))-sum(xi1)-sum(xi3)+m.thh2[j]*(sum(xi6)+sum(xi8)))
      m.thh2[j+1]=n/(sum(z1^m.thh1[j+1])+sum(xi2)+sum(xi4))
      m.thh3[j+1]=n2/(sum(S1-z1)+sum(xi5))
      m.th1=m.thh1[j+1];m.th2=m.thh2[j+1];m.th3=m.thh3[j+1]
    }
    m.thh1;m.thh2;m.thh3
    #convergence checking
    # par(mfrow=c(2,2))
    # plot(m.thh1,type='l');plot(m.thh2,type='l');plot(m.thh3,type='l')
    ml=c(m.thh1[it],m.thh2[it],m.thh3[it]);ml #ML estimates
    
    # Mean and Median duration (based on ML estimates)
    mean1=aa*(ml[2]^(-1/ml[1]))*gamma(1+1/ml[1]);mean1
    md1=aa*ml[2]^(-1/ml[1])*(log(2))^(1/ml[1]);md1 
    
    
    #===================================================================#
    #                                                                   #
    #                 Asymptotic Confidence Intervals (ACIs)            #
    #                                                                   #
    #===================================================================# 
    
    #Calculating the Information Matrix using Missing Information Principle
    #Calculate the expressions as in EM steps on Maximum likelihood estimates (MLEs)   
    
    ddd1=exp(-ml[2]*(lower1^ml[1]))-exp(-ml[2]*(S2^ml[1]));ddd1  # Denominator 1
    ddd2=exp(-ml[2]*(S3^ml[1]));ddd2                             # Denominator 2
    
    int00=function(u){ml[1]*ml[2]*(u^ml[1])*exp(-ml[2]*(u^ml[1]))}
    intt00=array(0,length(S2))
    for (i in 1:length(S2)) {
      l=lower1
      u=S2[i]
      intt00[i]=integrate(int00, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt00
    xi00=intt00/ddd1;xi00
    
    int11=function(u){log(u)*ml[1]*ml[2]*(u^(ml[1]-1))*exp(-ml[2]*(u^ml[1]))}
    intt11=array(0,length(S2))
    for (i in 1:length(S2)) {
      l=lower1
      u=S2[i]
      intt11[i]=integrate(int11, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt11
    xi11=intt11/ddd1;xi11   
    
    int22=function(u){ml[1]*ml[2]*(u^(2*ml[1]-1))*exp(-ml[2]*(u^ml[1]))}
    intt22=array(0,length(S2))
    for (i in 1:length(S2)) {
      l=lower1
      u=S2[i]
      intt22[i]=integrate(int22, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt22
    xi22=intt22/ddd1;xi22   
    
    intt33=array(0,length(S3))
    for (i in 1:length(S3)) {
      l=S3[i]
      u=Inf
      intt33[i]=integrate(int11, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt33
    xi33=intt33/ddd2;xi33   
    
    intt44=array(0,length(S3))
    for (i in 1:length(S3)) {
      l=S3[i]
      u=Inf
      intt44[i]=integrate(int22, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt44
    xi44=intt44/ddd2;xi44   
    
    S41=(S2-xi00);S41
    ddd3=exp(-ml[3]*lower1)-exp(-ml[3]*S41);ddd3
    
    int55=function(u){u*ml[3]*exp(-ml[3]*u)}
    intt55=array(0,length(S41))
    for (i in 1:length(S41)) {
      l=lower1
      u=S41[i]
      intt55[i]=integrate(int55, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt55
    xi55=intt55/ddd3;xi55   
    
    int66=function(u){log(u)*ml[1]*ml[2]*(u^(2*ml[1]-1))*exp(-ml[2]*(u^ml[1]))}
    intt66=array(0,length(S2))
    for (i in 1:length(S2)) {
      l=lower1
      u=S2[i]
      intt66[i]=integrate(int66, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt66
    xi66=intt66/ddd1;xi66
    
    int7=function(u){(log(u)^2)*ml[1]*ml[2]*(u^(2*ml[1]-1))*exp(-ml[2]*(u^ml[1]))}
    intt77=array(0,length(S2))
    for (i in 1:length(S2)) {
      l=lower1
      u=S2[i]
      intt77[i]=integrate(int7, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt77
    xi77=intt77/ddd1;xi77
    
    int88=function(u){log(u)*ml[1]*ml[2]*(u^(2*ml[1]-1))*exp(-ml[2]*(u^ml[1]))}
    intt88=array(0,length(S3))
    for (i in 1:length(S3)) {
      l=S3[i]
      u=Inf
      intt88[i]=integrate(int88, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt88
    xi88=intt88/ddd2;xi88
    
    intt99=array(0,length(S3))
    for (i in 1:length(S3)) {
      l=S3[i]
      u=Inf
      intt99[i]=integrate(int7, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    intt99
    xi99=intt99/ddd2;xi99
    
    #Elements of Information Matrix
    
    #Elements of I1 is obtained by taking the expectation of the negative 
    #second derivative of the log-likelihood
    a11= -n/(ml[1]^2)-ml[2]*sum((z1^ml[1])*(log(z1))^2)-ml[2]*(sum(xi77)+sum(xi99));a11
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
    alpha=0.05
    ci_l=ci_u=rep()
    for (i in 1:3) {
      ci_l[i]=ml[i]-qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
      ci_u[i]=ml[i]+qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
    }
    
    lmt=c(ci_l[1],ci_u[1],ci_l[2],ci_u[2],ci_l[3],ci_u[3]);lmt #CIs limits
    se=c(sqrt(diag(var_cov)[1]),sqrt(diag(var_cov)[2]),sqrt(diag(var_cov)[3])) #standard error of estimates
    
    ML=c(ml,lmt,se,md1,mean1);round(ML,4) #ML Estimates
    
    
    #===================================================================#
    #                                                                   #
    #                        Bayesian Estimation                        #
    #                                                                   #
    #===================================================================#  
    aaa=ccc=eee=1;bbb=ddd=fff=0 #setting the hyper-parameters for making non-informative
    b.thh1=b.thh2=b.thh3=rep()
    b.th1=b.th2=b.th3=rep()
    b.th1=ml[1];b.th2=ml[2];b.th3=ml[3]  #initialize values as MLE
    it1=250000                           #iterations 
    for(j in 1:it1){
      
      ww=runif(n2);ww1=runif(n3)#uniform variates
      tt1=(-(1/b.th2)*log(exp(-b.th2*(lower1^b.th1))-ww*(exp(-b.th2*(lower1^b.th1))-exp(-b.th2*(S2^b.th1)))))^(1/b.th1);tt1  #observation from non-recall
      tt2=(S3^b.th1-(1/b.th2)*log(1-ww1))^(1/b.th1);tt2  #observation right censored
      w_i=-(1/b.th3)*log(exp(-b.th3*lower1)-ww*(exp(-b.th3*lower1)-exp(-b.th3*(S2-tt1))));w_i  # for non-recall prob
      
      #posterior sample from alpha
      post=function(tth1,tth2){
        yy1=(tth1^(n+aaa-1))*prod(z1^(tth1-1))*prod(tt1^(tth1-1))*prod(tt2^(tth1-1))*exp(-tth2*(sum(z1^tth1)+sum(tt1^tth1)+sum(tt2^tth1)+bbb))
        return(yy1)}
      
      #MH Steps
      last=b.th1
      cand=abs(rnorm(1,ml[1],0.05));cand
      r1=ifelse(is.finite(post(cand,b.th2)/post(last,b.th2))==TRUE,post(cand,b.th2)/post(last,b.th2),0.01)
      if (runif(1)<min(r1,1)) last<-cand
      b.thh1[j+1]=last
      
      #posterior sample from beta
      b.thh2[j+1]<-rgamma(1,shape=n+ccc,scale=1/(sum(z1^b.thh1[j+1])+sum(tt1^b.thh1[j+1])+sum(tt2^b.thh1[j+1])+ddd))
      #posterior sample from lambda
      b.thh3[j+1]<-rgamma(1,shape=n2+eee,scale=1/(sum(S1-z1)+sum(w_i)+fff))
      b.th1=b.thh1[j+1];b.th2=b.thh2[j+1];b.th3=b.thh3[j+1]    
    }
    b.thh1;b.thh2;b.thh3
    length(b.thh1);length(b.thh2);length(b.thh3)
    
    #trace plot of generated chains 
    par(mfrow=c(2,2))
    traceplot(mcmc(b.thh1))
    traceplot(mcmc(b.thh2))
    traceplot(mcmc(b.thh3))
    
    #discard first 35000 to attain the stationary chain
    ch1=b.thh1[35001:it1];ch2=b.thh2[35001:it1];ch3=b.thh3[35001:it1]
    
    #check lag through ACF Plot
    par(mfrow=c(2,2))
    acf(ch1);acf(ch2);acf(ch3)
    zz2=seq(70,length(ch2),70);zz3=seq(20,length(ch3),20) #Taking the values at proper thin 
    ch111=ch1[zz2];ch222=ch2[zz2];ch333=ch3[zz3]
    length(ch111);length(ch222);length(ch333)
    ch11=ch111[1:3000];ch22=ch222[1:3000];ch33=ch333[1:3000] #final Markov Chain
    
    #quantile plots for chains
    par(mfrow=c(2,2))
    gx1=as.vector(ch11);gx1
    ii=1:length(gx1)
    q1.25=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.25)));q1.25
    q1.50=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.5)));q1.50
    q1.75=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.75)));q1.75
    
    gx2=na.omit(as.vector(ch22));gx2
    ii=1:length(gx2)
    q2.25=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.25),na.rm=TRUE));q2.25
    q2.50=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.5),na.rm=TRUE));q2.50
    q2.75=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.75),na.rm=TRUE));q2.75
    
    gx3=as.vector(ch33)
    ii=1:length(gx3)
    q3.25=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.25)));q3.25
    q3.50=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.5)));q3.50
    q3.75=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.75)));q3.75
    
    #Trace Plots
    matplot(ii,cbind(q1.25,q1.50,q1.75),main="",type="l",col=1,xlab="iteration",ylab=expression(alpha))
    matplot(ii,cbind(q2.25,q2.50,q2.75),main="",type="l",col=1,xlab="iteration",ylab=expression(beta))
    matplot(ii,cbind(q3.25,q3.50,q3.75),main="",type="l",col=1,xlab="iteration",ylab=expression(lambda))
    
    #ACF Plots
    acf(ch11,main="",col="black",xlab="Lag",ylab=expression("ACF"(alpha[1])))
    acf(ch22,main="",col="black",xlab="Lag",ylab=expression("ACF"(beta[1])))
    acf(ch33,main="",col="black",xlab="Lag",ylab=expression("ACF"(lambda[1])))
    
    #MAT Plots
    matplot(ch11,type="l",col="black",xlab="iterations",ylab=expression(alpha[1]),main="")
    matplot(ch22,type="l",col="black",xlab="iterations",ylab=expression(beta[1]),main="")
    matplot(ch33,type="l",col="black",xlab="iterations",ylab=expression(lambda[1]),main="")
    
    #density plots
    plot(density(ch11),col="black",type="l",ylab=expression("Density"(alpha[1])),main="")
    plot(density(ch22),col="black",type="l",ylab=expression("Density"(beta[1])),main="")
    plot(density(ch33),col="black",type="l",ylab=expression("Density"(lambda[1])),main="")
    
    #Bayes estimates under SELF 
    b.tth1=b.tth2=b.tth3=rep()
    b.tth1=mean(ch11);b.tth2=mean(ch22);b.tth3=mean(ch33)
    b.th=c(b.tth1,b.tth2,b.tth3);b.th  
    
    Mean1=mean(aa*(ch22^(-1/ch11))*gamma(1+1/ch11));Mean1  #mean 
    Med1=mean(aa*(ch22^(-1/ch11))*(log(2)^(1/ch11)));Med1  #median 
    
    #HPD Intervals using coda package
    h.th1=HPDinterval(mcmc(ch11));h.th2=HPDinterval(mcmc(ch22));h.th3=HPDinterval(mcmc(ch33))
    HPD_lim=c(h.th1[,1],h.th1[,2],h.th2[,1],h.th2[,2],h.th3[,1],h.th3[,2]);HPD_lim
    
    
    ML=c(ml,lmt,se,mean1,md1);ML           #ML Estimates along with confidence limit, standard error, mean and median 
    Bayes=c(b.th,HPD_lim,Mean1,Med1);Bayes #Bayes Estimates along with confidence limit, mean and median
 
    #===================================================================#
    #                                                                   #
    #                        End of R code                              #
    #                                                                   #
    #===================================================================# 
    