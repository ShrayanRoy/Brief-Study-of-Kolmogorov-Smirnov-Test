#What happenned in Multivariate Case????

rm(list=ls(all=T))
library(ggplot2)
library(mvtnorm)
library(VGAM)
library(copula)
library(LaplacesDemon)
library(pracma)
library(fMultivar)

defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))



set.seed(1234)
gk_check=function(r_CDF,F0,epsilon,ini.par,lower,upper){

    our_sample = matrix(0,nrow=1,ncol=2)
    n = 0
    Dn = 2*epsilon
    while(Dn[length(Dn)]>epsilon)
    {
      n=n+1
      our_sample=rbind(our_sample,r_CDF(1))
      Bi_empi_cdf=function(par)
      {
        x=par[1]
        y=par[2]
        abs(mean((our_sample[,1]<=x) & (our_sample[,2]<=y))-F0(x,y))
        
      }
      Dn[n]=optim(ini.par,Bi_empi_cdf,lower=lower,upper=upper,
                  method="L-BFGS-B",control=list(fnscale=-1))$value
    }
    
    my.data=data.frame(Index=1:n,Dn)
    graph.1=ggplot(my.data,aes(x = Index,y = Dn)) + geom_line(col="darkgoldenrod1",size = 2) +
      ggtitle("Plot of Sample Size(n) vs Dn ") + geom_point(col="red",size=1.5)+
      labs(x="Sample Size (n)",y="Value of Dn",
           subtitle=paste("n0=",n))+
      theme_bw(14) + defined_theme
    graph.1
}

#Multivariate Normal Distribution
g1=gk_check(function(n){rmvnorm(n,c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))},
         function(x,y){pmvnorm(lower=-Inf,upper=c(x,y),mean=c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))[1]},0.01,
         ini.par=c(0,0),lower=-Inf,upper=Inf)
g1=g1+ggtitle("Plot of Dn vs n for \n Bivariate Normal distribution")

g2=gk_check(function(n){rmvnorm(n,c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))},
            function(x,y){pmvnorm(lower=-Inf,upper=c(x,y),mean=c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))[1]},0.01,
            ini.par=c(0,0),lower=-Inf,upper=Inf)
g2=g2+ggtitle("Plot of Dn vs n for \n Bivariate Normal distribution")

g3=gk_check(function(n){rmvnorm(n,c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))},
            function(x,y){pmvnorm(lower=-Inf,upper=c(x,y),mean=c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))[1]},0.01,
            ini.par=c(0,0),lower=-Inf,upper=Inf)  
g3=g3+ggtitle("Plot of Dn vs n for \n Bivariate Normal distribution")

#Multivariate Gambel Distribution
g4=gk_check(function(n){rbifgmcop(n,apar=0.5)},
         function(x,y){pbifgmcop(x,y,apar=0.5)},0.001,ini.par=c(1,1),
         lower=c(0,0),upper=Inf)
g4=g4+ggtitle("Plot of Dn vs n for \n Bivariate Gambel distribution")
g5=gk_check(function(n){rbifgmcop(n,apar=0.5)},
         function(x,y){pbifgmcop(x,y,apar=0.5)},0.001,ini.par=c(1,1),
         lower=c(0,0),upper=Inf)
g5=g5+ggtitle("Plot of Dn vs n for \n Bivariate Gambel distribution")
g6=gk_check(function(n){rbifgmcop(n,apar=0.5)},
         function(x,y){pbifgmcop(x,y,apar=0.5)},0.001,ini.par=c(1,1),
         lower=c(0,0),upper=Inf)
g6=g6+ggtitle("Plot of Dn vs n for \n Bivariate Gambel distribution")
#Multivariate Student's t Distribution
g7=gk_check(function(n){rt2d(n,rho=0.333,nu=3)},
         function(x,y){pt2d(x,y,rho=0.333,nu=3)},0.01,
         ini.par=c(1,2),lower=-Inf,upper=Inf)
g7=g7+ggtitle("Plot of Dn vs n for \n Bivariate t distribution")
g8=gk_check(function(n){rt2d(n,rho=0.333,nu=3)},
         function(x,y){pt2d(x,y,rho=0.333,nu=3)},0.01,
         ini.par=c(1,2),lower=-Inf,upper=Inf)
g8=g8+ggtitle("Plot of Dn vs n for \n Bivariate t distribution")
g9=gk_check(function(n){rt2d(n,rho=0.333,nu=3)},
         function(x,y){pt2d(x,y,rho=0.333,nu=3)},0.01,
         ini.par=c(1,2),lower=-Inf,upper=Inf)
g9=g9+ggtitle("Plot of Dn vs n for \n Bivariate t distribution")
grid.arrange(g1,g2,g3,
     g4,g5,g6,
     g7,g8,g9,nrow=3)




#So, Glivenko Cantelli holds in Multivariate Case,
#Now, let's check what about KS-test
#============================ Distribution Under H0 ================================
emp_pdf=function(z)
{
  epsilon=0.001
  a=8*z*exp(-2*z*z)
  b=a+2*epsilon
  i=2
  while(abs(b-a)>epsilon)
  {
    a=b
    b=a+8*((-1)^(i-1))*i*i*z*exp(-2*i*i*z*z)
    i=i+1
  }
  return(b)
}
Asymptotic_distribution_H0 <- function(n,F0,r_CDF,tag = " ",ini.par,lower,upper){
  set.seed(seed = 1234)
  asymp_test_stat_val <- replicate(1000,{
    our_sample <- r_CDF(n)
    Bi_empi_cdf=function(par)
    {
      x=par[1]
      y=par[2]
      abs(mean((our_sample[,1]<=x) & (our_sample[,2]<=y))-F0(x,y))
      
    }
    performed_test =optim(ini.par,Bi_empi_cdf,lower=lower,upper=upper,
                method="L-BFGS-B",control=list(fnscale=-1))$value
    sqrt(n)*performed_test
  })
  print(mean(asymp_test_stat_val>=(1.35810)))
  ggplot(my.df <- data.frame(asymp_test_stat_val),aes(asymp_test_stat_val)) + 
    geom_histogram(aes(y = ..density..),breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
                                                     by = 1/sqrt(100)),col = 'black',fill = "aquamarine")+ 
  ylim(0,1.5)+ labs(x = 'sqrt(n)*Dn',y = 'Density',
                                                                                                                 title = "Asymptotic Distribution Under H0 \n for Kolmogorov-Smirnov Test Statistic Dn",
                                                                                                                 subtitle = paste(tag,"and n = ",n)) + geom_function(fun =function(x){vapply(x,FUN= emp_pdf,FUN.VALUE = 2)} ,col = "red") + 
    theme_bw(14) + defined_theme
  
}
Asymptotic_distribution_H0(100,function(x,y){pmvnorm(lower=-Inf,upper=c(x,y),mean=c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))[1]},
                           function(n){rmvnorm(n,c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))},
                          tag="For BVN((0,0),matrix(3,1,1,3)",
                          ini.par=c(0,0),lower=-Inf,upper=Inf)


Asymptotic_distribution_H0(100,function(x,y){pt2d(x,y,rho=0.333,nu=3)},
                           function(n){rt2d(n,rho=0.333,nu=3)},
                           tag="For Bivariate t distribution(rho=0.333,nu=3)",
                           ini.par=c(0,0),lower=-Inf,upper=Inf)

#Power Curve==================================================================
#================================================================================

Power_curve.n <- function(n,F0,r_CDF_H1,alpha.level,nullhypo,althypo,crit.val,lower,upper)
{     set.seed(seed = 1234)
      my.power = 0
      for(i in 1:length(n)){
        asymp_test_stat <- replicate(1000,{
          our_sample <- r_CDF_H1(n[i])
          Bi_empi_cdf=function(par)
          {
            x=par[1]
            y=par[2]
            abs(mean((our_sample[,1]<=x) & (our_sample[,2]<=y))-F0(x,y))
         
          }
          performed_test =optim(c(mean(our_sample[,1]),mean(our_sample[,2])),Bi_empi_cdf,lower=lower,upper=upper,
                                method="L-BFGS-B",control=list(fnscale=-1))$value
          performed_test
          print(performed_test)
        })
        my.power[i]=mean(asymp_test_stat >= crit.val[i])
        print(paste(my.power[i],"======================="))
      }
      power_curve <- data.frame(n,power = my.power)
      graph.1 <- ggplot(power_curve,aes(x = n,y = power)) + geom_line(col="darkgoldenrod1",size = 2) +
        ggtitle(paste("Simulated Power vs sample size for testing \n H0: ",nullhypo," \n vs H1: ",althypo)) + 
        geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
        theme_bw(14) + defined_theme
      print(graph.1)
      return(power_curve)
}

Power_curve.n(c(10,30,50,75,100,150,200),function(x,y){pt2d(x,y,rho=0.333,nu=3)},
              function(n){rmvnorm(n,c(0,0),sigma=matrix(c(3,1,1,3),nrow=2))},
              0.05,"Bivariate t Distribution(rho=0.333,nu=3)","Bivariate Normal with same mean and Dispersion",
              c(0.40925,0.24170,0.18845,1.35810/sqrt(75),1.35810/sqrt(100),
                1.35810/sqrt(150),1.35810/sqrt(200)),-Inf,Inf)
