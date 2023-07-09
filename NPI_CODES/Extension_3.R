######ALternatives to KS Test: Cramer-Von Misses Test and Anderson-Darling Test
library(goftest)
library(ggplot2)
library(VGAM)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))

#Distribution under H0
Asymptotic_distribution_H0 <- function(n,True_CDF,r_CDF,tag = " "){
  set.seed(seed = 1234)
  asymp_test_stat_val <- replicate(10000,{
    our_sample <- r_CDF(n)
    cvm_test <- cvm.test(our_sample,True_CDF)
    ad_test <- ad.test(our_sample,True_CDF)
    c(cvm_test$statistic,ad_test$statistic)
  })
  
  ggplot(my.df <- data.frame(test_stat = c(asymp_test_stat_val[1,],asymp_test_stat_val[2,]),
          Index = rep(c("CVM test","AD Test"),each = ncol(asymp_test_stat_val))),aes(x = test_stat,fill = Index)) + 
    geom_histogram(aes(y = ..density..),breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
    by = 1/sqrt(100)),col = 'black') + labs(x = 'sqrt(n)*Dn',y = 'Density',
    title = "Asymptotic Distribution Under H0 \n for One Sample CVM and AD test",
    subtitle = paste(tag,"and n = ",n)) + 
    theme_bw(14) + defined_theme + facet_wrap(.~Index)
  
}
graph01=Asymptotic_distribution_H0(200,function(x){punif(x,0,3)},function(n){runif(n,0,3)},
                           tag = "For U(0,3) ")
graph02=Asymptotic_distribution_H0(200,function(x){pnorm(x,0,1)},function(n){rnorm(n,0,1)},
                                   tag = "For N(0,1) ")
graph03=Asymptotic_distribution_H0(200,function(x){plaplace(x,0,1)},function(n){rlaplace(n,0,1)},
                                   tag = "For Laplace(0,1) ")
graph04=Asymptotic_distribution_H0(200,function(x){pcauchy(x,0,1)},function(n){rcauchy(n,0,1)},
                                   tag = "For C(0,1) ")
graph05=Asymptotic_distribution_H0(200,function(x){plogis(x,0,1)},function(n){rlogis(n,0,1)},
                                   tag = "For Logistic(0,1) ")
ggpubr::ggarrange(graph01,graph02,nrow=2)
ggpubr::ggarrange(graph03,graph04,graph05,nrow=3,ncol=1)

graph01=Asymptotic_distribution_H0(20,function(x){punif(x,0,3)},function(n){runif(n,0,3)},
                                   tag = "For U(0,3) ")
graph02=Asymptotic_distribution_H0(20,function(x){pnorm(x,0,1)},function(n){rnorm(n,0,1)},
                                   tag = "For N(0,1) ")
graph03=Asymptotic_distribution_H0(20,function(x){plaplace(x,0,1)},function(n){rlaplace(n,0,1)},
                                   tag = "For Laplace(0,1) ")
graph04=Asymptotic_distribution_H0(20,function(x){pcauchy(x,0,1)},function(n){rcauchy(n,0,1)},
                                   tag = "For C(0,1) ")
graph05=Asymptotic_distribution_H0(20,function(x){plogis(x,0,1)},function(n){rlogis(n,0,1)},
                                   tag = "For Logistic(0,1) ")
ggpubr::ggarrange(graph01,graph02,nrow=2)
ggpubr::ggarrange(graph03,graph04,graph05,nrow=3,ncol=1)

Power_curve.n <- function(n,H0_CDF,r_CDF_H1,H0_Pdf,H1_Pdf,nullhypo,althypo)
{
      set.seed(seed = 1234)
      my.power = matrix(0,nrow=length(n),ncol=3)
      for(i in 1:length(n)){
        asymp_test_stat_p.val <- replicate(100,{
          our_sample <- r_CDF_H1(n[i])
          ks_test_p.val <- ks.test(our_sample,H0_CDF,alternative = "two.sided")$p.value
          cvm_test_p.val<- cvm.test(our_sample,H0_CDF)$p.value
          ad_test_p.val<- ad.test(our_sample,H0_CDF)$p.value
          c(ks_test_p.val,cvm_test_p.val,ad_test_p.val)
        })
        my.power[i,]=apply(asymp_test_stat_p.val,1,function(x){mean(x<=0.05)})
      }
      power_curve <- data.frame(n=c(n,n,n),power=as.vector(my.power),Index=factor(rep(c("Kolmogorov Test",
                                                            "Cramer-Von Mises Test",
                                                            "Anderson-Darling Test"),each=length(n)),
                                                           levels= c("Kolmogorov Test",
                                                                     "Cramer-Von Mises Test",
                                                                     "Anderson-Darling Test")))
      graph.1 <- ggplot(power_curve,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
        ggtitle(paste("Simulated Power vs sample size for testing \n H0: ",nullhypo," vs H1: ",althypo)) + 
        geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
        theme_bw(14) + defined_theme

      graph.2 <- ggplot() + xlim(-5,5) + geom_function(fun = H1_Pdf,aes(colour = "H1 PDF"),size = 1.2) + 
        geom_function(fun = H0_Pdf,aes(colour = "H0 PDF"),size = 1.2) + labs(x = "x","PDF",
                                                                             title = "Plot of PDF's considered under H0 and H1",subtitle = "") + theme_bw(14) + defined_theme 
      
      ggpubr::ggarrange(graph.1,graph.2)
}

Power_curve.n(seq(10,200,by = 10),function(x){pnorm(x,0,1)},
              function(n){ qnorm(runif(n,pnorm(-1),pnorm(1)))},
              function(x){dnorm(x,0,1)},
              function(x){ifelse(abs(x) < 1,dnorm(x)/(pnorm(1) - pnorm(-1)),0)}," X ~ N(0,1)","X~Truncated Normal(-3,3)")

Power_curve.n(seq(10,200,by = 10),function(x){pnorm(x,0,1)},
              function(n){ rlaplace(n,0,1)},
              function(x){dnorm(x,0,1)},
              function(x){dlaplace(x,0,1)}," X ~ N(0,1)","X~Laplace(0,1)")

Power_curve.n(seq(10,200,by = 10),function(x){pcauchy(x,0,1)},
              function(n){ rlogis(n,0,1)},
              function(x){dcauchy(x,0,1)},
              function(x){dlogis(x,0,1)}," X ~ C(0,1)","X~Logistic(0,1)")


##Partial Specification=============================================
Power_Function_Normal <- function(n,r_CDF_H1)
{
  set.seed(seed = 1234)
  asymp_p.val <- replicate(1000,{
    my_sample = r_CDF_H1(n)
    Fx= function(x){pnorm(x,mean=mean(my_sample),sd = 1)}
    ks_test_p.val<- ks.test(my_sample, Fx,alternative = "two.sided")$p.value
    cvm_test_p.val<- cvm.test(my_sample,Fx)$p.value
    ad_test_p.val<- ad.test(my_sample,Fx)$p.value
    c(ks_test_p.val,cvm_test_p.val,ad_test_p.val)
  })
  apply(asymp_p.val,1,function(x){mean(x<=0.05)})
}
n=seq(10,200,by = 10)
pow_c.1=matrix(0,nrow=length(n),ncol=3)
for(i in 1:length(n))
{
  pow_c.1[i,]=Power_Function_Normal(n[i],function(n){rcauchy(n,2,1)})
}
power_curve.1 <- data.frame(n=c(n,n,n),power=as.vector(pow_c.1),Index=factor(rep(c("Kolmogorov Test",
                                                                                  "Cramer-Von Mises Test",
                                                                                  "Anderson-Darling Test"),each=length(n)),
                                                                            levels= c("Kolmogorov Test",
                                                                                      "Cramer-Von Mises Test",
                                                                                      "Anderson-Darling Test")))
graph.1 <- ggplot(power_curve.1,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
  ggtitle(paste("Simulated Power vs sample size for testing \n H0: ","N(mu,1)"," vs H1: ","C(mu,1)")) + 
  geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
  theme_bw(14) + defined_theme
graph.1

pow_c.2=matrix(0,nrow=length(n),ncol=3)
for(i in 1:length(n))
{
  pow_c.2[i,]=Power_Function_Normal(n[i],function(n){rlaplace(n,2,1)})
}
power_curve.2 <- data.frame(n=c(n,n,n),power=as.vector(pow_c.2),Index=factor(rep(c("Kolmogorov Test",
                                                                               "Cramer-Von Mises Test",
                                                                               "Anderson-Darling Test"),each=length(n)),
                                                                         levels= c("Kolmogorov Test",
                                                                                   "Cramer-Von Mises Test",
                                                                                   "Anderson-Darling Test")))
graph.2 <- ggplot(power_curve.2,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
  ggtitle(paste("Simulated Power vs sample size for testing \n H0: ","N(mu,1)"," vs H1: ","Laplace(mu,1)")) + 
  geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
  theme_bw(14) + defined_theme
graph.2

#power curve for Partially specified null
Power_Function1 <- function(n,r_CDF_H1)
{
  set.seed(seed = 1234)
  asymp_p_val <- replicate(10000,{
    my_sample <- r_CDF_H1(n)
    Fx=function(x){pnorm(x,mean(my_sample),sd(my_sample))}
    ks_test_p.val<- ks.test(my_sample, Fx,alternative = "two.sided")$p.value
    cvm_test_p.val<- cvm.test(my_sample,Fx)$p.value
    ad_test_p.val<- ad.test(my_sample,Fx)$p.value
    c(ks_test_p.val,cvm_test_p.val,ad_test_p.val)
  })
  apply(asymp_p_val,1,function(x){mean(x<=0.05)})
}

#taking the alternative as Cauchy with scale and location
n=seq(10,200,by = 10)
pow_c.3=matrix(0,nrow=length(n),ncol=3)
for(i in 1:length(n))
{
  pow_c.3[i,]=Power_Function1(n[i],function(n){rcauchy(n,0,1)})
}
power_curve.3 <- data.frame(n=c(n,n,n),power=as.vector(pow_c.3),Index=factor(rep(c("Kolmogorov Test",
                                                                                   "Cramer-Von Mises Test",
                                                                                   "Anderson-Darling Test"),each=length(n)),
                                                                             levels= c("Kolmogorov Test",
                                                                                       "Cramer-Von Mises Test",
                                                                                       "Anderson-Darling Test")))
graph.3 <- ggplot(power_curve.3,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
  ggtitle(paste("Simulated Power vs sample size for testing \n H0: ","N(mu,sigma^2)"," vs H1: ","C(mu,sigma)")) + 
  geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
  theme_bw(14) + defined_theme
graph.3
#taking the alternative as Laplace with scale and location
n=seq(10,200,by = 10)
pow_c.4=matrix(0,nrow=length(n),ncol=3)
for(i in 1:length(n))
{
  pow_c.4[i,]=Power_Function1(n[i],function(n){rlaplace(n,0,1)})
}
power_curve.4 <- data.frame(n=c(n,n,n),power=as.vector(pow_c.4),Index=factor(rep(c("Kolmogorov Test",
                                                                                   "Cramer-Von Mises Test",
                                                                                   "Anderson-Darling Test"),each=length(n)),
                                                                             levels= c("Kolmogorov Test",
                                                                                       "Cramer-Von Mises Test",
                                                                                       "Anderson-Darling Test")))
graph.4 <- ggplot(power_curve.4,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
  ggtitle(paste("Simulated Power vs sample size for testing \n H0: ","N(mu,sigma^2)"," vs H1: ","Laplace(mu,sigma)")) + 
  geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
  theme_bw(14) + defined_theme
graph.4


#power curve for Partially specified null
Power_Function2 <- function(n,r_CDF_H1)
{
  set.seed(seed = 1234)
  asymp_p_val <- replicate(10000,{
    my_sample <- r_CDF_H1(n)
    Fx=function(x){pcauchy(x,median(my_sample),as.vector((quantile(my_sample,0.75)-quantile(my_sample,0.25))/2))}
    ks_test_p.val<- ks.test(my_sample, Fx,alternative = "two.sided")$p.value
    cvm_test_p.val<- cvm.test(my_sample,Fx)$p.value
    ad_test_p.val<- ad.test(my_sample,Fx)$p.value
    c(ks_test_p.val,cvm_test_p.val,ad_test_p.val)
  })
  apply(asymp_p_val,1,function(x){mean(x<=0.05)})
}

#taking the alternative as normal with scale and location
n=seq(50,200,by = 10)
pow_c.5=matrix(0,nrow=length(n),ncol=3)
for(i in 1:length(n))
{
  pow_c.5[i,]=Power_Function2(n[i],function(n){rnorm(n,0,1)})
}
power_curve.5 <- data.frame(n=c(n,n,n),power=as.vector(pow_c.5),Index=factor(rep(c("Kolmogorov Test",
                                                                                   "Cramer-Von Mises Test",
                                                                                   "Anderson-Darling Test"),each=length(n)),
                                                                             levels= c("Kolmogorov Test",
                                                                                       "Cramer-Von Mises Test",
                                                                                       "Anderson-Darling Test")))
graph.5 <- ggplot(power_curve.5,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
  ggtitle(paste("Simulated Power vs sample size for testing \n H0: ","C(mu,sigma)"," vs H1: ","N(mu,sigma^2)")) + 
  geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
  theme_bw(14) + defined_theme
graph.5
#taking the alternative as Laplace with scale and location
n=seq(20,200,by = 10)
pow_c.6=matrix(0,nrow=length(n),ncol=3)
for(i in 1:length(n))
{
  pow_c.6[i,]=Power_Function2(n[i],function(n){rlaplace(n,0,1)})
}
power_curve.6 <- data.frame(n=c(n,n,n),power=as.vector(pow_c.6),Index=factor(rep(c("Kolmogorov Test",
                                                                                   "Cramer-Von Mises Test",
                                                                                   "Anderson-Darling Test"),each=length(n)),
                                                                             levels= c("Kolmogorov Test",
                                                                                       "Cramer-Von Mises Test",
                                                                                       "Anderson-Darling Test")))
graph.6 <- ggplot(power_curve.6,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
  ggtitle(paste("Simulated Power vs sample size for testing \n H0: ","C(mu,sigma)"," vs H1: ","Laplace(mu,sigma)")) + 
  geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
  theme_bw(14) + defined_theme
graph.6
