rm(list = ls(all = T))  #removes all objects
library(ggplot2)
library(VGAM)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))

#============================ Distribution Under H0 ================================

Exact_distribution_H0 <- function(n,True_CDF,r_CDF,tag = " "){
  set.seed(seed = 1234)
  test_stat_val <- replicate(10000,{
    our_sample <- r_CDF(n)
    performed_test <- ks.test(our_sample,True_CDF,alternative = "two.sided")
    performed_test$statistic
  })
  
  ggplot(my.df <- data.frame(test_stat_val),aes(test_stat_val)) + 
    geom_histogram(aes(y = ..density..),breaks = seq(0,1,by = 1/sqrt(10000)),
    col = 'black',fill = "cyan") + labs(x = 'Dn',y = 'Density',
    title = "Simulated Exact Distribution Under H0 of \n One Sample Kolmogorov-Smirnov Test Statistic Dn",
    subtitle = paste(tag,"and n = ",n)) + theme_bw(14) + defined_theme
}

#========================= Asymptotic Distribution Under H0 and H1 ==================
emp_pdf=function(z)
{
  epsilon=0.01
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
Asymptotic_distribution_H0 <- function(n,True_CDF,r_CDF,tag = " "){
  set.seed(seed = 1234)
  asymp_test_stat_val <- replicate(10000,{
    our_sample <- r_CDF(n)
    performed_test <- ks.test(our_sample,True_CDF,alternative = "two.sided")
    sqrt(n)*performed_test$statistic
  })
  
  ggplot(my.df <- data.frame(asymp_test_stat_val),aes(asymp_test_stat_val)) + 
    geom_histogram(aes(y = ..density..),breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
    by = 1/sqrt(100)),col = 'black',fill = "aquamarine") + labs(x = 'sqrt(n)*Dn',y = 'Density',
    title = "Asymptotic Distribution Under H0 \n for One Sample Kolmogorov-Smirnov Test Statistic Dn",
    subtitle = paste(tag,"and n = ",n)) + geom_function(fun =function(x){vapply(x,FUN= emp_pdf,FUN.VALUE = 2)} ,col = "red")+
    theme_bw(14) + defined_theme
  
}

Asymptotic_distribution_H1 <- function(n,True_CDF,r_CDF,H0_CDF,tag = " "){
  set.seed(seed = 1234)
  asymp_test_stat_val <- replicate(10000,{
    our_sample <- r_CDF(n)
    performed_test <- ks.test(our_sample,H0_CDF,alternative = "two.sided")
    sqrt(n)*performed_test$statistic
  })
  
  graph.1 <- ggplot(my.df <- data.frame(asymp_test_stat_val),aes(asymp_test_stat_val)) + 
    geom_histogram(aes(y = ..density..),breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
    by = 1/sqrt(100)),col = 'red',fill = "yellow") + labs(x = 'sqrt(n)*Dn',y = 'Density',
    title = "Asymptotic Distribution Under H1 for \n One Sample Kolmogorov-Smirnov Test Statistic Dn",
    subtitle = paste(tag,"and n = ",n)) + theme_bw(14) + defined_theme
  
  graph.1
  
}

#============= Distribution of P-Value ===================== 
PVal_Distribution <- function(n,rCDF,H0_CDF,tag = " "){
  set.seed(seed = 1234)
  p_val <- replicate(10000,{
    our_sample <- rCDF(n)
    ks.test(our_sample,H0_CDF,alternative = "two.sided")$p.value
  })
  
  ggplot(my.df <- data.frame(p_val),aes(p_val)) + 
    geom_histogram(aes(y = ..density..),breaks = seq(0,1 + 1/sqrt(100),
                                                     by = 1/sqrt(100)),col = 'black',fill = "cyan") + labs(x = 'sqrt(n)*Dn',y = 'Density',
                                                                                                           title = "Distribution of P-Value of \n One Sample Kolmogorov-Smirnov Test Statistic Dn",
                                                                                                           subtitle = paste(tag,"and n = ",n)) +
    theme_bw(14) + defined_theme
}

PVal_Distribution(50,function(n){rnorm(n,0,1)},function(x){pnorm(x,0,1)},"H0: X ~ N(0,1)")
PVal_Distribution(50,function(n){rnorm(n,0,1)},function(x){pnorm(x,0,2)},"H0: X ~ N(0,1) vs. H1 : X ~ N(0,2), H1 True")


#i) Uniform Distribution===================================================================================

Asymptotic_distribution_H0(200,function(x){punif(x,0,3)},function(n){runif(n,0,3)},
                           tag = "For U(0,3) ")
Asymptotic_distribution_H1(200,function(x){punif(x,-1,5)},function(n){runif(n,-1,5)},
                           function(x){punif(x,0,3)},tag = "H0: U(0,3), Truely U(-1,5)")

Asymptotic_distribution_H1(200,function(x){pnorm(x,1.5,sqrt(0.75))},function(n){rnorm(n,1.5,sqrt(0.75))},
                           function(x){punif(x,0,3)},tag = "H0: U(0,3), Truely N(1.5,sqrt(0.75))")

#ii)Normal Distribution ================================================================================
Asymptotic_distribution_H0(200,function(x){pnorm(x,4,3)},function(n){rnorm(n,4,3)},
                           tag = "For N(4,3) ")

Asymptotic_distribution_H1(200,function(x){pnorm(x,3,5)},function(n){rnorm(n,3,5)},
                           function(x){pnorm(x,4,3)},tag = "H0: N(4,3), Truely N(3,5)")

Asymptotic_distribution_H1(200,function(x){pnorm(x,4,1)},function(n){rnorm(n,4,1)},
                           function(x){pnorm(x,4,3)},tag = "H0: N(4,3), Truely N(4,1)")

Asymptotic_distribution_H1(200,function(x){pcauchy(x,4,3)},function(n){rcauchy(n,4,3)},
                           function(x){pnorm(x,4,3)},tag = "H0: N(4,3), Truely C(4,3)")
Asymptotic_distribution_H1(200,function(x){plogis(x,4,3)},function(n){rlogis(n,4,3)},
                           function(x){pnorm(x,4,3)},tag = "H0: N(4,3), Truely Logistic(4,3)")
Asymptotic_distribution_H1(200,function(x){plaplace(x,4,3)},function(n){rlaplace(n,4,3)},
                           function(x){pnorm(x,4,3)},tag = "H0: N(4,3), Truely Laplace(4,3)")

#iii)Cauchy Distribution=======================================================================
Asymptotic_distribution_H0(200,function(x){pcauchy(x,4,3)},function(n){rcauchy(n,4,3)},
                           tag = "For C(4,3) ")

Asymptotic_distribution_H1(200,function(x){pcauchy(x,3,5)},function(n){rcauchy(n,3,5)},
                           function(x){pcauchy(x,4,3)},tag = "H0: C(4,3), Truely C(3,5)")


Asymptotic_distribution_H1(200,function(x){pcauchy(x,4,1)},function(n){rcauchy(n,4,1)},
                           function(x){pcauchy(x,4,3)},tag = "H0: C(4,3), Truely C(4,1)")

Asymptotic_distribution_H1(200,function(x){pnorm(x,4,3)},function(n){rnorm(n,4,3)},
                           function(x){pcauchy(x,4,3)},tag = "H0: C(4,3), Truely N(4,3)")
Asymptotic_distribution_H1(200,function(x){plaplace(x,4,3)},function(n){rlaplace(n,4,3)},
                           function(x){pcauchy(x,4,3)},tag = "H0: C(4,3), Truely Laplace(4,3)")
Asymptotic_distribution_H1(200,function(x){plogis(x,4,3)},function(n){rlogis(n,4,3)},
                           function(x){pcauchy(x,4,3)},tag = "H0: C(4,3), Truely Logistic(4,3)")
#iv)Double Exponential distribution===============================================================================================

Asymptotic_distribution_H0(200,function(x){plaplace(x,4,3)},function(n){rlaplace(n,4,3)},
                           tag = "For Laplace(4,3) ")

Asymptotic_distribution_H1(200,function(x){plaplace(x,3,5)},function(n){rlaplace(n,3,5)},
                           function(x){plaplace(x,4,3)},tag = "H0: Laplace(4,3), Truely Laplace(3,5)")

Asymptotic_distribution_H1(200,function(x){plaplace(x,4,1)},function(n){rlaplace(n,4,1)},
                           function(x){plaplace(x,4,3)},tag = "H0: Laplace(4,3), Truely Laplace(4,1)")
Asymptotic_distribution_H1(200,function(x){pnorm(x,4,3)},function(n){rnorm(n,4,3)},
                           function(x){plaplace(x,4,3)},tag = "H0: Laplace(4,3), Truely N(4,3)")
Asymptotic_distribution_H1(200,function(x){pcauchy(x,4,3)},function(n){rcauchy(n,4,3)},
                           function(x){plaplace(x,4,3)},tag = "H0: Laplace(4,3), Truely C(4,3)")
Asymptotic_distribution_H1(200,function(x){plogis(x,4,3)},function(n){rlogis(n,4,3)},
                           function(x){plaplace(x,4,3)},tag = "H0: Laplace(4,3), Truely Logistic(4,3)")
#v)Logistic Distribution==========================================================

Asymptotic_distribution_H0(200,function(x){plogis(x,4,3)},function(n){rlogis(n,4,3)},
                           tag = "For Logistic(4,3) ")

Asymptotic_distribution_H1(200,function(x){plogis(x,3,5)},function(n){rlogis(n,3,5)},
                           function(x){plogis(x,4,3)},tag = "H0: Logistic(4,3), Truely Logistic(3,5)")

Asymptotic_distribution_H1(200,function(x){plogis(x,4,1)},function(n){rlogis(n,4,1)},
                           function(x){plogis(x,4,3)},tag = "H0: Logistic(4,3), Truely Logistic(4,1)")
Asymptotic_distribution_H1(200,function(x){pnorm(x,4,3)},function(n){rnorm(n,4,3)},
                           function(x){plogis(x,4,3)},tag = "H0: Logistic(4,3), Truely N(4,3)")
Asymptotic_distribution_H1(200,function(x){pcauchy(x,4,3)},function(n){rcauchy(n,4,3)},
                           function(x){plogis(x,4,3)},tag = "H0: Logistic(4,3), Truely C(4,3)")
Asymptotic_distribution_H1(200,function(x){plaplace(x,4,3)},function(n){rlaplace(n,4,3)},
                           function(x){plogis(x,4,3)},tag = "H0: Logistic(4,3), Truely Laplace(4,3)")

#vi)Weibull Distribution======================================================================================

Asymptotic_distribution_H0(200,function(x){pweibull(x,shape=1,scale=3)},function(n){rweibull(n,shape=1,scale=3)},
                           tag = "For Weibull(1,3) ")

Asymptotic_distribution_H1(200,function(x){pweibull(x,shape=3,scale=3)},function(n){rweibull(n,shape=3,scale=3)},
                           function(x){pweibull(x,shape=1,scale=3)},tag = "H0: Weibull(1,3), Truely Weibull(3,3)")
Asymptotic_distribution_H1(200,function(x){pweibull(x,shape=2,scale=2)},function(n){rweibull(n,shape=2,scale=2)},
                           function(x){pweibull(x,shape=1,scale=3)},tag = "H0: Weibull(1,3), Truely Weibull(2,2)")

