rm(list = ls(all = T))  #removes all objects
library(ggplot2)
library(VGAM)
library(dplyr)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))
Power_Function1 <- function(n,r_CDF_H1)
{
  set.seed(seed = 1234)
  asymp_p_val <- replicate(10000,{
    my_sample <- r_CDF_H1(n)
    Fx=function(x){pnorm(x,mean(my_sample),sd(my_sample))}
    ks_test_p.val<- ks.test(my_sample, Fx,alternative = "two.sided")$p.value
    sw_test_p.val<-shapiro.test(my_sample)$p.value
    c(ks_test_p.val,sw_test_p.val)
  })
  apply(asymp_p_val,1,function(x){mean(x<=0.05)})
}

#taking the alternative as Laplace with scale and location
n=seq(30,200,by = 10)
pow_c.3=matrix(0,nrow=length(n),ncol=2)
for(i in 1:length(n))
{
  pow_c.3[i,]=Power_Function1(n[i],function(n){rlaplace(n,0,1)})
}
power_curve.3 <- data.frame(n=c(n,n),power=as.vector(pow_c.3),Index=factor(rep(c("Kolmogorov Test",
                                                                                 "Shapiro-Wilk Test"),each=length(n)),
                                                                           levels= c("Kolmogorov Test",
                                                                                     "Shapiro-Wilk Test")))
graph.3 <- ggplot(power_curve.3,aes(x = n,y = power,col = Index)) + geom_line(size = 1.2) +
  ggtitle(paste("Simulated Power vs sample size for testing \n H0: ","N(mu,sigma^2)"," vs H1: ","Laplace(mu,sigma)")) + 
  geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
  theme_bw(14) + defined_theme
graph.3                                                                                   
                                                                                 