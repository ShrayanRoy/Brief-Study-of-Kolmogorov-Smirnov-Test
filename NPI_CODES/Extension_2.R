#What happens in discrete case???
rm(list=ls())
library(dgof)
library(ggplot2)

defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))


###Does Usual Kolmogorov Fails?
#Binomial Case===============
set.seed(seed = 1234)
my.power=0
p=seq(0.1,0.9,by=0.05)
n=10
f <- function(x){pbinom(x,5,0.5)}
for(i in 1:length(p)){
  asymp_test_stat_pval <- replicate(10000,{
    our_sample <- rbinom(n,5,p[i])
    performed_test <- ks.test(our_sample,f,alternative = "two.sided",exact=TRUE)
    
    performed_test$p.value
  })
  my.power[i]=mean(asymp_test_stat_pval<0.05)
}
power_curve <- data.frame(p,power = my.power)
ggplot(power_curve,aes(x = p,y = power)) + geom_line(col="darkgoldenrod1",size = 2) +
  ggtitle(paste("Simulated Power Curve for testing \n H0: ","X~Bin(5,0.5)"," vs H1: ","X~Bin(5,p),0<p<1")) + geom_point(col="red",size=1.5)+
  labs(p,y="Simulated Power")+geom_hline(yintercept = 0.05,col="blue",size=1.2,linetype="dashed")+
  theme_bw(14) + defined_theme

###For Poisson Case
set.seed(seed = 1234)
my.power=0
lambda=1:5
n=10
f <- function(x){ppois(x,1)}
for(i in 1:length(lambda)){
  asymp_test_stat_pval <- replicate(10000,{
    our_sample <- rpois(n,lambda[i])
    performed_test <- ks.test(our_sample,f,alternative = "two.sided",exact=TRUE)
    
    performed_test$p.value
  })
  my.power[i]=mean(asymp_test_stat_pval<0.05)
}
power_curve <- data.frame(lambda,power = my.power)
ggplot(power_curve,aes(x = lambda,y = power)) + geom_line(col="darkgoldenrod1",size = 2) +
  ggtitle(paste("Simulated Power Curve for testing \n H0: ","X~Poisson(1)"," vs H1: ",
                "X~Poisson(",expression(lambda),")",expression(lambda),">1")) + 
  geom_point(col="red",size=1.5)+ labs(x=expression(lambda),y="Simulated Power") +
  geom_hline(yintercept = 0.05,col="blue",size=1.2,linetype="dashed")+
  theme_bw(14) + defined_theme  
####For Binomial Distribution==============================
set.seed(seed = 1234)
my.power=0
p=seq(0.1,0.9,by=0.05)
n=10
f <- stepfun(0:5,c(0,pbinom(0:5,5,0.5)))
      for(i in 1:length(p)){
        asymp_test_stat_pval <- replicate(10000,{
          our_sample <- rbinom(n,5,p[i])
          performed_test <- dgof::ks.test(our_sample,f,alternative = "two.sided",exact=TRUE)
        
  performed_test$p.value
        })
        my.power[i]=mean(asymp_test_stat_pval<0.05)
      }
      power_curve <- data.frame(p,power = my.power)
      ggplot(power_curve,aes(x = p,y = power)) + geom_line(col="darkgoldenrod1",size = 2) +
        ggtitle(paste("Simulated Power Curve for testing \n H0: ","X~Bin(5,0.5)"," vs H1: ","X~Bin(5,p),0<p<1")) + geom_point(col="red",size=1.5)+
        labs(p,y="Simulated Power")+geom_hline(yintercept = 0.05,col="blue",size=1.2,linetype="dashed")+
        theme_bw(14) + defined_theme
      
####For Poisson Distribution==============================
set.seed(seed = 1234)
      my.power=0
      lambda=1:5
      n=10
      f <- stepfun(0:10,c(0,ppois(0:10,1)))
      for(i in 1:length(lambda)){
        asymp_test_stat_pval <- replicate(10000,{
          our_sample <- rpois(n,lambda[i])
          performed_test <- dgof::ks.test(our_sample,f,alternative = "two.sided",exact=TRUE)
          
          performed_test$p.value
        })
        my.power[i]=mean(asymp_test_stat_pval<0.05)
      }
power_curve <- data.frame(lambda,power = my.power)
ggplot(power_curve,aes(x = lambda,y = power)) + geom_line(col="darkgoldenrod1",size = 2) +
        ggtitle(paste("Simulated Power Curve for testing \n H0: ","X~Poisson(1)"," vs H1: ",
                      "X~Poisson(",expression(lambda),")",expression(lambda),">1")) +
                  geom_point(col="red",size=1.5)+
        labs(x=expression(lambda),y="Simulated Power")+geom_hline(yintercept = 0.05,col="blue",size=1.2,linetype="dashed")+
        theme_bw(14) + defined_theme      
      
      
    
    




