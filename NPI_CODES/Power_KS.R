rm(list = ls(all = T))  #removes all objects
library(ggplot2)
library(VGAM)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))

#============================ Power Function ====================================

Power_Function <- function(n,H1_CDF,r_CDF_H1,H0_CDF,alpha.level){
#  set.seed(seed = 987654321)
  asymp_test_stat_val <- replicate(10000,{
    our_sample <- r_CDF_H1(n)
    performed_test <- ks.test(our_sample,H0_CDF,alternative = "two.sided",exact=TRUE)
    performed_test$p.value
  })
  mean(asymp_test_stat_val<0.05)
}

Power_Function(5,function(x){pnorm(x,0,1)},function(n){rnorm(n,0,1)},
               function(x){pnorm(x,0,1)},alpha.level = 0.05)
Power_Function(20,function(x){pnorm(x,0,1)},function(n){rnorm(n,0,1)},
               function(x){pnorm(x,0,1)},alpha.level = 0.05)
Power_Function(50,function(x){pnorm(x,0,1)},function(n){rnorm(n,0,1)},
               function(x){pnorm(x,0,1)},alpha.level = 0.05)
Power_Function(100,function(x){pnorm(x,0,1)},function(n){rnorm(n,0,1)},
               function(x){pnorm(x,0,1)},alpha.level = 0.05)

Power_Function(5,function(x){plogis(x,0,1)},function(n){rlogis(n,0,1)},
               function(x){plogis(x,0,1)},alpha.level = 0.05)
Power_Function(20,function(x){plogis(x,0,1)},function(n){rlogis(n,0,1)},
               function(x){plogis(x,0,1)},alpha.level = 0.05)
Power_Function(50,function(x){plogis(x,0,1)},function(n){rlogis(n,0,1)},
               function(x){plogis(x,0,1)},alpha.level = 0.05)
Power_Function(100,function(x){plogis(x,0,1)},function(n){rlogis(n,0,1)},
               function(x){plogis(x,0,1)},alpha.level = 0.05)

library(VGAM)
Power_Function(5,function(x){prayleigh(x,2)},function(n){rrayleigh(n,2)},
               function(x){prayleigh(x,2)},alpha.level = 0.05)
Power_Function(20,function(x){prayleigh(x,2)},function(n){rrayleigh(n,2)},
               function(x){prayleigh(x,2)},alpha.level = 0.05)
Power_Function(50,function(x){prayleigh(x,2)},function(n){rrayleigh(n,2)},
               function(x){prayleigh(x,2)},alpha.level = 0.05)
Power_Function(100,function(x){prayleigh(x,2)},function(n){rrayleigh(n,2)},
               function(x){prayleigh(x,2)},alpha.level = 0.05)


#Illustration=============================================================
#i)Uniform Distribution===================================================
a <- seq(2,3,length.out = 20)
my.power <- 0
for(i in 1:length(a)){
 my.power[i] = Power_Function(50,function(x){punif(x,0,a[i])},function(n){runif(n,0,a[i])},
                    function(x){punif(x,0,3)},alpha.level = 0.05)
}

power_curve <- data.frame(a,power = my.power)
ggplot(power_curve,aes(x = a,y = power)) + geom_line(col="red",size = 1.5) +
  ggtitle("Simulated Power Curve for testing \n H0 : X~U(0,3) vs H1 : X~U(0,a),a<3") + 
  theme_bw(14) + defined_theme

