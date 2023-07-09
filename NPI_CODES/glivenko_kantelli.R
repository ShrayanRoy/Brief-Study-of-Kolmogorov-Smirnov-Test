rm(list=ls())
set.seed(1234)
library(ggplot2)
library(VGAM)
library(gridExtra)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))


#===============================================================================

gk_check=function(r_CDF,F0,epsilon,R){
  n0 = replicate(R,{
  y = NULL
  n = 0
  Dn = 2*epsilon
  while(Dn[length(Dn)]>epsilon)
  {
   n=n+1
   y=c(y,r_CDF(1))
   Dn[n]=ks.test(y,function(x){F0(x)},alternative = "two.sided")$statistic
  }

  my.data=data.frame(Index=1:n,Dn)
  graph.1=ggplot(my.data,aes(x = Index,y = Dn)) + geom_line(col="darkgoldenrod1",size = 2) +
    ggtitle("Plot of Sample Size(n) vs Dn ") + geom_point(col="red",size=1.5)+
    labs(x="Sample Size (n)",y="Value of Dn")+
    theme_bw(14) + defined_theme
  print(graph.1)
  return(n)
 })
 return(n0)
}
arrangeGrob()
gk_check(r_CDF=function(n){rcauchy(n,0,1)},
         F0=function(x){pcauchy(x,0,1)},epsilon=0.01,R=10)

gk_check(r_CDF=function(n){rexp(n,1)},
         F0=function(x){pexp(x,1)},epsilon=0.01,R=10)

r_CDF=function(n){
  u=runif(n,pnorm(-2),pnorm(2))
  qnorm(u)
}
my_CDF = function(x){
  ifelse(x < -2,0,
   ifelse(x < 2,(pnorm(x) - pnorm(-2))/(pnorm(2) - pnorm(-2)),1))
}

gk_check(r_CDF,my_CDF,epsilon = 0.01,R = 10)



