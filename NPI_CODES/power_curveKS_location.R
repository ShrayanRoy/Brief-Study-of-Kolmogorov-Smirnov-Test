rm(list = ls(all = T))  #removes all objects
library(ggplot2)
library(VGAM)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))

#=======================================================================================================

Power_Function.location <- function(n,H0_CDF,r_CDF_H0,parameter,alpha.level,nullhypo,althypo,is.exact=FALSE,crit.val= FALSE,my.col){
  if(is.exact== TRUE & crit.val== FALSE)
  {
    print(noquote("WARNING: Please input critical value"))
  }else{
  set.seed(seed = 1234)
  my.power=0
  if(is.exact== FALSE){
  for(i in 1:length(parameter)){
  asymp_test_stat_val <- replicate(10000,{
    our_sample <- r_CDF_H0(n)+parameter[i]
    performed_test <- ks.test(our_sample,H0_CDF,alternative = "greater",exact=FALSE)
    sqrt(n)*performed_test$statistic
  })
  my.power[i]=mean(4*(asymp_test_stat_val)^2 >= qchisq(alpha.level,df = 2,lower.tail = F,))
  }
  power_curve <- data.frame(parameter,power = my.power)
  ggplot(power_curve,aes(x = parameter,y = power)) + geom_line(col=my.col,size = 2) +
    ggtitle(paste("Simulated Power Curve for testing H0: ",nullhypo," vs H1: ",althypo)) + geom_point(col="red",size=1.2)+
    labs(x="a",y="Simulated Power")+
    theme_bw(14) + defined_theme
  }else{
    for(i in 1:length(parameter)){
      exact_test_stat_val <- replicate(10000,{
        our_sample <- r_CDF_H0(n)+parameter[i]
        performed_test <- ks.test(our_sample,H0_CDF,alternative = "greater",exact=TRUE)
        performed_test$statistic
      })
      my.power[i]=mean(exact_test_stat_val >= crit.val)
    }
    power_curve <- data.frame(parameter,power = my.power)
    ggplot(power_curve,aes(x = parameter,y = power)) + geom_line(col=my.col,size = 2) +
      ggtitle(paste("Simulated Power Curve for testing H0: ",nullhypo," vs H1: ",althypo)) + geom_point(col="red",size=1.2)+
      labs(x = "a",y="Simulated Power")+
      theme_bw(14) + defined_theme
    
   }
 }
}
  
Power_Function.location(50,function(x){pnorm(x,0,3)},function(n){rnorm(n,0,3)},
                        seq(-1.5,0,length.out=20),0.05,"X~N(0,3)","X~N(a,3),-1.5<a<0",my.col = "yellow")
Power_Function.location(50,function(x){pnorm(x,0,1)},function(n){rnorm(n,0,1)},
               seq(-1.5,0,length.out=20),0.05,"X~N(0,1)","X~N(a,1),-1.5<a<0",my.col = "yellow")

Power_Function.location(50,function(x){pnorm(x,0,3)},function(n){rcauchy(n,0,3)},
                        seq(-1.5,0,length.out=20),0.05,"X~Cauchy(0,3)","X~Cauchy(a,3),-1.5<a<0",my.col = "aquamarine")
Power_Function.location(50,function(x){pnorm(x,0,1)},function(n){rcauchy(n,0,1)},
                        seq(-1.5,0,length.out=20),0.05,"X~Cauchy(0,1)","X~Cauchy(a,1),-1.5<a<0",my.col = "aquamarine")

Power_Function.location(50,function(x){punif(x,0,3)},function(n){runif(n,0,3)},
                        seq(-1.5,0,length.out=20),0.05,"X~Unif(0,3)","X~Unif(a,a+3),-1.5<a<0",my.col = "seagreen2")
Power_Function.location(50,function(x){punif(x,0,1)},function(n){runif(n,0,1)},
                        seq(-1.5,0,length.out=20),0.05,"X~Unif(0,1)","X~Unif(a,a+1),-1.5<a<0",my.col = "seagreen2")

Power_Function.location(50,function(x){pexp(x,3)},function(n){rexp(n,3)},
                        seq(-1.5,0,length.out=20),0.05,"X~Shifted Exp(0,3)","X~Shifted Exp(a,3),-1.5<a<0",my.col = "cyan1")
Power_Function.location(50,function(x){pexp(x,1)},function(n){rexp(n,1)},
                        seq(-1.5,0,length.out=20),0.05,"X~Shifted (0,1)","X~Shifted Exp(a,1),-1.5<a<0",my.col = "cyan1")


#Power_Function.location(10,function(x){pnorm(x,0,3)},function(n){rnorm(n,0,3)},
#                        seq(-2,0,length.out=20),0.05,"X~N(0,3)","X~N(a,3),-2<a<0",is.exact = TRUE)

#==============================================================================================

Power_Function.scale <- function(n,H0_CDF,r_CDF_H0,parameter,alpha.level,nullhypo,althypo,is.exact= FALSE,crit.val=FALSE,my.col){
  if(is.exact==TRUE & crit.val == FALSE){
    print(noquote("WARNING: Please input critical value"))
  }else{
    if(is.exact==FALSE){
  set.seed(seed = 1234)
  my.power=0
  for(i in 1:length(parameter)){
    asymp_test_stat_val <- replicate(10000,{
      our_sample <- r_CDF_H0(n)/parameter[i]
      performed_test <- ks.test(our_sample,H0_CDF,alternative = "greater")
      sqrt(n)*performed_test$statistic
    })
    my.power[i]=mean(4*(asymp_test_stat_val)^2 >= qchisq(alpha.level,df = 2,lower.tail = F))
  }
  power_curve <- data.frame(parameter,power = my.power)
  ggplot(power_curve,aes(x = parameter,y = power)) + geom_line(col=my.col,size = 1.5) +
    ggtitle(paste("Simulated Power Curve for testing \n H0: ",nullhypo," vs H1: ",althypo)) + 
    geom_point(col="red",size=1.5) + labs(x=expression(theta),y="Simulated Power")+
    theme_bw(14) + defined_theme
    }else{set.seed(seed = 1234)
      my.power=0
      for(i in 1:length(parameter)){
        exact_test_stat_val <- replicate(10000,{
          our_sample <- r_CDF_H0(n)/parameter[i]
          performed_test <- ks.test(our_sample,H0_CDF,alternative = "greater",exact=TRUE)
          performed_test$statistic
        })
        my.power[i]=mean(exact_test_stat_val >= crit.val)
      }
      power_curve <- data.frame(parameter,power = my.power)
      ggplot(power_curve,aes(x = parameter,y = power)) + geom_line(col="darkgoldenrod1",size = 1.5) +
        ggtitle(paste("Simulated Power Curve for testing \n H0: ",nullhypo," vs H1: ",althypo)) + 
        geom_point(col="red",size=1.5) + labs(x=expression(theta),y="Simulated Power")+
        theme_bw(14) + defined_theme + geom_hline(yintercept = 0.05,col = "black",linewidth = 1.1,linetype="dashed")
    
    }
  }
}
Power_Function.scale(40,function(x){punif(x,0,3)},function(n){runif(n,0,2)},
                        seq(1,1.5,length.out=20),0.05,"X~U(0,2)","X~U(0,2/a),1<a<1.5",
                     "yellow")
Power_Function.scale(40,function(x){punif(x,0,3)},function(n){runif(n,1,3)},
                     seq(1,1.5,length.out=20),0.05,"X~U(1,3)","X~U(1,1 + 2/a),1<a<1.5","yellow")

Power_Function.scale(40,function(x){pexp(x,1)},function(n){rexp(n,1)},
                     seq(1,1.5,length.out = 20),0.05,"X~Shifted Exp(0,1)","X~Shifted Exp(0,1/a),1 < a < 1.5",
                     "aquamarine")
Power_Function.scale(40,function(x){pexp(x,1)},function(n){rexp(n,1)},
                     seq(1,1.5,length.out = 20),0.05,"X~Shifted Exp(2,1)","X~1/a),1 < a < 1.5",'aquamarine')

Power_Function.scale(40,function(x){ppareto(x,scale=1,shape=2)},function(n){rpareto(n,scale=1,shape=2)},
                     seq(1,1.5,length.out = 20),0.05,"X~Pareto(scale=1,shape=2)","X~Pareto(scale=1/a,shape=2),1 < a < 1.5",
                     "seagreen2")

Power_Function.scale(40,function(x){prayleigh(x,1)},function(n){rrayleigh(n,1)},
                     seq(1,1.5,length.out = 20),0.05,"X~Rayleigh(0,1)","X~Rayleigh(1/a),1 < a < 1.5",
"seagreen2")

#================================================================================

Power_curve.n <- function(n,H0_CDF,r_CDF_H1,H0_Pdf,H1_Pdf,alpha.level,nullhypo,althypo,is.exact=FALSE,crit.val=FALSE)
{
  if(is.exact==TRUE & max(crit.val)==FALSE)
   {print(noquote("WARNING:Please input critical value"))}
  else{
    if(is.exact==FALSE){
      set.seed(seed = 1234)
      my.power = 0
      for(i in 1:length(n)){
         asymp_test_stat_p.val <- replicate(10000,{
         our_sample <- r_CDF_H1(n[i])
         test_p.val <- ks.test(our_sample,H0_CDF,alternative = "two.sided")$p.value
         test_p.val
      })
      my.power[i]=mean(asymp_test_stat_p.val <= 0.05)
     }
      power_curve <- data.frame(n,power = my.power)
      graph.1 <- ggplot(power_curve,aes(x = n,y = power)) + geom_line(col="darkgoldenrod1",size = 2) +
        ggtitle(paste("Simulated Power vs sample size for testing \n H0: ",nullhypo," vs H1: ",althypo)) + 
        geom_point(col="red",size=1.5) + labs(x = "Sample Size(n)",y = "Simulated Power") +
        theme_bw(14) + defined_theme
    }else{
       set.seed(seed = 1234)
       my.power = 0
       for(i in 1:length(n)){
         exact_test_stat_val <- replicate(10000,{
           our_sample <- r_CDF_H1(n[i])
           test.val <- ks.test(our_sample,H0_CDF,alternative = "two.sided",exact=TRUE)$statistic
           test.val
         })
         my.power[i]=mean(exact_test_stat_val >=crit.val[i])
       }
       power_curve <- data.frame(n,power = my.power)
       graph.1 <- ggplot(power_curve,aes(x = n,y = power)) + geom_line(col="darkgoldenrod1",size = 2) +
         ggtitle(paste("Simulated Power vs sample size for testing \n H0: ",nullhypo," vs H1: ",althypo)) + 
         geom_point(col="red",size=1.5) + labs(x="Sample Size(n)",y="Simulated Power")+
         theme_bw(14) + defined_theme
     }
    graph.2 <- ggplot() + xlim(-30,30) + geom_function(fun = H1_Pdf,aes(colour = "H1 CDF"),size = 1.2) + 
      geom_function(fun = H0_Pdf,aes(colour = "H0 PDF"),size = 1.2) + labs(x = "x","PDF",
      title = "Plot of PDF's considered under H0 and H1",subtitle = "") + theme_bw(14) + defined_theme 
    
    ggpubr::ggarrange(graph.1,graph.2)
  }
}

library(VGAM)
Power_curve.n(seq(50,200,by=10),function(x){pnorm(x,0,1)},function(n){rcauchy(n,0,1)},
              function(x){dnorm(x,0,1)},function(x){dcauchy(x,0,1)},0.05,"X~N(0,1)","X~C(0,1)")
Power_curve.n(seq(50,200,by=10),function(x){pnorm(x,0,1)},function(n){rlaplace(n,0,1)},
              function(x){dnorm(x,0,1)},function(x){dlaplace(x,0,1)},0.05,"X~N(0,1)","X~Laplace(0,1)")

Power_curve.n(seq(50,200,by=10),function(x){pcauchy(x,0,1)},function(n){rlogis(n,0,1)},
              function(x){dcauchy(x,0,1)},
              function(x){dlogis(x,0,1)},0.05,"X~C(0,1)","X~Logistic(0,1)")
Power_curve.n(seq(50,200,by=10),function(x){pcauchy(x,0,1)},
              function(n){rnorm(n,0,1)},
              function(x){dcauchy(x,0,1)},
              function(x){dnorm(x,0,1)},
              0.05,"X~C(0,1)","X~N(0,1)")
Power_curve.n(c(5,30,55,80,100),function(x){pnorm(x,0,1)},
              function(n){rcauchy(n,0,1)},
              function(x){dnorm(x,0,1)},
              function(x){dcauchy(x,0,1)},
              0.05,"X~N(0,1)","X~C(0,1)",is.exact = TRUE,crit.val = c(0.56327,0.40925,0.33760,0.29407,0.26402))

Power_curve.n(seq(10,200,by = 10),function(x){pnorm(x,0,1)},
              function(n){y = rnorm(n,0,1);ifelse(y < 3,y,3)},
              function(x){qnorm(x,0,1)},
              function(x){ifelse(x < 3,dnorm(x)/pnorm(3),0)},0.05," X ~ N(0,1)","X~F")

r_CDF=function(n){
  u=runif(n,pnorm(-3),pnorm(3))
  qnorm(u)
}
hist(r_CDF(10000),freq = F,breaks = sqrt(10000),)
Power_curve.n(seq(20,200,by=10),function(x){pnorm(x,0,1)},
              function(n){r_CDF(n)},function(x){dnorm(x,0,1)},
              function(x){ifelse(abs(x) < 3,dnorm(x)/(pnorm(3) - pnorm(-3)),0)},
              0.05,"X~N(0,1)","X~F1")
#======================================================================================
