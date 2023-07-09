rm(list = ls(all = T))  #removes all objects
library(ggplot2)
library(VGAM)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))


#======

power_animation2 <- function(n){
  set.seed(1234)
  
  graph.data <- NULL
  
  for(i in 1:length(n)){
    
    test.stat <- replicate(3000,{
      my.s <- rnorm(n[i],0,1)
      sqrt(n[i])*ks.test(my.s,"pnorm",alternative = "two.sided")$statistic
    })
   graph.data <- rbind(graph.data,data.frame(test_stat = test.stat,Sample_size = c(n[i]),
                        Dist = c("Normal(0,1)"),crit.val = c(1.35810)))
  }
  
  for(i in 1:length(n)){
    
    test.stat <- replicate(3000,{
      my.s <- rcauchy(n[i],0,1)
      sqrt(n[i])*ks.test(my.s,"pnorm",alternative = "two.sided")$statistic
    })
    graph.data <- rbind(graph.data,data.frame(test_stat = test.stat,Sample_size = c(n[i]),
                                              Dist = c("Cauchy(0,1)"),crit.val = c(1.35810)))
  }
  
  for(i in 1:length(n)){
    
    test.stat <- replicate(3000,{
      my.s <- rlaplace(n[i],0,1)
      sqrt(n[i])*ks.test(my.s,"pnorm",alternative = "two.sided")$statistic
    })
    graph.data <- rbind(graph.data,data.frame(test_stat = test.stat,Sample_size = c(n[i]),
                                              Dist = c("Laplace(0,1)"),crit.val = c(1.35810)))
  }
  
  for(i in 1:length(n)){
    
    test.stat <- replicate(3000,{
      my.s <- rlogis(n[i],0,1)
      sqrt(n[i])*ks.test(my.s,"pnorm",alternative = "two.sided")$statistic
    })
    graph.data <- rbind(graph.data,data.frame(test_stat = test.stat,Sample_size = c(n[i]),
                                              Dist = c("Logistic(0,1)"),crit.val = c(1.35810)))
  }
  
  graph.data$Dist <- factor(graph.data$Dist,levels = c("Normal(0,1)","Cauchy(0,1)","Laplace(0,1)","Logistic(0,1)"))
  
  graph.1 <- ggplot(graph.data,aes(x = test_stat,frame = Sample_size,fill =  Dist)) +
    geom_histogram(aes(y = ..density..),breaks = seq(0,ceiling(sqrt(max(n))) + 1/sqrt(100),
        by = 1/sqrt(100)),position = "identity") + geom_vline(xintercept = 1.35810,col = "black",linetype = "dashed",
        size = 0.9) + labs(x = 'sqrt(n)*Dn',y = 'Density') + 
        theme_bw(14) + defined_theme + facet_wrap(.~Dist)
  
  plotly::ggplotly(graph.1) 
}


power_animation2(c(10,20,50,70,100,150,200,250))

#========================================================

power_animation1 <- function(n){
  
  power.df <- NULL
  a <- seq(0.5,3,by = 0.2)
  
  for(i in 1:length(a)){
    
    my.cdf <- function(x){
      ifelse(x < -a[i],ifelse(x < a[i], 
            (pnorm(x) - pnorm(-a[i]))/(pnorm(a[i]) - pnorm(-a[i])),1),0)}
    
    for(j in 1:length(n)){
      
     asymp.p_val <- replicate(1000,{
       my.s <- qnorm(runif(n[j],pnorm(-a[i]),pnorm(a[i])))
       ks.test(my.s,"pnorm",alternative = "two.sided")$p.value
     })
    
     power.df <- rbind(power.df,data.frame(Power_val = mean(asymp.p_val <= 0.05),
                  Sample_size = c(n[j]),a_val = c(a[i])))  
      
    }
  }
  
  graph.1 <- ggplot(power.df,aes(x = Sample_size,frame = a_val)) +
        geom_line(aes(y = Power_val),col = "aquamarine",size = 1.2) + 
        geom_point(aes(y =  Power_val),col = "red") + 
        labs(x = "Sample Size",y = "Empirical Power",
        title = "Empirical Power Curve for testing H0: X~N(0,1) against Truncated Normals") +
        theme_bw(14) + defined_theme
  
  library(plotly)
  plotly::ggplotly(graph.1) %>% animation_opts(
    1000, easing = "elastic", redraw = FALSE
  )
  
}

power_animation1(c(10,20,30,50,70,100,150,200,250))

#==================================================================================

Confidence_band <- function(n,x,H0_CDF,rCDF_H0,D_na,tag){
 
  my.s <- rCDF_H0(n) 
  empircial_cdf <- ecdf(my.s)
  Ln.x <- vapply(empircial_cdf(x) - D_na,FUN = function(z){max(0,z)},FUN.VALUE = 2)
  Un.x <- vapply(empircial_cdf(x) + D_na,FUN = function(z){min(1,z)},FUN.VALUE = 2)
  
  graph.data <- data.frame(x = x,val = c(Ln.x,H0_CDF(x),Un.x),
                  Index = factor(rep(c("Ln (Lower)","Fx","Un (Upper)"),each =  length(x)),
                  levels = c("Ln (Lower)","Fx","Un (Upper)")),
                  Lower = Ln.x,Upper = Un.x)
  ggplot(graph.data,aes(x,val,col = Index)) +  
    geom_ribbon(data = graph.data,aes(ymin = Lower, ymax = Upper),alpha = 0.1,col = "snow") + 
    geom_line(linewidth = 1.2) + labs(x = "x",y = "Fn and Fx",title = "Confidence Band for F(x)",
         subtitle = tag) + theme_bw(14) + defined_theme
    
}

g1 <- Confidence_band(10,seq(-5,5,by = 0.1),function(x){pnorm(x,0,1)},
                function(n){rnorm(n,0,1)},0.40925,"When H0: X ~ N(0,1) and n = 10")
g2 <- Confidence_band(10,seq(-5,5,by = 0.1),function(x){pcauchy(x,0,1)},
                function(n){rcauchy(n,0,1)},0.40925,"When H0: X ~ C(0,1) and n = 10")

g3 <- Confidence_band(20,seq(-5,5,by = 0.1),function(x){plaplace(x,0,1)},
                function(n){rlaplace(n,0,1)},0.29407,"When H0: X ~ Laplace(0,1) and n = 20")
g4 <- Confidence_band(20,seq(-5,5,by = 0.1),function(x){ifelse(x < -2,0,ifelse(x < 2,(pnorm(x) - pnorm(-2))/(pnorm(2) - pnorm(-2)),1))},
                function(n){qnorm(runif(n,pnorm(-2),pnorm(2)))},0.29407,"When H0: X ~ Trunc N(0,1) with (-2,2)and n = 20")

gridExtra::grid.arrange(g1,g2,g3,g4,ncol = 2)

#======================== for H1 cdf ====

Confidence_band_H1 <- function(n,x,H1_CDF,rCDF_H0,D_na,tag){
  
  my.s <- rCDF_H0(n) 
  empircial_cdf <- ecdf(my.s)
  Ln.x <- vapply(empircial_cdf(x) - D_na,FUN = function(z){max(0,z)},FUN.VALUE = 2)
  Un.x <- vapply(empircial_cdf(x) + D_na,FUN = function(z){min(1,z)},FUN.VALUE = 2)
  
  graph.data <- data.frame(x = x,val = c(Ln.x,H1_CDF(x),Un.x),
                           Index = factor(rep(c("Ln (Lower)","F1(x)","Un (Upper)"),each =  length(x)),
                                          levels = c("Ln (Lower)","F1(x)","Un (Upper)")),
                           Lower = Ln.x,Upper = Un.x)
  ggplot(graph.data,aes(x,val,col = Index)) +  
    geom_ribbon(data = graph.data,aes(ymin = Lower, ymax = Upper),alpha = 0.1,col = "snow") + 
    geom_line(linewidth = 1.2) + labs(x = "x",y = "Fn and F1(x)",title = "Confidence Band for F1(x)",
       subtitle = tag) + theme_bw(14) + defined_theme
  
}

g1 <- Confidence_band_H1(30,seq(-5,5,by = 0.1),function(x){plaplace(x,1,1)},
                function(n){rnorm(n,0,1)},0.24170,"When H0: X ~ Normal(0,1) and F1 = Laplace(1,1),n = 30")
g2 <- Confidence_band_H1(30,seq(-5,5,by = 0.1),function(x){plogis(x,0,1)},
                function(n){rcauchy(n,0,1)},0.24170,"When H0: X ~ Cauchy(0,1) and F1 = Logistic(0,1),n = 30")
g3 <- Confidence_band_H1(30,seq(-5,5,by = 0.1),function(x){plogis(x,1,1)},
                function(n){rnorm(n,0,1)},0.24170,"When H0: X ~ Normal(0,1) and F1 = Logistic(1,1),n = 30")
g4 <- Confidence_band_H1(30,seq(0,5,by = 0.1),function(x){pexp(x,1)},
                   function(n){rweibull(n,4,1)},0.24170,"When H0: X ~ Weibull(4,1) and F1 = Exp(1),n = 30")

gridExtra::grid.arrange(g1,g2,g3,g4,ncol = 2)

#===========================================

Coverage <- function(n,x,H1_CDF,rCDF_H0,D_na){
  set.seed(12345)
  s<-replicate(1000,{
       my.s <- rCDF_H0(n) 
      empircial_cdf <- ecdf(my.s)
      Ln.x <- vapply(empircial_cdf(x) - D_na,FUN = function(z){max(0,z)},FUN.VALUE = 2)
      Un.x <- vapply(empircial_cdf(x) + D_na,FUN = function(z){min(1,z)},FUN.VALUE = 2)
      (all(Ln.x < H1_CDF(x) & Un.x > H1_CDF(x))) 
   })
  mean(s)
}

Coverage.matrix <- matrix(
c(Coverage(20,seq(-5,5,by = 0.1),function(x){plaplace(x,0,1)},
           function(n){rlaplace(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pnorm(x,0,1)},
         function(n){rlaplace(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){plogis(x,0,1)},
         function(n){rlaplace(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pcauchy(x,0,1)},
         function(n){rlaplace(n,0,1)},0.29407),

Coverage(20,seq(-5,5,by = 0.1),function(x){plaplace(x,0,1)},
         function(n){rnorm(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pnorm(x,0,1)},
         function(n){rnorm(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){plogis(x,0,1)},
         function(n){rnorm(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pcauchy(x,0,1)},
         function(n){rlaplace(n,0,1)},0.29407),

Coverage(20,seq(-5,5,by = 0.1),function(x){pnorm(x,0,1)},
         function(n){rlogis(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pnorm(x,0,1)},
         function(n){rlogis(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){plogis(x,0,1)},
         function(n){rlogis(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pcauchy(x,0,1)},
         function(n){rlaplace(n,0,1)},0.29407),

Coverage(20,seq(-5,5,by = 0.1),function(x){pnorm(x,0,1)},
         function(n){rcauchy(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pnorm(x,0,1)},
         function(n){rcauchy(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){plogis(x,0,1)},
         function(n){rcauchy(n,0,1)},0.29407),
Coverage(20,seq(-5,5,by = 0.1),function(x){pcauchy(x,0,1)},
         function(n){rcauchy(n,0,1)},0.29407)),
ncol = 4,byrow = T)


colnames(Coverage.matrix) = rownames(Coverage.matrix) = c("Laplace(0,1)","Normal(0,1)","Logis(0,1)","Cauchy(0,1)")
as.data.frame(Coverage.matrix)


#===============================================================================

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

#Data is generated from Exp(1)
n <- 20
exp.ks <- replicate(3000,{
  ini.s <- rexp(n,1)
  my.s <- cummin(ini.s)*(1:n)
  sqrt(n)*ks.test(my.s,"pexp",alternative = "two.sided")$statistic
})

g1 <- ggplot(g.data <- data.frame(test_stat = exp.ks),aes(x = exp.ks)) +
  geom_histogram(aes(y = ..density..),col = "black",fill = "seagreen2",breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
  by = 1/sqrt(100)),position = "identity") +
  labs(x = 'sqrt(n)*Dn',y = 'Density',title = "Distribution of Kolmogorov-Smirnov test",
       subtitle = "When Data is dependent Exponential(1) and n = 20") + 
  geom_function(fun = function(x){vapply(x, FUN = emp_pdf, FUN.VALUE = 2)},col = "red",size = 0.9) + theme_bw(14) + 
  defined_theme

n <- 30
exp.ks <- replicate(3000,{
  ini.s <- rexp(n,1)
  my.s <- cummin(ini.s)*(1:n)
  sqrt(n)*ks.test(my.s,"pexp",alternative = "two.sided")$statistic
})

g2 <- ggplot(g.data <- data.frame(test_stat = exp.ks),aes(x = exp.ks)) +
  geom_histogram(aes(y = ..density..),col = "black",fill = "aquamarine",breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
                                                                      by = 1/sqrt(100)),position = "identity") +
  labs(x = 'sqrt(n)*Dn',y = 'Density',title = "Distribution of Kolmogorov-Smirnov test",
       subtitle = "When Data is dependent Exponential(1) and n = 30") + 
  geom_function(fun = function(x){vapply(x, FUN = emp_pdf, FUN.VALUE = 2)},col = "red",size = 0.9) + theme_bw(14) + 
  defined_theme

gridExtra::grid.arrange(g1,g2,ncol = 2)

#Data is generated from N(0,1)
n <- 20
exp.norm <- replicate(3000,{
  ini.s <- rnorm(n,1)
  my.s <- cumsum(ini.s)/sqrt(1:n)
  sqrt(n)*ks.test(my.s,"pnorm",alternative = "two.sided")$statistic
})

g1 <- ggplot(g.data <- data.frame(test_stat = exp.norm),aes(x = exp.norm)) +
  geom_histogram(aes(y = ..density..),col = "black",fill = "hotpink",breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
                                                                      by = 1/sqrt(100)),position = "identity") +
  labs(x = 'sqrt(n)*Dn',y = 'Density',title = "Distribution of Kolmogorov-Smirnov test",
       subtitle = "When Data is dependent Normal(0,1) and n = 20") + 
  geom_function(fun = function(x){vapply(x, FUN = emp_pdf, FUN.VALUE = 2)},col = "red",size = 0.9) + theme_bw(14) + 
  defined_theme

n <- 30
exp.norm <- replicate(3000,{
  ini.s <- rnorm(n,1)
  my.s <- cumsum(ini.s)/sqrt(1:n)
  sqrt(n)*ks.test(my.s,"pnorm",alternative = "two.sided")$statistic
})

g2 <- ggplot(g.data <- data.frame(test_stat = exp.norm),aes(x = exp.norm)) +
  geom_histogram(aes(y = ..density..),col = "black",fill = "coral1",breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
                                                                      by = 1/sqrt(100)),position = "identity") +
  labs(x = 'sqrt(n)*Dn',y = 'Density',title = "Distribution of Kolmogorov-Smirnov test",
       subtitle = "When Data is dependent Normal(0,1) and n = 30") + 
  geom_function(fun = function(x){vapply(x, FUN = emp_pdf, FUN.VALUE = 2)},col = "red",size = 0.9) + theme_bw(14) + 
  defined_theme

gridExtra::grid.arrange(g1,g2,ncol = 2)

#Data is generated from C(0,1)
n <- 20
exp.c <- replicate(3000,{
  ini.s <- rcauchy(n,1)
  my.s <- cumsum(ini.s)/(1:n)
  sqrt(n)*ks.test(my.s,"pcauchy",alternative = "two.sided")$statistic
})

g1 <- ggplot(g.data <- data.frame(test_stat = exp.c),aes(x = exp.c)) +
  geom_histogram(aes(y = ..density..),col = "black",fill = "seagreen2",breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
                                                                      by = 1/sqrt(100)),position = "identity") +
  labs(x = 'sqrt(n)*Dn',y = 'Density',title = "Distribution of Kolmogorov-Smirnov test",
       subtitle = "When Data is dependent Cauchy(0,1) and n = 20") + 
  geom_function(fun = function(x){vapply(x, FUN = emp_pdf, FUN.VALUE = 2)},col = "red",size = 0.9) + theme_bw(14) + 
  defined_theme

n <- 30
exp.c <- replicate(3000,{
  ini.s <- rcauchy(n,1)
  my.s <- cumsum(ini.s)/(1:n)
  sqrt(n)*ks.test(my.s,"pcauchy",alternative = "two.sided")$statistic
})

g2 <- ggplot(g.data <- data.frame(test_stat = exp.c),aes(x = exp.c)) +
  geom_histogram(aes(y = ..density..),fill = "aquamarine",col = "black",breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
                                                                      by = 1/sqrt(100)),position = "identity") +
  labs(x = 'sqrt(n)*Dn',y = 'Density',title = "Distribution of Kolmogorov-Smirnov test",
       subtitle = "When Data is dependent Cauchy(0,1) and n = 30") + 
  geom_function(fun = function(x){vapply(x, FUN = emp_pdf, FUN.VALUE = 2)},col = "red",size = 0.9) + theme_bw(14) + 
  defined_theme

gridExtra::grid.arrange(g1,g2,ncol = 2)

n <- 10000
ini.s <- rexp(n,1)
my.s <- cummin(ini.s)*(1:n)
my.test_stat <- NULL
for(i in 1:n){
 my.test_stat <- c(my.test_stat,ks.test(my.s[1:i],"pexp",alternative = "two.sided")$statistic)  
}

ggplot(g.data <- data.frame(Dn = my.test_stat,n = 1:n),aes(x = n,y = my.test_stat))+
  geom_line(col = "red") + labs(x = 'Dn',y = '',title = "Convergence of Kolmogorov-Smirnov test",
                               subtitle = "When Data is generated from dependent Exp(1) ") + 
  theme_bw(14) + 
  defined_theme

n <- 10000
ini.s <- rnorm(n,0,1)
my.s <- cummin(ini.s)/sqrt(1:n)
my.test_stat <- NULL
for(i in 1:n){
  my.test_stat <- c(my.test_stat,ks.test(my.s[1:i],"pnorm",alternative = "two.sided")$statistic)  
}

ggplot(g.data <- data.frame(Dn = my.test_stat,n = 1:n),aes(x = n,y = my.test_stat))+
  geom_line(col = "red") + labs(x = 'Dn',y = '',title = "Convergence of Kolmogorov-Smirnov test",
                                subtitle = "When Data is generated from dependent Normal(0,1) ") + 
  theme_bw(14) + 
  defined_theme

#=====================================

rm(list = ls(all = T))
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))











