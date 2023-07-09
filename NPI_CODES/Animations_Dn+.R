rm(list = ls(all = T))  #removes all objects
library(ggplot2)
library(VGAM)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 11,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 10),plot.title = element_text(family = "serif",
                                                                                     colour = "red", hjust = -0.01),legend.text = element_text(size = 10,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))

#===============================================================================


Asymptotic_distribution_H0 <- function(n,True_CDF,r_CDF,parameter,tag = " "){
  set.seed(seed = 1234)
  Distribution_data <- NULL
  for(i in 1:length(parameter)){
   asymp_test_stat_val <- replicate(10000,{
     our_sample <- r_CDF(n) + parameter[i]
     performed_test <- ks.test(our_sample,True_CDF,alternative = "greater")
     sqrt(n)*performed_test$statistic
   })
   Distribution_data <- rbind(Distribution_data,
                        data.frame(Test_stat = asymp_test_stat_val,theta = parameter[i]))
  }
  
  graph.1 <- ggplot(Distribution_data,aes(Test_stat,frame = theta)) + 
    geom_histogram(aes(y = ..density..),position = "identity",breaks = seq(0,ceiling(sqrt(n)) + 1/sqrt(100),
    by = 1/sqrt(100)),col = 'black',fill = "cyan") + labs(x = 'sqrt(n)*Dn+',y = 'Density',
    title = "Asymptotic Distribution Under H0 of Sample Kolmogorov-Smirnov Test Statistic Dn+",
    subtitle = tag) + theme_bw(14) + defined_theme
  plotly::ggplotly(graph.1)
}  
