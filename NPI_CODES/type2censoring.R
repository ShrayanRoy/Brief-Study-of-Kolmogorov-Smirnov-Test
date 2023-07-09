#=======================================================
rm(list = ls(all = T)) #removes all objects

#For type two censoring Data
Ks.Cens2Test <- function(data,H0_CDF,r){
  censored_data <- sort(data)[1:r]
  Z_censored_data <- H0_CDF(censored_data)
  test.statistic <- max(abs((((1:r) - 0.5)/length(data))  - Z_censored_data)) +  (0.5/length(data))
  return(test.statistic)
}

#=============================================

emp_pdf=function(x){
  pdf1 <-  function(z){
    epsilon = 0.01; a = 8*z*exp(-2*z*z);b = a + 2*epsilon;i = 2
    while(abs(b-a)>epsilon)
    {
      a = b;b = a+8*((-1)^(i-1))*i*i*z*exp(-2*i*i*z*z);i = i+1
    }
    return(b)
  }
  return(vapply(x,FUN = pdf1,FUN.VALUE = 2))
}


par(mfrow = c(2,2))

Ks.data_weibull <- replicate(10000,{
  sampled_data <- rweibull(30,scale = 2,shape = 1)
  Ks.Cens2Test(sampled_data,function(x){ pweibull(x,scale = 2,shape = 1) },10)
})
hist(sqrt(30)*Ks.data_weibull,col = "yellow",border = "red",freq = F,
     xlab = "Modified KS Statistic",
     main = "Asymptotic Distribution of Modified Dn for Type 2 Censoring and r = 10",
     breaks = seq(0,sqrt(30) + 1/sqrt(100),by = 1/sqrt(1000)),col.main = "#FF3333")
curve(emp_pdf,add = T,col = "blue")
mtext("Data from Weibull(2,1) Distribution",side = 3)
box()


Ks.data_exp <- replicate(10000,{
  sampled_data <- rexp(30,rate = 0.5)
  Ks.Cens2Test(sampled_data,function(x){ pexp(x,rate = 0.5) },10)
})
hist(sqrt(30)*Ks.data_exp,col = "yellow",border = "red",freq = F,
     xlab = "Modified KS Statistic",
     main = "Asymptotic Distribution of Modified Dn for Type 2 Censoring and r = 10",
     breaks = seq(0,sqrt(30) + 1/sqrt(100),by = 1/sqrt(1000)),col.main = "#FF3333")
curve(emp_pdf,add = T,col = "blue")
mtext("Data from Exp(1) Distribution",side = 3)
box()


Ks.data_exp1 <- replicate(10000,{
  sampled_data <- rweibull(30,2,1)
  Ks.Cens2Test(sampled_data,function(x){ pweibull(x,2,1)},20)
})
hist(sqrt(30)*Ks.data_exp1,col = "yellow",border = "red",freq = F,
     xlab = "Modified KS Statistic",
     main = "Asymptotic Distribution of Modified Dn for Type 2 Censoring and r = 20",
     breaks = seq(0,sqrt(30) + 1/sqrt(100),by = 1/sqrt(1000)),col.main = "#FF3333")
curve(emp_pdf,add = T,col = "blue")
mtext("Data from Weibull(2,1) Distribution",side = 3)
box()


Ks.data_exp1 <- replicate(10000,{
  sampled_data <- rexp(30,rate = 0.5)
  Ks.Cens2Test(sampled_data,function(x){ pexp(x,rate = 0.5) },20)
})
hist(sqrt(30)*Ks.data_exp1,col = "yellow",border = "red",freq = F,
     xlab = "Modified KS Statistic",
     main = "Asymptotic Distribution of Modified Dn for Type 2 Censoring and r = 20",
     breaks = seq(0,sqrt(30) + 1/sqrt(100),by = 1/sqrt(1000)),col.main = "#FF3333")
curve(emp_pdf,add = T,col = "blue")
mtext("Data from Exp(1) Distribution",side = 3)
box()


