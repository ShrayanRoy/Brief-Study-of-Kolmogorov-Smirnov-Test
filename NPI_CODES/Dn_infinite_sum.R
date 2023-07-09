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
emp_pdf(1)
z=seq(0,5,by=0.01)
pdf=vapply(z,FUN=emp_pdf,FUN.VALUE = 2)
plot(z,pdf,type="l")
