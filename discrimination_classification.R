library(tidyverse)


## Ex 11.7
tibble(
    x1=seq(-1,1,length=100),
    x2=seq(-0.5,1.5,length=100),
    y1=1-abs(x),
    y2=(1-abs(x2-0.5))
) %>% ggplot() +
    geom_point(aes(x=x1,y=y1)) +
    geom_point(aes(x=x2,y=y2))

## Ex 11.8
tibble(
    x1=seq(-1,1,length=100),
    x2=seq(-1.5,2.5,length=100),
    y1=1-abs(x),
    y2=(0.25*(2-abs(x2-0.5)))
) %>% ggplot() +
    geom_point(aes(x=x1,y=y1)) +
    geom_point(aes(x=x2,y=y2))




spooled<-matrix(c(1,-1/3,-1/3,4),nrow=2,ncol=2,byrow = T)
# 
x<-matrix(c(.3571-.9556, .4667, .0714, .9000-.9556),nrow=2,ncol=2,byrow = T)
ref(x)
install.packages('matlib')
library('matlib')
matlib::echelon(x,verbose=T)
a <- c(0.7797828, 1) %>% as.matrix()
t(a)%*%spooled%*%a
b<-a/sqrt(4.088206)
t(b)%*%spooled%*%b


ybar <- matrix(c(-8.644,34,-5.89,30.44,-10.44,29.94),nrow=3,ncol=2,byrow = T )
ybar
yhat<-c(-5.865,32.025) %>% as.matrix()
yhat<-c(-5.865,32.025) %>% as.matrix()
yhat<-c(-5.865,32.025) %>% as.matrix()

source('utility_functions.R')
y<-ones(3)%*%t(yhat)
res<-y-ybar
res^2 %>% rowSums()
res%*%t(res)



a<- matrix(c(0.1239,-0.0076,-0.0799,0.0092),byrow=T,nrow=2,ncol=2)
a
eigen(a)
