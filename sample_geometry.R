library(tidyverse)
source("utility_functions.R")
#ch3

#ex3.1
X<-matrix(c(9,5,1,1,3,2),ncol=2,nrow = 3,byrow = F)

X %>% as_tibble() %>% ggplot(aes(x=V1,y=V2)) + 
    geom_point() +
    geom_point(aes(x=mean(V1),y=mean(V2)),color="red")  #mean coordinate

# p=2 vector in n=3 space
library(plotly)
plot_ly(data = as.data.frame(t(X)), x=~V1,y=~V2, z=~V3, type="scatter3d") 
?plot_ly






# 3.9

X<-matrix(c(12,18,14,20,16,17,20,16,18,19,29,38,30,38,35),ncol=3,nrow = 5,byrow = F)
X
xbar <- colMeans(X)
xbar
xbar<-ones(nrow(X))%*%xbar #as n by p matrix
Z <- X-xbar #mean corrected
Z
#notice that col 3 is col1+col2
# therefore col1+col2 - col3 = 0
# ie columns are linearly dependent
a<-matrix(c(1,1,-1), nrow = 3,ncol=1)
a
Z%*%a # = 0 

S <- my_covmatrix(Z)
det(S) # the generalized sample variance is basically zero 
S%*%a # = 0 , ie a can be rescaled to be an eigenvector corresponding to eigenvalue zero


# 3.14
X<-matrix(c(9,5,1,1,3,2),ncol=2,nrow = 3,byrow = F)
c <- c(-1,2) %>% as.matrix()
b <- c(2,3) %>% as.matrix()
b
X
S<-my_covmatrix(X)
X
t(b)%*%X[1,]
t(b)%*%X[2,]
t(b)%*%X[3,]

xbar<- colMeans(X) %>% as.matrix()
t(b)%*%xbar
t(c)%*%xbar
t(b)%*%S%*%b
t(c)%*%S%*%c
t(b)%*%S%*%c


#  3.15
X<-matrix(c(1,6,8,4,2,3,3,6,3),ncol=3,nrow = 3,byrow = F)
n<-nrow(X)
xbar <- colMeans(X)
XBAR<-ones(nrow(X))%*%xbar #as n by p matrix
I<- diag(n)
J <- ones(n)%*%t(ones(n))
S <- (1/(n-1))*t(X)%*%(I-(1/n)*J)%*%X # or:
S <- (1/(n-1))*t(X)%*%(X-XBAR)
my_covmatrix(X)

#linear combinartions
b <- c(1,1,1) %>% as.matrix()
c <- c(1,2,-3) %>% as.matrix()

t(b)%*%xbar #sample mean of b'X
t(c)%*%xbar  #sample mean of c'X
t(b)%*%S%*%b #sample variance of b'X
t(c)%*%S%*%c #sample variance of c'X
t(b)%*%S%*%c  #sample covariance of b'X an c'X




cc <- matrix( c(6, -10/4, -10/4, 1.5),nrow=2,ncol=2)
cc <- matrix( c(48.68,-26,-26,14),nrow=2,ncol=2)
cc<-cc/2
solve(cc)
det(cc)
4/27


Y<-matrix(c(-3,4,5,15,16,11),nrow=3,ncol=2, byrow = F)
Y
cov(Y)
xbar<-colMeans(Y)
xbar<-ones(3)%*%t(xbar)
S<-(1/2)*t(Y)%*%(Y-xbar)
solve(cov(Y))
i<-solve(S)
my_covmatrix(Y)
mu <- c(4,14) %>% as.matrix()
ones(n)%*%t(mu)
xbar <- as.matrix(xbar)
xbar
n<-3
n*(t(xbar-mu)%*%i%*%(xbar-mu))
7/9
