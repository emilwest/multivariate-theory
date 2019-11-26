library(tidyverse)
source("utility_functions.R")
#ch10

#example 10.1
Z <- matrix( c(1,.4,.5,.6,
               .4,1,.3,.4,
               .5,.3,1,.2,
               .6,.4,.2,1), nrow=4,ncol=4 )
Z
# PARTITION 
p<-2
q<-2
n<-q+p
P11 <- Z[1:p,1:q] 
P12 <- Z[1:p,(q+1):n]
P21 <- t(P12)
P22 <- Z[(q+1):n,(q+1):n]

P11inv<- get_square_root(P11) %>% solve()
P22inv<- P22 %>% solve()

A <- P11inv%*%P12%*%P22inv%*%P21%*%P11inv
eigen(A)
e1<-eigen(A)$vectors[,1]*(-1)
e1 <- e1 %>% as.matrix()

a1 <- P11inv%*%e1
b1 <- P22inv%*%P21%*%a1
varb1 <- t(b1)%*%P22%*%b1
sdd<- (1/sqrt(varb1)) %>% as.vector() #make it vector to behave like scalar
#now scale b1 by 1/sqrt(variance of b1) = standard deviation of b1 to get b1 to have variance=1
b1<- sdd*b1

# the first pair of canonical variates are:
a1 # 0.86*Z1 + 0.28*Z2 where Z1 and Z2 in first partirion 
b1 # 0.54*Z1 + 0.74*Z2 where Z1 and Z2 in second partirion 
# and their canon corr is 
sqrt(eigen(A)$values[1])

# onoly works for 2 partitons currently
partition_2d <- function(mat, p, q){
    n <- p+q
    S11 <- mat[1:p,1:p] 
    S12 <- mat[1:p,(p+1):n]
    S21 <- t(S12)
    S22 <- mat[(p+1):n,(p+1):n]
    return(list(
        mat = mat,
        S11 = S11,
        S12 = S12,
        S21 = S21,
        S22 = S22
    ))
    
}



# exercise 10.2
S <- matrix( c(8,2,3,1,
               2,5,-1,3,
               3,-1,6,-2,
               1,3,-2,7), nrow=4,ncol=4 )
S
p<-2
q<-2
n<-q+p
S11 <- S[1:p,1:q] 
S12 <- S[1:p,(q+1):n]
S21 <- t(S12)
S22 <- S[(q+1):n,(q+1):n]

S11inv <- get_square_root(S11) %>% solve()
S22inv <- S22 %>% solve()

lambda <- sqrt(eigen(S11)$values) 
lambda <-1/lambda
lambda<- lambda %>% diag()

eigen(S11)$vectors%*%lambda%*%t(eigen(S11)$vectors )
0.000409

A <- S11inv%*%S12%*%S22inv%*%S21%*%S11inv
A
eigen(A)
sqrt(.3)
sqrt(.239)

e1<-eigen(A)$vectors[,1] %>% as.matrix()
e2<-eigen(A)$vectors[,2] %>% as.matrix()
e1%*%t(e2)
a1 <- S11inv%*%e1
a2 <- S11inv%*%e2
b1 <- S22inv%*%S21%*%a1
b2 <- S22inv%*%S21%*%a2


varb1 <- t(b1)%*%S22%*%b1
sdd<- (1/sqrt(varb1)) %>% as.vector() #make it vector to behave like scalar
varb2 <- t(b2)%*%S22%*%b2
sdd2<- (1/sqrt(varb2)) %>% as.vector() #make it vector to behave like scalar

#now scale b1 by 1/sqrt(variance of b1) = standard deviation of b1 to get b1 to have variance=1
b1<- sdd*b1
b2<- sdd2*b2
a1
a2
b1
b2
det(ones(1))
ones(1)*2 %>% solve()
ones(1)*5 %>% solve()



#install.packages("CCA")
library(CCA)


# 10.9

R <- matrix(c(1,.6328,.2412,.0586,
              .6328,1,-.0553,.0655,
              .2412,-.0553,1,.4248,
              .0586,.0655,.4248,1), nrow=4,ncol=4)
R
#::cc(X = R[,1:2],Y=R[,3:4])

partition_2d(R, 3, 1)


my_CCA <- function(covmat, X, Y, population=T){
    partitions <- partition_2d(covmat, X, Y)
    S <- partitions$mat
    S11 <- partitions$S11
    S11_sqrt_inv <- get_square_root(partitions$S11) %>% solve()
    S22_sqrt_inv <- get_square_root(partitions$S22) %>% solve()
    S11_inv <- partitions$S11 %>% solve()
    S22_inv <- partitions$S22 %>% solve()
    S22 <- partitions$S22
    S12 <- partitions$S12
    S21 <- partitions$S21
    
    
    A <- S11_sqrt_inv%*%S12%*%S22_inv%*%S21%*%S11_sqrt_inv
    B <- S22_sqrt_inv%*%S21%*%S11_inv%*%S12%*%S22_sqrt_inv
    
    e <- eigen(A)$vectors*(-1)
    f <- eigen(B)$vectors*(-1)
    
    a <- matrix(NA, ncol=ncol(e),nrow = nrow(e))
    b <- matrix(NA, ncol=ncol(e),nrow = nrow(e))
    
    for (i in 1:ncol(e)){
        #print(i)
        a[,i] <- as.matrix(S11_sqrt_inv%*%e[,i])
        if (population==F){
            b[,i] <- S22_sqrt_inv%*%f[,i] %>% as.matrix()
        } else {
            b_temp <- S22inv%*%S21%*%a[,i] %>% as.matrix
            varb <- t(b_temp)%*%S22%*%b_temp
            sdd<- (1/sqrt(varb)) %>% as.vector()
            b[,i] <- sdd*b_temp
        }
    }
    
    return(
        list(
            U = t(a),
            V = t(b),
            canonical_correlation = sqrt(eigen(A)$values),

        )
    )
}

# p= num variables in partition 1, q=num variables in partition 2
test_if_partitions_are_uncorrelated <- function(alpha=0.05, covmat, p, q, n){
    partitions <- partition_2d(covmat, p, q)
    S11 <- partitions$S11
    S22 <- partitions$S22
    S <- partitions$mat
    
    # H0: covariance S12=0 ie partition 1 is uncorrelated with partition 2
    # H1: covariance S12!=0 ie they are correlated 
    LRT_statistic <- -n*log(( det(S11)*det(S22) )/ det(S))
    critical <- qchisq(alpha, df = (p*q),lower.tail = F)
    p <- pchisq(critical, df = (p*q),lower.tail = F)
    Barletts_corrected <- -(n-1-0.5*(p+q+1))*log(( det(S11)*det(S22) )/ det(S))
    uncorrelated <- Barletts_corrected > critical
    
    return(
        list(
            LRT_statistic = LRT_statistic,
            Barletts_corrected = Barletts_corrected,
            critical = critical,
            p = p,
            uncorrelated = uncorrelated
        )
    )
}

test_if_partitions_are_uncorrelated(alpha=0.05,covmat = R,p = 2,q = 2, n=140)

# sequential_test <- function(alpha=0.05, canonical_correlations, p, q, n){
#     #H0: canonical correlations 1,..,k are not zero and k+1,...,p = 0
#     #H1: i:th canonical correlation is not zero for some i >= k+1
#     for (i in 1:p){
#         i <- i+1
#         
#         
#     }
#     Barletts_corrected <- -(n-1-0.5*(p+q+1))*log(prod(1-canonical_correlations[-i]^2))
#     return(Barletts_corrected)    
# }
# 
# sequential_test(alpha = 0.05, canonical_correlations = my_CCA(R,2,2,population = F,n=140)$canonical,
#                 p=2,q=2,n=140)


my_CCA(S,2,2)
my_CCA(R,2,2,population = F,n=140)

#install.packages("CCP")
library(CCP)
rho <- my_CCA(R,2,2,population = F)$canonical
p.asym(rho = rho,N = 140,p = 2,q = 2) # SEQUENTIAL TEST
# 