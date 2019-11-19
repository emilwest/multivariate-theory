library(tidyverse)
source("utility_functions.R")

lambda1<-1.96
lambda2<-0.68
lambda3<-0.36
e1<-c(.625, .593, .507)
e2<-c(-.219,-.491,.843)
e3<-c(.749, -.638, -.177)


R<-matrix(c(.81,.63,.45,
            .63,.49,.35,
            .45,.35,.25), ncol=3,nrow=3)
R
eigen(R)


###############
# ex 9.19

df<-read_table2(file="data/T9-12.txt", col_names = F)
df
n<-nrow(df)
# assume orthogonal factor model

# a) standardize variables and obtain either PCA or ML solutions 
# for m=2 and m=3
X <- df %>% as.matrix()
S <- my_covmatrix(X)
R <- convert_covmatrix_to_corrmatrix(S)
S_diag <- S %>% diag() %>% diag() #get variances as diagonal matrix
D <- S_diag %>% solve() %>% sqrt() #get one over standard deviation 
xbar<-colMeans(df) %>% as.matrix()
xbar_np<-ones(n)%*%t(xbar) #creates nxp matrix of means
dim(X)

# STANDARDIZED X
Z <- (X-xbar_np)%*%D

# use R since S has some variables with unproportionately large varances
ei <- eigen(R)
lambdas <- ei$values
e <- ei$vectors

# m= number of common factors
m<-2
m<-3
#diag(lambdas) makes the scalar be multiplied to each column
L<- e[,1:m]%*%diag(sqrt(lambdas[1:m])) 
L #loading matrix
LL <- L%*%t(L) 
LL
communalities<- rowSums(L^2)
psi<- (1 - communalities) %>% diag()
psi 
LL+psi # fitted correlation matrix
R
residual_matrix <- R-(LL+psi)
residual_matrix %>% round(3)
# residual matrix for m=3 is smaller, 
# ie 3 factors represents the observations better 

###
# using maximum likelihood estimation

#?factanal()
# DEFAULT: rotated
# varimax is default
mle2 <- factanal(Z, factors = 2)
mle3 <- factanal(Z, factors = 3,rotation="varimax")

mle2$loadings
#x2 (sales profitability) and x7 (math ability) are high for factor 1,
#x4 (creativity) and x6 (abstract reasoning) for factor 2 and 3
mle3$loadings 
mle2$uniquenesses #psi values
mle3$uniquenesses #psi values

LL <- mle2$loadings[1:7,]%*%t(mle2$loadings[1:7,])
fitted_mle2<- LL+diag(mle2$uniquenesses)
residual_matrix_mle2 <- R-(fitted_mle2)

LL <- mle3$loadings[1:7,]%*%t(mle3$loadings[1:7,])
fitted_mle3<- LL+diag(mle3$uniquenesses)
residual_matrix_mle3 <- R-(fitted_mle3)

residual_matrix_mle2 %>% round(3)
residual_matrix_mle3 %>% round(3)
#mle with 3 factors has smaller values in residual matrix

#unrotated
factanal(Z, factors = 2,rotation = "none")
factanal(Z, factors = 3,rotation = "none",covmat = R,n.obs = n)


# d) test of model fit
# H0: good model fit
# H1: not good model fit
#maximum likelihood
alpha<-0.01
p<-7
m<-2
lambda <- det(fitted_mle2)/det(R)
#barletts correction
q <- (n-1-(2*p + 4*m + 5)/6)*log(lambda)
d <- ((p-m)^2 - p - m)/2
critial <- qchisq(p=alpha,df=d,lower.tail =F) 
q > critial
# Reject H0

m<-3
lambda <- det(fitted_mle3)/det(R)
#barletts correction
q <- (n-1-(2*p + 4*m + 5)/6)*log(lambda)
d <- ((p-m)^2 - p - m)/2
critial <- qchisq(p=alpha,df=d,lower.tail =F) 
q > critial
# Reject H0

# by the chisaure criterion, niether m=2 or m=3 models appear to fit








#############################################################
# ex 9.31

df <- read_table2(file = "data/T8-6.txt",col_names = F)
df
colnames(df) <- c(100,200,400,800,1500,5000,10000,"marathon", "country")
colnames(df)

df_meters_per_second <- df %>% 
    select(everything()) %>% 
    mutate(hundrams = sort(100/`100`) ,
           twohndrams = sort(200/`200`),
           fourhundams = sort(400/`400`),
           eightms = sort(800/(`800`*60)),
           femtonms = sort(1500/(`1500`*60)),
           femms = sort(5000/(`5000`*60)),
           tenms = sort(10000/(`10000`*60)),
           marathonms = sort(42195/(marathon*60))
    ) %>%
    select( hundrams, twohndrams, fourhundams, 
            eightms, femtonms, femms, tenms, marathonms)

df_meters_per_second
S<-my_covmatrix( as.matrix(df_meters_per_second) )
R<-convert_covmatrix_to_corrmatrix(S)
S
R
eigen(S)
X <- df_meters_per_second %>% as.matrix()
X
factanal(X, factors = 1, covmat = S)
factanal(X, factors = 2, covmat = S)
# factor 1: higher loadings on longer distnace, may represent endurance
# factor 2: higher loadings on shorter distances, may repr. speed
plot(X)

# check outliers
library(reshape2)
melt(df_meters_per_second) %>% 
    ggplot(aes(x=variable,y=value)) +
    geom_point(mapping = aes(color=variable)) 




#############################################################
# ex 9.33
df <- read_table2(file = "data/T4-6.txt",col_names = F)
df
colnames(df) <-c("lndep" ,"Supp", "Benev",
                 "Conform" ,"Leader", "Gender", "Socio")
df
X <- df[,1:5] %>% as.matrix()
X
S<-my_covmatrix(X)
R<-convert_covmatrix_to_corrmatrix(S)
R
princomp(R)
eigen(R)
factanal(X, factors = 2, covmat = R)
factanal(X, factors = 2, covmat = R,rotation = "none")
factanal(X, factors = 3, covmat = R) # error: 3 factors are too many for 5 variables
# probably the heywood case for m=3
scree_plot(eigenvalue = eigen(R)$values)


# Principal component solution of factor model.
# if rotated: varimax rotation is default
# varimax selects the orthogonal transformation T 
# that maximizes V which contains scaled rotated coefficients.
# An orth. rotation simply rotates the coordinate axes.
# maximizing V corresponds to "spreading out" 
# the squares of the loadings on each factor as much as possible.
# Rotating the factors is similar to sharpening the focus 
# on a microscope in order to see details more clearly. 
# It makes the factors more interpretable. 
PCF <- function( factors, covmat, rotate=T){
    ei <- eigen(covmat)
    lambdas <- ei$values
    n <- nrow(covmat)
    e <- ei$vectors
    L <- e[,1:factors]%*%diag(sqrt(lambdas[1:factors])) 
    if (rotate==T){
        rotated <- varimax(as.matrix(L)) 
        L <- rotated$loadings[1:n,] %>% as.matrix()
    }
    LL <- L%*%t(L) 
    communalities<- rowSums(L^2)
    psi <- (1 - communalities) %>% diag()
    fitted <- LL+psi
    residual_matrix <- covmat-fitted
    variance_communality <- sum(communalities)
    variance_per_factor <- lambdas[1:factors]
    percent_variance_per_factor <- variance_per_factor/sum(lambdas)
    
    return(
        list(L=L, #loading matrix
             LL=LL,
             communalities = communalities,
             psi = psi,
             fitted =  fitted,
             residual_matrix = residual_matrix,
             variance_communality = variance_communality,
             variance_per_factor = variance_per_factor,
             percent_variance_per_factor = percent_variance_per_factor,
             if (rotate==T){rotated = rotated$loadings}
             )
    )
}

PCF(2,R)
PCF(3,R, rotate = F)
PCF(3,R, rotate = T)

