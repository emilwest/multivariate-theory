m<-matrix(c(7476.45,303.62,303.62,26.19),ncol=2,nrow=2)
m
m<-matrix(c(1,0.2015,0.2015,1),ncol=2,nrow=2)
m
eigen(m)


m<-matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),ncol=3,nrow=3)
m
eigen(m)
1/sqrt(3)


m<-matrix(c(1,0.75,0.63,0.64,
            0.75,1,0.69,0.74,
            0.63,0.69,1,0.66,
            0.64,0.74,0.66,1),
          ncol=4,nrow=4
)


################


# ex. 8.13
library(tidyverse)

df<-read_table2(file="data/T1-7.txt", col_names = F)
df
cor(df) #R
cov(df) #S
eigen(cor(df))
df

my_covmatrix <- function(mat){
    X <- as.matrix(mat) # n by p matrix
    n <- nrow(X)
    dimensions <- dim(X)
    I <- diag(dimensions[1])
    ones <- matrix(1L, nrow = dimensions[1], ncol = 1)
    J <- ones%*%t(ones)
    
    # p by n * ... * n by p = p by p matrix 
    S <- (1/(n-1))*t(X)%*%(I-(1/n)*J)%*%X
    return(S)
}
convert_covmatrix_to_corrmatrix <- function(S){
    # diag diag makes it a diagonal matrix
    # solve gets inverse
    S_inv <- S %>% diag() %>% diag() %>% solve() %>% sqrt()
    R <- S_inv%*%S%*%S_inv
    return(R)
}
convert_corrmatrix_to_covmatrix <- function(R,S){
    # diag diag makes it a diagonal matrix
    D <- S %>% diag() %>% diag() %>%  sqrt()
    S <- D%*%R%*%D
    return(S)
}
ones <- function(int){
    return(matrix(1L,nrow = int, ncol=1))
}
get_quantiles <- function(df){
    i = 1:nrow(df)
    p = (i-0.5)/nrow(df)
    q = qnorm(p)
    return(q)
}


S <- my_covmatrix(df)
R <- convert_covmatrix_to_corrmatrix(S)
convert_corrmatrix_to_covmatrix(R,S)

#notice x1 has higer variance than other x
#therefore its more stable to use the corr matrix R
diag(S)

R

eigen(R)

eigens<-eigen(R)$values %>% enframe() 

eigens %>% select(value) %>% mutate(varexplained = value/n(),
                                    cumvarexplained = cumsum(value)/n() )

# 4 variables needed to summarize 
neweigens <- eigens %>% filter(name<=4) %>% select(value) %>% as_vector()

neww<- neweigens %>% diag() %>% sqrt()
# correrlation between x:s and first four principal components
eigen(R)$vectors[,1:4]%*%neww




###################################
# ex 8.19

df <- read_table2(file = "data/T1-9.txt",col_names = F)
df

colnames(df) <- c(100,200,400,800,1500,3000,"marathon", "country")
colnames(df)

df_meters_per_second <- df %>% 
    select(everything()) %>% 
    mutate(hundrams = sort(100/`100`) ,
           twohndrams = sort(200/`200`),
           fourhundams = sort(400/`400`),
           eightms = sort(800/(`800`*60)),
           femtonms = sort(1500/(`1500`*60)),
           threems = sort(3000/(`3000`*60)),
           marathonms = sort(42195/(marathon*60))
    ) %>%
    select( hundrams, twohndrams, fourhundams, 
           eightms, femtonms, threems, marathonms)

S<-my_covmatrix( as.matrix(df_meters_per_second) )
S
eigen(S)
cov(df_meters_per_second)



#############
# ex 8.20 men
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
eigen(S)



#############################################################
# ex 8.22

df <- read_table2(file = "data/T1-10.txt",col_names = F)
df

colnames(df) <- c("Breed","SalePr", "YrHgt" ,"FtFrBody", "PrctFFB" ,"Frame" 
  ,"BkFat" ,"SaleHt", "SaleWt")
df

df <- df %>% select(-Breed,-SalePr)
df

S <- my_covmatrix(df)
R <- convert_covmatrix_to_corrmatrix(S)
R

#PCA with S:
eigens<-eigen(S)$values %>% enframe() 

eigens %>% select(value) %>% mutate(varexplained = value/n(),
                                    cumvarexplained = cumsum(value)/n(),
                                    propexplained = value/sum(value),
                                    cumprop = cumsum(value)/sum(value))


#SCREE PLOT!!!!!!!!!
eigens %>% ggplot(aes(y=value,x=name)) +
    geom_point(fill="black", colour="white", shape=21, size=2.5) +
    geom_path()


# 1 variable explan 80% of the variation , 1+2 explain almost 100%
neweigens <- eigens %>% filter(name<=2) %>% select(value) %>% as_vector()
#EIGENVECTORS:

eigen(S)$vectors %>% round(6)

eigen(S)$vectors %>% as_tibble() %>%
    ggplot(aes(x=V2,y=V1)) +
    geom_point()


#p.460 construct sample principal components
colMeans(df)
xbar <- colMeans(df) %>% as.matrix()
e1 <-eigen(S)$vectors[,1]  %>% as.matrix()
e1<-e1*(-1)
e2 = eigen(S)$vectors[,2]  %>% as.matrix()
X <- df %>% as.matrix()
X

one<-ones(nrow(df))
xbar_np<-one%*%t(xbar) #creates nxp matrix of means

(X-xbar_np)%*%e1
(X-xbar_np)%*%e2

df2 <- df %>% mutate(
    #y1 = (X-xbar_np)%*%e1,
    #y2 = (X-xbar_np)%*%e2,
    y1 = (X)%*%e1,
    y2 = (X)%*%e2,
    q  = get_quantiles(X)
) %>% select(y1,y2,q)

df2
df2 %>% ggplot(aes(y=y1,x=y2)) +
                   geom_point() 

### QQPLOT
qqplot(y = df2$y1, x = get_quantiles(df) )
#with ggplot
df2 %>%
    ggplot(aes(x=q,y=sort(y1))) +
    geom_point()



qnorm(df$YrHgt)

#PCA with R:
eigens<-eigen(R)$values %>% enframe() 

eigens %>% select(value) %>% mutate(varexplained = value/n(),
                                    cumvarexplained = cumsum(value)/n(),
                                    propexplained = value/sum(value),
                                    cumprop = cumsum(value)/sum(value))

#SCREE PLOT!!!!!!!!!
eigens %>% ggplot(aes(y=value,x=name)) +
    geom_point(fill="black", colour="white", shape=21, size=2.5) +
    geom_path()


# 3 variables explain almost 90% of variiation
neweigens <- eigens %>% filter(name<=2) %>% select(value) %>% as_vector()
#EIGENVECTORS:
eigen(R)$vectors %>% round(6) 


e1 <-eigen(R)$vectors[,1]  %>% as.matrix()
e1 <- e1*(-1)
e2 = eigen(R)$vectors[,2]  %>% as.matrix()
df2 <- df %>% mutate(
    #y1 = (X-xbar_np)%*%e1,
    #y2 = (X-xbar_np)%*%e2,
    y1 = (X)%*%e1,
    y2 = (X)%*%e2,
    q  = get_quantiles(X)
) %>% select(y1,y2,q)

df2 %>% ggplot(aes(y=y1,x=y2)) +
    geom_point() 

qqplot(y = df2$y1, x = get_quantiles(df))

#with ggplot
df2 %>%
    ggplot(aes(x=q,y=sort(y1))) +
    geom_point()







#######################################################
# ex 8.27
df <- read_table2(file = "data/T7-7.txt",col_names = T)
df
colnames(df) <- c("BL","EM", "SF" ,"BS", "AFL" ,"LFF" 
                  ,"FFF","ZST")
df

df <- df %>% select(BL,EM,SF,BS)

is.na(df)
df<-na.omit(df)
S <- my_covmatrix(df)
R <- convert_covmatrix_to_corrmatrix(S)
R
S
#PCA with S:
eigens<-eigen(S)$values %>% enframe() 

eigens %>% select(value) %>% mutate(varexplained = value/n(),
                                    cumvarexplained = cumsum(value)/n(),
                                    propexplained = value/sum(value),
                                    cumprop = cumsum(value)/sum(value))


#SCREE PLOT!!!!!!!!!
eigens %>% ggplot(aes(y=value,x=name)) +
    geom_point(fill="black", colour="white", shape=21, size=2.5) +
    geom_path()

# 1 comp summarize all variation
#EIGENVECTORS:
eigen(S)$vectors %>% round(6)


#p.460 construct sample principal components
colMeans(df)
xbar <- colMeans(df) %>% as.matrix()
e1 <-eigen(S)$vectors[,1]  %>% as.matrix()
e1<-e1
e2 = eigen(S)$vectors[,2]  %>% as.matrix()
e2 = e2*(-1)
X <- df %>% as.matrix()
X

one<-ones(nrow(df))
xbar_np<-one%*%t(xbar) #creates nxp matrix of means

(X-xbar_np)%*%e1
(X-xbar_np)%*%e2

df2 <- df %>% mutate(
    #y1 = (X-xbar_np)%*%e1,
    #y2 = (X-xbar_np)%*%e2,
    y1 = (X)%*%e1,
    y2 = (X)%*%e2,
    q  = get_quantiles(X)
) %>% select(y1,y2,q)

df2
df2 %>% ggplot(aes(y=y1,x=y2)) +
    geom_point() 

### QQPLOT
qqplot(y = df2$y1, x = get_quantiles(df) )
#with ggplot
df2 %>%
    ggplot(aes(x=q,y=sort(y1))) +
    geom_point()





#############
# FOR R

eigens<-eigen(R)$values %>% enframe() 

eigens %>% select(value) %>% mutate(varexplained = value/n(),
                                    cumvarexplained = cumsum(value)/n(),
                                    propexplained = value/sum(value),
                                    cumprop = cumsum(value)/sum(value))


#SCREE PLOT!!!!!!!!!
eigens %>% ggplot(aes(y=value,x=name)) +
    geom_point(fill="black", colour="white", shape=21, size=2.5) +
    geom_path()

# 1 comp summarize all variation
#EIGENVECTORS:
eigen(R)$vectors %>% round(6)


#p.460 construct sample principal components
colMeans(df)
xbar <- colMeans(df) %>% as.matrix()
e1 <-eigen(R)$vectors[,1]  %>% as.matrix()
e1<-e1
e2 = eigen(R)$vectors[,2]  %>% as.matrix()
e2 = e2*(-1)
X <- df %>% as.matrix()
X

one<-ones(nrow(df))
xbar_np<-one%*%t(xbar) #creates nxp matrix of means

(X-xbar_np)%*%e1
(X-xbar_np)%*%e2

df2 <- df %>% mutate(
    y1 = (X-xbar_np)%*%e1,
    y2 = (X-xbar_np)%*%e2,
    #y1 = (X)%*%e1,
    #y2 = (X)%*%e2,
    q  = get_quantiles(X)
) %>% select(y1,y2,q)

df2
df2 %>% ggplot(aes(y=y1,x=y2)) +
    geom_point() 

### QQPLOT
qqplot(y = df2$y1, x = get_quantiles(df) )
#with ggplot
df2 %>%
    ggplot(aes(x=q,y=sort(y1))) +
    geom_point()






