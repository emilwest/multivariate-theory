
ones <- function(int){
    return(matrix(1L,nrow = int, ncol=1))
}
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
get_quantiles <- function(df){
    i = 1:nrow(df)
    p = (i-0.5)/nrow(df)
    q = qnorm(p)
    return(q)
}

scree_plot <- function(eigenvalue){
    eigens <- eigenvalue %>% enframe() 
    eigens %>% ggplot(aes(y=value,x=name)) +
        geom_point(fill="black", colour="white", shape=21, size=2.5) +
        geom_path() +
        labs(title = "Scree plot",
             x= "Component Number",
             y= "Eigenvalue")
    
}