#library(plyr) #for rbind.fill.matrix
library(tidyverse)

df<-read_table2(file="data/T11-7.txt", col_names = F)

df <- tibble(X1= c(6,5,8,4,7,3,1,2,2,5,3,2), X2= c(7,9,6,9,9,3,6,3,3,1,1,3),X3=c("A","A","A","A","A","B","B","B","C","C","C","C") )
df
# exercise 6.25
#Construct a one-way MANOVA of the crude-oil data listed in Table 11.7 on page 662.
#Construct 95% simultaneous confidence intervals to determine which mean components differ among the populations. 
#(You may want to consider transformations of the data to make them more closely conform to the usual MANOVA assumptions.)

# groups
group_column <- "X6"
group_column <- "X3"


groupnames<-df %>% select(group_column) %>% unique()
factors <- df %>% select(-group_column) %>% colnames()
n_factors <- length(factors)

X <- df %>% select(-group_column) %>% as.matrix()
X[,1]
groupnames[[1]][1]
g<-nrow(groupnames)

df %>% filter(group_column == groupnames[[1]][1])
df %>% select(group_column) %>% filter( group_column =="Wilhelm")

df 



#means per group
df %>% group_by(X3) %>% summarise_all(.funs = mean)

#overall means regardless of group
df %>% select(-group_column) %>% summarise_all(.funs = mean)


#########

#each col represents a group 
group_means <- df %>% group_by(X) %>% summarise_all(.funs = mean) %>% select(-X6) %>% as.matrix() %>% t()

#overall means regardless of group
overall_sample_mean <-df %>% select(-X6) %>% summarise_all(.funs = mean) %>% as.matrix() %>% t()

colSums(df$X1 %>% as.matrix())/56

n_j<-df %>% group_by(X6) %>% tally()
n_j$n
max(n_j$n) #max group number


get_n_from_g <- function(g){
    #groupnames[[1]][g]
    return( n_j %>% filter(X6 == groupnames[[1]][g] ) %>% select(n) %>% as.numeric() )
}


get_n_from_g(1)
get_n_from_g(2)
get_n_from_g(3)

 

r1<- df$X1[which(df$X6 %in% groupnames[[1]][1]) ] %>% as.matrix() %>% t()
r2<- df$X1[which(df$X6 %in% groupnames[[1]][2]) ] %>% as.matrix() %>% t()
r3<- df$X1[which(df$X6 %in% groupnames[[1]][3]) ] %>% as.matrix() %>% t()

rbind.fill.matrix(r1,r2)

#rows=groups, columns=observations
X_ij <- do.call("rbind.fill.matrix", list(r1,r2,r3) )
X_ij_template <-  replace(X_ij, X_ij != is.na(X_ij), 1)

overall_mean <- overall_sample_mean["X1",][[1]]
overall_mean <- replace(X_ij, X_ij != is.na(X_ij), overall_mean)
#estimated treatment effect: group mean - overall mean
group_mean <- group_means["X1",] %>% as.matrix() 

group_mean <- do.call("rbind", list(
replace(X_ij[1,] , X_ij[1,] != is.na(X_ij[1,]), group_mean[1,] ),
replace(X_ij[2,] , X_ij[2,] != is.na(X_ij[2,]), group_mean[2,] ),
replace(X_ij[3,] , X_ij[3,] != is.na(X_ij[3,]), group_mean[3,] )
))

#treatment effect:
est_treatment_effect <- group_mean-overall_mean

#residuals
residuals<-X_ij-group_mean


X_ij == overall_mean+est_treatment_effect+residuals

SS_obs<-sum(X_ij^2,na.rm = T)
SS_mean<-sum(overall_mean^2,na.rm = T)
SS_tr<-sum(est_treatment_effect^2,na.rm = T)
SS_res<-sum(residuals^2,na.rm = T)
SS_corrected <- SS_obs-SS_mean

sum(X_ij,na.rm = T)

#https://stats.idre.ucla.edu/spss/library/spss-library-my-sums-of-squares-dont-add-up/ 
SS_mean+SS_tr+SS_res
SS_obs
SS_corrected



X[,1] 
df$X1[which(df$X6 %in% groupnames[[1]][1]) ] %>% as.character() %>% paste(collapse = " ")


df[,"X6"]



get_indices_per_group <- function(df){
    indices_list <- list()
    for (group in 1:g){
        #convert to character since each list element can't have different lengths
        indices_list[group] <- which(df$X6 %in% groupnames[[1]][group]) %>% as.character() %>% paste(collapse = " ")
    }
    return(indices_list)
}

group_ind<-get_indices_per_group(df)

convert_string_to_numeric <- function(group_ind,g){
    test<-group_ind[[g]] %>% strsplit(split = " ") # group_ind[[g]] gets the indices in group g
    test<-lapply(test,as.numeric)[[1]] # [[1]] gets the vector only
    return(test)
}



test<-group_ind[[2]] %>% strsplit(split = " ") 
test<-lapply(test,as.numeric)[[1]]

as.numeric(test)
X[1:7]
convert_string_to_numeric(group_ind, g=3)



tmp_matrix <- matrix(NA, nrow=g, ncol=max(n_j$n) )
tmp_matrix

replace(tmp_matrix[1,], tmp_matrix[1,] != is.na(tmp_matrix[1,]) , convert_string_to_numeric(group_ind, g=1)  )


tmp_matrix[1, 1:get_n_from_g(1)] <- convert_string_to_numeric(group_ind, g=1)

X[convert_string_to_numeric(group_ind, g=1)  ,1]
X[1, convert_string_to_numeric(group_ind, g=1) ]


x_f <- X[,1]
x_f[convert_string_to_numeric(group_ind, g=3) ]

for (group in 1:g){
    #print(convert_string_to_numeric(group_ind, g=group))
    values_at_g <- X[convert_string_to_numeric(group_ind, g=group)  ,1] #change to X_f[convert_string_to_numeric(group_ind, g=group) ]
    tmp_matrix[group, 1:get_n_from_g(group) ] <- values_at_g
}
#str(tmp_matrix)



fill_observation_matrix <- function(tmp_matrix,f){
    
    for (group in 1:g){
        #print(convert_string_to_numeric(group_ind, g=group))
        values_at_g <- X[convert_string_to_numeric(group_ind, g=group), f] 
        tmp_matrix[group, 1:get_n_from_g(group) ] <- values_at_g
    }
    return(tmp_matrix)
}


obs_matrix<-fill_observation_matrix(tmp_matrix,f)
obs_matrix

#fill mean matrix
mean_value <-  overall_sample_mean[1,][[1]] #replace overall_sample_mean[1,] with overall_sample_mean[f,]

overall_sample_mean_matrix <- replace(obs_matrix, obs_matrix != is.na(obs_matrix), mean_value) # replace obs_matrix with its overall mean

#fill estimated treatment effect matrix
#group mean - overall mean
is.matrix(obs_matrix)
is.matrix(group_mean_matrix)
is.matrix(estimated_treatment_effect_matrix)
is.matrix(residual_matrix)
str(tmp_matrix)

fill_group_mean_matrix <- function(tmp_matrix,f){
    for (group in 1:g){
        tmp_matrix[group, 1:get_n_from_g(group) ] <- group_means[f,group ] #change to group_means[f,group]
    }
    return(tmp_matrix)
}
group_mean_matrix <- fill_group_mean_matrix(tmp_matrix)
estimated_treatment_effect_matrix <- group_mean_matrix - overall_sample_mean_matrix

residual_matrix <- obs_matrix-group_mean_matrix


obs_matrix == overall_sample_mean_matrix + estimated_treatment_effect_matrix + residual_matrix #TRUE

SS_obs<-sum(X_ij^2,na.rm = T)
SS_mean<-sum(overall_mean^2,na.rm = T)
SS_tr<-sum(est_treatment_effect^2,na.rm = T)
SS_res<-sum(residuals^2,na.rm = T)
SS_corrected <- SS_obs-SS_mean

str(obs_matrix)
str(overall_mean)

SS_obs<- get_sum_of_squares(obs_matrix)
SS_mean<- get_sum_of_squares(overall_sample_mean_matrix)
SS_tr<-get_sum_of_squares(estimated_treatment_effect_matrix)
SS_res<-get_sum_of_squares(residual_matrix)
SS_corrected <- SS_obs-SS_mean

get_sum_of_squares <- function(mat){
    return(sum( mat[!is.na(mat)]^2 ))
}










get_SS_decomposition_for_factor_f <- function(f){
    tmp_matrix <- matrix(NA, nrow=g, ncol=max(n_j$n))
    obs_matrix<-fill_observation_matrix(tmp_matrix,f)
    mean_value <-  overall_sample_mean[f,][[1]] 
    overall_sample_mean_matrix <- replace(obs_matrix, obs_matrix != is.na(obs_matrix), mean_value) # replaces obs_matrix with its overall mean
    group_mean_matrix <- fill_group_mean_matrix(tmp_matrix,f)
    estimated_treatment_effect_matrix <- group_mean_matrix - overall_sample_mean_matrix
    residual_matrix <- obs_matrix-group_mean_matrix
    
    SS_obs<- get_sum_of_squares(obs_matrix)
    SS_mean<- get_sum_of_squares(overall_sample_mean_matrix)
    SS_tr<-get_sum_of_squares(estimated_treatment_effect_matrix)
    SS_res<-get_sum_of_squares(residual_matrix)
    SS_corrected <- SS_obs-SS_mean
    

    return(
        list(
            SS_obs = SS_obs,
            SS_mean = SS_mean,
            SS_tr = SS_tr, 
            SS_res = SS_res,
            SS_corrected = SS_corrected,
            obs_matrix=obs_matrix,
            overall_sample_mean_matrix = overall_sample_mean_matrix,
            group_mean_matrix = group_mean_matrix,
            estimated_treatment_effect_matrix = estimated_treatment_effect_matrix,
            residual_matrix = residual_matrix
        )
    )
    
}


test<-get_SS_decomposition_for_factor_f(1)

test
B


B[f,f] <- SS_tr
W[f,f] <- SS_res 
T[f,f] <- SS_corrected







################### 



# TREATMENT MATRIX
B <- matrix(0,nrow = n_factors, ncol = n_factors)
# RESIDUAL MATRIX
W <- matrix(0,nrow = n_factors, ncol = n_factors)
# TOTAL MATRIX 
T <- matrix(0,nrow = n_factors, ncol = n_factors)


# treatment
#sum_l to g * n_l * (xbar_l - xbar) * (xbar_l - xbar)'

results_list <- list()
for (f in 1:n_factors){
    #print(f)
    #X_f <- X[,f] 
    
    
    # step 1: calculate SS for each variable
    SS <- get_SS_decomposition_for_factor_f(f)
    results_list[[f]] <- SS
    
    B[f,f] <- SS$SS_tr
    W[f,f] <- SS$SS_res 
    T[f,f] <- SS$SS_corrected
    
}


B
W
T

results_list[[1]]$obs_matrix%*%results_list[[2]]$obs_matrix
convert_NA_to_zero <- function(mat){
    mat[is.na(mat)] <- 0
    return(mat)
}
convert_NA_to_zero(results_list[[1]]$obs_matrix)%*%convert_NA_to_zero(results_list[[2]]$obs_matrix)

multiply_and_sum <- function(i,j, mat){
    
    S<-sum(
        crossprod(
            convert_NA_to_zero(results_list[[i]][mat]),
            convert_NA_to_zero(results_list[[j]][mat])
        )
    )
    return(S)
}

crossprod(results_list[[4]]$estimated_treatment_effect_matrix[,1] , results_list[[5]]$estimated_treatment_effect_matrix[,1])

unique_tr <- results_list[[4]]$estimated_treatment_effect_matrix[,1]

n_j

get_treatment_crossproduct <- function(i,j){
    tot <- 0
    for (group in 1:g){
        group_mean_one <- results_list[[i]]$estimated_treatment_effect_matrix[g,1] 
        group_mean_two <- results_list[[j]]$estimated_treatment_effect_matrix[g,1] 
        n<-get_n_from_g(group)
        tot = tot + n*group_mean_one*group_mean_two
    }
    return(tot)
}

get_treatment_crossproduct(1,4)


sum(
    crossprod(
        convert_NA_to_zero(results_list[[2]]$obs_matrix),
        convert_NA_to_zero(results_list[[3]]$obs_matrix)
    )
)

sum()

sum(crossprod(convert_NA_to_zero(results_list[[f]]$obs_matrix),convert_NA_to_zero(results_list[[f_col]]$obs_matrix)))

multiply_and_sum(2,3, "obs_matrix")

# CROSS PRODUCTS
for (f in 1:n_factors){
    SS_row <- results_list[[f]]
    for(f_col in 1:n_factors){
        if (f!=f_col){
            SS_col <- results_list[[f_col]]
            total<- sum(crossprod(convert_NA_to_zero(results_list[[f]]$obs_matrix),convert_NA_to_zero(results_list[[f_col]]$obs_matrix)))
            
            mean<-sum(crossprod(convert_NA_to_zero(results_list[[f]]$overall_sample_mean_matrix),convert_NA_to_zero(results_list[[f_col]]$overall_sample_mean_matrix)))
            treatment <- get_treatment_crossproduct(f,f_col)
            residual<-sum(crossprod(convert_NA_to_zero(results_list[[f]]$residual_matrix),convert_NA_to_zero(results_list[[f_col]]$residual_matrix)))
            obs_corrected <- total-mean
            B[f,f_col] <- treatment
            W[f,f_col] <- residual
            T[f,f_col] <- obs_corrected
        }
    }
}
W
T
B

det(W)/det(B+W)









