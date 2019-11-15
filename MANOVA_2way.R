library(tidyverse)

df <- tibble(X1= c(6,4,8,2,3,-3,4,-4,-3,-4,3,-4), X2 = c(8,6,12,6,8,2,3,3,2,-5,-3,-6), l=c(1,1,1,1,2,2,2,2,3,3,3,3) ,k=rep(1:4,3) )
df

rowmean <- df %>% select(-k) %>% group_by(l) %>% summarise_all(.funs = "mean")

colmean <- df %>% select(-l) %>% group_by(k) %>% summarise_all(.funs = "mean")

overallmean <- df %>% select(-l,-k) %>% summarise_all(.funs = "mean")
groupnames<-df %>% select(group_column) %>% unique()


group_means <- rowmean %>% select(X1,X2) %>% as.matrix() %>% t() #Xbar_l
row_means <- colmean %>% select(X1,X2) %>% as.matrix() %>% t() #Xbar_k
overall_sample_mean <-  overallmean %>% select(X1,X2) %>% as.matrix() %>% t() #Xbar

g<-nrow(groupnames)
b<-ncol(row_means)
X <- df %>% select(-group_column,-k) %>% as.matrix()
X
group_column<-"l"

n_j<-df %>% group_by_at(group_column) %>% tally()
n_j
colnames(n_j) <- c("group","n")
n_j$n
max(n_j$n) #max group number

get_n_from_g <- function(g){
    #groupnames[[1]][g]
    return( n_j %>% filter(group == groupnames[[1]][g] ) %>% select(n) %>% as.numeric() )
}

get_n_from_g(1)

get_indices_per_group <- function(df){
    indices_list <- list()
    for (group in 1:g){
        #convert to character since each list element can't have different lengths
        indices_list[group] <- which(df[group_column][[1]] %in% groupnames[[1]][group]) %>% as.character() %>% paste(collapse = " ")
    }
    return(indices_list)
}

group_ind<-get_indices_per_group(df)

convert_string_to_numeric <- function(group_ind,g){
    test<-group_ind[[g]] %>% strsplit(split = " ") # group_ind[[g]] gets the indices in group g
    test<-lapply(test,as.numeric)[[1]] # [[1]] gets the vector only
    return(test)
}


tmp_matrix <- matrix(NA, nrow=g, ncol=max(n_j$n) )
fill_observation_matrix <- function(tmp_matrix,f){
    
    for (group in 1:g){
        #print(convert_string_to_numeric(group_ind, g=group))
        values_at_g <- X[convert_string_to_numeric(group_ind, g=group), f] 
        tmp_matrix[group, 1:get_n_from_g(group) ] <- values_at_g
    }
    return(tmp_matrix)
}
obs_matrix<-fill_observation_matrix(tmp_matrix,1)
obs_matrix


#fill mean matrix
mean_value <-  overall_sample_mean[1,][[1]] #replace overall_sample_mean[1,] with overall_sample_mean[f,]
overall_sample_mean_matrix <- replace(obs_matrix, obs_matrix != is.na(obs_matrix), mean_value) # replace obs_matrix with its overall mean

row_means[1,] 

fill_group_mean_matrix <- function(tmp_matrix,f){
    for (group in 1:g){
        tmp_matrix[group, 1:get_n_from_g(group) ] <- group_means[f,group ] 
    }
    return(tmp_matrix)
}

group_mean_matrix <- fill_group_mean_matrix(tmp_matrix,1)
estimated_treatment_effect_matrix <- group_mean_matrix - overall_sample_mean_matrix


fill_row_mean_matrix <- function(tmp_matrix,f){
    for (row in 1:b){
        tmp_matrix[1:3, row ] <- row_means[f,row] 
    }
    return(tmp_matrix)
}


row_mean_matrix <- fill_row_mean_matrix(tmp_matrix,1)
estimated_treatment_effect_factor2_matrix <- row_mean_matrix - overall_sample_mean_matrix



residual_matrix <- obs_matrix-group_mean_matrix-row_mean_matrix+overall_sample_mean_matrix







residual_matrix <- obs_matrix-group_mean_matrix

obs_matrix == overall_sample_mean_matrix + estimated_treatment_effect_matrix + residual_matrix #TRUE


get_sum_of_squares <- function(mat){
    return(sum( mat[!is.na(mat)]^2 ))
}

SS_obs<- get_sum_of_squares(obs_matrix)
SS_mean<- get_sum_of_squares(overall_sample_mean_matrix)
SS_tr<-get_sum_of_squares(estimated_treatment_effect_matrix)
SS_res<-get_sum_of_squares(residual_matrix)
SS_corrected <- SS_obs-SS_mean






get_SS_decomposition_for_factor_f <- function(f){
    tmp_matrix <- matrix(NA, nrow=g, ncol=max(n_j$n))
    obs_matrix<-fill_observation_matrix(tmp_matrix,f)
    mean_value <-  overall_sample_mean[f,][[1]] 
    overall_sample_mean_matrix <- replace(obs_matrix, obs_matrix != is.na(obs_matrix), mean_value) # replaces obs_matrix with its overall mean
    group_mean_matrix <- fill_group_mean_matrix(tmp_matrix,f)
    estimated_treatment_effect_matrix <- group_mean_matrix - overall_sample_mean_matrix
    row_mean_matrix <- fill_row_mean_matrix(tmp_matrix,f)
    estimated_treatment_effect_factor2_matrix <- row_mean_matrix - overall_sample_mean_matrix
    
    
    
    residual_matrix <- obs_matrix-group_mean_matrix-row_mean_matrix+overall_sample_mean_matrix
    
    #residual_matrix <- obs_matrix-group_mean_matrix
    
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

test<-get_SS_decomposition_for_factor_f(2)
test


################### 



# TREATMENT MATRIX
B <- matrix(0,nrow = g, ncol = b)
# RESIDUAL MATRIX
W <- matrix(0,nrow = g, ncol = b)
# TOTAL MATRIX 
T <- matrix(0,nrow = g, ncol = b)


results_list <- list()
for (f in 1:n_factors){
    # step 1: calculate SS for each variable
    SS <- get_SS_decomposition_for_factor_f(f)
    results_list[[f]] <- SS
    
    B[f,f] <- SS$SS_tr
    W[f,f] <- SS$SS_res 
    T[f,f] <- SS$SS_corrected
}
results_list




