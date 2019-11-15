library(tidyverse)
df <- tibble(X1= c(6,5,8,4,7,3,1,2,2,5,3,2), X2= c(7,9,6,9,9,3,6,3,3,1,1,3),X3=c("A","A","A","A","A","B","B","B","C","C","C","C") )
df<-read_table2(file="data/T11-7.txt", col_names = F)
df
group_column <- "X3"
group_column <- "X6"
groupnames<-df %>% select(group_column) %>% unique()
factors <- df %>% select(-group_column) %>% colnames()
n_factors <- length(factors)

X <- df %>% select(-group_column) %>% as.matrix()
X[,1]
groupnames[[1]][1]
g<-nrow(groupnames)

df %>% filter(group_column == groupnames[[1]][1])
df %>% select(group_column) %>% filter( group_column =="Wilhelm")



df<-as.data.frame(df)
#means per group

df[[group_column]]

#overall means regardless of group
df %>% select(-group_column) %>% summarise_all(.funs = mean)

#each col represents a group 
group_means <- df %>% group_by_at(group_column) %>% summarise_all(.funs = mean) %>% select(-group_column) %>% as.matrix() %>% t()

#overall means regardless of group
overall_sample_mean <-df %>% select(-group_column) %>% summarise_all(.funs = mean) %>% as.matrix() %>% t()


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
obs


#fill mean matrix
mean_value <-  overall_sample_mean[1,][[1]] #replace overall_sample_mean[1,] with overall_sample_mean[f,]
overall_sample_mean_matrix <- replace(obs_matrix, obs_matrix != is.na(obs_matrix), mean_value) # replace obs_matrix with its overall mean



fill_group_mean_matrix <- function(tmp_matrix,f){
    for (group in 1:g){
        tmp_matrix[group, 1:get_n_from_g(group) ] <- group_means[f,group ] 
    }
    return(tmp_matrix)
}

group_mean_matrix <- fill_group_mean_matrix(tmp_matrix,1)
estimated_treatment_effect_matrix <- group_mean_matrix - overall_sample_mean_matrix

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

test<-get_SS_decomposition_for_factor_f(2)
test


################### 



# TREATMENT MATRIX
B <- matrix(0,nrow = n_factors, ncol = n_factors)
# RESIDUAL MATRIX
W <- matrix(0,nrow = n_factors, ncol = n_factors)
# TOTAL MATRIX 
T <- matrix(0,nrow = n_factors, ncol = n_factors)


results_list <- list()
for (f in 1:n_factors){
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



p<-as.vector(results_list[[1]]$obs_matrix)*as.vector(results_list[[2]]$obs_matrix) 
sum(as.vector(results_list[[1]]$obs_matrix)*as.vector(results_list[[2]]$obs_matrix) ,na.rm = T)
sum(p,na.rm = T)


results_list[[1]]$estimated_treatment_effect_matrix[,1]




get_treatment_crossproduct <- function(i,j){
    tot <- 0
    for (group in 1:g){
        group_mean_one <- results_list[[i]]$estimated_treatment_effect_matrix[group,1] 
        group_mean_two <- results_list[[j]]$estimated_treatment_effect_matrix[group,1] 
        n<-get_n_from_g(group)
        tot = tot + n*group_mean_one*group_mean_two
        cat(group ,"\n")
    }
    return(tot)
}

get_n_from_g(1)*results_list[[1]]$estimated_treatment_effect_matrix[1,1]*results_list[[2]]$estimated_treatment_effect_matrix[1,1] +
    get_n_from_g(2)*results_list[[1]]$estimated_treatment_effect_matrix[2,1]*results_list[[2]]$estimated_treatment_effect_matrix[2,1] +
    get_n_from_g(3)*results_list[[1]]$estimated_treatment_effect_matrix[3,1]*results_list[[2]]$estimated_treatment_effect_matrix[3,1]


get_treatment_crossproduct(1,2)





# CROSS PRODUCTS
for (f in 1:n_factors){
    for(f_col in 1:n_factors){
        #if (f!=f_col){
            total <- sum(as.vector(results_list[[f]]$obs_matrix)*as.vector(results_list[[f_col]]$obs_matrix) ,na.rm = T)
            
            mean <- sum(as.vector(results_list[[f]]$overall_sample_mean_matrix)*as.vector(results_list[[f_col]]$overall_sample_mean_matrix) ,na.rm = T)
            treatment <- get_treatment_crossproduct(f,f_col)
            residual<- sum(as.vector(results_list[[f]]$residual_matrix)*as.vector(results_list[[f_col]]$residual_matrix) ,na.rm = T)
            obs_corrected <- total-mean
            B[f,f_col] <- treatment
            W[f,f_col] <- residual
            T[f,f_col] <- obs_corrected
        #}
    }
}

W
T
B

det(W)/det(B+W)














nl <- c(271,138,107)
n<- sum(nl)

x1 <- matrix(c(2.066,0.480,0.082,0.360), nrow=1, ncol=4)
x2 <- matrix(c(2.167,0.596,0.124,0.418), nrow=1, ncol=4)
x3 <- matrix(c(2.273,0.521,0.125,0.383), nrow=1, ncol=4)
xbar <- round((nl[1]*x1 + nl[2]*x2 + nl[3]*x3 )/ n ,digits = 3)
xbar
x1
x2
x3


?manova
manova()
#X3 = group
m1<-manova(cbind(X1, X2) ~ X3, data = df)
summary(m1,test = "Wilks")



aov(data = df, formula = )
