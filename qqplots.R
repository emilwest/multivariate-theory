x1<-c(62.97,65.45,92.01,95.04,108.28,152.36,165.68,263.99,265.19,285.06)
x1

x2<-c(17.05,16.59,10.91,14.14,9.52,25.33,18.54,15.73,8.10,11.13)
x2 <- sort(x2)
x2
n<-10
mean(x1)
mean(x2)

i<-1:10

p<-(i-0.5)/n
p #pvalues
q<-qnorm(p) #critvalues
plot(q,x1) 
plot(q,x2)


cor(q,x1)
cor(q,x2)






########################
# Ex 4.37
library(tidyverse)
df<-read.table(file = "T1-9.txt", sep = "")
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
         marathonms = sort(42195/(marathon*60)),
         i = 1:nrow(df),
         p = (i-0.5)/nrow(df),
         q = qnorm(p),
         x = ( mat[i,] - t(x_bar)[1,] ) %>% as.matrix(),
         d = t(x)%*%solve(cov(mat))%*%x
         
         ) %>%
  select(country, hundrams, twohndrams, fourhundams, eightms, femtonms, threems, marathonms,
         i, p, q, d)

df_meters_per_second

## QQ plots
plot(df_meters_per_second$q, df_meters_per_second$hundrams)
plot(df_meters_per_second$q, df_meters_per_second$twohndrams)
plot(df_meters_per_second$q, df_meters_per_second$fourhundams)
plot(df_meters_per_second$q, df_meters_per_second$eightms)
plot(df_meters_per_second$q, df_meters_per_second$femtonms)
plot(df_meters_per_second$q, df_meters_per_second$threems)
plot(df_meters_per_second$q, df_meters_per_second$marathonms)


### TESTING MARGINAL NORMALITY WITH Q-Q PLOTS
library(reshape2)

melt(df_meters_per_second %>% select(-country,-p,-i) , id.vars = "q" ) %>%
  ggplot(aes(x=q,y=value)) +
  geom_point( mapping = aes(color=variable)) +
  labs(title = "Q-Q plots of meters per second at different running distances",
       x= "Q=standardized N(0,1) values",
       y= "Meters per second")
  
#marathon -> not a straight line -> probably not normal


# alternativt:
melt(df_meters_per_second %>% select(-country,-p,-i) , id.vars = "q" ) %>%
  ggplot(aes(x=q,y=value)) +
  geom_point(mapping = aes(color=variable)) +
  geom_smooth(method = "lm",aes(q,value), se=F) + #GEOM_SOOOTH FUNKAR PÅ FACET 
  facet_wrap(~variable, scales="free_y") + #free_y makes it easier to see plots
  labs(title = "Q-Q plots of meters per second at different running distances",
       x= "Q=standardized N(0,1) values",
       y= "Meters per second")




# CORRELATIONS
cor(df_meters_per_second %>%
      select(hundrams:marathonms) , df_meters_per_second$q)


                                                        
### TESTING JOINT NORMALITY WITH CHI SQUARE PLOTS

x_bar <- df_meters_per_second %>%
  select(hundrams, twohndrams, fourhundams, eightms, femtonms, threems, marathonms) %>%
  summarise_all(mean) %>% as.matrix() %>% t()

mat<-df_meters_per_second %>%
  select(hundrams, twohndrams, fourhundams, eightms, femtonms, threems, marathonms) %>% 
  as.matrix() 

x_bar

S<-colSums( (mat-x_bar[,1])^2 )/nrow(mat)

d_j<-t(x)%*%solve(cov(mat))%*%x
d_j

x<-(mat[44,]-x_bar) %>% as.matrix()
x
t(x)
i<-1:nrow(df)
x = ( mat[i,] - t(x_bar)[1,] ) %>% as.matrix()
d = t(x)%*%solve(cov(mat))%*%x
