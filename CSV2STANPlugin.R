# replace "~/Desktop/Shared Code for Soon" by "folder-location"

# load libraries 
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
n_chains <- 4


# load data  
input <- function(inputfile) {
dd <<- read_csv(file=str_c(inputfile) ) %>% 
                select(Antibiotic,Setting,
		       Rcgc_AMR_ALL,
		       Rcgc_AMR_DEF,
		       Rcgc_AMR_DEF_zero_replaced,
		       Rtax_e4,
		       Rtax_enterobacteriaceae,
		       Rtax_enterobacterales,
		       Observed_Resistant_Infections,
		       Total_Observed_Infections_With_AST_Data)
}
run <- function() {}
output <- function(outputfile) {
#  split column Antibiotic
dd1 <- dd %>% separate(Antibiotic, into=c("Antibiotic","count"),sep = " " ) %>% select(-count)
# rename cols
names(dd1)[3:8] <- c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv")

n_atbs <- dd1 %>% group_by(Antibiotic) %>% n_groups() # 16
saveRDS(n_atbs, paste(outputfile, "n_atbs.rds", sep="/"))
#  extract infection data
df_transf2 <- dd1 %>% select(Setting,Antibiotic,Total_Observed_Infections_With_AST_Data) %>% 
                spread(key=Antibiotic ,value=Total_Observed_Infections_With_AST_Data)
df_transf3 <- dd1 %>% select(Setting,Antibiotic,Observed_Resistant_Infections) %>% 
                spread(key=Antibiotic ,value=Observed_Resistant_Infections)
#
setting_v <- df_transf2 %>% .$Setting
count_s <- df_transf2  %>% select(-Setting) %>% as.matrix() 
count_sr <- df_transf3  %>% select(-Setting) %>% as.matrix()
saveRDS(count_s, paste(outputfile, "count_s.rds", sep="/"))
saveRDS(count_sr, paste(outputfile, "count_sr.rds", sep="/"))
saveRDS(setting_v, paste(outputfile, "setting_v.rds", sep="/"))
data_exists <- replace(count_sr, !is.na(count_sr), 1 ) # make helper that says where data is
data_exists <- replace(data_exists , is.na(data_exists), 0 ) # and where it isn't
count_s_hna <- replace(count_s,is.na(count_s),0) # hide NAs with 0
count_sr_hna <- replace(count_sr,is.na(count_sr),0) # hide NAs with 0

# taxonomic predictors (standadise)
df_transf <- dd1 %>% select(Setting,Antibiotic,b4div) %>% 
                spread(key=Antibiotic ,value=b4div)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_1a <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,entdiv) %>% 
                spread(key=Antibiotic ,value=entdiv)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_1b <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,entbacdiv) %>% 
                spread(key=Antibiotic ,value=entbacdiv)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_1c <- (z - mean(z)) / sd(z) # standardise

# CARD predictors (standadise)
df_transf <- dd1 %>% select(Setting,Antibiotic,ar12) %>% 
                spread(key=Antibiotic ,value=ar12)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_2a <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,ar2) %>% 
                spread(key=Antibiotic ,value=ar2)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_2b <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,ar1o2) %>% 
                spread(key=Antibiotic ,value=ar1o2)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_2c <- (z - mean(z)) / sd(z) # standardise

# impute counts where counts was 0 based on imputing from two other settings
count_s_for_pred <- count_s
saveRDS(count_s_for_pred, paste(outputfile, "count_s_for_pred.rds", sep="/"))
for (i in 1:ncol(count_s)) {
                icol <-  count_s[,i]
                if (sum(is.na(icol)) > 0) {
                                icol_nonna <- icol[!is.na(icol)] 
                                count_s_for_pred[which(is.na(icol)),i] <- round(mean(icol_nonna))
                }
}
}

