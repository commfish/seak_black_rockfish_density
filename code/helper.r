# helper file
# libraries and functions

library(TMB)
library(tidyverse)
library(fngr)
theme_set(theme_report())

# scale data to 1
range01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}


clean_up_II <- function(data, variable){
  
  brf %>% 
    group_by(year) %>%
    tally %>% 
    ungroup() %>% 
    mutate(tot = sum(n)) -> smp_size
  
  brf %>% 
    dplyr::select(X = {{variable}}, year) %>% 
    drop_na(X) %>% 
    group_by(X, year) %>% 
    summarise(nn = n()) %>% 
    ungroup() %>% 
    mutate(prop = nn / sum(nn)) %>% 
    left_join(smp_size) %>% 
    mutate(prop = prop * nn) %>% 
    group_by(X) %>% 
    summarise(prop = sum(prop) / mean(tot)) %>% 
    ungroup %>% 
    mutate(prop = range01(prop / sum(prop))) 
  
}

clean_up <- function(data, variable){
  
  data %>% 
    dplyr::select(X = {{variable}}, year) %>% 
    drop_na(X) %>% 
    group_by(X) %>% 
    summarise(n = n()) %>% 
    # group_by(year) %>% 
    mutate(prop = range01(n / sum(n))) %>% 
    dplyr::select(-n)
  
}
