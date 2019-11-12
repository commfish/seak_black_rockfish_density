# black rockfish selectivity at age
# ben.williams@alaska.gov

# load ----
source(here::here("code/helper.r"))

setwd(here::here("code"))

# data ----

# note if the lengths are not rounded there are occasionaly more precise estimates
# these will not allow the model to converge
# combine sport and comm age data

read_csv(here::here("data/br_bio.csv"), guess = 50000) %>% 
  rename_all(tolower) %>% 
  dplyr::select(year, Area = g_management_area_code, 
                length = length_millimeters, age, sex = sex_code) %>% 
  filter(sex == 2) %>%
  dplyr::select(-sex) %>%
  mutate(length = round(length / 10)) %>% 
  bind_rows(read_csv(here::here("data/sport_brf_bio_se.csv"), guess = 50000) %>% 
              rename_all(tolower) %>% 
              filter(sex=="F") %>%
              dplyr::select(year, Area = area, 
                            length, age, sex ) %>% 
              dplyr::select(-sex) %>%
              mutate(length = round(length / 10))) %>% 
  drop_na() %>% 
  filter(Area=="CSEO", !is.na(age), !is.na(length)) -> brf


# choose data ----
# length based
# 
# dat <- clean_up(brf, length) %>%
#   left_join(data.frame(X = 20:65), .) %>%
#   mutate(prop = replace_na(prop, 0)) -> dat

# age based 

clean_up(brf, age) %>% 
  left_join(data.frame(X = 1:30), .) %>%
  mutate(prop = replace_na(prop, 0)) -> dat

# pre-flight check that the data are in a good form and starting values are reasonable 
# von b
Linf <- 55
kappa <- 0.15
t0 <- -0.5

brf %>% 
  mutate(fit = Linf * (1 - exp(-kappa * (age - t0)))) %>% 
  ggplot(aes(age, length)) + 
  geom_point() + 
  geom_line(aes(y = fit)) +
  expand_limits(y = 0)

# selectivity

alpha = 9
beta = -1.0
M = 0.123
F = 0.1

dat %>% 
  mutate(sel =  1 / (1 + exp(alpha + beta * X ))) -> dd

dd$N <- 1

for(i in 2:nrow(dd)){
  dd$N[i] = dd$N[i-1] * exp(-M -F * dd$sel[i-1])
}

dd %>% 
  mutate(C = N * (1 - exp(-M -F * sel)) * F * sel/(M + F * sel),
         fit_prop = C / max(C)) %>% 
  ggplot(aes(X, prop)) + 
  geom_bar(stat = "identity", alpha = 0.3) +
  geom_line(aes(y = fit_prop))

# TMB model ----
compile("brf_spr.cpp")
dyn.load(dynlib("brf_spr"))

# data 
data = list(age = brf$age,
            length = brf$length,
            X = dat$X, 
            prop = dat$prop)

# starting parameters
params <- list(Linf = 50, kappa = 0.15, t0 = 1.0, log_von_sigma = 0.001,
               Fcur = 0.10, M = 0.123,
               alpha = 9, beta = -1.0, log_sel_sigma = 0.0001)

# add parameter bounds 
map = list(M = factor(NA), alpha = factor(NA)) #M = factor(NA)
L = list(Linf = 20, kappa = 0.05, t0 = -10.0, log_von_sigma = 0.0001,
         Fcur = 0.01, 
         # M = 0.123,
         #alpha = 5, 
         beta = -3, 
         log_sel_sigma = 0.0001)
U = list(Linf = 70, kappa = 0.35, t0 = 5.0, log_von_sigma = 10,
         Fcur = 0.15, 
         # M = 0.123,
         #alpha = 15, 
         beta = -0.05, 
         log_sel_sigma = 10)

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL="brf_spr",
                   map = map)

# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr,
              upper = U,
              lower = L)

best <- model$env$last.par.best
rep <- sdreport(model)

best
rep


yfit <- model$report()$yfit
fit_prop <- model$report()$fit_prop
a <- model$report()$alpha
b <- model$report()$beta
N <- model$report()$Na
S <- model$report()$Sa
C <- model$report()$Ca

# von b

brf %>% 
  mutate(fit = yfit) %>% 
  ggplot(aes(age, length)) + 
  geom_point(alpha = 0.1) + 
  geom_line(aes(y = yfit)) + 
  expand_limits(y = 0)

# selectivity 
dat %>% 
  mutate(sel = S,
         catch = C,
         numbers = N,
         fit = fit_prop) %>% 
  ggplot(aes(X, prop)) + 
  geom_bar(stat = "identity", alpha = 0.3) +
  geom_line(aes(y = fit), color = 1)

fem_A50 <- -a / b
female <- rep
male_A50 <- -a / b
male <- rep
