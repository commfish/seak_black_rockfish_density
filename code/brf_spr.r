# black rockfish selectivity at age
# ben.williams@alaska.gov

# load ----
source("code/helper.r")

setwd(here::here("code"))

# data ----

# note if the lengths are not rounded there are occasionaly more precise estimates
# these will not allow the model to converge
# combine sport and comm age data

# Female only data for sport and commercial fisheries
read_csv(here::here("data/br_bio.csv"), guess = 50000) %>% 
  rename_all(tolower) %>% 
  dplyr::select(year, Area = g_management_area_code, 
                length = length_millimeters, age, sex = sex_code,
                weight = weight_kilograms) %>% 
  filter(sex == 2) %>%
  dplyr::select(-sex) %>%
  mutate(length = round(length / 10),
         fish = "comm") %>% 
  bind_rows(read_csv(here::here("data/sport_brf_bio_se.csv"), guess = 50000) %>% 
              rename_all(tolower) %>% 
              filter(sex=="F") %>%
              dplyr::select(year, Area = area, 
                            length, age, sex , weight = wt_kg) %>% 
              dplyr::select(-sex) %>%
              mutate(length = round(length / 10),
                     fish = "sport")) %>% 
  filter(Area %in% c("CSEO"), !is.na(age), !is.na(length), 
         year!=2001) -> brf

# read_csv(here::here("data/kodiak_brf.csv"), guess = 50000) %>%
#   filter(age>0 & age<26, sex==2, !is.na(age), !is.na(length),
#          year>2001) -> brf

# estimate length/weight relationship

lw <- unname(lm(log(weight) ~ log(length), data = brf)$coef)
lw_int = (lw[1])
lw_slope = (lw[2])

data.frame(length = 1:70) %>% 
  mutate(weight = exp(lw_int + lw_slope * log(length)))  %>% 
  ggplot(aes(length, weight)) + 
  geom_line() + 
  geom_point(data = brf, alpha = 0.2)


# maturity at age from Kodiak study

mat_int = -7.521637
mat_slope = 0.717806

data.frame(age = 1:25) %>% 
  mutate(maturity = exp(mat_int + mat_slope * age) / (1 + exp(mat_int + mat_slope * age)))  %>% 
  ggplot(aes(age, maturity)) + 
  geom_line()

# choose data ----
# length based
# 
# dat <- clean_up(brf, length) %>%
#   left_join(data.frame(X = 20:65), .) %>%
#   mutate(prop = replace_na(prop, 0)) -> dat

# age based 

clean_up(brf, age) %>% 
  # group_by(year) %>% 
  left_join(expand.grid(X = 1:30), .) %>%
  mutate(prop = replace_na(prop, 0)) -> dat

# pre-flight check that the data are in a good form and starting values are reasonable 
# von b
Linf <- 54
kappa <- 0.15
t0 <- -0.5

brf %>% 
  mutate(fit = Linf * (1 - exp(-kappa * (age - t0)))) %>% 
  ggplot(aes(age, length)) + 
  geom_point() + 
  geom_line(aes(y = fit)) +
  expand_limits(y = 0)

# selectivity

alpha = 10
beta = -1.0
M = 0.123
F = 0.07

dat %>% 
  mutate(sel =  1 / (1 + exp(alpha + beta * X ))) -> dd

dd$N <- 1

for(i in 2:nrow(dd)){
  dd$N[i] = dd$N[i-1] * exp(-M -F * dd$sel[i-1])
}

dd %>%
  # group_by(year) %>% 
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
            prop = dat$prop,
            lw_int = lw_int,
            lw_slope = lw_slope,
            mat_b0 = mat_int,
            mat_b1 = mat_slope,
            mu_M = 0.13, # based on barefoot ecologist tool
            sigma_M = 0.01)

# starting parameters
params <- list(Linf = 55, kappa = 0.15, t0 = -0.5, log_von_sigma = 0.001,
               Fcur = 0.07, M = 0.13,
               alpha = 10, beta = -1.0, log_sel_sigma = 0.0001)

# add parameter bounds 
# map = list(M = factor(NA))
           # alpha = factor(NA), beta = factor(NA)) #M = factor(NA)
map = list()

L = list(Linf = 20, kappa = 0.05, t0 = -10.0, log_von_sigma = 0.0001,
         Fcur = 0.03, 
         M = 0.05,
         alpha = 5, 
         beta = -2.5, 
         log_sel_sigma = 0.0001)

U = list(Linf = 70, kappa = 0.35, t0 = 5.0, log_von_sigma = 10,
         Fcur = 0.2, 
         M = 0.3,
         alpha = 17, 
         beta = -0.9, 
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
rep

best
summary(rep, select = "all", p.value = T)


yfit <- model$report()$yfit
fit_prop <- model$report()$fit_prop
a <- model$report()$alpha
b <- model$report()$beta
N <- model$report()$Na
unN <- model$report()$UnNa
S <- model$report()$Sa
C <- model$report()$Ca
M <- model$report()$M
Fcur <- model$report()$Fcur
fished <- model$report()$fished 
unfished <- model$report()$unfished 
model$report()$exploitB

ifelse(is.nan(fished), 0, fished) -> fished
ifelse(is.nan(unfished), 0, unfished) -> unfished

data.frame(age = dat$X,
           fished = fished, 
           unfished = unfished) %>% 
  ggplot(aes(age, unfished)) + 
  geom_line() + 
  geom_line(aes(y = fished), color = 4) + 
  expand_limits(y = c(0, 0.03))
sum(fished) / sum(unfished) # SPR

# von b

brf %>% 
  mutate(fit = yfit) %>% 
  ggplot(aes(age, length)) + 
  geom_point(alpha = 0.1) + 
  geom_line(aes(y = yfit)) + 
  expand_limits(x = 0, y = 0)

# age comp 
dat %>% 
  mutate(sel = S,
         catch = C,
         numbers = N,
         fit = fit_prop) %>% 
  ggplot(aes(X, prop)) + 
  geom_bar(stat = "identity", alpha = 0.3) +
  geom_line(aes(y = fit), color = 1)

fem_A50 <- -a / b # 50% selected
# female <- rep
# fem_M <- M
# male_A50 <- -a / b
# male <- rep
# male_M <- M
