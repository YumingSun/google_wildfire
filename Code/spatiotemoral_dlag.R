library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(gtrendsR)
scale2 <- function(x, na.rm = TRUE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
num_date_region = function(s_t_data){
  date_to_numeric_ea = setNames(seq(1:length(unique(s_t_data$week))),
                                unique(s_t_data$week))
  metro_to_numeric_ea = setNames(seq(1:length(unique(s_t_data$Metro))),
                                 unique(s_t_data$Metro))
  
  s_t_data$Metro_num = 1
  s_t_data$date_num_str = 1
  for (i in 1:dim(s_t_data)[1]) {
    s_t_data$Metro_num[i] = metro_to_numeric_ea[[s_t_data$Metro[i]]]
    s_t_data$date_num_str[i] = date_to_numeric_ea[[s_t_data$week[i]]]
  }
  s_t_data$date_num_unstr = s_t_data$date_num_str
  return(s_t_data)
}

get_sample = function(name,data){
  return(lapply(data, function(s) s$latent[paste0(name,":1"), 1]))
}

data_loc = "../Data/"
res_loc = "../Results/"
spatio_tempo_data = read_csv(paste0(data_loc,
                                    "master_data.csv"))
## Sensitivity
spatio_tempo_data = read_csv(paste0(data_loc,
                                    "master_data_sens.csv"))

spatio_tempo_data = spatio_tempo_data %>%
  filter(Metro != "Yuma, AZ-El Centro, CA") %>%
  filter(Metro != "Reno, NV") %>%
  filter(Metro != "Medford-Klamath Falls, OR")


spatio_tempo_data = spatio_tempo_data %>%
  mutate(
    across(`air pollution`:`air filter`, ~ pmax(.x, 1e-6))
  )

##sensitivity
spatio_tempo_data = spatio_tempo_data %>%
  mutate(
    across(`air filter`:`air purifier`, ~ pmax(.x, 1e-6))
  )

##INLA
spatio_tempo_data = spatio_tempo_data %>%
  mutate(Metro = factor(Metro)) %>%
  mutate_at(vars(starts_with("PM25")), scale2) %>%
  mutate_at(vars(starts_with("smokePM")), scale2) %>%
  mutate_at(
    c("temperature", "median_household_income", "population",
      "incident_counts", "MaxTemp", "MinTemp", "fireScale",
      "AvgTemp", "MaxHumidity", "MinHumidity", "AvgHumidity",
      "MaxWindSpeed", "MinWindSpeed", "AvgWindSpeed",
      "MaxPressure", "MinPressure", "AvgPressure",
      "MaxDewPoint", "MinDewPoint", "AvgDewPoint"
    ),
    scale2
  )%>%
  mutate_at(vars(contains("RPL")), scale2) %>%
  mutate_at(vars(contains("NOINT")), scale2) %>%
  select(-temperature, -MaxTemp, -MinTemp, -MaxHumidity,
         -MinHumidity, -incident_acres_burned,
         -MaxWindSpeed, -MinWindSpeed, -MaxPressure,
         -MinPressure, -MaxDewPoint, -MinDewPoint,
         -RPL_THEME1_median, -RPL_THEME2_median, -RPL_THEME3_median,
         -RPL_THEME4_median, -RPL_THEMES_median, -EP_NOINT_median, -E_NOINT_median,
         
  ) %>%
  rename(air_pollution = "air pollution",
         air_purifier = "air purifier",
         air_quality = "air quality",
         air_filter = "air filter")%>%
  # Sensitivity
  # rename(air_purifier = "air purifier",
  #        air_filter = "air filter")%>%
  group_by(Metro) %>%
  mutate(
    smokePM_lag1 = lag(smokePM, n = 1, order_by = week, default = 0),
    smokePM_lag2 = lag(smokePM, n = 2, order_by = week, default = 0),
    smokePM_lag3 = lag(smokePM, n = 3, order_by = week, default = 0),
    smokePM_lag4 = lag(smokePM, n = 4, order_by = week, default = 0),
    
    pm25_lag1 = lag(pm25, n = 1, order_by = week, default = 0),
    pm25_lag2 = lag(pm25, n = 2, order_by = week, default = 0),
    pm25_lag3 = lag(pm25, n = 3, order_by = week, default = 0),
    pm25_lag4 = lag(pm25, n = 4, order_by = week, default = 0),
    
    fireScale_lag1 = lag(fireScale, n = 1, order_by = week, default = 0),
    fireScale_lag2 = lag(fireScale, n = 2, order_by = week, default = 0),
    fireScale_lag3 = lag(fireScale, n = 3, order_by = week, default = 0),
    fireScale_lag4 = lag(fireScale, n = 4, order_by = week, default = 0),
    
    AvgTemp_lag1 = lag(AvgTemp, n = 1, order_by = week, default = 0),
    AvgTemp_lag2 = lag(AvgTemp, n = 2, order_by = week, default = 0),
    AvgTemp_lag3 = lag(AvgTemp, n = 3, order_by = week, default = 0),
    AvgTemp_lag4 = lag(AvgTemp, n = 4, order_by = week, default = 0),
    
    AvgHumidity_lag1 = lag(AvgHumidity, n = 1, order_by = week, default = 0),
    AvgHumidity_lag2 = lag(AvgHumidity, n = 2, order_by = week, default = 0),
    AvgHumidity_lag3 = lag(AvgHumidity, n = 3, order_by = week, default = 0),
    AvgHumidity_lag4 = lag(AvgHumidity, n = 4, order_by = week, default = 0),
    
    AvgWindSpeed_lag1 = lag(AvgWindSpeed, n = 1, order_by = week, default = 0),
    AvgWindSpeed_lag2 = lag(AvgWindSpeed, n = 2, order_by = week, default = 0),
    AvgWindSpeed_lag3 = lag(AvgWindSpeed, n = 3, order_by = week, default = 0),
    AvgWindSpeed_lag4 = lag(AvgWindSpeed, n = 4, order_by = week, default = 0),
    
    AvgPressure_lag1 = lag(AvgPressure, n = 1, order_by = week, default = 0),
    AvgPressure_lag2 = lag(AvgPressure, n = 2, order_by = week, default = 0),
    AvgPressure_lag3 = lag(AvgPressure, n = 3, order_by = week, default = 0),
    AvgPressure_lag4 = lag(AvgPressure, n = 4, order_by = week, default = 0),
    
    AvgDewPoint_lag1 = lag(AvgDewPoint, n = 1, order_by = week, default = 0),
    AvgDewPoint_lag2 = lag(AvgDewPoint, n = 2, order_by = week, default = 0),
    AvgDewPoint_lag3 = lag(AvgDewPoint, n = 3, order_by = week, default = 0),
    AvgDewPoint_lag4 = lag(AvgDewPoint, n = 4, order_by = week, default = 0),
  ) %>%
  ungroup()


dma_map = st_read(paste0(data_loc, 'NatDMA/', 'NatDMA.shp'))
cal_map = dma_map %>%
  mutate(NAME = gsub(" ([A-Z]{2})$", ", \\1", NAME)) %>%
  mutate(NAME = case_when(
    NAME == "Yuma AZ-El Centro, CA" ~ "Yuma, AZ-El Centro, CA",
    TRUE ~ NAME
  )) %>%
  filter(NAME %in% levels(spatio_tempo_data$Metro))

cal_map = st_transform(cal_map, crs = 3857)
cal_nb = poly2nb(cal_map)
adj_matrix = nb2mat(cal_nb, style = "B", zero.policy = TRUE)

rownames(adj_matrix) = cal_map$NAME 
colnames(adj_matrix) = cal_map$NAME


spatio_tempo_data = spatio_tempo_data %>%
  mutate(Metro = as.character(Metro),
         week = as.character(week))

spatio_tempo_data = num_date_region(spatio_tempo_data)
spatio_tempo_data$st_id = seq(1:2871)
#sensitivity
spatio_tempo_data$st_id = seq(1:2310)


nlags = 3
s_t_data = spatio_tempo_data
exposure = "smokePM"
outcome_vec = c("air_pollution","air_quality", "air_purifier", "air_filter")
# Sensitivity
outcome_vec = c("air_purifier", "air_filter")

time_dependent_interacts = c("AvgTemp", "AvgHumidity", "AvgWindSpeed",
                             "AvgPressure","fireScale")
time_independent_interacts = c("RPL_THEME1_mean", "RPL_THEME2_mean",
                               "RPL_THEME3_mean", "RPL_THEME4_mean")

lags = c("",paste0("_lag", seq(1,nlags)))
exposure_lags_vec = paste0(exposure,lags)
exposure_lags = paste(exposure_lags_vec, collapse = " + ")
time_dependent_interacts_lst = list()
time_dependent_posterior_name = list()
for (v in time_dependent_interacts) {
  v_interact_vec = paste0(exposure_lags_vec, ":", paste0(v,lags))
  v_interact = paste(v_interact_vec, collapse = " + ")
  time_dependent_posterior_name[[v]] = v_interact_vec
  time_dependent_interacts_lst[[v]] = v_interact
}
time_dependent_posterior_name[["smokePM2.5"]] = exposure_lags_vec
time_independent_interacts_lst = list()
time_independent_posterior_name = list()
for (v in time_independent_interacts){
  v_interact_vec = paste0(exposure_lags_vec, ":", v)
  v_interact = paste(v_interact_vec, collapse = " + ")
  time_independent_posterior_name[[v]] = v_interact_vec
  time_independent_interacts_lst[[v]] = v_interact
}

post_samp_lst = list()
fixed_effect_lst = list()
for (outcome in outcome_vec) {
  f = paste0(outcome, " ~ ", exposure_lags, " + ",
             paste0(time_dependent_interacts_lst, collapse = " + "), " + ",
             paste0(time_independent_interacts_lst, collapse = " + "),
             " + f(Metro_num, model = 'bym', graph = adj_matrix, constr = TRUE) + f(date_num_str, model='rw2', constr = TRUE) + f(st_id, model='iid', constr = TRUE)")
  model_formula = as.formula(f)
  
  mod = inla(
    formula = model_formula,
    family  = "gamma",                      
    data    = s_t_data,
    control.family = list(link = "log"),    
    control.predictor = list(compute = TRUE),
    control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.fixed     = list(
      mean = list(default = 0, intercept = 0),
      prec = list(default = 1e-4, intercept = 1e-4)
    ),
    quantiles = c(0.025, 0.5, 0.975),
    verbose=FALSE
  )
  fixed_effect = mod$summary.fixed
  fixed_effect = tibble::rownames_to_column(fixed_effect, "Coefficient") %>%
    filter(Coefficient != "(Intercept)") %>%
    dplyr::select(Coefficient, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`) %>%
    mutate(outcome = outcome,
           exposure = exposure)
  
  fixed_effect_lst[[paste0(outcome,'_',exposure)]] = fixed_effect
  
  samp = inla.posterior.sample(1000, mod)
  
  intercept_effect = lapply(samp, function(s) s$latent["(Intercept):1", 1])
  intercept_effect = do.call(cbind,intercept_effect)
  
  loc_random_effect = lapply(samp, function(s) s$latent[paste0("Metro_num:",seq(1:11)), 1])
  
  loc_random_effect = do.call(cbind,loc_random_effect)
  rownames(loc_random_effect) = unique(s_t_data$Metro)
  
  date_str_effect = lapply(samp, function(s) s$latent[paste0("date_num_str:",seq(1:261)), 1])
  # Sensitivity
  # date_str_effect = lapply(samp, function(s) s$latent[paste0("date_num_str:",seq(1:210)), 1])
  date_str_effect = do.call(cbind,date_str_effect)
  rownames(date_str_effect) = paste0("str:",unique(s_t_data$week))
  
  error_effect = lapply(samp, function(s) s$latent[paste0("st_id:",seq(1:2871)), 1])
  # Sensitivity
  # error_effect = lapply(samp, function(s) s$latent[paste0("st_id:",seq(1:2310)), 1])
  error_effect = do.call(cbind,error_effect)
  
  time_independent_samples = lapply(time_independent_posterior_name,
                                    get_sample,
                                    data = samp)
  time_independent_samples = lapply(time_independent_samples,
                                    function(s) do.call(cbind,s))
  
  time_independent_samples = do.call(rbind,time_independent_samples)
  
  time_dependent_samples = lapply(time_dependent_posterior_name,
                                  get_sample,
                                  data = samp)
  
  
  time_dependent_samples = lapply(time_dependent_samples,
                                  function(s) do.call(cbind,s))
  time_dependent_samples = do.call(rbind,time_dependent_samples)
  
  
  posterior_samp = rbind(time_dependent_samples,
                         time_independent_samples,
                         intercept_effect,
                         loc_random_effect,
                         date_str_effect,
                         error_effect)
  post_samp_lst[[paste0(outcome,'_',exposure)]] = posterior_samp
}

save(post_samp_lst, file = paste0(res_loc,'post_samp_lst_new_run.RData'))
save(fixed_effect_lst, file = paste0(res_loc,'fixed_effect_lst_new_run.RData'))


