library(tidyverse)
library(mapdata)
library(maps)
library(viridis)
library(gtools)
library(igraph)
library(ggpubr)
library(viridis)
library(INLA)
library(sf)

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

read_svi = function(data_loc, year){
  svi_data = read_csv((paste0(data_loc, sprintf("svi_%d_by_DMA.csv",year))),skip = 2)
  colnames(svi_data) = c("Metro",
                         "RPL_THEME1_mean", "RPL_THEME1_median",
                         "RPL_THEME2_mean", "RPL_THEME2_median",
                         "RPL_THEME3_mean", "RPL_THEME3_median",
                         "RPL_THEME4_mean", "RPL_THEME4_median",
                         "RPL_THEMES_mean", "RPL_THEMES_median",
                         "EP_NOINT_mean", "EP_NOINT_median",
                         "E_NOINT_mean", "E_NOINT_median"
  )
  svi_data = svi_data %>%
    mutate(Metro = trimws(gsub("\\."," ",Metro))) %>%
    mutate(Metro = case_when(
      Metro == "Chico Redding CA" ~ "Chico-Redding CA",
      Metro == "Fresno Visalia CA" ~ "Fresno-Visalia CA",
      Metro == "Monterey Salinas CA" ~ "Monterey-Salinas CA",
      Metro == "Sacramento Stockton Modesto CA" ~ "Sacramento-Stockton-Modesto CA",
      Metro == "San Francisco Oakland San Jose CA" ~ "San Francisco-Oakland-San Jose CA",
      Metro == "Santa Barbara Santa Maria San Luis Obispo CA" ~ "Santa Barbara-Santa Maria-San Luis Obispo CA",
      Metro == "Yuma AZ El Centro CA" ~ "Yuma AZ-El Centro CA",
      TRUE ~ Metro
    ))
  return(svi_data)
}

get_sample = function(name,data){
  return(lapply(data, function(s) s$latent[paste0(name,":1"), 1]))
}

get_sum_sample = function(name_lst,data){
  sample_lst =lapply(name_lst, 
                     get_sample,data = data)
  return(Reduce(`+`, sample_lst))
}

summary_samples = function(samp){
  return_vec = c(mean(samp), 
                 quantile(samp,c(0.025,0.975))
  )
  names(return_vec)[1] = "mean"
  return(
    return_vec
  )
}

summary_sample_lst = function(sample_lst){
  summary_lst = lapply(sample_lst, summary_samples)
  return(do.call(rbind, summary_lst))
}


fit_dlag_mod = function(s_t_data,adj_mat,outcome,exposure,
                        time_dependent_interacts,
                        time_independent_interacts,
                        nlags = 4){
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
  
  f = paste0(outcome, " ~ 1 + ", exposure_lags, " + ",
             paste0(time_dependent_interacts_lst, collapse = " + "), " + ",
             paste0(time_independent_interacts_lst, collapse = " + "),
             " + f(Metro_num, model='bym',graph=adj_mat) + f(date_num_str, model='rw2') + f(date_num_unstr, model='iid')"
             )
  model_formula = as.formula(f)
  mod = inla(model_formula, family="poisson", control.compute = list(dic=T,config = T),
             quantiles=c(0.025, 0.5,0.975), data=s_t_data
  )
  samp = inla.posterior.sample(1000, mod)
  
  intercept_effect = lapply(samp, function(s) s$latent["(Intercept):1", 1])
  intercept_effect = do.call(cbind,intercept_effect)
  
  loc_random_effect = lapply(samp, function(s) s$latent[paste0("Metro_num:",seq(1:12)), 1])
  
  loc_random_effect = do.call(cbind,loc_random_effect)
  rownames(loc_random_effect) = unique(s_t_data$Metro)
  
  date_str_effect = lapply(samp, function(s) s$latent[paste0("date_num_str:",seq(1:261)), 1])
  date_str_effect = do.call(cbind,date_str_effect)
  rownames(date_str_effect) = paste0("str:",unique(s_t_data$week))
  
  date_unstr_effect = lapply(samp, function(s) s$latent[paste0("date_num_unstr:",seq(1:261)), 1])
  date_unstr_effect = do.call(cbind,date_unstr_effect)
  rownames(date_unstr_effect) = paste0("unstr:",unique(s_t_data$week))
  
  time_independent_samples = lapply(time_independent_posterior_name,
                                    get_sample,data = samp)
 
  time_independent_samples = lapply(time_independent_samples,
                                    function(s) do.call(cbind,s))
  time_independent_samples = do.call(rbind,time_independent_samples)
  
  
  time_dependent_samples = lapply(time_dependent_posterior_name,
                                  get_sample,data = samp)

  time_dependent_samples = lapply(time_dependent_samples,
                                  function(s) do.call(cbind,s))
  time_dependent_samples = do.call(rbind,time_dependent_samples)
  
  fixed_effect = mod$summary.fixed
  fixed_effect = tibble::rownames_to_column(fixed_effect, "Coefficient") %>%
    filter(Coefficient != "(Intercept)") %>%
    dplyr::select(Coefficient, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`) %>%
    mutate(outcome = outcome,
           exposure = exposure)
  time_effect = exp(mod$summary.random$date_num_str[,"mean"] + mod$summary.random$date_num_unstr[,"mean"])
  time_effect_df = data.frame(
    week = as.Date(unique(s_t_data$week)),
    search_rate = time_effect
  )%>%
    mutate(outcome = outcome,
           exposure = exposure)
  
  region_effect = exp(mod$summary.random$Metro_num[1:length(unique(s_t_data$Metro)),"mean"])
  region_effect_df = data.frame(
    DMA = unique(s_t_data$Metro),
    region_effect = region_effect
  )%>%
    mutate(outcome = outcome,
           exposure = exposure)
  
  return(list("posterior_sample" = rbind(time_dependent_samples,
                                         time_independent_samples,
                                         intercept_effect,
                                         loc_random_effect,
                                         date_str_effect,
                                         date_unstr_effect),
              "fixed_effect" = fixed_effect,
              "time_effect" = time_effect_df,
              "region_effect" = region_effect_df))
}



data_loc = "../Data/"
res_loc = "../Result/"
master.data = read_csv(paste0(data_loc,"master_data.csv"))

spatio_tempo_data = master.data %>%
  select(-incident_names) %>%
  mutate(
    incident_acres_burned = case_when(
      is.na(incident_acres_burned) ~ 0,
      TRUE ~ incident_acres_burned
    ),
    incident_counts = case_when(
      is.na(incident_counts) ~ 0,
      TRUE ~ incident_counts
    )
  )

spatio_tempo_data = spatio_tempo_data %>%
  mutate(Metro = factor(Metro)) %>%
  mutate_at(vars(starts_with("PM25")), scale2) %>%
  mutate_at(vars(starts_with("smokePM")), scale2) %>%
  mutate_at(
    c("temperature", "median_household_income", "population",
      "incident_acres_burned", "incident_counts", "MaxTemp", "MinTemp",
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
         -MinHumidity, 
         -MaxWindSpeed, -MinWindSpeed, -MaxPressure,
         -MinPressure, -MaxDewPoint, -MinDewPoint,
         -RPL_THEME1_median, -RPL_THEME2_median, -RPL_THEME3_median,
         -RPL_THEME4_median, -RPL_THEMES_median, -EP_NOINT_median, -E_NOINT_median
  ) %>%
  rename(fireScale = "incident_acres_burned",
         air_filter = "air filter",
         air_pollution = "air pollution",
         air_purifier = "air purifier",
         air_quality = "air quality") %>%
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


adj_matrix = matrix(0,nrow = 12, ncol = 12)
rownames(adj_matrix) = levels(spatio_tempo_data$Metro)
colnames(adj_matrix) = levels(spatio_tempo_data$Metro)

# Eureka CA
adj_matrix['Eureka, CA','San Francisco-Oakland-San Jose, CA'] = 1.0
adj_matrix['Eureka, CA','Chico-Redding, CA'] = 1.0
# Chico-Redding CA
adj_matrix['Chico-Redding, CA','Eureka, CA'] = 1.0
adj_matrix['Chico-Redding, CA','San Francisco-Oakland-San Jose, CA'] = 1.0
adj_matrix['Chico-Redding, CA','Sacramento-Stockton-Modesto, CA'] = 1.0
# San Francisco-Oakland-San Jose CA
adj_matrix['San Francisco-Oakland-San Jose, CA','Eureka, CA'] = 1.0
adj_matrix['San Francisco-Oakland-San Jose, CA','Chico-Redding, CA'] = 1.0
adj_matrix['San Francisco-Oakland-San Jose, CA','Sacramento-Stockton-Modesto, CA'] = 1.0
adj_matrix['San Francisco-Oakland-San Jose, CA','Monterey-Salinas, CA'] = 1.0
adj_matrix['San Francisco-Oakland-San Jose, CA','Fresno-Visalia, CA'] = 1.0
# Sacramento-Stockton-Modesto CA
adj_matrix['Sacramento-Stockton-Modesto, CA','Chico-Redding, CA'] = 1.0
adj_matrix['Sacramento-Stockton-Modesto, CA','San Francisco-Oakland-San Jose, CA'] = 1.0
adj_matrix['Sacramento-Stockton-Modesto, CA','Fresno-Visalia, CA'] = 1.0
# Monterey-Salinas CA
adj_matrix['Monterey-Salinas, CA','San Francisco-Oakland-San Jose, CA'] = 1.0
adj_matrix['Monterey-Salinas, CA','Fresno-Visalia, CA'] = 1.0
adj_matrix['Monterey-Salinas, CA','Bakersfield, CA'] = 1.0
adj_matrix['Monterey-Salinas, CA','Santa Barbara-Santa Maria-San Luis Obispo, CA'] = 1.0
# Fresno-Visalia CA
adj_matrix['Fresno-Visalia, CA','Monterey-Salinas, CA'] = 1.0
adj_matrix['Fresno-Visalia, CA','Santa Barbara-Santa Maria-San Luis Obispo, CA'] = 1.0
adj_matrix['Fresno-Visalia, CA','Bakersfield, CA'] = 1.0
adj_matrix['Fresno-Visalia, CA','Sacramento-Stockton-Modesto, CA'] = 1.0
adj_matrix['Fresno-Visalia, CA','San Francisco-Oakland-San Jose, CA'] = 1.0
adj_matrix['Fresno-Visalia, CA','Los Angeles, CA'] = 1.0
# Santa Barbara-Santa Maria-San Luis Obispo CA
adj_matrix['Santa Barbara-Santa Maria-San Luis Obispo, CA','Monterey-Salinas, CA'] = 1.0
adj_matrix['Santa Barbara-Santa Maria-San Luis Obispo, CA','Fresno-Visalia, CA'] = 1.0
adj_matrix['Santa Barbara-Santa Maria-San Luis Obispo, CA','Los Angeles, CA'] = 1.0
adj_matrix['Santa Barbara-Santa Maria-San Luis Obispo, CA','Bakersfield, CA'] = 1.0
# Bakersfield CA
adj_matrix['Bakersfield, CA','Monterey-Salinas, CA'] = 1.0
adj_matrix['Bakersfield, CA','Fresno-Visalia, CA'] = 1.0
adj_matrix['Bakersfield, CA','Santa Barbara-Santa Maria-San Luis Obispo, CA'] = 1.0
adj_matrix['Bakersfield, CA','Los Angeles, CA'] = 1.0
# Los Angeles CA
adj_matrix['Los Angeles, CA','Santa Barbara-Santa Maria-San Luis Obispo, CA'] = 1.0
adj_matrix['Los Angeles, CA','Bakersfield, CA'] = 1.0
adj_matrix['Los Angeles, CA','Fresno-Visalia, CA'] = 1.0
adj_matrix['Los Angeles, CA','Palm Springs, CA'] = 1.0
adj_matrix['Los Angeles, CA','San Diego, CA'] = 1.0
adj_matrix['Los Angeles, CA','Yuma, AZ-El Centro, CA'] = 1.0
# Palm Springs CA
adj_matrix['Palm Springs, CA','Los Angeles, CA'] = 1.0
adj_matrix['Palm Springs, CA','San Diego, CA'] = 1.0
adj_matrix['Palm Springs, CA','Yuma, AZ-El Centro, CA'] = 1.0
# San Diego CA
adj_matrix['San Diego, CA','Los Angeles, CA'] = 1.0
adj_matrix['San Diego, CA','Palm Springs, CA'] = 1.0
adj_matrix['San Diego, CA','Yuma, AZ-El Centro, CA'] = 1.0
# Yuma AZ-El Centro CA
adj_matrix['Yuma, AZ-El Centro, CA','Los Angeles, CA'] = 1.0
adj_matrix['Yuma, AZ-El Centro, CA','San Diego, CA'] = 1.0
adj_matrix['Yuma, AZ-El Centro, CA','Palm Springs, CA'] = 1.0
adj_matrix_new = t(apply(adj_matrix,1, function(x) x/sum(x)))

spatio_tempo_data = spatio_tempo_data %>%
  mutate(Metro = as.character(Metro),
         week = as.character(week))

spatio_tempo_data = num_date_region(spatio_tempo_data)




exposure_vec = c("smokePM")
outcome_vec = c("air_pollution", "air_quality", "air_filter", "air_purifier")
time_dependent_vec = c("AvgTemp", "AvgHumidity", "AvgWindSpeed",
                       "AvgPressure","fireScale")
time_independent_vec = c("RPL_THEME1_mean", "RPL_THEME2_mean",
                         "RPL_THEME3_mean", "RPL_THEME4_mean")
nlags = 3
sample_lst = list()
fixed_effect_lst = list()
time_effect_lst = list()
region_effect_lst = list()
for (outcome in outcome_vec) {
  for (exposure in exposure_vec) {
    mod_data = spatio_tempo_data
    adj_mat_data = adj_matrix_new
    res = fit_dlag_mod(s_t_data = mod_data,
                       adj_mat = adj_mat_data,
                       outcome = outcome, exposure = exposure,
                       time_dependent_interacts = time_dependent_vec,
                       time_independent_interacts = time_independent_vec,
                       nlags = nlags)
    sample_lst[[paste0(outcome,'_',exposure)]] = res[["posterior_sample"]]
    fixed_effect_lst[[paste0(outcome,'_',exposure)]] = res[["fixed_effect"]]
    time_effect_lst[[paste0(outcome,'_',exposure)]] = res[["time_effect"]]
    region_effect_lst[[paste0(outcome,'_',exposure)]] = res[["region_effect"]]
    print(paste0(outcome,'_',exposure))
  }
}
save(sample_lst, file = "posterior_coefficients_dlag.RData")
save(fixed_effect_lst, file = "fixed_effect_dlag.RData")
save(time_effect_lst, file = "time_effect_dlag.RData")
save(region_effect_lst, file = "region_effect_dlag.RData")
