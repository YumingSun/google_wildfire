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
library(latex2exp)
library(xtable)

create_coef_for_ci = function(data_lst, data_name, confounder_name,
                              time_dependent = T){
  if (time_dependent) {
    intercept = data_lst[[data_name]][c("smokePM:1","smokePM_lag1:1",
                                        "smokePM_lag2:1","smokePM_lag3:1"),]
    confounder = data_lst[[data_name]][c(sprintf("smokePM:%s:1",confounder_name),
                                         sprintf("smokePM_lag1:%s_lag1:1",confounder_name),
                                         sprintf("smokePM_lag2:%s_lag2:1",confounder_name),
                                         sprintf("smokePM_lag3:%s_lag3:1",confounder_name)),
    ]
  } else {
    intercept = data_lst[[data_name]][c("smokePM:1","smokePM_lag1:1",
                                        "smokePM_lag2:1","smokePM_lag3:1"),]
    confounder = data_lst[[data_name]][c(sprintf("smokePM:%s:1",confounder_name),
                                         sprintf("smokePM_lag1:%s:1",confounder_name),
                                         sprintf("smokePM_lag2:%s:1",confounder_name),
                                         sprintf("smokePM_lag3:%s:1",confounder_name)),
    ]
  }
  
  coef = rbind(intercept,confounder)
  return(coef)
}

mean_ci_from_post = function(data_vec){
  avg = mean(data_vec)
  ci = quantile(data_vec,p = c(0.025, 0.975))
  res_num = c(avg, ci)
  res_out = sprintf("%.3f (%.3f, %.3f)", 
                    res_num[1],res_num[2],res_num[3])
  return(res_out)
}

create_data_for_ci = function(var, population_data, lb = 0.1, ub = 0.9){
  df = data.frame(
    v_0 = seq(from = quantile(population_data[[var]],lb),
              to = quantile(population_data[[var]],ub),
              length.out = 103)
  ) %>%
    mutate(v_1 = lag(v_0, n = 1, default = NA),
           v_2 = lag(v_0, n = 2, default = NA),
           v_3 = lag(v_0, n = 3, default = NA)) %>%
    drop_na() %>%
    mutate(v_0_scale = (v_0 - mean(population_data[[var]]))/sd(population_data[[var]]),
           v_1_scale = (v_1 - mean(population_data[[var]]))/sd(population_data[[var]]),
           v_2_scale = (v_2 - mean(population_data[[var]]))/sd(population_data[[var]]),
           v_3_scale = (v_3 - mean(population_data[[var]]))/sd(population_data[[var]]),
           i_0 = rep(1.0, 100), i_1 = rep(1.0, 100), i_2 = rep(1.0, 100),
           i_3 = rep(1.0, 100))
  
  return(df)
}

calculate_rrsi = function(df, coef, var_name){
  rrsi_sample = data.matrix(df[,c("i_0","i_1","i_2","i_3", 
                                  "v_0_scale","v_1_scale",
                                  "v_2_scale","v_3_scale")]) %*%coef
  rrsi_sample = exp(rrsi_sample)
  rrsi_sample_summary = apply(rrsi_sample, 1, quantile,c(0.025,0.5,0.975))
  
  rrsi_df = data.frame(cbind(df$v_0, t(rrsi_sample_summary)))
  rrsi_df[['var_name']] = var_name
  colnames(rrsi_df) = c("var_value", "2.5%CI", "50%CI", "97.5%CI","var_name")
  return(rrsi_df)
  
}

data_loc = "../Data/"
img_loc = "../Results/"
master.data = read_csv(paste0(data_loc,"master_data.csv"))
dma_map = st_read(paste0(data_loc,"NatDMA/NatDMA.shp"))
dma_map$NAME = gsub("-", " ", dma_map$NAME)
colnames(dma_map)[4] = "DMA"

load("posterior_coefficients_dlag.RData")
load("fixed_effect_dlag.RData")
load("time_effect_dlag.RData")
load("region_effect_dlag.RData")
# Descriptive Analysis
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

spatio_tempo_data$DMA = gsub("-", " ", spatio_tempo_data$Metro)
spatio_tempo_data$DMA = gsub(",", "", spatio_tempo_data$DMA)

wild_fire_temporal = spatio_tempo_data %>%
  group_by(week) %>%
  summarise(incident_counts_sum = sum(incident_counts),
            acre_burned_sum = sum(incident_acres_burned),
            smokePM_mean = mean(smokePM),
            air_filter_mean = mean(`air filter`),
            air_purifier_mean = mean(`air purifier`),
            air_pollution_mean = mean(`air pollution`),
            air_quality_mean = mean(`air quality`)) %>%
  ungroup()

google_search_temporal_long = wild_fire_temporal %>%
  pivot_longer(cols = c("air_purifier_mean","air_pollution_mean" ),
               names_to = "Search Terms Original",
               values_to = "RSI(%)") %>%
  mutate(`Search Terms` = factor(`Search Terms Original`, labels = c("Air pollution","Air purifier")))

smokePM_temporal_long = wild_fire_temporal %>%
  pivot_longer(cols = c("smokePM_mean"),
               names_to = "Average Smoke PM2.5",
               values_to = "RSI(%)") 
smokePM_temporal = ggplot(wild_fire_temporal) + 
  geom_line(aes(x = week, y = smokePM_mean),linewidth=0.8) +
  scale_x_date(date_labels = "%b-%Y", date_breaks  ="3 month") +
  xlab("Time")+
  ylab(TeX("\\textbf{Smoke PM2.5} ($\\mu g/m^3$)")) + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,face = "bold",
                                   angle = 90, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 10,face = "bold"))

search_temporal = ggplot(google_search_temporal_long) + 
  geom_line(aes(x = week, y = `RSI(%)`,color = `Search Terms`),linewidth=0.8) +
  scale_x_date(date_labels = "%b-%Y", date_breaks  ="3 month") +
  xlab("Time")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,face = "bold",
                                   angle = 90, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 10,face = "bold"))


fire_scale = ggplot(wild_fire_temporal) + 
  geom_bar(aes(x = week, y= acre_burned_sum/(1E5)),stat="identity") +
  scale_x_date(date_labels = "%b-%Y", date_breaks  ="3 month") +
  ylab(TeX("\\textbf{Acres Burned} ($\\times 10^5$ acres)"))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.x  = element_blank(),
        axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank())

fire_counts = ggplot(wild_fire_temporal) + 
  geom_bar(aes(x = week, y= incident_counts_sum),stat="identity") +
  scale_x_date(date_labels = "%b-%Y", date_breaks  ="3 month") +
  ylab(TeX("\\textbf{Wildfire Counts}"))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.x  = element_blank(),
        axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank())

p1 = ggarrange(fire_counts, search_temporal, heights = c(0.7, 2),
               ncol = 1, nrow = 2,align = "v")
p2 = ggarrange(fire_scale, smokePM_temporal, heights = c(0.7, 2),
               ncol = 1, nrow = 2,align = "v")

wild_fire_loc = spatio_tempo_data %>%
  group_by(Metro) %>%
  summarise(incident_counts_sum = sum(incident_counts),
            acre_burned_sum = sum(incident_acres_burned),
            smokePM_mean = mean(smokePM),
            air_filter_mean = mean(`air filter`),
            air_purifier_mean = mean(`air purifier`),
            air_pollution_mean = mean(`air pollution`),
            air_quality_mean = mean(`air quality`)) 
wild_fire_loc$DMA = gsub("-", " ", wild_fire_loc$Metro)
wild_fire_loc$DMA = gsub(",", "", wild_fire_loc$DMA)


plot_google_search = dma_map %>%
  inner_join(wild_fire_loc, by = "DMA") %>%
  pivot_longer(
    cols = c("air_purifier_mean","air_pollution_mean"
    ),
    names_to = "Category",
    values_to = "RSI (%)"
  ) %>% 
  mutate(
    search_terms = case_when(
      Category == "air_purifier_mean" ~ "Air purifier",
      Category == "air_pollution_mean" ~ "Air pollution",
    )
  )

plot_smoke = dma_map %>%
  inner_join(wild_fire_loc, by = "DMA")


google_search_loc = ggplot(data = plot_google_search) +
  geom_sf(aes(fill = `RSI (%)`)) +
  facet_wrap(~ search_terms, ncol = 2) + 
  scale_fill_viridis_c() +
  theme_minimal() +
  theme_void() +
  theme(legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold",size=10),
        strip.text = element_text(face="bold",size=10))

smoke_loc = ggplot(data = plot_smoke) +
  geom_sf(aes(fill = smokePM_mean)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme_void() + labs(fill=TeX('\\textbf{Smoke PM}$_{2.5}$',bold=TRUE)) + 
  theme(legend.title = element_text(face = "bold"))

p3 = ggarrange(google_search_loc, smoke_loc, widths = c(2, 1.2),
               ncol = 2, nrow = 1)

# main effects
fixed_effect_df = bind_rows(fixed_effect_lst, .id = "column_label")

main_effect_df = fixed_effect_df %>%
  filter(Coefficient %in% c(
    "pm25", "pm25_lag1", "pm25_lag2", "pm25_lag3", "pm25_lag4",
    "smokePM", "smokePM_lag1", "smokePM_lag2", "smokePM_lag3", "smokePM_lag4"
  )) %>%
  mutate(
    lag = case_when(
      is.na(word(Coefficient, 2, sep = "_")) ~ "lag0",
      TRUE ~ word(Coefficient, 2, sep = "_")
    ),
    lag_num = case_when(
      lag == "lag0" ~ 0,
      lag == "lag1" ~ 1,
      lag == "lag2" ~ 2,
      lag == "lag3" ~ 3,
    )
  ) %>%
  mutate(mean = exp(mean),`0.025quant` = exp(`0.025quant`),
         `0.975quant` = exp(`0.975quant`))

p4 = main_effect_df %>%
  filter(outcome %in% c("air_pollution", "air_purifier")) %>%
  mutate(outcome = factor(outcome, labels = c("Air pollution", "Air purifier"))) %>%
  ggplot(aes(x = lag_num, y = mean)) +
  geom_pointrange(aes(ymin=`0.025quant`, ymax = `0.975quant`)) +
  geom_line() + 
  xlab("Lag (week)") +
  ylab("Relative RSI") + 
  facet_wrap(~ outcome, ncol = 2)+#,  scales = "free") + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=15),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x = element_text(size = 15,face = "bold"))

appendix_p1 = main_effect_df %>%
  filter(outcome %in% c("air_quality", "air_filter")) %>%
  mutate(outcome = factor(outcome, labels = c("Air quality", "Air filter"))) %>%
  ggplot(aes(x = lag_num, y = mean)) +
  geom_pointrange(aes(ymin=`0.025quant`, ymax = `0.975quant`)) +
  geom_line() + 
  xlab("Lag (week)") +
  ylab("Relative RSI") + 
  facet_wrap(~ outcome, ncol = 2)+#,  scales = "free") + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=15),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x = element_text(size = 15,face = "bold"))


pollution_temp_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                         "AvgTemp")
pollution_humid_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                          "AvgHumidity")
pollution_wind_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                         "AvgWindSpeed")
pollution_pressure_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                             "AvgPressure")
pollution_fire_scale_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                               "fireScale")
pollution_socio_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                          "RPL_THEME1_mean", time_dependent = F)
pollution_house_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                          "RPL_THEME2_mean", time_dependent = F)
pollution_race_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                         "RPL_THEME3_mean", time_dependent = F)
pollution_trans_coef = create_coef_for_ci(sample_lst,"air_pollution_smokePM",
                                          "RPL_THEME4_mean", time_dependent = F)


purifier_temp_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                        "AvgTemp")
purifier_humid_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                         "AvgHumidity")
purifier_wind_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                        "AvgWindSpeed")
purifier_pressure_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                            "AvgPressure")
purifier_fire_scale_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                              "fireScale")
purifier_socio_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                         "RPL_THEME1_mean", time_dependent = F)
purifier_house_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                         "RPL_THEME2_mean", time_dependent = F)
purifier_race_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                        "RPL_THEME3_mean", time_dependent = F)
purifier_trans_coef = create_coef_for_ci(sample_lst,"air_purifier_smokePM",
                                         "RPL_THEME4_mean", time_dependent = F)

quality_temp_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                       "AvgTemp")

quality_humid_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                        "AvgHumidity")

quality_wind_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                       "AvgWindSpeed")

quality_pressure_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                           "AvgPressure")
quality_fire_scale_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                             "fireScale")
quality_socio_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                        "RPL_THEME1_mean", time_dependent = F)
quality_house_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                        "RPL_THEME2_mean", time_dependent = F)
quality_race_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                       "RPL_THEME3_mean", time_dependent = F)
quality_trans_coef = create_coef_for_ci(sample_lst,"air_quality_smokePM",
                                        "RPL_THEME4_mean", time_dependent = F)


filter_temp_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                      "AvgTemp")
filter_humid_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                       "AvgHumidity")
filter_wind_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                      "AvgWindSpeed")
filter_pressure_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                          "AvgPressure")
filter_fire_scale_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                            "fireScale")
filter_socio_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                       "RPL_THEME1_mean", time_dependent = F)
filter_house_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                       "RPL_THEME2_mean", time_dependent = F)
filter_race_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                      "RPL_THEME3_mean", time_dependent = F)
filter_trans_coef = create_coef_for_ci(sample_lst,"air_filter_smokePM",
                                       "RPL_THEME4_mean", time_dependent = F)

temp_df = create_data_for_ci("AvgTemp", spatio_tempo_data)
humid_df = create_data_for_ci("AvgHumidity", spatio_tempo_data)
wind_df = create_data_for_ci("AvgWindSpeed", spatio_tempo_data)
pressure_df = create_data_for_ci("AvgPressure", spatio_tempo_data)
fire_scale_df = create_data_for_ci("incident_acres_burned", spatio_tempo_data,
                                   lb = 0.9, ub = 0.99)
socio_df = create_data_for_ci("RPL_THEME1_mean", spatio_tempo_data)
house_df = create_data_for_ci("RPL_THEME2_mean", spatio_tempo_data)
race_df = create_data_for_ci("RPL_THEME3_mean", spatio_tempo_data)
trans_df = create_data_for_ci("RPL_THEME4_mean", spatio_tempo_data)


pollution_temp_rrsi = calculate_rrsi(temp_df,pollution_temp_coef,"Temperature (°F)")
pollution_humid_rrsi = calculate_rrsi(humid_df,pollution_humid_coef,"Humidity (%)")
pollution_wind_rrsi = calculate_rrsi(wind_df,pollution_wind_coef,"Wind Speed (mph)")
pollution_pressure_rrsi = calculate_rrsi(pressure_df,pollution_pressure_coef,
                                         "Atmospheric Pressure (inHg)")
pollution_fire_scale_rrsi = calculate_rrsi(fire_scale_df,pollution_fire_scale_coef,
                                           "Fire Scale (acres)")
pollution_socio_rrsi = calculate_rrsi(socio_df,pollution_socio_coef,
                                      "Socioeconomic Status")
pollution_house_rrsi = calculate_rrsi(house_df,pollution_house_coef,
                                      "Household Characteristics")
pollution_race_rrsi = calculate_rrsi(race_df,pollution_race_coef,
                                     "Racial & Ethnic Minority Status")
pollution_trans_rrsi = calculate_rrsi(trans_df,pollution_trans_coef,
                                      "Housing Type & Transportation")

purifier_temp_rrsi = calculate_rrsi(temp_df,purifier_temp_coef,"Temperature (°F)")
purifier_humid_rrsi = calculate_rrsi(humid_df,purifier_humid_coef,"Humidity (%)")
purifier_wind_rrsi = calculate_rrsi(wind_df,purifier_wind_coef,"Wind Speed (mph)")
purifier_pressure_rrsi = calculate_rrsi(pressure_df,purifier_pressure_coef,
                                        "Atmospheric Pressure (inHg)")
purifier_fire_scale_rrsi = calculate_rrsi(fire_scale_df,purifier_fire_scale_coef,
                                          "Fire Scale (acres)")
purifier_socio_rrsi = calculate_rrsi(socio_df,purifier_socio_coef,
                                     "Socioeconomic Status")
purifier_house_rrsi = calculate_rrsi(house_df,purifier_house_coef,
                                     "Household Characteristics")
purifier_race_rrsi = calculate_rrsi(race_df,purifier_race_coef,
                                    "Racial & Ethnic Minority Status")
purifier_trans_rrsi = calculate_rrsi(trans_df,purifier_trans_coef,
                                     "Housing Type & Transportation")

pollution_rrsi = rbind.data.frame(pollution_temp_rrsi, pollution_humid_rrsi,
                                  pollution_wind_rrsi, pollution_pressure_rrsi,
                                  pollution_fire_scale_rrsi,
                                  pollution_socio_rrsi, pollution_house_rrsi,
                                  pollution_race_rrsi, pollution_trans_rrsi) %>%
  mutate(var_name = factor(var_name, 
                           levels = c("Fire Scale (acres)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))


purifier_rrsi = rbind.data.frame(purifier_temp_rrsi, purifier_humid_rrsi,
                                 purifier_wind_rrsi, purifier_pressure_rrsi,
                                 purifier_fire_scale_rrsi,
                                 purifier_socio_rrsi, purifier_house_rrsi,
                                 purifier_race_rrsi, purifier_trans_rrsi) %>%
  mutate(var_name = factor(var_name, 
                           levels = c("Fire Scale (acres)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))

quality_temp_rrsi = calculate_rrsi(temp_df,quality_temp_coef,"Temperature (°F)")
quality_humid_rrsi = calculate_rrsi(humid_df,quality_humid_coef,"Humidity (%)")
quality_wind_rrsi = calculate_rrsi(wind_df,quality_wind_coef,"Wind Speed (mph)")
quality_pressure_rrsi = calculate_rrsi(pressure_df,quality_pressure_coef,
                                       "Atmospheric Pressure (inHg)")
quality_fire_scale_rrsi = calculate_rrsi(fire_scale_df,quality_fire_scale_coef,
                                         "Fire Scale (acres)")
quality_socio_rrsi = calculate_rrsi(socio_df,quality_socio_coef,
                                    "Socioeconomic Status")
quality_house_rrsi = calculate_rrsi(house_df,quality_house_coef,
                                    "Household Characteristics")
quality_race_rrsi = calculate_rrsi(race_df,quality_race_coef,
                                   "Racial & Ethnic Minority Status")
quality_trans_rrsi = calculate_rrsi(trans_df,quality_trans_coef,
                                    "Housing Type & Transportation")

filter_temp_rrsi = calculate_rrsi(temp_df,filter_temp_coef,"Temperature (°F)")
filter_humid_rrsi = calculate_rrsi(humid_df,filter_humid_coef,"Humidity (%)")
filter_wind_rrsi = calculate_rrsi(wind_df,filter_wind_coef,"Wind Speed (mph)")
filter_pressure_rrsi = calculate_rrsi(pressure_df,filter_pressure_coef,
                                      "Atmospheric Pressure (inHg)")
filter_fire_scale_rrsi = calculate_rrsi(fire_scale_df,filter_fire_scale_coef,
                                        "Fire Scale (acres)")
filter_socio_rrsi = calculate_rrsi(socio_df,filter_socio_coef,
                                   "Socioeconomic Status")
filter_house_rrsi = calculate_rrsi(house_df,filter_house_coef,
                                   "Household Characteristics")
filter_race_rrsi = calculate_rrsi(race_df,filter_race_coef,
                                  "Racial & Ethnic Minority Status")
filter_trans_rrsi = calculate_rrsi(trans_df,filter_trans_coef,
                                   "Housing Type & Transportation")

quality_rrsi = rbind.data.frame(quality_temp_rrsi, quality_humid_rrsi,
                                quality_wind_rrsi, quality_pressure_rrsi,
                                quality_fire_scale_rrsi,
                                quality_socio_rrsi, quality_house_rrsi,
                                quality_race_rrsi, quality_trans_rrsi) %>%
  mutate(var_name = factor(var_name, 
                           levels = c("Fire Scale (acres)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))


filter_rrsi = rbind.data.frame(filter_temp_rrsi, filter_humid_rrsi,
                               filter_wind_rrsi, filter_pressure_rrsi,
                               filter_fire_scale_rrsi,
                               filter_socio_rrsi, filter_house_rrsi,
                               filter_race_rrsi, filter_trans_rrsi) %>%
  mutate(var_name = factor(var_name, 
                           levels = c("Fire Scale (acres)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))

p5 = ggplot(data = pollution_rrsi, aes(x = var_value)) + 
  geom_line(aes(y = `50%CI`)) +
  geom_ribbon(aes(ymin=`2.5%CI`, ymax=`97.5%CI`), linetype=2, alpha=0.1) +
  facet_wrap(~ var_name, ncol = 3, scales = "free") + 
  xlab("") +
  ylab("cumulative Relative RSI") + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=12),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.title.x = element_text(size = 12,face = "bold"))

p6 = ggplot(data = purifier_rrsi, aes(x = var_value)) + 
  geom_line(aes(y = `50%CI`)) +
  geom_ribbon(aes(ymin=`2.5%CI`, ymax=`97.5%CI`), linetype=2, alpha=0.1) +
  facet_wrap(~ var_name, ncol = 3, scales = "free") + 
  xlab("") +
  ylab("cumulative Relative RSI") + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=12),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.title.x = element_text(size = 12,face = "bold"))

loc_data = spatio_tempo_data %>%
  select(week, Metro, AvgTemp, AvgHumidity, AvgWindSpeed, AvgPressure, incident_acres_burned,
         smokePM, RPL_THEME1_mean, RPL_THEME2_mean, RPL_THEME3_mean,
         RPL_THEME4_mean) %>%
  group_by(Metro) %>%
  summarise(AvgTemp = mean(AvgTemp,na.rm = T), 
            AvgHumidity = mean(AvgHumidity,na.rm = T),
            AvgWindSpeed = mean(AvgWindSpeed,na.rm = T), 
            AvgPressure = mean(AvgPressure,na.rm = T),
            incident_acres_burned = mean(incident_acres_burned,na.rm = T), 
            smokePM = mean(smokePM,na.rm = T),
            RPL_THEME1_mean = mean(RPL_THEME1_mean,na.rm = T), 
            RPL_THEME2_mean = mean(RPL_THEME2_mean,na.rm = T),
            RPL_THEME3_mean = mean(RPL_THEME3_mean,na.rm = T),
            RPL_THEME4_mean = mean(RPL_THEME4_mean,na.rm = T)
  ) %>%
  mutate(AvgTemp_scale = (AvgTemp - mean(spatio_tempo_data$AvgTemp,na.rm = T))/sd(spatio_tempo_data$AvgTemp,na.rm = T), 
         AvgTemp_scale_lag1 = AvgTemp_scale,
         AvgTemp_scale_lag2 = AvgTemp_scale,
         AvgTemp_scale_lag3 = AvgTemp_scale,
         AvgHumidity_scale = (AvgHumidity - mean(spatio_tempo_data$AvgHumidity,na.rm = T))/sd(spatio_tempo_data$AvgHumidity,na.rm = T), 
         AvgHumidity_scale_lag1 = AvgHumidity_scale,
         AvgHumidity_scale_lag2 = AvgHumidity_scale,
         AvgHumidity_scale_lag3 = AvgHumidity_scale,
         AvgWindSpeed_scale = (AvgWindSpeed - mean(spatio_tempo_data$AvgWindSpeed,na.rm = T))/sd(spatio_tempo_data$AvgWindSpeed,na.rm = T),
         AvgWindSpeed_scale_lag1 = AvgWindSpeed_scale,
         AvgWindSpeed_scale_lag2 = AvgWindSpeed_scale,
         AvgWindSpeed_scale_lag3 = AvgWindSpeed_scale,
         AvgPressure_scale = (AvgPressure - mean(spatio_tempo_data$AvgPressure,na.rm = T))/sd(spatio_tempo_data$AvgPressure,na.rm = T), 
         AvgPressure_scale_lag1 = AvgPressure_scale,
         AvgPressure_scale_lag2 = AvgPressure_scale,
         AvgPressure_scale_lag3 = AvgPressure_scale,
         incident_acres_burned_scale = (incident_acres_burned - mean(spatio_tempo_data$incident_acres_burned,na.rm = T))/sd(spatio_tempo_data$incident_acres_burned,na.rm = T),
         incident_acres_burned_scale_lag1 = incident_acres_burned_scale,
         incident_acres_burned_scale_lag2 = incident_acres_burned_scale,
         incident_acres_burned_scale_lag3 = incident_acres_burned_scale,
         smokePM_scale = (smokePM - mean(spatio_tempo_data$smokePM,na.rm = T))/sd(spatio_tempo_data$smokePM,na.rm = T), 
         smokePM_scale_lag1 = smokePM_scale,
         smokePM_scale_lag2 = smokePM_scale,
         smokePM_scale_lag3 = smokePM_scale,
         RPL_THEME1_mean_scale = (RPL_THEME1_mean - mean(spatio_tempo_data$RPL_THEME1_mean,na.rm = T))/sd(spatio_tempo_data$RPL_THEME1_mean,na.rm = T), 
         RPL_THEME1_mean_scale_lag1 = RPL_THEME1_mean_scale,
         RPL_THEME1_mean_scale_lag2 = RPL_THEME1_mean_scale,
         RPL_THEME1_mean_scale_lag3 = RPL_THEME1_mean_scale,
         RPL_THEME2_mean_scale = (RPL_THEME2_mean - mean(spatio_tempo_data$RPL_THEME2_mean,na.rm = T))/sd(spatio_tempo_data$RPL_THEME2_mean,na.rm = T), 
         RPL_THEME2_mean_scale_lag1 = RPL_THEME2_mean_scale,
         RPL_THEME2_mean_scale_lag2 = RPL_THEME2_mean_scale,
         RPL_THEME2_mean_scale_lag3 = RPL_THEME2_mean_scale,
         RPL_THEME3_mean_scale = (RPL_THEME3_mean - mean(spatio_tempo_data$RPL_THEME3_mean,na.rm = T))/sd(spatio_tempo_data$RPL_THEME3_mean,na.rm = T), 
         RPL_THEME3_mean_scale_lag1 = RPL_THEME3_mean_scale,
         RPL_THEME3_mean_scale_lag2 = RPL_THEME3_mean_scale,
         RPL_THEME3_mean_scale_lag3 = RPL_THEME3_mean_scale,
         RPL_THEME4_mean_scale = (RPL_THEME4_mean - mean(spatio_tempo_data$RPL_THEME4_mean,na.rm = T))/sd(spatio_tempo_data$RPL_THEME4_mean,na.rm = T),
         RPL_THEME4_mean_scale_lag1 = RPL_THEME4_mean_scale,
         RPL_THEME4_mean_scale_lag2 = RPL_THEME4_mean_scale,
         RPL_THEME4_mean_scale_lag3 = RPL_THEME4_mean_scale,
  ) %>%
  rename(
    `smokePM:AvgTemp:1` = AvgTemp_scale,
    `smokePM_lag1:AvgTemp_lag1:1` = AvgTemp_scale_lag1,
    `smokePM_lag2:AvgTemp_lag2:1` = AvgTemp_scale_lag2,
    `smokePM_lag3:AvgTemp_lag3:1` = AvgTemp_scale_lag3,
    `smokePM:AvgHumidity:1` = AvgHumidity_scale,
    `smokePM_lag1:AvgHumidity_lag1:1` = AvgHumidity_scale_lag1,
    `smokePM_lag2:AvgHumidity_lag2:1` = AvgHumidity_scale_lag2,
    `smokePM_lag3:AvgHumidity_lag3:1` = AvgHumidity_scale_lag3,
    `smokePM:AvgWindSpeed:1` = AvgWindSpeed_scale,
    `smokePM_lag1:AvgWindSpeed_lag1:1` = AvgWindSpeed_scale_lag1,
    `smokePM_lag2:AvgWindSpeed_lag2:1` = AvgWindSpeed_scale_lag2,
    `smokePM_lag3:AvgWindSpeed_lag3:1` = AvgWindSpeed_scale_lag3,
    `smokePM:AvgPressure:1` = AvgPressure_scale,
    `smokePM_lag1:AvgPressure_lag1:1` = AvgPressure_scale_lag1,
    `smokePM_lag2:AvgPressure_lag2:1` = AvgPressure_scale_lag2,
    `smokePM_lag3:AvgPressure_lag3:1` = AvgPressure_scale_lag3,
    `smokePM:fireScale:1` = incident_acres_burned_scale,
    `smokePM_lag1:fireScale_lag1:1` = incident_acres_burned_scale_lag1,
    `smokePM_lag2:fireScale_lag2:1` = incident_acres_burned_scale_lag2,
    `smokePM_lag3:fireScale_lag3:1` = incident_acres_burned_scale_lag3,
    `smokePM:RPL_THEME1_mean:1` = RPL_THEME1_mean_scale,
    `smokePM_lag1:RPL_THEME1_mean:1` = RPL_THEME1_mean_scale_lag1,
    `smokePM_lag2:RPL_THEME1_mean:1` = RPL_THEME1_mean_scale_lag2,
    `smokePM_lag3:RPL_THEME1_mean:1` = RPL_THEME1_mean_scale_lag3,
    `smokePM:RPL_THEME2_mean:1` = RPL_THEME2_mean_scale,
    `smokePM_lag1:RPL_THEME2_mean:1` = RPL_THEME2_mean_scale_lag1,
    `smokePM_lag2:RPL_THEME2_mean:1` = RPL_THEME2_mean_scale_lag2,
    `smokePM_lag3:RPL_THEME2_mean:1` = RPL_THEME2_mean_scale_lag3,
    `smokePM:RPL_THEME3_mean:1` = RPL_THEME3_mean_scale,
    `smokePM_lag1:RPL_THEME3_mean:1` = RPL_THEME3_mean_scale_lag1,
    `smokePM_lag2:RPL_THEME3_mean:1` = RPL_THEME3_mean_scale_lag2,
    `smokePM_lag3:RPL_THEME3_mean:1` = RPL_THEME3_mean_scale_lag3,
    `smokePM:RPL_THEME4_mean:1` = RPL_THEME4_mean_scale,
    `smokePM_lag1:RPL_THEME4_mean:1` = RPL_THEME4_mean_scale_lag1,
    `smokePM_lag2:RPL_THEME4_mean:1` = RPL_THEME4_mean_scale_lag2,
    `smokePM_lag3:RPL_THEME4_mean:1` = RPL_THEME4_mean_scale_lag3
  )


loc_data$`smokePM:1` = 1.0
loc_data$`smokePM_lag1:1` = 1.0
loc_data$`smokePM_lag2:1` = 1.0
loc_data$`smokePM_lag3:1` = 1.0

interact_name = rownames(sample_lst$air_pollution_smokePM)[1:40]

pollution_rrsi = data.matrix(loc_data[,interact_name]) %*% data.matrix(sample_lst$air_pollution_smokePM[interact_name,])
rownames(pollution_rrsi) = loc_data$Metro
pollution_rrsi_mean = exp(apply(
  pollution_rrsi,
  1,
  mean
))

pollution_rrsi_df = data.frame(DMA = names(pollution_rrsi_mean), 
                               RRSI = pollution_rrsi_mean, row.names = NULL) 

pollution_rrsi_df$DMA = gsub(",", "", pollution_rrsi_df$DMA)
pollution_rrsi_df$DMA = gsub("-", " ", pollution_rrsi_df$DMA)
pollution_rrsi_map = dma_map %>%
  inner_join(pollution_rrsi_df, by = "DMA") 
pollution_rrsi_map$Search = "Air pollution"

quality_rrsi = data.matrix(loc_data[,interact_name]) %*% data.matrix(sample_lst$air_quality_smokePM[interact_name,])
rownames(quality_rrsi) = loc_data$Metro
quality_rrsi_mean = exp(apply(
  quality_rrsi,
  1,
  mean
))

quality_rrsi_df = data.frame(DMA = names(quality_rrsi_mean), 
                             RRSI = quality_rrsi_mean, row.names = NULL) 

quality_rrsi_df$DMA = gsub(",", "", quality_rrsi_df$DMA)
quality_rrsi_df$DMA = gsub("-", " ", quality_rrsi_df$DMA)
quality_rrsi_map = dma_map %>%
  inner_join(quality_rrsi_df, by = "DMA") 
quality_rrsi_map$Search = "Air quality"

purifier_rrsi = data.matrix(loc_data[,interact_name]) %*% data.matrix(sample_lst$air_purifier_smokePM[interact_name,])
rownames(purifier_rrsi) = loc_data$Metro
purifier_rrsi_mean = exp(apply(
  purifier_rrsi,
  1,
  mean
))

purifier_rrsi_df = data.frame(DMA = names(purifier_rrsi_mean), 
                              RRSI = purifier_rrsi_mean, row.names = NULL) 

purifier_rrsi_df$DMA = gsub(",", "", purifier_rrsi_df$DMA)
purifier_rrsi_df$DMA = gsub("-", " ", purifier_rrsi_df$DMA)

purifier_rrsi_map = dma_map %>%
  inner_join(purifier_rrsi_df, by = "DMA") 

purifier_rrsi_map$Search = "Air purifier"

filter_rrsi = data.matrix(loc_data[,interact_name]) %*% data.matrix(sample_lst$air_filter_smokePM[interact_name,])
rownames(filter_rrsi) = loc_data$Metro
filter_rrsi_mean = exp(apply(
  filter_rrsi,
  1,
  mean
))

filter_rrsi_df = data.frame(DMA = names(filter_rrsi_mean), 
                            RRSI = filter_rrsi_mean, row.names = NULL) 

filter_rrsi_df$DMA = gsub(",", "", filter_rrsi_df$DMA)
filter_rrsi_df$DMA = gsub("-", " ", filter_rrsi_df$DMA)

filter_rrsi_map = dma_map %>%
  inner_join(filter_rrsi_df, by = "DMA") 

filter_rrsi_map$Search = "Air filter"

rrsi_map = rbind(pollution_rrsi_map, purifier_rrsi_map)
rrsi_map_appendix = rbind(quality_rrsi_map, filter_rrsi_map)

p7 = ggplot(data = rrsi_map) +
  geom_sf(aes(fill = RRSI)) + 
  scale_fill_viridis_c() +
  facet_wrap(~ Search, ncol = 2) +
  theme_minimal() +
  theme_void() +
  theme(legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold",size=10),
        strip.text = element_text(face="bold",size=10))

appendix_p2 = ggplot(data = rrsi_map_appendix) +
  geom_sf(aes(fill = RRSI)) + 
  scale_fill_viridis_c() +
  facet_wrap(~ Search, ncol = 2) +
  theme_minimal() +
  theme_void() +
  theme(legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold",size=10),
        strip.text = element_text(face="bold",size=10))

