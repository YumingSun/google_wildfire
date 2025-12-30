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

data_loc = "../Data/"
res_loc = "../Results/"

spatio_tempo_data = read_csv(paste0(data_loc,"master_data.csv"))
spatio_tempo_data = spatio_tempo_data %>%
  filter(Metro != "Yuma, AZ-El Centro, CA") %>%
  filter(Metro != "Reno, NV") %>%
  filter(Metro != "Medford-Klamath Falls, OR")

gtrend_overall = read_csv(paste0(data_loc,"gtrend_overall.csv"))

spatio_tempo_data_des = spatio_tempo_data %>%
  select(Metro, fireScale, AvgTemp, AvgHumidity, AvgWindSpeed,
         AvgPressure, RPL_THEME1_mean, RPL_THEME2_mean, 
         RPL_THEME3_mean, RPL_THEME4_mean,
         smokePM)

# create table 1 and table 2
des_tab = spatio_tempo_data_des %>%
  group_by(Metro) %>%
  summarise(
    fireScale_mean = mean(fireScale * 10^3, na.rm = TRUE),
    fireScale_sd   = sd(fireScale* 10^3, na.rm = TRUE),
    
    AvgTemp_mean = mean(AvgTemp, na.rm = TRUE),
    AvgTemp_sd   = sd(AvgTemp, na.rm = TRUE),
    
    AvgHumidity_mean = mean(AvgHumidity, na.rm = TRUE),
    AvgHumidity_sd   = sd(AvgHumidity, na.rm = TRUE),
    
    AvgWindSpeed_mean = mean(AvgWindSpeed, na.rm = TRUE),
    AvgWindSpeed_sd   = sd(AvgWindSpeed, na.rm = TRUE),
    
    AvgPressure_mean = mean(AvgPressure, na.rm = TRUE),
    AvgPressure_sd   = sd(AvgPressure, na.rm = TRUE),
    
    RPL_THEME1_mean = mean(RPL_THEME1_mean, na.rm = TRUE),
    
    RPL_THEME2_mean = mean(RPL_THEME2_mean, na.rm = TRUE),
    
    RPL_THEME3_mean = mean(RPL_THEME3_mean, na.rm = TRUE),
    
    RPL_THEME4_mean = mean(RPL_THEME4_mean, na.rm = TRUE),
    
    smokePM_mean = mean(smokePM, na.rm = TRUE),
    smokePM_sd   = sd(smokePM, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    fireScale = sprintf("%.2f (%.2f)", fireScale_mean, fireScale_sd),
    Temp = sprintf("%.2f (%.2f)", AvgTemp_mean, AvgTemp_sd),
    Humidity = sprintf("%.2f (%.2f)", AvgHumidity_mean, AvgHumidity_sd),
    Wind = sprintf("%.2f (%.2f)", AvgWindSpeed_mean, AvgWindSpeed_sd),
    Pressure = sprintf("%.2f (%.2f)", AvgPressure_mean, AvgPressure_sd),
  ) %>%
  select(Metro, fireScale, Temp, Humidity, Wind, Pressure,
         RPL_THEME1_mean, RPL_THEME2_mean,
         RPL_THEME3_mean, RPL_THEME4_mean)

des_tab1 = des_tab %>%
  select(Metro, fireScale, Temp, Humidity, Wind, Pressure)

des_tab2 = des_tab %>%
  select(Metro, RPL_THEME1_mean, RPL_THEME2_mean, 
         RPL_THEME3_mean, RPL_THEME4_mean)



spatio_tempo_data$DMA = gsub("-", " ", spatio_tempo_data$Metro)
spatio_tempo_data$DMA = gsub(",", "", spatio_tempo_data$DMA)

wild_fire_temporal = spatio_tempo_data %>%
  group_by(week) %>%
  summarise(acre_burned_sum = sum(incident_acres_burned),
            smokePM_mean = mean(smokePM)) %>%
  ungroup() %>%
  left_join(gtrend_overall, by = join_by(week == Week)) %>%
  select(week, acre_burned_sum, smokePM_mean, `air pollution`, `air purifier`)


google_search_temporal_long = wild_fire_temporal %>%
  pivot_longer(cols = c("air pollution","air purifier" ),
               names_to = "Search Terms Original",
               values_to = "RSI(%)") %>%
  mutate(`Search Terms` = factor(`Search Terms Original`, 
                                 levels = c("air pollution", "air purifier"),
                                 labels = c("air pollution","air purifier")))

smokePM_temporal_long = wild_fire_temporal %>%
  pivot_longer(cols = c("smokePM_mean"),
               names_to = "Average Smoke PM2.5",
               values_to = "RSI(%)") 

google_search_smokePM_temporal_long = wild_fire_temporal %>%
  pivot_longer(cols = c("air purifier","air pollution","smokePM_mean"),
               names_to = "Search Terms and PM2.5",
               values_to = "Count") %>%
  mutate(
    count_scale = case_when(
      `Search Terms and PM2.5` == "air purifier" ~ Count,
      `Search Terms and PM2.5` == "air pollution" ~ Count,
      `Search Terms and PM2.5` == "smokePM_mean" ~ Count *2
    ),
    "Type" = factor(`Search Terms and PM2.5`, 
                    levels = c("air pollution", "air purifier", "smokePM_mean"),
                    labels = c("air pollution","air purifier", "Smoke PM2.5"))
  )

search_smokePM_temporal = ggplot(google_search_smokePM_temporal_long) + 
  geom_line(aes(x = week, y = count_scale,color = Type),linewidth=0.4)+
  scale_y_continuous(
    name = "Scaled RSI",
    sec.axis = sec_axis(~.*0.5,name=TeX("\\textbf{Smoke PM2.5} ($\\mu g/m^3$)"))
  ) + 
  scale_x_date(date_labels = "%b-%Y", date_breaks  ="4 month") +
  xlab("Time")+
  theme_classic()+
  theme(plot.title = element_text(hjust = unit(7, 'mm')),
        axis.text.x = element_text(size = unit(7, 'mm'),face = "bold",
                                   angle = 90, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = unit(7, 'mm'),face = "bold"),
        axis.title = element_text(size = unit(7, 'mm'),face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = unit(7, 'mm'),face = "bold"),
        legend.position.inside=c(0.2,.85),
        plot.margin = margin(0, 2, 0, 2, "pt"))

fire_scale = ggplot(wild_fire_temporal) + 
  geom_bar(aes(x = week, y= acre_burned_sum/(1E5)),stat="identity") +
  ylab(TeX("\\textbf{Acres Burned} ($\\times 10^5$)"))+
  theme_classic()+
  theme(axis.text.y = element_text(size = unit(7, 'mm'),face = "bold"),
        axis.title.y = element_text(size = unit(7, 'mm'),face = "bold"),
        axis.text.x  = element_blank(),
        axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank(),
        plot.margin = margin(2, 0, 0, 2, "pt"))


figure_1 = ggarrange(fire_scale, search_smokePM_temporal, 
               heights = c(0.7, 2),
               ncol = 1, nrow = 2,align = "v")


wild_fire_loc = spatio_tempo_data %>%
  group_by(Metro) %>%
  summarise(acre_burned_sum = sum(incident_acres_burned),
            smokePM_mean = mean(smokePM),
            air_purifier_mean = mean(`air purifier`),
            air_pollution_mean = mean(`air pollution`),
  ) 
wild_fire_loc$DMA = gsub("-", " ", wild_fire_loc$Metro)
wild_fire_loc$DMA = gsub(",", "", wild_fire_loc$DMA)

dma_map = st_read(paste0(data_loc,"NatDMA/NatDMA.shp"))
dma_map$NAME = gsub("-", " ", dma_map$NAME)
colnames(dma_map)[4] = "DMA"

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
      Category == "air_purifier_mean" ~ "air purifier",
      Category == "air_pollution_mean" ~ "air pollution",
    )
  )

plot_smoke = dma_map %>%
  inner_join(wild_fire_loc, by = "DMA")

google_search_loc_pollution = plot_google_search %>%
  filter(search_terms == "air pollution") %>%
  ggplot() +
  geom_sf(aes(fill = `RSI (%)`)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme_void() + labs(title="air pollution",fill = "Scaled RSI") +
  theme(legend.title = element_text(face = "bold",size=unit(3, 'mm')),
        legend.text = element_text(size=unit(3, 'mm')),
        plot.title = element_text(face = "bold",size=unit(3, 'mm'),
                                  hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text = element_text(face="bold",size=unit(3, 'mm')),
        legend.key.size = unit(0.1, 'cm'),
        legend.position.inside=c(0.9,0.7),
        panel.spacing = unit(0, "cm"))

google_search_loc_purifier = plot_google_search %>%
  filter(search_terms == "air purifier") %>%
  ggplot() +
  geom_sf(aes(fill = `RSI (%)`)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme_void() + labs(title="air purifier",fill = "Scaled RSI") +
  theme(legend.title = element_text(face = "bold",size=unit(3, 'mm')),
        legend.text = element_text(size=unit(3, 'mm')),
        plot.title = element_text(face = "bold",size=unit(3, 'mm'),
                                  hjust=0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text = element_text(face="bold",size=unit(3, 'mm')),
        legend.key.size = unit(0.1, 'cm'),
        legend.position.inside=c(0.9,0.7),
        panel.spacing = unit(0, "cm"))

smoke_loc = ggplot(data = plot_smoke) +
  geom_sf(aes(fill = smokePM_mean)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme_void() + labs(fill=TeX('\\textbf{Smoke PM}$_{2.5}$',bold=TRUE),
                      title=TeX('\\textbf{Smoke PM}$_{2.5}$',bold=TRUE)) + 
  theme(legend.title = element_text(face = "bold",size=unit(3, 'mm')),
        legend.text = element_text(size=unit(3, 'mm')),
        plot.title = element_text(face = "bold",size=unit(3, 'mm'),hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text = element_text(face="bold",size=unit(3, 'mm')),
        legend.key.size = unit(0.1, 'cm'),
        legend.position.inside=c(0.9,0.7),
        panel.spacing = unit(0, "cm"))

figure_2 = ggarrange(google_search_loc_pollution,
               google_search_loc_purifier,
               smoke_loc, widths = c(1, 1, 1),
               ncol = 3, nrow = 1)


load(paste0(res_loc,"post_samp_lst.RData"))
load(paste0(res_loc,"fixed_effect_lst.RData"))

# Sensitivity for pre-pandemic period
load(paste0(res_loc,"post_samp_lst_sens.RData"))
load(paste0(res_loc,"fixed_effect_lst_sens.RData"))


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

figure_3 = main_effect_df %>%
  filter(outcome %in% c("air_pollution", "air_purifier")) %>%
  mutate(outcome = factor(outcome, labels = c("air pollution", "air purifier"))) %>%
  ggplot(aes(x = lag_num, y = mean)) +
  geom_pointrange(aes(ymin=`0.025quant`, ymax = `0.975quant`), size = 0.05) +
  geom_line() + 
  xlab("Lag (week)") +
  ylab("RSIR") + 
  scale_y_continuous(breaks = c(1, 2, 3, 4))+
  facet_wrap(~ outcome, ncol = 2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=unit(4, 'mm')),
        axis.text.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.text.x = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.x = element_text(size = unit(4, 'mm'),face = "bold"))

supplement_figure_b_1 = main_effect_df %>%
  filter(outcome %in% c("air_quality", "air_filter")) %>%
  mutate(outcome = factor(outcome, labels = c("air quality", "air filter"))) %>%
  ggplot(aes(x = lag_num, y = mean)) +
  geom_pointrange(aes(ymin=`0.025quant`, ymax = `0.975quant`), size = 0.05) +
  geom_line() +
  xlab("Lag (week)") +
  ylab("RSIR") + 
  facet_wrap(~ outcome, ncol = 2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=unit(4, 'mm')),
        axis.text.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.text.x = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.x = element_text(size = unit(4, 'mm'),face = "bold"))

# sensitivity analysis for pre-pandamic period
supplement_figure_b_8 = main_effect_df %>%
  filter(outcome %in% c("air_filter", "air_purifier")) %>%
  mutate(outcome = factor(outcome, 
                          levels = c("air_filter", "air_purifier"),
                          labels = c("Air filter", "Air purifier"))) %>%
  ggplot(aes(x = lag_num, y = mean)) +
  geom_pointrange(aes(ymin=`0.025quant`, ymax = `0.975quant`), size = 0.05) +
  geom_line() + 
  xlab("Lag (week)") +
  ylab("RSIR") + 
  facet_wrap(~ outcome, ncol = 2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=unit(4, 'mm')),
        axis.text.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.text.x = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.x = element_text(size = unit(4, 'mm'),face = "bold"))

# alternative pm2.5 estimates
spatio_tempo_data_newpm = read_csv(paste0(data_loc,"master_data_newpm.csv"))
spatio_tempo_data_newpm = spatio_tempo_data_newpm %>%
  filter(Metro != "Yuma, AZ-El Centro, CA") %>%
  filter(Metro != "Reno, NV") %>%
  filter(Metro != "Medford-Klamath Falls, OR")

wild_fire_temporal_newpm = spatio_tempo_data_newpm %>%
  group_by(week) %>%
  summarise(acre_burned_sum = sum(incident_acres_burned),
            smokePM_mean = mean(smokePM)) %>%
  ungroup() %>%
  left_join(gtrend_overall, by = join_by(week == Week)) %>%
  select(week, acre_burned_sum, smokePM_mean, `air quality`, `air purifier`)

google_search_smokePM_temporal_long_pm = wild_fire_temporal_newpm %>%
  pivot_longer(cols = c("air purifier","air quality","smokePM_mean"),
               names_to = "Search Terms and PM2.5",
               values_to = "Count") %>%
  mutate(
    count_scale = case_when(
      `Search Terms and PM2.5` == "air purifier" ~ Count,
      `Search Terms and PM2.5` == "air quality" ~ Count,
      `Search Terms and PM2.5` == "smokePM_mean" ~ Count *2
    ),
    "Type" = factor(`Search Terms and PM2.5`, 
                    levels = c("air quality", "air purifier", "smokePM_mean"),
                    labels = c("Air quality","Air purifier", "Smoke PM2.5"))
  )

compare_pm = rbind(
  google_search_smokePM_temporal_long %>%
    filter(Type == "Smoke PM2.5") %>%
    select(week, Count) %>%
    mutate(data = "Burke et al."),
  google_search_smokePM_temporal_long_pm %>%
    filter(Type == "Smoke PM2.5") %>%
    select(week, Count) %>%
    mutate(data = "Aguilera et al.")
)

supplement_figure_b_3 = ggplot(compare_pm) + 
  geom_line(aes(x = week, y = Count,color = data),linewidth=0.4)+
  scale_y_continuous(
    name = TeX("\\textbf{Smoke PM2.5} ($\\mu g/m^3$)")
  ) + 
  scale_x_date(date_labels = "%b-%Y", date_breaks  ="4 month") +
  xlab("Time")+
  theme_classic()+
  theme(plot.title = element_text(hjust = unit(7, 'mm')),
        axis.text.x = element_text(size = unit(7, 'mm'),face = "bold",
                                   angle = 90, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = unit(7, 'mm'),face = "bold"),
        axis.title = element_text(size = unit(7, 'mm'),face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = unit(7, 'mm'),face = "bold"),
        legend.position.inside=c(0.2,.85),
        plot.margin = margin(0, 2, 0, 2, "pt"))

wild_fire_loc_newpm = spatio_tempo_data_newpm %>%
  group_by(Metro) %>%
  summarise(acre_burned_sum = sum(incident_acres_burned),
            smokePM_mean = mean(smokePM),
            air_purifier_mean = mean(`air purifier`),
            air_pollution_mean = mean(`air pollution`)
  ) 
wild_fire_loc_newpm$DMA = gsub("-", " ", wild_fire_loc_newpm$Metro)
wild_fire_loc_newpm$DMA = gsub(",", "", wild_fire_loc_newpm$DMA)

plot_smoke_newpm = dma_map %>%
  inner_join(wild_fire_loc_newpm, by = "DMA")

plot_smoke_compare = rbind(
  plot_smoke_newpm %>% mutate(data = "Aguilera et al."),
  plot_smoke %>% mutate(data = "Burke et al.")
)

supplement_figure_b_4 = ggplot(data = plot_smoke_compare) +
  geom_sf(aes(fill = smokePM_mean)) +
  scale_fill_viridis_c() +
  facet_wrap(~ data, ncol = 2)+
  theme_minimal() +
  theme_void() + labs(fill=TeX('\\textbf{Smoke PM}$_{2.5}$',bold=TRUE)) + 
  theme(legend.title = element_text(face = "bold",size=unit(3, 'mm')),
        legend.text = element_text(size=unit(3, 'mm')),
        plot.title = element_text(face = "bold",size=unit(3, 'mm'),hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text = element_text(face="bold",size=unit(3, 'mm')),
        legend.key.size = unit(0.1, 'cm'),
        legend.position.inside=c(0.9,0.7),
        panel.spacing = unit(0, "cm"))


# load for alternative PM2.5 estimates
# alternative PM2.5 estimates main effect
load(paste0(res_loc,"post_samp_lst_newpm.RData"))
load(paste0(res_loc,"fixed_effect_lst_newpm.RData"))

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

supplement_figure_b_5 = main_effect_df %>%
  filter(outcome %in% c("air_pollution", "air_purifier")) %>%
  mutate(outcome = factor(outcome, labels = c("air pollution", "air purifier"))) %>%
  ggplot(aes(x = lag_num, y = mean)) +
  geom_pointrange(aes(ymin=`0.025quant`, ymax = `0.975quant`), size = 0.05) +
  geom_line() + 
  xlab("Lag (week)") +
  ylab("RSIR") + 
  scale_y_continuous(breaks = c(1, 2, 3, 4))+
  facet_wrap(~ outcome, ncol = 2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=unit(4, 'mm')),
        axis.text.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.y = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.text.x = element_text(size = unit(4, 'mm'),face = "bold"),
        axis.title.x = element_text(size = unit(4, 'mm'),face = "bold"))
#end alternative PM2.5 main effect

# effect modifiers
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
load(paste0(res_loc,"post_samp_lst.RData"))
load(paste0(res_loc,"fixed_effect_lst.RData"))
sample_lst = post_samp_lst
spatio_tempo_data = read_csv(paste0(data_loc, 
                                    "master_data.csv"))
spatio_tempo_data = spatio_tempo_data %>%
  filter(Metro != "Reno, NV") %>%
  filter(Metro != "Medford-Klamath Falls, OR") %>%
  filter(Metro != "Yuma, AZ-El Centro, CA")

# Load for pre-pandemic sensitivity analysis
load(paste0(res_loc,"post_samp_lst_sens.RData"))
load(paste0(res_loc,"fixed_effect_lst_sens.RData"))
spatio_tempo_data = read_csv(paste0(data_loc, 
                                    "master_data_sens.csv"))
spatio_tempo_data = spatio_tempo_data %>%
  filter(Metro != "Reno, NV") %>%
  filter(Metro != "Medford-Klamath Falls, OR") %>%
  filter(Metro != "Yuma, AZ-El Centro, CA")

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

# run for supplement table b1 and b2
mean_ci_from_post = function(data_vec){
  avg = mean(data_vec)
  ci = quantile(data_vec,p = c(0.025, 0.975))
  res_num = c(avg, ci)
  res_out = sprintf("%.3f (%.3f, %.3f)", 
                    res_num[1],res_num[2],res_num[3])
  return(res_out)
}


modifiers = c("fireScale", "Temp", "Humidity", "WindSpeed", "Pressure",
              "THEME1", "THEME2", "THEME3", "THEME4")
modifier_quality_lst = list()

for (m in modifiers) {
  loc_id = str_detect(rownames(sample_lst$air_quality_smokePM), m)
  res_ = apply(sample_lst$air_quality_smokePM[loc_id,],1,mean_ci_from_post)
  names(res_) = c("Lag0", "Lag1", "Lag2", "Lag3")
  modifier_quality_lst[[m]] = res_
  print(m)
}

supplement_tab_b_1 = as.data.frame(do.call(rbind, modifier_quality_lst))

modifier_filter_lst = list()

for (m in modifiers) {
  loc_id = str_detect(rownames(sample_lst$air_filter_smokePM), m)
  res_ = apply(sample_lst$air_filter_smokePM[loc_id,],1,mean_ci_from_post)
  names(res_) = c("Lag0", "Lag1", "Lag2", "Lag3")
  modifier_filter_lst[[m]] = res_
  print(m)
}

supplement_tab_b_1 = as.data.frame(do.call(rbind, modifier_filter_lst))

for (m in modifiers) {
  loc_id = str_detect(rownames(sample_lst$air_filter_smokePM), m)
  res_ = apply(sample_lst$air_filter_smokePM[loc_id,],1,mean_ci_from_post)
  names(res_) = c("Lag0", "Lag1", "Lag2", "Lag3")
  modifier_filter_lst[[m]] = res_
  print(m)
}

supplement_tab_b_1 = as.data.frame(do.call(rbind, modifier_filter_lst))
#end for supplement b1 and b2

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


create_data_for_ci = function(var, population_data, lb = 0.025, ub = 0.975){
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
temp_df = create_data_for_ci("AvgTemp", spatio_tempo_data)
humid_df = create_data_for_ci("AvgHumidity", spatio_tempo_data)
wind_df = create_data_for_ci("AvgWindSpeed", spatio_tempo_data)
pressure_df = create_data_for_ci("AvgPressure", spatio_tempo_data)
fire_scale_df = create_data_for_ci("fireScale", spatio_tempo_data,
                                   lb = 0.8, ub = 0.99)
socio_df = create_data_for_ci("RPL_THEME1_mean", spatio_tempo_data,
                              lb = 0.2, ub = 0.8)
house_df = create_data_for_ci("RPL_THEME2_mean", spatio_tempo_data,
                              lb = 0.2, ub = 0.8)
race_df = create_data_for_ci("RPL_THEME3_mean", spatio_tempo_data,
                             lb = 0.2, ub = 0.8)
trans_df = create_data_for_ci("RPL_THEME4_mean", spatio_tempo_data,
                              lb = 0.2, ub = 0.8)

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

pollution_temp_rrsi = calculate_rrsi(temp_df,pollution_temp_coef,"Temperature (°F)")
pollution_humid_rrsi = calculate_rrsi(humid_df,pollution_humid_coef,"Humidity (%)")
pollution_wind_rrsi = calculate_rrsi(wind_df,pollution_wind_coef,"Wind Speed (mph)")
pollution_pressure_rrsi = calculate_rrsi(pressure_df,pollution_pressure_coef,
                                         "Atmospheric Pressure (inHg)")
pollution_fire_scale_rrsi = calculate_rrsi(fire_scale_df,pollution_fire_scale_coef,
                                           "Fire Scale (per mille)")
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
                                          "Fire Scale (per mille)")
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
                           levels = c("Fire Scale (per mille)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))%>%
  mutate(var_value = case_when(
    var_name == "Fire Scale (per mille)" ~ var_value * 1000,
    TRUE ~ var_value
  ))


purifier_rrsi = rbind.data.frame(purifier_temp_rrsi, purifier_humid_rrsi,
                                 purifier_wind_rrsi, purifier_pressure_rrsi,
                                 purifier_fire_scale_rrsi,
                                 purifier_socio_rrsi, purifier_house_rrsi,
                                 purifier_race_rrsi, purifier_trans_rrsi) %>%
  mutate(var_name = factor(var_name, 
                           levels = c("Fire Scale (per mille)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))%>%
  mutate(var_value = case_when(
    var_name == "Fire Scale (per mille)" ~ var_value * 1000,
    TRUE ~ var_value
  ))

quality_temp_rrsi = calculate_rrsi(temp_df,quality_temp_coef,"Temperature (°F)")
quality_humid_rrsi = calculate_rrsi(humid_df,quality_humid_coef,"Humidity (%)")
quality_wind_rrsi = calculate_rrsi(wind_df,quality_wind_coef,"Wind Speed (mph)")
quality_pressure_rrsi = calculate_rrsi(pressure_df,quality_pressure_coef,
                                       "Atmospheric Pressure (inHg)")
quality_fire_scale_rrsi = calculate_rrsi(fire_scale_df,quality_fire_scale_coef,
                                         "Fire Scale (per mille)")
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
                                        "Fire Scale (per mille)")
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
                           levels = c("Fire Scale (per mille)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))%>%
  mutate(var_value = case_when(
    var_name == "Fire Scale (per mille)" ~ var_value * 1000,
    TRUE ~ var_value
  ))

filter_rrsi = rbind.data.frame(filter_temp_rrsi, filter_humid_rrsi,
                               filter_wind_rrsi, filter_pressure_rrsi,
                               filter_fire_scale_rrsi,
                               filter_socio_rrsi, filter_house_rrsi,
                               filter_race_rrsi, filter_trans_rrsi) %>%
  mutate(var_name = factor(var_name, 
                           levels = c("Fire Scale (per mille)", "Temperature (°F)", "Humidity (%)", 
                                      "Wind Speed (mph)", 
                                      "Atmospheric Pressure (inHg)",
                                      "Socioeconomic Status",
                                      "Household Characteristics",
                                      "Racial & Ethnic Minority Status",
                                      "Housing Type & Transportation"
                           )))%>%
  mutate(var_value = case_when(
    var_name == "Fire Scale (per mille)" ~ var_value * 1000,
    TRUE ~ var_value
  ))

figure_4 = ggplot(data = pollution_rrsi, aes(x = var_value)) + 
  geom_line(aes(y = `50%CI`)) +
  geom_ribbon(aes(ymin=`2.5%CI`, ymax=`97.5%CI`), linetype=2, alpha=0.1) +
  facet_wrap(~ var_name, ncol = 3, scales = "free") + 
  xlab("") +
  ylab("cumulative RSIR") + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=unit(6, 'mm')),
        axis.text.y = element_text(size = unit(6, 'mm'),face = "bold"),
        axis.title.y = element_text(size = unit(6, 'mm'),face = "bold"),
        axis.text.x = element_text(size = unit(6, 'mm'),face = "bold"),
        axis.title.x = element_text(size = unit(6, 'mm'),face = "bold")) 


figure_5 = ggplot(data = purifier_rrsi, aes(x = var_value)) + 
  geom_line(aes(y = `50%CI`)) +
  geom_ribbon(aes(ymin=`2.5%CI`, ymax=`97.5%CI`), linetype=2, alpha=0.1) +
  facet_wrap(~ var_name, ncol = 3, scales = "free") + 
  xlab("") +
  ylab("cumulative RSIR") + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold",size=unit(6, 'mm')),
        axis.text.y = element_text(size = unit(6, 'mm'),face = "bold"),
        axis.title.y = element_text(size = unit(6, 'mm'),face = "bold"),
        axis.text.x = element_text(size = unit(6, 'mm'),face = "bold"),
        axis.title.x = element_text(size = unit(6, 'mm'),face = "bold"))



loc_data = spatio_tempo_data %>%
  filter(Metro != "Reno, NV") %>%
  filter(Metro != "Medford-Klamath Falls, OR") %>%
  filter(Metro != "Yuma, AZ-El Centro, CA") %>%
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

pollution_rrsi_ci = exp(apply(
  pollution_rrsi,
  1,
  quantile,
  c(0.025,0.975)
))

pollution_rrsi_ci_txt = sprintf("(%.2f, %.2f)", pollution_rrsi_ci["2.5%",], 
                                pollution_rrsi_ci["97.5%",])
pollution_rrsi_mean_txt = sprintf("%.2f", pollution_rrsi_mean)

dma_names = names(pollution_rrsi_mean)
pollution_rrsi_tab_txt = cbind(dma_names,pollution_rrsi_mean_txt,pollution_rrsi_ci_txt)


pollution_rrsi_df = data.frame(DMA = names(pollution_rrsi_mean), 
                               RRSI = pollution_rrsi_mean, row.names = NULL) 

pollution_rrsi_df$DMA = gsub(",", "", pollution_rrsi_df$DMA)
pollution_rrsi_df$DMA = gsub("-", " ", pollution_rrsi_df$DMA)
pollution_rrsi_map = dma_map %>%
  inner_join(pollution_rrsi_df, by = "DMA") 
pollution_rrsi_map$Search = "air pollution"

quality_rrsi = data.matrix(loc_data[,interact_name]) %*% data.matrix(sample_lst$air_quality_smokePM[interact_name,])
rownames(quality_rrsi) = loc_data$Metro
quality_rrsi_mean = exp(apply(
  quality_rrsi,
  1,
  mean
))

quality_rrsi_ci = exp(apply(
  quality_rrsi,
  1,
  quantile,
  c(0.025,0.975)
))

quality_rrsi_ci_txt = sprintf("(%.2f, %.2f)", quality_rrsi_ci["2.5%",], 
                              quality_rrsi_ci["97.5%",])
quality_rrsi_mean_txt = sprintf("%.2f", quality_rrsi_mean)

dma_names = names(quality_rrsi_mean)
quality_rrsi_tab_txt = cbind(dma_names,quality_rrsi_mean_txt,quality_rrsi_ci_txt)


quality_rrsi_df = data.frame(DMA = names(quality_rrsi_mean), 
                             RRSI = quality_rrsi_mean, row.names = NULL) 

quality_rrsi_df$DMA = gsub(",", "", quality_rrsi_df$DMA)
quality_rrsi_df$DMA = gsub("-", " ", quality_rrsi_df$DMA)
quality_rrsi_map = dma_map %>%
  inner_join(quality_rrsi_df, by = "DMA") 
quality_rrsi_map$Search = "air quality"

purifier_rrsi = data.matrix(loc_data[,interact_name]) %*% data.matrix(sample_lst$air_purifier_smokePM[interact_name,])
rownames(purifier_rrsi) = loc_data$Metro
purifier_rrsi_mean = exp(apply(
  purifier_rrsi,
  1,
  mean
))

purifier_rrsi_ci = exp(apply(
  purifier_rrsi,
  1,
  quantile,
  c(0.025,0.975)
))

purifier_rrsi_ci_txt = sprintf("(%.2f, %.2f)", purifier_rrsi_ci["2.5%",], 
                               purifier_rrsi_ci["97.5%",])
purifier_rrsi_mean_txt = sprintf("%.2f", purifier_rrsi_mean)

dma_names = names(purifier_rrsi_mean)
purifier_rrsi_tab_txt = cbind(dma_names,purifier_rrsi_mean_txt,purifier_rrsi_ci_txt)


purifier_rrsi_df = data.frame(DMA = names(purifier_rrsi_mean), 
                              RRSI = purifier_rrsi_mean, row.names = NULL) 

purifier_rrsi_df$DMA = gsub(",", "", purifier_rrsi_df$DMA)
purifier_rrsi_df$DMA = gsub("-", " ", purifier_rrsi_df$DMA)

purifier_rrsi_map = dma_map %>%
  inner_join(purifier_rrsi_df, by = "DMA") 

purifier_rrsi_map$Search = "air purifier"

filter_rrsi = data.matrix(loc_data[,interact_name]) %*% data.matrix(sample_lst$air_filter_smokePM[interact_name,])
rownames(filter_rrsi) = loc_data$Metro
filter_rrsi_mean = exp(apply(
  filter_rrsi,
  1,
  mean
))

filter_rrsi_ci = exp(apply(
  filter_rrsi,
  1,
  quantile,
  c(0.025,0.975)
))

filter_rrsi_ci_txt = sprintf("(%.2f, %.2f)", filter_rrsi_ci["2.5%",], 
                             filter_rrsi_ci["97.5%",])
filter_rrsi_mean_txt = sprintf("%.2f", filter_rrsi_mean)

dma_names = names(filter_rrsi_mean)
filter_rrsi_tab_txt = cbind(dma_names,filter_rrsi_mean_txt,filter_rrsi_ci_txt)


filter_rrsi_df = data.frame(DMA = names(filter_rrsi_mean), 
                            RRSI = filter_rrsi_mean, row.names = NULL) 

filter_rrsi_df$DMA = gsub(",", "", filter_rrsi_df$DMA)
filter_rrsi_df$DMA = gsub("-", " ", filter_rrsi_df$DMA)

filter_rrsi_map = dma_map %>%
  inner_join(filter_rrsi_df, by = "DMA") 

filter_rrsi_map$Search = "air filter"

rrsi_map = rbind(pollution_rrsi_map, purifier_rrsi_map)
rrsi_map_appendix = rbind(quality_rrsi_map, filter_rrsi_map)

pollution_rrsi_tab_df = data.frame(pollution_rrsi_tab_txt)
colnames(pollution_rrsi_tab_df) = c("DMA", "RRSI", "95% CI")

purifier_rrsi_tab_df = data.frame(purifier_rrsi_tab_txt)
colnames(purifier_rrsi_tab_df) = c("DMA", "RRSI", "95% CI")

quality_rrsi_tab_df = data.frame(quality_rrsi_tab_txt)
colnames(quality_rrsi_tab_df) = c("DMA", "RRSI", "95% CI")

filter_rrsi_tab_df = data.frame(filter_rrsi_tab_txt)
colnames(filter_rrsi_tab_df) = c("DMA", "RRSI", "95% CI")

table_3 = inner_join(pollution_rrsi_tab_df, purifier_rrsi_tab_df, by = "DMA",
                    suffix = c(".pollution", ".purifier"))

colnames(rrsi_map)[14] = "RSIR"
figure_6 = rrsi_map %>%
  ggplot() +
  geom_sf(aes(fill = RSIR)) + 
  scale_fill_viridis_c(
    trans = "log",                    
    breaks = c(0.5, 1, 2, 5, 10, 20, 30),
    labels = scales::label_number(accuracy = 0.1),
    name = "RSIR"
  ) +
  facet_wrap(~ Search, ncol = 2) + 
  theme_minimal() +
  theme_void() +
  theme(legend.title = element_text(face = "bold",size=unit(4, 'mm')),
        plot.title = element_text(face = "bold",size=unit(5, 'mm')),
        strip.text = element_text(face="bold",size=unit(5, 'mm')),
        legend.text = element_text(size=unit(4, 'mm')),
        legend.key.size = unit(0.3, 'cm'),
        legend.position=c(0.95,0.6))

colnames(rrsi_map_appendix)[14] = "RSIR"
appendix_p4 = rrsi_map_appendix %>%
  ggplot() + 
  geom_sf(aes(fill = RSIR)) +
  scale_fill_viridis_c(
    trans = "log",                  
    breaks = c(0.5, 1, 2, 5, 10, 20, 30),
    labels = scales::label_number(accuracy = 0.1),
    name = "RSIR"
  ) +
  facet_wrap(~ Search, ncol = 2) + 
  theme_minimal() +
  theme_void() +
  theme(legend.title = element_text(face = "bold",size=unit(4, 'mm')),
        plot.title = element_text(face = "bold",size=unit(5, 'mm')),
        strip.text = element_text(face="bold",size=unit(5, 'mm')),
        legend.text = element_text(size=unit(4, 'mm')),
        legend.key.size = unit(0.3, 'cm'),
        legend.position=c(0.95,0.6))

