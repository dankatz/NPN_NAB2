### General information ####################################################################
# Project goals: 
# -assess correlations between NPN and NAB data 
# -provide proof of concept for future continental-scale airborne pollen models
#
# current/potential project participants:
# Dan Katz, Theresa Crimmins, Liz Vogt, Arie Managan, Claudia Brown, Ellen Denny, 
# Shubhayu Saha, Ambarish (Rish) Vaidyanathan
#
# This script is for data exploration and analyses
# Note: One of the goals is to keep this as simple as possible;
# full phenological models are a future project


### set up working environment #############################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(imputeTS)
library(zoo)
library(ggpmisc)
library(viridis)
library(ggthemes)



### load in and prepare nab data ###############################################################
nab_raw <- read_csv("C:/Users/dsk856/Box/texas/NAB/NAB2009_2019_pollen_200508.csv", guess_max = 92013) %>% 
  filter(site == "Armonk" | site == "Marietta (Atlanta)")

#expand to include missing dates
date_station_grid <- expand_grid(seq(min(nab_raw$date), max(nab_raw$date), by = '1 day'), unique(nab_raw$site)) %>% 
  `colnames<-`(c("Date", "site")) %>%
  filter(!is.na(site)) %>% 
  ungroup()

nab <- left_join(date_station_grid, nab_raw) %>% 
  dplyr::select(Date, site, 
                Acer, Alnus, Ambrosia, Betula, Cupressaceae, Fraxinus, 
                Pinaceae, Populus, Quercus, Ulmus) %>% 
  pivot_longer(cols = c(Acer, Alnus, Ambrosia, Betula, Cupressaceae, Fraxinus, Pinaceae, Populus, Quercus, Ulmus),
               names_to = "taxon", values_to = "pol") %>% 
  mutate(s_year = year(Date),
         doy = yday(Date),
         date_noyr = format(Date, format="%m-%d")) %>% 
  arrange(taxon, site, Date) %>%
  group_by(taxon, site) %>%
  mutate(pol_m =  round(na_interpolation(pol), 1), #linear interpolation of missing data #maxgap
         pol_m_ma = round(rollmean(pol_m, 7, na.pad=TRUE),1),
         pol_na = case_when(is.na(pol) ~ 1)) %>%  #1 week moving average
  ungroup()

#heatmaps of numner of observations taxon x year x NAB station
nab %>% 
  group_by(taxon, site, s_year) %>%
  filter(pol > 25) %>% #filter observations where pollen concentration was > x grains/m3
  filter(!is.na(pol)) %>% 
  summarize(obs_per_year = n()) %>% 
  mutate(obs_per_year_d = #cut(obs_per_year, breaks = c(0, 4, 10), labels = c("<5 days", "5-10 days", ">10 days")))
           case_when(obs_per_year < 5.1 ~ "<5",
                     obs_per_year > 5.1 &  obs_per_year < 10.1 ~ "5 - 10",
                     obs_per_year > 10.2 ~ ">10")) %>% 
  #visualize number of obs per year per station
  ggplot(aes(x = taxon, y = as.factor(s_year), fill = obs_per_year_d)) + geom_tile() + facet_wrap(~site) +
  scale_fill_viridis_d(name = "n/yr \n(days where pollen\n >25 grains/m3)") +
  theme_bw() + ylab("year")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#visualizing the NAB data as a time series
nab %>% 
  #filter(site == "Marietta (Atlanta)") %>%  #filter(site == "Armonk") %>%
  #filter(s_year == 2020) %>% 
  filter(taxon == "Cupressaceae") %>% 
  ggplot(aes(x = doy, y = pol + .5)) + theme_bw() + scale_y_log10() + 
  geom_point(aes(x= doy, y = pol_m + .5), col = "red", size = 0.3) + geom_point() + 
  ylab(pollen~grains~per~m^3) + scale_x_continuous(limits = c(1, 135)) + facet_grid(site~s_year) 

#calculating the maximum observed pollen concentration per taxon per station per year
pol_yearly_max <- nab %>%
  group_by(s_year, site, taxon) %>%
  slice_max(pol) %>%
  rename(max_pol = pol,
         max_pol_date = Date) %>%
  mutate(max_pol_doy = yday(max_pol_date)) %>%
  dplyr::select(site, s_year, taxon, max_pol, max_pol_date, max_pol_doy)

nab_join <- left_join(nab, pol_yearly_max) %>% 
  mutate(site = case_when(site == "Marietta (Atlanta)" ~ "Atlanta",
                          site == "Armonk" ~ "Armonk",
                          site == "Waterbury" ~ "Waterbury")) %>% 
  rename(date = Date)




### load in and prepare NPN data ###############################################################
npn_raw <- read_csv("C:/Users/dsk856/Box/texas/pheno/npn_common_anemo_active_flow.csv", guess_max = 672676)

#creating summary metrics for npn data
hist(npn_raw$tmean_dif)
hist(npn_raw$NAB_min_dist, breaks = 50)

filt_tmean_dif <- 2 #filter NPN observations that are within X degrees celsius of the nearest NAB station
filt_dist <- 322 #filter NPN observations that are within X km of the nearest NAB station #322 km = 200 miles


npn <- npn_raw %>% 
  filter(tmean_dif > -filt_tmean_dif & tmean_dif < filt_tmean_dif) %>% #MAT filtering
  filter(NAB_min_dist < filt_dist) %>% #distance filtering
  mutate(s_year = year(observation_date),
         doy = yday(observation_date),
         date_noyr = format(observation_date, format="%m-%d"),
         taxon = genus,
         taxon = case_when(genus == "Juniperus" ~ "Cupressaceae",
                           genus == "Taxodium" ~ "Cupressaceae",
                           genus == "Thuja" ~ "Cupressaceae",
                           genus == "Chamaecyparis" ~ "Cupressaceae",
                           genus == "Metasequoia" ~ "Cupressaceae",
                           genus == "Calocedrus" ~ "Cupressaceae",
                           TRUE ~ taxon)) %>% 
  rename(site = NAB_station, 
         date = observation_date) 

#expand to include missing dates
date_station_grid_npn <- expand_grid(seq(min(npn$date), max(npn$date), by = '1 day'), unique(npn$site),
                                     unique(npn$taxon)) %>% 
  `colnames<-`(c("date", "site","taxon")) %>%
  filter(!is.na(site)) %>% 
  ungroup()

npn_summary <- npn %>% 
  group_by(taxon, date, doy, date_noyr, s_year, site) %>% 
  summarize(mean_flow = mean(phenophase_status),
            mean_prop_flow = mean(flow_prop),
            n_obs = as.numeric(n())) %>% 
  ungroup() %>% 
  left_join(date_station_grid_npn, .) %>% 
  mutate(n_obs = case_when( is.na(n_obs) ~ 0,
                            n_obs > 0 ~ n_obs))%>% 
  arrange(taxon, site, s_year, date)

npn_season_summary_nobs <- npn_summary %>% 
  group_by(taxon, s_year, site) %>% 
  summarize(nobs_per_season = sum(n_obs))
npn_season_summary_nobs_open_flow <- npn%>% 
  filter(phenophase_status == 1) %>% 
  group_by(taxon, s_year, site) %>% 
  summarize(nobs_yes_per_season = n())

npn_join <-  left_join(npn_summary, npn_season_summary_nobs) %>% 
  left_join(., npn_season_summary_nobs_open_flow) %>% 
  filter(nobs_yes_per_season > 0) %>% 
  mutate(mean_flow_m = round(na_interpolation(mean_flow),1),
         mean_flow_m_ma = round(rollmean(mean_flow_m, 7, na.pad=TRUE),2),
         mean_prop_flow_m = round(na_interpolation(mean_prop_flow), 1),
         mean_prop_flow_m_ma = round(rollmean(mean_prop_flow_m, 7, na.pad=TRUE),2),
         mean_nobs_ma = round(rollmean(n_obs, 7, na.pad=TRUE),2))

#heatmaps of numner of observations taxon x year x NAB station
npn_season_summary_nobs_open_flow %>% 
  group_by(taxon, site, s_year) %>%
  #visualize number of obs per year per station
  ggplot(aes(x = taxon, y = as.factor(s_year), fill = nobs_yes_per_season)) + geom_tile() + facet_wrap(~site) +
  scale_fill_viridis(name = "count", trans = "log10", breaks = c(1,10, 50, 100, 250)) +
  #breaks = my_breaks, labels = my_breaks)
  #scale_fill_viridis_c(name = "flowering events observed/yr (n)") +
  theme_bw() + ylab("year")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#visualizing the NPN data as a time series
npn_join %>% 
  filter(site != "Waterbury") %>%
  filter(s_year > 2013 & s_year < 2021) %>% 
  filter(taxon == "Quercus") %>% 
  ggplot(aes(x = doy, y = mean_flow_m_ma, color = n_obs)) + theme_few() + scale_y_log10() + 
  geom_point() + geom_line(color = "gray")+
  ylab("mean proportion of individuals flowering (smoothed)") + #scale_x_continuous(limits = c(50, 250)) + 
  scale_color_viridis_c(option = "magma", direction = -1, name = "observations per day (n)")+ facet_grid(site~s_year) 

taxon_npn_yearly_max <- npn_join %>% 
  group_by(taxon, site, s_year) %>% 
  slice_max(mean_prop_flow_m_ma, with_ties = FALSE) %>% 
  mutate(date_prop_flow_m_ma_max = date,
         doy_prop_flow_m_ma_max = yday(date_prop_flow_m_ma_max))



### combine NAB and NPN data #########################################################################
nab_npn <- left_join(nab_join, npn_join)

#directly compare data
formula <- y ~ x 
nab_npn %>% 
  filter(nobs_yes_per_season > 50) %>% #only include taxa x seasons where there's decent NPN data
  ggplot(aes(y = pol + 1, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab(pollen~grains~per~m^3)


nab_npn %>% 
  filter(nobs_yes_per_season > 50) %>% #only include taxa x seasons where there's decent NPN data
  filter(doy > 70 & doy < 150) %>% 
  filter(taxon == "Quercus") %>% 
  ggplot(aes(y = pol + 1)) + geom_histogram() + theme_bw() +
  facet_wrap(site~taxon) + scale_y_log10() 

#visualize time series side by side
nab_npn %>% 
  filter(nobs_yes_per_season > 50) %>% #only include taxa x seasons where there's decent NPN data
  filter(taxon == "Quercus") %>% 
  ggplot(aes(x = doy, y = pol/max_pol)) + geom_point() + theme_bw() + facet_wrap(site~s_year) +
  geom_line(aes(x = doy, y = mean_prop_flow_m_ma, col = mean_nobs_ma)) +
  scale_x_continuous(limits  = c(70, 150)) +
  scale_color_viridis_c(name = "NPN:\nmean observations/day", direction = -1) + 
  ylab("proportion (airborne pollen & flowers open)")

#example of NPN records from Acer species
npn %>% filter( taxon == "Acer") %>% group_by(species) %>% 
  #filter(phenophase_status == 1) %>% 
  summarize(n_obs = n())

