# Original script by Dan Katz
# modifications by T Crimmins
# 8-26-21
#
### General information ####################################################################
# Project goals: 
# -assess correlations between NPN and NAB data 
# -provide proof of concept for future continental-scale airborne pollen models
#
# current project participants:
# Dan Katz, Theresa Crimmins, Liz Vogt, Shubhayu Saha
# potential project participants: 
# Arie Managan, Claudia Brown, Ellen Denny, Ambarish (Rish) Vaidyanathan
#
# This script is for data exploration and analyses
# Note: One of the goals is to keep this as simple as possible;
# full phenological models are a future project
# comparison to health outcomes data is a future project
###########################################################################################

### set up working environment #############################################################
#rm(list=ls())
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
library(ggpubr)

#library(MASS)
#library(vcd) #install.packages('vcd')
#library(AeRobiology)

### load in and prepare NAB data ###############################################################
nab_raw <- read_csv("data/NAB_pollen_210621.csv", guess_max = 92013) 

### remove Denver rows since no meaningful data in there #####
nab_raw <- filter(nab_raw, site != "Denver") %>% 
  rename(dates = Date)

#expand to include missing dates
# date_station_grid is a dataframe with a row for each date in the record for each station
date_station_grid <- expand_grid(seq(min(nab_raw$date), max(nab_raw$date), by = '1 day'), unique(nab_raw$site)) %>% 
  `colnames<-`(c("dates", "site")) %>%
  filter(!is.na(site)) %>% 
  ungroup()

# EXPANDED LIST OF GENERA TO INCLUDE ALL TREES IN NAB and NN
nab <- left_join(date_station_grid, nab_raw) %>% 
  dplyr::select(dates, site, 
                Acer, Alnus, Betula, Carpinus.Ostrya, Carya, Celtis, Corylus, Cupressaceae, Fagus, Fraxinus, 
                Juglans, Liquidambar, Myrica, Pinaceae, Populus, Platanus, Quercus, Salix, Tilia, Tsuga, Ulmus) %>% 
  pivot_longer(cols = c(Acer, Alnus, Betula, Carpinus.Ostrya, Carya, Celtis, Corylus, Cupressaceae, Fagus, Fraxinus, 
                        Juglans, Liquidambar, Myrica, Pinaceae, Populus, Platanus, Quercus, Salix, Tilia, Tsuga, Ulmus),
               names_to = "taxon", values_to = "pol") %>% 
  arrange(site, taxon, dates) %>% 
  mutate(years = year(dates),
         ydays = yday(dates))

# rescale pollen counts to 0-1 
nab <- nab %>%
  group_by(site, taxon, years) %>%
  mutate(polpct = scales::rescale(pol, to=c(0,1))) %>%  #for each site*taxon*year
  group_by(site, taxon) %>% 
  mutate(polpct_allyrs = scales::rescale(pol, to=c(0,1))) #for each site*taxon (across years)
  #ggplot(nab, aes(x = polpct_allyrs)) + geom_histogram() + theme_bw() + facet_grid(taxon~site) #graphical check

## create season definitions based on pollen integral =================================================
#total pollen measured per site/taxon/year
#nab_focal_season_max <- nab %>% 
#  group_by(site, taxon, years) %>% 
#  summarize(sum_pol = sum(pol, na.rm = TRUE))

# #total pollen measured per site/taxon/year - using scaled pollen values
# nab_focal_season_max <- nab %>% 
#   group_by(site, taxon, years) %>% 
#   summarize(sum_pctpol = sum(polpct, na.rm = TRUE))

#total pollen measured per site/taxon
nab_focal_season_max <- nab %>% 
  group_by(site, taxon) %>% 
  summarize(sum_pol = sum(pol, na.rm = TRUE))

#position in pollen season
#nab_seasons <- nab %>% left_join(., nab_focal_season_max) %>% 
#  group_by(site, taxon, years) %>% 
#  mutate(cumu_pol = cumsum(replace_na(pol, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
#         cumu_pol_r = cumu_pol/sum_pol,          #relative sum of pollen
#         in_95season = case_when(cumu_pol_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
#                                 cumu_pol_r >= 0.025 & cumu_pol_r <= 0.975~ "in 95% season",
#                                 cumu_pol_r > 0.975 ~ "not in 95% season"))  

# #position in pollen season - using scaled pollen values
# nab_seasons <- nab %>% left_join(., nab_focal_season_max) %>% 
#   group_by(site, taxon, years) %>% 
#   mutate(cumu_pol = cumsum(replace_na(polpct, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
#          cumu_pol_r = cumu_pol/sum_pctpol,          #relative sum of pollen
#          in_95season = case_when(cumu_pol_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
#                                  cumu_pol_r >= 0.025 & cumu_pol_r <= 0.975~ "in 95% season",
#                                  cumu_pol_r > 0.975 ~ "not in 95% season"))  

#position in pollen season - using sum of pollen /taxon/site across all years
nab_seasons <- nab %>% left_join(., nab_focal_season_max) %>% 
  group_by(site, taxon) %>% 
  arrange(site, taxon, ydays, years ) %>% 
  mutate(cumu_pol = cumsum(replace_na(pol, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
         cumu_pol_r = cumu_pol/sum_pol,          #relative sum of pollen
         in_95season = case_when(cumu_pol_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
                                 cumu_pol_r >= 0.025 & cumu_pol_r <= 0.975~ "in 95% season",
                                 cumu_pol_r > 0.975 ~ "not in 95% season"))  


####### load in and prepare NPN data ###############################################################
#npn_raw <- read_csv("data/200mibuffer_8-12-21.csv", guess_max = 672676)
npn_raw <- read_csv("data/200mibuffer-inclusive.csv", guess_max = 672676)

filt_tmean_dif <- 2 #filter NPN observations that are within X degrees celsius of the nearest NAB station

npn <- npn_raw %>% 
  filter(tmean_dif > -filt_tmean_dif & tmean_dif < filt_tmean_dif) %>% #MAT filtering
  mutate(years = year(observation_date),
         doy = yday(observation_date),
         dates_noyr = format(observation_date, format="%m-%d"),
         taxon = genus,
         taxon = case_when(genus == "Juniperus" ~ "Cupressaceae",
                           genus == "Taxodium" ~ "Cupressaceae",
                           genus == "Thuja" ~ "Cupressaceae",
                           genus == "Chamaecyparis" ~ "Cupressaceae",
                           genus == "Metasequoia" ~ "Cupressaceae",
                           genus == "Calocedrus" ~ "Cupressaceae",
                           genus == "Carpinus" ~ "Carpinus.Ostrya",
                           genus == "Ostrya" ~ "Carpinus.Ostrya",
                           genus == "Pinus" ~ "Pinaceae",
                           genus == "Picea" ~ "Pinaceae",
                           genus == "Larix" ~ "Pinaceae", 
                           genus == "Pseudotsuga" ~ "Pinaceae",
                           genus == "Tsuga" ~ "Pinaceae",
                           genus == "Abies" ~ "Pinaceae",
                           TRUE ~ taxon)) %>% 
  dplyr::rename(dates = observation_date,
                site = NABStn
  ) 

#expand to include missing dates
date_station_grid_npn <- expand_grid(seq(min(npn$dates), max(npn$dates), by = '1 day'), unique(npn$site),
                                     unique(npn$taxon)) %>% 
  `colnames<-`(c("dates", "site","taxon")) %>%
  filter(!is.na(site)) %>% 
  ungroup() %>% 
  mutate(years = year(dates),
         phenophase_status = NA,
         flow_prop = NA)

# CALCULATE PROPORTION OF "YES" RECORDS FOR FLOWERING (npn_summary$mean_prop_flow)
npn_summary <- 
  bind_rows(date_station_grid_npn, npn) %>% 
  arrange(taxon, site, years, dates) %>% 
  group_by(taxon, site, years, dates, doy) %>% 
  dplyr::summarize(mean_flow = mean(phenophase_status),
                   mean_prop_flow = mean(flow_prop),
                   n_obs = sum(!is.na(observation_id))) %>% #do not include NA values in n() calculation
  ungroup() %>% 
  mutate(doy = yday(dates), #include some derived date variables that were missing earlier
         date_noyr = format(dates, format="%m-%d"))


# CALCULATES NUMBER OF NN OBS PER YEAR - used for filtering out taxa w/few records later on in script
npn_season_summary_nobs <- npn_summary %>% 
  group_by(taxon, years, site) %>% 
  dplyr::summarize(nobs_per_season = sum(n_obs))
npn_season_summary_nobs_open_flow <- npn %>% 
  filter(phenophase_status == 1) %>% 
  group_by(taxon, years, site) %>% 
  dplyr::summarize(nobs_yes_per_season = n())

# CALCULATES 7-day MOVING AVERAGE FOR NN OPEN FLOWERS
npn_join <-  left_join(npn_summary, npn_season_summary_nobs) %>% 
  left_join(., npn_season_summary_nobs_open_flow) %>% 
  filter(nobs_yes_per_season > 0) %>% 
  mutate(mean_flow_m = round(na_interpolation(mean_flow),1),
         mean_flow_m_ma = round(rollmean(mean_flow_m, 7, na.pad=TRUE),2),
         mean_prop_flow_m = round(na_interpolation(mean_prop_flow), 1),
         mean_prop_flow_m_ma = round(rollmean(mean_prop_flow_m, 7, na.pad=TRUE),2),
         mean_nobs_ma = round(rollmean(n_obs, 7, na.pad=TRUE),2))

##### figure out DOY for 2.5% and 97.% of "open flowers" #############################################
# fill in any NAs in mean_prop_flow_m_ma with 0.0
npn_join[c("mean_prop_flow_m_ma")][is.na(npn_join[c("mean_prop_flow_m_ma")])] <- 0 

npn_focal_season_integral <- npn_join %>% 
  group_by(site, taxon) %>% 
  summarize(sum_mean_prop_flow_m_ma = sum(mean_prop_flow_m_ma, na.rm = TRUE))

#add cumulative sum field (by site*taxon*year)
npn_seasons <- npn_join %>% left_join(., npn_focal_season_integral) %>% 
  group_by(site, taxon) %>% 
  arrange(site, taxon, doy) %>% 
  mutate(cumu_flow = cumsum(replace_na(mean_prop_flow_m_ma, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
         cumu_flow_r = cumu_flow/sum_mean_prop_flow_m_ma,          #relative sum of pollen
         in_npn_95season = case_when(cumu_flow_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
                                     cumu_flow_r >= 0.025 & cumu_flow_r <= 0.975~ "in 95% season",
                                     cumu_flow_r > 0.975 ~ "not in 95% season"))  

### merge NAB & NPN ####################################################################

#prepare NAB data for joining
nab_seasons_join <- nab_seasons %>% rename(in_pol95season = in_95season)

nabnpn <- left_join(nab_seasons_join, npn_seasons)


### data exploration ##################################################################

#season summaries
test <- nabnpn %>% dplyr::select(site, taxon, years, sum_pol) %>% distinct() 

#nab data richness
nabnpn %>% dplyr::select(site, taxon, years, sum_pol) %>% distinct() %>% 
  ggplot(aes(x= years, y = sum_pol))+ geom_point() + facet_grid(site~taxon) + theme_bw() + scale_y_log10()

#npn data richness
nabnpn %>% dplyr::select(site, taxon, years, nobs_yes_per_season) %>% distinct() %>% 
  ggplot(aes(x= years, y = nobs_yes_per_season))+ geom_point() + facet_grid(site~taxon) + theme_bw() + scale_y_log10()

#visualize npn data
nabnpn %>% 
  #filter(site == "Armonk") %>% 
  filter(nobs_yes_per_season > 50) %>% 
  #filter(in_npn_95season == "in 95% season") %>% 
  ggplot(aes(x = doy, y = mean_prop_flow_m_ma, group = as.factor(years),
             color = in_npn_95season)) + geom_line() + facet_grid(site~taxon) + theme_bw() 

#visualize nab data
nabnpn %>% 
  filter(site == "Armonk") %>% 
  filter(sum_pol > 200) %>% 
  ggplot(aes(x = doy, y = pol + 1, group = as.factor(years),
             color = in_pol95season)) + geom_point() + facet_grid(site~taxon) + theme_bw()  + scale_y_log10()

#comparing nab and npn data - Pearson's
formula <- y ~ x 
nabnpn %>% 
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma, y = pol + 1)) + 
  geom_point(alpha = 0.5) + facet_grid(site~taxon) + theme_bw()  + scale_y_log10() +
  geom_smooth(method = "lm") +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black")

#comparing nab and npn data - Pearson's - using scaled pollen values
formula <- y ~ x 
nabnpn %>% 
  # filter(sum_pctpol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma, y = polpct + 1)) + 
  geom_point(alpha = 0.5) + facet_grid(site~taxon) + theme_bw() # + scale_y_log10() +
geom_smooth(method = "lm") +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black")

#comparing nab and npn data - Spearman's  - BUT - the "scale_y_log10" is still in here, do we need to remove that??
formula <- y ~ x 
nabnpn %>% 
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma, y = pol + 1)) + 
  geom_point(alpha = 0.5) + facet_grid(site~taxon) + theme_bw()  + scale_y_log10() +
  geom_smooth(method = "lm") + 
  stat_cor(method = "spearman")

#comparing nab and npn data - Spearman's - using scaled pollen values
formula <- y ~ x 
nabnpn %>% 
  #  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma, y = polpct + 1)) + 
  geom_point(alpha = 0.5) + facet_grid(site~taxon) + theme_bw()  + scale_y_log10() +
  geom_smooth(method = "lm") + 
  stat_cor(method = "spearman")