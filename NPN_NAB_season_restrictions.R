# Original script by Dan Katz
# modifications by T Crimmins
# 8-13-21
#
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
###########################################################################################
# Changes in this version 
# - calculate DOYs 2.5% and 97.% of total pollen count is reached; DOYs 2.5% and 97.5% of "open flowers" is reached
# - use these dates to determine the pollen season (exclude dates outside of this window from further analysis)
# - scale pollen counts from 0-1 by site*taxon*year
#
### set up working environment #############################################################
rm(list=ls())
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
library(MASS)
library(vcd) #install.packages('vcd')
library(AeRobiology)
library(reshape2) #dk: I'm probably going to update instances of reshape2 to dplyr/tidyr when I see them

### load in and prepare NAB data ###############################################################
#nab_raw <- read_csv("~/RProjects/DanK_analyses/NAB_data/NAB_pollen_210621.csv", guess_max = 92013) 
nab_raw <- read_csv("C:/Users/dsk856/Box/things for other people/NAB_NPN/NAB_pollen_210621.csv", guess_max = 92013) 
#dk: we should probably switch over to the 'here' package for this; in the meantime, we can just comment/uncomment this line

### remove Denver rows since no meaningful data in there #####
nab_raw <- filter(nab_raw, site != "Denver") 

### SUBSET TO FLOWER MOUND [or other site] HERE ######
#nab_raw <- subset(nab_raw, site == "Flower Mound")
#nab_raw <- subset(nab_raw, site == "Minneapolis")

#expand to include missing dates
# date_station_grid is a dataframe with a row for each date in the record for each station
date_station_grid <- expand_grid(seq(min(nab_raw$date), max(nab_raw$date), by = '1 day'), unique(nab_raw$site)) %>% 
  `colnames<-`(c("Date", "site")) %>%
  filter(!is.na(site)) %>% 
  ungroup()

# EXPANDED LIST OF GENERA TO INCLUDE ALL TREES IN NAB and NN
nab <- left_join(date_station_grid, nab_raw) %>% 
  dplyr::select(Date, site, 
                Acer, Alnus, Betula, Carpinus.Ostrya, Carya, Celtis, Corylus, Cupressaceae, Fagus, Fraxinus, 
                Juglans, Liquidambar, Myrica, Pinaceae, Populus, Platanus, Quercus, Salix, Tilia, Tsuga, Ulmus) %>% 
  pivot_longer(cols = c(Acer, Alnus, Betula, Carpinus.Ostrya, Carya, Celtis, Corylus, Cupressaceae, Fagus, Fraxinus, 
                        Juglans, Liquidambar, Myrica, Pinaceae, Populus, Platanus, Quercus, Salix, Tilia, Tsuga, Ulmus),
               names_to = "taxon", values_to = "pol") %>% 
  arrange(site, taxon, Date) %>% 
  mutate(Year = year(Date))


## create season definitions based on pollen integral =================================================
#total pollen measured per site/taxon/year
nab_focal_season_max <- nab %>% 
              group_by(site, taxon, Year) %>% 
              summarize(sum_pol = sum(pol, na.rm = TRUE))
            
#position in pollen season
nab_seasons <- nab %>% left_join(., nab_focal_season_max) %>% 
              group_by(site, taxon, Year) %>% 
              mutate(cumu_pol = cumsum(replace_na(pol, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
                     cumu_pol_r = cumu_pol/sum_pol,          #relative sum of pollen
                     in_95season = case_when(cumu_pol_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
                                             cumu_pol_r >= 0.025 & cumu_pol_r <= 0.975~ "in 95% season",
                                             cumu_pol_r > 0.975 ~ "not in 95% season"))  
  
# # rescale pollen counts to 0-1 for each site*taxon*year #dk: I think that any re-scaling should be done at the very end, not here
# nab$Year<-year(nab$Date)
# 
# nab <- nab %>%
#   group_by(site, taxon, Year) %>%
#   mutate(polpct = scales::rescale(pol, to=c(0,1)))

# dk: I'm just going to do it manually since I'd rather know exactly what's happening and doing it for each station is 
# a little awkward and doesn't scale as well, especially because this is running pretty slow 



# # reshape table from long to wide - need to do this for AeRobiology package to work
# nab1 <- dcast(nab, Date + site ~ taxon, value.var = "polpct")
# 
# # AeRobiology only seems to work for one site at a time?? Run through NAB data for each site, generate PS (pollen season) stats one at a time
# Armonk <- subset(nab1, site == "Armonk")
# Armonk$Date <- as.Date(Armonk$Date, format = "%Y-%m-%d")
# Armonk <- subset(Armonk, select = -site)
# ArmonkPS <- calculate_ps(Armonk, method = "percentage", perc = 95)
# ArmonkPS$site <- "Armonk"
# 
# Carrolton <- subset(nab1, site == "Carrolton")
# Carrolton$Date <- as.Date(Carrolton$Date, format = "%Y-%m-%d")
# Carrolton <- subset(Carrolton, select = -site)
# CarroltonPS <- calculate_ps(Carrolton, method = "percentage", perc = 95)
# CarroltonPS$site <- "Carrolton"
# 
# FlowerMound <- subset(nab1, site == "Flower Mound")
# FlowerMound$Date <- as.Date(FlowerMound$Date, format = "%Y-%m-%d")
# FlowerMound <- subset(FlowerMound, select = -site)
# FlowerMoundPS <- calculate_ps(FlowerMound, method = "percentage", perc = 95)
# FlowerMoundPS$site <- "Flower Mound"
# 
# Minneapolis <- subset(nab1, site == "Minneapolis")
# Minneapolis$Date <- as.Date(Minneapolis$Date, format = "%Y-%m-%d")
# Minneapolis <- subset(Minneapolis, select = -site)
# MinneapolisPS <- calculate_ps(Minneapolis, method = "percentage", perc = 95)
# MinneapolisPS$site <- "Minneapolis"
# 
# Springfield <- subset(nab1, site == "Springfield")
# Springfield$Date <- as.Date(Springfield$Date, format = "%Y-%m-%d")
# Springfield <- subset(Springfield, select = -site)
# SpringfieldPS <- calculate_ps(Springfield, method = "percentage", perc = 95)
# SpringfieldPS$site <- "Springfield"
# 
# Waterbury <- subset(nab1, site == "Waterbury")
# Waterbury$Date <- as.Date(Waterbury$Date, format = "%Y-%m-%d")
# Waterbury <- subset(Waterbury, select = -site)
# WaterburyPS <- calculate_ps(Waterbury, method = "percentage", perc = 95)
# WaterburyPS$site <- "Waterbury"
# 
# # paste all of the "sitePS" dfs back together into one df
# nabPS <- rbind(ArmonkPS, CarroltonPS, FlowerMoundPS, MinneapolisPS, SpringfieldPS, WaterburyPS)
# 
# # drop unneeded columns
# nabPS <- subset(nabPS, select = -(ln.ps:daysth))
# 
# # clean up unneeded files
# rm(Armonk)
# rm(ArmonkPS)
# rm(Carrolton)
# rm(CarroltonPS)
# rm(FlowerMound)
# rm(FlowerMoundPS)
# rm(Minneapolis)
# rm(MinneapolisPS)
# rm(Springfield)
# rm(SpringfieldPS)
# rm(Waterbury)
# rm(WaterburyPS)

####### load in and prepare NPN data ###############################################################
#npn_raw <- read_csv("~/RProjects/DanK_analyses/200mibuffer_NNrecords_alltaxa_7-24-21.csv", guess_max = 672676)
npn_raw <- read_csv("C:/Users/dsk856/Box/things for other people/NAB_NPN/200mibuffer_NNrecords_alltaxa_7-24-21.csv", guess_max = 672676)

# smaller test file
#npn_raw <- read_csv("~/RProjects/DanK_analyses/200mibuffer_practice_8-12-21.csv", guess_max = 672676)

filt_tmean_dif <- 2 #filter NPN observations that are within X degrees celsius of the nearest NAB station

npn <- npn_raw %>% 
  filter(tmean_dif > -filt_tmean_dif & tmean_dif < filt_tmean_dif) %>% #MAT filtering
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
                           genus == "Carpinus" ~ "Carpinus.Ostrya",
                           genus == "Ostrya" ~ "Carpinus.Ostrya",
                           genus == "Pinus" ~ "Pinaceae",
                           genus == "Picea" ~ "Pinaceae",
                           genus == "Larix" ~ "Pinaceae", 
                           genus == "Pseudotsuga" ~ "Pinaceae",
                           genus == "Tsuga" ~ "Pinaceae",
                           genus == "Abies" ~ "Pinaceae",
                           TRUE ~ taxon)) %>% 
  dplyr::rename(site = NABStn, 
                date = observation_date) 

#expand to include missing dates
date_station_grid_npn <- expand_grid(seq(min(npn$date), max(npn$date), by = '1 day'), unique(npn$site),
                                     unique(npn$taxon)) %>% 
  `colnames<-`(c("date", "site","taxon")) %>%
  filter(!is.na(site)) %>% 
  ungroup() %>% 
  mutate(phenophase_status = NA,
         flow_prop = NA)

# CALCULATE PROPORTION OF "YES" RECORDS FOR FLOWERING (npn_summary$mean_prop_flow)
#dk : leaving off on this section for the moment. need to figure out how to bring along no-observation days
npn_summary <- 
  bind_rows(date_station_grid_npn, npn) %>% 
  arrange(taxon, site, s_year, date) %>% 
  group_by(taxon, site, s_year, date, doy, date_noyr) %>% 
  dplyr::summarize(mean_flow = mean(phenophase_status),
                   mean_prop_flow = mean(flow_prop),
                   n_obs = as.numeric(n())) %>% 
  ungroup() %>% 
  #left_join(date_station_grid_npn, .) %>% 
  mutate(n_obs = case_when( is.na(n_obs) ~ 0,
                            n_obs > 0 ~ n_obs))


# CALCULATES NUMBER OF NN OBS PER SEASON (YEAR?) - used for filtering out taxa w/few records later on in script
npn_season_summary_nobs <- npn_summary %>% 
  group_by(taxon, s_year, site) %>% 
  dplyr::summarize(nobs_per_season = sum(n_obs))
npn_season_summary_nobs_open_flow <- npn %>% 
  filter(phenophase_status == 1) %>% 
  group_by(taxon, s_year, site) %>% 
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

#add cumulative sum field (by site*taxon*year)
npn_join <- npn_join %>% 
  group_by(site, taxon, s_year) %>% 
  arrange(site, taxon, date) %>% 
  mutate(cum_mean_prop_flow_m_ma = cumsum(mean_prop_flow_m_ma))

#sort by site, taxon, date - maybe not actually necessary
npn_join <- npn_join %>%
  arrange(site, taxon, date)

# figure out DOY that 2.5% and 97.% of max value in cum_mean_prop_flow_m_ma is met for each site*taxon*year
npn_flow_season <- npn_join %>%
  group_by(site, taxon, s_year) %>%
  dplyr::summarize(NPNDOY2.5 = doy[min(which(cum_mean_prop_flow_m_ma >= 0.025*max(cum_mean_prop_flow_m_ma)))],
                   NPNDOY97.5 = doy[min(which(cum_mean_prop_flow_m_ma >= 0.975*max(cum_mean_prop_flow_m_ma)))])

# join season length DOYs for NAB and NN 
NABNPNseason <- merge(npn_flow_season, nabPS, by.x=c("site", "taxon", "s_year"), by.y=c("site", "type", "seasons"))

NABNPNseason <- NABNPNseason %>% 
  rename(
    NABDOY2.5 = st.jd,
    NABDOY97.5 = en.jd
  )

###########################################################################################################
# trim NAB dataframe (nab) using min/max of 2.5 and 97.5DOYs for further analyses
# for now, generate nab dfs using both "inclusive" (bigger window) and "constrained" (smaller window) time frames as a form of sensitivity analyis
# might just pick one route based on these results for generating final results

# add DOY column to nab
nab$DOY <- yday(nab$Date)

# subset using more inclusive criteria (min DOY of NPNDOY2.5 and NABDOY2.5; max DOY of NPNDOY97.5 and NABDOY97.5)
nabmerge <- merge(nab, NABNPNseason, by.x=c("site", "taxon", "Year"), by.y=c("site", "taxon", "s_year"), all.x = TRUE)

# trim using inclusive dates
nab_inclusive_sub <- nabmerge %>%
  group_by(site, taxon, Year) %>%
  filter(DOY > (min(NPNDOY2.5, NABDOY2.5)) & DOY < (max(NPNDOY97.5, NABDOY97.5)))              

# create dataframe using "inclusive" dates
nab_inclusive_sub <- subset(nab_inclusive_sub, select = -c(NPNDOY2.5, NPNDOY97.5, NABDOY2.5, NABDOY97.5, st.dt, en.dt))

# repeat for the inverse - generating dataframe with less inclusive date range for analysis
nab_constrained_sub <- nabmerge %>%
  group_by(site, taxon, Year) %>%
  filter(DOY > (max(NPNDOY2.5, NABDOY2.5)) & DOY < (min(NPNDOY97.5, NABDOY97.5)))   

nab_constrained_sub <- subset(nab_constrained_sub, select = -c(NPNDOY2.5, NPNDOY97.5, NABDOY2.5, NABDOY97.5, st.dt, en.dt))

# trim NPN dataframe (nab) using min/max of 2.5 and 97.5DOYs for further analyses
npnmerge <- merge(npn_join, NABNPNseason, by.x=c("site", "taxon", "s_year"), by.y=c("site", "taxon", "s_year"), all.x = TRUE)

# create dataframe using "inclusive" dates
npn_inclusive_sub <- npnmerge %>%
  group_by(site, taxon, s_year) %>%
  filter(doy > (min(NPNDOY2.5, NABDOY2.5)) & doy < (max(NPNDOY97.5, NABDOY97.5)))              

npn_inclusive_sub <- subset(npn_inclusive_sub, select = -c(NPNDOY2.5, NPNDOY97.5, NABDOY2.5, NABDOY97.5, st.dt, en.dt))

# repeat for the inverse - generating dataframe with less inclusive date range for analysis
npn_constrained_sub <- npnmerge %>%
  group_by(site, taxon, s_year) %>%
  filter(doy > (max(NPNDOY2.5, NABDOY2.5)) & doy < (min(NPNDOY97.5, NABDOY97.5)))   

npn_constrained_sub <- subset(npn_constrained_sub, select = -c(NPNDOY2.5, NPNDOY97.5, NABDOY2.5, NABDOY97.5, st.dt, en.dt))

#####################################
# clean up files 
rm(nabmerge)
rm(nab1)
rm(npn_raw)
rm(npnmerge)
rm(date_station_grid)
rm(date_station_grid_npn)
rm(nab_raw)
rm(nab)
rm(npn)
rm(npnPS)
rm(npn_season_summary_nobs)
rm(npn_season_summary_nobs_open_flow)
rm(npn_summary)
rm(npn_join)
rm(NABNPNseason)
rm(npn_flow_season)

########################### merge NAB & NPN ##########
# inclusive dates
nabnpn_incl <- merge(nab_inclusive_sub, npn_inclusive_sub, by.x=c("site", "taxon", "Date"), by.y=c("site", "taxon", "date"))

# constrained dates
nabnpn_const <- merge(nab_constrained_sub, npn_constrained_sub, by.x=c("site", "taxon", "Date"), by.y=c("site", "taxon", "date"))

# clean up files
rm(nab_constrained_sub)
rm(npn_constrained_sub)
rm(nab_inclusive_sub)
rm(npn_inclusive_sub)

############### histograms, correlations by site & taxon! ###################
# INCLUSIVE DATES RESULTS

Armonkincl <- subset(nabnpn_incl, site == "Armonk")

formula <- y ~ x 
Armonkincl %>% 
  ggplot(aes(y = polpct + 0.01, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Armonkincl %>% 
  filter(taxon == "Acer") %>% 
  ggplot(aes(x = doy, y = mean_prop_flow_m_ma, col = as.factor(Year))) + geom_line() + theme_bw() +facet_wrap(Year~taxon) +
  geom_point(aes(x= doy, y = polpct))


Carroltonincl <- subset(nabnpn_incl, site == "Carrolton")

formula <- y ~ x 
Carroltonincl %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

FlowerMoundincl <- subset(nabnpn_incl, site == "Flower Mound")

formula <- y ~ x 
FlowerMoundincl %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Minneapolisincl <- subset(nabnpn_incl, site == "Minneapolis")

formula <- y ~ x 
Minneapolisincl %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Springfieldincl <- subset(nabnpn_incl, site == "Springfield")

formula <- y ~ x 
Springfieldincl %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Waterburyincl <- subset(nabnpn_incl, site == "Waterbury")

formula <- y ~ x 
Waterburyincl %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

# CONSTRAINED DATES RESULTS
Armonkconst <- subset(nabnpn_const, site == "Armonk")

formula <- y ~ x 
Armonkconst %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Carroltonconst <- subset(nabnpn_const, site == "Carrolton")

formula <- y ~ x 
Carroltonconst %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

FlowerMoundconst <- subset(nabnpn_const, site == "FlowerMound")

formula <- y ~ x 
FlowerMoundconst %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Minneapolisconst <- subset(nabnpn_const, site == "Minneapolis")

formula <- y ~ x 
Minneapolisconst %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Springfieldconst <- subset(nabnpn_const, site == "Springfield")

formula <- y ~ x 
Springfieldconst %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

Waterburyconst <- subset(nabnpn_const, site == "Waterbury")

formula <- y ~ x 
Waterburyconst %>% 
  ggplot(aes(y = polpct, x = mean_prop_flow_m_ma)) + geom_point(alpha = 0.3) + theme_bw() +
  facet_wrap(site~taxon) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10() +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black") +
  xlab("plants flowering \n(proportion; 7 day moving average)") + ylab("scaled pollen count")

