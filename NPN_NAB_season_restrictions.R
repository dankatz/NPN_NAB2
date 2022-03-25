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
#NAB data were assembled in this script: #C:/Users/danka/Box/texas/NAB/extract_pollen_data_from_NPNdata220308.R
nab_raw <- read_csv("C:/Users/danka/Box/Cornell/1 national pollen model/NAB_NPN/NAB_data_request_220308e.csv") %>% 
  mutate(site = NAB_station)%>% 
  rename(dates = Date) 

#nab_raw <- read_csv("data/NAB_pollen_220128c.csv", guess_max = 92013) 

#names(nab_raw)
#test <- filter(nab_raw, site == "NYC")

unique(nab_raw$site)
#remove Denver rows since the data there isn't useful for us
nab_raw <- filter(nab_raw, site != "Denver") 

#expand to include missing dates
# date_station_grid is a dataframe with a row for each date in the record for each station
date_station_grid <- expand_grid(seq(min(nab_raw$dates), max(nab_raw$dates), by = '1 day'), unique(nab_raw$site)) %>% 
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
#npn_raw <- read_csv("data/200mibuffer-inclusive_220128.csv", guess_max = 672676)
npn_raw <- read_csv("C:/Users/danka/Box/Cornell/1 national pollen model/NAB_NPN/npn_200mi_2c_220321.csv")

filt_tmean_dif <- 2 #filter NPN observations that are within X degrees Celsius of the nearest NAB station

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
                site = NAB_station
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
npn_join[c("mean_prop_flow")][is.na(npn_join[c("mean_prop_flow")])] <- 0 

npn_focal_season_integral <- npn_join %>% 
  group_by(site, taxon) %>% 
  summarize(sum_mean_prop_flow = sum(mean_prop_flow, na.rm = TRUE))

#add cumulative sum field (by site*taxon*year) and define 95% season
#switched over to using the proportion of flowers that were open to define season; 
#not using interpolated data or NAs in season definitions
npn_seasons <- npn_join %>% left_join(., npn_focal_season_integral) %>% 
  group_by(site, taxon) %>% 
  arrange(site, taxon, doy) %>% 
  #filter(nobs_yes_per_season > 50) %>% 
  mutate(cumu_flow = cumsum(replace_na(mean_prop_flow, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
         cumu_flow_r = cumu_flow/sum_mean_prop_flow,          #relative sum of pollen
         in_npn_95season = case_when(cumu_flow_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
                                     cumu_flow_r >= 0.025 & cumu_flow_r <= 0.975~ "in 95% season",
                                     cumu_flow_r > 0.975 ~ "not in 95% season")) %>% 
  #manual correction for Cupressaceae season straddling the end of the year
  mutate(in_npn_95season = case_when(taxon == "Cupressaceae" & site == "Carrolton" & doy < 145 ~ "in 95% season", 
                                     taxon == "Cupressaceae" & site == "Carrolton" & doy >= 145 & doy < 345~ "not in 95% season",
                                     taxon == "Cupressaceae" & site == "Carrolton" & doy >= 345 ~ "in 95% season",
                                     TRUE ~ in_npn_95season))

### merge NAB & NPN ####################################################################

#prepare NAB data for joining
nab_seasons_join <- nab_seasons %>% rename(in_pol95season = in_95season)

nabnpn <- left_join(nab_seasons_join, npn_seasons)

unique(nabnpn$taxon)
### data exploration ##################################################################

# #season summaries
# nabnpn %>% dplyr::select(site, taxon, years, sum_pol) %>% distinct() 
# 
# #nab data richness
# nabnpn %>% dplyr::select(site, taxon, years, sum_pol) %>% distinct() %>% 
#   ggplot(aes(x= years, y = sum_pol))+ geom_point() + facet_grid(site~taxon) + theme_bw() + scale_y_log10()
# 
# #npn data richness
# nabnpn %>% dplyr::select(site, taxon, years, nobs_yes_per_season) %>% distinct() %>% 
#   ggplot(aes(x= years, y = nobs_yes_per_season))+ geom_point() + facet_grid(site~taxon) + theme_bw() + scale_y_log10()

#visualize npn data
nabnpn %>% 
  #filter(site == "Armonk") %>% 
  filter(nobs_yes_per_season > 50) %>% 
  #filter(in_npn_95season == "in 95% season") %>% 
  ggplot(aes(x = doy, y = mean_prop_flow_m_ma, group = as.factor(years),
             color = in_npn_95season)) + geom_point() + facet_grid(site~taxon) + theme_bw() 

#visualize nab data
nabnpn %>% 
  filter(site == "Armonk") %>% 
  filter(sum_pol > 200) %>% 
  ggplot(aes(x = doy, y = pol + 1, group = as.factor(years),
             color = in_pol95season)) + geom_point() + facet_grid(site~taxon) + theme_bw()  + scale_y_log10()

## visualize a few time series as examples
nabnpn %>% 
  filter(site == "Waterbury") %>%  #unique(nabnpn$site)
  filter(taxon == "Acer") %>%  #unique(nabnpn$taxon)
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(doy > 80 & doy < 155) %>% 
  ggplot(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = mean_prop_flow_m_ma * 100)) + geom_line(col = "blue") + theme_few() + facet_wrap(~years) +
  geom_point(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = polpct * 100)) + xlab("date") + 
  scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
  theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
  scale_x_date(date_labels = "%b")


nabnpn %>% 
  filter(site == "Armonk") %>%  #unique(nabnpn$site)
  filter(taxon == "Quercus") %>%  #unique(nabnpn$taxon)
  filter(sum_pol > 200) %>% 
  filter(years == 2018 | years == 2019) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(doy > 95 & doy < 175) %>% 
  ggplot(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = mean_prop_flow_m_ma * 100)) + geom_line(col = "blue") + 
  theme_few() + facet_wrap(~years, ncol = 1) +
  geom_point(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = polpct * 100)) + xlab("date") + 
  scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
  theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
  scale_x_date(date_labels = "%b %d")





### Fig. 3: overall comparisons ###############################
#comparing nab and npn data - Pearson's
formula <- y ~ x 
nabnpn %>% 
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma, y = pol + 1)) + 
  geom_point(alpha = 0.5) + facet_grid(site~taxon) + theme_bw()  + scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE, label.x = .9, color = "black")

#comparing nab and npn data - Pearson's - using scaled pollen values
formula <- y ~ x 
nabnpn %>% 
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma, y = polpct)) + theme_bw() + #  scale_y_log10() +
  xlab("observed in flower (%)") + ylab("airborne pollen (% of maximum)") + 
  geom_point(alpha = 0.5) + facet_wrap (site ~ taxon) +#facet_grid(site~taxon) + 
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson") 
  # stat_poly_eq(aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
  #              formula = formula, parse = TRUE, label.x = .9, color = "black")

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
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma * 100, y = polpct * 100)) + 
  geom_point(alpha = 0.3) + facet_grid(site~taxon) + ggthemes::theme_few()  + #scale_y_log10() +
  xlab("observed in flower (%)") + ylab("airborne pollen (% of maximum)") + 
  #geom_smooth(method = "lm") + 
  stat_cor(method = "spearman", cor.coef.name = "rho")


### creating a table of correlations by taxon x site
cor_spear <- nabnpn %>%  
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  #filter(site == "Armonk" & taxon == "Quercus") %>% 
  filter(!is.na(mean_prop_flow_m_ma)) %>% 
  filter(!is.na(polpct)) %>% 
  group_by(site, taxon, years) %>% 
  summarize(n_obs = n(),
            cor_spear = cor(mean_prop_flow_m_ma, polpct, method = "spearman"),
            cor_p_value = cor.test(mean_prop_flow_m_ma, polpct, method = "spearman")$p.value) %>% 
  mutate(cor_spear = round(cor_spear, 2)) %>% 
  arrange(taxon)
cor_spear #unique(cor_spear$taxon)
write_csv(cor_spear, "C:/Users/danka/Box/things for other people/NAB_NPN/spearman_taxon_site_year_220322.csv")
#write.table(cor_spear, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
#dir()

ggplot(cor_spear, aes(x = taxon, y = cor_spear)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = site), width = 0.2) + ggthemes::theme_few() +
  ylab("Spearman correlation between airborne pollen and flowering") +
  theme(axis.text.x = element_text(face = "italic"))

# 
# # creating a similar table without the temperature restriction 
# # (i.e., manually re-running script after changing filt_tmean_dif)
# cor_spear_notemp <- nabnpn %>%  
#   filter(sum_pol > 200) %>% 
#   filter(nobs_yes_per_season > 50) %>% 
#   filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
#   #filter(site == "Armonk" & taxon == "Quercus") %>% 
#   filter(!is.na(mean_prop_flow_m_ma)) %>% 
#   filter(!is.na(polpct)) %>% 
#   group_by(site, taxon, years) %>% 
#   summarize(cor_spear = cor(mean_prop_flow_m_ma, polpct, method = "spearman"),
#             n = n()) %>% 
#   mutate(cor_spear = round(cor_spear, 2)) %>% 
#   arrange(taxon)
# cor_spear_notemp
# write.table(cor_spear_notemp, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# 
# 




### figure for different Acer species over time ###################################
# #exploring species composition for Acer for NPN
# npn_raw %>% 
#   filter(genus == "Acer") %>% 
#   group_by(species) %>% 
#   summarize(n = n())
# 
# head(npn_raw)

acer_ny <- npn_raw %>% 
  mutate(years = year(observation_date)) %>% 
  filter(years == 2016) %>% 
  filter(genus == "Acer") %>% 
  # filter(species == "rubrum" | #species == "saccharum" | #species == "platanoides" | species == "pensylvanicum" |
  #          species == "saccharinum" | 
  #          species == "negundo"
  #        ) %>% 
  filter(NABStn == "Springfield") %>% 
  filter(distNAB < 321869 * 1) %>% #321869 = 200 miles
  filter(day_of_year > 50) %>% 
  filter(day_of_year < 150) %>% 
  arrange(species,  observation_date) %>% 
  group_by(species, day_of_year) %>% 
  #group_by() %>% 
  # mutate(mean_flow = mean(phenophase_status), 
  #        mean_flow_m = round(na_interpolation(mean_flow),1),
  #        mean_flow_m_ma = round(rollmean(mean_flow_m, 7, na.pad=TRUE),2)) 
  dplyr::summarize(mean_flow = mean(phenophase_status),
                   mean_prop_flow = mean(flow_prop),
                   n_obs = sum(!is.na(observation_id)))  %>% #do not include NA values in n() calculation
  filter(n_obs > 3) %>% 
  ungroup() 

acer_ny %>% 
  filter(species != "platanoides") %>% 
  ggplot( aes(x = as.Date(day_of_year, origin = as.Date("2018-01-01")), y = mean_flow, color = species)) + #geom_point() + 
  ggthemes::theme_few() + ylab("flowering (% of observations)") + xlab("") + 
  scale_x_date(limits = c(ymd("2018-03-01"), ymd("2018-05-24"))) + 
  geom_line(aes(x = as.Date(day_of_year, origin = as.Date("2018-01-01")), 
                y=rollmean(mean_flow, 14, na.pad=TRUE))) +#+ facet_wrap(~species)
  scale_color_discrete(labels = c(#expression(italic("negundo")),
                                  expression(italic("Acer rubrum")),
                                  expression(italic("Acer saccharum")))) +
  ggtitle("Maples near Springfield")


#comparing nab and npn data - Spearman's - using scaled pollen values
formula <- y ~ x 
nabnpn %>% 
  filter(sum_pol > 200) %>% 
  filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  filter(taxon == "Acer") %>% 
  ggplot(aes(x = mean_prop_flow_m_ma * 100, y = polpct * 100)) + 
  geom_point(alpha = 0.3) + facet_grid(site~years) + ggthemes::theme_few()  + #scale_y_log10() +
  xlab("observed in flower (%)") + ylab("airborne pollen (% of maximum)") + 
  #geom_smooth(method = "lm") + 
  stat_cor(method = "spearman", cor.coef.name = "rho")


nabnpn %>% 
  filter(site == "Springfield") %>%  #unique(nabnpn$site)
  filter(taxon == "Acer") %>%  #unique(nabnpn$taxon)
  filter(sum_pol > 200) %>% 
  filter(years == 2016 ) %>% 
  filter(nobs_yes_per_season > 30) %>% 
  filter(doy > 50 & doy < 175) %>% 
  ggplot(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = mean_prop_flow_m_ma * 100)) + geom_line(col = "blue") + 
  theme_few() + facet_grid(site~years) +
  geom_point(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = polpct * 100)) + xlab("date") + 
  scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
  theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
  scale_x_date(date_labels = "%b %d")
