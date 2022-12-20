### Observations from the USA National Phenology Network can be leveraged to model airborne pollen  ###########
# Aerobiologia, 2022 
# Daniel Katz, Elizabeth Vogt, Arie Manangan, Claudia Brown, Dan Dalan, Kai Zhu, Yiluan Song, and Theresa Crimmins
# 
# For questions about the manuscript or script, please contact Dan Katz: dankatz@cornell.edu
#
# This script includes data assembly, analysis, and visualization. A separate branch ('direct_NAB_NPN_comparison') was created for an offshoot of the 
# analysis that is now being developed as a separate analysis.
# Note: This analysis includes proprietary data from the National Allergy Bureau; only the NPN portion is fully reproducible without access to that data


### set up working environment #######################################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(rnpn)
library(sf)
library(geosphere)
library(prism)
library(readr)
library(here) 
library(purrr)
library(ggplot2)
library(imputeTS)
library(zoo)
library(ggpmisc)
library(viridis)
library(ggthemes)
library(ggpubr)

here::i_am("NPN_NAB_combined_script.R") #using the 'here' package for consistent file structure system
# This script was developed in:
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)

sessionInfo(package = NULL)
### Prepare NPN data #################################################################################
### select top anemophilous taxa from NPN -----------------------------------------
npn_plants <- npn_species(kingdom = "Plantae")


#top woody anemophilous angiosperms
acer_species_list <- c(777,1843,59,778,1,2,1591,60,779,780,3,781,61,1199)
alnus_species_list <- c(62,63,319)
betula_species_list <- c(97, 1439, 98, 430, 1850, 1339, 1851, 99, 1805)
fraxinus_species_list <- c(74,872,873,75,1350)
populus_species_list <- c(1361,320,976,977,1188,27,1481)
quercus_species_list <- c(705,100,1365,757,1870,987,1690,1484,988,316,297,1485,1190,765,1486,
                          301,704,101,1691,1212,989,1366,102,1756,1213,1755,1487,1159,305)
ulmus_species_list <- c(1192,1048,1049,1215,1216)
misc_species_list <- c(935, #Myrica
                       823, #Carpinus
                       1176, 1177, 824, 67, 68, #Carya
                       829, 1342, 1924, 1605, #Celtis
                       71, 72, #Corylus
                       1353, 80, #Juglans
                       81, #Liquidamber
                       2007, #Morus
                       1360, #Olea
                       970, 1211, #Platanus
                       1007, 717, 1875, 1494, 2066, 293, 322, 1493, 1006, 1163, 77, 1371, #salix
                       1009, 1008, 1372, 1010, 1876, #Salix
                       93, 1775, 1776, 1777 #Tilia
)

#herbaceous angiosperms 
herbaceous_species_list <-c(
  441, 1885, 435, 1889, 1005, 1618, #Amaranthaceae 
  145,788,146, #Ambrosia
  1436, 1902, 105, 796, 797, 1900, 798, 1901, #Artemisia
  #Chenopodiaceae
  969, #Plantago
  #Rumex
  1986, 1050)#Urticaceae

#grasses
poaceae_species_list <- filter(npn_plants, family_name == "Poaceae") %>%  
  pull(species_id) 


#pollen cones
pinus_species_list <- c(1629,1686,965,762,50,295,220,1480,219,51,966,967,52,25,968,1687,53,54) #Pinaceae
cupressaceae_species_list <- c(43,1743,1354,289,902,291,290,44, #junipers
                               1723,831,1356,296,1040) #other

list_all_focal_taxa <- c(acer_species_list, alnus_species_list, betula_species_list, fraxinus_species_list, 
                         populus_species_list, quercus_species_list, ulmus_species_list,
                         pinus_species_list, cupressaceae_species_list, 
                         misc_species_list,
                         herbaceous_species_list, poaceae_species_list)

### download and process data -------------------------------------------------------------
npn_direct <- read_csv( here("data", "NPN_220620.csv")) #try reading in data if it's already downloaded
if(exists("npn_direct") == FALSE){ #does the npn_direct object exist? If not, download it:
  npn_direct <- npn_download_status_data(
    request_source = 'Daniel Katz, Cornell and/or Theresa Crimmins',
    species_ids = list_all_focal_taxa,
    years = c(as.character(2009:2021)), #years to include
    phenophase_ids = c(501, 502, 495, 503), #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" 
                                            #conifers: 495 ==  open pollen cones, 503 == pollen release from cones
    additional_fields = c("Observed_Status_Conflict_Flag", "partner_group")
  )
  write_csv(npn_direct, here("data", "NPN_220620.csv"))
}


# get total records per site
sitetotals <- npn_direct %>% 
  group_by(site_id) %>%
  dplyr::summarise(n=n())
sitetotals <- as.data.frame(sitetotals)

# filter out observations with conflict flags
conflict <- npn_direct[which(npn_direct$observed_status_conflict_flag %in% c("MultiObserver-StatusConflict", "OneObserver-StatusConflict")),]  

# use that to get total conflicts per site
conflict_totals <- conflict %>% 
  group_by(site_id) %>%
  dplyr::summarise(nflag=n())
conflict_totals <- as.data.frame(conflict_totals)

# merge the total records and total conflicts (replace NAs with zeroes for sites
# with no conflicts)
hi_conflict <- merge(conflict_totals, sitetotals, by="site_id", all=TRUE)
hi_conflict[is.na(hi_conflict)] <- 0

# calculate percentage of conflict records
hi_conflict$percentage <- hi_conflict$nflag/hi_conflict$n*100
# un-comment the next two lines for a histogram of conflicts per site
# hist(hi_conflict$percentage, xlab="Percentage conflict records", 
#     ylab="Number of sites", main="Frequency of conflicts per site")

# identify sites for which conflict records make up more than 5 percent of total
# records
lowconflict_sites <- hi_conflict[which(hi_conflict$percentage<=5),]
lowconflict_sites <- as.data.frame(lowconflict_sites[,1])
names(lowconflict_sites) <- "site_id"

# remove the high-conflict sites from the previously filtered data 
# (with single conflict flags removed) to see the difference... 
npn_direct_flag <- merge(npn_direct, lowconflict_sites, by="site_id")

npn_flow <- filter(npn_direct_flag, phenophase_id == 501 | phenophase_id == 495)


# flowering intensity value (which is entered about 82% of the time)
npn_active_flow <- npn_flow %>%
  filter(phenophase_status != -1) %>% #removing observations where the observer was unsure whether the phenophase was occurring
  #filter(intensity_value != -9999) %>% #removing observations where the intensity wasn't recorded
  mutate(
    yr = year(observation_date), 
    flow_prop = case_when(
      phenophase_status == 0  ~ 0,
      intensity_value == "Less than 5%" ~ 0.025,
      intensity_value ==  "5-24%"~ (0.05+0.24)/2,
      intensity_value == "25-49%" ~ (0.25+0.49)/2,
      intensity_value == "50-74%" ~ (0.5+0.74)/2,
      intensity_value == "75-94%" ~ (0.75+0.94)/2,
      intensity_value == "95% or more" ~ 0.97,
      intensity_value == "-9999" ~ -9999)) 
hist(npn_active_flow$flow_prop)


#number of observations over time
#npn_active_flow %>% group_by(yr) %>% summarize(n = n())

#number of observations over time by genus
#npn_active_flow %>% group_by(genus, yr) %>% summarize(n = n())

#number of observations by genus
npn_active_flow %>% #filter(flow_prop != 0) %>% 
  group_by(genus) %>% summarize(n = n())

#number of observations by genus for open flowers/pollen cones
npn_active_flow %>% filter(flow_prop != 0) %>% 
  group_by(genus) %>% summarize(n = n())



#observations for pollen release
npn_pol <- filter(npn_direct_flag, phenophase_id == 502 | phenophase_id == 503) %>% 
  filter(phenophase_status != -1)

npn_pol %>% 
  group_by(phenophase_status) %>% 
  dplyr::summarize(n_obs = n())

### SI 1: total observations per taxa ----------------------------------------------------------
SI_1_phenophase_observed <- npn_flow %>% group_by(genus, species) %>% 
  summarize(n_obs_flow_looked_for = n())

SI_1_active_flow_observed <- npn_active_flow %>% 
  filter(intensity_value > 0) %>% 
  group_by(genus, species) %>% 
  summarize(n_obs_flow_active = n())

SI_1 <- left_join(SI_1_phenophase_observed, SI_1_active_flow_observed) %>% 
  mutate(n_obs_flow_active = replace_na(n_obs_flow_active, 0)) %>% 
  mutate(taxon = paste(genus, species, sep = " "))
write_csv(SI_1, "SI_1_220622.csv")


### average temperature of January - April in year of obs for each NPN obs site------------------------------------------------
prism_set_dl_dir("C:/Users/dsk273/Documents/prism")

npn_active_flow_yr_sf <- npn_active_flow %>% 
  dplyr::select(longitude, latitude, site_id, yr) %>% 
  distinct() %>% 
  st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 

#prism_set_dl_dir("~/prism")
#get_prism_monthlys(type = "tmean", years = 2009:2021, keepZip = TRUE, mon = 1:4)

#function to extract time period for each year of data
extract_temp_yr <- function(focal_yr){ #focal_yr <- 2009
  tmean_rast_yr_mo <- prism_archive_subset(temp_period = "monthly", type = "tmean", years = focal_yr, mon = 1:4)
  tmean_rast2_yr_mo <- pd_stack(tmean_rast_yr_mo)
  r_mean <- raster::calc(tmean_rast2_yr_mo, mean) #raster::plot(r_mean)
  
  npn_active_flow_yr_sf_focal_yr <- filter(npn_active_flow_yr_sf, yr == focal_yr)
  tmean_data <- unlist(raster::extract(x = r_mean, 
                                       y = npn_active_flow_yr_sf_focal_yr)) %>% as.data.frame()  
  
  npn_active_flow_yr_sf_focal_yr2 <- npn_active_flow_yr_sf_focal_yr %>%  
    dplyr::select(site_id, yr) %>% 
    mutate(tmean_jan_apr = as.numeric(unlist(tmean_data)))
  
  npn_active_flow_yr_sf_focal_yr2$geometry <- NULL
  print(focal_yr)
  return(npn_active_flow_yr_sf_focal_yr2)
}

#extract_temp_yr(focal_yr = 2010, mo_start = 1, mo_end = 4)
spring_temp_site_yr <- purrr::map_dfr(.x = 2009:2021, .f = extract_temp_yr)
npn_active_flow <- left_join(npn_active_flow, spring_temp_site_yr)


### load in and prepare NAB data ###############################################################
#NAB data were assembled in this script: #C:/Users/danka/Box/texas/NAB/extract_pollen_data_from_NPNdata220308.R
nab_raw <- read_csv(here("data", "NAB_data_request_220308e.csv")) %>% 
  mutate(site = NAB_station)%>% 
  rename(dates = Date) 
#names(nab_raw)
#test <- filter(nab_raw, site == "NYC")

unique(nab_raw$site)
length(unique(nab_raw$site))

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

nab_season_obs_summary <- nab %>% group_by(site, taxon, years) %>% 
        filter(!is.na(pol)) %>% 
        filter(pol != 0) %>% 
        summarize(n_season_non_zero = n())

nab <- left_join(nab, nab_season_obs_summary)

#extract temperature for the NAB stations and add it to active flowers dataframe
NAB_coords <- nab_raw %>% dplyr::select(NAB_station, Lat, Long) %>%
  distinct() %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

NAB_coords_not_sf <- nab_raw %>% dplyr::select(NAB_station, Lat, Long) %>% distinct() %>% rename(site = NAB_station)

#function to extract time period for each year of data
extract_temp_yr_NAB <- function(focal_yr){
  tmean_rast_yr_mo <- prism_archive_subset(temp_period = "monthly", type = "tmean", years = focal_yr, mon = 1:4)
  tmean_rast2_yr_mo <- pd_stack(tmean_rast_yr_mo)
  r_mean <- raster::calc(tmean_rast2_yr_mo, mean)
  
  #npn_active_flow_yr_sf_focal_yr <- filter(npn_active_flow_yr_sf, yr == focal_yr)
  tmean_data <- unlist(raster::extract(x = r_mean, 
                                       y = NAB_coords)) %>% as.data.frame()  
  
  nab_active_flow_yr_sf_focal_yr2 <- NAB_coords %>%  
    mutate(tmean_jan_apr = as.numeric(unlist(tmean_data)),
           yr = focal_yr)
  
  nab_active_flow_yr_sf_focal_yr2$geometry <- NULL
  print(focal_yr)
  return(nab_active_flow_yr_sf_focal_yr2)
}

#extract_temp_yr(focal_yr = 2010, mo_start = 1, mo_end = 4)
spring_temp_site_yr_nab <- purrr::map_dfr(.x = 2009:2021, .f = extract_temp_yr_NAB)%>% 
  rename(years = yr,
         site = NAB_station)

#add in spring air temperature
nab <- left_join(nab, spring_temp_site_yr_nab)



### Fig 1: Nature's Notebook interface and map of observations #####################################################
# This figure was made by Theresa Crimmins and is not included in this script. 
# The re-made publication version of the figure is saved here:
# "C:\Users\dsk273\Box\writing\NAB_NPN short communication\Fig1_photopea_221219.psd"

### Fig 2: spring temperature vs days open Quercus flowers were recorded #####################################################
xlab <- "winter—spring temperature (°C)"
npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubra" | species == "velutina" | species == "alba" | species == "palustris") %>% 
  filter(longitude < -60 & longitude > -90) %>% 
  filter(day_of_year < 181) %>% 
  #filter(flow_prop > 0) %>% 
#  filter(day_of_year < 173) %>% 
 # filter(day_of_year > 81) %>% 
  filter(flow_prop != 0) %>% #change to >0 to exclude the observations where flowering intensity wasn't recorded (i.e., 2009 - 2011)
  ggplot(aes(x = tmean_jan_apr, y = as.Date(day_of_year, origin = as.Date("2018-01-01")))) + geom_point(alpha = 0.1) + ggthemes::theme_few() +
  #facet_wrap(~genus) + 
  #scale_color_viridis_c() + 
  #facet_wrap(~species) +
  geom_smooth(method = "lm", se = FALSE) + 
  xlab(xlab) + ylab("flowering observed (date)") + scale_y_date(date_breaks = "1 month", date_labels =  "%b") +
  #stat_regline_equation(label.y = as.Date(180, origin = as.Date("2018-01-01")), aes(label = ..eq.label..), label.x = 15) +
  
  #the r2 and p value are now hardcoded for greater control, the original dynamic version is here:
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = as.Date(185, origin = as.Date("2018-01-01")), label.x = 15)
  annotate(geom = "text",x = 20, y =  as.Date(200, origin = as.Date("2018-01-01")), 
             label= "italic(r)^2~'='~0.66", parse=TRUE) +
  annotate(geom = "text",x = 20, y =  as.Date(190, origin = as.Date("2018-01-01")), 
           label= "p~'<'~'0.0001'", parse=TRUE) 
  
ggsave(filename = "C:/Users/dsk273/Box/Cornell/writing/NAB_NPN short communication/revision to be submitted October 22/fig2.tiff", dpi = 600, width = 119, height = 119, units = "mm")

#Quercus regression
npn_active_flow_Qu <- npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubrum") %>% 
  filter(longitude < -60 & longitude > -90) %>% 
  filter(day_of_year < 181) %>% 
  filter(flow_prop != 0)  #change to >0 to exclude the observations where flowering intensity wasn't recorded (i.e., 2009 - 2011)
qu_fit_lm <- lm( day_of_year ~ tmean_jan_apr, data = npn_active_flow_Qu)
qu_fit <- summary(qu_fit_lm)
qu_fit
mean(abs(qu_fit_lm$residuals)) #mean absolute error

# qu_tmean_predicted_doy <- tmean_rast2 * qu_fit$coefficients[2] + qu_fit$coefficients[1]
# raster::plot(qu_tmean_predicted_doy)

### Fig SI 5: For each oak species, spring temperature vs flowering observations #####################################################
xlab <- "winter—spring temperature (°C)"
npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubra" | species == "velutina" | species == "alba" | species == "palustris") %>% 
  filter(longitude < -60 & longitude > -90) %>% 
  filter(day_of_year < 181) %>% 
  filter(flow_prop != 0) %>% #change to >0 to exclude the observations where flowering intensity wasn't recorded (i.e., 2009 - 2011)
  ggplot(aes(x = tmean_jan_apr, y = as.Date(day_of_year, origin = as.Date("2018-01-01")), color = species)) + geom_point(alpha = 0.1) + ggthemes::theme_few() +
  geom_smooth(method = "lm", se = FALSE) + 
  xlab(xlab) + ylab("flowering observed (date)") + scale_y_date(date_breaks = "1 month", date_labels =  "%b") +
  theme(legend.text = element_text(face = "italic"))
  #facet_wrap(~species) +
  #stat_regline_equation(label.y = as.Date(180, origin = as.Date("2018-01-01")), aes(label = ..eq.label..), label.x = 15) +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = as.Date(185, origin = as.Date("2018-01-01")), label.x = 15)

#Quercus regression
npn_active_flow_Qu <- npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubrum") %>% 
  filter(longitude < -60 & longitude > -90) %>% 
  filter(day_of_year < 181) %>% 
  filter(flow_prop != 0)  #change to >0 to exclude the observations where flowering intensity wasn't recorded (i.e., 2009 - 2011)
qu_fit_lm <- lm( day_of_year ~ tmean_jan_apr, data = npn_active_flow_Qu)
qu_fit <- summary(qu_fit_lm)

mean(abs(qu_fit_lm$residuals)) #mean absolute error

# qu_tmean_predicted_doy <- tmean_rast2 * qu_fit$coefficients[2] + qu_fit$coefficients[1]
# raster::plot(qu_tmean_predicted_doy)


###SI 2: map residuals for Quercus flowering ~ spring temperature ================================================================================
npn_active_flow_Qu <- npn_active_flow_Qu  %>% 
  mutate(residuals = qu_fit$residuals,
         resid_cat = case_when(residuals < - 28 ~ -28,
                               residuals > 28 ~ 28,
                               TRUE ~ residuals))
ggplot(npn_active_flow_Qu, aes(x = longitude, y = latitude, color = residuals)) + geom_point() + scale_color_viridis_c()

us_boundary <- sf::read_sf("C:/Users/dsk273/Documents/NPN_NAB2/s_22mr22.shp")
us_boundary <- filter(us_boundary, LON < -60 & LON > -92 & LAT > 20)
ggplot(us_boundary) +   geom_sf(data = us_boundary, colour = "black", fill = NA) +
  geom_jitter(aes(x = longitude, y = latitude , col = resid_cat),# size = pollen),#col = hilo2), pollen / max_p
              data = npn_active_flow_Qu, alpha = .7, size = 3, width = 0.2, height = 0.2)  + 
   scale_color_distiller(palette = "Spectral", name = "residuals (days)", breaks = c(-28, -21, -14, -7, 7, 14, 21, 28), labels = c("< -28", "-21", "-14", "-7", "7", "14", "21", ">28")) +
  xlab("") + ylab("") + #theme_few() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  coord_sf(datum=NA) #removes sf induced gridlines



#oak observation sample size
npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubra" | species == "velutina" | species == "alba" | species == "palustris") %>% 
  filter(longitude < -60 & longitude > -90) %>%  filter(day_of_year < 181) %>% filter(flow_prop != 0) %>% summarize(n = n())

npn_active_flow %>% filter(genus == "Quercus") %>%
  filter(longitude < -60 & longitude > -90) %>%  filter(day_of_year < 181) %>% filter(flow_prop != 0) %>% group_by(species) %>% summarize(n = n())

#is the strength of the association different in the north vs south?
npn_active_flow_Qu <- npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubrum") %>% 
  filter(latitude < 40) %>% #switch this to check north vs south
  filter(longitude < -60 & longitude > -90) %>% 
  filter(day_of_year < 181) %>% 
  filter(flow_prop != 0)  #change to >0 to exclude the observations where flowering intensity wasn't recorded (i.e., 2009 - 2011)
summary(lm( day_of_year ~ tmean_jan_apr, data = npn_active_flow_Qu))





### Fig 3: Quercus peak doy NPN vs Quercus peak doy NAB #####################################################
# predicted Quru peak open flow day at each NAB site

peak_qu_npn_pred_nab_sites_data <- nab %>% ungroup() %>% filter(taxon == "Quercus") %>% 
  select(site, taxon, years, tmean_jan_apr) %>% distinct() 
peak_qu_npn_pred_nab_sites_pred_results <- predict.lm(object = qu_fit_lm, newdata = peak_qu_npn_pred_nab_sites_data, se.fit = TRUE, interval = "prediction", level = 0.34)

peak_qu_npn_pred_nab_sites <- nab %>% ungroup() %>% filter(taxon == "Quercus") %>% 
  select(site, taxon, years, tmean_jan_apr) %>% distinct() %>% 
  mutate(quru_peak_day_pred_fit = peak_qu_npn_pred_nab_sites_pred_results$fit[,1],
         quru_peak_day_pred_lwr = peak_qu_npn_pred_nab_sites_pred_results$fit[,2],
         quru_peak_day_pred_upr = peak_qu_npn_pred_nab_sites_pred_results$fit[,3]
         )

# peak_qu_npn_pred_nab_sites <- nab %>% ungroup() %>% filter(taxon == "Quercus") %>% 
#   select(site, taxon, years, tmean_jan_apr) %>% distinct() %>% 
#   mutate(quru_peak_day_pred = tmean_jan_apr * qu_fit$coefficients[2] + qu_fit$coefficients[1])

# empirical Quru peak at each NAB site
NAB_quru_peak_day <- left_join(nab, NAB_coords_not_sf) %>% 
  filter(taxon == "Quercus") %>% filter(years > 2008) %>% 
  filter(Long < -60 & Long > -90) %>% 
  filter(n_season_non_zero  > 14) %>% 
  arrange(site, dates) %>% 
  group_by(site, years)%>% 
  mutate(pol_m = round(na_interpolation(pol),1),
         pol_m_ma = round(rollmean(pol_m, 21, na.pad=TRUE),2),
         NAs_in_window = case_when(is.na(pol) ~ 1,
                                   !is.na(pol) ~ 0),
         NAs_in_window_mean = round(rollmean(NAs_in_window, 21, na.pad=TRUE),2)
  )%>% 
  slice_max(pol_m_ma, with_ties = FALSE) %>% 
  filter(pol_m_ma >100)
NAB_NPN_quru_peak_day <- left_join(NAB_quru_peak_day, peak_qu_npn_pred_nab_sites)

ggplot(NAB_NPN_quru_peak_day, aes(x = as.Date(quru_peak_day_pred_fit, origin = as.Date("2018-01-01")), 
                                  xmin = as.Date(quru_peak_day_pred_lwr, origin = as.Date("2018-01-01")),
                                  xmax = as.Date(quru_peak_day_pred_upr, origin = as.Date("2018-01-01")),
                                  y = as.Date(ydays, origin = as.Date("2018-01-01")))) + 
  geom_smooth(method = "lm", se = FALSE,  lwd =0.6, col = "skyblue3") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_errorbarh(aes(color = as.factor(site)), alpha = 0.2, height = 0) +
  geom_point(aes(color = as.factor(site))) + ggthemes::theme_few() + 
  xlab("modeled peak flowering (date)")+ ylab("peak pollen (date)") +
  scale_color_discrete(guide="none") +
  
  #the r2 and p value are now hardcoded for greater formatting control
  annotate(geom = "text",x =  as.Date(88, origin = as.Date("2018-01-01")), y =  as.Date(143, origin = as.Date("2018-01-01")), 
           label= "italic(r)^2~'='~0.81", parse=TRUE) +
  annotate(geom = "text",x =  as.Date(89, origin = as.Date("2018-01-01")), y =  as.Date(139, origin = as.Date("2018-01-01")), 
           label= "p<~'0.0001'", parse=TRUE) +
  coord_cartesian(xlim = c(as.Date(85, origin = as.Date("2018-01-01")), as.Date(145, origin = as.Date("2018-01-01"))),
                  ylim = c(as.Date(85, origin = as.Date("2018-01-01")), as.Date(145, origin = as.Date("2018-01-01"))))
  

ggsave(filename = "C:/Users/dsk273/Box/Cornell/writing/NAB_NPN short communication/revision to be submitted October 22/fig3.tiff", dpi = 600, width = 119, height = 119, units = "mm")

fit_lm <- lm( ydays ~ quru_peak_day_pred_fit, data = NAB_NPN_quru_peak_day)
fit <- summary(fit_lm)

mean(abs(fit_lm$residuals)) #mean absolute error

#number of station-years
dplyr::select(NAB_NPN_quru_peak_day, site, years) %>% distinct()

#number of stations in case study
length(unique(NAB_NPN_quru_peak_day$site))

#number of stations in NAB dataset
length(unique(nab$site))

#are both variables normally distributed
hist(NAB_NPN_quru_peak_day$quru_peak_day_pred_fit)
shapiro.test(NAB_NPN_quru_peak_day$quru_peak_day_pred_fit)

hist(NAB_NPN_quru_peak_day$ydays)
shapiro.test(NAB_NPN_quru_peak_day$ydays)

car::qqPlot(NAB_NPN_quru_peak_day$ydays, ylab = "test")
car::qqPlot(NAB_NPN_quru_peak_day$quru_peak_day_pred_fit)
qqplot(x = NAB_NPN_quru_peak_day$ydays, y = NAB_NPN_quru_peak_day$quru_peak_day_pred_fit, xlab = "observed flowering date", ylab = "predicted pollen date")


### SI 3: version of Fig 3 by site  #############################################################################################
ggplot(NAB_NPN_quru_peak_day, aes(x = as.Date(quru_peak_day_pred_fit, origin = as.Date("2018-01-01")), 
                                  y = as.Date(ydays, origin = as.Date("2018-01-01"))   )) + #, group = site
  geom_point(aes(color = as.factor(site))) + ggthemes::theme_few() + facet_wrap(~site) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  xlab("peak flowering (date)")+ ylab("peak pollen (date)") +
  scale_color_discrete(guide="none") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.y = as.Date(142, origin = as.Date("2018-01-01")), 
           label.x = as.Date(83, origin = as.Date("2018-01-01")), cor.coef.name = "r", digits = 3)


###SI 4: map residuals for Quercus flowering ~ peak pollen (for fig. 3) ================================================================================
NAB_NPN_quru_peak_day <- NAB_NPN_quru_peak_day %>% ungroup() %>% mutate(resids = fit_lm$residuals)

npn_active_flow_Qu <- npn_active_flow_Qu  %>% 
  mutate(residuals = qu_fit$residuals,
         resid_cat = case_when(residuals < - 28 ~ -28,
                               residuals > 28 ~ 28,
                               TRUE ~ residuals))
ggplot(npn_active_flow_Qu, aes(x = longitude, y = latitude, color = residuals)) + geom_point() + scale_color_viridis_c()

us_boundary <- sf::read_sf("C:/Users/dsk273/Documents/NPN_NAB2/s_22mr22.shp")
us_boundary <- filter(us_boundary, LON < -60 & LON > -92 & LAT > 20)
ggplot(us_boundary) +   geom_sf(data = us_boundary, colour = "black", fill = NA) +
  geom_jitter(aes(x = Long, y = Lat , col = resids),# size = pollen),#col = hilo2), pollen / max_p
              data = NAB_NPN_quru_peak_day, alpha = .7, size = 3, width = 0.3, height = 0.3)  + 
  scale_color_distiller(palette = "Spectral", name = "residuals (days)") +
  xlab("") + ylab("") + #theme_few() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  coord_sf(datum=NA) #removes sf induced gridlines