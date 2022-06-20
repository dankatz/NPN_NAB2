### NPN-NAB comparison manuscript  ####################################################################
# Current project participants:
# Liz Vogt, Theresa Crimmins, Arie Managan, Claudia Brown, Dan Dalan, Kai Zhu, Yiluan Song, and Dan Katz
# This script includes data assembly, analysis, and visualization and is the compilation of several 
# previous versions, also stored in this repo

### set up working environment #######################################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(rnpn)
library(sf)
library(geosphere)
library(prism)
library(readr)
library(here) #install.packages("Rtools")
library(purrr)
library(ggplot2)
library(imputeTS)
library(zoo)
library(ggpmisc)
library(viridis)
library(ggthemes)
library(ggpubr)

#rm(list=ls())

here::i_am("NPN_NAB_combined_script.R") #using the 'here' package for consistent file structure system

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
#filter(npn_plants, genus == "Tilia")
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
  phenophase_ids = c(501, 502,495, 503), #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" #conifers: 495 ==  503 ==
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
# un-comment the next two lines for a histogram of conflicts per site... 
#hist(hi_conflict$percentage, xlab="Percentage conflict records", 
#     ylab="Number of sites", main="Frequency of conflicts per site")

# identify sites for which conflict records make up more than 5 percent of total
# records
lowconflict_sites <- hi_conflict[which(hi_conflict$percentage<=5),]
lowconflict_sites <- as.data.frame(lowconflict_sites[,1])
names(lowconflict_sites) <- "site_id"

# remove the high-conflict sites from the previously filtered data 
# (with single conflict flags removed) to see the difference... 
npn_direct_flag <- merge(npn_direct, lowconflict_sites, by="site_id")

#write_csv(npn_direct, "data/npn_direct_220128.csv") #keeping a local copy to avoid having to re-download data
#npn_direct <- read_csv("data/npn_direct_220128.csv")
npn_flow <- filter(npn_direct_flag, phenophase_id == 501 | phenophase_id == 495)


# flowering intensity value (which is entered about 82% of the time)
npn_active_flow <- npn_flow %>%
  filter(phenophase_status != -1) %>% #removing observations where the observer was unsure whether the phenophase was occurring
  filter(intensity_value != -9999) %>% 
  mutate(
    yr = year(observation_date), 
    flow_prop = case_when(
      phenophase_status == 0  ~ 0,
      intensity_value == "Less than 5%" ~ 0.025,
      intensity_value ==  "5-24%"~ (0.05+0.24)/2,
      intensity_value == "25-49%" ~ (0.25+0.49)/2,
      intensity_value == "50-74%" ~ (0.5+0.74)/2,
      intensity_value == "75-94%" ~ (0.75+0.94)/2,
      intensity_value == "95% or more" ~ 0.97)) 
hist(npn_active_flow$flow_prop)

#observations for pollen release
npn_pol <- filter(npn_direct_flag, phenophase_id == 502 | phenophase_id == 503) %>% 
  filter(phenophase_status != -1)

npn_pol %>% 
  group_by(phenophase_status) %>% 
  dplyr::summarize(n_obs = n())


### extract mean annual air temperature for each NPN observation site: flowers -----------------------------
### Mean annual temperature 
npn_active_flow_sf <- npn_active_flow %>% 
  dplyr::select(longitude, latitude, site_id) %>% 
  distinct() %>% 
  st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 


#prism_set_dl_dir("~/prism") 
prism_set_dl_dir("C:/Users/dsk273/Documents/prism")
get_prism_normals(type = "tmean", resolution = "4km", annual = TRUE)

#choosing which rasters need to be included
tmean_rast <- prism_archive_subset(type = "tmean", temp_period = "annual normals", resolution = "4km")
tmean_rast2 <- pd_stack(tmean_rast)

#extract the data from prism raster(s) #note: allowing the raster package to handle the CRS transformation
tmean_data <- unlist(raster::extract(x = tmean_rast2, 
                                     y = npn_active_flow_sf)) %>% as.data.frame() 

npn_active_flow3 <- npn_active_flow_sf %>%  
  dplyr::select(site_id) %>% 
  mutate(tmean = as.numeric(unlist(tmean_data)))

npn_active_flow3$geometry <- NULL
npn_active_flow <- left_join(npn_active_flow, npn_active_flow3)

### average temperature of January - April ------------------------------------------------
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
spring_temp_site_yr <- purrr::map_dfr(.x = 2011:2021, .f = extract_temp_yr)
npn_active_flow <- left_join(npn_active_flow, spring_temp_site_yr)





### extract mean annual air temperature for each NPN observation site: pollen -----------------------------
# npn_pol_sf <- npn_pol %>% 
#   dplyr::select(longitude, latitude, site_id) %>% 
#   distinct() %>% 
#   st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 
# 
# 
# prism_set_dl_dir("~/prism")
# get_prism_normals(type = "tmean", resolution = "4km", annual = TRUE)
# 
# #choosing which rasters need to be included
# tmean_rast <- prism_archive_subset(type = "tmean", temp_period = "annual normals", resolution = "4km")
# tmean_rast2 <- pd_stack(tmean_rast)
# 
# #extract the data from prism raster(s) #note: allowing the raster package to handle the CRS transformation
# tmean_data <- unlist(raster::extract(x = tmean_rast2, 
#                                      y = npn_pol_sf)) %>% as.data.frame() 
# 
# npn_pol_sf3 <- npn_pol_sf %>%  
#   dplyr::select(site_id) %>% 
#   mutate(tmean = as.numeric(unlist(tmean_data)))
# 
# npn_pol_sf3$geometry <- NULL
# npn_pol <- left_join(npn_pol, npn_pol_sf3)

### extract mean annual air temperature for each selected NAB site --------------------------------
#NAB data were assembled in this script: #C:/Users/danka/Box/texas/NAB/extract_pollen_data_from_NPNdata220308.R
NAB <- read_csv(here("data", "NAB_data_request_220308e.csv")) #contact Dan Katz if you have any questions about this dataset

#extract temperature for the NAB stations and add it to active flowers dataframe
NAB_coords <- NAB %>% dplyr::select(NAB_station, Lat, Long) %>% 
  distinct() %>% 
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) 

tmean_data_NAB <- unlist(raster::extract(x = tmean_rast2, #matrix(c(NAB_tx_coords$long, NAB_tx_coords$lat), ncol = 2), 
                                         y = NAB_coords )) %>% as.data.frame() 
NAB_coords_tmean <- NAB_coords %>% mutate(tmean_NAB = as.numeric(unlist(tmean_data_NAB))) 
NAB_coords_tmean$geometry <- NULL




### extract spring mean air temperature for each selected NAB site --------------------------------
#prism_set_dl_dir("~/prism")
#get_prism_monthlys(type = "tmean", years = 2009:2021, keepZip = TRUE, mon = 1:4)

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
spring_temp_site_yr_nab <- purrr::map_dfr(.x = 2009:2021, .f = extract_temp_yr_NAB)





### calculate geographic distance from each NPN site to each NAB site -----------------------------
NAB_coords_notsf <- NAB %>% dplyr::select(NAB_station, NAB_station_lat = Lat, NAB_station_long = Long) %>% distinct() 
#NAB_coords_notsf <- NAB_coords_notsf[1:2,] #for testing with a subset of stations 
npn_active_flow_coords_only <- npn_active_flow %>% dplyr::select(longitude, latitude) %>% as.matrix(.)
npn_coords_all <- as.matrix(npn_active_flow[,c(6,5)])

#function for extracting distance
dist_calc_fun <- function(NAB_station, NAB_station_long, NAB_station_lat){
  #NAB_coords_row <- NAB_coords_notsf[NAB_station_row,] #NAB_coords_row <- NAB_coords_notsf[1,]
  #NAB_coords_focal <- c(NAB_station_long, NAB_station_lat) #NAB_coords_row[,3:2] #make sure to have it in the order of long, lat
  npn_active_flow_dist_i <- distm( x = npn_active_flow_coords_only, 
                                   y = c(NAB_station_long, NAB_station_lat), 
                                   fun = distHaversine)
  
  npn_active_flow_i <- npn_active_flow %>% mutate(distNAB = as.numeric(npn_active_flow_dist_i),
                                                  NAB_station = NAB_station)
  npn_active_flow_i_200m <- filter(npn_active_flow_i, distNAB <= 400000) #filter observations within x meters to keep size down
  return(npn_active_flow_i_200m)
}

#Sys.time()
NPN_near_NAB <- pmap_dfr(NAB_coords_notsf, dist_calc_fun) #takes ~ 30 min on laptop
#Sys.time()

#ggplot(NPN_near_NAB, aes(x = distNAB)) + geom_histogram() + facet_wrap(~NAB_station)
#ggplot(NPN_near_NAB, aes(x = longitude, y = latitude)) + geom_point() 


### add in average annual temperature for NAB stations and calculate difference ---------------------
NPN_near_NAB2 <- left_join(NPN_near_NAB, NAB_coords_tmean) %>% 
  mutate(tmean_dif = as.numeric(tmean_NAB - tmean))



### export data to file (data exploration is next script) --------------------------------------------
write_csv(NPN_near_NAB2, here("data", "NPN_near_NAB_220427.csv"))

# readr::write_csv(NPN_near_NAB2, "data/200mibuffer-inclusive_220321.csv")
# write.csv(NPN_near_NAB2, "data/200mibuffer-inclusive_220321b.csv")

#double check that everything added up well
NPN_near_NAB2 %>% 
  group_by(NAB_station) %>% 
  filter(genus == "Quercus") %>%   
  summarize(n = n())




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

nab <- nab %>% #filter(taxon == "Ulmus") %>% 
  mutate(taxon = case_when(taxon == "Ulmus" & ydays > 175 & site != "Santa Barbara" ~ "Ulmus_fall", 
         TRUE ~ taxon))


# rescale pollen counts to 0-1 
nab <- nab %>%
  group_by(site, taxon, years) %>%
  mutate(polpct = scales::rescale(pol, to=c(0,1))) %>%  #for each site*taxon*year
  group_by(site, taxon) %>% 
  mutate(polpct_allyrs = scales::rescale(pol, to=c(0,1))) #for each site*taxon (across years)
#ggplot(nab, aes(x = polpct_allyrs)) + geom_histogram() + theme_bw() + facet_grid(taxon~site) #graphical check


#add in spring air temperature
spring_temp_site_yr_nab <- spring_temp_site_yr_nab %>% rename(years = yr,
                                                              site = NAB_station)
nab <- left_join(nab, spring_temp_site_yr_nab)

## create season definitions based on pollen integral =================================================

#total pollen measured per site/taxon (across all years)
nab_focal_all_years <- nab %>% 
  group_by(site, taxon) %>% 
  summarize(sum_pol_all_yrs = sum(pol, na.rm = TRUE))

#total pollen measured per site/taxon/year
nab_focal_season_sum <- nab %>% 
  group_by(site, taxon, years) %>% 
  summarize(sum_pol_season = sum(pol, na.rm = TRUE))

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

#Temporary solution for Cupressaceae season definitions
#for custom fixes of Cupressaceae seasons that straddle the calendar year
TX_sites <- c("Austin", "Dallas", "Flower Mound", "Houston", "San Antonio B", "Waco A", "Waco B", "Oklahoma City A", "Phoenix",
              "Tulsa A")
west_coast_sites <- c("Eugene", "Pleasanton", "San Diego", "San Jose")


#position in pollen season - using sum of pollen /taxon/site across all years
nab_seasons <- nab %>% 
  left_join(., nab_focal_all_years) %>% 
  left_join(., nab_focal_season_sum) %>% 
  group_by(site, taxon) %>% 
  arrange(site, taxon, ydays, years ) %>% 
  mutate(cumu_pol = cumsum(replace_na(pol, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
         cumu_pol_r = cumu_pol/sum_pol_all_yrs,          #relative sum of pollen
         in_95season = case_when(cumu_pol_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
                                 cumu_pol_r >= 0.025 & cumu_pol_r <= 0.975~ "in 95% season",
                                 cumu_pol_r > 0.975 ~ "not in 95% season"),
         in_99polseason = case_when(cumu_pol_r < 0.005 ~ "not in 99% season", #is the observation in the 95% pollen season?
                                    cumu_pol_r >= 0.005 & cumu_pol_r <= 0.995~ "in 99% season",
                                    cumu_pol_r > 0.995 ~ "not in 99% season"))%>% 
  #manual correction for Cupressaceae season straddling the end of the year: Texas/SW sites
  mutate(in_95season = case_when(
    taxon == "Cupressaceae" & site %in% TX_sites & ydays < 125 ~ "in 95% season", 
    taxon == "Cupressaceae" & site %in% TX_sites & ydays >= 125 & ydays < 275 ~ "not in 95% season",
    taxon == "Cupressaceae" & site %in% TX_sites & ydays >= 275 ~ "in 95% season",
    TRUE ~ in_95season )) %>% 
  #manual correction for Cupressaceae season straddling the end of the year: west coast sites
  mutate(in_95season = case_when(
    taxon == "Cupressaceae" & site %in% west_coast_sites & ydays < 200 ~ "in 95% season", 
    taxon == "Cupressaceae" & site %in% west_coast_sites & ydays >= 200 & ydays < 300 ~ "not in 95% season",
    taxon == "Cupressaceae" & site %in% west_coast_sites & ydays >= 300 ~ "in 95% season",
    TRUE ~ in_95season )) %>% 
  #same but for 99% season
  mutate(in_99polseason = case_when(
    taxon == "Cupressaceae" & site %in% TX_sites & ydays < 115 ~ "in 99% season", 
    taxon == "Cupressaceae" & site %in% TX_sites & ydays >= 115 & ydays < 255 ~ "not in 99% season",
    taxon == "Cupressaceae" & site %in% TX_sites & ydays >= 255 ~ "in 99% season",
    TRUE ~ in_99polseason )) %>% 
  #manual correction for Cupressaceae season straddling the end of the year: west coast sites
  mutate(in_99polseason = case_when(
    taxon == "Cupressaceae" & site %in% west_coast_sites & ydays < 180 ~ "in 99% season", 
    taxon == "Cupressaceae" & site %in% west_coast_sites & ydays >= 180 & ydays < 290 ~ "not in 99% season",
    taxon == "Cupressaceae" & site %in% west_coast_sites & ydays >= 290 ~ "in 99% season",
    TRUE ~ in_99polseason ))



# #visual checks of NAB season definitions
nab_seasons %>%
  filter(taxon == "Cupressaceae") %>%
  #filter(taxon == "Ulmus") %>%
  #filter(years > 2011) %>%
  #filter(sum_pol_season > 200) %>%
  #filter(nobs_yes_per_season > 50) %>%
  #filter(in_npn_95season == "in 95% season") %>%
  ggplot(aes(x = ydays, y = polpct, group = as.factor(years),
             color = in_99polseason)) + geom_point() + facet_wrap(~site) + theme_bw()
# 
# #visual checks of NAB season definitions: 99%
# nab_seasons %>% 
#   filter(taxon == "Quercus") %>% 
#   #filter(years > 2011) %>% 
#   filter(sum_pol_season > 200) %>% 
#   #filter(nobs_yes_per_season > 50) %>% 
#   #filter(in_npn_95season == "in 95% season") %>% 
#   ggplot(aes(x = ydays, y = polpct, group = as.factor(years),
#              color = in_99polseason)) + geom_point() + facet_wrap(~site) + theme_bw() 


####### load in and prepare NPN data ###############################################################
npn_raw <- read_csv(here("data", "NPN_near_NAB_220427.csv"))
  
filt_tmean_dif <- 2 #filter NPN observations that are within X degrees Celsius of the nearest NAB station

npn <- npn_raw %>% 
  filter(tmean_dif > -filt_tmean_dif & tmean_dif < filt_tmean_dif) %>% #MAT filtering
  filter(distNAB < (200 * 1000)) %>% #filter by distance from NAB; needs it in meters #160934 = 100 miles, 321869 = 200 miles, 482803 = 300 miles
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

#hist(npn$distNAB)


#expand to include missing dates
date_station_grid_npn <- expand_grid(seq(min(npn$dates), max(npn$dates), by = '1 day'), unique(npn$site),
                                     unique(npn$taxon)) %>% 
  `colnames<-`(c("dates", "site","taxon")) %>%
  filter(!is.na(site)) %>% 
  ungroup() %>% 
  mutate(years = year(dates),
         phenophase_status = NA,
         flow_prop = NA)

# CALCULATE PROPORTION OF "YES" RECORDS FOR FLOWERING for each taxon x site x date (npn_summary$mean_prop_flow)
npn_summary <- 
  bind_rows(date_station_grid_npn, npn) %>% 
  arrange(taxon, site, years, dates) %>% 
  group_by(taxon, site, years, dates) %>% 
  dplyr::summarize(mean_flow = mean(phenophase_status, na.rm = TRUE),
                   mean_prop_flow = mean(flow_prop, na.rm = TRUE),
                   n_obs = sum(!is.na(observation_id)),
                   tmean = mean(tmean, na.rm = TRUE),
                   tmean_NAB = mean(tmean_NAB, na.rm = TRUE),
                   distNAB_mean = mean(distNAB, na.rm = TRUE),
                   lat_mean = mean(latitude, na.rm = TRUE),
                   long_mean = mean(longitude, na.rm = TRUE)) %>% #do not include NA values in n() calculation
  ungroup() %>% 
  mutate(doy = yday(dates), #include some derived date variables that were missing earlier
         date_noyr = format(dates, format="%m-%d"))


# CALCULATE NUMBER OF NN OBS PER YEAR - used for filtering out taxa w/few records later on in script
npn_season_summary_nobs <- npn_summary %>% 
  group_by(taxon, years, site) %>% 
  dplyr::summarize(nobs_per_season = sum(n_obs))

npn_season_summary_nobs_open_flow <- npn %>% 
  filter(phenophase_status == 1) %>% 
  group_by(taxon, years, site) %>% 
  dplyr::summarize(nobs_yes_per_season = n())

# CALCULATE 7-day MOVING AVERAGE FOR NN OPEN FLOWERS
npn_join <-  left_join(npn_summary, npn_season_summary_nobs) %>% 
  left_join(., npn_season_summary_nobs_open_flow) %>% 
  filter(nobs_yes_per_season > 0) %>% 
  mutate(mean_flow_m = round(na_interpolation(mean_flow),1),
         mean_flow_m_ma = round(rollmean(mean_flow_m, 7, na.pad=TRUE),2),
         mean_prop_flow_m = round(na_interpolation(mean_prop_flow), 1),
         mean_prop_flow_m_ma = round(rollmean(mean_prop_flow_m, 7, na.pad=TRUE),2),
         mean_prop_flow_ma = round(rollapply(mean_prop_flow, width=7, 
                                             FUN=function(x) mean(x, na.rm=TRUE), by=1, partial=TRUE, fill=NA),2),
         mean_nobs_ma = round(rollmean(n_obs, 7, na.pad=TRUE),2))


#figure out the unique number of observers per taxon x site x year
npn_nobs_season <- npn %>% 
  dplyr::select(site, taxon, years, individual_id) %>% 
  distinct()%>% 
  group_by(site, taxon, years) %>% 
  summarize(unique_observers = n()) 



##### season definitions for NPN  #############################################
# fill in any NAs in mean_prop_flow_m_ma with 0.0
npn_join[c("mean_prop_flow")][is.na(npn_join[c("mean_prop_flow")])] <- 0 


# #define the pollen season within a year
# npn_focal_season_integral <- npn_join %>% 
#   group_by(site, taxon, years) %>% 
#   summarize(sum_mean_prop_flow = sum(mean_prop_flow, na.rm = TRUE))

#for defining the pollen season across years
npn_focal_season_integral <- npn_join %>% 
  filter(nobs_yes_per_season > 10) %>% 
  group_by(site, taxon) %>% 
  summarize(sum_mean_prop_flow = sum(mean_prop_flow, na.rm = TRUE))


#add cumulative sum field (by site*taxon) and define 95% & 99%season across all years
#switched over to using the proportion of flowers that were open to define season; 
#not using interpolated data or NAs in season definitions
TX_sites <- c("Austin", "Dallas", "Flower Mound", "Houston", "San Antonio B", "Waco A", "Waco B")

npn_seasons <- npn_join %>% left_join(., npn_focal_season_integral) %>% 
  left_join(., npn_nobs_season) %>% 
  group_by(site, taxon) %>%  #Could add years here to convert to season definition within year
  arrange(site, taxon, doy, years) %>% 
  filter(nobs_yes_per_season > 10) %>% 
  mutate(cumu_flow = cumsum(replace_na(mean_prop_flow, 0)),  #cumulative sum of pollen #putting NAs as 0s for calculating seasonal sums 
         cumu_flow_r = cumu_flow/sum_mean_prop_flow,          #relative sum of pollen
         in_npn_95season = case_when(cumu_flow_r < 0.025 ~ "not in 95% season", #is the observation in the 95% pollen season?
                                     cumu_flow_r >= 0.025 & cumu_flow_r <= 0.975~ "in 95% season",
                                     cumu_flow_r > 0.975 ~ "not in 95% season"),
         in_npn_99season = case_when(cumu_flow_r < 0.005 ~ "not in 99% season", #is the observation in the 99% pollen season?
                                     cumu_flow_r >= 0.005 & cumu_flow_r <= 0.995~ "in 99% season",
                                     cumu_flow_r > 0.995 ~ "not in 99% season")) %>% 
  #manual correction for Cupressaceae season straddling the end of the year
  mutate(in_npn_95season = case_when(taxon == "Cupressaceae" & site %in% TX_sites & doy < 125 ~ "in 95% season", 
                                     taxon == "Cupressaceae" & site %in% TX_sites & doy >= 125 & doy < 330 ~ "not in 95% season",
                                     taxon == "Cupressaceae" & site %in% TX_sites & doy >= 330 ~ "in 95% season",
                                     TRUE ~ in_npn_95season),
         in_npn_99season = case_when(taxon == "Cupressaceae" & site %in% TX_sites & doy < 120 ~ "in 99% season", 
                                     taxon == "Cupressaceae" & site %in% TX_sites & doy >= 120 & doy < 310 ~ "not in 99% season",
                                     taxon == "Cupressaceae" & site %in% TX_sites & doy >= 310 ~ "in 99% season",
                                     TRUE ~ in_npn_99season))

#visual checks of NPN season definitions
# npn_seasons %>%
#   filter(taxon == "Quercus") %>%
#   #filter(nobs_yes_per_season > 50) %>%
#   #filter(in_npn_95season == "in 95% season") %>%
#   ggplot(aes(x = doy, y = mean_prop_flow_ma, group = as.factor(years),
#              color = in_npn_95season)) + geom_point() + facet_grid(site~years) + theme_bw()
# 
npn_seasons %>%
  filter(taxon == "Ulmus") %>%
  #filter(nobs_yes_per_season > 50) %>%
  #filter(in_npn_95season == "in 95% season") %>%
  ggplot(aes(x = doy, y = mean_prop_flow_m_ma, group = as.factor(years),
             color = in_npn_99season)) + geom_point() + facet_grid(site~years) + theme_bw()




### merge NAB & NPN ####################################################################

#prepare NAB data for joining
nab_seasons_join <- nab_seasons %>% rename(in_pol95season = in_95season)
nabnpn <- left_join(nab_seasons_join, npn_seasons) 

#create labels that will easily allow genus (but not families!) to be italicized using the ggtext package 
nabnpn <- nabnpn %>% mutate(taxon_labs = case_when(taxon == "Cupressaceae" ~ "Cupressaceae",
                                                   taxon == "Pinaceae" ~ "Pinaceae",
                                                   TRUE ~ paste0("<i>", taxon,"</i>")),
                            taxon_labs = as.factor(taxon_labs))

 # check italicizing
# nabnpn %>% sample_n(1) %>% 
#   ggplot(aes(x = taxon, y = years)) + geom_point() +
#   scale_x_discrete(labels = levels(nabnpn$taxon_labs)) +
#   theme(axis.text.x = ggtext::element_markdown())



### misc data exploration ##################################################################
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
# 
# #visualize npn data
# nabnpn %>% 
#   #filter(site == "Armonk") %>% 
#   filter(nobs_yes_per_season > 50) %>% 
#   #filter(in_npn_95season == "in 95% season") %>% 
#   ggplot(aes(x = doy, y = mean_prop_flow_m_ma, group = as.factor(years),
#              color = in_npn_95season)) + geom_point() + facet_grid(site~taxon) + theme_bw() 
# 
# nabnpn %>% 
#   filter(taxon == "Quercus") %>% 
#   filter(nobs_yes_per_season > 20) %>% 
#   #filter(in_npn_95season == "in 95% season") %>% 
#   ggplot(aes(x = doy, y = mean_prop_flow_m_ma, group = as.factor(years),
#              color = in_npn_95season)) + geom_point() + facet_grid(site~years) + theme_bw() 
# 
# #visualize nab data
# nabnpn %>% 
#   filter(site == "Armonk") %>% 
#   filter(sum_pol > 200) %>% 
#   ggplot(aes(x = doy, y = pol + 1, group = as.factor(years),
#              color = in_pol95season)) + geom_point() + facet_grid(site~taxon) + theme_bw()  + scale_y_log10()
# 
# ## visualize a few time series as examples
# nabnpn %>% 
#   filter(site == "Waterbury") %>%  #unique(nabnpn$site)
#   filter(taxon == "Quercus") %>%  #unique(nabnpn$taxon)
#   filter(sum_pol_season > 200) %>% 
#   filter(nobs_yes_per_season > 50) %>% 
#   filter(doy > 80 & doy < 155) %>% 
#   ggplot(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = mean_prop_flow_ma * 100)) + geom_line(col = "blue") + theme_few() + facet_wrap(~years) +
#   geom_point(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = polpct * 100)) + xlab("date") + 
#   scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
#   theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
#   scale_x_date(date_labels = "%b%d")



### Fig 2: examples of time series and correlation #############################################################

# panel A: Armonk Quercus time series 2018
#add another example here
fig_site <- "Armonk"
fig_taxon <- "Quercus"
fig_year <- 2018
fig_seasons <- nabnpn %>% filter(site == fig_site & taxon == fig_taxon & years == fig_year) %>%  
  arrange(taxon, site, years, dates)
fig_season_pol_start <- fig_seasons$dates[min(which (fig_seasons$in_pol95season == "in 95% season"))]
fig_season_pol_end   <- fig_seasons$dates[max(which (fig_seasons$in_pol95season == "in 95% season"))]
# fig_season_npn_start <- fig_seasons$dates[min(which (fig_seasons$in_npn_95season == "in 95% season"))]
# fig_season_npn_end   <- fig_seasons$dates[max(which (fig_seasons$in_npn_95season == "in 95% season"))]
fig_season_pol_start <- fig_seasons$dates[min(which (fig_seasons$in_99polseason == "in 99% season"))]
fig_season_pol_end   <- fig_seasons$dates[max(which (fig_seasons$in_99polseason == "in 99% season"))]
fig_season_npn_start <- fig_seasons$dates[min(which (fig_seasons$in_npn_99season == "in 99% season"))]
fig_season_npn_end   <- fig_seasons$dates[max(which (fig_seasons$in_npn_99season == "in 99% season"))]


# correlation
panel_a_cor <- 
  nabnpn %>% 
  filter(site == fig_site) %>%  #unique(nabnpn$site)
  filter(taxon == fig_taxon) %>%  #unique(nabnpn$taxon)
  filter(years == fig_year) %>% 
  #filter(in_pol95season == "in 95% season" & in_npn_95season == "in 95% season") %>% 
  filter(in_99polseason == "in 99% season" & in_npn_99season == "in 99% season") %>% 
  #filter(nobs_yes_per_season > 50) %>% 
  #filter(doy > 95 & doy < 175) %>% 
  summarise(out = cor(polpct, mean_prop_flow_ma, use = "complete.obs", method = "spearman")) 
# ggplot(aes(y = polpct * 100, x = mean_prop_flow_ma * 100)) + geom_point(alpha = 0.7) + 
#   theme_few() + ylab("airborne pollen (% of maximum)") + xlab("flowering (% of observations)")


panel_a <- nabnpn %>% 
  filter(site == fig_site) %>%  #unique(nabnpn$site)
  filter(taxon == fig_taxon) %>%  #unique(nabnpn$taxon)
  filter(years == fig_year) %>% 
  # filter(nobs_yes_per_season > 50) %>% 
  # filter(sum_pol_season > 200) %>% 
  filter(doy > 95 & doy < 160) %>% 
  ggplot(aes(x = dates, y = mean_prop_flow_ma * 100)) + geom_line(col = "blue") + 
  theme_few() + #facet_wrap(~years) +
  geom_point(aes(x = dates, y = polpct * 100), alpha = 0.3) + xlab("date") + 
  scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
  annotate("text", x= ymd_hms("2018/6/05 00:00:00"), y= 98, label= paste0("\U03C1 = ", round(panel_a_cor[3], 2))) + 
  theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
  #scale_x_date(date_labels = "%b %d") +
  # geom_segment(x = fig_season_pol_start, xend = fig_season_pol_end,#pollen 95% season line
  #              y = -2, yend = -2, col = "black", lwd = 2)  +
  # geom_segment(x = fig_season_npn_start, xend = fig_season_npn_end,#npn 95% season line
  #              y = -3, yend = -3, col = "blue", lwd = 2)  +
  geom_line(aes(y= zoo::rollmean(polpct * 100, 7, na.pad=TRUE)), col = "gray20") 




# panel B: Springfield Acer time series 2016
fig_site <- "Springfield"
fig_taxon <- "Acer"
fig_year <- 2016
fig_seasons <- nabnpn %>% filter(site == fig_site & taxon == fig_taxon & years == fig_year) %>%  
  arrange(taxon, site, years, dates)
#fig_season_pol_start <- fig_seasons$dates[min(which (fig_seasons$in_pol95season == "in 95% season"))]
#fig_season_pol_end   <- fig_seasons$dates[max(which (fig_seasons$in_pol95season == "in 95% season"))]
# fig_season_npn_start <- fig_seasons$dates[min(which (fig_seasons$in_npn_95season == "in 95% season"))]
# fig_season_npn_end   <- fig_seasons$dates[max(which (fig_seasons$in_npn_95season == "in 95% season"))]
fig_season_pol_start <- fig_seasons$dates[min(which (fig_seasons$in_99polseason == "in 99% season"))]
fig_season_pol_end   <- fig_seasons$dates[max(which (fig_seasons$in_99polseason == "in 99% season"))]
fig_season_npn_start <- fig_seasons$dates[min(which (fig_seasons$in_npn_99season == "in 99% season"))]
fig_season_npn_end   <- fig_seasons$dates[max(which (fig_seasons$in_npn_99season == "in 99% season"))]

# correlation
panel_b_cor <- 
  nabnpn %>% 
  filter(site == fig_site) %>%  #unique(nabnpn$site)
  filter(taxon == fig_taxon) %>%  #unique(nabnpn$taxon)
  filter(years == fig_year) %>% 
  #filter(in_pol95season == "in 95% season" & in_npn_95season == "in 95% season") %>% 
  filter(in_99polseason == "in 99% season" & in_npn_99season == "in 99% season") %>% 
  #filter(nobs_yes_per_season > 50) %>% 
  #filter(doy > 95 & doy < 175) %>% 
  summarise(out = cor(polpct, mean_prop_flow_ma, use = "complete.obs", method = "spearman")) 
# ggplot(aes(y = polpct * 100, x = mean_prop_flow_ma * 100)) + geom_point(alpha = 0.7) + 
#   theme_few() + ylab("airborne pollen (% of maximum)") + xlab("flowering (% of observations)")

panel_b <- nabnpn %>% 
  filter(site == fig_site) %>%  #unique(nabnpn$site)
  filter(taxon == fig_taxon) %>%  #unique(nabnpn$taxon)
  filter(years == fig_year) %>% 
  # filter(nobs_yes_per_season > 50) %>% 
  # filter(sum_pol_season > 200) %>% 
  filter(doy > 50 & doy < 170) %>% 
  ggplot(aes(x = dates, y = mean_prop_flow_ma * 100)) + geom_line(col = "blue") + 
  theme_few() + facet_wrap(~years) +
  geom_point(aes(x = dates, y = polpct * 100), alpha = 0.3) + xlab("date") + 
  scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
  theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
  #scale_x_date(date_labels = "%b %d") +
  annotate("text", x= ymd_hms("2016/6/11 00:00:00"), y= 98, label= paste0("\U03C1 = ", round(panel_b_cor[3], 2))) + 
  # 
  # geom_segment(x = fig_season_pol_start, xend = fig_season_pol_end,#pollen 95% season line
  #              y = -2, yend = -2, col = "black", lwd = 2)  +
  # geom_segment(x = fig_season_npn_start, xend = fig_season_npn_end,#npn 95% season line
  #              y = -3, yend = -3, col = "blue", lwd = 2)  +
  #geom_line(aes(y= zoo::rollmean(polpct * 100, 7, na.pad=TRUE)), col = "gray20") 
  geom_line(aes(y= rollapply(polpct * 100, width=7, FUN=function(x) mean(x, na.rm=TRUE), by=1, partial=TRUE, fill=NA)), col = "gray10") 




# panel C: Atlanta Quercus time series 2014
#add another example here
fig_site <- "Springfield"
fig_taxon <- "Acer"
fig_year <- 2018
fig_seasons <- nabnpn %>% filter(site == fig_site & taxon == fig_taxon & years == fig_year) %>%  
  arrange(taxon, site, years, dates)
#fig_season_pol_start <- fig_seasons$dates[min(which (fig_seasons$in_pol95season == "in 95% season"))]
#fig_season_pol_end   <- fig_seasons$dates[max(which (fig_seasons$in_pol95season == "in 95% season"))]
# fig_season_npn_start <- fig_seasons$dates[min(which (fig_seasons$in_npn_95season == "in 95% season"))]
# fig_season_npn_end   <- fig_seasons$dates[max(which (fig_seasons$in_npn_95season == "in 95% season"))]
fig_season_pol_start <- fig_seasons$dates[min(which (fig_seasons$in_99polseason == "in 99% season"))]
fig_season_pol_end   <- fig_seasons$dates[max(which (fig_seasons$in_99polseason == "in 99% season"))]
fig_season_npn_start <- fig_seasons$dates[min(which (fig_seasons$in_npn_99season == "in 99% season"))]
fig_season_npn_end   <- fig_seasons$dates[max(which (fig_seasons$in_npn_99season == "in 99% season"))]

# correlation
panel_c_cor <- 
  nabnpn %>% 
  filter(site == fig_site) %>%  #unique(nabnpn$site)
  filter(taxon == fig_taxon) %>%  #unique(nabnpn$taxon)
  filter(years == fig_year) %>% 
  #filter(in_pol95season == "in 95% season" & in_npn_95season == "in 95% season") %>% 
  filter(in_99polseason == "in 99% season" & in_npn_99season == "in 99% season") %>% 
  #filter(nobs_yes_per_season > 50) %>% 
  #filter(doy > 95 & doy < 175) %>% 
  summarise(out = cor(polpct, mean_prop_flow_ma, use = "complete.obs", method = "spearman")) 
# ggplot(aes(y = polpct * 100, x = mean_prop_flow_ma * 100)) + geom_point(alpha = 0.7) + 
#   theme_few() + ylab("airborne pollen (% of maximum)") + xlab("flowering (% of observations)")

panel_c <- nabnpn %>% 
  filter(site == fig_site) %>%  #unique(nabnpn$site)
  filter(taxon == fig_taxon) %>%  #unique(nabnpn$taxon)
  filter(years == fig_year) %>% 
  # filter(nobs_yes_per_season > 50) %>% 
  # filter(sum_pol_season > 200) %>% 
  filter(doy > 70 & doy < 120) %>% 
  ggplot(aes(x = dates, y = mean_prop_flow_ma * 100)) + geom_line(col = "blue") + 
  theme_few() + #facet_wrap(~years) +
  geom_point(aes(x = dates, y = polpct * 100), alpha = 0.3) + xlab("date") + 
  scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
  theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
  annotate("text", x= ymd_hms("2017/4/26 00:00:00"), y= 98, label= paste0("\U03C1 = ", round(panel_c_cor[3], 2))) +
  #scale_x_date(date_labels = "%b %d") +
  # geom_segment(x = fig_season_pol_start, xend = fig_season_pol_end,#pollen 95% season line
  #              y = -2, yend = -2, col = "black", lwd = 2)  +
  # geom_segment(x = fig_season_npn_start, xend = fig_season_npn_end,#npn 95% season line
  #              y = -3, yend = -3, col = "blue", lwd = 2)  
  geom_line(aes(y= rollapply(polpct * 100, width=7, FUN=function(x) mean(x, na.rm=TRUE), by=1, partial=TRUE, fill=NA)), col = "gray10") 


cowplot::plot_grid(panel_a, panel_b, panel_c, nrow = 3)



### Fig. 3: Acer version and mean annual temperature ###########################

## Panel A: DOY of season peak for a taxon as a function of mean annual temperature
test <- slice_sample(npn, n = 10000)

nabnpn %>% 
  filter(taxon == "Acer") %>% 
  group_by(site, years) %>% 
  slice(cumu_flow_r)


summary(nabnpn$cumu_flow_r)


npn_pol %>% filter(genus == "Acer") %>% 
  filter(phenophase_status == 1) %>% 
  ggplot(aes(x = tmean, y = day_of_year, color = longitude)) + geom_point(alpha = 0.9) + theme_bw() +
  facet_wrap(~species) + scale_color_viridis_c()
  
npn_active_flow %>% filter(genus == "Acer") %>% filter(species == "rubrum" | species == "saccharinum" | species == "saccharum" | 
                                                         species == "negundo") %>% 
  filter(longitude < -60 & longitude > -120) %>% 
  filter(day_of_year < 175) %>% 
  filter(flow_prop > 0) %>% 
  ggplot(aes(x = tmean_jan_apr, y = day_of_year, color = species)) + geom_point(alpha = 0.1) + theme_bw() +
  #scale_color_viridis_c() + 
  #facet_wrap(~species) +
  geom_smooth(method = "lm")

#Acer rubrum regression
npn_active_flow_Acru <- npn_active_flow %>% filter(genus == "Acer") %>% filter(species == "rubrum") %>% 
  filter(longitude < -60 & longitude > -130) %>% 
  filter(day_of_year < 175) %>% 
  filter(flow_prop > 0)
acru_fit <- summary(lm( day_of_year ~ tmean_jan_apr, data = npn_active_flow_Acru))
acru_tmean_predicted_doy <- tmean_rast2 * acru_fit$coefficients[2] + acru_fit$coefficients[1]
raster::plot(acru_tmean_predicted_doy)

#Acer saccharinum regression
npn_active_flow_Acsaccharinum <- npn_active_flow %>% filter(genus == "Acer") %>% filter(species == "saccharinum") %>% 
  filter(longitude < -60 & longitude > -130) %>% 
  filter(day_of_year < 175) %>% 
  filter(flow_prop > 0)
Acsaccharinum_fit <- summary(lm( day_of_year ~ tmean_jan_apr, data = npn_active_flow_Acsaccharinum))
Acsaccharinum_tmean_predicted_doy <- tmean_rast2 * Acsaccharinum_fit$coefficients[2] + Acsaccharinum_fit$coefficients[1]
raster::plot(Acsaccharinum_tmean_predicted_doy)

#Acer negundo regression
npn_active_flow_Acne <- npn_active_flow %>% filter(genus == "Acer") %>% filter(species == "negundo") %>% 
  filter(longitude < -60 & longitude > -130) %>% 
  filter(day_of_year < 175) %>% 
  filter(flow_prop > 0)
acne_fit <- summary(lm( day_of_year ~ tmean_jan_apr, data = npn_active_flow_Acne))
acne_tmean_predicted_doy <- tmean_rast2 * acne_fit$coefficients[2] + acne_fit$coefficients[1]
raster::plot(acne_tmean_predicted_doy)



### Fig 3: time series of observed pollen and time series of observed maple flowering and predicted maple flowering

# empirical flowering curve for Acru
#difference from predicted peak day
npn_flow_empir_Acru <- 
  npn_active_flow %>% filter(genus == "Acer") %>% filter(species == "rubrum") %>% 
  filter(longitude < -60 & longitude > -130) %>% 
  filter(day_of_year < 175) %>% 
  mutate(doy_peak_pred = tmean_jan_apr * acru_fit$coefficients[2] + acru_fit$coefficients[1],
         dif_from_peak = day_of_year - round(doy_peak_pred, 0)) %>% 
  group_by(dif_from_peak) %>% 
  summarize(flow_prop_mean = mean(flow_prop))
site_mean_dif_df <- data.frame(dif_from_peak = -30:30)
npn_flow_empir_Acru <- left_join(npn_flow_empir_Acru, site_mean_dif_df) %>% 
    ungroup() %>% 
    mutate(flow_prop_mean_m = zoo::rollapply(flow_prop_mean, 7, #one week moving average
                                              mean, na.rm = TRUE, partial = TRUE, align='center'))

ggplot(npn_flow_empir_Acru, aes(x = dif_from_peak, y = flow_prop_mean_m)) + geom_line()


# empirical flowering curve for Acne
#difference from predicted peak day
npn_flow_empir_Acne <- 
  npn_active_flow %>% filter(genus == "Acer") %>% filter(species == "negundo") %>% 
  filter(longitude < -60 & longitude > -130) %>% 
  filter(day_of_year < 175) %>% 
  mutate(doy_peak_pred = tmean_jan_apr * acne_fit$coefficients[2] + acne_fit$coefficients[1],
         dif_from_peak = day_of_year - round(doy_peak_pred, 0)) %>% 
  group_by(dif_from_peak) %>% 
  summarize(flow_prop_mean = mean(flow_prop))
site_mean_dif_df <- data.frame(dif_from_peak = -30:30)
npn_flow_empir_Acne <- left_join(npn_flow_empir_Acne, site_mean_dif_df) %>% 
  ungroup() %>% 
  mutate(flow_prop_mean_m = zoo::rollapply(flow_prop_mean, 7, #one week moving average
                                           mean, na.rm = TRUE, partial = TRUE, align='center'))

ggplot(npn_flow_empir_Acne, aes(x = dif_from_peak, y = flow_prop_mean_m)) + geom_line()


# predicted Acru open flower curve at an NAB site
head(nabnpn)
focal_peak_acru <- nab %>% ungroup() %>% filter(site == "Springfield" & years == 2018) %>% 
  select(tmean_jan_apr) %>% distinct() %>% 
  mutate(acru_peak_day_pred = tmean_jan_apr * acru_fit$coefficients[2] + acru_fit$coefficients[1])

npn_flow_empir_acru_focal_site <- npn_flow_empir_Acru %>% 
  mutate(doy = dif_from_peak + round(focal_peak$acru_peak_day_pred[1],0))

# predicted Acne open flower curve at an NAB site
focal_peak_acne <- nab %>% ungroup() %>% filter(site == "Springfield" & years == 2018) %>% 
  select(tmean_jan_apr) %>% distinct() %>% 
  mutate(acne_peak_day_pred = tmean_jan_apr * acne_fit$coefficients[2] + acne_fit$coefficients[1])

npn_flow_empir_acne_focal_site <- npn_flow_empir_Acne %>% 
  mutate(doy = dif_from_peak + round(focal_peak_acne$acne_peak_day_pred[1],0))



acer_ny <- npn_raw %>% 
  mutate(years = year(observation_date)) %>% 
  filter(years == 2016) %>% 
  filter(genus == "Acer") %>% 
  # filter(species == "rubrum" | #species == "saccharum" | #species == "platanoides" | species == "pensylvanicum" |
  #          species == "saccharinum" | 
  #          species == "negundo"
  #        ) %>% 
  filter(NAB_station == "Springfield") %>% 
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
  ggtitle("Maples near Springfield") +
  geom_line(data = npn_flow_empir_acru_focal_site, aes(x = as.Date(doy, origin = as.Date("2018-01-01")),
                                                       y = flow_prop_mean_m
                                                       ), color = "black")+
  geom_line(data = npn_flow_empir_acne_focal_site, aes(x = as.Date(doy, origin = as.Date("2018-01-01")),
                                                       y = flow_prop_mean_m), color = "green")



ggplot(npn_flow_empir_acru_focal_site, aes(x = doy, y = flow_prop_mean_m)) + geom_point()



### Fig 3: Quercus peak doy NPN vs Quercus peak doy NAB #####################################################


xlab <- "spring temperature (°C)"
npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubra" | species == "velutina" | species == "alba" | species == "palustris") %>% 
  filter(longitude < -60 & longitude > -90) %>% 
  filter(day_of_year < 181) %>% 
  filter(flow_prop > 0) %>% 
  ggplot(aes(x = tmean_jan_apr, y = as.Date(day_of_year, origin = as.Date("2018-01-01")))) + geom_point(alpha = 0.1) + ggthemes::theme_few() +
  facet_wrap(~genus) + 
  #scale_color_viridis_c() + 
  #facet_wrap(~species) +
  geom_smooth(method = "lm") + 
  xlab(xlab) + ylab("flowering observed (date)") + scale_y_date(date_breaks = "1 month", date_labels =  "%b") +
  #stat_regline_equation(label.y = as.Date(180, origin = as.Date("2018-01-01")), aes(label = ..eq.label..), label.x = 15) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = as.Date(185, origin = as.Date("2018-01-01")), label.x = 15)

#Quercus regression
npn_active_flow_Qu <- npn_active_flow %>% filter(genus == "Quercus") %>% #filter(species == "rubrum") %>% 
  filter(longitude < -60 & longitude > -90) %>% 
  filter(day_of_year < 181) %>% 
  filter(flow_prop > 0) 
qu_fit <- summary(lm( day_of_year ~ tmean_jan_apr, data = npn_active_flow_Qu))
qu_tmean_predicted_doy <- tmean_rast2 * qu_fit$coefficients[2] + qu_fit$coefficients[1]
raster::plot(qu_tmean_predicted_doy)


# predicted Quru peak open flow day at each NAB site
head(nabnpn)
peak_qu_npn_pred_nab_sites <- nab %>% ungroup() %>% filter(taxon == "Quercus") %>% 
  select(site, taxon, years, tmean_jan_apr) %>% distinct() %>% 
  mutate(quru_peak_day_pred = tmean_jan_apr * qu_fit$coefficients[2] + qu_fit$coefficients[1])

# empirical Quru peak at each NAB site
NAB_quru_peak_day <- nabnpn %>% filter(taxon == "Quercus") %>% #filter(years > 2009) %>% 
  #filter(sum_pol_season > 100) %>% 
  arrange(site, dates) %>% 
  group_by(site, years)%>% 
  mutate(pol_m = round(na_interpolation(pol),1),
         pol_m_ma = round(rollmean(pol_m, 7, na.pad=TRUE),2),
         #n_obs_real = 
         #mean_nobs_ma = round(rollmean(n_obs, 7, na.pad=TRUE),2)
         )%>% 
    
         #filter(mean_nobs_ma > 3) %>% 
         filter(pol_m_ma > 50) %>% 
         slice_max(pol_m_ma)
NAB_NPN_quru_peak_day <- left_join(NAB_quru_peak_day, peak_qu_npn_pred_nab_sites)

ggplot(NAB_NPN_quru_peak_day, aes(x = as.Date(quru_peak_day_pred, origin = as.Date("2018-01-01")), 
                                  y = as.Date(doy, origin = as.Date("2018-01-01")))) + geom_point(aes(color = as.factor(site))) + ggthemes::theme_few() + geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  xlab("predicted peak flowering (date)")+ ylab("observed peak pollen (date)") +
  scale_color_discrete(guide="none") +
  # stat_regline_equation(label.y = as.Date(155, origin = as.Date("2018-01-01")), aes(label = ..eq.label..), 
  #                       label.x = as.Date(70, origin = as.Date("2018-01-01"))) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.y = as.Date(158, origin = as.Date("2018-01-01")), 
           label.x = as.Date(70, origin = as.Date("2018-01-01")))




###########################################################

# ## Panel B: map of US with color of predicted peak DOY
# test <- npn_active_flow %>% filter(genus == "Quercus") %>% 
#   filter(flow_prop > 0) %>% 
#   filter(species == "rubra" | species == "palustris") 
# lm(day_of_year ~ tmean, data = test)
# 
# tmean_rast2
# plot(tmean_rast2)
# raster::plot(tmean_rast2)
# 
# tmean_predicted_doy <- tmean_rast2 * -2.93 + 158.95
# raster::plot(tmean_predicted_doy)
# 
# npn_active_flow %>% filter(genus == "Quercus") %>% 
#   filter(species == "rubra" | species == "palustris") %>% 
#   filter(day_of_year < 150 & day_of_year > 50) %>% 
#   filter(flow_prop > 0) %>% 
#   #filter(phenophase_status == 1) %>%
#   ggplot(aes(x = longitude, y = latitude, color = day_of_year)) + geom_point(size = 3, alpha =0.5) + scale_color_viridis_c() + theme_bw()



## Panel C: predicted peak DOY based on NPN-MAT compared to observed NAB peak
#extract predicted DOY for red oak

head(nab)


### Fig. 3: overall comparisons ###############################
### creating a table of correlations by taxon x site
cor_spear_nobs <- nabnpn %>%  
  #filter(sum_pol_season > 100) %>% 
  #filter(nobs_yes_per_season > 50) %>% 
  filter(in_npn_99season == "in 99% season" & in_99polseason == "in 99% season") %>% 
  filter(!is.na(mean_prop_flow_ma)) %>% 
  filter(!is.na(polpct)) %>% 
  group_by(site, taxon, years) %>% 
  summarize(n_obs_comparison = n())


cor_spear <- nabnpn %>%  
  left_join(., cor_spear_nobs) %>% 
  filter(sum_pol_season > 100) %>% 
  filter(nobs_yes_per_season > 10) %>% 
  filter(n_obs_comparison > 10) %>% 
  #filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
  filter(in_npn_99season == "in 99% season" & in_99polseason == "in 99% season") %>% 
  filter(!is.na(mean_prop_flow_ma)) %>% 
  filter(!is.na(polpct)) %>% 
  group_by(site, taxon, taxon_labs, years) %>% 
  mutate(tmean_dif = tmean - tmean_NAB) %>% 
  summarize(n_obs = n(),
            cor_spear = cor(mean_prop_flow_ma, polpct, method = "spearman", use="complete.obs"),
            cor_p_value = cor.test(mean_prop_flow_ma, polpct, method = "spearman", use="complete.obs")$p.value,
            unique_observers = mean(unique_observers),
            tmean_dif = mean(tmean_dif, na.rm = TRUE),
            tmean = mean(tmean, na.rm = TRUE),
            lat_mean = mean(lat_mean, na.rm = TRUE),
            long_mean = mean(long_mean, na.rm = TRUE),
            distNAB_mean = mean(distNAB_mean, na.rm = TRUE)
  ) %>% 
  mutate(cor_spear = round(cor_spear, 2),
         cor_p_value_discrete = case_when(cor_p_value >= 0.05 ~ "ns",
                                          cor_p_value < 0.05 &
                                            cor_p_value >= 0.01 ~ "p < 0.05",
                                          cor_p_value < 0.01 ~ "p < 0.01")) %>% 
  mutate(cor_p_value_discrete = factor(cor_p_value_discrete, levels = c("ns", "p < 0.05", "p < 0.01"))) %>% 
  arrange(taxon) %>% 
  filter(!is.na(cor_spear)) %>% 
  ungroup()

#str(cor_spear$taxon_labs)

  cor_spear <- cor_spear %>% mutate(taxon_labs2 = forcats::fct_drop(taxon_labs))
  
cor_spear #unique(cor_spear$taxon) 
str(cor_spear$taxon_labs2)

ggplot(cor_spear, aes(x = taxon_labs2, y = cor_spear)) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(group = site,
                                                     color = cor_p_value_discrete), width = 0.1, alpha = .15) + ggthemes::theme_few() +
  ylab("Spearman correlation between airborne pollen and flowering") +
  scale_color_manual(values = c("gray30", "dodgerblue4", "blue4"), name = "") +
  scale_x_discrete(labels = levels(cor_spear$taxon_labs2), name = "taxa") +
  theme(axis.text.x = ggtext::element_markdown())

ggsave(filename = "Fig_2.jpg", width = 20, height = 15, units = "cm", dpi = 300, scale = 1.25)

#some stats for results section
cor_spear %>%  
  summarize(spear_mean = mean(cor_spear, na.rm = TRUE),
            spear_sd = sd(cor_spear, na.rm = TRUE))

cor_spear %>%  
  group_by(taxon) %>% 
  summarize(spear_mean = mean(cor_spear, na.rm = TRUE),
            spear_sd = sd(cor_spear, na.rm = TRUE))


### SI 1: sensitivity analyses for temperature and distance cutoffs ##################################






######## Scraps #######################################################################################

# 
# # #latitude
# # ggplot(cor_spear, aes(x = taxon, y = cor_spear)) + 
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_boxplot(outlier.shape = NA) + geom_jitter(aes(group = site, color = lat_mean), width = 0.1, ) + ggthemes::theme_few() +
# #   ylab("Spearman correlation between airborne pollen and flowering") +
# #   theme(axis.text.x = element_text(face = "italic")) + 
# #   scale_color_viridis_c()
# 
# 
# #latitude
# ggplot(cor_spear, aes(x = as.factor(site), y = cor_spear)) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_jitter(aes(group = as.factor(years)), width = 0.1) + 
#   ggthemes::theme_few() +
#   ylab("Spearman correlation between airborne pollen and flowering") +
#   theme(axis.text.x = element_text(face = "italic")) +
#   facet_wrap(~taxon)
# 
# #geographic
# ggplot(cor_spear, aes(x = long_mean, y = lat_mean, color = cor_spear)) +
#   geom_jitter(aes(group = as.factor(years)), width = 0.01) + 
#   ggthemes::theme_few() +
#   scale_color_viridis_c() +
#   facet_wrap(~taxon)
# 
# 
# 
# length(unique(nabnpn$site))




# # creating a similar table without the temperature restriction 
# # (i.e., manually re-running script after changing filt_tmean_dif)
# cor_spear_notemp <- nabnpn %>%  
#   filter(sum_pol > 200) %>% 
#   filter(nobs_yes_per_season > 50) %>% 
#   filter(in_npn_95season == "in 95% season" & in_pol95season == "in 95% season") %>% 
#   #filter(site == "Armonk" & taxon == "Quercus") %>% 
#   filter(!is.na(mean_prop_flow_ma)) %>% 
#   filter(!is.na(polpct)) %>% 
#   group_by(site, taxon, years) %>% 
#   summarize(cor_spear = cor(mean_prop_flow_ma, polpct, method = "spearman"),
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
  filter(NAB_station == "Springfield") %>% 
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
  ggplot(aes(x = mean_prop_flow_ma * 100, y = polpct * 100)) + 
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
  ggplot(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = mean_prop_flow_ma * 100)) + geom_line(col = "blue") + 
  theme_few() + facet_grid(site~years) +
  geom_point(aes(x = as.Date(doy, origin = as.Date("2018-01-01")), y = polpct * 100)) + xlab("date") + 
  scale_y_continuous(name="flowering (% of observations)", sec.axis=sec_axis(~., name="airborne pollen (% of maximum)")) +
  theme(axis.title.y.left=element_text(color="blue"), axis.text.y.left=element_text(color="blue")) +
  scale_x_date(date_labels = "%b %d")


