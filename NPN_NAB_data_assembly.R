### General information ####################################################################
# Project goals: 
# -assess correlations between NPN and NAB data 
# -provide proof of concept for future continental-scale airborne pollen models
#
# current/potential project participants:
# Dan Katz, Theresa Crimmins, Liz Vogt, Arie Managan, Claudia Brown, Ellen Denny, 
# Shubhayu Saha, Ambarish (Rish) Vaidyanathan
#
# This script is for assembling relevant datasets and filtering NPN observations
# Note: One of the goals is to keep this as simple as possible;
# full phenological models are a future project


### set up working environment #############################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(rnpn)
library(sf)
library(prism)

### select a few top anemophilous taxa from NPN #####################################
#anemophilous angiosperms
acer_species_list <- c(777,1843,59,778,1,2,1591,60,779,780,3,781,61,1199)
alnus_species_list <- c(62,63,319)
betula_species_list <- c(97, 1439, 98, 430, 1850, 1339, 1851, 99, 1805)
fraxinus_species_list <- c(74,872,873,75,1350)
populus_species_list <- c(1361,320,976,977,1188,27,1481)
quercus_species_list <-c(705,100,1365,757,1870,987,1690,1484,988,316,297,1485,1190,765,1486,
                         301,704,101,1691,1212,989,1366,102,1756,1213,1755,1487,1159,305)
ulmus_species_list <- c(1192,1048,1049,1215,1216)

#herbaceous angiosperms
ambrosia_species_list <- c(145,788,146)

#pollen cones
pinus_species_list <- c(1629,1686,965,762,50,295,220,1480,219,51,966,967,52,25,968,1687,53,54)
cupressaceae_species_list <- c(43,1743,1354,289,902,291,290,44, #junipers
                               1723,831,1356,296,1040) #other

list_all_focal_taxa <- c(acer_species_list, alnus_species_list, betula_species_list, fraxinus_species_list, 
                         populus_species_list, quercus_species_list, ulmus_species_list,
                         ambrosia_species_list, pinus_species_list, cupressaceae_species_list)

###download and process data #####################################
npn_direct <- npn_download_status_data(
  request_source = 'Daniel Katz, UT Austin',
  species_ids = list_all_focal_taxa,
  years = c(as.character(2008:2021)), #years to include
  phenophase_ids = c(501, 502,495, 503) #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" #conifers: 495 ==  503 ==
)
#names(npn_direct)

npn_flow <- filter(npn_direct, phenophase_id == 501 | phenophase_id == 495)

# looks like flowering intensity value was only entered very rarely
npn_active_flow <- npn_flow %>%
  filter(phenophase_status != -1) %>% #removing observations where the observer was unsure whether the phenophase was occurring
  mutate(flow_prop = case_when(
    phenophase_status == 0  ~ 0,
    intensity_value == "Less than 5%" ~ 0.025,
    intensity_value ==  "5-24%"~ (0.05+0.24)/2,
    intensity_value == "25-49%" ~ (0.25+0.49)/2,
    intensity_value == "50-74%" ~ (0.5+0.74)/2,
    intensity_value == "75-94%" ~ (0.75+0.94)/2,
    intensity_value == "95% or more" ~ 0.97,
    phenophase_status == 1  ~ 0.5)) #assuming that when intensity value isn't given, a tree is halfway through flowering
hist(npn_active_flow$flow_prop)

#how many observations are for pollen release?                    
filter(npn_direct, phenophase_id == 502 | phenophase_id == 503) %>% group_by(phenophase_status) %>% summarize(n_obs = n())

### calculate geographic distance from each NPN site to the nearest focal NAB site #########
NAB_coords <- data.frame(station = c("Armonk", "Atlanta", "Waterbury"),
                         lat = c(41.1299814, 33.974,	41.5493),	
                         long = c(-73.7310037, -84.5493, -73.0659)) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) 

npn_active_flow_sf <- npn_active_flow %>% 
  select(longitude, latitude, site_id) %>% 
  distinct() %>% 
  st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 


distances <- st_distance(npn_active_flow_sf, NAB_coords, by_element = FALSE) /1000 #calculate distances and convert to km
distances_df <- as.data.frame(distances) 
distances_min <- apply(distances_df, 1, FUN = min) #minimum distance to a NAB station
which_station_closest <- apply(distances_df, 1, function(x) which(x == min(x, na.rm = TRUE))) #which station is closest
NAB_station_lookup <- data.frame(NAB_station = NAB_coords$station, n_lookup = 1:length(unique(NAB_coords$station)))
station_looked_up <- left_join(data.frame(n_lookup = which_station_closest), NAB_station_lookup)
npn_active_flow_sf <- mutate(npn_active_flow_sf, NAB_min_dist = distances_min, NAB_station = station_looked_up$NAB_station)

npn_active_flow2 <- npn_active_flow_sf 
npn_active_flow2$geometry <- NULL
npn_active_flow <- left_join(npn_active_flow, npn_active_flow2)

### extract mean annual air temperature for each NPN observation site ##################################
prism_set_dl_dir("C:/Users/dsk856/Documents/prismtmp")
get_prism_normals(type = "tmean", resolution = "4km", annual = TRUE)

#choosing which rasters need to be included
tmean_rast <- prism_archive_subset(type = "tmean", temp_period = "annual normals", resolution = "4km")
tmean_rast2 <- pd_stack(tmean_rast)

#extract the data from prism raster(s) #note: allowing the raster package to handle the CRS transformation
tmean_data <- unlist(raster::extract(x = tmean_rast2, #matrix(c(NAB_tx_coords$long, NAB_tx_coords$lat), ncol = 2), 
                                     y = npn_active_flow_sf #, buffer = 10
)) %>% as.data.frame() 

npn_active_flow3 <- npn_active_flow_sf %>%  
  dplyr::select(site_id) %>% 
  mutate(tmean = unlist(tmean_data)) 
npn_active_flow3$geometry <- NULL
npn_active_flow <- left_join(npn_active_flow, npn_active_flow3)



### extract mean annual air temperature for each selected NAB site #####################################
# selected NAB sites: Atlanta, Armonk, Waterbury
# note: if we wanted to do this for Cupressaceae, that's covered by my most recent NAB data
# request

#extract temperature for the NAB stations and add it to active flowers dataframe
tmean_data_NAB <- unlist(raster::extract(x = tmean_rast2, #matrix(c(NAB_tx_coords$long, NAB_tx_coords$lat), ncol = 2), 
                                         y = NAB_coords )) %>% as.data.frame() 
NAB_coords_tmean <- NAB_coords %>% mutate(tmean_NAB = unlist(tmean_data_NAB)) %>% 
  rename(NAB_station = station)
NAB_coords_tmean$geometry <- NULL

npn_active_flow <- left_join(npn_active_flow, NAB_coords_tmean) %>% 
  mutate(tmean_dif = tmean_NAB - tmean)



### export data to file (data exploration is next script) ##################################
readr::write_csv(npn_active_flow, "C:/Users/dsk856/Box/texas/pheno/npn_common_anemo_active_flow.csv")

