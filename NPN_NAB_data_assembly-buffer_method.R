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
#
# **This script uses Theresa's method of including each NN obs within a specified buffer aroudn a NAB station to that station
#
### set up working environment #############################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(rnpn)
library(sf)
library(geosphere)
library(prism)

### select a few top anemophilous taxa from NPN #####################################
#anemophilous angiosperms
acer_species_list <- c(777,1843,59,778,1,2,1591,60,779,780,3,781,61,1199)
alnus_species_list <- c(62,63,319)
betula_species_list <- c(97, 1439, 98, 430, 1850, 1339, 1851, 99, 1805)
fraxinus_species_list <- c(74,872,873,75,1350)
populus_species_list <- c(1361,320,976,977,1188,27,1481)
quercus_species_list <- c(705,100,1365,757,1870,987,1690,1484,988,316,297,1485,1190,765,1486,
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
  request_source = 'Daniel Katz, UT Austin and/or Theresa Crimmins',
  species_ids = list_all_focal_taxa,
  years = c(as.character(2008:2021)), #years to include
  phenophase_ids = c(501, 502,495, 503), #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" #conifers: 495 ==  503 ==
  additional_fields = c("Observed_Status_Conflict_Flag", "partner_group")
)
#names(npn_direct)

npn_direct <- subset(npn_direct, observed_status_conflict_flag == "-9999")

npn_direct <- filter(npn_direct, !(partner_group %in% c("CSU Chico NSCI 102", "SCMF Naturalists", 
                                                        "Sycamore Canyon Wilderness Park - Riverside",
                                                        "Marist College", "Sam Hughes Neighborhood Association",
                                                        "UNCO BIO 111", "Maricopa Cooperative Extension",
                                                        "Pima County Extension", "Lasell College",
                                                        "UofL campus", "Ursinus College", "U of A Campus Arboretum",
                                                        "RMC Campus Phenology", "SUNY Geneseo", "AZ Project WET")))
                                      
npn_flow <- filter(npn_direct, phenophase_id == 501 | phenophase_id == 495)

# looks like flowering intensity value was only entered very rarely - TMC: actually it's around 82%
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
filter(npn_direct, phenophase_id == 502 | phenophase_id == 503) %>% 
  group_by(phenophase_status) %>% 
  dplyr::summarize(n_obs = n())

### extract mean annual air temperature for each NPN observation site ##################################
npn_active_flow_sf <- npn_active_flow %>% 
  dplyr::select(longitude, latitude, site_id) %>% 
  distinct() %>% 
  st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 


prism_set_dl_dir("~/RProjects/DanK_analyses/prism")
get_prism_normals(type = "tmean", resolution = "4km", annual = TRUE)

#choosing which rasters need to be included
tmean_rast <- prism_archive_subset(type = "tmean", temp_period = "annual normals", resolution = "4km")
tmean_rast2 <- pd_stack(tmean_rast)

#extract the data from prism raster(s) #note: allowing the raster package to handle the CRS transformation
tmean_data <- unlist(raster::extract(x = tmean_rast2, 
                                     y = npn_active_flow_sf 
)) %>% as.data.frame() 

npn_active_flow3 <- npn_active_flow_sf %>%  
  dplyr::select(site_id) %>% 
  mutate(tmean = unlist(tmean_data)) 

npn_active_flow3$geometry <- NULL

npn_active_flow <- left_join(npn_active_flow, npn_active_flow3)

### extract mean annual air temperature for each selected NAB site #####################################

#extract temperature for the NAB stations [and add it to active flowers dataframe - NOT SURE ABOUT THIS STEP]
########## IF WE CAN GET MORE NAB STATIONS - NEED TO ADD THEIR LAT-LONG INFO HERE ################
NAB_coords <- data.frame(station = c("Armonk", "Waterbury", "Flower Mound", "Minneapolis", "Springfield", "Carrolton"),
                         lat = c(41.1299814, 41.5493, 33.0439926,  44.9749718, 40.7002184, 33.0439971),	
                         long = c(-73.7310037, -73.068176, -97.0780927, -93.2756491, -74.3244188, -96.8341116)) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) 

tmean_data_NAB <- unlist(raster::extract(x = tmean_rast2, 
                                         y = NAB_coords )) %>% as.data.frame() 

NAB_coords_tmean <- NAB_coords %>% mutate(tmean_NAB = unlist(tmean_data_NAB)) 

NAB_coords_tmean$geometry <- NULL

### calculate geographic distance from each NPN site to each NAB site #########
########## IF WE CAN GET MORE NAB STATIONS - NEED TO ADD THEIR STUFF INFO HERE ################

# calculate distance from each NAB stn and each NN obs station one at a time, filter by buffer, and append

npn_active_flow$distNAB <- distm(as.matrix(npn_active_flow[,c(6,5)]), c(-73.7310037, 41.1299814), fun = distHaversine) #Armonk
npn_Armonk <- npn_active_flow[(npn_active_flow$distNAB <= 321869),] #200 miles
npn_Armonk$NABStn <- "Armonk"

npn_active_flow <- subset(npn_active_flow, select= -distNAB)
npn_active_flow$distNAB <- distm(as.matrix(npn_active_flow[,c(6,5)]), c(-97.0780927, 33.0439926), fun = distHaversine) #FlowerMound
npn_FlowerMnd <- npn_active_flow[(npn_active_flow$distNAB <= 321869),] #200 miles
npn_FlowerMnd$NABStn <- "FlowerMound"

npn_buffer_ok <- rbind(npn_Armonk, npn_FlowerMnd) #append FlowerMnd to Armonk in new dataframe
rm(npn_Armonk)
rm(npn_FlowerMnd)

npn_active_flow <- subset(npn_active_flow, select= -distNAB)
npn_active_flow$distNAB <- distm(as.matrix(npn_active_flow[,c(6,5)]), c(-93.2756491, 44.9749718), fun = distHaversine) #Minneapolis
npn_Minneapolis <- npn_active_flow[(npn_active_flow$distNAB <= 321869),] #200 miles
npn_Minneapolis$NABStn <- "Minneapolis"

npn_buffer_ok <- rbind(npn_buffer_ok, npn_Minneapolis) #append Minneapolis to nn_buffer_ok rm(npn_Armonk)
rm(npn_Minneapolis)

npn_active_flow <- subset(npn_active_flow, select= -distNAB)
npn_active_flow$distNAB <- distm(as.matrix(npn_active_flow[,c(6,5)]), c(-74.3244188, 40.7002184), fun = distHaversine) #Springfield
npn_Springfield <- npn_active_flow[(npn_active_flow$distNAB <= 321869),] #200 miles
npn_Springfield$NABStn <- "Springfield"

npn_buffer_ok <- rbind(npn_buffer_ok, npn_Springfield) #append Springfield
rm(npn_Springfield)

npn_active_flow <- subset(npn_active_flow, select= -distNAB)
npn_active_flow$distNAB <- distm(as.matrix(npn_active_flow[,c(6,5)]), c(-73.068176, 41.5493), fun = distHaversine) #Waterbury
npn_Waterbury <- npn_active_flow[(npn_active_flow$distNAB <= 321869),] #200 miles
npn_Waterbury$NABStn <- "Waterbury"

npn_buffer_ok <- rbind(npn_buffer_ok, npn_Waterbury) #append Waterbury
rm(npn_Waterbury)

npn_active_flow <- subset(npn_active_flow, select= -distNAB)
npn_active_flow$distNAB <- distm(as.matrix(npn_active_flow[,c(6,5)]), c(-96.8341116, 33.0439971), fun = distHaversine) #Carrolton
npn_Carrolton <- npn_active_flow[(npn_active_flow$distNAB <= 321869),] #200 miles
npn_Carrolton$NABStn <- "Carrolton"

npn_buffer_ok <- rbind(npn_buffer_ok, npn_Carrolton) #append Carrolton
rm(npn_Carrolton)

##### ADD IN TMEAN VALUE FOR NAB STATIONS, CALCULATE tmean_dif

npn_buffer_ok <- left_join(npn_buffer_ok, NAB_coords_tmean, by = c("NABStn" = "station"))
npn_buffer_ok$tmean_dif <- (npn_buffer_ok$tmean_NAB - npn_buffer_ok$tmean)


### export data to file (data exploration is next script) ##################################
readr::write_csv(npn_active_flow, "data/200mibuffer-inclusive.csv")

