### NPN-NAB comparison manuscript  ####################################################################
# Current project participants:
# Liz Vogt, Theresa Crimmins, Arie Managan, Claudia Brown, Dan Dalan, and Dan Katz
# This script includes data assembly, analysis, and visualization and is the compilation of several 
# previous versions, also stored in this repo

# Note: DK is  working on assembling and cleaning up this script (April 27; check back in a day or two for a prettier version)

### set up working environment #############################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(rnpn)
library(sf)
library(geosphere)
library(prism)
library(readr)

### select a few top anemophilous taxa from NPN #####################################

#top woody anemophilous angiosperms
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
pinus_species_list <- c(1629,1686,965,762,50,295,220,1480,219,51,966,967,52,25,968,1687,53,54) #Pinaceae
cupressaceae_species_list <- c(43,1743,1354,289,902,291,290,44, #junipers
                               1723,831,1356,296,1040) #other

list_all_focal_taxa <- c(acer_species_list, alnus_species_list, betula_species_list, fraxinus_species_list, 
                         populus_species_list, quercus_species_list, ulmus_species_list,
                         ambrosia_species_list, pinus_species_list, cupressaceae_species_list)

###download and process data #####################################
npn_direct <- npn_download_status_data(
  request_source = 'Daniel Katz, Cornell and/or Theresa Crimmins',
  species_ids = list_all_focal_taxa,
  years = c(as.character(2009:2021)), #years to include
  phenophase_ids = c(501, 502,495, 503), #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" #conifers: 495 ==  503 ==
  additional_fields = c("Observed_Status_Conflict_Flag", "partner_group")
)

#write_csv(npn_direct, "C:/Users/danka/Box/Cornell/1 national pollen model/NAB_NPN/NPN_220308.csv")
#npn_direct <- read_csv("C:/Users/danka/Box/Cornell/1 national pollen model/NAB_NPN/NPN_220308.csv") 
#names(npn_direct)

#Theresa - is this QA/QC section depreciated? If so, let's remove these next 8 lines
# npn_direct <- subset(npn_direct, observed_status_conflict_flag == "-9999")
# npn_direct <- filter(npn_direct, !(partner_group %in% c("CSU Chico NSCI 102", "SCMF Naturalists",
#                                                        "Sycamore Canyon Wilderness Park - Riverside",
#                                                        "Marist College", "Sam Hughes Neighborhood Association",
#                                                        "UNCO BIO 111", "Maricopa Cooperative Extension",
#                                                        "Pima County Extension", "Lasell College",
#                                                       "UofL campus", "Ursinus College", "U of A Campus Arboretum",
#                                                        "RMC Campus Phenology", "SUNY Geneseo", "AZ Project WET")))

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
#prism_set_dl_dir("C:/Users/dsk856/Documents/prismtmp")
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
  mutate(tmean = as.numeric(unlist(tmean_data)))

npn_active_flow3$geometry <- NULL

npn_active_flow <- left_join(npn_active_flow, npn_active_flow3)

### extract mean annual air temperature for each selected NAB site #####################################
#NAB data were assembled in this script: #C:/Users/danka/Box/texas/NAB/extract_pollen_data_from_NPNdata220308.R
NAB <- read_csv("C:/Users/danka/Box/Cornell/1 national pollen model/NAB_NPN/NAB_data_request_220308e.csv")

#extract temperature for the NAB stations and add it to active flowers dataframe
NAB_coords <- NAB %>% dplyr::select(NAB_station, Lat, Long) %>% 
  distinct() %>% 
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) 

tmean_data_NAB <- unlist(raster::extract(x = tmean_rast2, #matrix(c(NAB_tx_coords$long, NAB_tx_coords$lat), ncol = 2), 
                                         y = NAB_coords )) %>% as.data.frame() 
NAB_coords_tmean <- NAB_coords %>% mutate(tmean_NAB = as.numeric(unlist(tmean_data_NAB))) 
NAB_coords_tmean$geometry <- NULL


### calculate geographic distance from each NPN site to each NAB site ################################
#NAB_stations <- unique(NAB_coords_tmean$NAB_station)
NAB_coords_notsf <- NAB %>% dplyr::select(NAB_station, NAB_station_lat = Lat, NAB_station_long = Long) %>% distinct() 
#NAB_coords_notsf <- NAB_coords_notsf[1:2,] #for testing with a subset of stations 
npn_active_flow_coords_only <- npn_active_flow %>% dplyr::select(longitude, latitude) %>% as.matrix(.)
npn_coords_all <- as.matrix(npn_active_flow[,c(6,5)])

#function for extracting 
dist_calc_fun <- function(NAB_station, NAB_station_long, NAB_station_lat){
  #NAB_coords_row <- NAB_coords_notsf[NAB_station_row,] #NAB_coords_row <- NAB_coords_notsf[1,]
  #NAB_coords_focal <- c(NAB_station_long, NAB_station_lat) #NAB_coords_row[,3:2] #make sure to have it in the order of long, lat
  npn_active_flow_dist_i <- distm( x = npn_active_flow_coords_only, 
                                   y = c(NAB_station_long, NAB_station_lat), 
                                   fun = distHaversine)
  
  npn_active_flow_i <- npn_active_flow %>% mutate(distNAB = as.numeric(npn_active_flow_dist_i),
                                                  NAB_station = NAB_station)
  npn_active_flow_i_200m <- filter(npn_active_flow_i, distNAB <= 321869) #filter observations within 200 miles
  return(npn_active_flow_i_200m)
}

Sys.time()
NPN_near_NAB <- pmap_dfr(NAB_coords_notsf, dist_calc_fun)
Sys.time()

#ggplot(NPN_near_NAB, aes(x = distNAB)) + geom_histogram() + facet_wrap(~NAB_station)
#ggplot(NPN_near_NAB, aes(x = longitude, y = latitude)) + geom_point() 


##### ADD IN TMEAN VALUE FOR NAB STATIONS, CALCULATE tmean_dif ################

NPN_near_NAB2 <- left_join(NPN_near_NAB, NAB_coords_tmean) %>% 
  mutate(tmean_dif = as.numeric(tmean_NAB - tmean))

str(NPN_near_NAB2)
str(NAB_coords_tmean)



### export data to file (data exploration is next script) ##################################
write_csv(NPN_near_NAB2, "C:/Users/danka/Box/Cornell/1 national pollen model/NAB_NPN/npn_200mi_2c_220321.csv")

# readr::write_csv(NPN_near_NAB2, "data/200mibuffer-inclusive_220321.csv")
# write.csv(NPN_near_NAB2, "data/200mibuffer-inclusive_220321b.csv")



#double check that everything added up well
NPN_near_NAB2 %>% 
  group_by(NAB_station) %>% 
  filter(genus == "Quercus") %>%   
  summarize(n = n())

