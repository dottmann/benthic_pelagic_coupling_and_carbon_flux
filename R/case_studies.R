
## Script name: Case studies
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: February 2024
## Last update:  July 2024
##
## ---------------------------
##
## Readme:
##
## This script analyzes each of the 5 case-study areas
##
## ---------------------------


#####################################################################
# Load libraries:
library(tidyverse)
library(ggthemes)
library(geosphere)
library(FEISTY)
library(marmap)
library(mapdata)
library(patchwork)

# Reset directory:
# Get the directory path of the currently running script and set directory:
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path.expand(script_dir))
setwd('../')

# Clear environment:
rm(list = ls())

# Create "exclude" function:
"%ni%" <- Negate("%in%")

# Load selected COBALT points:
load(file = "data/cobalt_NAtlantic_slopes.Rdata")

cobalt_NAtlantic_slopes <- cobalt_NAtlantic_slopes %>% 
  filter(!is.na(photic), !is.na(lz_prod), !is.na(mz_prod))

# Load FEISTY simulations:
load(file = "data/simulations_NAtlantic_slopes.Rdata")

my_sims$Biomass <- my_sims$Biomass %>%
  filter(funGroup %ni% c("Zoo", "benthos"))

# Load trawling data of the Celtic sea:
load(file = "data/celtic_trawl_data.Rdata")

# Load trawling data of the Mid Atlantic Bight:
load(file = "data/mid_trawl_data.Rdata")


#-----------------------------------------------
# Edit data:
data <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$Biomass, by = "id") %>%
  filter(region != "Others")

#-------------------------------------------------
# Get color palette:
my_params <- setupVertical2()
my_palette <- my_params$my_palette
region_pallet <- c("#6699CC", "#117733", "#AA4499", "#332288", "grey85", "#A65B00")

#----------------------------------------------------------------------------
# Plot relative biomas:
df <- data %>%
  mutate(depth_round = round(depth / 100) * 100) %>%
  group_by(depth_round, coast, region, funGroup) %>%
  reframe(mean_totBiomass = mean(totBiomass),
            mean_propBiomass = mean(propBiomass),
            sd_totBiomass = sd(totBiomass),
            se_totBiomass = sd_totBiomass/sqrt(n()),
            sd_propBiomass = sd(propBiomass),
            se_propBiomass = propBiomass/sqrt(n()))

# Plot diet as a function of 4 flux pathways:
data2 <- cobalt_NAtlantic_slopes %>%
  filter(region != "Others") %>%
  mutate(coast = "East") %>%
  left_join(my_sims$diet, by = "id") %>%
  group_by(id, region, coast, depth, flux_pathway) %>%
  summarise(ingestion = sum(ingestion)) %>%
  mutate(rel_ingestion = ingestion / sum(ingestion))

df2 <- data2 %>%
  mutate(depth_round = round(depth / 100) * 100) %>%
  group_by(depth_round, region, coast, flux_pathway) %>%
  reframe(mean_ingestion = mean(ingestion),
          sd_ingestion = sd(ingestion),
          mean_rel_ingestion = mean(rel_ingestion),
          sd_rel_ingestion = sd(rel_ingestion))


data_full2 <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$diet, by = "id") %>%
  mutate(region = "a) NE & NW Atlantic") %>%
  group_by(id, region, coast, depth, flux_pathway) %>%
  summarise(ingestion = sum(ingestion)) %>%
  mutate(rel_ingestion = ingestion / sum(ingestion))

df_full2 <- data_full2 %>%
  mutate(depth_round = round(depth / 100) * 100) %>%
  group_by(depth_round, region, coast, flux_pathway) %>%
  reframe(mean_ingestion = mean(ingestion),
          sd_ingestion = sd(ingestion),
          mean_rel_ingestion = mean(rel_ingestion),
          sd_rel_ingestion = sd(rel_ingestion))

df_full2 <- rbind(df_full2, df2)

df_full2 <- df_full2 %>%
   mutate(region = case_when(region == "Western Ireland" ~ "b) Western Ireland",
                             region == "Gulf of Mexico" ~ "c) Gulf of Mexico",
                             region == "Mid Atlantic Bight" ~ "d) Mid Atlantic Bight",
                             region == "North Caribbean" ~ "e) North Caribbean",
                             region == "Sahelian upwelling" ~ "f) Sahelian upwelling",
                             T ~ region))

p <- ggplot() +
  geom_line(data = df_full2, 
            aes(x = depth_round, y = 100 * mean_rel_ingestion, color = flux_pathway, linetype = coast)) +
  geom_ribbon(data = df_full2, 
              aes(x = depth_round, ymin = 100 * ifelse((mean_rel_ingestion - sd_rel_ingestion) < 0, 0, (mean_rel_ingestion - sd_rel_ingestion)),
                  ymax = 100 * (mean_rel_ingestion + sd_rel_ingestion), fill = flux_pathway,
                  linetype = coast), 
              alpha = 0.1) +
  scale_fill_manual(values = c("#A65B00", as.vector(my_palette)[c(8, 6, 5)])) +
  scale_color_manual(values = c("#A65B00", as.vector(my_palette)[c(8, 6, 5)])) +
  geom_hline(yintercept = 0, alpha = .3) +
  xlab("Bottom depth (m)") +
  ylab(expression(paste("% flux of C to demersal fishes"))) +
  xlim(c(NA, 2500)) +
  facet_wrap(vars(region)) +
  theme_base() +
  theme(plot.background = element_blank(),
        legend.title = element_blank())

p
ggsave("plots/prop_flux_to_demersals_cases_all.png", p, height = 60 , width = 80, units = "mm", scale = 3)


#-----------------------------------------------------------------------------------------------------
# Compare detrital flux vs active carbon transport:
# Calculate total active flux:
data4 <- cobalt_NAtlantic_slopes %>%
  filter(
    # depth >500 & depth <1800,
    region != "Others") %>%
  left_join(my_sims$diet, by = "id") %>%
  mutate(flux_pathway2 = case_when(flux_pathway %in% c("Midwater", "Epipelagic") ~ "Pelagic",
                                   flux_pathway == "Demersal" ~  "Demersal",
                                   T ~ flux_pathway),
         detrital_flux = dflux / 9) %>% # Convert to C equivalents: 1 g C = 9 g WW
  group_by(id, long, lat, region, depth, detrital_flux, flux_pathway2) %>%
  summarise(ingestion = sum(ingestion) / 9) %>% # Convert to C equivalents: 1 g C = 9 g WW
  ungroup()

temp <- data4 %>%
  filter(flux_pathway2 != "Demersal") %>%
  group_by(id) %>% 
  mutate(sum_ingestion = sum(ingestion),# If we remove demesal, how much of the remaining ingestion is benthic vs (midwater + pelagic)?
         rel_ingestion = ingestion / sum_ingestion) %>%
  filter(flux_pathway2 != "Benthic") %>% # We want to keep the ingested proportion of pelagic food
  dplyr::select(id, rel_ingestion) %>%
  right_join(data4, by = "id") %>%
  filter(flux_pathway2 == "Demersal") %>% # To calculate how much of the demersal flux comes originated in the pelagic system, we multiply it by the proportion of ingested pelagic food
  mutate(pelagic_from_demersal = rel_ingestion * ingestion) %>%
  dplyr::select(id, pelagic_from_demersal) %>%
  right_join(data4, by = "id") %>%
  filter(flux_pathway2 == "Pelagic") %>%
  mutate(flux = ingestion + pelagic_from_demersal, # Now we add the pelagic flux + the demersal flux that was originated in the pelagic system
         flux_pathway = "Active") %>%
  ungroup() %>%
  dplyr::select(id, long, lat,  region, depth, flux, flux_pathway)


data4 <- data4 %>%
  filter((row_number() - 1) %% 3 == 0) %>%
  dplyr::select(id, long, lat, region, depth, detrital_flux) %>%
  rename(flux = detrital_flux) %>%
  mutate(flux_pathway = "Detritus") %>%
  rbind(temp) %>%
  mutate(flux = case_when(flux < 0 ~ 0,
                   T ~ flux))

# Calculate proportion of active flux:
data4 %>%
  filter(depth >= 500, depth <= 1800) %>%
  group_by(region, flux_pathway) %>%
  summarise(flux = sum(flux)) %>%
  pivot_wider(names_from = flux_pathway, values_from = flux) %>%
  mutate(flux_ratio = 100 * Active / Detritus)
  
  
df4 <- data4 %>%
  filter(flux_pathway ==  "Active") 

grid_points <- df4 %>%
  dplyr::select(lat, long) %>%
  rename(latitude = lat,
         longitude = long)


# Function to calculate the corners of each grid cell
get_corners <- function(lat, lon) {
  half_res <- 0.25 / 2
  lat_min <- lat - half_res
  lat_max <- lat + half_res
  lon_min <- lon - half_res
  lon_max <- lon + half_res
  return(matrix(c(lon_min, lat_min, 
                  lon_min, lat_max, 
                  lon_max, lat_max, 
                  lon_max, lat_min, 
                  lon_min, lat_min), 
                ncol = 2, byrow = TRUE))
}

# Create a list to store the corners of each grid cell
corners_list <- mapply(get_corners, df4$lat, df4$long, SIMPLIFY = FALSE)

# Calculate the area for each grid cell in m^2
grid_areas <- sapply(corners_list, function(corners) {
  areaPolygon(corners)
})

# Add the area to the original data frame
df4 <- df4 %>%
  mutate(area = grid_areas,
         c_area = area * flux) # m^2 * g m^-2 = g

# Calculate total active C and CO2 flux:
df4 %>%
  filter(depth >= 500, depth <= 1800) %>%
  group_by(region) %>%
  summarise(area_km2 = sum(area) / 1e6, # m^-2 to km^2
    c_region = sum(c_area) / 1e6,   # g to tonnes of C year^-1
    # co2_region = (44/12) * c_region, # Million tonnes of CO2 year^-1
    c_region_km2 = c_region * 1000 / area_km2) # tones km^-2 to kg km^-2
    # CO2_km2 = co2_region * 1e6 / area_km2) 


#------------------------------------------------------------------------------------------------
# Special case of Celtic sea:
load("data/ICESsurveys.RData")

# Get all countries and transform to the same CRS
cntries <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")

# Get bathymetry data:
b = getNOAA.bathy(lon1 = -16, lon2 = -8, lat1 = 51, lat2 = 55, 
                  resolution = 1)
# Querying NOAA database ...
# This may take seconds to minutes, depending on grid size

# Convert bathymetry to data frame:
bf = fortify.bathy(b)

# Get regional polygons:
reg = map_data("world2Hires")

# Convert lat longs
reg$long = (360 - reg$long)*-1

lons = c(-16, -8)
lats = c(51, 55)

survey <- survey3 %>%
  filter(Depth <= 750, 
         ShootLong > -16, ShootLong < -8, 
         ShootLat > 51, ShootLat < 55, 
         class != "Cephalopoda") 

p1 <- survey %>%
  ggplot() +
  
  # add bathymetry contour
  geom_contour(data = bf, 
               aes(x = x, y = y, z = z),
               # breaks = c(-1000),
               breaks = seq(from = -200, to = -2600, by = -200),
               size = c(.3),
               colour = c("grey90")) +
  
  ylab("Longitude") + xlab("Latitude") +
  
  geom_sf(data = cntries, fill = "grey", colour = NA) +
  
  
  geom_point(aes(x = ShootLong, y = ShootLat, color = Depth))+
  scale_color_continuous(trans = 'reverse') +
  labs(color = "Bottom depth (m)") +
  
  # configure projection and plot domain
  coord_sf(xlim = lons, ylim = lats) +
  theme_base() +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = c(.26, .7)) +
  ggtitle("a")


p1


data_celtic_feisty <- data %>%
  filter(region == "Western Ireland") %>%
  dplyr::select(depth, funGroup, totBiomass) %>%
  mutate(source = "FEISTY")

df_celtic_trawl <- df_celtic_trawl %>%
  ungroup(HaulID) %>%
  dplyr::select(-HaulID) %>%
  rename(depth = Depth,
         totBiomass = wgtlencpue_q) %>%
  mutate(funGroup = case_when(funGroup == "small_pelagic" ~ "smallPel",
                              funGroup == "mesopelagic" ~ "mesoPel",
                              funGroup == "large_pelagic" ~ "largePel",
                              funGroup == "bathypelagic" ~ "bathyPel",
                              T ~ "demersals"),
         source = "trawl")

df_celtic_trawl %>%
  mutate(depth_round = round(depth / 100) * 100) %>%
  group_by(depth_round) %>%
  summarise(count = n())

df_celtic <- rbind(data_celtic_feisty, df_celtic_trawl)

df_celtic_2 <- df_celtic %>%
  mutate(depth_round = round(depth / 100) * 100) %>%
  group_by(depth_round, funGroup, source) %>%
  reframe(sd_totBiomass = sd(totBiomass),
          se_totBiomass = sd_totBiomass/sqrt(n()),
          totBiomass = mean(totBiomass)) %>%
  filter(funGroup == "demersals")


p2 <- ggplot() +
  geom_line(data = df_celtic_2,
             aes(x = depth_round, y = totBiomass, color = source)) +
  geom_ribbon(data = df_celtic_2, 
              aes(x = depth_round, ymin = ifelse((totBiomass - sd_totBiomass) < 0, 0, (totBiomass - sd_totBiomass)), ymax = totBiomass + sd_totBiomass, fill = source), 
              alpha = 0.2) +
  scale_fill_manual(values = c("#FF00FF", "gray50")) +
  scale_color_manual(values = c("#FF00FF", "gray50")) +
  geom_hline(yintercept = 0, alpha = .3) +
  ylab(expression(paste("Biomass (g ww ",m^2,")"))) +  
  xlab("Bottom depth (m)") +
  xlim(c(NA, 700)) +
  theme_base() +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = c(.8, .8),
        legend.title = element_blank()) +
  ggtitle("b")

p2


p <- p1 + p2  &
  theme(plot.background = element_blank()) &
  theme(plot.tag = element_text(size = 10))

p
ggsave("plots/totBiomass_celtic_map.png", p, height = 36 , width = 80, units = "mm", scale = 3)


# Mid Atlantic Bight:
data_mid_feisty <- data %>%
  filter(region == "Mid Atlantic Bight") %>%
  dplyr::select(depth, funGroup, totBiomass) %>%
  mutate(source = "FEISTY") %>%
  filter(funGroup == "demersals")


df_mid_trawl %>%
  mutate(depth_round = round(depth / 100) * 100) %>%
  group_by(depth_round) %>%
  summarise(count = n()) 

data_mid_trawl_unique <- df_mid_trawl %>%
  group_by(haulid) %>%
  slice(1)

# Set map limits:
lons = c(-78, -65)
lats = c(34, 42)

p3 <- data_mid_trawl_unique %>%
  ggplot() +
  
  # add bathymetry contour
  geom_contour(data = bf, 
               aes(x = x, y = y, z = z),
               # breaks = c(-1000),
               breaks = seq(from = -200, to = -2600, by = -200),
               size = c(.3),
               colour = c("grey90")) +
  
  ylab("Longitude") + xlab("Latitude") +
  
  geom_sf(data = cntries, fill = "grey", colour = NA) +
  
  
  geom_point(aes(x = lon, y = lat, color = depth)) +
  scale_color_continuous(trans = 'reverse', limits = c(750, 0)) +
  labs(color = "Depth (m)") +
  
  # configure projection and plot domain
  coord_sf(xlim = lons, ylim = lats) +
  theme_base() +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "none") +
  ggtitle("c")

p3


df_mid_trawl <- df_mid_trawl %>%
  ungroup() %>%
  dplyr::select(depth, funGroup, totBiomass, source)


df_mid <- rbind(data_mid_feisty, df_mid_trawl)

df_mid_2 <- df_mid %>%
  mutate(depth_round = round(depth / 100) * 100) %>%
  group_by(depth_round, funGroup, source) %>%
  reframe(sd_totBiomass = sd(totBiomass),
          se_totBiomass = sd_totBiomass/sqrt(n()),
          totBiomass = mean(totBiomass)) 


p4 <- ggplot() +
  geom_line(data = df_mid_2,
            aes(x = depth_round, y = totBiomass, color = source)) +
  geom_ribbon(data = df_mid_2, 
              aes(x = depth_round, ymin = ifelse((totBiomass - sd_totBiomass) < 0, 0, (totBiomass - sd_totBiomass)), ymax = totBiomass + sd_totBiomass, fill = source), 
              alpha = 0.2) +
  scale_fill_manual(values = c("#FF00FF", "gray50")) +
  scale_color_manual(values = c("#FF00FF", "gray50")) +
  geom_hline(yintercept = 0, alpha = .3) +
  ylab(expression(paste("Biomass (g ww ",m^2,")"))) +  
  xlab("Bottom depth (m)") +
  xlim(c(NA, 700)) +
  theme_base() +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  ggtitle("d") 

p4


p <- p3 + p4  &
  theme(plot.background = element_blank()) &
  theme(plot.tag = element_text(size = 10))

p
ggsave("plots/totBiomass_mid_map.png", p, height = 36 , width = 80, units = "mm", scale = 3)


p <- p1 + p2 + p3 + p4 +
  plot_layout(ncol = 2) &
  theme(plot.background = element_blank()) &
  theme(plot.tag = element_text(size = 10))

p
ggsave("plots/totBiomass_map.png", p, height = 70 , width = 80, units = "mm", scale = 3)


#-------------------------------------------------

# creates feeding flux from prey to predator
getFeeding = function(sim) {
  p <- sim$p
  u <- sim$u
  
  # get last 40% of timeseries
  etaTime <- 0.4 
  ixTime  <- which(sim$t>=((1-etaTime)*sim$t[sim$nTime]))
  
  # Number of groups and biomass:
  ngroup <- p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])] #resources (4) + fish 
  biomass <- u
  
  #Average of the biomass : 
  prey     <- colMeans(biomass[ixTime,]) # mean value of the last 40% time 
  predator <- colMeans(biomass[ixTime,]) # mean value of the last 40% time 
  
  # estimate encounters
  Enc <- p$V[1] * p$theta[1,] * prey 
  for (i in 2:length(predator)){
    Enc <- rbind(Enc, p$V[i] * p$theta[i,] * prey)
  }
  
  # estimate encounters per predator
  Encspecies = rowSums(Enc)
  
  # estimate the mortality generated
  mortpr <-  p$Cmax[1] * p$V[1] * p$theta[1,] / (p$Cmax[1]+ Encspecies[1]) * predator[1]
  for (i in 2:length(predator)){
    mortpr <- rbind(mortpr,p$Cmax[i] * p$V[i] * p$theta[i,] / (p$Cmax[i]+ Encspecies[i]) * predator[i])
  }
  
  # estimate the flux from prey to predator
  mortpr <- ifelse(is.na(mortpr),0,mortpr)
  flux_prey_to_pred <- t(t(mortpr)*prey)
  rownames(flux_prey_to_pred) <- paste("pred",colnames(p$theta),sep="_")
  colnames(flux_prey_to_pred) <- paste("prey",colnames(p$theta),sep="_")
  return(flux_prey_to_pred)
}


# Load FEISTY simulations:
load(file = "data/simulations_NAtlantic_slopes.Rdata")

# Load selected COBALT points:
load(file = "data/cobalt_NAtlantic_slopes.Rdata")

cobalt_NAtlantic_slopes <- cobalt_NAtlantic_slopes %>% 
  filter(!is.na(photic), !is.na(lz_prod), !is.na(mz_prod))

# Edit data:
data <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$Biomass, by = "id") %>%
  filter(region != "Others")

temp <- cobalt_NAtlantic_slopes %>% 
  filter(!is.na(photic), !is.na(lz_prod), !is.na(mz_prod), region != "Others") %>%
  group_by(region) %>%
  summarise(mz_prod = mean(mz_prod),
            lz_prod = mean(lz_prod),
            dflux = mean(dflux),
            depth = mean(depth),
            photic = mean(photic),
            Temp_C = mean(Temp_C),
            tob = mean(tob),
            Temp_m_C = mean(Temp_m_C, na.rm = T))

out <- data.frame()
out2 <- data.frame()
region <- temp$region
SpId <- c("Zooplankton", "Benthos", "Small pelagics", "Large pelagics", "Demersals", "Midwater fish")

my_setup <-  for(i in 1:nrow(temp)) { # 
  
  # Customize FEISTY setup:
  my_setup <- setupVertical2(szprod = temp[i, 2], lzprod = temp[i, 3], # Pelagic productivities
                             dfbot = temp[i, 4], # detrital flux reaching the bottom
                             nStages = 10, # No. of size groups
                             Tp = as.numeric(temp[i, 7]), # Temperature of the surface
                             Tb = as.numeric(temp[i, 8]), # Temperature of the bottom
                             Tm = as.numeric(temp[i, 9]), # Temperature of the midwater
                             depth = as.numeric(temp[i, 5]), # Bottom depth
                             photic = as.numeric(temp[i, 6])) # photic depth
  
  
  sim <- simulateFEISTY(bCust    = T,
                        p      = my_setup, 
                        tEnd   = 200,
                        USEdll = TRUE) 
  
  
  resources <- sim$R %>%
    as.data.frame() %>%
    tail(80) %>%
    summarise(across(everything(), mean)) %>%
    mutate(Zoo = smallZoo + largeZoo) %>%
    dplyr::select(Zoo, benthos) %>%
    t() %>%
    as.data.frame() %>%
    rename(totBiomass = V1) %>%
    rownames_to_column(var = "funGroup") %>%
    mutate(region = region[i])
  
  totBiomass <- sim$totBiomass %>%
    as.data.frame() %>%
    tail(80) %>%
    summarise(across(everything(), mean)) %>%
    rename(smallPel = totBiomass.smallPel,
           mesoPel = totBiomass.mesoPel,
           largePel = totBiomass.largePel,
           midwPred = totBiomass.midwPred,
           demersals = totBiomass.demersals) %>%
    mutate(midwater = mesoPel + midwPred) %>%
    dplyr::select(-mesoPel, -midwPred) %>%
    t() %>%
    as.data.frame() %>%
    rename(totBiomass = V1) %>%
    rownames_to_column(var = "funGroup") %>%
    # mutate(propBiomass = totBiomass / sum(totBiomass)) %>%
    mutate(region = region[i])
  
  totBiomass <- rbind(resources, totBiomass)
  
  out <- rbind(out, totBiomass)
  
  # Create line width: 
  Flux <- getFeeding(sim)
  # Flux <- c(Flux)
  
  # Identify columns that contain the strings for each functional group:
  zoo_cols <- grep("Zoo", colnames(Flux))
  benthos_cols <- grep("benthos", colnames(Flux))
  smallPel_cols <- grep("smallPel", colnames(Flux))
  mesoPel_cols <- grep("mesoPel", colnames(Flux))
  largePel_cols <- grep("largePel", colnames(Flux))
  midwPred_cols <- grep("midwPred", colnames(Flux))
  demersals_cols <- grep("demersals", colnames(Flux))
  midwater_cols <- c(mesoPel_cols, midwPred_cols)
  
  # Sum the columns from each functional group:
  prey_zoo <- rowSums(Flux[, zoo_cols])
  prey_benthos <- Flux[, benthos_cols]
  prey_smallPel <- rowSums(Flux[, smallPel_cols])
  prey_mesoPel <- rowSums(Flux[, mesoPel_cols])
  prey_largePel <- rowSums(Flux[, largePel_cols])
  prey_midwPred <- rowSums(Flux[, midwPred_cols])
  prey_demersals <- rowSums(Flux[, demersals_cols])
  prey_midwater <- rowSums(Flux[, midwater_cols])
  
  # Create a new matrix with the summed column:
  new_Flux <- cbind(prey_zoo, prey_benthos, prey_smallPel, prey_largePel, prey_demersals, prey_midwater)
  
  # Convert negative values to 0
  new_Flux <- pmax(new_Flux, 0 )
  
  # Identify rows that contain the strings for each functional group:
  zoo_rows <- grep("Zoo", rownames(new_Flux))
  benthos_rows <- grep("benthos", rownames(new_Flux))
  smallPel_rows <- grep("smallPel", rownames(new_Flux))
  mesoPel_rows <- grep("mesoPel", rownames(new_Flux))
  largePel_rows <- grep("largePel", rownames(new_Flux))
  midwPred_rows <- grep("midwPred", rownames(new_Flux))
  demersals_rows <- grep("demersals", rownames(new_Flux))
  midwater_rows <- c(mesoPel_rows, midwPred_rows)
  
  
  # Sum the rows from each functional group:
  pred_zoo <- colSums(new_Flux[zoo_rows, ])
  pred_benthos <- new_Flux[benthos_rows, ]
  pred_smallPel <- colSums(new_Flux[smallPel_rows, ])
  pred_mesoPel <- colSums(new_Flux[mesoPel_rows, ])
  pred_largePel <- colSums(new_Flux[largePel_rows, ])
  pred_midwPred <- colSums(new_Flux[midwPred_rows, ])
  pred_demersals <- colSums(new_Flux[demersals_rows, ])
  pred_midwater <- colSums(new_Flux[midwater_rows, ])
  

  
  # Create a new matrix with the summed rows:
  Flux <- rbind(pred_zoo, pred_benthos, pred_smallPel, pred_largePel, pred_demersals, pred_midwater)
  Flux[3, 1:6] <- 0 # We don't care about small pelagics feeding on others
  Flux[1:6, 5] <- 0 # We don't care about who feeds on demersals
  diag(Flux) <- 0
  Flux[6, 4] <- 0 # We don't care about midwaters feeding on large pelagics
  Flux <- c(t(Flux)) 
  
  # Marker size depends on biomass following a cubic transformation
  # Get fish biomass predicted by FEISTY:
  fgroup_order <- c("Zoo", "benthos", "smallPel", "largePel", "demersals", "midwater")
  
  xx <- data %>% 
    mutate(funGroup = case_when(funGroup %in% c("mesoPel", "midwPred") ~ "midwater",
                                T ~ funGroup)) %>%
    group_by(id, region, funGroup) %>%
    reframe(totBiomass = sum(totBiomass)) %>%
    ungroup() %>%
    group_by(region, funGroup) %>%
    reframe(mean_totBiomass = mean(totBiomass)) %>%
    mutate(funGroup = factor(funGroup, levels = fgroup_order)) %>%
    arrange(region, funGroup)
  
 
  
  Msize <- subset(xx, region == temp$region[i])$mean_totBiomass / max(xx$mean_totBiomass) #max(totBiomass$totBiomass)
  Msize[Msize == 0] <- NA
  Msize <- Msize^(1/2)  
  
  # Set values of each coordinate (i.e. size and water column position) and put together:
    coord_1 <- data.frame(index = 1:length(Flux),
                        mc = rep(c(0, 0, 1, 2, 2.5, 2), length(Msize)), 
                        depth = rep(c(1, 0, 1, 1, 0, .5), length(Msize)), 
                        SpId = rep(SpId, length(Msize)),
                        Msize = Msize, 
                        LineWdth = (Flux/93)^(1/3), 
                        Alpha = (Flux/93)^(1/3)) 
  
  
  coord_2 <- data.frame(index = 1:length(Flux), # Notice repetition grouped by "each" to change order
                        mc = rep(c(0, 0, 1, 2, 2.5, 2), each = length(Msize)), 
                        depth = rep(c(1, 0, 1, 1, 0, .5), each = length(Msize)), 
                        SpId = rep(SpId, each = length(Msize)),
                        Msize = rep(Msize, each = length(Msize)),
                        LineWdth = (Flux/93)^(1/3), 
                        Alpha = (Flux/93)^(1/3)) 
  
  # Combine in a data frame:
  df <- rbind(coord_1, coord_2)
  
  df <- df %>%
    filter(LineWdth > 0)
  
  df <- df[order(df$index),]
  
  df <- df %>%
    mutate(region = region[i])
  
  out2 <- rbind(out2, df)
} 

  
plotdf <- subset(out2, region == "Western Ireland")
msize <- max(plotdf$Msize)

p <- ggplot() +
  geom_line(data = plotdf, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
            show.legend = F, linewidth = subset(out2, region == "Western Ireland")$LineWdth * 1) +
  geom_point(data = subset(out2, region == "Western Ireland"), aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
  scale_color_manual(values = c("#A65B00", "#228833", "#33BBEE", "#AA4499", "#EE6677", "#999933")) +
  scale_radius(limits = c(0, NA), range = c(0, 10 * msize)) +
  guides(size = "none") +
  theme_base() +
  theme(legend.position = "none", 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = region_pallet[1], size = 1.3),
        panel.background = element_rect(fill = "transparent")) +
  lims(x = c(-.25, 2.75), y = c(-.25, 1.25))

p
ggsave("plots/network_Western_Ireland.png", p, height = 23 , width = 30, units = "mm", scale = .9, dpi = 1000)


plotdf <- subset(out2, region == "North Caribbean")
msize <- max(plotdf$Msize)

p <- ggplot() +
  geom_line(data = plotdf, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
            show.legend = F, linewidth = subset(out2, region == "North Caribbean")$LineWdth * 1) +
  geom_point(data = subset(out2, region == "North Caribbean"), aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
  scale_color_manual(values = c("#A65B00", "#228833", "#33BBEE", "#AA4499", "#EE6677", "#999933")) +
  scale_radius(limits = c(0, NA), range = c(0, 10 * msize)) +
  guides(size = "none")+
  theme_base() +
  theme(legend.position = "none", 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),        
        panel.border = element_rect(color = region_pallet[4], size = 1.3),
        panel.background = element_rect(fill = "transparent")) +
  lims(x = c(-.25, 2.75), y = c(-.25, 1.25))

p
ggsave("plots/network_caribbean.png", p, height = 23 , width = 30, units = "mm", scale = .9, dpi = 1000)


plotdf <- subset(out2, region == "Gulf of Mexico")
msize <- max(plotdf$Msize)

p <- ggplot() +
  geom_line(data = plotdf, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
            show.legend = F, linewidth = subset(out2, region == "Gulf of Mexico")$LineWdth * 1) +
  geom_point(data = subset(out2, region == "Gulf of Mexico"), aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
  scale_color_manual(values = c("#A65B00", "#228833", "#33BBEE", "#AA4499", "#EE6677", "#999933")) +
  scale_radius(limits = c(0, NA), range = c(0, 10 * msize)) +
  guides(size = "none")+
  theme_base() +
  theme(legend.position = "none", 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = region_pallet[2], size = 1.3),
        panel.background = element_rect(fill = "transparent")) +
  lims(x = c(-.25, 2.75), y = c(-.25, 1.25))

p
ggsave("plots/network_gom.png", p, height = 23 , width = 30, units = "mm", scale = .9, dpi = 1000)


plotdf <- subset(out2, region == "Mid Atlantic Bight")
msize <- max(plotdf$Msize)

p <- ggplot() +
  geom_line(data = plotdf, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
            show.legend = F, linewidth = subset(out2, region == "Mid Atlantic Bight")$LineWdth * 1) +
  geom_point(data = subset(out2, region == "Mid Atlantic Bight"), aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
  scale_color_manual(values = c("#A65B00", "#228833", "#33BBEE", "#AA4499", "#EE6677", "#999933")) +
  scale_radius(limits = c(0, NA), range = c(0, 10 * msize)) +
  guides(size = "none")+
  theme_base() +
  theme(legend.position = "none", 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = region_pallet[3], size = 1.3),
        panel.background = element_rect(fill = "transparent")) +
  lims(x = c(-.25, 2.75), y = c(-.25, 1.25))

p
ggsave("plots/network_mabight.png", p, height = 23 , width = 30, units = "mm", scale = .9, dpi = 1000)



plotdf <- subset(out2, region == "Sahelian upwelling")
msize <- max(plotdf$Msize)

p <- ggplot() +
  geom_line(data = plotdf, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
            show.legend = F, linewidth = subset(out2, region == "Sahelian upwelling")$LineWdth * 1) +
  geom_point(data = subset(out2, region == "Sahelian upwelling"), aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
  scale_color_manual(values = c("#A65B00", "#228833", "#33BBEE", "#AA4499", "#EE6677", "#999933")) +
  scale_radius(limits = c(0, NA), range = c(0, 10 * msize)) +
  guides(size = "none")+
  theme_base() +
  theme(legend.position = "none", 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = region_pallet[6], size = 1.3),
        panel.background = element_rect(fill = "transparent")) +
  lims(x = c(-.25, 2.75), y = c(-.25, 1.25))

p
ggsave("plots/network_sahelian.png", p, height = 23 , width = 30, units = "mm", scale = .9, dpi = 2000)

p <- ggplot() +
  geom_line(data = out2, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
            show.legend = F, linewidth = out2$LineWdth*3) +
  geom_point(data = out2, aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
  scale_color_manual(values = c("#A65B00", "#228833", "#33BBEE", "#AA4499", "#EE6677", "#999933")) +
  scale_radius(limits = c(0, NA), range = c(0, 20)) +
  guides(size = "none") +
  theme_base() +
  theme(legend.position = "none", 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank()) +
  lims(x = c(-.15, 3.15), y = c(-.15, 1.15)) +
  facet_wrap(~region)

p
# ggsave("plots/network_all.png", p, height = 50 , width = 80, units = "mm", scale = 1.8)



# Get reference values:

Msize <- c(50, 5, .5, .05) / max(xx$mean_totBiomass)
Flux <- c(50, 10, 1, 0)

Msize <- Msize^(1/2) 

# Set values of each coordinate (i.e. size and water column position) and put together:
coord_1 <- data.frame(index = 1:4,
                      mc = c(0.1, 1.4, 2.3, 3), 
                      depth = c(0, 0, 0, 0), 
                      SpId = c("A", "B", "C", "D"),
                      Msize = Msize, 
                      # LineWdth = (Flux/93)^(1/3), 
                      Alpha = 1) 


coord_2 <- data.frame(index = 1:4, # Notice repetition grouped by "each" to change order
                      mc = c(1, 2, 3, 3), 
                      depth = c(0, 0, 0, 0), 
                      SpId = c("A", "B", "C", "D"),
                      Msize = NA,
                      # LineWdth = (Flux/93)^(1/3), 
                      Alpha = 1) 

# Combine in a data frame:
refdf <- rbind(coord_1, coord_2)

p <- ggplot() +
  # geom_line(data = refdf, aes(x = mc, y = depth, group = index),
  #           show.legend = F, linewidth = refdf$LineWdth) +
  geom_point(data = refdf, aes(x = mc, y = depth, size = Msize), stroke = 0, shape = 16, color = "gray60") +
  scale_radius(limits = c(0, NA), range = c(0, 10 * max(Msize))) +
  guides(size = "none") +
  theme_base() +
  theme(legend.position = "none", 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(size = 1.3),
        panel.background = element_rect(fill = "transparent")) +
  lims(x = c(-.6, 3.2)) 

p
ggsave("plots/network_reference.png", p, height = 23 , width = 30, units = "mm", scale = .9, dpi = 2000)


#              END OF SCRIPT