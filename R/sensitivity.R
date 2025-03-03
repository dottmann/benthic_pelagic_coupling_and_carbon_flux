
## Script name: Sensitivity
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: January 2024
## Last update:  March 2025
##
## ---------------------------
##
## Readme:
##
## This script tuns a sensitivity analysis to multiple FEISTY parameters
##
## ---------------------------


#####################################################################
# Load libraries:
library(tidyverse)
library(ggthemes)
library(viridis)
library(sf)
library(marmap)
library(mapdata)
library(geosphere)

# Reset directory:
# Get the directory path of the currently running script and set directory:
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path.expand(script_dir))
setwd('../')

# Clear environment:
rm(list = ls())

# Create "exclude" function:
"%ni%" <- Negate("%in%")

# Load functions:
source("R/pathway_functions.R")


#-------------------------------------------------------------------------------------------------
# Run FEISTY ----------------------> this can be run in parallel, or save as a separate file

# Load selected COBALT data:
load(file = "data/cobalt_NAtlantic_slopes.Rdata")

# Filter data to remove NAs:  
cobalt_NAtlantic_slopes <- cobalt_NAtlantic_slopes %>%
  filter(!is.na(photic), !is.na(lz_prod), !is.na(mz_prod))

time0 <- Sys.time()
# my_sims <- run_feisty_sensitivity(cobalt_NAtlantic_slopes)
Sys.time() -time0

# Outfile data:
# save(my_sims, file = "data/simulations_NAtlantic_slopes_sensitivity_am0.9.Rdata")

# Load simulations:
# Basic run:
load(file = "data/simulations_NAtlantic_slopes.Rdata")

benchmark <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$Biomass, by = "id") %>%
  filter(region == "Western Ireland", funGroup == "demersals", depth >500 & depth <1800) %>%
  summarise(mean = mean(totBiomass))

# Sensitivity runs:
load(file = "data/simulations_NAtlantic_slopes_sensitivity.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_Cmax0.9.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_Cmax1.1.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_dfbot0.9.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_dfbot1.1.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_dfbot2.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_epsAssim0.9.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_epsAssim1.1.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_metabolism0.9.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_metabolism1.1.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_bc1.1.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_bc0.9.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_bm1.1.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_bm0.9.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_am1.1.Rdata")

load(file = "data/simulations_NAtlantic_slopes_sensitivity_am0.9.Rdata")



my_sims$Biomass <- my_sims$Biomass %>%
  filter(funGroup %ni% c("Zoo", "benthos"))

#------------
my_params <- setupVertical2()
my_palette <- my_params$my_palette


#-----------------------------------------------
# Edit data:
data <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$Biomass, by = "id") %>%
  filter(region != "Others")


# Calculate Biomass:
sensitivity_test <- data %>%
  filter(region == "Western Ireland", funGroup == "demersals", depth >500 & depth <1800) %>%
  summarise(mean = mean(totBiomass))

100 * (sensitivity_test - benchmark) / benchmark

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

# Load trawling data of the Celtic sea:
load(file = "data/celtic_trawl_data.Rdata")


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
  xlab("Depth (m)") +
  xlim(c(NA, 700)) +
  theme_base() +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = c(.8, .8),
        legend.title = element_blank()) 

p2


p
ggsave("plots/totBiomass_western_ireland.png", p2, height = 36 , width = 40, units = "mm", scale = 3)


#     END OF SCRIPT
############################