
## Script name: Run FEISTY in N Atlantic slopes
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: January 2024
## Last update:  March 2025 2024
##
## ---------------------------
##
## Readme:
##
## This script loads and manipulates COBALT data
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

# Load functions:
source("R/pathway_functions.R")

# Create "exclude" function:
"%ni%" <- Negate("%in%")

# Load forcing data from COBALT:
load(file = "data/forcings_1991_1995_15_averaged.Rdata")

data_cobalt <- as.data.frame(forcings) %>%
  rename(mz_prod = medium_zoo, 
         lz_prod = large_zoo)


#-----------------------------------------------
# Edit data COBALT:

data_cobalt <- data_cobalt %>%
  mutate(id = 1:nrow(data_cobalt))


#---------------------
# Define areas of study:

# West coast:
# Build a polygon:
ply = data.frame(
  lon = c( -16, -10, -4, 2.5, 5, 14,
           20, 12, 8, 5, -6,-9, -9, -6, 0, -1, -8, -8, -5, -13.5, -16, -16,
           -19, -19, -18, -9.5, -10.5, -10.5, -10.5, -4, -13, -16, -16),
  lat = c(54, 59, 62, 65, 68, 70,
          70, 68, 65, 62, 59, 54, 50, 48, 45, 43, 43, 38, 36, 25, 20, 10,
          10, 20, 25, 36, 38, 40, 44.5, 44.5, 48,
          51, 54)
)

# Check polygon:
ggplot() +
  geom_polygon(data = ply, aes(x = lon, y = lat), fill = "lightblue", color = "black") +
  coord_fixed(ratio = 1) +  # Adjust the aspect ratio if needed
  theme_minimal()


# Convert to sf geometry list:
ply_sf <- st_sf(st_sfc(geometry = st_polygon(list(as.matrix(ply)))))
st_crs(ply_sf) = 4326


# Get coordinates of cobalt data:
coords <- data_cobalt %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# Overlay shape with COBALT coordinates:
overlap <- st_intersects(ply_sf, coords)

overlap <- overlap %>% 
  as.data.frame() %>%
  rename(COBALT_id = col.id)


# Filter COBALT data for only point swithin the polyong:
grid_match <- data_cobalt[unique(overlap$COBALT_id), ]

# Check
ggplot() +
  geom_point(data = grid_match, aes(x = long, y = lat))

# Ready cobalt data:
cobalt_NEAtlantic_slopes <- grid_match %>%
  dplyr::select(c("id", "lat", "long", "mz_prod", "lz_prod", "dflux", "depth", "photic", "Temp_C", "tob", "Temp_m_C")) %>%
  mutate(coast = "East")

ggplot() +
  geom_point(data = cobalt_NEAtlantic_slopes,
             aes(x = depth, y = dflux, color = lat)) +
  theme_base()
 
# Outfile data:
save(cobalt_NEAtlantic_slopes, file = "data/cobalt_NEAtlantic_slopes.Rdata")


#-------------------------------------
# West coast:

# Build polygon 1 for the west coast:
ply1 = data.frame(
  lon = c(-40, 
          -45, -55, -62, -62, -56, -55, -74, -81, -81, -76, -71, -65, -65,
          -56, -60, -60, -67, -75, -76, -72, -55, -51, -41, -52, -57, -55, -45, -40,
          -40),
  lat = c(62, 
          62, 65, 65, 56, 53, 47, 40, 30, 25, 20, 18, 17, 10,
          10, 17, 19, 20, 25, 30, 37, 44, 41, 47, 56, 60, 61, 58, 58,
          62)
)

# Check polygon:
ggplot() +
  geom_polygon(data = ply1, aes(x = lon, y = lat), fill = "lightblue", color = "black") +
  coord_fixed(ratio = 1) +  # Adjust the aspect ratio if needed
  theme_minimal()


# Convert to sf geometry list:
ply_sf1 <- st_sfc(geometry = st_polygon(list(as.matrix(ply1))))
st_crs(ply_sf1) = 4326


# Polygon 2 for the Gulf of Mexico:
# Repeat for the GOM:
ply2 = data.frame(
  lon = c(-98, 
          -85,-82,-89,
          -92, -98,
          -98),
  lat = c(28, 
          31, 23, 21, 
          18, 18,
          28)
)

# Check polygon:
ggplot() +
  geom_polygon(data = ply2, aes(x = lon, y = lat), fill = "lightblue", color = "black") +
  coord_fixed(ratio = 1) +  # Adjust the aspect ratio if needed
  theme_minimal()

# Convert to sf geometry list:
ply_sf2 <- st_sfc(geometry = st_polygon(list(as.matrix(ply2))))
st_crs(ply_sf2) = 4326

combined_poly <- st_union(ply_sf1, ply_sf2)

# Get coordinates of cobalt data:
coords <- data_cobalt %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# Overlay shape with COBALT coordinates:
overlap <- st_intersects(combined_poly, coords)

overlap <- overlap %>% 
  as.data.frame() %>%
  rename(COBALT_id = col.id)


# Filter COBALT data for only point swithin the polyong:
grid_match <- data_cobalt[unique(overlap$COBALT_id), ]

# Check:
ggplot() +
  geom_point(data = grid_match, aes(x = long, y = lat))


# Ready cobalt data:
cobalt_NWAtlantic_slopes <- grid_match %>%
  dplyr::select(c("id", "lat", "long", "mz_prod", "lz_prod", "dflux", "depth", "photic", "Temp_C", "tob", "Temp_m_C")) %>%
  mutate(coast = "West")

ggplot() +
  geom_point(data = cobalt_NWAtlantic_slopes, 
             aes(x = depth, y = dflux, color = lat)) +
  theme_base()


# Outfile data:
save(cobalt_NWAtlantic_slopes, file = "data/cobalt_NWAtlantic_slopes.Rdata")


#------------------------------------------------------------------------------
# Merge East and West files in one:
cobalt_NAtlantic_slopes <- cobalt_NEAtlantic_slopes %>%
  rbind(cobalt_NWAtlantic_slopes) %>%
  mutate(region = case_when(long < -82 ~ "Gulf of Mexico",
                            lat >51 & lat <55 & long >-16 ~ "Western Ireland",
                            lat < 20 & long > -25 ~ "Sahelian upwelling",
                            lat >35 & lat <41 & long <(-30) ~ "Mid Atlantic Bight",
                            lat <22 & lat >15 & long <(-30) ~ "North Caribbean",
                            T ~ "Others"))


# Outfile data:
save(cobalt_NAtlantic_slopes, file = "data/cobalt_NAtlantic_slopes.Rdata")

# Using marmap and mapdata
# Get bathymetry data:
b = getNOAA.bathy(lon1 = -100, lon2 = 20, lat1 = 10, lat2 = 70, 
                  resolution = 2)
# Querying NOAA database ...
# This may take seconds to minutes, depending on grid size

# Convert bathymetry to data frame:
bf = fortify.bathy(b)

# Get regional polygons:
reg = map_data("world2Hires")

# Convert lat longs
reg$long = (360 - reg$long)*-1

# Set map limits
lons = c(-94, 11)
lats = c(12, 68)

# Get all countries and transform to the same CRS
cntries <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")

region_pallet <- c("#117733", "#AA4499", "#332288", "grey85", "#661100", "#6699CC")

p <- ggplot() +
  
  geom_sf(data = cntries, fill = "grey55", colour = NA) +

  # COBALT points:
  geom_point(data = cobalt_NAtlantic_slopes, 
             aes(x = long, y = lat, color = region), size = .05) +
  scale_color_manual(values = region_pallet) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  
  # configure projection and plot domain
  coord_sf(xlim = lons, ylim = lats) +
  
  # formatting
  ylab("Lat") + xlab("Lon") +
  theme_base() +
  theme(plot.background = element_blank(),
        legend.position = c(.52, .33),
        legend.title = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(),
        legend.key.width = unit(0, "pt"))

p
ggsave("plots/cobalt_points_N_Atlantic.png", p, height = 50 , width = 80, units = "mm", scale = 2)



#-------------------------------------------------------------------------------------------------
# Run FEISTY ----------------------> this can be run in parallel, or save as a separate file

# Load selected COBALT data:
load(file = "data/cobalt_NAtlantic_slopes.Rdata")

# Filter data to remove NAs:  
cobalt_NAtlantic_slopes <- cobalt_NAtlantic_slopes %>%
  filter(!is.na(photic), !is.na(lz_prod), !is.na(mz_prod))

time0 <- Sys.time()
# my_sims <- run_feisty(cobalt_NAtlantic_slopes)
Sys.time() -time0

# Outfile data:
# save(my_sims, file = "data/simulations_NAtlantic_slopes.Rdata")

# Load simulations:
load(file = "data/simulations_NAtlantic_slopes.Rdata")

# # Correct label of pelagic fish:
my_sims$Biomass <- my_sims$Biomass %>%
  filter(funGroup %ni% c("Zoo", "benthos"))

#------------
my_params <- setupVertical2()
my_palette <- my_params$my_palette


# Plot biomasses:
df1 <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$Biomass, by = "id") %>%
  mutate(propBiomass = case_when(funGroup %in% c("midwPred", "mesoPel") & propBiomass == 0 & depth <= 250 ~ NA,
                                T ~ propBiomass)) %>%
  mutate(funGroup = ifelse(funGroup %in% c("midwPred", "mesoPel"), "midwater", funGroup)) %>%
  group_by(id, lat, long, depth, coast, funGroup) %>%
  summarize(totBiomass = sum(totBiomass),
            propBiomass = sum(propBiomass)) %>%
  ungroup()

df1_2 <- df1 %>%
    mutate(propBiomass = case_when(is.na(propBiomass) ~ 0,
                                   T ~ 100 * propBiomass),
         totBiomass = case_when(is.na(totBiomass) ~ 0,
                                 T ~ totBiomass)) 

# Plot it in a map:
lons = c(-94, 11)
lats = c(12, 68)

df1_2 <- df1_2 %>%
  mutate(funGroup = case_when(funGroup == "demersals" ~ "a)",
                              funGroup == "smallPel" ~ "b)",
                              funGroup == "largePel" ~ "c)",
                              funGroup == "midwater" ~ "d)",
                              T ~ funGroup))

p <- ggplot() +
  
  geom_sf(data = cntries, fill = "grey", colour = NA) +
  
  # Flux:
  geom_point(data = subset(df1_2),
             aes(x = long, y = lat, color = propBiomass), size = .1) + 
  scale_color_viridis() +
  labs(color = "%") +
  
  # configure projection and plot domain
  coord_sf(xlim = lons, ylim = lats) +
  
  # formatting
  ylab("Lat") + xlab("Lon") +
  # ggtitle(label = "Biomass of functional groups") +
  theme_base() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10),
        plot.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = "bottom") +
  facet_wrap(vars(funGroup)) 

p
# ggsave("plots/propBiomass_NAtlantic_map.png", p, height = 70 , width = 80, units = "mm", scale = 3)


# Distinguishing by coast:
df1_2 <- df1 %>%
  mutate(propBiomass = case_when(is.na(propBiomass) ~ 0,
                                 T ~ propBiomass),
         totBiomass = case_when(is.na(totBiomass) ~ 0,
                                 T ~ totBiomass),
         depth = round(depth / 100) * 100) %>%
  group_by(depth, coast, funGroup) %>%
  reframe(sd_totBiomass = sd(totBiomass),
          totBiomass = mean(totBiomass),
          sd_propBiomass = sd(propBiomass),
          propBiomass = mean(propBiomass))

#------------

# Plot total biomass:
df1_3 <- df1 %>%
  mutate(propBiomass = case_when(is.na(propBiomass) ~ 0,
                                 T ~ propBiomass),
         totBiomass = case_when(is.na(totBiomass) ~ 0,
                                T ~ totBiomass),
         depth = round(depth / 100) * 100) %>%
  group_by(depth, funGroup) %>%
  reframe(sd_totBiomass = sd(totBiomass),
          totBiomass = mean(totBiomass),
          sd_propBiomass = sd(propBiomass),
          propBiomass = mean(propBiomass))

combined_groups <- df1_3 %>%
  group_by(depth) %>%
  summarise(
    funGroup = "Total",
    sd_totBiomass = NA,
    totBiomass = sum(totBiomass),
    sd_propBiomass = NA,
    propBiomass = sum(propBiomass)
  )

df1_3 <- rbind(df1_3, combined_groups)

df1_3 <-df1_3 %>%
  mutate(totBiomass = case_when(totBiomass < .1 ~ 0.1,
                                 T ~ totBiomass))


p <- ggplot() +
  geom_line(data = df1_3,
            aes(x = depth, y = totBiomass, color = funGroup, linetype = funGroup)) +
  scale_linetype_manual(values = c("demersals" = "solid", "largePel" = "solid", "midwater" = "solid", "smallPel" = "longdash", "Total" = "solid"),
                        labels = c("demersals" = "Demersals", "largePel" = "Large pelagics", "midwater" = "Midwater fish", "smallPel" = "Small pelagics", "Total" = "Total")) +
  # geom_ribbon(data = df1_3,
  #             aes(x = depth, ymin = ifelse((totBiomass - sd_totBiomass) < 0.1, 0.1, (totBiomass - sd_totBiomass)), ymax = totBiomass + sd_totBiomass, fill = funGroup),
  #             alpha = 0.05) +
  scale_fill_manual(values = c(as.vector(my_palette)[c(8, 6, 5, 4)], "black")) +
  scale_color_manual(values = c(as.vector(my_palette)[c(8, 6, 5, 4)], "black"),
                     labels = c("demersals" = "Demersals", "largePel" = "Large pelagics", "midwater" = "Midwater fish", "smallPel" = "Small pelagics", "Total" = "Total")) +
  geom_hline(yintercept = 0.1, alpha = 1, color = "white", size = 1) +
  xlab("Depth (m)") +
  ylab(expression(paste("Biomass (g ww ", m^-2, ")"))) +
  xlim(0, 2500) +
  scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30)) +
  theme_base() +
  theme(#legend.position = c(0.7, 0.7),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 3))

p
ggsave("plots/totBiomass_NAtlantic_log.png", p, height = 60 , width = 55, units = "mm", scale = 2)


#------------
# Plot diet as a function of 4 flux pathways:
df2 <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$diet, by = "id") %>%
  ungroup() %>%
  group_by(id, long, lat, coast, region, depth, flux_pathway) %>%
  summarise(ingestion = sum(ingestion)) %>%
  mutate(rel_ingestion = ingestion / sum(ingestion),
         rel_ingestion = case_when(flux_pathway == "Midwater" & rel_ingestion == 0 & depth <= 250 ~ NA,
                                 T ~ rel_ingestion))

df2_2 <- df2 %>%
  mutate(rel_ingestion = case_when(is.na(rel_ingestion) ~ 0,
                                   T ~ 100 * rel_ingestion)) 

df2_2 <- df2_2 %>%
  mutate(flux_pathway = case_when(flux_pathway == "Benthic" ~ "a)",
                                  flux_pathway == "Demersal" ~ "b)",
                                  flux_pathway == "Epipelagic" ~ "c)",
                                  flux_pathway == "Midwater" ~ "d)",
                                  T ~ flux_pathway)
    
  )

lons = c(-94, 11)
lats = c(12, 68)

p <- ggplot() +
  
  geom_sf(data = cntries, fill = "grey", colour = NA) +
  
  # Flux:
  geom_point(data = subset(df2_2),
             aes(x = long, y = lat, color = rel_ingestion), size = .1) + 
  scale_color_viridis() +
  labs(color = "%") +
  
  # configure projection and plot domain
  coord_sf(xlim = lons, ylim = lats) +
  
  # formatting
  ylab("Lat") + xlab("Lon") +
  facet_wrap(vars(flux_pathway)) +
  # ggtitle(label = "Source of C to demersal community") +
  theme_base() + 
  theme(axis.title = element_blank(),
    # axis.text = element_blank(),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.background = element_blank(),
    strip.text = element_text(hjust = 0)) 

p
ggsave("plots/prop_flux_to_demersals_NAtlantic_map.png", p, height = 70 , width = 80, units = "mm", scale = 3)


# Separating coasts:
df2_2 <- df2 %>%
  mutate(rel_ingestion = case_when(is.na(rel_ingestion) ~ 0,
                                 T ~ rel_ingestion),
         depth = round(depth / 100) * 100) %>%
  group_by(depth, coast, flux_pathway) %>%
  reframe(sd_rel_ingestion = sd(rel_ingestion),
          rel_ingestion = mean(rel_ingestion))


p <- ggplot() +
  geom_line(data = df2_2,
            aes(x = depth, y = rel_ingestion, color = flux_pathway, linetype = coast)) +
  geom_ribbon(data = df2_2, 
              aes(x = depth, ymin = ifelse((rel_ingestion - sd_rel_ingestion) < 0, 0, (rel_ingestion - sd_rel_ingestion)), ymax = rel_ingestion + sd_rel_ingestion, fill = flux_pathway, linetype = coast), 
              alpha = 0.05) +
  scale_fill_manual(values = c("#A65B00", as.vector(my_palette)[c(8, 5, 6)])) +
  scale_color_manual(values = c("#A65B00", as.vector(my_palette)[c(8, 5, 6)])) +
  geom_hline(yintercept = 0, alpha = .3) +
  xlim(0, 2500) +
  theme_base() +
  theme(plot.background = element_blank())

p
ggsave("plots/prop_flux_to_demersals_coasts.png", p, height = 60 , width = 80, units = "mm", scale = 3)


# Plot diet as a function of 2 flux pathways:
df4 <- cobalt_NAtlantic_slopes %>%
  left_join(my_sims$diet, by = "id") %>%
  mutate(flux_pathway2 = case_when(flux_pathway %in% c("Midwater", "Epipelagic") ~ "Pelagic",
                                   flux_pathway == "Demersal" ~  "Demersal",
                                   T ~ flux_pathway),
         detrital_flux = dflux / 9) %>% # Convert to C equivalents: 1 g C = 9 g WW
  group_by(id, long, lat, coast, depth, detrital_flux, flux_pathway2) %>%
  summarise(ingestion = sum(ingestion) / 9) %>% # Convert to C equivalents: 1 g C = 9 g WW
  ungroup()

temp <- df4 %>%
  filter(flux_pathway2 != "Demersal") %>%
  group_by(id) %>% 
  mutate(sum_ingestion = sum(ingestion), # If we remove demesal, how much of the remaining ingestion is benthic vs (midwater + pelagic)?
         rel_ingestion = ingestion / sum_ingestion) %>%
  filter(flux_pathway2 != "Benthic") %>% # We want to keep the ingested proportion of pelagic food
  dplyr::select(id, rel_ingestion) %>%
  right_join(df4, by = "id") %>%
  filter(flux_pathway2 == "Demersal") %>% # To calculate how much of the demersal flux comes originated in the pelagic system, we multiply it by the proportion of ingested pelagic food
  mutate(pelagic_from_demersal = rel_ingestion * ingestion) %>%
  dplyr::select(id, pelagic_from_demersal) %>%
  right_join(df4, by = "id") %>%
  filter(flux_pathway2 == "Pelagic") %>%
  mutate(flux = ingestion + pelagic_from_demersal, # Now we add the pelagic flux + the demersal flux that was originated in the pelagic system
         flux_pathway = "Active") %>%
  ungroup() %>%
  dplyr::select(id, long, lat, coast, depth, flux, flux_pathway)


df4 <- df4 %>%
  filter((row_number() - 1) %% 3 == 0) %>% # Filter only benthic source
  dplyr::select(id, long, lat, coast, depth, detrital_flux) %>% # Keep detrital flux and other characteristics of the gridpoint (no need of ingestion)
  rename(flux = detrital_flux) %>%
  mutate(flux_pathway = "Detritus") %>%
  rbind(temp) %>%
  mutate(flux = case_when(flux < 0 ~ 0,
                          T ~ flux))


df4_2 <- df4 %>%
  mutate(flux = case_when(is.na(flux) ~ 0,
                                   T ~ flux),
         depth = round(depth / 100) * 100) %>%
  group_by(depth, coast, flux_pathway) %>%
  reframe(sd_flux = sd(flux),
          flux = mean(flux))


p <- ggplot() +
  geom_line(data = df4_2,
            aes(x = depth, y = flux, color = flux_pathway, linetype = coast)) +
  geom_ribbon(data = df4_2, 
              aes(x = depth, ymin = ifelse((flux - sd_flux) < 0, 0, (flux - sd_flux)), ymax = flux + sd_flux, fill = flux_pathway, linetype = coast), 
              alpha = 0.05) +
  scale_fill_manual(values = c("#FF00FF", "gray50")) +
  scale_color_manual(values = c("#FF00FF", "gray50")) +
  geom_hline(yintercept = 0, alpha = .3) +
  xlim(0, 2500) +
  theme_base() 

p
#  ggsave("plots/flux_to_seafloor_coasts.png", p, height = 60 , width = 80, units = "mm", scale = 3)


# Make a map of the proportion of actife flux relative to detrital flux:
df5 <- df4 %>%
  pivot_wider(names_from = flux_pathway, values_from = flux) %>%
  mutate(flux_ratio = 100 * Active / Detritus)
  

# Set map limits
lons = c(-94, 11)
lats = c(12, 68)

p <- ggplot() +
  
  geom_sf(data = cntries, fill = "grey", colour = NA) +
  
  # Flux:
  geom_point(data = df5,
             aes(x = long, y = lat, color = flux_ratio),
             size = .25) + 
  scale_color_viridis() +
  labs(color = "%") +
  
  # configure projection and plot domain
  coord_sf(xlim = lons, ylim = lats) +
  
  # formatting
  ylab("Lat") + xlab("Lon") +
  # ggtitle(label = "Active flux to sea floor") +
  theme_base() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10),
        plot.background = element_blank(),
        legend.position = "bottom") 

p
ggsave("plots/flux_to_seafloor_map.png", p, height = 35 , width = 45, units = "mm", scale = 4)




#-------------------------------------------------------------------------------------------------------------
# Calculate proportion of active flux:
df4 %>%
  filter(depth >= 500, depth <= 1800) %>%
  group_by(coast, flux_pathway) %>%
  summarise(flux = sum(flux)) %>% 
  pivot_wider(names_from = flux_pathway, values_from = flux) %>%
  mutate(flux_ratio = 100 * Active / Detritus) 
  

# Calculate total active flux:
df4 <- df4 %>%
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
         c_area = area * flux)

# Calculate total active C and CO2 flux:
df4 %>%
  # filter(depth >= 500, depth <= 1800) %>%
  group_by(coast) %>%
  summarise(area_km2 = sum(area) / 1e6, # m^-2 to km^-2
            c_by_coast = sum(c_area) / 1e6,   # g to tonnes of C year^-1
            # CO2_by_coast = (44/12) * c_by_coast, # Million tonnes of CO2
            c_coast_km2 = c_by_coast * 1000 / area_km2) # tones km^-2 to kg km^-2
            # CO2_km2 = CO2_by_coast * 1e6 / area_km2) 

df4 %>%
  summarise(area_km2 = sum(area) / 1e6, # m^-2 to km^-2
            c_by_coast = sum(c_area) / 1e6,   # g to tonnes of C year^-1
            # CO2_by_coast = (44/12) * c_by_coast, # Million tonnes of CO2
            c_coast_km2 = c_by_coast * 1000 / area_km2) # tones km^-2 to kg km^-2
# CO2_km2 = CO2_by_coast * 1e6 / area_km2) 


#     END OF SCRIPT
############################

