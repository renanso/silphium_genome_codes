################################################################################
# Geographic and Climate Data Visualization for Silphium Samples
################################################################################
# Purpose: Visualizes sample collection locations on maps with climate data layers
# Author: Renan Souza
# Created: R 4.4.2
# Output: Maps showing sample locations, hardiness zones, and climate variables
################################################################################

# Clear workspace
rm(list=ls())

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
# Spatial and mapping visualization
library("sf")                   # Simple features for spatial data handling
library("rnaturalearth")        # Natural Earth map data
library("rnaturalearthdata")    # Additional Natural Earth datasets
library("maps")                 # Map drawing utilities
library("ggspatial")            # Spatial extensions for ggplot2

# Data manipulation and visualization
library("ggplot2")              # Grammar of graphics for plotting
library("ggrepel")              # Better label placement to avoid overlaps
library("tidyverse")            # Data manipulation toolkit
library("data.table")           # Fast data frame operations

# Geographic data handling
library("geodata")              # Download and manage geographic data
library("terra")                # Raster and vector data manipulation
library("raster")               # Raster data analysis

# US-specific utilities
library("USAboundaries")        # US state and county boundaries (v0.3.0 required)
library("tools")                # Utility functions
library("conflicted")           # Manage function name conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")


# ============================================================================
# DATA INPUT
# ============================================================================
# Load sample coordinates
data <- read.csv("coordinates3.csv", header = TRUE)
# Expected columns: clone ID, longitude (long), latitude (lat)

# Set plotting theme to black and white
theme_set(theme_bw())

# ============================================================================
# SECTION 1: BASIC WORLD MAP WITH SAMPLE LOCATIONS
# ============================================================================
# Download world country boundaries at medium resolution
world <- ne_countries(scale = "medium", returnclass = "sf")

# Download rivers and lakes layer for geographic context
rivers50 <- ne_download(scale = 50, type = 'rivers_lake_centerlines', 
                        category = 'physical')
# Crop rivers to study region (central North America)
rivers_cropped <- st_crop(st_as_sf(rivers50), xmin = -95, xmax = -87,
                          ymin = 25, ymax = 48)


# Map 1: World map with sample points and labels
ggplot(data = world) +
  geom_sf() +
  geom_point(data = data, aes(x = long, y = lat), size = 2, shape = 23, fill = "darkred") +
  # Label points with clone ID, avoiding overlaps
  geom_label_repel(aes(x = long, y = lat, label = clone), size = 3, color = "black", 
                   data = data, max.overlaps = 200) +
  geom_sf(data = rivers_cropped, col = 'blue') +
  coord_sf(xlim = c(-110, -70), ylim = c(20, 55), expand = FALSE)

# ============================================================================
# SECTION 2: US STATE BOUNDARIES MAP
# ============================================================================
# Load US state boundaries for geographic reference
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

# Convert state names to title case for consistency
states$ID <- toTitleCase(states$ID)


# Map 2: US state map with sample points
ggplot(data = world) +
  geom_sf(data = states) +
  geom_point(data = data, aes(x = long, y = lat), size = 1, shape = 23, fill = "darkred") +
  geom_sf(data = rivers_cropped, col = 'blue') +
  coord_sf(xlim = c(-110, -70), ylim = c(20, 55), expand = FALSE)

# ============================================================================
# SECTION 3: USDA HARDINESS ZONES VISUALIZATION
# ============================================================================
# Load USDA Plant Hardiness Zone shapefile (pre-downloaded, 2012 edition)
# Download source: https://prism.oregonstate.edu/projects/public/phm/2012/
shp_hardness <- read_sf('phm_us_shp.shp')

# Filter to zones suitable for Silphium (9b and warmer)
# Silphium species have high cold tolerance up to zone 3-4
shp_hardness_subset <- shp_hardness %>%
  filter(str_detect(ZONE, '9b|10a|10b|11a|11b'))

# Load US state boundaries (lower 48 states only for reference)
usa <- us_boundaries(type = "state", resolution = "low") %>% 
  filter(!state_abbr %in% c("PR", "AK", "HI"))


# Map 3: Hardiness zone visualization with sample locations
ggplot() +
  geom_sf(data = shp_hardness, aes(fill = ZONE), color = NA) +
  geom_sf(data = usa, color = 'black', fill = NA) +
  theme_void() +  # Remove lat/long grid lines for cleaner map
  geom_point(data = data, aes(x = long, y = lat), size = 1.5, 
             shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-107, -70), ylim = c(20, 55), expand = FALSE)

# ============================================================================
# SECTION 4: CLIMATE DATA LAYERS FROM WORLDCLIM v2.1
# ============================================================================
# WorldClim is a free climate database with 10-minute resolution data
# Reference: https://www.worldclim.org/
# Variables obtained for the continental USA:
#   - tavg: Average temperature (°C)
#   - tmin: Minimum temperature (°C)
#   - tmax: Maximum temperature (°C)
#   - prec: Precipitation (mm)
#   - srad: Solar radiation (kJ m-2 day-1)
#   - elev: Elevation (m)

# NOTE: Update paths to match your local directory where WorldClim data is stored
data_tavg <- worldclim_country("USA", var = "tavg", res = 10, version = "2.1", 
                               path = "~/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim")
data_tmin <- worldclim_country("USA", var = "tmin", res = 10, version = "2.1", 
                               path = "~/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim")
data_tmax <- worldclim_country("USA", var = "tmax", res = 10, version = "2.1", 
                               path = "~/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim")
data_prec <- worldclim_country("USA", var = "prec", res = 10, version = "2.1", 
                               path = "~/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim")
data_srad <- worldclim_country("USA", var = "srad", res = 10, version = "2.1", 
                               path = "~/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim")
data_elev <- worldclim_country("USA", var = "elev", res = 10, version = "2.1", 
                               path = "~/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim")

vapr.files <- list.files("~/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim/climate/wc2.1_30s_vapr", ".tif", full.names=TRUE)
vapr <- stack(vapr.files)
data_vapr <-rast(vapr)

points <- vect(data, geom=c("long", "lat"), crs = "EPSG:4326")
conflicts_prefer(terra::extract)
values_tavg <- extract(data_tavg, points)
values_tmin <- extract(data_tmin, points)
values_tmax <- extract(data_tmax, points)
values_prec <- extract(data_prec, points)
values_srad <- extract(data_srad, points)
values_elev <- extract(data_elev, points)
values_vapr <- extract(data_vapr, points)

final_data<- cbind(data,values_tavg,values_tmin,values_tmax,values_prec,values_srad,values_elev,values_vapr)
write.csv(final_data, "climate_variable.csv")

#Heat Index calculation
#Thornthwaite 1948 on mean temperature for months 6,7,8
final_data<- read.csv("climate_variable.csv")
final_data$hit_og<-((final_data[,11]^1.514)/5)+((final_data[,12]^1.514)/5)+((final_data[,13]^1.514)/5)

#Thornthwaite 1948 on highest temperature for months 6,7,8
final_data$hit_hi<-((final_data[,37]^1.514)/5)+((final_data[,38]^1.514)/5)+((final_data[,39]^1.514)/5)

## Aridity index (ai) (Global-AI geodataset have been multiplied by a factor of 10,000)
##change the path to your Global-AI geodataset
ai.files <- list.files("/Users/renansouza/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/worldclim/climate/Global-AI_v3_monthly", ".tif", full.names=TRUE)
ai <- stack(ai.files)
data_ai <-rast(ai)
values_ai <- raster::extract(data_ai, points)
values_ai$summer_total<-rowSums(values_ai[,c(7:9)])
values_ai$summer_mean<-rowMeans(values_ai[,c(7:9)])

final_data2<-cbind(final_data[,c(2,3,4,85,86)],values_ai[,c(14,15)])

write.csv(final_data2, "final_gea2.csv")

#read.csv("final_gea_index2.csv")
###plot aridity index

#map <- stack("/My Drive/projects/007_silphium50K/genome_scans/gea/Global-AI_ET0_v3_annual/ai_v3_yr.tif")
map <- stack("/Users/renansouza/Library/CloudStorage/GoogleDrive-rsouza@hudsonalpha.org/My Drive/projects/007_silphium50K/genome_scans/gea/Global-AI_ET0_v3_annual/ai_v3_yr.tif")
#e <- as(extent(-130, -60, 20, 55), 'SpatialPolygons')
e <- as(extent(-100, -87, 29, 45), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
r <- crop(map, e)
tmax_Jan_09_rs<- stack(r)

#--- convert to data.frame ---#
tmax_Jan_09_df <-
  as.data.frame(tmax_Jan_09_rs, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit()

#remove 0
tmax_Jan_09_df2<- tmax_Jan_09_df[tmax_Jan_09_df$awi_pm_sr_yr >0, ]

colnames(tmax_Jan_09_df2)<- c("x","y","ai")
tmax_Jan_09_df2$ai2<-tmax_Jan_09_df2$ai/10000
#--- take a look ---#
head(tmax_Jan_09_df2)

png(filename="plot_ai_v3.png", width = 600, height = 750, units = "px")

### focus on midwest

g_tmax_map <- ggplot(data = tmax_Jan_09_df2) +
    geom_raster(aes(x = x, y = y, fill = `ai2`)) +
    #scale_fill_manual(values = mycolors) +
    scale_fill_viridis_c(option = "turbo", direction = -1, alpha = 0.9, 
                         begin = 0, end = 0.8, name = "Aridity", breaks = c(0.2,0.4,0.6,0.8,1,1.2)) +
    #theme_void() +
    theme( legend.position.inside =  c(0.8, 0.2),
           panel.background = element_blank(), 
           axis.line = element_line(colour = "black"),
           axis.title = element_text(size=22),
           axis.text = element_text(size=14),
           legend.text = element_text(size=17),
           legend.title = element_text(size=20)) +
   borders("state", xlim=c(-100, -88), ylim=c(30,45), colour = 'black') +
  geom_point(data = data, aes(x = long, y = lat), size = 1, shape = 23, fill = "black") +
  labs(x = "Longitude", y = "Latitude") +
  coord_cartesian(xlim=c(-99, -88), ylim=c(30,44))
g_tmax_map

dev.off()


###plot heat index mean

points2 <- vect(tmax_Jan_09_df2, geom=c("x", "y"), crs = "EPSG:4326")

values_tavg <- extract(data_tavg, points2)
values_tmin <- extract(data_tmin, points2)
values_tmax <- extract(data_tmax, points2)
values_prec <- extract(data_prec, points2)
values_srad <- extract(data_srad, points2)
values_elev <- extract(data_elev, points2)
values_vapr <- extract(data_vapr, points2)

plot_data<- cbind(tmax_Jan_09_df2,values_tavg,values_tmin,values_tmax,values_prec,values_srad,values_elev,values_vapr)

plot_data[,11]

#Heat Index calculation
#Thornthwaite 1948 on mean temperature for months 6,7,8
plot_data$hit_og<-((plot_data[,11]^1.514)/5)+((plot_data[,12]^1.514)/5)+((plot_data[,13]^1.514)/5)

#Thornthwaite 1948 on highest temperature for months 6,7,8
plot_data$hit_hi<-((plot_data[,37]^1.514)/5)+((plot_data[,38]^1.514)/5)+((plot_data[,39]^1.514)/5)


plot_data2<-data.frame(plot_data$x,plot_data$y, plot_data$hit_og, plot_data$hit_hi)

##Heat index mean
png(filename="plot_him_v3.png", width = 600, height = 750, units = "px")

g_tmax_map <- ggplot(data = plot_data2) +
  geom_raster(aes(x = plot_data.x, y = plot_data.y, fill = `plot_data.hit_og`)) +
  #scale_fill_manual(values = mycolors) +
  scale_fill_viridis_c(option = "turbo", direction = 1, alpha = 0.9, 
                       begin = 0, end = 0.8, name = "Heat\nIndex\nMean", breaks = c(0,20,40,60,80,99)) +
  #theme_void() +
  theme( legend.position.inside =  c(0.8, 0.2),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.title = element_text(size=22),
         axis.text = element_text(size=14),
         legend.text = element_text(size=17),
         legend.title = element_text(size=18)) +
  borders("state", xlim=c(-100, -88), ylim=c(30,45), colour = 'black') 

### focus on midwest
g_tmax_map +
  geom_point(data = data, aes(x = long, y = lat), size = 1, shape = 23, fill = "black") +
  labs(x = "Longitude", y = "Latitude") +
  coord_cartesian(xlim=c(-99, -88), ylim=c(30,44))

dev.off()

##Heat index high

png(filename="plot_hih_v3.png", width = 650, height = 750, units = "px")

g_tmax_map <- ggplot(data = plot_data2) +
  geom_raster(aes(x = plot_data.x, y = plot_data.y, fill = `plot_data.hit_hi`)) +
  #scale_fill_manual(values = mycolors) +
  scale_fill_viridis_c(option = "turbo", direction = 1, alpha = 0.9, 
                       begin = 0, end = 0.8, name = "Heat\nIndex\nHigh", breaks = c(0,20,40,60,80,100,120,140)) +
  #theme_void() +
  theme( legend.position.inside =  c(0.8, 0.2),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.title = element_text(size=22),
         axis.text = element_text(size=14),
         legend.text = element_text(size=17),
         legend.title = element_text(size=18)) +
  borders("state", xlim=c(-100, -88), ylim=c(30,45), colour = 'black') 

### focus on midwest
g_tmax_map +
  geom_point(data = data, aes(x = long, y = lat), size = 1, shape = 23, fill = "black") +
  labs(x = "Longitude", y = "Latitude") +
  coord_cartesian(xlim=c(-99, -88), ylim=c(30,44))

dev.off()

