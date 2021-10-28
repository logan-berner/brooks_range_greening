# This R script takes the Pretty Fork GPS track from Romian Dial and 
# extracts Landsat data for these points using GEE. 
# A helpful resource for working with spatial data in R: https://www.maths.lancs.ac.uk/~rowlings/Teaching/UseR2012/cheatsheet.html
# Author: Logan Berner, NAU
# Date: 2021-09-27

# SET UP -------------------------------------------------------------------------------------------
rm(list=ls())
require(rgdal)
require(sf)
require(lsatTS)
require(rgee)
require(dplyr)
require(raster)

# require(devtools)
# install_github('https://github.com/logan-berner/lsatTS/', auth_token = '')

# setwd('C:/Users/Logan/My Drive/research/side_projects/dial_brooks_range/')
ak.aea.proj <- '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

# Read in transect shapefile
transect.pts <- st_read('data/pixel_walking/PrettyFork.shp')

# reproject to aea
transect.pts <- transect.pts %>% st_transform(crs = ak.aea.proj)

# create raster over region
grid.resolution.m = 100 
coords <- st_coordinates(transect.pts)
aoi.r <- raster(xmn = min(coords[,1]-100), xmx = max(coords[,1]+grid.resolution.m),
       ymn = min(coords[,2]-100), ymx = max(coords[,2]+grid.resolution.m), resolution = grid.resolution.m, crs = ak.aea.proj)
aoi.r[] <- 1:ncell(aoi.r)

plot(aoi.r)
plot(transect.pts, add=T)

# subsample points based on raster grid cells 
transect.pts$grid.cell <- raster::extract(aoi.r, transect.pts)
transect.pts <- transect.pts %>% group_by(grid.cell) %>% do(sample_n(.,1)) %>% as.data.frame %>% st_as_sf
plot(transect.pts, add=T, color = 'black')

# Add unique ID to each sample point
transect.pts <- transect.pts %>% st_sf %>% mutate(sample_id = paste0('pt', 1:nrow(st_coordinates(transect.pts))), 
                                           region = c("pretty_fork"))

# transect.pts <- transect.pts[,c(6,7)]

# write out 
outname <- paste0('output/pixel_walking_subset_tracks/pretty_fork_sample_',grid.resolution.m,'m_grid.shp')
st_write(transect.pts, dsn = outname, append = F)

# Extract Landsat time-series for points along transect
ee_Initialize()
lsat_export_ts(transect.pts, file_prefix = 'lsat_pretty_fork_test_', drive_export_dir = 'earth_engine', 
               startJulian = 142, endJulian = 243)

## NOTE-- ONCE DATA EXTRACTION IS COMPELET, MOVE THE CSV FILES OUT OF YOUR GOOGLE DRIVE AND INTO THE /data/lsat_extracs/ folder 

# END SCRIPT ----------------------------------------------------------------------------------------