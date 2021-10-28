# This R script develops Landsat NDVI time series for sampling sites along a transect across the Brooks Range.
# Author: Logan Berner, NAU
# Date: 2021-10-28

# SET UP WORK SPACE ---------------------------------------------------------------------------
rm(list=ls())
require(lsatTS)
require(data.table)
require(rgdal)
require(sf)
require(ggplot2)
require(ggmap)
require(ggridges)
require(R.utils)
# setwd('C:/Users/Logan/My Drive/research/side_projects/dial_brooks_range/')

# PROCESS LANDSAT DATA -----------------------------------------------------------------------
# Load files
temp.files <- list.files('data/lsat_extracts/', pattern = 'pretty_fork', full.names = T)
lsat.dt <- do.call("rbind", lapply(temp.files, fread))

pts.sf <- st_read('output/pixel_walking_subset_tracks/pretty_fork_sample_100m_grid.shp')


# Parse data, filter to clear-sky observations, compute mean surface reflectance among pxls w/in each window around a site
lsat.dt <- lsat_general_prep(dt = lsat.dt)

# Clean the data, filtering out clouds, snow, water, radiometric and geometric errors
lsat.dt <- lsat_clean_data(lsat.dt, geom.max = 15, cloud.max = 80, sza.max = 60, 
                           filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)

# Check availability of Landsat data
data.summary.dt <- lsat_summarize_data_avail(lsat.dt)
data.summary.dt

# Compute NDVI
lsat.dt <- lsat_calc_spec_index(lsat.dt, 'ndvi')

# Cross-calibrate NDVI among sensors using RF models
lsat.dt <- lsat_calibrate_rf(lsat.dt, band = 'ndvi', doy.rng = 151:242, min.obs = 5, frac.train = 0.80, overwrite.col = T, outdir = 'output/ndvi_xcal/')

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, vi = 'ndvi', window.yrs = 5, window.min.obs = 10, vi.min = 0.1, spl.fit.outfile = F, progress = T)

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, vi = 'ndvi', min.frac.of.max = 0.75)

# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, vi = 'ndvi', min.obs = 10, reps = 2, min.frac.of.max = 0.75, outdir = NA)

# append vegetation class to each point
lsat.gs.dt$vegstack <- pts.sf$Vegstack[match(lsat.gs.dt$site, pts.sf$site)]

# Write out data.table with growing season summaries
fwrite(lsat.gs.dt, 'output/pretty_fork_lsat_ndvi_max_timeseries.csv')

# Compute temporal trend in annual maximum NDVI
lsat.trnd.dt <- lsat_calc_trend(lsat.gs.dt, 'ndvi.max', yrs = 2000:2020, sig = 0.10)
lsat.trnd.dt$vegstack <- pts.sf$Vegstack[match(lsat.trnd.dt$site, pts.sf$site)]

# write out csv
fwrite(lsat.trnd.dt, 'output/pretty_fork_lsat_ndvi_max_trends.csv')

# write out shapefile
lsat.trnd.sf = st_as_sf(lsat.trnd.dt, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
mkdirs('output/pixel_walking_lsat_trend_tracks/')
st_write(lsat.trnd.sf, dsn = 'output/pixel_walking_lsat_trend_tracks/lsat_/pretty_fork_lsat_ndvi_trends_wgs84', driver = 'ESRI Shapefile', append = F)


# FIGURES ================================================================================

# HISTROGRAM OF OVERALL % CHANGE
ggplot(lsat.trnd.dt, aes(total.change.pcnt, fill=..x..)) +
  geom_histogram(bins = 50, size = 0.25, color = 'gray20') +
  scale_fill_gradient2(low="darkorange4", mid='white', high="darkgreen", limits = c(-75,75), midpoint = 0) +
  labs(x = expression("Change in Landsat NDVI"['max']~"from 2000 to 2020 (%)"), y = 'Number of sampling sites') +
  theme_bw() + theme(legend.position = 'none') + xlim(-75, 75)

ggsave('figures/pretty_fork_lsat_ndvi_trend_hist.jpg', width = 5, height = 3, units = 'in')

# RIDGES PLOTS
ggplot(lsat.trnd.dt, aes(x=total.change.pcnt, y=vegstack, fill=..x..)) +
  geom_density_ridges_gradient() +
  scale_fill_gradient2(low="darkorange4", mid='white', high="darkgreen", midpoint = 0, name = '% change') +
  xlim(-50,50) + 
  labs(x = expression("Change in Landsat NDVI"['max']~"from 2000 to 2020 (%)"), y = 'Density') +
  theme_bw()

ggsave('figures/pretty_fork_lsat_ndvi_trend_ridges.jpg', width = 10, height = 8, units = 'in')


# 
# # histogram of percent change in Landsat NDVI from 2000 to 2020
# ggplot(lsat.trnd.dt, aes(total.change.pcnt, fill=..x..)) +
#   geom_histogram(bins = 30, size = 0.25, color = 'gray20') + 
#   scale_fill_gradientn(colours = c('brown4','brown1','white','green','green4'),
#                        breaks = c(-50,-20, -5, 0, 5, 20, 50), limits = c(-50,50)) + 
#   labs(x = expression("Change in Landsat NDVI"['max']~"from 2000 to 2020 (%)"), y = 'Number of sampling sites') + 
#   theme_bw() + theme(legend.position = 'none') + xlim(-50, 50)
# 
# ggsave('figures/pretty_fork_lsat_ndvi_trend_hist.jpg', width = 5, height = 3, units = 'in')
# 
# # transparent background
# ggplot(lsat.trnd.dt, aes(total.change.pcnt, fill=..x..)) +
#   geom_histogram(bins = 30, size = 0.25, color = 'gray20') + xlim(-50, 100) +
#   scale_fill_gradient2(low="darkorange4", mid='ivory', high="darkgreen", limits = c(-50,100), midpoint = 0) +
#   labs(x = "Change in Landsat NDVI from 2000 to 2020 (%)", y = 'Number of sampling sites') +
#   theme(legend.position = 'none', panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA))

# # density plot
# lsat.trnd.dens.dt <- density(lsat.trnd.dt$total.change.pcnt, n=100000)
# lsat.trnd.dens.dt <- data.table(x=lsat.trnd.dens.dt$x, y=lsat.trnd.dens.dt$y)
# lsat.trnd.dens.dt <- lsat.trnd.dens.dt[x >= -50 & x <= 100]
# ggplot(lsat.trnd.dens.dt, aes(x, y)) +
#   geom_segment(aes(xend = x, yend = 0, color = x)) +
#   geom_line() +
#   scale_color_gradient2(low="brown4", mid='white', high="green", limits = c(-50,100)) +
#   labs(x = "Percent change in NDVI from 2000 to 2020", y = 'Density') +
#   theme_bw()