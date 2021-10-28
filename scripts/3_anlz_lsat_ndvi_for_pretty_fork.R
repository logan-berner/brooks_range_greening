# This R script develops Landsat NDVI time series for sampling sites along a transect across the Brooks Range.
# Author: Logan Berner, NAU
# Date: 2021-04-08

# SET UP WORK SPACE ---------------------------------------------------------------------------
rm(list=ls())
require(data.table)
require(ggplot2)
setwd('C:/Users/Logan/My Drive/research/side_projects/dial_brooks_range/')

lsat.dt <- fread('output/pretty_fork_lsat_ndvi_max_timeseries.csv')

lsat.avg.dt <- lsat.dt[, .(ndvi.max.avg = mean(ndvi.max)), by = c('site','vegstack')]

lsat.smry.dt <- lsat.avg.dt[, .(ndvi.max.med = median(ndvi.max.avg)), by = c('vegstack')]
setorder(lsat.smry.dt, ndvi.max.med)

lsat.avg.dt[, vegstack := factor(vegstack, levels = lsat.smry.dt$vegstack)]

ggplot(lsat.avg.dt, aes(x = vegstack, y = ndvi.max.avg)) + 
  geom_boxplot() + coord_flip() + 
  labs(y = expression("Average Landsat NDVI"['max']~"from 2000 to 2020 (unitless)"))

ggsave('figures/pretty_fork_lsat_ndvi_averages_bwplot.jpg', width = 10, height = 8, units = 'in')  
