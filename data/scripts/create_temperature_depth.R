## Making temperature and depth data ###
## VAST needs temperature and depth information for every
## "knot" in every year. We make those data here
## by extrapolating from survey data

library(tidyverse)
library(sf)
library(gstat)
library(raster)
library(rasterVis)
library(rgdal)
library(viridis)
library(here)

# ggplot theme
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14,color="black"),
        # axis.text=element_blank(),
        axis.text=element_text(family="sans",size=8,color="black"),
        panel.grid.major = element_blank(),
        panel.border=element_rect(color='black',fill=NA))
theme_set(plot_theme)

# Load survey data
opi_dat_long <- read_rds(here::here('data','processed','longform_opilio.rds'))
cod_dat <- read_rds(here::here('data','processed','cod_dat_clean.rds'))

# Pull out combined temperature/depth information
tdepth1 <- cod_dat %>% 
  dplyr::select(year,midlat,midlon,depth) %>% 
  rename(lat=midlat,lon=midlon)
tdepth2 <- opi_dat_long %>% 
  dplyr::select(year,lon,lat,depth,temp)

tdepth <- bind_rows(tdepth1,tdepth2) %>% 
  distinct() %>% 
  filter(!is.na(lat),!is.na(lon),year<2018)

tdepth_map <- tdepth %>% 
  ggplot(aes(lon,lat,col=temp))+
  geom_point()+
  facet_wrap(~year)
tdepth_map

# Use a convex hull outlining the survey area to extrapolate
# flip around the dateline
tdepth2 <- tdepth %>% mutate(lon_flipped=180 + ifelse(lon>0, lon-360, lon))
dat.sf <- st_as_sf(tdepth2,coords=c('lon_flipped','lat'),crs=4326) %>% st_transform("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")
dat.sp <- dat.sf %>% as_Spatial()

# convex hull raster of entire survey area
dat.ch <- st_convex_hull(dat.sf %>% st_union())
r <- raster(dat.ch %>% as_Spatial())
res(r) <- 10 #1 km resolution

## function to interpolate temperature data for a given year
interpolate_data <- function(df,yr,variable) {
  #subset data
  dat.sp <- df %>% 
    filter(year==yr,!is.na(temp)) %>% 
    as_Spatial()
  
  #interpolate using inverse distance weighting
  idm<-gstat(formula=as.formula(paste0(variable,'~','1')),locations=dat.sp)
  dat.idw <- interpolate(r,idm) %>% 
    mask(as_Spatial(dat.ch)) %>% 
    as.data.frame(xy=TRUE) %>% 
    mutate(year=yr)
  
  dat.idw
}

# interpolate temperature
purrr::map(1982:2017, ~{
  interpolate_data(df=dat.sf,yr=.x,variable='temp')
}) %>% bind_rows() -> nbt.interpolated
nbt.interpolated <- nbt.interpolated %>% 
  rename(temp=var1.pred) %>% 
  filter(!is.na(temp))

# interpolate depth
purrr::map(1982:2017, ~{
  interpolate_data(df=dat.sf,yr=.x,variable='depth')
}) %>% bind_rows() -> depth.interpolated
depth.interpolated <- depth.interpolated %>% 
  rename(depth=var1.pred) %>% 
  filter(!is.na(depth))
# we just want one constant depth (depth shouldn't change)
depth.interpolated <- depth.interpolated %>% 
  group_by(x,y) %>% 
  summarise(depth=mean(depth,na.rm=T)) %>% 
  mutate(year='all')

# one year's plot function (for checking)
# state outline
# ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
#   st_union() %>% 
#   st_transform(st_crs(dat.sf))

plot_data <- function(dat,variable) {
  dat %>% 
    ggplot()+
    geom_raster(aes_string("x","y",fill=variable),na.rm=T,alpha=0.8,interpolate = FALSE)+
    # geom_sf(data=ak,fill='gray80')+
    labs(x='',y='',fill='',title='')+
    coord_equal()+
    facet_wrap(~year)+
    scale_fill_viridis(na.value=NA,option="C")
}
plot_data(nbt.interpolated,variable='temp')
plot_data(depth.interpolated,variable='depth')

# save
write_rds(depth.interpolated,here::here('data','processed','depth_interpolated.rds'))
write_rds(nbt.interpolated,here::here('data','processed','nbt_interpolated.rds'))
