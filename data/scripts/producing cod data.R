# Compile cod data
library(tidyverse)
library(here)

cod_dat_raw <- read_csv(here::here('data','raw','race_length_by_haul.csv'),skip=7)

# remove unnecessary columns and rename
cod_dat_rn <- cod_dat_raw %>%
  select(Year,`Haul Join ID`,`Starting Latitude (dd)`,`Starting Longitude (dd)`,`Ending Latitude (dd)`,`Ending Longitude (dd)`, Stratum,
         `Satisfactory Gear Performance`,`Bottom Depth`,`Length (mm)`,`Frequency`,`Sample Type`,`Length Type`) %>%
  rename(hauljoin=`Haul Join ID`,year=Year,startlat=`Starting Latitude (dd)`,startlon=`Starting Longitude (dd)`,
          endlat=`Ending Latitude (dd)`,endlon=`Ending Longitude (dd)`,gear_satisfactory=`Satisfactory Gear Performance`,depth=`Bottom Depth`,length=`Length (mm)`,
         frequency=`Frequency`)

# Join appropriate GIS stations and empty hauls (cross-referenced from snow crab data)
haul_join_key <- read_rds(here::here('data','processed','haul_join_key.rds'))

cod_dat <- cod_dat_rn %>% 
  full_join(haul_join_key,by=c("hauljoin","year")) %>% 
  # fill in zeroes for unrepresented hauls, and use 0.01 for area swept
  mutate(frequency=coalesce(frequency,0)) %>% 
  rename(area_km2=AreaSwept_km2) %>% mutate(area_km2=coalesce(area_km2,0.01)) %>% 
  # for those hauls missing midlat/midlon, fill with startlat/startlon
  mutate(midlon=coalesce(midlon,startlon),midlat=coalesce(midlat,startlat)) %>% 

  filter(year>1981,year<2018) %>% 
  filter(!is.na(length))
  
print(paste("The number of missing latitudes is",sum(is.na(cod_dat$midlat))))
print(paste("The number of missing longitudes is",sum(is.na(cod_dat$midlon))))
print(paste("The number of missing stations IDs is",sum(is.na(cod_dat$station))))
print(paste("The number of missing years is",sum(is.na(cod_dat$year))))


# Aggregate by station and year
cod_dat_clean <- cod_dat %>%
  
  select(towid,station,year,area_km2,midlon,midlat,depth,length,frequency)

write_rds(cod_dat_clean,here::here('data','processed','cod_dat_clean.rds'))
