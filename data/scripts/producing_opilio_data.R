# Produce snow crab data for use in VAST

library(tidyverse)
library(here)

t<-proc.time()
data<- read_csv(here::here('data','raw','opilio.csv'),skip = 5,col_types = 'dddddddcccddddcddcddddddddddcddddddddddddddddc')
mat_df<-read_csv(here::here('data','raw','maturity.csv'),col_types = cols())

# haul join key (for comparing to other species, e.g. cod, collected in the same survey)
haul_join_key <- data %>% select(HAULJOIN,AKFIN_SURVEY_YEAR,GIS_STATION,HAUL,MID_LATITUDE,MID_LONGITUDE,AREA_SWEPT_VARIABLE) %>% 
  distinct() %>%
  rename(hauljoin=HAULJOIN,year=AKFIN_SURVEY_YEAR,midlat=MID_LATITUDE,
         midlon=MID_LONGITUDE,station=GIS_STATION,haul_number=HAUL,
         AreaSwept_km2=AREA_SWEPT_VARIABLE) %>% 
  #add another id column, combining year, station, and haul number
  unite(towid,year,station,haul_number,sep="_",remove=F)
write_rds(haul_join_key,here::here('data','processed','haul_join_key.rds'))

# data before 1982 unreliable
use_yr<-seq(1982,2017)

data <- data %>% filter(AKFIN_SURVEY_YEAR>1981)

binsize=5
sizes		 <-seq(27.5,132.5,binsize)
# size bin matching key
sizesdf <- tibble(dn=seq(25,130,5),up=seq(30,135,5),size=sizes)

#=====================PROCESS RAW DATA====================#

opi.dat <- data %>% 
  
  # remove obs with no sex
  filter(!is.na(SEX)) %>% 
  
  # select variables of interest
  select(AKFIN_SURVEY_YEAR,GIS_STATION,HAUL,MID_LATITUDE,MID_LONGITUDE, AREA_SWEPT_VARIABLE,
         GEAR_TEMPERATURE,GEAR_DEPTH,BOTTOM_DEPTH,SEX,WIDTH,CLUTCH_SIZE,SHELL_CONDITION,
         CALCULATED_WEIGHT,SAMPLING_FACTOR) %>% 
  
  #rename some variables
  rename(year=AKFIN_SURVEY_YEAR,station=GIS_STATION,haul_number=HAUL,lat=MID_LATITUDE,lon=MID_LONGITUDE,
         area_km2=AREA_SWEPT_VARIABLE,temp=GEAR_TEMPERATURE,depth_gear=GEAR_DEPTH,depth=BOTTOM_DEPTH) %>% 
  #add another id column, combining year, station, and haul number
  unite(towid,year,station,haul_number,sep="_",remove=F) %>% 
  
  # sex as character
  mutate(sex=ifelse(SEX==2,"Female","Male")) %>% 

  # add size class bins by matching to the "key" above 
  # by rounding down and up to the nearest multiple of 5
  mutate(dn=5*floor((WIDTH+0.1)/5),up=5*ceiling((WIDTH+0.1)/5)) %>% 
  left_join(sizesdf,by=c("dn","up")) %>% 
  filter(!is.na(size)) %>% 
  select(-dn,-up) %>% 
  
  # join the male maturity dataset by size bin for use in calculating probability of male maturity below
  left_join(mat_df,by=c('size'='Size')) %>% 
  
  # add probability of maturity for both males and females
  mutate(prob_mature=case_when(
    SEX==2 & CLUTCH_SIZE>0 ~ 1,
    SEX==2 & CLUTCH_SIZE==0 ~ 0,
    # probs of male maturity come from 'mat_df' dataframe
    SEX==1 & SHELL_CONDITION<=2 ~ New,
    SEX==1 & SHELL_CONDITION >2 ~ Old
  )
  ) %>% 
  mutate(prob_immature=1-prob_mature) %>% 
  
  # summarize values of interest
  group_by(towid,station,lat,lon,year,sex,size) %>% 
  summarise(
    # number and weight of mature and immature crabs
    num_mature=sum(SAMPLING_FACTOR*prob_mature,na.rm=T),
    num_immature=sum(SAMPLING_FACTOR*prob_immature,na.rm=T),
    weight_mature=sum(CALCULATED_WEIGHT*SAMPLING_FACTOR*prob_mature,na.rm=T),
    weight_immature=sum(CALCULATED_WEIGHT*SAMPLING_FACTOR*prob_immature,na.rm=T),
    
    #other useful descriptives
    area_km2 = first(area_km2),
    depth_gear=first(depth_gear),
    depth=first(depth),
    temp=first(temp)
  ) %>% 
  ungroup()

## convert to long form
opi_dat_long <- opi.dat %>%
  pivot_longer(names_to = "category",values_to="value",num_mature:weight_immature) %>% 
  
  #rename some stuff and add units
  mutate(maturity=ifelse(grepl('immature',category),"Immature","Mature"),
         units=ifelse(grepl('weight',category),"kilos","numbers")) %>% 
  
  #turn implicit missing values into explicit missing values
  #e.g., for every station/year that is in the data (i.e., was surveyed), but
  #did not collect any immature females, put in a zero for that observation
  complete(sex,size,maturity,units,nesting(towid,station,year,area_km2,lon,lat,depth,depth_gear,temp),fill=list(value=0)) %>% 
  ungroup() %>% 
  
  #select final variables of interest
  select(towid,station,year,area_km2,lon,lat,depth,depth_gear,temp,sex,size,maturity,value,units) %>% 
  arrange(station,year,towid)

## NOTE: We arrive at a different total number of observations than in the previous interation:
## this seems to be because of the way we're counting NAs; the previous analysis counted year/station combinations
## that are present in the raw data, but without any sample info (sex, weight, etc.), as valid observations
## THIS new analysis does not do that (i.e., it throws out station/year combinations for which there is no useful data,
## and does not assume that those are zeros)

## SAVE
write_rds(opi_dat_long,here::here('data','processed','longform_opilio.rds'))
proc.time()-t