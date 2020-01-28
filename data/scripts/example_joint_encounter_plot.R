## Example time series output of joint encounter rates

# Defined as product of probability of occurrence of 2 categories
# Required data for this script to run:
## VAST outputs: 
## Report (TMB Obj$report())
## dat (Data that was originally passed to VAST)

library(tidyverse)

## ggplot theme (just for aesthetic look)
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14,color="black"),
        # axis.text=element_blank(),
        axis.text=element_text(family="sans",size=8,color="black"),
        panel.grid.major = element_blank(),
        panel.border=element_rect(color='black',fill=NA),
        strip.text = element_text(size=8))
theme_set(plot_theme)


# load(here::here('data','VAST output',"Save.RData"))
# Report <- Save$Report
# dat <- Save$Data

# Estimated encounter probability by knot/category/year
est <- Report$R1_gcy

# Choose species to calculate joint occurence
# Must be names of categories from original VAST data object
species_1 <- "Opilio Immature"
species_2 <- "Medium Cod"

knots = nrow(Report$R1_gcy)
cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
years <- seq(min(dat$year),max(dat$year))

# Organize species data into long-form data frame instead of 3-D array
spp1_dat <- est[,which(cats==species_1),] %>% as_tibble() %>% set_names(as.character(years)) %>% 
  mutate(knot=row_number()) %>% pivot_longer(-knot,names_to = "year",values_to="enc_1") %>% 
  mutate(year=as.numeric(year))
spp2_dat <- est[,which(cats==species_2),] %>% as_tibble() %>% set_names(as.character(years)) %>% 
  mutate(knot=row_number()) %>% pivot_longer(-knot,names_to = "year",values_to="enc_2")%>% 
  mutate(year=as.numeric(year))

# Join the two species' data by year and knot
df <- spp1_dat %>% left_join(spp2_dat,by=c('year','knot')) %>% 
  # calc joint encounter
  mutate(joint_enc=enc_1*enc_2) %>% 
  ungroup()%>% 
  # calculate means by year
  group_by(year) %>% 
  summarise(mean_joint_enc=mean(joint_enc,na.rm=T)) %>% 
  ungroup()

# Plot time series of mean joint encounter rate across knots for each year
out <- df %>% 
  ggplot(aes(year,mean_joint_enc))+
  geom_line()+
  geom_point()+
  labs(x="Year",y="Mean Joint Encounter Rate",title="")+
  ylim(c(0.8,1))+
  scale_x_continuous(breaks=seq(1980,2020,by=5))+
  theme(axis.text.x = element_text(angle=90))
out

## Implementation with PESC data
load(here::here('data','processed','PESCs_PCod.RData'))
pesc <- Save

pesc_knots <- pesc$PESC_Knots
# convert PESC to long form
pesc_knots %<>% mutate(knot=row_number()) %>% 
  pivot_longer(names_to='year',values_to = 'stomach_contents',1:31) %>% 
  mutate(year=as.numeric(year))%>% 
  # calculate means by year
  group_by(year) %>% 
  summarise(mean_pesc=mean(stomach_contents,na.rm=T)) %>% 
  ungroup() 

# Plot time series of mean PESC across knots for each year
out2 <- pesc_knots %>% 
  ggplot(aes(year,mean_pesc))+
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks=seq(1980,2020,by=5))+
  labs(x="Year",y="Mean PESC",title="")+
  theme(axis.text.x = element_text(angle=90))
out2

## And putting them together with two y-axes
combined_df <- df %>% 
  left_join(pesc_knots,by=c('year'))

combined_plot <- combined_df %>% 
  ggplot()+
  geom_line(aes(year,mean_joint_enc))+
  geom_point(aes(year,mean_joint_enc))+
  # weird fudging to get the secondary axis to look right
  geom_line(aes(year,mean_pesc*0.2/350+0.75),col='red')+
  geom_point(aes(year,mean_pesc*0.2/350+0.75),col='red')+
  scale_y_continuous(limits=c(0.75,1),sec.axis=sec_axis(~(.-0.75)*1750,name="Mean PESC"))+
  scale_x_continuous(breaks=seq(1980,2020,by=5))+
  labs(x="Year",y="Mean Joint Encounter Rate",title="Joint Encounter Rate and\nMean Estimated PESC")+
  theme(axis.text.x = element_text(angle=90),
        axis.title.y.right = element_text(color='red'))
combined_plot

ggsave(here::here('data','VAST output','plots','joint_enc_vs_pesc.png'),h=6,w=6)

## Split by depth zone
depth.interpolated <- read_rds(here::here('data','processed','depth_interpolated.rds'))

Spatial_List <- read_rds(here::here('data','VAST output',"Spatial_List.rds"))

find_depth_x <- function(Spatial_List) {
  # locations of knots from Spatial_List object
  loc_query <- Spatial_List[['loc_x']]
  # locations of interpolated depth info
  loc_covars <- depth.interpolated %>% distinct(x,y)
  # nearest neighbors (knots and interpolated covariates)
  which_nn <- nn2(data=loc_covars,query=loc_query,k=1)$nn.idx[,1]
  # for each knot in each year, find nearest neighbor from the interpolated datasets
  X_x <- tibble(knot=1:nrow(loc_query),depth=depth.interpolated$depth[which_nn])

  return(X_x)
}
depth_x <- find_depth_x(Spatial_List=Spatial_List) %>% 
  # add zonal definition
  mutate(depth_zone=case_when(
    depth<50 ~ "Inner",
    depth>=50 & depth<=100 ~ "Middle",
    depth>100 ~"Outer"
  ))

# add zone information to data
df <- spp1_dat %>% left_join(spp2_dat,by=c('year','knot')) %>% 
  # calc joint encounter
  mutate(joint_enc=enc_1*enc_2) %>% 
  ungroup()%>% 
  left_join(depth_x,by='knot') %>% 
  # calculate means by year
  group_by(year,depth_zone) %>% 
  summarise(mean_joint_enc=mean(joint_enc,na.rm=T)) %>% 
  ungroup()

pesc_knots <- pesc$PESC_Knots %>% 
  mutate(knot=row_number()) %>% 
  pivot_longer(names_to='year',values_to = 'stomach_contents',1:31) %>% 
  mutate(year=as.numeric(year))%>% 
  left_join(depth_x,by='knot') %>% 
  # calculate means by year
  group_by(year,depth_zone) %>% 
  summarise(mean_pesc=mean(stomach_contents,na.rm=T)) %>% 
  ungroup() 

combined_df <- df %>% 
  left_join(pesc_knots,by=c('year','depth_zone'))

combined_plot <- combined_df %>% 
  ggplot()+
  geom_line(aes(year,mean_joint_enc,linetype=depth_zone))+
  geom_point(aes(year,mean_joint_enc))+
  # weird fudging to get the secondary axis to look right
  geom_line(aes(year,mean_pesc*0.2/350+0.75,linetype=depth_zone),col='seagreen')+
  geom_point(aes(year,mean_pesc*0.2/350+0.75),col='seagreen')+
  scale_y_continuous(limits=c(0.75,1),sec.axis=sec_axis(~(.-0.75)*1750,name="Mean PESC"))+
  scale_x_continuous(breaks=seq(1980,2020,by=5))+
  scale_linetype(name="Depth Zone")+
  labs(x="Year",y="Mean Joint Encounter Rate",title="Joint Encounter Rate and\nMean Estimated PESC")+
  theme(axis.text.x = element_text(angle=90),
        axis.title.y.right = element_text(color='seagreen'))
combined_plot

ggsave(here::here('data','VAST output','plots','joint_enc_vs_pesc_zones.png'),combined_plot,h=6,w=8)

## Map knot numbers (for reference)
library(sf)
ak <- read_sf(here('data','spatial','cb_2017_02_anrc_500k.shp')) %>% 
  st_union() %>% 
  st_transform(26904)

loc_knots <- Spatial_List$loc_x %>% as_tibble() %>% set_names(c("x","y")) %>% 
  mutate(knot=row_number()) %>% 
  st_as_sf(coords=c('x','y'),crs=st_crs(ak))
bbox <- st_bbox(loc_knots)

knots_plot <- ggplot()+ 
  geom_sf(data=ak,fill='gray30',col=NA)+
  geom_sf_text(data=loc_knots,aes(label=knot),col='gray30',size=2)+
  coord_sf(crs = st_crs(ak), datum = NA) +
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  theme(axis.text = element_blank())+
  labs(x="",y="")

ggsave(here::here('data','VAST output','plots','knot_locs.png'),knots_plot,w=4,h=4)

# testing- low occurrence places of immature crab
# number of years <0.1 prob of occurrence
test <- spp1_dat %>% mutate(is_low=enc_1<0.1) %>% group_by(knot) %>% 
  summarise(is_low=sum(is_low)) %>% ungroup() %>% 
  mutate(is_low_alot=is_low>5) %>% 
  right_join(loc_knots) %>% st_as_sf()
library(viridis)
ggplot()+ 
  geom_sf(data=ak,fill='gray30',col=NA)+
  geom_sf_text(data=test,aes(label=knot,col=is_low_alot),size=2)+
  coord_sf(crs = st_crs(ak), datum = NA) +
  # scale_color_viridis()+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  theme(axis.text = element_blank())+
  labs(x="",y="")
