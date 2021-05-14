# Produce manuscript figures
# 1. Factor Maps
# 2. Factor loadings (PCA)
# 3. 2-dimensional factors
# 4. Correlations between classes
# 5. Lagged correlations (time series)
# 6. Time series indices

library(tidyverse)
library(here)
library(sf)
library(VAST)
library(abind)
library(viridis)
library(ggsci)
library(RANN)
library(cowplot)
library(marmap)
library(raster)
select <- dplyr::select

# ggplot theme
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

# basemap
# library(rnaturalearth)
# library(rnaturalearthdata)
# ak <- ne_states(country='United States of America',geounit="Alaska",returnclass = 'sf')

ak <- read_sf(here('data','spatial','cb_2017_02_anrc_500k.shp')) %>% 
  st_union() %>% 
  st_transform(26904)
          
## data for plotting
fp = here::here('data','VAST output')
load(here::here('data','VAST output','plots',"derived_quantities.Rdata"))
Extrapolation_List <- read_rds(paste0(fp,"/Extrapolation_List.rds"))
Spatial_List <- read_rds(paste0(fp,"/Spatial_List.rds"))
load(paste0(fp,"/Save.RData"))
load(paste0(fp,"/TMBData.RData"))
Report <- Save$Report
SDr <- Save$SDReport
dat=Save$Data
Year_Set = seq(min(dat$year),max(dat$year))
Years2Include = which( Year_Set %in% sort(unique(dat$year)))
Index = plot_biomass_index( DirName=fp, TmbData=TmbData, Sdreport=SDr, Year_Set=Year_Set, 
                            category_names=levels(dat$spp),Years2Include=Year_Set, use_biascorr=TRUE )

MapDetails_List = make_map_info( "Region"="Eastern_Bering_Sea",spatial_list = Spatial_List, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# factors <- plot_factors( Report=Save$Report, ParHat=Save$ParHat, Data=TmbData, SD=SDr, mapdetails_list=MapDetails_List, 
#                          Year_Set=Year_Set, category_names=levels(dat$spp), plotdir=paste0(here::here('data','VAST output','plots'),'/'))
# write_rds(factors,here::here('data','VAST output','factors.rds'))
factors <- read_rds(here::here('data','VAST output','factors.rds'))

depth.interpolated <- read_rds(here::here('data','processed','depth_interpolated.rds'))
nbt.interpolated <- read_rds(here::here('data','processed','nbt_interpolated.rds'))

find_temperature <- function(Spatial_List,dat) {
  # locations of knots from Spatial_List object
  loc_query <- Spatial_List[['loc_i']]
  # temperature summarized across yrs
  mbt <- nbt.interpolated %>% 
    group_by(x,y) %>% 
    summarise(meantemp=mean(temp,na.rm=T)) %>% 
    ungroup()
  # locations of interpolated temperature and depth info
  loc_covars <- nbt.interpolated %>% distinct(x,y)
  # nearest neighbors (knots and interpolated covariates)
  which_nn <- nn2(data=loc_covars,query=loc_query,k=1)$nn.idx[,1]
  # for each knot in each year, find nearest neighbor from the interpolated datasets
  X_x <- tibble(idx=1:nrow(loc_query),meantemp=mbt$meantemp[which_nn])
  # spatial
  X_i <- dat %>% select(Lon,Lat) %>% bind_cols(select(X_x,meantemp))
  return(X_i)
}

# bathymetry and mean temperature
bering_bathy_rast <- getNOAA.bathy(lon1 = -178,lon2 = -158,54,63,resolution = 10) %>% 
  marmap::as.raster()
bering_bathy <- bering_bathy_rast %>% rasterToContour(levels=c(-50,-100)) %>% st_as_sf() %>% 
  st_transform(26904)
lims.bathy <- st_bbox(bering_bathy)
# mean_nbt <- find_temperature(Spatial_List,dat) %>% st_as_sf(coords=c("lon","lat"),crs=4326) %>% 
#   st_transform(26904)
# mean_nbt_rast <- rasterize(mean_nbt,
#                            y=projectRaster(bering_bathy_rast,crs= "+proj=utm +zone=4 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"),
#                            field="meantemp",
#                            fun=mean)
# mean_nbt_xy <- mean_nbt_rast %>% as.data.frame(xy=T)
dat_i <- dat %>% select(lon,lat) %>% st_as_sf(coords=c('lon','lat'),crs=4326) %>% st_transform(26904) %>% st_coordinates() %>% 
  as.data.frame()

bathy_samples <- ggplot()+  
  geom_sf(data=ak,fill='gray30',col=NA)+
  geom_point(data=dat_i,aes(X,Y),col='gray30',size=1)+
  geom_sf(data = bering_bathy, col='darkgreen',size=1)+
  coord_sf(crs = st_crs(ak), datum = NA) +
  xlim(lims.bathy[1],lims.bathy[3])+ylim(lims.bathy[2],lims.bathy[4])+
  scale_fill_viridis(discrete=FALSE,na.value="transparent")+
  # scale_color_viridis(option="D")+
  theme(axis.text = element_blank(),
        panel.spacing.x = unit(0,"pt"),
        panel.spacing.y=unit(0,"pt"))+
  labs(x="",y="")
bathy_samples
ggsave(bathy_samples,file=paste0(fp,"bathy_samples.png"),width = 4,h=4)

## Raw data size frequency distributions
#### Import Data ####
load("data/longform_opilio2.Rdata")
load("data/cod_dat_size_number.Rdata")

# size frequency of immature crabs
opi_imm_n <- opi_dat_long %>% filter(maturity=="Immature",units=="numbers",size<58) %>% 
  group_by(size) %>% 
  summarise(totn=sum(value)) %>% 
  ungroup()
opi_imm_sf<-opi_imm_n %>% ggplot(aes(size,totn))+
  geom_bar(stat='identity')+labs(x="Carapace Width (mm)",y="Total Number",title="Immature Snow Crab\nSize Frequency")+
  theme(plot.title = element_text(size=10))

opi_spawn <- opi_dat_long %>% filter(maturity=="Mature", sex=="Female",units=="numbers") %>% 
  group_by(size) %>% 
  summarise(totn=sum(value)) %>% 
  ungroup()
opi_spawn_sf<-opi_spawn %>% ggplot(aes(size,totn))+
  geom_bar(stat='identity')+
  lims(x=c(0,80))+
  labs(x="Carapace Width (mm)",y="Total Number",title="Spawner Snow Crab\nSize Frequency")+
  theme(plot.title = element_text(size=10))
gadus <- cod_dat %>% 
  mutate(size_class=case_when(
    length<=200 ~ "small cod",
    length>200 & length<=800 ~ "medium cod",
    length>800 ~ "large cod"
  )) %>% 
  group_by(length,size_class) %>% 
  summarise(abun=sum(frequency,na.rm=T)) %>% 
  ungroup()
gadus_sf <- gadus %>%
  ggplot(aes(length,abun,group=size_class,fill=size_class))+
  geom_col(alpha=0.5)+
  scale_fill_npg()+
  theme(legend.position = c(0.6,0.8),
        legend.text=element_text(size=8))+
  theme(plot.title = element_text(size=10))+
  labs(x="Fork Length (mm)",y="Total Number",title="Pacific Cod\nSize Frequency",fill="")
gadus_sf

size_freq_all <- plot_grid(gadus_sf,opi_imm_sf,opi_spawn_sf,nrow=1)
ggsave(size_freq_all,file=paste0(fp,"size_freq_all.png"),width = 8,h=4)

### Plots derived from model outputs ###

map_points <- function(Extrapolation_List) {
  pts <- Extrapolation_List$Data_Extrap[,c('Lon','Lat','Area_in_survey_km2')] %>% 
    st_as_sf(coords=c('Lon','Lat'),crs=4326) %>% 
    #convert to AK UTM zone
    st_transform(26904)
  pts[,c("E_km","N_km")] <- st_coordinates(pts)
  pts$idx <-as.integer(Spatial_List$PolygonList$NN_Extrap$nn.idx)
  pts <- filter(pts,Area_in_survey_km2>0) %>% select(-Area_in_survey_km2)
  
  return(pts)
}

map_density <- function(dat,Report,Extrapolation_list,category,years='all') {
  
  dens <- Report$D_gcy
  
  cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  pts <- map_points(Extrapolation_List)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$year),max(dat$year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$year)))]
  knots <- dim(dens)[1]
  
  cats <- tools::toTitleCase(levels(dat$spp))

  df <- tibble(knot=rep(1:knots,length(Years2Include)),year=rep(Years2Include,each=knots),dens=log(as.numeric(dens[,category,])))
  sp_df<-pts %>% left_join(df,by=c('idx'='knot')) %>% 
    group_by(year) %>% sample_n(15000,replace=T)
  
  if(any(years!='all')){
    sp_df <- pts %>% left_join(df,by=c('idx'='knot')) %>% 
      filter(year %in% years)
  }
  # make the facetted plot by year
  # cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  out<-ggplot()+
    geom_point(data=sp_df,aes(E_km,N_km,col=dens),shape=".")+
    scale_color_gradientn(colors=cols)+
    facet_wrap(~year)+
    labs(title="",x="Eastings",y="Northings",col="Log Numbers")+
    geom_sf(data=ak,fill='gray50',col=NA)+
    coord_sf(crs = st_crs(ak), datum = NA) +
    xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
    # scale_color_viridis()+
    theme(axis.text = element_blank(),
          panel.spacing.x = unit(0,"pt"),
          panel.spacing.y=unit(0,"pt"))

  return(out)
}

# Plot factor maps and factor loadings by category (species)
plot_fct_loadings <- function(factors,dat,rotated=TRUE) {
  
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  # collect factor loadings for the different predictors
  if(rotated){
    loadings <- factors$Rotated_loadings
  }
  else{
    loadings <- factors$Loadings
  }
  
  loadings_lng <- map_df(names(loadings),.f = function(x){
    m <- loadings[[x]]
    if(is.matrix(m)){
      tibble(load=as.numeric(m),spp=rep(cats,ncol(m)),fct=x,fct_num=rep(1:ncol(m),each=nrow(m)))
    }
  })
  # # proportion of variance explained
  # round(100*sum(L_pj[,whichfactor]^2)/sum(L_pj^2),1)
  fct_labeller <- c(
    "Omega1" = "Spatial Variation\nEncounter",
    "Omega2" = "Spatial Variation\nPositive Abundance",
    "Epsilon1" = "Spatio-temporal Variation\nEncounter",
    "Epsilon2" = "Spatio-temporal Variation\nPositive Abundance"
  )
  # indicator for pos/neg
  loadings_lng <- loadings_lng %>% 
    mutate(fct=factor(fct,levels=c("Omega1","Omega2","Epsilon1","Epsilon2"))) %>% 
    # proportion of variance explained
    group_by(fct,fct_num) %>% 
    mutate(loadsq=sum(load^2)) %>%
    ungroup() %>% group_by(fct) %>% 
    mutate(tot=sum(load^2),prop=round(loadsq/tot*100,1)) %>% 
    ungroup() %>% 
    mutate(is_pos=ifelse(load>0,"yes","no"))
  pl <- rev(pal_npg()(2))
  out <-ggplot(loadings_lng,aes(spp,load,fill=is_pos))+
    geom_bar(stat='identity')+
    geom_hline(yintercept=0,col='black')+
    geom_text(aes(x=4.5,y=1.5,label=prop),check_overlap = T)+
    scale_fill_manual(values=pl)+
    facet_grid(fct_num~fct,labeller = labeller(fct=fct_labeller))+
    guides(fill='none')+
    labs(x="Class (Species/Size Combination)",y="Loading",title="")+
    theme(axis.text.x = element_text(angle=60,hjust=1))
  out
}
# plot 2-dimensional loadings
plot_2d_fct_loadings <- function(dat,factors,rotated=TRUE,which_lpred,tit=""){
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  # collect factor loadings for the different predictors
  if(rotated){
    loadings <- factors$Rotated_loadings
  }
  else{
    loadings <- factors$Loadings
  }
  
  loadings_lng <- map_df(names(loadings),.f = function(x){
    m <- loadings[[x]]
    if(is.matrix(m)){
      tibble(load=as.numeric(m),spp=rep(cats,ncol(m)),fct=x,fct_num=rep(1:ncol(m),each=nrow(m)))
    }
  })
  # # proportion of variance explained
  # round(100*sum(L_pj[,whichfactor]^2)/sum(L_pj^2),1)
  
  # indicator for pos/neg
  loadings_lng <- loadings_lng %>% 
    # proportion of variance explained
    group_by(fct,fct_num) %>% 
    mutate(loadsq=sum(load^2)) %>%
    ungroup() %>% group_by(fct) %>% 
    mutate(tot=sum(load^2),prop=round(loadsq/tot*100,1)) %>% 
    ungroup() %>% 
    mutate(is_pos=ifelse(load>0,"yes","no"))

  pca <- loadings_lng %>% 
    filter(fct_num < 3) %>% 
    select(spp,fct,fct_num,load) %>%
    unite("fct",fct,fct_num) %>% 
    group_by(spp) %>% 
    spread(fct,load) %>% 
    mutate(spp_sym= ifelse(spp %in% c('Opilio Immature','Opilio Spawner'),'crab','cod')) %>% 
    select(spp,spp_sym,contains(which_lpred)) %>% 
    ungroup()
  xmax <- ceiling(max(pca[,3]))
  xmin <- floor(min(pca[,3]))
  ymax <- ceiling(max(pca[,4]))
  ymin <- floor(min(pca[,4]))
  
  varnames <- names(pca)
  
  pl <- pal_npg()(5)
  pca_plot <- pca %>% 
    ggplot(aes_string(varnames[3],varnames[4],col="spp",shape="spp"))+
    geom_point(size=5)+
    coord_equal()+
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    xlim(xmin,xmax)+ylim(ymin,ymax)+
    scale_color_manual(values=pl,breaks = c("Opilio Immature","Opilio Spawner","Small Cod","Medium Cod","Large Cod"))+
    scale_shape_manual(values=c(17,17,16,16,16),breaks = c("Opilio Immature","Opilio Spawner","Small Cod","Medium Cod","Large Cod"))+
    # guides(shape='none')+
    labs(title=tit,x="Factor 1",y="Factor 2",col='',shape='')+
    theme(legend.text = element_text(size=8),
          plot.title = element_text(size=12),
          legend.position = c(0.15,0.25))
  pca_plot
}
plot_fct_maps <- function(Extrapolation_List,Report,dat,which_lpred,which_f,average=FALSE){
  
  pts <- map_points(Extrapolation_List)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$year),max(dat$year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$year)))]
  
  #pull out factor map data, for enc prob and pos values
  fcts <- factors$Rotated_factors %>% pluck(which_lpred)
  n_knots <- dim(fcts)[1]
  fct <- fcts[,which_f,]
  n_knots<-ifelse(grepl('Omega',which_lpred),length(fct),dim(fct)[1])
  n_yrs <- ifelse(grepl('Omega',which_lpred),1,dim(fct)[2])

  df <- tibble(lpred=which_lpred,which_f=which_f,knot=rep(1:n_knots,n_yrs),value=as.numeric(fct))
  if(grepl('Epsilon',which_lpred)) {df <- df %>% mutate(year=rep(Years2Include,each=n_knots))}
  
  # cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  sp_df<-pts %>% left_join(df,by=c('idx'='knot')) %>% group_by(idx) %>% sample_frac(0.5)
  
  if(average){
    sp_df_mean <- pts %>% left_join(df,by=c('idx'='knot')) %>% group_by(idx,E_km,N_km) %>%
      filter(!is.na(value)) %>% 
      summarise(meanval=mean(value)) %>% 
      ungroup() %>% 
      group_by(idx) %>% sample_frac(0.5)
  }
  
  if(grepl('Omega',which_lpred)){
    out<-ggplot()+
      geom_point(data=sp_df,aes(E_km,N_km,col=value),size=0.75)+
      # geom_contour(data = bering_bathy,
      #              aes(x=x, y=y, z=z),
      #              breaks=c(-100),
      #              size=c(0.3),
      #              colour="black")+
      # geom_contour(data = bering_bathy, 
      #              aes(x=x, y=y, z=z),
      #              breaks=c(-50),
      #              size=c(0.3),
      #              colour="black")+
      scale_color_gradient2(low="darkblue",mid="lightgreen",high="red",midpoint=0)+
      labs(title="",x="Eastings",y="Northings",col="")+
      geom_sf(data=ak,fill='gray50',col=NA)+
      coord_sf(crs = st_crs(ak), datum = NA) +
      xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
      theme(axis.text = element_blank(),
            panel.spacing.x = unit(0,"pt"),
            panel.spacing.y=unit(0,"pt"))
  }
  if(grepl('Epsilon',which_lpred)){
    if(average){
      out<-ggplot()+
        geom_point(data=sp_df_mean,aes(E_km,N_km,col=meanval))+
        scale_color_gradient2(low="darkblue",mid="lightgreen",high="red",midpoint=0)+
        labs(title="",x="Eastings",y="Northings",col="")+
        geom_sf(data=ak,fill='gray50',col=NA)+
        coord_sf(crs = st_crs(ak), datum = NA) +
        xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
        theme(axis.text = element_blank(),
              panel.spacing.x = unit(0,"pt"),
              panel.spacing.y=unit(0,"pt"))
    }
    else{
      out<-ggplot()+
        geom_point(data=sp_df,aes(E_km,N_km,col=value),shape=".")+
        scale_color_gradient2(low="darkblue",mid="lightgreen",high="red",midpoint=0)+
        facet_wrap(~year)+
        labs(title="",x="Eastings",y="Northings",col="")+
        geom_sf(data=ak,fill='gray50',col=NA)+
        coord_sf(crs = st_crs(ak), datum = NA) +
        xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
        theme(axis.text = element_blank(),
              panel.spacing.x = unit(0,"pt"),
              panel.spacing.y=unit(0,"pt"))  
    }
  }
  # make the facetted plot by year

    out
}

# Plot species correlations (density correlations)
cor.sig <- function(x,y){
  corr <- cor(x,y,method='pearson')
  test <- suppressWarnings(cor.test(x,y,method='pearson'))
  out <- ifelse(test$p.value<0.05,corr,NA)
  out
}

plot_category_correlations <- function(Report, Extrapolation_List,dat,type="density"){
  
  if(!(type %in% c('density','enc','pos'))){
    stop("type parameter must be one of 'density','enc','pos'")
  }
  
  if(type=='density'){
    est <- log(Report$D_xcy) 
    lab <- "Density"
  }
  if(type=='enc'){
    est <- Report$R1_xcy
    lab <- "Encounter Probability"
  }
  if(type=='pos'){
    est <- Report$R2_xcy
    lab <- "Positive Catch Rate"
  }
  
  # row bind all the years for each species
  report_dims <- dim(Report$D_xcy)
  dim(est) <- c(report_dims[1]*report_dims[3],report_dims[2])
  
  # get correlations
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  spp_cor <- numeric()
  
  for(i in 1:length(cats)){
    for(j in 1:length(cats)){
      spp_cor[(length(cats)*(i-1))+j] = cor.sig(est[,i],est[,j])
    }
  }
  spp_cor[spp_cor==1] <- NA
  
  df <- tibble(spp1=rep(cats,each=length(cats)),spp2=rep(cats,length(cats)),corr=spp_cor)
  
  out <- df %>% 
    ggplot(aes(spp1,spp2,fill=corr))+
    geom_tile()+
    scale_fill_gradient2(low='darkblue',high='darkred',mid='white',limits=c(-0.2,0.1),na.value='gray80')+
    geom_text(aes(label=round(corr,2)))+
    coord_equal()+
    guides(fill='none')+
    labs(x="",y="",title="")+
    theme(axis.text = element_text(angle=60,hjust=1),
          legend.title = element_text(size=14),
          panel.border = element_blank())
  out
}

map_category_correlations <- function(Report, Extrapolation_List,dat,type="density",species_1,species_2,scl.lims=c(-1,1)){
  
  if(!(type %in% c('density','enc','pos'))){
    stop("type parameter must be one of 'density','enc','pos'")
  }
  
  if(type=='density'){
    est <- log(Report$D_xcy) 
    lab <- "Density"
  }
  if(type=='enc'){
    est <- Report$R1_xcy
    lab <- "Encounter Probability"
  }
  if(type=='pos'){
    est <- Report$R2_xcy
    lab <- "Positive Catch Rate"
  }
  
  # row bind all the years for each species
  # report_dims <- dim(Report$D_xcy)
  # dim(est) <- c(report_dims[1]*report_dims[3],report_dims[2])
  
  # get correlations for each knot
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  # spp_cor <- numeric()
  
  spp_corr_by_knot <- array(dim=c(dim(est)[1],dim(est)[2],dim(est)[2]))
  for(k in 1:dim(est)[1]){
    corrmat <- numeric()
    for(i in 1:length(cats)){
      for(j in 1:length(cats)){
        corrmat[(length(cats)*(i-1))+j] <- cor.sig(est[k,i,],est[k,j,])
      }
    }
    spp_corr_by_knot[k,,] <- corrmat
  }
  knots = nrow(Report$D_xcy)
  df <- tibble(knot=rep(1:knots,length(cats)*length(cats)),spp1=rep(rep(cats,each=knots),length(cats)),spp2=rep(cats,each=knots*length(cats)),corr=as.numeric(spp_corr_by_knot))
  
  df_cat <- df %>% 
    filter(spp1==species_1,spp2==species_2)
  
  pts <- map_points(Extrapolation_List)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  
  sp_df<-pts %>% left_join(df_cat,by=c('idx'='knot')) %>% group_by(idx) %>% sample_frac(0.5) %>% ungroup()
  
  out<-ggplot()+
    geom_point(data=sp_df,aes(E_km,N_km,col=corr))+
    scale_color_gradient2(low='darkblue',mid='lightgreen',high='darkred',midpoint=0,na.value = 'gray80',limits=scl.lims)+
    labs(title="",x="Eastings",y="Northings",col="")+
    geom_sf(data=ak,fill='gray50',col=NA)+
    coord_sf(crs = st_crs(ak), datum = NA) +
    xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
    # scale_color_viridis()+
    theme(axis.text = element_blank(),
          panel.spacing.x = unit(0,"pt"),
          panel.spacing.y=unit(0,"pt"),
          legend.position = c(0.1,0.25),
          legend.key.size = unit(0.5,'cm'),
          legend.text = element_text(size=10))
  out
}

plot_abundance <- function(Index_table,dat){
  abun <- Index_table %>% as.data.frame() %>% 
    mutate(log_abun=log(Estimate_metric_tons),upper=log_abun+SD_log,lower=log_abun-SD_log) %>% 
    mutate(Category=tools::toTitleCase(as.character(Category))) %>% 
    mutate(Category=factor(Category,levels=c("Opilio Immature","Opilio Spawner","Small Cod","Medium Cod","Large Cod")))
  
  ggplot(abun,aes(Year,log_abun,ymax=upper,ymin=lower,fill=Category))+
    geom_ribbon()+
    geom_line()+
    scale_fill_npg()+
    guides(fill='none')+
    facet_wrap(~Category,ncol=1,scales="free_y")+
    labs(x="Year",y="Log Abundance")

}

plot_cog <- function(Report,dat){
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  cog_arr <- Report$mean_Z_cym
  ncats <- dim(cog_arr)[1]
  nyrs <- dim(cog_arr)[2]
  years <- seq(min(dat$year),max(dat$year))
  dim(cog_arr) <- c(ncats*nyrs,dim(cog_arr)[3])
  
  cog_tbl <- tibble(category=rep(cats,nyrs),year=rep(years,each=ncats),East=cog_arr[,1],North=cog_arr[,2]) %>% 
    gather(direction,value,East,North)
  
  ggplot(cog_tbl,aes(year,value,col=category))+
    geom_line(size=1)+
    scale_color_npg()+
    guides(color='none')+
    facet_grid(direction~category,scales="free")+
    labs(x="Year",y="Kilometers")+
    theme(strip.text=element_text(size=10))
}

plot_2d_cog <- function(Report,dat) {
  
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  cog_arr <- Report$mean_Z_cym
  ncats <- dim(cog_arr)[1]
  nyrs <- dim(cog_arr)[2]
  years <- seq(min(dat$year),max(dat$year))
  dim(cog_arr) <- c(ncats*nyrs,dim(cog_arr)[3])
  
  cog_tbl <- tibble(category=rep(cats,nyrs),year=rep(years,each=ncats),east=cog_arr[,1],north=cog_arr[,2]) %>% 
    group_by(category) %>% 
    mutate(leadeast=lead(east,1),leadnorth=lead(north,1))
  out <- ggplot(cog_tbl,aes(east,north,col=year,xend=leadeast,yend=leadnorth))+
    geom_point()+
    geom_segment(arrow=arrow(angle=20,length = unit(0.2, "cm"),type='closed'))+
    labs(col='Year',x="Eastings",y="Northings")+
    coord_equal()+
    scale_color_viridis_c()+
    scale_size_continuous(range=c(1,4))+
    facet_wrap(~category)
  out
}

# Factor correlations with temperature

find_density_covariates <- function(Spatial_List,dat) {
  # locations of knots from Spatial_List object
  loc_query <- Spatial_List[['loc_x']]
  # years
  yrs <- seq(min(dat$year),max(dat$year))
  # locations of interpolated temperature and depth info
  loc_covars <- nbt.interpolated %>% distinct(x,y)
  # nearest neighbors (knots and interpolated covariates)
  which_nn <- nn2(data=loc_covars,query=loc_query,k=1)$nn.idx[,1]
  # for each knot in each year, find nearest neighbor from the interpolated datasets
  X_xtp <- array(data=NA,dim=c(nrow(loc_query),length(yrs),2))
  # fill in covariates for each knot in each year
  # depth is constant across years
  for(j in 1:length(yrs)) {
    yr_temp <- nbt.interpolated %>% filter(year==yrs[j])
    X_xtp[,j,1] <- yr_temp$temp[which_nn]
    X_xtp[,j,2] <- depth.interpolated$depth[which_nn]
  }
  return(X_xtp)
}

map_temperature <- function(Spatial_List,dat,years='all') {
  
  covars <- find_density_covariates(Spatial_List,dat)
  cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  pts <- map_points(Extrapolation_List)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$year),max(dat$year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$year)))]
  knots <- dim(covars)[1]
  
  cats <- tools::toTitleCase(levels(dat$spp))
  
  df <- tibble(knot=rep(1:knots,length(Years2Include)),year=rep(Years2Include,each=knots),meantemp=as.numeric(covars[,,1]))
  sp_df<-pts %>% left_join(df,by=c('idx'='knot')) %>% 
    group_by(year) %>% sample_n(15000,replace=T)
  
  if(any(years!='all')){
    sp_df <- pts %>% left_join(df,by=c('idx'='knot')) %>% 
      filter(year %in% years)
  }
  # make the facetted plot by year
  # cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  out<-ggplot()+
    geom_point(data=sp_df,aes(E_km,N_km,col=meantemp),shape=".")+
    scale_color_gradientn(colors=cols)+
    facet_wrap(~year)+
    labs(title="",x="Eastings",y="Northings",col="Temperature")+
    geom_sf(data=ak,fill='gray50',col=NA)+
    coord_sf(crs = st_crs(ak), datum = NA) +
    xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
    # scale_color_viridis()+
    theme(axis.text = element_blank(),
          panel.spacing.x = unit(0,"pt"),
          panel.spacing.y=unit(0,"pt"))
  
  return(out)
}


fct_covariate_correlation <- function(Spatial_List,dat,factors,which_lpred,which_f,average=FALSE){
  covars <- find_density_covariates(Spatial_List,dat)
  knots <- 1:dim(covars)[1]
  years <- seq(min(dat$year),max(dat$year))
  dim(covars) <- c(dim(covars)[1]*dim(covars)[2],2)
  covars <- as.data.frame(covars) %>% set_names("temp","depth") %>% mutate(knot=rep(knots,length(years)),year=rep(years,each=length(knots)))

  fcts <- factors$Rotated_factors %>% pluck(which_lpred)
  n_knots <- dim(fcts)[1]
  fct <- fcts[,which_f,]
  n_knots<-ifelse(grepl('Omega',which_lpred),length(fct),dim(fct)[1])
  n_yrs <- ifelse(grepl('Omega',which_lpred),1,dim(fct)[2])
  
  if(grepl('Epsilon',which_lpred)) {
    # if looking at spatio-temporal correlations, use temperature anomalies
    covars_anom <- covars %>%
      select(knot,year,temp) %>% 
      group_by(knot) %>% 
      mutate(temp_anom=temp-mean(temp)) %>% 
      ungroup()
    df <- tibble(lpred=which_lpred,which_f=which_f,knot=rep(1:n_knots,n_yrs),value=as.numeric(fct)) %>% 
      mutate(year=rep(Years2Include,each=n_knots)) %>% 
      left_join(covars_anom,by=c('knot','year')) %>% 
      filter(!is.na(temp_anom))
    if(average){
      df <- df %>% group_by(knot,lpred,which_f) %>% summarize(meanval=mean(value),meananom=mean(temp_anom)) %>% ungroup()
      out <- df %>% 
        group_by(lpred,which_f) %>% 
        summarize(temp_anom_corr=cor.sig(meanval,meananom))
    }
    else{
      out <- df %>% 
        group_by(lpred,which_f,year) %>% 
        summarize(temp_anom_corr=cor.sig(value,temp_anom))
    }
  }
  if(grepl('Omega',which_lpred)) {
    df <- tibble(lpred=which_lpred,which_f=which_f,knot=rep(1:n_knots,n_yrs),value=as.numeric(fct)) %>% 
      left_join(covars,by=c('knot'))%>% 
      filter(!is.na(temp),!is.na(depth)) %>% 
      distinct()
    out <- df %>% 
      group_by(lpred,which_f) %>% 
      summarize(tempcorr=cor.sig(value,temp),depthcorr=cor.sig(value,depth))
  }
  out
}

abun_covariate_correlation <- function(Spatial_List,Extrapolation_List,dat,Report,category,map_result=FALSE,scl.lims=c(-1,1)){
  covars <- find_density_covariates(Spatial_List,dat)
  knots <- 1:dim(covars)[1]
  years <- seq(min(dat$year),max(dat$year))
  dim(covars) <- c(dim(covars)[1]*dim(covars)[2],2)
  covars <- as.data.frame(covars) %>% set_names("temp","depth") %>% mutate(knot=rep(knots,length(years)),year=rep(years,each=length(knots)))
  
  covars_anom <- covars %>%
    select(knot,year,temp) %>% 
    group_by(knot) %>% 
    mutate(temp_anom=temp-mean(temp)) %>% 
    ungroup()
  
  dens <- Report$D_xcy
  Year_Set = seq(min(dat$year),max(dat$year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$year)))]
  knots <- dim(dens)[1]
  
  cats <- tools::toTitleCase(levels(dat$spp))
  
  out <- tibble(knot=rep(1:knots,length(Years2Include)),year=rep(Years2Include,each=knots),dens=log(as.numeric(dens[,category,]))) %>% 
    left_join(covars_anom,by=c('knot','year')) %>% 
    group_by(knot) %>% 
    summarize(temp_anom_corr=cor.sig(dens,temp_anom))
  
  if(map_result){
    pts <- map_points(Extrapolation_List)
    lims <- st_bbox(pts)
    st_geometry(pts) <- NULL
    Year_Set = seq(min(dat$year),max(dat$year))
    Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$year)))]
    knots <- dim(dens)[1]
    
    cats <- tools::toTitleCase(levels(dat$spp))
    
    sp_df<-pts %>% left_join(out,by=c('idx'='knot')) %>% 
      group_by(idx) %>% sample_frac(0.5)
    
    col_grad <- function(sp_df){
      if(max(sp_df$temp_anom_corr,na.rm=T)>0&min(sp_df$temp_anom_corr,na.rm = T)<0){
        scale_color_gradient2(low='darkblue',mid='lightgreen',high='darkred',na.value='gray80',midpoint = 0,limits=scl.lims)
      }
      else{
        scale_color_gradient(low='darkblue',high='lightgreen',na.value='gray80',limits=scl.lims)
      }
    }
    out_map<-ggplot()+
      geom_point(data=sp_df,aes(E_km,N_km,col=temp_anom_corr))+
      col_grad(sp_df)+
      labs(title="",x="Eastings",y="Northings",col="")+
      geom_sf(data=ak,fill='gray50',col=NA)+
      coord_sf(crs = st_crs(ak), datum = NA) +
      xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
      labs(x="Eastings",y="Northings",col="")+
      theme(axis.text = element_blank(),
            panel.spacing.x = unit(0,"pt"),
            panel.spacing.y=unit(0,"pt"),
            legend.position = c(0.1,0.25),
            legend.key.size=unit(0.5,'cm'),
            legend.text=element_text(size=10))
    out_map
  }
  else{
    out
  }
}
