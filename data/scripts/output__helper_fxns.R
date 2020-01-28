### Helper functions to visualize VAST outputs ####
### specific to Eastern Bering Sea ###
## these are mostly to make nicer versions of the standard plots ##

library(tidyverse)
library(sf)
library(VAST)
library(abind)
library(viridis)
library(ggsci)
library(here)

# basemap
# library(rnaturalearth)
# library(rnaturalearthdata)
# ak <- ne_states(country='United States of America',geounit="Alaska",returnclass = 'sf')

ak <- read_sf(here::here('data','spatial','cb_2017_02_anrc_500k.shp')) %>% 
  st_union() %>% 
  st_transform(26904)

## data for plotting
fp = here::here('data','VAST output','plots') %>% paste0("/")
load(paste0(fp,"derived_quantities.Rdata"))

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

# Plot predicted density
plot_density <- function(Region,Spatial_List, Report,Extrapolation_List,dat,fp,saveplots=TRUE) {
  
  dens <- log(Report$D_gcy)
  cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  pts <- Extrapolation_List$Data_Extrap[,c('Lon','Lat','Area_in_survey_km2')] %>% 
    st_as_sf(coords=c('Lon','Lat'),crs=4326) %>% 
    #convert to AK UTM zone
    st_transform(26904)
  pts[,c("E_km","N_km")] <- st_coordinates(pts)
  pts$idx <-as.integer(Spatial_List$PolygonList$NN_Extrap$nn.idx)
  pts <- filter(pts,Area_in_survey_km2>0) %>% select(-Area_in_survey_km2)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$Year),max(dat$Year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$Year)))]
  knots <- dim(dens)[1]
  
  cats <- tools::toTitleCase(levels(dat$spp))
  
  plotlist <- list()
  for(i in 1:dim(dens)[2]){
    df <- tibble(knot=rep(1:knots,length(Years2Include)),year=rep(Years2Include,each=knots),dens=as.numeric(dens[,i,]))
    sp_df<-pts %>% left_join(df,by=c('idx'='knot'))
    # make the facetted plot by year
    cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
    out<-ggplot()+
      geom_point(data=sp_df,aes(E_km,N_km,col=dens),size=0.25)+
      scale_color_gradientn(colors=cols)+
      facet_wrap(~year)+
      labs(title=cats[i],x="Eastings",y="Northings",col="")+
      geom_sf(data=ak,fill='gray50',col=NA)+
      coord_sf(crs = st_crs(ak), datum = NA) +
      xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
      theme(axis.text = element_blank(),
            panel.spacing.x = unit(0,"pt"),
            panel.spacing.y=unit(0,"pt"))
    # save the plot
    if(saveplots) ggsave(plot=out,filename = paste0(fp,"dens_",cats[i],".png"),w=6,h=6)
    plotlist[[i]] <- out
  }
  return(plotlist)
}

# Plot factor maps and factor loadings by category (species)
plot_fct_loadings <- function(factors,dat,rotated=TRUE,fp,saveplots=TRUE) {
  
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
  out <-ggplot(loadings_lng,aes(spp,load,fill=is_pos))+
    geom_bar(stat='identity')+
    geom_hline(yintercept=0,col='black')+
    geom_text(aes(x=4,y=1.5,label=prop),check_overlap = T)+
    facet_grid(fct_num~fct)+
    guides(fill='none')+
    labs(x="Species",y="Loading",title="Factor Loadings")+
    theme(axis.text.x = element_text(angle=60,hjust=1))
  out
  if(saveplots){
    r <- ifelse(rotated,"rotated","unrotated")
    ggsave(plot=out,filename = paste0(fp,"factor_loadings_",r,".png"),w=6,h=6)
  }
  return(loadings_lng)
}

plot_fct_maps <- function(Region,MapDetails_List, Report,dat,fp,saveplots=TRUE){
  
  #pull out factor map data, for enc prob and pos values
  fcts <- factors$Rotated_factors
  epsilon <-fcts[c("Epsilon1","Epsilon2")]
  n_factors <-   dim(epsilon[[1]])[2]
  
  df<-map_df(names(epsilon),.f=function(x) {
    arr <- epsilon[[x]]
    n_knots<-dim(arr)[1]
    n_yrs <- dim(arr)[3]
    dim(arr) <- c(n_knots*n_yrs,dim(arr)[2])
    arr <- as.data.frame(arr)
    names(arr) <- paste0(c("Factor_"),1:ncol(arr))
    arr$knot=rep(1:n_knots,each=n_yrs)
    arr$year=rep(Years2Include,n_knots)
    arr$lpred <- x
    arr
  })

  cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  pts <- MapDetails_List$PlotDF[,c('Lon','Lat','x2i','Include')] %>% 
    st_as_sf(coords=c('Lon','Lat'),crs=4326) %>% 
    #convert to AK UTM zone
    st_transform(26904)
  pts[,c("E_km","N_km")] <- st_coordinates(pts)
  pts$idx <-as.integer(pts$x2i)
  pts <- filter(pts,Include>0) %>% select(-x2i,-Include)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$Year),max(dat$Year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$Year)))]
  knots <- unique(df$knot)
  
  sp_df<-pts %>% left_join(df,by=c('idx'='knot'))
  
  plotlist <- list()
  for(i in 1:n_factors){
    sp_dfx <- sp_df %>% dplyr::select(E_km,N_km,year,lpred,3+i) %>% 
    fctnam <- tail(names(sp_dfx),1)
    sp_dfx$value <- sp_dfx[,fctnam]
    # make the facetted plot by year
    out<-ggplot()+
      geom_point(data=sp_dfx,aes(E_km,N_km,col=value))+
      scale_color_gradientn(colors=cols)+
      facet_wrap(~year)+
      labs(title=cats[i],x="Eastings",y="Northings",col="")+
      geom_sf(data=ak,fill='gray50',col=NA)+
      coord_sf(crs = st_crs(ak), datum = NA) +
      xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
      theme(axis.text = element_blank(),
            panel.spacing.x = unit(0,"pt"),
            panel.spacing.y=unit(0,"pt"))
    # save the plot
    if(saveplots) ggsave(plot=out,filename = paste0(fp,"dens_",cats[i],".png"),w=6,h=6)
    plotlist[[i]] <- out
  }
  return(plotlist)
}
# Plot species correlations (density correlations)
plot_category_correlations <- function(Report,Spatial_List, Extrapolation_List,dat,type="density"){
  
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
  cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  # row bind all the years for each species
  report_dims <- dim(Report$D_xcy)
  dim(est) <- c(report_dims[1]*report_dims[3],report_dims[2])
  
  # get correlations
  spp_cor <- cor(est,est)
  
  # make graph
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  df <- tibble(corr=as.numeric(spp_cor),spp1=rep(cats,each=length(cats)),spp2=rep(cats,length(cats))) %>% 
    mutate(corr=ifelse(round(corr,2)==1.00,NA,corr))
  
  df %>% 
    ggplot(aes(spp1,spp2,fill=corr))+
    geom_tile()+
    scale_fill_gradient2()+
    geom_text(aes(label=round(corr,2)))+
    coord_equal()+
    labs(x="",y="",title=paste("Correlation in Predicted",lab),fill=expression(rho))+
    theme(axis.text.x = element_text(angle=60,hjust=1),
          legend.title = element_text(size=14))

}

plot_2d_cog <- function(Spatial_List,Report,Sdreport,dat,dens.tbl,use_biascorr=F) {
  
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  CogName <- 'mean_Z_cym'
  
  cog_arr <- Save$Report$mean_Z_cym
  ncats <- dim(cog_arr)[1]
  nyrs <- dim(cog_arr)[2]
  years <- seq(min(dat$Year),max(dat$Year))
  dim(cog_arr) <- c(ncats*nyrs,dim(cog_arr)[3])
  cog_tbl <- tibble(category=rep(cats,nyrs),year=rep(years,each=ncats),east=cog_arr[,1],north=cog_arr[,2])
  
  # index
  dens_tbl <- dens.tbl %>% rename(category=Category,year=Year,est_mt=Estimate_metric_tons) %>% 
    select(category,year,est_mt) %>% 
    mutate(category=tools::toTitleCase(category))
  
  cog_dens <- left_join(cog_tbl,dens_tbl,by=c('category','year'))
  
  out <- ggplot(cog_dens,aes(east,north,col=year,size=est_mt))+
    geom_point()+
    coord_equal()+
    facet_wrap(~category)
  out
}
