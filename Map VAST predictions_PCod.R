##############################################################################################################################
## 
## This script maps VAST output for the PCod diet study
## 
##############################################################################################################################
rm( list = ls() )

######## Define directories
InputFiles = 'E:/Diet_composition_study/Input'
OutputFiles = 'E:/Diet_composition_study/Output'
MyMaps = 'E:/Diet_composition_study/Maps_Pcod'

######## Load required R libraries
library( VAST )

######## Load information and define some settings
setwd( InputFiles )
load( "Save_orig_1.RData" )
Report_1 = Save_orig$Report
load( "Save_orig_2.RData" )
Report_2 = Save_orig$Report
Y_gs_1 = log(Report_1$Index_xcyl + Report_2$Index_xcyl + 1e-15)
Year_Set = 1984:2015
Years2Include = c( 1:20, 22:32 )
Y_gs_2 = Y_gs_1[,2,Years2Include,1,drop=TRUE]
load( "MapDetails_List_PCod.RData" )
MapDetails_List <- MapDetails_List
PlotDF <- MapDetails_List$PlotDF
Y_gs_3 = Y_gs_2[PlotDF[,'x2i'],,drop=FALSE]
Y_gs = Y_gs_3[PlotDF[,'Include']>0,] 
setwd( OutputFiles )
Log_PESC_Extrapolation_grid_cells <- Y_gs
save( Log_PESC_Extrapolation_grid_cells, file = "Log_PESC_Extrapolation_grid_cells_PCod.RData" )
category_names = 1

######## Produce a map for 1984
tI <- 1
file_name = "Log_PESC_Pacific_cod_1984"
Y_gt = Y_gs[,tI]
zlim = range( Y_gs, na.rm = TRUE )
map_list = MapDetails_List
panel_labels = category_names
projargs = '+proj=longlat'
map_resolution = "medium"
setwd( MyMaps )
working_dir = paste0( getwd(),"/" )
Format = "png" 
Res = 300
add = TRUE
outermargintext = c( "          Longitude (°E)", "          Latitude (°N)" )
mar = c( 0, 0, 2, 0 )
oma = c( 4, 4, 0, 0 )
land_color = "grey48"
cex.legend = 1.25
country = NULL
Y_gt = matrix( Y_gt, ncol = 1 )
MapSizeRatio = map_list$MapSizeRatio
n_cells = nrow( Y_gt )
mfrow = ceiling( sqrt( ncol( Y_gt ) ) )
mfrow = c( mfrow, ceiling( ncol( Y_gt ) / mfrow ) )
col = colorRampPalette( colors = c( "darkblue", "blue", "lightblue", "lightgreen", "yellow", "orange", "red" ) )
col = col( 1000 )
legend_x = c( 0.05, 0.10 )
legend_y = c( 0.04, 0.34 )
loc_g = map_list$PlotDF[which( map_list$PlotDF[,'Include'] > 0 ), c( 'Lon', 'Lat' )]
CRS_orig = sp::CRS( '+proj=longlat' )
CRS_proj = sp::CRS( projargs )
map_data = rnaturalearth::ne_countries( scale = switch( map_resolution, "low" = 110, "medium" = 50, 
	"high" = 10, 50 ), country = country )
map_data = sp::spTransform( map_data, CRSobj = CRS_proj )
Par = list( mfrow = mfrow, mar = mar, oma = oma )
png( file = paste0( working_dir, file_name, ".png" ), width = Par$mfrow[2]*MapSizeRatio[2],
	height = Par$mfrow[1]*MapSizeRatio[1], res = Res, units = 'in')
		#par( Par )
		par( mfrow = c( 1, 1 ) )
		par( mar = c( 4.5, 4.5, 0.5, 0.5 ) )
		Points_orig = sp::SpatialPointsDataFrame( coords = loc_g, data = data.frame( "y" = Y_gt[,1] ), 
			proj4string = CRS_orig )
		Points_LongLat = sp::spTransform( Points_orig, sp::CRS( '+proj=longlat' ) )
		Points_proj = sp::spTransform( Points_orig, CRS_proj )
		cell.size = mean( diff( Points_proj@bbox[1,] ), diff( Points_proj@bbox[2,] ) ) / floor( sqrt( n_cells ) )
    		Raster_proj = plotKML::vect2rast( Points_proj, cell.size = cell.size )
		xlim = Raster_proj@bbox[1,]
		ylim = Raster_proj@bbox[2,]
    		image( Raster_proj, col = col, zlim = zlim, xlim = xlim, ylim = ylim )
		sp::plot( map_data, col = land_color, add = TRUE )
		box()
    		xl = ( 1 - legend_x[1] ) * par( 'usr' )[1] + ( legend_x[1] ) * par( 'usr' )[2]
    		xr = ( 1 - legend_x[2] ) * par( 'usr' )[1] + ( legend_x[2] ) * par( 'usr' )[2]
    		yb = ( 1 - legend_y[1] ) * par( 'usr' )[3] + ( legend_y[1] ) * par( 'usr' )[4]
    		yt = ( 1 - legend_y[2] ) * par( 'usr' )[3] + ( legend_y[2] ) * par( 'usr' )[4]
		align = c( "lt","rb" )[2]
     	 	gradient = c( "x", "y" )[2]
		plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
			legend = round( seq( zlim[1], zlim[2], length = 4 ), 1 ), 
			rect.col = col, cex = cex.legend, align = align, gradient = gradient )
		#axis( 1, c( 121, 123, 125 ), c( "121", "123", "125" ) )
		axis( 2 )
		#mtext( side = 1, outer = TRUE, outermargintext[1], cex = 1.25, line = par()$oma[1]/2, padj = -1.8 )
		mtext( side = 2, outer = TRUE, outermargintext[2], cex = 1.25, line = par()$oma[2]/2, padj = 1.8 )
dev.off() 

