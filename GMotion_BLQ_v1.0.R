#
#------- GMotion -----------------------------------
#
# Version:      1.0
# Project Name: RMC-Alliance
# Author:       i.marzan@csic.es
# Date:         June 2023
# License:      CC BY-NC
# Description:  This R script is dedicated to postprocess the ground motion data from the Copernicus EGMS platform.
#               The code reads the EGMS CSV files, subsets the data according to specified region of interest, and
#               generates interactive maps and plots to conduct timeseries analysis.
#               The main products are timeseries multilayer maps and timeseries cross-sections.
# 
# Comment:      Although this code can be applied to any zip file downloaded from the Copernicus EGMS platform, 
#               this code is parameterized for the subsidence experienced in the Bologna area, Italy, file EGMS_L3_E44N23_100km_U.zip
#               It serves as a practical exercise in the context of the advanced course RawMatCop-Alliance, 
#               focusing on Earth Observation applied to Raw Materials. 
#               The EGMS data can be found in the working directory .../Data
#
# cite:         Use the Zenodo doi available in the Github readme
# -----------------------------------------------------
#
# wroking directory
setwd("..../Data") 
getwd()
#
##############
### 1. PACKAGES
###### 

# Installing packages from CRAN respository
if (!require(tidyverse)) {install.packages("tidyverse")}
if (!require(terra)) {install.packages("terra")}
if (!require(raster)) {install.packages("raster")}
if (!require(akima)) {install.packages("akima")}
if (!require(plotly)) {install.packages("htmlwidgets")}
if (!require(leaflet.extras)) {install.packages("leaflet.extras")}   

#' Installing packages from the Github repository because development versions are needed, 
#' as advanced functionalities are required
if (!require(devtools)) {install.packages("devtools")} # devtools allows access to development packages
devtools::install_github("rstudio/leaflet") 
devtools::install_github("rspatial/terra")


# loading installed packages into the session
library(tidyverse)   # For data manipulation and visualization
library(akima)       # For interpolation of irregularly and regularly spaced data
library(raster)      # For spatial raster data manipulation
library(terra)       # Next-generation interface to spatial data in R, based on the raster package
library(leaflet)     # For creating interactive maps within R
library(leaflet.extras)
library(htmlwidgets)      # To transform graphic product into a HTML file


###################
### 2. SUBSETTING
########

# these are the 4 EGMS tiles around Bologna
EGMS_zip = "EGMS_L3_E44N23_100km_U.zip"

# ROI
# The EGMS projection is in EPSG:3035 ETRS89-extended / LAEA Europe / Lambert azimuthal equal-area projection
# limits to the ROI in meters 
W_limit = 4365000 ; E_limit = 4450000
S_limit = 2375000 ; N_limit = 2420000


# Subsetting without unzip 

# csv file inside the zip file
EGMS_csv = gsub(".zip", ".csv", EGMS_zip)

# Subset using R unz
egms = subset(read_csv(unz(EGMS_zip, EGMS_csv)), 
             easting  < E_limit & easting  > W_limit & 
             northing < N_limit & northing > S_limit)
dim(egms)
head(egms, 2)
summary(egms[, 0:12]) # Analysis of twelve first columns because the timeseries begins there.
class(egms)

#' Data distribution
boxplot(egms$mean_velocity) 


#'##############
### 3. DATA INTERPOLATION 
#'##############

# 3.1. Build the grid data to nest the interpolation

# Set the grid spacing in meters
grid_step <- 400  

# Define x and y values of the grid
x_cells <- with(egms, seq(min(easting), max(easting), by = grid_step))
y_cells <- with(egms, seq(min(northing), max(northing), by = grid_step))

# 3.2. Perform spline interpolation on mean_velocity
egms_vel <- with(egms, interp(x = easting, y = northing, z = mean_velocity, 
                 xo = x_cells, yo = y_cells,
                 duplicate = "mean", linear=TRUE)) # if linear=False cubic spline is computed
class(egms_vel)

# 3.3. Convert the interpolation result to a raster
egms_vel_r <- raster(egms_vel) 

   dim(egms_vel_r)    # raster dimension

   plot(egms_vel_r)   # Plot the raster 

   boxplot(egms_vel_r) #' Data distribution

#'
# 3.4. In some areas, interpolation can produce mathematical artifacts (outliers), especially at dataset boundaries,
#'    skewing analysis and visualization. Apply this to identify and remove these extreme values in the raster.
#'
#'    Establish outlier boundaries for mean velocity
      vel_lower_limit = -40
      vel_upper_limit = 40

#     Replace values outside these limits with NA.
      egms_vel_r2 = egms_vel_r          # preserve the original raster, work on the copy
      egms_vel_r2[egms_vel_r2 < vel_lower_limit | egms_vel_r2 > vel_upper_limit] <- NA

#     Plot the raster again
      plot(egms_vel_r2)
      boxplot(egms_vel_r2)

#' 3.5. The interpolation results are not valid in areas with no data. 
#'    This areas need to be masked.  

#    3.5.1 Get egms coordinates file with the LAEA Europe coordinates reference system (CRS). Coordinates are in meters.  
      egms_vel_m <- vect(cbind(egms$easting, egms$northing), crs = "EPSG:3035") # SpatVector object of terra package

#     3.5.2 Convert the velocity raster to SaptRaster object (terra), more efficient       
      egms_vel_rm = rast(egms_vel_r2)
      # Apply LAEA Europe CRS, then coordinates are in meters
      crs(egms_vel_rm) <- "EPSG:3035" 
        
      plot(egms_vel_rm) 
      points(egms_vel_m)

#     3.5.3 Now compute distances between for each raster cell to the near egms point 
      distance_rm <- distance(egms_vel_rm, egms_vel_m)

#     3.5.4 Masking
           egms_vel_rm2 = egms_vel_rm # preserve the original raster, work on the copy

      #    Determine the maximum distance for interpolation confidence
           offset = 700
           
      #    Change all values in egms_vel_rm2 to NA if the corresponding distance in distance_rm is greater than x
           egms_vel_rm2[distance_rm > offset] <- NA # distance in meters, choose the offset.

      plot(egms_vel_rm2) 
      points(egms_vel_m)


#'##############
### 4. PLOTTING MEAN VELOCITY OVER DYNAMIC MAP
#'##############

# 4.1 Project the raster to a geographic reference system, coordinates are in degrees
egms_vel_rd = project(egms_vel_rm2, "EPSG:4326")
      

# 4.2 Color palette 
#  Define a color ramp with the InSAR color convention and scale limits 
   egms_colors = c("darkred", "red", "orange", "yellow", "green","cyan", "dodgerblue", "blue", "darkblue")
   zlim_v <- c(-25, 25) # mm/y velocity limits 
#
#' Create a color palette with colorNumeric function 
#' using the previously defined variables egms_colors and zlim_v 
   mypal_vel <- colorNumeric(palette = egms_colors, # color gradient to use
                          domain = zlim_v,       # range of the numeric variable
                          na.color = NA          # no color for NA values
                          )

   
# 4.3 Create a dynamic map with the lonlat (degrees) raster and a custom color palette
   p_vel <- plet(egms_vel_rd,            # raster file
              alpha=0.7,              # opacity
              tiles="Stamen.Terrain", # basemap provider, more options: https://leaflet-extras.github.io/leaflet-providers/preview/
              col = mypal_vel,        # color palette
              legend = NULL           # no legend
             )  %>% 
         addLegend(
              position = "bottomleft",  # add a legend
              pal = mypal_vel,          # legend colors
              values = zlim_v,          # legend values
              title = "mean-vel mm/y",  # legend title
              opacity = 0.7             # legend opacity
             )  %>% 
         addScaleBar("bottomright"      # add scalebar
             ) %>%              
         addDrawToolbar(                # add toolbar 
                        polylineOptions = TRUE,   # Enable polyline draw
                        polygonOptions      = FALSE, 
                        circleOptions       = FALSE, 
                        circleMarkerOptions = FALSE, 
                        markerOptions       = FALSE, 
                        rectangleOptions    = FALSE,
             )%>% 
         addReverseSearchOSM(           # add location to polyline nodes
                        showSearchLocation = TRUE, 
                        showFeature = FALSE, 
                        fitBounds   = FALSE, 
                        displayText = FALSE
             ) 

p_vel

# Transform the plot into a HTML file
 htmlwidgets::saveWidget(p_vel, "mean_velocity.html")


#'##############
### 5. MAPPING GROUND MOTION TIMESERIES 
#'##############

# 5.1. Define the timeseries monitoring period
#
#   Specifying first, last date columns and the sequential increment
    first_col = 12
    last_col  = ncol(egms)
    leap_col  = 5

    timeseries = colnames(egms)[seq(first_col, last_col, by = leap_col)] 

# 5.2. Creat the function GMOTION to apply linear interpolation to the ground motion timeseries.
#    GMOTION takes the selected monitoring period, performs a linear interpolation for each data column, 
#    and returns a list of rasters.
#
   GMotion <- function(timeseries) {
      # Perform the spline interpolation using the same grid define in section 3
      grd <- with (egms, akima::interp(x = easting, y = northing, z = timeseries, 
                                   xo = x_cells, yo = y_cells,
                                   duplicate = "mean", linear = TRUE)) # if linear=False cubic spline is computed
      r = raster(grd)  
      return(r)
      }

# 5.3. Apply the function to each column and store the results in a list. 
# 
#    Attention! 
#    This operation can take several minutes depending on the parameters: grid_step & leap_col
     egms_ts <- lapply(egms[timeseries], GMotion)

#    The result is a list of rasters
     paste("The result is a", class(egms_ts), "of", length(egms_ts) , "rasters")

# 5.4. Stack the raster list into a multi-layer raster object
     egms_ts_r <- stack(egms_ts)
     
     dim(egms_ts_r)      # raster dimension

     boxplot(egms_ts_r)  #' Data distribution              
     
# 5.5. In some areas, interpolation can produce mathematical artifacts (outliers), especially at dataset boundaries,
     #'    skewing analysis and visualization. Apply this to identify and remove these extreme values in the raster.
     #'
     #'    Establish outlier boundaries for displacement
     disp_lower_limit = -200
     disp_upper_limit = 200
     
     #     Replace values outside these limits with NA.
     egms_ts_r2 = egms_ts_r          # preserve the original raster, work on the copy
     egms_ts_r2[egms_ts_r2 < disp_lower_limit | egms_ts_r2 > disp_upper_limit] <- NA
     
     #     Plot the raster again
     boxplot(egms_ts_r2)
     

# 5.6. Convert raster to SaptRaster object (terra), more efficient 
     egms_ts_rm = rast(egms_ts_r2)
#         Apply LAEA Europe CRS, then coordinates are in meters
          crs(egms_ts_rm) <- "EPSG:3035"


# 5.7. Masking areas with no data
     egms_ts_rm2 = egms_ts_rm # copy raster

#    offset variable (3.5.4 Masking) is the maximum distance for interpolation confidence
#    Change all values in egms_ts_rm2 to NA if the corresponding distance in distance_rm is greater than x
     egms_ts_rm2[distance_rm > offset] <- NA


#'##############
### 6. PLOTTING TIME SERIES OVER DYNAMIC MAP
#'##############

#  Project the raster to a geographic reference system, coordinates are in degrees.  
egms_ts_rd = project(egms_ts_rm2, "EPSG:4326")


# Define the limits of the color palette for displacement data
zlim_disp <- c(-150, 150)

# Create a colorNumeric object "mypal_disp"
mypal_disp <- colorNumeric(palette = egms_colors, domain = zlim_disp, na.color = NA) 


# Plot timeseries on a dynamic base map
p_ts <- plet(egms_ts_rd,                   # raster file
              y = c(1:nlyr(egms_ts_rd)),   # raster layers 
              alpha=0.7,                   # opacity
              tiles="Stamen.Terrain",      # basemap provider, more options: https://leaflet-extras.github.io/leaflet-providers/preview/
              col = mypal_disp,            # color palette
              legend = NULL,               # no legend
              ) %>%           
        addLegend(position = "bottomleft", # Add a legend to the map
              pal = mypal_disp,            # Legend colors
              values = zlim_disp,          # Legend values
              title = "Disp. mm",          # Legend title
              opacity = 0.7                # Legend opacity
             ) %>%             
        addDrawToolbar(       # add toolbar 
              polylineOptions     = TRUE,   # Enable polyline draw
              polygonOptions      = FALSE, 
              circleOptions       = FALSE, 
              circleMarkerOptions = FALSE, 
              markerOptions       = FALSE, 
              rectangleOptions    = FALSE,
             ) %>% 
        addReverseSearchOSM(  # add location to polyline nodes
              showSearchLocation = TRUE, 
                     showFeature = FALSE, 
                     fitBounds   = FALSE, 
                     displayText = FALSE
             ) 

p_ts

### export as a web page
htmlwidgets::saveWidget(p_ts, "Timeseries.html")


#'##############
### 7. TIMESERIES CROSS-SECTION
#'##############

#' 7.1. Create the cross-section to project the timeseries displacements

#'   Use option of draw a polyline in the dinamic map to get maker defining the cross-section coordinates 

#    cross-section coords and spatvector format
     line_coords <- rbind(c(11.58334, 44.52499), c(11.08664, 44.63011))

#    geographic projection
     ts_line <- vect(line_coords, "lines", crs = "EPSG:4326")

# 7.2. Extract displacement values from each raster cell that the cross-section intersects
     disp <- extract(egms_ts_rd, ts_line, ID = FALSE, xy = TRUE)   # use the lonlat raster

#    split disp file in coordinates and values
     disp_values = subset(disp, select = -c(x, y))
     disp_coord  = subset(disp, select =  c(x, y))                                      
#    round coordinates to 5 decimals
     disp_coord = round(disp_coord, 5) 

# 7.3. Plot the cross-section 
     line_length <- perim(ts_line) # Calculate the total length of the line

#    Create a sequence of distances along the line using the line length and the number of intersected cells
     distant <- seq(from = 0, to = line_length, length.out = nrow(disp))
   # round distant to 0 decimals
     distant = round(distant, 0)

#    Create a data frame with the distances and displacement values
     cross_sect <- data.frame(distance = distant, disp_values)
     dim(cross_sect)

#    Reshape the data frame to a long format using the pivot_longer function
#    This simplifies the process of plotting time series with differentiated colors 
#    by concatenating all displacement columns and assigning labels for each date 
     cross_sect_long = pivot_longer(cross_sect, cols = -1, names_to = "date", values_to = "displacement")

#    Remove 'X' character from date column
     cross_sect_long$date <- gsub("X", "", cross_sect_long$date)

#    Create the plot
     ggplot(cross_sect_long, aes(distance, displacement, color = date)) +
            geom_line() +
            theme(legend.position = "bottom") + 
            xlab("distance (m)") +
            ylab("displacement (mm)")

# 
# Attention!
#    The line is plotted using the first quadrant convention in GIS, where the origin (0, 0) is located at the most northwestern point.

#    Values and coordinates of the start and end of the cross-section
cat( 
  "Distance", head(cross_sect$distance, 1), ", coordinates", paste(head(disp_coord, 1), collapse=" "),
  "\nDistance", tail(cross_sect$distance, 1), ", coordinates", paste(tail(disp_coord, 1), collapse=" ")    
     )


# Plot the last line over the map
p_vel_pts_line = p_vel %>% lines(ts_line, col="white", lwd=2)
p_vel_pts_line

### export as a web page
htmlwidgets::saveWidget(p_vel_pts_line, "mean_velocity_points_line.html")


#'##############
### 8. EXPORT
#'##############

# Write the line to kml
writeVector(ts_line, filename = "ts_line.kml", filetype = "KML")

# write raster to geotiff
writeRaster(egms_vel_rd, filename = "egms_vel_rd.tif", filetype = "GTiff")
writeRaster(egms_ts_rd, filename = "egms_ts_rd.tif", filetype = "GTiff")

#'##############
### END
#'##############
#'







