# GMotion
This R script is dedicated to postprocess the ground motion data from the Copernicus EGMS platform. 
The code reads the EGMS CSV files, subsets the data according to specified region of interest, 
and generates interactive maps and plots to conduct timeseries analysis. 
The main products are timeseries multilayer maps and timeseries cross-sections.

This is the RStudio version, which uses the subsidence in Bologna as a case study. 
It serves as a practical exercise in the context of the advanced course RawMatCop-Alliance, focusing on 
Earth Observation applied to Raw Materials. The EGMS data can be found in the working directory .../Data

Although this code can be applied to any zip file downloaded from the Copernicus EGMS platform, 
this code is parameterized for the subsidence experienced in the Bologna area, Italy, file EGMS_L3_E44N23_100km_U.zip
Files do not need to be unziped. 

Main parameters that need to be modifed for appling to a new area are: 
1. The ROI.................. W-E & S-N limits
2. grid_step.................for the raster spacing. Use large values when running the first time to reduce time consuming.
3. vel_lower_limit, vel_upper_limit, disp_lower_limit, disp_upper_limit............lower & upper limits to remove outliers
4. offset....................maximum distance for interpolation confidence
5. zlim_v and zlim_disp......limits for the color palettes
6. first_col, last_col, leap_col.........for the timeseries monitoring period
7. line_coords...............cross-section coordinates 
