## Purpose: Create a raster of the number of fires observed within a given time period
## Author: Zack Steel

nfires_surface <- function(
    landscape, # feature(s) that represent the landscape of interest
    fires, #shapefile of fire perimeters including a year attribute
    year_lab = "Year", # label of the fire year column
    years = NULL, # a two value vector giving the range of fire years to consider, if null all years in periemter dataset will be included
    out_res = 30, # resolution of raster output, defaults to 30m
    out_raster = NULL, #path if saving raster, if return
    return_rast = F #return the raster to user?
)
{
  library(tidyverse)
  library(sf)
  library(terra)
  
  ## assign fire year column name
  fires <- rename(fires, Fire_Year = !! (sym(year_lab)))
  
  ## Limit by landscape of interest
  landscape <- st_transform(landscape, crs = st_crs(fires))
  pers <- fires[landscape,] 
  
  ## and years of interest if specified
  if(!is.null(years)) {
    pers <- filter(pers, Fire_Year >= min(years), Fire_Year <= max(years))
  }
  
  # Create raster template 
  r <- rast(ext = ext(landscape), res = out_res, crs = crs(pers))
  
  # Initialize raster values to 0  
  r[] <- 0
  
  # Loop through polygon years
  yrs <- unique(pers$Fire_Year) %>% 
    sort()
  
  for(yr in yrs){
    
    ## pull out and dissolve year
    pers.yr <- filter(pers, Fire_Year == yr) %>% 
      summarise() %>% 
      mutate(value = 1)
    
    # Convert polygon to raster
    pi <- rasterize(pers.yr, r, field = "value") 
    
    # Add 1 to cells where polygon raster == 1
    r[pi == 1] <- r[pi == 1] + 1
    
  }
  
  ## Save resulting raster and return to user if specified, return if no raster path is provided
  if(!is.null(out_raster)) {
    writeRaster(r, filename = out_raster)
  }
  
  if(return_rast | !is.null(out_raster)){
    return(r)
  } 
  
}