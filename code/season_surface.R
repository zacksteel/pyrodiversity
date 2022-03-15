## Purpose: Calculate a fire order-weighted seasonality (ignition date) surface across a landscape
## Author: Zack Steel

season_surface <- function(landscape, # feature(s) that represent the landscape of interest
                           fires, # features representing fire severity
                           fire_years = "Year", # label of the fire year column
                           fire_day = "jday", # label of the fire ignition Julian date column
                           decay_rate = 0.5, # Between [0,1); a rate of .5 means each subsequent value recieves 1/2 of the previous
                           out_raster = NULL, #path if saving raster, if return
                           raster_template = "data/spatial/CBI_template.tif"
                        ) 
  {
  library(tidyverse)
  library(sf)
  library(terra)
  library(here)
  
  ## Pull in severity raster to use as template
  r_template <- rast(here(raster_template))
  
  ## Conform CRS of features to the template
  landscape <- st_transform(landscape, crs = st_crs(r_template))
  fires <- st_transform(fires, crs = st_crs(r_template))
  
  ## Which fires intersect with the landscape of interest?
  keep <- suppressMessages(st_intersects(fires, landscape)) %>%
    apply(1, any) 
  fires <- fires[keep,]
  
  ## assign fire year column name
  fires <- rename(fires, Fire_Year = !! (sym(fire_years)))
  # fires <- rename(fires, Fire_Year = Year)
  
  noburn_r <- vect(landscape) %>% 
    rast(resolution = res(r_template), vals = NA)
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    warning("No fires intersect landscape, returning NA landscape")
    
    if(is.null(out_raster)) {
      return(noburn_r)
    } else {
      writeRaster(noburn_r, filename = out_raster, overwrite = T)
    }
  }
    ## if there are fires
  else {
    ## Make sure day is between 1 and 365 (ie. that is looks like a Julian day)
    tmp <- dplyr::select(fires, fire_day) %>%
      pull(fire_day)
    if(min(tmp) < 1 | max(tmp) > 365) {
      stop("fire_day attribute must be Julian day format, but at least one value is less than 1 or greater than 365. \nPlease check format")
    }
    
    ## add NA day to base landscape
    landscape <- mutate(landscape,
                        jday = NA)
    
    ## simplify to just year and day attributes, "dissolve" and convert to cosine of radians to get circular date
    f_day <- dplyr::select(fires, year = Fire_Year, jday = all_of(fire_day)) %>%
      ## add landscape feature as a starting year
      ## make unique year and day
      mutate(yrdy = as.integer(paste0(year, jday))) %>% 
      group_by(year, jday, yrdy) %>%
      summarise(geometry = suppressMessages(st_union(geometry)),
                .groups = "drop_last") %>%
      ## Convert jday to degrees then radians
      mutate(deg_day = (jday/366) * 360,
             rad_day = deg_day*pi / 180,
             crad_day = cos(rad_day)) %>% 
      ungroup() 
    
    ## get all years
    years <- pull(f_day, year) %>%
      unique() %>% 
      sort()
    
    ## Created a raster for each year there was at least one fire
    sea_yrs <- lapply(years, function(x) {
      r <- filter(f_day, year == x) %>% 
        st_collection_extract("POLYGON") %>% 
        vect() %>% 
        terra::rasterize(y = noburn_r, fun = "last", field = "jday") %>% 
        suppressWarnings()
      
      ## Align with the landscape raster
      return(raster::resample(r, noburn_r))
      
    })
    
    ## Assign burned areas as 1
    burned_l <- lapply(sea_yrs, function(x) {
      x[!is.na(x)] <- 1
      x[is.na(x)] <- 0
      return(x)
    })
    ## Add them up to get the cummulative number of fires
    fire_cnt <- Reduce("+", burned_l, accumulate = T)
    
    ## Convert back to NA if a pixel didn't burn in a given year
    fire_order <- lapply(1:length(burned_l), function(x) {
      x <- burned_l[[x]] * fire_cnt[[x]]
      x[x == 0] <- NA
      return(x)
    })
    
    ## Created a stack of weight rasters 
    ## create raster of max fires so we are flipping the weight and 
    ## making recent burns more important
    stck <- rast(fire_order)
    maxorder <- if(nlyr(stck) == 1) {1} else {
      max(stck, na.rm = T)
    }
    
    w <- lapply(fire_order, function(x) {
      ## exponentiating to zero gives a weight of 1, higher numbers get lower weights
      (1 - decay_rate) ^ (maxorder - x)
    }) %>%
      rast()
    
    ## stack severity layers
    sea_stack <- rast(sea_yrs)
    
    ## Get weighted average
    weighted_sea <- if(nlyr(sea_stack) == 1) {sea_yrs[[1]]} else {
      weighted.mean(sea_stack, w, na.rm = T)}
    
    ## Keep just the area within the buffered landscape
    land_sea <- crop(weighted_sea, vect(landscape)) %>% 
      mask(vect(landscape))
    
    ## Return raster to user if path is not provided for writing raster
    if(is.null(out_raster)) {
      return(land_sea)
    } else {
      writeRaster(land_sea, filename = out_raster, overwrite = T)
    }
    
    
  }
  
}
