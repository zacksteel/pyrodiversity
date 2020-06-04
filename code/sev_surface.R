## Purpose: Calculate a fire order-weighted burn severity surface across a landscape
## Author: Zack Steel

sev_surface <- function(landscape, # feature(s) that represent the landscape of interest
                        fires, # features representing fire severity
                        ID, #label of the unique landscape identifier
                        severity_dir, # directory where fire rasters are held
                        fire_years = "Year", # label of the fire year column
                        decay_rate = 0.5, # Between [0,1); a rate of .5 means each subsequent value recieves 1/2 of the previous
                        out_dir #path to hold output rasters
                     ) {
  library(tidyverse)
  library(sf)
  library(raster)
  library(lwgeom)
  
  ## filename of output
  fn_r <- paste0(out_dir, "/sev_", as.data.frame(landscape)[1,ID], ".tif")
  
  ## Which fires intersect with the landscape of interest?
  ## Fires should also be utm
  keep <- suppressMessages(st_intersects(fires, landscape)) %>%
    apply(1, any) 
  fires <- fires[keep,] 
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    warning("No fires intersect landscape, returning a zero severity landscape")
    
    ## Pull in severity raster to use as template
    r_template <- raster("data/spatial/CBI_template.tif")
    noburn_r <- raster(landscape, resolution = res(r_template), vals = 0) %>% 
      mask(landscape)
    crs(noburn_r) <- crs(r_template)

    writeRaster(noburn_r, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r) }
  
  ## otherwise calculate a weighted severity landscape
  else {
    
    ## Make a list of severity rasters
    fire_files <- as.character(fires$Fire_ID) %>%
      sapply( function(x) list.files(severity_dir, pattern = x, full.names = T))
    
    ## loop through each
    rasters <- lapply(1:length(fire_files), function(i) {
      
      ## pull out fire and info
      fire <- fires[i,]
      id <- as.character(fire$Fire_ID)
      year <- pull(fire, {{fire_years}})
      
      ## Get cbi raster
      r <- raster(fire_files[[i]], layer = 1) 
      
      names(r) <- as.character(year)
      return(r)
    })
    
    ## Create a landscape-wide empty raster 
    land_r <- raster(landscape, resolution = res(rasters[[1]]), vals = NA)
    crs(land_r) <- crs(rasters[[1]])
    
    ## get vector of burn years
    years <- pull(fires, {{fire_years}}) %>%
      unique()
    
    ## Created a raster for each year there was at least one fire
    sev_yrs <- lapply(years, function(x) {
      # print(x)
      ## Which match the year in question
      keep <- sapply(rasters, names) == paste0("X", x)
      fs <- rasters[keep]
      
      ## Mosaic if more than one, otherwise unlist and return
      if(length(fs) == 1) {m <- fs[[1]]} else
      {
        ## Set up awkward do.call for list of rasters
        # names(fs)[1:2] <- c('x','y')
        fs$fun <- max # if overlapping fires use max severity
        fs$na.rm <- TRUE
        
        m <- do.call(mosaic, fs)
      }
      
      ## Align with the landscape raster
      return(raster::resample(m, land_r))
      
    })
    
    ## make sure there are no negative values
    sev_yrs <- lapply(sev_yrs, function(x) {
      x[x<0] <- 0
      return(x)
    })
    
    ## Assign burned areas as 1
    burned_l <- lapply(sev_yrs, function(x) {
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
    stck <- stack(fire_order)
    maxorder <- if(nlayers(stck) == 1) {1} else {
      max(stck, na.rm = T)
    }
    
    w <- lapply(fire_order, function(x) {
      ## exponentiating to zero gives a weight of 1, higher numbers get lower weights
      (1 - decay_rate) ^ (maxorder - x)
    }) %>%
      stack()
    
    ## stack severity layers
    sev_stack <- stack(sev_yrs)
    
    ## Get weighted average
    weighted_sev <- if(nlayers(sev_stack) == 1) {sev_yrs[[1]]} else {
      weighted.mean(sev_stack, w, na.rm = T)}
    
    ## Keep just the area within the buffered landscape
    land_sev <- crop(weighted_sev, landscape)
    
    ## Any pixels that haven't burned assign zero
    land_sev[is.na(land_sev)] <- 0
    
    ## make everything outside the landscape NA again
    land_sev <- mask(land_sev, landscape)
    
    ## save
    writeRaster(land_sev, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r)
  }
  
  return(out)
  
}