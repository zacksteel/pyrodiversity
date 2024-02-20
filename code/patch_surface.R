## Purpose: Calculate a fire order-weighted patch size (as defined by severity) surface across a landscape
## Author: Zack Steel

patch_surface <- function(landscape, # feature(s) that represent the landscape of interest
                          fires, # features representing fire severity
                          severity_dir, # directory where fire rasters are held
                          fire_years = "Year", # label of the fire year column
                          decay_rate = 0.5, # Between [0,1); a rate of .5 means each subsequent value recieves 1/4 of the previous
                          out_raster = NULL, #path if saving raster, if return
                          raster_template = "data/spatial/CBI_template.tif"
                     ) {
  library(terra)
  # library(stars)
  library(tidyverse)
  library(sf)
  library(units)
  library(here)
  
  ## some checks
  if(decay_rate > 1) stop("Cannot have a decay rate greater than 1")
  
  ## Pull in severity raster to use as template
  if(is.character(raster_template)) {
    r_template <- rast(here(raster_template))
  } else {
    r_template <- raster_template
  }
  
  ## Conform CRS of features to the template
  landscape <- st_transform(landscape, crs = st_crs(r_template))
  fires <- st_transform(fires, crs = st_crs(r_template))
  
  ## Which fires intersect with the landscape of interest?
  fires <- fires[landscape,]
  
  ## Create a landscape-wide empty raster 
  land_r <- vect(landscape) %>% 
    rast(resolution = res(r_template), vals = NA) %>% 
    mask(vect(landscape))
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    warning("No fires intersect landscape, returning a single patch")
    ## Calculate area
    landscape$log_area <- log(st_area(landscape))
    
    noburn_r <- terra::rasterize(vect(landscape), land_r, 
                          field = "log_area")
    
    if(is.null(out_raster)) {
      return(noburn_r)
    } else {
      writeRaster(noburn_r, filename = out_raster, overwrite = T)
    }
    
    ## if fires do more stuff
    } else {
      ## Make a list of severity rasters
      fire_files <- as.character(fires$Fire_ID) %>%
        sapply( function(x) list.files(severity_dir, pattern = x, full.names = T))
      
      miss <- sapply(fire_files, length) %>% table() %>% as.data.frame() %>% 
        filter(`.` == 0) %>% pull(Freq)
      
      if(length(miss) > 0) {stop(cat(miss, 'missing fires in severity library\n'))}
    
      ## loop through each
      rasters <- lapply(1:length(fire_files), function(i) {
        
        ## pull out fire and info
        fire <- fires[i,]
        id <- as.character(fire$Fire_ID)
        year <- pull(fire, {{fire_years}})
        # year <- pull(fire, Fire_Year)
        path <- fire_files[[i]]
        
        ## Get cbi raster
        r <- rast(path)
        names(r) <- as.character(year)
        return(r)
      })
      
      ## And remove any nulls
      rasters[sapply(rasters, is.null)] <- NULL
    
      ## get vector of burn years
      years <- sapply(rasters, names) %>% 
        as.integer() %>% 
        unique() %>% 
        sort()
    
    ## Created a raster for each year there was at least one fire
    patch_yrs <- list()
    
    #### This is a slow step
    for(x in years) {
      ## Which match the year in question
      keep <- sapply(rasters, names) == x
      fs <- rasters[keep]
      
      ## Merge (mosaic sometimes fails) if more than one, otherwise unlist and return
      if(length(fs) == 1) {m <- fs[[1]]} else
      {
        ## set up SpatRasterCollection avoids awkward do.call
        rsrc <- sprc(fs)
        ## use try catch around the following function to capture warnings
        tryCatch({
          # m <- mosaic(rsrc, fun = 'max')
          m <- merge(rsrc)
        },
        warning = function(w) {
          # Print custom message and stop the script
          print("A warning occurred during the merge operation. Stopping the script.")
          # Optionally, print the original warning message
          print(conditionMessage(w))
          # Stop execution
          stop("Script stopped due to a warning in merge operation.")
        })
      }
      
      ## Convert to features with standard cut-offs from Miller and Thode 2007
      breaks <- c(0, .1, 1.25, 2.25, 3)
      rc <- classify(m, rcl = breaks, include.lowest = T)
      
      ## Convert to vector
      fc <- as.polygons(rc, dissolve = T) %>% 
        disagg()
      ## Calculate log-area
      fc$log_area <- log(expanse(fc, unit = 'ha'))
      ## Back to raster
      r_area <- terra::rasterize(fc, rc, field = 'log_area')
      
      ## Align with the landscape raster
      out <- terra::resample(r_area, land_r)
      
      patch_yrs <- c(patch_yrs, out)
    }

    
    ## Manipulate rasters to set up order weighting
    ## Assign burned areas as 1
    burned_l <- lapply(patch_yrs, function(x) {
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
    
    ## If decay rate is one (i.e., only interested in last fire/maxorder) then assign weighting directly. zero raised to anything returns a 1 here for some reason
    if(decay_rate == 1) {
      
      w <- lapply(fire_order, function(x) {
        ## which pixels were the last fire?
        tmp.r <- maxorder == x
        ## assign a 1
        tmp.r[] <- ifelse(tmp.r[] == TRUE, 1, 0)
        return(tmp.r)
      }) %>% 
        rast()
      
    } else 
      ## for decay rates less than 1
    {
      w <- lapply(fire_order, function(x) {
        ## exponentiating to zero gives a weight of 1, higher numbers get lower weights
        ## exponent term makes the fire order relative the max (e.g., if there have been 3 fires and we are on the third it gets full weight)
        (1 - decay_rate) ^ (maxorder - x)
      }) %>%
        ## 'stack' the list of rasters
        rast()
    }
    
    ## stack severity layers
    patch_stack <- rast(patch_yrs)
    
    ## Get weighted average
    ## Take out garbage, the next step is a big memory suck
    gc()
    weighted_patch <- if(nlyr(patch_stack) == 1) {patch_yrs[[1]]} else {
      weighted.mean(patch_stack, w, na.rm = T)}
    
    ## Keep just the area within the buffered landscape
    land_patch <- mask(weighted_patch, vect(landscape))
    
    ## Calculate patches of unburned area; pretty slow
    ub_poly <- suppressMessages(st_difference(landscape, fires[,"geometry"])) %>% 
      vect() %>% 
      disagg()
    ## Calculate log-area
    ub_poly$log_area <- log(expanse(ub_poly, unit = 'ha'))
    ## Back to raster
    ub_r <- terra::rasterize(ub_poly, land_patch, field = 'log_area')
    
    ## Combine weighted burned and unburned
    land_patch2 <- terra::merge(land_patch, ub_r)
    
    ## Return raster to user if path is not provided for writing raster
    if(is.null(out_raster)) {
      return(land_patch2)
    } else {
      writeRaster(land_patch2, filename = out_raster, overwrite = T)
    }
  }
  
}