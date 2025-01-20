## Purpose: Calculate a fire order-weighted burn severity surface across a landscape
## Author: Zack Steel

sev_surface <- function(landscape, # feature(s) that represent the landscape of interest
                        fires, # features representing fire severity
                        severity_dir, # directory where fire rasters are held
                        fire_years = "Year", # label of the fire year column
                        decay_rate = 0.5, # Between [0,1); a rate of .5 means each subsequent value recieves 1/2 of the previous
                        out_raster = NULL, #path if saving raster, if return
                        raster_template = "data/spatial/CBI_template.tif"
                     ) {
  library(tidyverse)
  library(sf)
  library(terra)
  library(lwgeom)
  
  ## some checks
  if(decay_rate > 1) stop("Cannot have a decay rate greater than 1")
  ## make sure fires has a Fire_ID column, stop if not
  if(!"Fire_ID" %in% names(fires)) stop("Fires must have a Fire_ID column to match severity rasters")
  message(paste0("Running with decay rate of ", decay_rate))
  
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
  fires1 <- fires[landscape,]
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires1) == 0) {
    warning("No fires intersect landscape, returning a zero severity landscape")
  
    if(is.null(out_raster)) {
      return(noburn_r)
    } else {
      writeRaster(noburn_r, filename = out_raster, overwrite = T)
    }
    } else {
  
  ## otherwise calculate a weighted severity landscape

    
    ## Make a list of severity rasters; they will be empty if not in the directory
    fire_files <- as.character(fires1$Fire_ID) %>%
      sapply( function(x) list.files(severity_dir, pattern = x, full.names = T))
    
    miss <- sapply(fire_files, length) %>% table() %>% as.data.frame() %>% 
      filter(`.` == 0) %>% pull(Freq)
    
    if(length(miss) > 0) {stop(cat(miss, 'missing fires in severity library\n'))}
    
    ## loop through each
    rasters <- lapply(1:length(fire_files), function(i) {
      # print(i)
      ## pull out fire and info
      fire <- fires1[i,]
      id <- as.character(fire$Fire_ID)
      year <- pull(fire, {{fire_years}})
      
      ## Get cbi raster
      r <- rast(fire_files[[i]]) 
      
      names(r) <- as.character(year)
      return(r)
    })
    
    ## Create a landscape-wide empty raster 
    land_r <- vect(landscape) %>% 
      rast(resolution = res(rasters[[1]]), vals = NA)
    
    ## get vector of burn years
    years <- pull(fires1, {{fire_years}}) %>%
      unique() %>% 
      ## sort them in ascending order
      sort()
    
    ## Created a raster for each year there was at least one fire
    sev_yrs <- lapply(years, function(x) {
      # print(x)
      ## Which match the year in question
      keep <- sapply(rasters, names) == x
      fs <- rasters[keep]
      
      ## Mosaic if more than one, otherwise unlist and return
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
      
      ## Align with the landscape raster
      return(terra::resample(m, land_r))
      
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
    sev_stack <- rast(sev_yrs)
    
    ## Clean up
    gc()
    
    ## Get weighted average
    weighted_sev <- if(nlyr(sev_stack) == 1) {sev_yrs[[1]]} else {
      weighted.mean(sev_stack, w, na.rm = T) 
      } ## unburned areas are NA so must not be included in the weighted mean
    
    ## Keep just the area within the buffered landscape
    land_sev <- crop(weighted_sev, vect(landscape))
    
    ## Any pixels that haven't burned assign zero
    land_sev[is.na(land_sev)] <- 0
    
    ## make everything outside the landscape NA again
    land_sev <- mask(land_sev, vect(landscape))
    
    ## Return raster to user if path is not provided for writing raster
    if(is.null(out_raster)) {
      return(land_sev)
    } else {
      writeRaster(land_sev, filename = out_raster, overwrite = T)
    }
  }
  
}