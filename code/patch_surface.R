## Purpose: Calculate a fire order-weighted patch size (as defined by severity) surface across a landscape
## Author: Zack Steel

patch_surface <- function(landscape, # feature(s) that represent the landscape of interest
                          fires, # features representing fire severity
                          ID, #label of the unique landscape identifier
                          severity_dir, # directory where fire rasters are held
                          fire_years = "Year", # label of the fire year column
                          decay_rate = 0.5, # Between [0,1); a rate of .5 means each subsequent value recieves 1/4 of the previous
                          buffer_landscape = 0, # Extension of landscape (in meters) if needed for point calculations near edges
                          ## Buffer could also be used to reduce edge effects on patch area calculations
                          out_dir #path to hold output rasters
                     ) {
  library(raster)
  library(stars)
  library(tidyverse)
  library(sf)
  library(fasterize)
  
  ## Simplify landscape if needed; Created a buffered landscape
  land_buf <- st_union(landscape) %>%
    ## Switch to utm for buffering
    st_transform(3310) %>%
    st_buffer(buffer_landscape) %>%
    st_transform(crs = st_crs(landscape)) %>%
    ## need to assign a geometry column
    st_sf() 
  
  ## filename of output
  fn_r <- paste0(out_dir, "/pat_", as.data.frame(landscape)[1,ID], ".tif")
  
  ## Which fires intersect with the landscape of interest?
  ## Fires should also be utm
  keep <- suppressMessages(st_intersects(fires, land_buf)) %>%
    apply(1, any) 
  fires <- fires[keep,] 
  
  ## Pull in severity raster to use as template
  r_template <- raster("data/spatial/CBI_template.tif")
  
  ## Create a landscape-wide empty raster 
  #### for some reason crs isn't carrying over automatically anymore (R upgrade?)
  land_r <- as_Spatial(land_buf) %>% 
    raster(resolution = res(r_template), vals = NA) %>% 
    mask(land_buf)
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    warning("No fires intersect landscape, returning a single patch")
    ## Calculate area
    land_buf$log_area <- log(st_area(land_buf))
    
    noburn_r <- fasterize(land_buf, land_r, 
                          field = "log_area")
    
    writeRaster(noburn_r, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r)
    } else {
      ## Make a list of severity rasters
      fire_files <- as.character(fires$Fire_ID) %>%
        sapply( function(x) list.files(severity_dir, pattern = x, full.names = T))
    
      ## loop through each
      rasters <- lapply(1:length(fire_files), function(i) {
        
        ## pull out fire and info
        fire <- fires[i,]
        id <- as.character(fire$Fire_ID)
        year <- pull(fire, {{fire_years}})
        path <- fire_files[[i]]
        
        ## Get cbi raster
        r <- raster(path, layer = 1)
        names(r) <- as.character(year)
        return(r)
      })
      
      ## And remove any nulls
      rasters[sapply(rasters, is.null)] <- NULL
    
      ## get vector of burn years
      years <- sapply(rasters, names) %>% 
        substring(2) %>% 
        as.integer()
    
    ## Created a raster for each year there was at least one fire
    patch_yrs <- list()
    
    #### This is a very slow step
    for(x in years) {
      ## Which match the year in question
      keep <- sapply(rasters, names) == paste0("X", x)
      fs <- rasters[keep]
      
      ## Mosaic if more than one, otherwise unlist and return
      if(length(fs) == 1) {m <- fs[[1]]} else
      {
        ## Set up awkward do.call for list of rasters
        fs2 <- fs
        fs2$fun <- max # if overlapping fires use max severity
        fs2$na.rm <- TRUE
        
        m <- do.call(mosaic, fs2)
      }
      
      ## Convert to features with standard cut-offs from Miller and Thode 2007
      breaks <- c(0, .1, 1.25, 2.25, 3)
      rc <- cut(m, breaks = breaks)
        ## Convert to vector
        ## stars version is about 3 times faster
        ## https://r-spatial.github.io/stars/articles/stars5.html
        
      fc <- st_as_stars(rc) %>% 
        st_as_sf() %>% 
        group_by(layer) %>%
        summarize()
      
      ## convert from multi to single polygons
      fc_polys <- st_cast(fc, "MULTIPOLYGON", warn = F) %>% # sometimes necessary but not always
        st_cast("POLYGON", warn = F) %>%
        ## Calculate log-area and convert back to raster
        mutate(log_area = log(st_area(.)))
      
      r_area <- fasterize(fc_polys, rc, field = "log_area")
      
      
      ## Align with the landscape raster
      out <- raster::resample(r_area, land_r)
      
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
    patch_stack <- stack(patch_yrs)
    
    ## Get weighted average
    weighted_patch <- if(nlayers(patch_stack) == 1) {patch_yrs[[1]]} else {
      weighted.mean(patch_stack, w, na.rm = T)}
    
    ## Keep just the area within the buffered landscape
    land_patch <- mask(weighted_patch, land_buf)
    
    ## Calculate patches of unburned area
    ub_poly <- suppressMessages(st_difference(land_buf, fires[,"geometry"])) %>%
      st_cast("MULTIPOLYGON", warn = F) %>% # sometimes necessary but not always
      st_cast("POLYGON", warn = F) %>%
      mutate(log_area = log(st_area(.)))
    ub_r <- fasterize(ub_poly, land_patch, field = "log_area")
    
    ## Combine weighted burned and unburned
    land_patch2 <- mosaic(land_patch, ub_r, fun=min)
    
    ## save
    writeRaster(land_patch2, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r)
  }
  
  return(out)
  
}