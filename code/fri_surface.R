## Purpose: Calculate a fire order-weighted mean return interval surface across a landscape
## Author: Zack Steel

fri_surface <- function(landscape, # feature(s) that represent the landscape of interest
                        fires, #shapefile of fire perimeters including a year attribute
                        fire_years = "Year", # label of the fire year column
                        start_year, # Year prior to start of dataset (for landsat: 1983)
                        drop_censored_int = F, # last fire before dataset begins is unknown, do we assume the first interval starts on the start year (F) or drop this first censored interval (T) 
                        end_year = NULL, # Year after end of dataset 
                        decay_rate = 0.5, # Importance decay rate of the "invisible mosaic", between [0,1)
                        out_raster = NULL, #path if saving raster, if return
                        raster_template = "data/spatial/CBI_template.tif" # raster or path to raster
  ) {
  library(tidyverse)
  library(sf)
  library(terra)
  # library(here)
  
  ## some checks
  if(decay_rate > 1) stop("Cannot have a decay rate greater than 1")
  
  ## Pull in severity raster to use as template
  if(is.character(raster_template)) {
    r_template <- rast(raster_template)
  } else {
    r_template <- raster_template
  }
  
  ## assign fire year column name
  fires <- rename(fires, Fire_Year = !! (sym(fire_years))) |> 
    ## make sure Fire_Year is numeric
    mutate(Fire_Year = as.numeric(Fire_Year))
  
  ## if end year is not supplied use dataset
  if(is.null(end_year)) {
    end_year <- max(fires$Fire_Year) + 1
  }
  
  ## Only keep fires that occurred between the start and end years
  fires <- filter(fires, Fire_Year > start_year, Fire_Year < end_year) 
  
  ## Conform CRS of features to the template
  landscape <- st_transform(landscape, crs = st_crs(r_template))
  fires <- st_transform(fires, crs = st_crs(r_template))
  
  ## Which fires intersect with the landscape of interest?
  fires <- fires[landscape,]
  
  ## landscape raster of start year
  land_r <- vect(landscape) %>% 
    rast(resolution = res(r_template), vals = start_year) %>% 
    mask(vect(landscape))
  ## and end year
  end_r <- vect(landscape) %>% 
    rast(resolution = res(r_template), vals = end_year) %>% 
    mask(vect(landscape))

  ## FRI raster without fires
  noburn_r <- vect(landscape) %>%
    rast(resolution = res(r_template), vals = end_year - start_year) %>% 
    mask(vect(landscape))
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    message("No fires intersect landscape")
    
    ## if dropping first interval we don't have any intervals
    if(drop_censored_int) {
      noburn_r <- setValues(noburn_r, NA)
    }
    
    if(is.null(out_raster)) {
      return(noburn_r)
    } else {
      writeRaster(noburn_r, filename = out_raster, overwrite = T)
    }
    ## if there are fires
    } else {
    
    ## get all fire years
    years <- pull(fires, Fire_Year) %>%
      unique() %>% 
      sort()
    
    ## Created a raster for each year there was at least one fire
    fri_yrs <- lapply(years, function(x) {
      r <- filter(fires, Fire_Year == x) %>% 
        st_collection_extract("POLYGON") %>% 
        vect() %>% 
        terra::rasterize(y = land_r, field = "Fire_Year") %>% 
        suppressWarnings()
      
      ## Align with the landscape raster
      terra::resample(r, land_r)
    })
    
    ## Add start & end year raster
    fri_yrs <- append(list(land_r), fri_yrs) %>% 
      append(list(end_r))
    
    
    ## if dropping first interval we want the first raster to be NA
    if(drop_censored_int) {
      fri_yrs[[1]] <- setValues(fri_yrs[[1]], NA)
    }
    
    ## Get years since last fire for each year
    ## set up empty raster
    ysf_yrs <- lapply(fri_yrs, function(x) setValues(x, NA) )
    for(i in 2:length(fri_yrs)) {
      ## memory intensive step, clean up
      gc()
      ## make raster of most recent event
      ## pull all preceeding rasters
      rc <- rast(fri_yrs[(i-1):1]) %>% 
        ## get max year (most recent) for each pixel
        app(fun = max, na.rm = T)

      ## years since most recent fire
      r <- fri_yrs[[i]] - rc
      ysf_yrs[[i]] <- r
    }
    
    ## If just calculating years since last fire (decay = 1)
    ## We can stop here and take the last ysf_yrs raster
    if(decay_rate == 1) {
      out <- ysf_yrs[[length(ysf_yrs)]] %>% 
        mask(vect(landscape))
    } else {
      ## Assign burned areas as 1
      burned_l <- lapply(ysf_yrs, function(x) {
        x[!is.na(x)] <- 1
        x[is.na(x)] <- 0
        return(x)
      })
      burned_l[[1]] <- setValues(ysf_yrs[[1]], NaN)
      
      ## Add them up to get the cumulative number of fires
      ## Skip first empty year and add back in
      fire_cnt <- Reduce("+", burned_l[2:length(burned_l)], accumulate = T)
      fire_cnt <- append(burned_l[1], fire_cnt)
      
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
        ## exponentiating to zero gives a weight of 1 (in the case of the most recent interval; btwn last year and last fire), higher numbers get lower weights
        ## exponent term makes the fire order relative the max (e.g., if there have been 3 fires and we are on the third it gets full weight)
        (1 - decay_rate) ^ (maxorder - x)
      }) %>%
        rast()
      
      ## stack ysf layers
      ysf_stack <- rast(ysf_yrs)
      
      ## Get weighted average
      weighted_fri <- if(nlyr(ysf_stack) == 1) {ysf_stack[[1]]} else {
        weighted.mean(ysf_stack, w, na.rm = T)} 
      
      ## mask one more time
      out <- mask(weighted_fri, vect(landscape))
    }
    
    
    
    ## Return raster to user if path is not provided for writing raster
    if(is.null(out_raster)) {
      return(out)
    } else {
      writeRaster(out, filename = out_raster, overwrite = T)
    }
    
  }
  
}
