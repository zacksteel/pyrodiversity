## Purpose: Calculate a fire order-weighted mean return interval surface across a landscape
## Author: Zack Steel

fri_surface <- function(landscape, # feature(s) that represent the landscape of interest
                        fires, #shapefile containing fire IDs associated with severity rasters
                        ID, #label of the unique landscape identifier
                        fire_years = "Year", # label of the fire year column
                        start_year, # Year prior to start of dataset (for landsat: 1983)
                        end_year, # Year after end of dataset (for MTBS: 2018)
                        decay_rate = 0.5, # Importance decay rate of the "invisible mosaic", between [0,1)
                        out_dir, #path to hold output rasters
                        raster_template = "data/spatial/CBI_template.tif"
  ) {
  library(tidyverse)
  library(sf)
  library(raster)
  library(stars)
  library(fasterize)
  
  ## filename of output
  fn_r <- paste0(out_dir, "/fri_", as.data.frame(landscape)[1,ID], ".tif")
  
  ## Pull in severity raster to use as template
  r_template <- raster(raster_template)
  
  ## Conform CRS of features to the template
  landscape <- st_transform(landscape, crs = st_crs(r_template))
  fires <- st_transform(fires, crs = st_crs(r_template))
  
  ## Which fires intersect with the landscape of interest?
  keep <- suppressMessages(st_intersects(fires, landscape)) %>%
    apply(1, any) 
  fires <- fires[keep,]
  
  ## assign fire year column name
  fires <- rename(fires, Fire_Year = !! (sym(fire_years)))

  ## landscape raster of start year
  land_r <- as_Spatial(landscape) %>% 
    raster(resolution = res(r_template), vals = start_year) %>% 
    mask(landscape)
  ## and end year
  end_r <- as_Spatial(landscape) %>% 
    raster(resolution = res(r_template), vals = end_year) %>% 
    mask(landscape)

  ## FRI raster without fires
  noburn_r <- as_Spatial(landscape) %>%
    raster(resolution = res(r_template), vals = end_year - start_year) %>% 
    mask(landscape)
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    message("No fires intersect landscape")
    
    writeRaster(noburn_r, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r) } else {

    ## limit to landscape
    fires_land <- suppressMessages(st_intersection(fires, landscape))
    
    ## get all years
    years <- pull(fires_land, Fire_Year) %>%
      unique() %>% 
      sort()
    
    ## Created a raster for each year there was at least one fire
    #### Slow step here, stars or terra versions of rasterize and resample might help
    fri_yrs <- lapply(years, function(x) {
      r <- filter(fires_land, Fire_Year == x) %>% 
        st_collection_extract("POLYGON") %>% 
        fasterize(raster = land_r, field = "Fire_Year") %>% 
        suppressWarnings()
      
      ## Align with the landscape raster
      return(raster::resample(r, land_r))
      
    })
    
    ## Add start & end year raster
    fri_yrs <- c(land_r, fri_yrs, end_r)
    
    ## Get years since last fire for each year
    ## set up empty raster
    ysf_yrs <- lapply(fri_yrs, function(x) x * NA)
    for(i in 2:length(fri_yrs)) {
      ## make raster of most recent event
      #### stars or terra version of calc may improve speed here as well
      rc <- calc(stack(fri_yrs[i-1:i]), fun = max)
      ## years since most resent fire
      r <- fri_yrs[[i]] - rc
      ysf_yrs[[i]] <- r
    }
    
    ## Assign burned areas as 1
    burned_l <- lapply(ysf_yrs, function(x) {
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
    
    ## stack ysf layers
    ysf_stack <- stack(ysf_yrs)
    
    ## Get weighted average
    weighted_fri <- if(nlayers(ysf_stack) == 1) {ysf_stack[[1]]} else {
      weighted.mean(ysf_stack, w, na.rm = T)}
    

    ## save
    writeRaster(weighted_fri, filename = fn_r, overwrite = T)

    out <- basename(fn_r)
  }
  
  return(out)
  
}
