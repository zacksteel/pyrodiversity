## Purpose: Calculate a fire order-weighted mean return interval surface across a landscape
## Author: Zack Steel

fri_surface <- function(landscape, # feature(s) that represent the landscape of interest
                        fires, #shapefile containing fire IDs associated with severity rasters
                        ID, #label of the unique landscape identifier
                        fire_years = "Year", # label of the fire year column
                        start_year, # Year prior to start of dataset (for landsat: 1983)
                        end_year, # Year after end of dataset (for MTBS: 2018)
                        decay_rate = 0.5, # Importance decay rate of the "invisible mosaic", between [0,1)
                        out_dir #path to hold output rasters
  ) {
  library(tidyverse)
  library(sf)
  library(raster)
  library(stars)
  library(fasterize)
  
  ## filename of output
  fn_r <- paste0(out_dir, "/fri_", as.data.frame(landscape)[1,ID], ".tif")
  
  ## Which fires intersect with the landscape of interest?
  keep <- suppressMessages(st_intersects(fires, landscape)) %>%
    apply(1, any) 
  fires <- fires[keep,]
  
  ## Pull in severity raster to use as template
  r_template <- raster("data/spatial/CBI_template.tif")
  land_r <- raster(landscape, resolution = res(r_template), vals = 1) %>% 
    mask(landscape)
  noburn_r <- raster(landscape, resolution = res(r_template), vals = end_year - start_year) %>% 
    mask(landscape)
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    message("No fires intersect landscape")
    
    writeRaster(noburn_r, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r) } 
  
  else {

    ## limit to landscape
    fires_land <- suppressMessages(st_intersection(fires, landscape))
    
    ## add starting year to base landscape
    landscape <- mutate(landscape,
                        year = start_year)
    land_year <- dplyr::select(landscape, year)
    ## simplify to just year attributes and "dissolve"
    f_year <- dplyr::select(fires_land, year = one_of(fire_years)) %>% 
      ## add landscape feature as a starting year
      rbind(f_year,land_year) %>%
      group_by(year) %>%
      summarise() 
    
    ## Attempt to run an intersection
    ## If topology errors crop up, rebuild from rasters, much slower but less error prone
    
    ## Make index table of year
    i_tab <- data.frame(index = 1:nrow(f_year), year = f_year$year)
    
    f_int <- try(st_intersection(f_year), silent = T)
    
    if("try-error" %in% class(f_int)) {
      message("Topology issues, rebuilding features (slower)")
    }
    
    if("try-error" %in% class(f_int)) {
      ## Convert to raster and then back again
      ## orders of magnitude faster with stars package than raster
      ## https://r-spatial.github.io/stars/articles/stars5.html
      f_year2  <- lapply(1:nrow(f_year), function(x) 
        st_rasterize(f_year[x,"year"])) %>% 
        lapply(function(x) st_as_sf(x, merge = T))
      f_year2 <- do.call(rbind, f_year2) %>% 
        ## make multipolygons
        group_by(year) %>% 
        summarise() %>% 
        st_buffer(dist = 0)
      
      f_int <- st_intersection(f_year2)
      
      ## Make index table of year
      i_tab <- data.frame(index = 1:nrow(f_year2), year = f_year2$year)
    }
    
    
    ## map years to origin index and subsequent calculations
    d <- mutate(f_int,
                years = map(origins, function(x) i_tab[x,"year"]),
                ## Add ending year
                years = map(years, function(x) c(x, end_year)),
                ## get intervals between years
                intervals = map(years, diff),
                ## Set order of events (intervals in this case)
                order = map(intervals, function(x) length(x):1),
                ## calculate mean return interval using specified decay rate = [0,1)
                ## exponentiating to zero gives a weight of 1, higher numbers get lower weights
                fri = map2(intervals, order, function(x, y) weighted.mean(x, (1 - decay_rate)^(y-1))),
                ## Also get years since last fire
                yslf = map(intervals, last),#function(x) x[[1]]),
                ## unlist 
                fri = unlist(fri),
                yslf = unlist(yslf))
    
    ## Pull out polygons (i.e. drop linesstrings from collections)
    # d2 <- st_collection_extract(d, "POLYGON")
    d2 <- d %>% 
      filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON"))
    
    fri_r <- fasterize(d2, noburn_r, field = "fri") 
    
    ## save
    writeRaster(fri_r, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r)
  }
  
  return(out)
  
}