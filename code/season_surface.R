## Purpose: Calculate a fire order-weighted seasonality (ignition date) surface across a landscape
## Author: Zack Steel

season_surface <- function(landscape, # feature(s) that represent the landscape of interest
                           fires, # features representing fire severity
                           ID, #label of the unique landscape identifier
                           fire_years = "Year", # label of the fire year column
                           fire_day = "jday", # label of the fire ignition Julian date column
                           decay_rate = 0.5, # Between [0,1); a rate of .5 means each subsequent value recieves 1/2 of the previous
                           out_dir #path to hold output rasters
                        ) 
  {
  library(tidyverse)
  library(sf)
  library(raster)
  library(stars)
  library(fasterize)
  
  ## filename of output
  fn_r <- paste0(out_dir, "/sea_", as.data.frame(landscape)[1,ID], ".tif")
  
  ## Which fires intersect with the landscape of interest?
  keep <- suppressMessages(st_intersects(fires, landscape)) %>%
    apply(1, any) 
  fires <- fires[keep,]
  
  ## Pull in severity raster to use as template
  r_template <- raster("data/spatial/CBI_template.tif")
  noburn_r <- as_Spatial(landscape) %>% 
    raster(resolution = res(r_template), vals = NA)
  
  ## If no fires within the landscape skip a lot and treat the landscape as a single patch
  if(nrow(fires) == 0) {
    warning("No fires intersect landscape, returning NA landscape")

    writeRaster(noburn_r, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r) } 
  
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
    f_day <- dplyr::select(fires, year = one_of(fire_years), jday = one_of(fire_day)) %>%
      ## add landscape feature as a starting year
      ## make unique year and day
      mutate(yrdy = as.integer(paste0(year, jday))) %>% 
      group_by(year, jday, yrdy) %>%
      summarise(geometry = st_union(geometry)) %>%
      ## Convert jday to degrees then radians
      mutate(deg_day = (jday/366) * 360,
             rad_day = deg_day*pi / 180,
             crad_day = cos(rad_day)) %>% 
      ungroup() 
    
    ## Attempt to run an intersection to set up order of fires and areas of overlap
    ## If topology errors crop up, rebuild from rasters, slower but less error prone
    
    ## Make index table of jday
    i_tab <- data.frame(index = 1:nrow(f_day),
                        yrdy = f_day$yrdy,
                        year = f_day$year,
                        jday = f_day$jday,
                        crad_day = f_day$crad_day)
    
    f_int <- try(st_intersection(f_day), silent = T)
    
    if("try-error" %in% class(f_int)) {
      warning("Topology issues, rebuilding features (slower)")
    }
    
    if("try-error" %in% class(f_int)) {
      
      ## Convert to raster and then back again
      ## orders of magnitude faster with stars package than raster
      ## https://r-spatial.github.io/stars/articles/stars5.html
      f_day2  <- lapply(1:nrow(f_day), function(x) 
        st_rasterize(f_day[x,"yrdy"])) %>% 
        lapply(function(x) st_as_sf(x, merge = T))
      f_day2 <- do.call(rbind, f_day2) %>% 
        ## make multipolygons
        group_by(yrdy) %>% 
        summarise() %>% 
        ungroup() %>% 
        ## Add year and day info back in
        merge(i_tab, by = "yrdy") 
      
      ## still getting some ring self-interseciont issues, buffer 0 first
      f_int <- st_transform(f_day2, 3310) %>% 
        st_buffer(0) %>% 
        st_transform(crs = st_crs(landscape)) %>% 
        st_set_precision(1e5) %>% 
        # st_make_valid() %>% 
        st_intersection()
      
      
    }
    
    ## map years to origin index and subsequent calculations
    d <- mutate(f_int,
                ## lists of cos radian days, julian days, and years
                crad_days = map(origins, function(x) i_tab[x,"crad_day"]),
                jdays = map(origins, function(x) i_tab[x,"jday"]),
                years = map(origins, function(x) i_tab[x,"year"]),
                ## Set order of events 
                order = map(crad_days, function(x) length(x):1),
                ## calculate mean return interval using specified decay rate = [0,1)
                season_w = map2(crad_days, order, function(x, y) weighted.mean(x, (1 - decay_rate)^(y-1))),
                ## unlist 
                season_w = unlist(season_w), ## Convert back to julian day
                jday_w = acos(season_w),
                jday_w = jday_w * 180 / pi,
                jday_w = jday_w / 360 * 366)
    
    ## Pull out polygons (i.e. drop linesstrings from collections)
    d2 <- st_collection_extract(d, "POLYGON")
    
    sea_r <- fasterize(d2, noburn_r, field = "season_w") %>%
      mosaic(noburn_r, fun = max) %>% 
      ## limit to landscape
      mask(landscape)
    
    writeRaster(sea_r, filename = fn_r, overwrite = T)
    
    out <- basename(fn_r)
    
  }
  
  return(out)
  
}