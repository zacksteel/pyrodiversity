## Purpose: Calculate functional dispersion of fire regime trait(s) 
## Project: Pyrodiversity

pyrodiv_calc <- function(
  traits, #vector of paths to trait rasters or rasters themselves
  frich = F, #logical, whether to also calculate function richness
  pca_axes = "max", #number of PC dimensions to use when calculating FRic
  mask = NULL #optional mask layer path (e.g. remove non-flammable areas)
) {
  library(tidyverse)
  library(FD)
  library(raster)
  
  ## Read in rasters
  ts <- stack(traits)
  if(!is.null(mask)) {
    m <- if(class(mask) == "RasterLayer") {m <- mask} else {
      m <- raster(mask)
    }

    ## projection mask if needed
    if(!(identical(crs(ts), crs(m)) & identical(extent(ts), extent(m)))) {
      m <- suppressWarnings(projectRaster(from = m, crs = crs(ts), res = res(ts), #gives a warning about missing values, maybe because are only == 1?, function appear to work
                             method = "ngb", alignOnly = F)) %>% 
        ## necessary to resample to align extents
        resample(y = ts, method = "ngb")
    }
    ## run mask
    ts <- mask(ts, m) 
  }
  
  ## get crosstabs as frequencies/abundances of unique trait combos
  if(nlayers(ts) == 1) {
    ctab <- freq(ts, digits = 2) %>% 
      as.data.frame()
    names(ctab) <- c(names(ts), "Freq")
  } else {
    ctab <- crosstab(ts, digits = 3, long = T, useNA = T) 
  }
  
  ## drop all NA row
  ctab <- ctab %>% 
    mutate(nacnt = rowSums(is.na(.))) %>% 
    filter(nacnt < length(.) - 2) %>% 
    dplyr::select(-nacnt)

  
  ## Warn user about number of unique species
  if(nrow(ctab) > 1000) {
    message(paste0(
      nrow(ctab), " unique 'species' (fire trait combinations).
  Large numbers of combinations (e.g. > 5000) can take a ridiculous amount of time to run. 
  Consider reducing precision of raster layers by rounding. 
  E.g. Frequency - 1 year increments; Patch size - 1 log HA increments; 
  Severity - 0.5 CBI increments; Seasonality - 0.1 increments."))
    }
  
  ## when only one value and NA (e.g. season layer), throws an error
  ## adding a single fudge cell works around this without inflating dispersion
  cmean <- colMeans(ctab, na.rm = T)
  add <- cmean + cmean/1000
  add['Freq'] <- 1

  ctab <- rbind(ctab, add)
  
  # add "species" names
  row.names(ctab) <- sapply(1:nrow(ctab), function(x) paste0("sp",x))
  
  ## convert to trait and abundance 
  abun <- rownames_to_column(ctab) %>% 
    pivot_wider(names_from = rowname, id_cols = Freq, values_from = Freq) %>% 
    as.matrix()
  traits <- dplyr::select(ctab, -Freq)
  
  ## run diversity function
  div <- dbFD(traits, abun, stand.x = T, 
              corr = "cailliez",
              calc.FRic = frich, 
              m = pca_axes,
              calc.FGR = F,
              calc.CWM = F,
              calc.FDiv = F,
              messages = F)
  
  ## Also calculate Simpsons
  div$simpson <- diversity(abun, index = "simpson")
  
  ## add to adataframe
  d <- as.data.frame(div)
  
  return(d)
  
}