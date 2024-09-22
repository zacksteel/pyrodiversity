## Purpose: Calculate functional dispersion of fire regime trait(s) 
## Project: Pyrodiversity

global_fd <- function(
  traits, #vector of paths to trait rasters or rast objects
  tr_wt = NULL, # relative weights for traits
  frich = F, #logical, whether to also calculate function richness
  pca_axes = "max", #number of PC dimensions to use when calculating FRic
  mask_lyr = NULL #optional mask layer path (e.g. remove non-flammable areas)
) {
  library(tidyverse, quietly = T, warn.conflicts = F, verbose = F)
  library(FD, quietly = T, warn.conflicts = F, verbose = F)
  library(terra, quietly = T, warn.conflicts = F, verbose = F)
  
  ## Read in rasters
  if(class(traits) != 'SpatRaster') {
    ts <- rast(traits)
  } else
  {
    ts <- traits
  }
  
  ## if null trait weights, force an equal weight
  if(is.null(tr_wt)) {tr_wt <- rep(1, nlyr(ts))}
  
  if(!is.null(mask_lyr)) {
    if(class(mask_lyr) == "SpatRaster") {m <- mask_lyr} else 
      {
      m <- rast(mask_lyr)
      }

    ## projection mask if needed
    if(!(identical(crs(ts), crs(m)) & identical(ext(ts), ext(m)))) {
      m <- suppressWarnings(project(m, ts, method = "near")) #%>%  #gives a warning about missing values, maybe because are only == 1?, function appear to work
        ## necessary to resample to align extents
        # resample(y = ts, method = "near")
    }
    ## run mask
    ts <- terra::mask(ts, m) 
  }
  
  ## get crosstabs as frequencies/abundances of unique trait combos
  if(nlyr(ts) == 1) {
    ctab <- freq(ts, digits = 3) %>% 
      as.data.frame() %>% 
      select(value, Freq = count)
    # names(ctab) <- c(names(ts), "Freq")
  } else {
    ctab <- crosstab(ts, digits = 3, long = T, useNA = T) %>% 
      rename(Freq = n)
  }
  
  ## dbFD cannot handle a trait with one or zero non-na values
  cnt.vals <- dplyr::select(ctab, -Freq) %>% 
    purrr::map( ~length(unique(na.omit(.x)))) < 2
  if(TRUE %in% cnt.vals) stop(paste0("Trait ", which(cnt.vals), " has one or zero unique non-NaN values. This breaks the functional diversity calculation."))
  
  ## remove instances with only NA traits (don't consider the last Freq column)
  ctab <- filter(ctab, if_any(1:last_col()-1, ~ !is.na(.)))
  
  # ctab <- ctab %>% 
  #   mutate(nacnt = rowSums(is.na(.))) %>% 
  #   filter(nacnt < length(.) - 2) %>% 
  #   dplyr::select(-nacnt)

  
  ## Warn user about number of unique species
  if(nrow(ctab) > 1000) {
    message(paste0(
      nrow(ctab), " unique 'species' (fire trait combinations).
  Large numbers of combinations (e.g. > 5000) can take a ridiculous amount of time to run. 
  Consider reducing precision of raster layers by rounding. 
  E.g. Frequency - 1 year increments; Patch size - 1 log HA increments; 
  Severity - 0.5 CBI increments; Seasonality - 0.1 increments."))
    }
  
  #### On small landscapes (e.g., around points) this fudge cell has a meaningful impact
  #### Need to find a better way to address this bug
  ## when only one value and NA (e.g. season layer), throws an error
  ## adding a single fudge cell works around this without inflating dispersion
  # cmean <- colMeans(ctab, na.rm = T)
  # add <- cmean + cmean/1000
  # add['Freq'] <- 1
  # 
  # ctab2 <- rbind(ctab, add)
  
  ctab2 <- ctab
  
  # add "species" names
  row.names(ctab2) <- sapply(1:nrow(ctab2), function(x) paste0("sp",x))
  
  ## convert to trait and abundance 
  abun <- rownames_to_column(ctab2) %>% 
    ## only have one 'site' when calculating for the full landscape
    dplyr::select(rowname, Freq) %>% 
    pivot_wider(names_from = rowname, values_from = Freq) %>% 
    as.matrix()
  traits2 <- dplyr::select(ctab2, -Freq)
  
  ## run diversity function
  div <- dbFD(traits2, abun, 
              w = tr_wt,
              stand.x = T, 
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