## Purpose: Calculate focal (with moving window) landscape functional diversity
## Project: pyrodiv_trends
## Upstream:
## Downstream:

focal_fd = function(traits, #a list of rasters (SpatRaster) with the same extent and resolution
                    # tr_names = NULL, #character vector of trait names
                    points, #optional point vector file to sample around, #### not yet implemented ####
                    w, #odd number window size (rows/columns), passed to terra::focalValues
                    metric = "FDis", #character vector of FD metrics to return. Options: 'nbsp', 'FRic', 'FEve', 'FDis'
                    pca_axes = "max" #number of PC dimensions to use when calculating FRic
                    )
{
  library(tidyverse)
  library(terra)
  library(FD)
  
  if(class(traits) != 'list') {
    traits = list(traits)
  }
  
  #### add if names are null some generic naming
  #### add option to extract using point layer instead of full scene
  #### add option to pass either a list or multilayer SpatRaster file
  #### add optional internal buffer to avoid dropping edge NA cells?
  
  ## add arbitrary names
  tr_names = sapply(1:length(traits), function(x) paste0('trait', x))
  
  ## slightly inelegant way to using multiple layers
  ctab = tibble(name = tr_names, rast = traits) %>% 
    mutate(ctab = purrr::map2(rast, name, ~focalValues(.x, w = w) %>% 
                                as.data.frame() %>% 
                                rownames_to_column(var = "com") %>% 
                                mutate(com = as.integer(com)) %>% 
                                pivot_longer(-com, values_to = "value") %>% 
                                mutate(trait = .y))) %>% 
    ## combine counts from each trait
    pull(ctab) %>% 
    bind_rows() %>% 
    ## somewhat slow step
    pivot_wider(id_cols = c(com, name), names_from = trait, values_from = value) %>% 
    ## remove the position ID
    dplyr::select(-name) %>%
    ## remove instances with only NA traits
    filter(if_any(2:last_col(), ~ !is.na(.))) %>% 
    ## Count up species/unique fire histories in each community
    group_by_all() %>% 
    count() %>% 
    ungroup() %>% 
    ## name communities and species
    mutate(com = paste0('com', com)) 
  
  ## Get species (unique combo of traits)
  spp = dplyr::select(ctab, -com, -n) %>% 
    unique() %>% 
    rowid_to_column() %>% 
    mutate(spp = paste0("sp", rowid)) %>% 
    dplyr::select(-rowid)
  
  ## add species number back to tibble
  ctab2 = left_join(ctab, spp, by = tr_names)
  
  ## Make the abundance matrix
  abun = pivot_wider(ctab2, id_cols = com, names_from = spp, values_from = n) %>% 
    column_to_rownames("com") %>% 
    ## replace missing (NA) values with zero
    replace(is.na(.), 0)
  
  ## make trait table
  tr = dplyr::select(ctab2, -com, -n) %>% 
    unique() %>% 
    column_to_rownames("spp")
  
  frich = ifelse('FRic' %in% c(metric), T, F)
  
  ## If rasters are very large = we have many 'communities' to consider run subsets in parallel
  abun_split = mutate(abun,
                     gp = ceiling(row_number()/10000)) %>% 
    group_split(gp, .keep = F) 
  
  library(foreach)
  library(doParallel)
  
  #### Some issues when splitting the data results in a species not occurring in any community
  #### How does this affect the outputs? Can we trick the function by adding in a dummy community then removing after?
  #### Does trimming the species list affect FDis and FRic? It would affect FEve.
  #### As we increase the number and resolution of traits this is likely to happen more and more
  
  registerDoParallel(1)  # use multicore, set to the number of our cores
  tic()
  div_l = foreach (i=1:length(abun_split)) %dopar% {
    
    library(FD)
    
    abun_sub = abun_split[i][[1]]
    
    div_sub = dbFD(tr, abun_sub, stand.x = T, 
                   corr = "cailliez",
                   calc.FRic = frich,
                   m = pca_axes,
                   calc.FGR = F,
                   calc.CWM = F,
                   calc.FDiv = F,
                   messages = F)
  }
  toc()
  
  div = purrr::map(div_l, data.frame) %>% 
    do.call(rbind, .)
  
  identical(div$FDis, is.numeric(div.full$FDis))
  
  ## remap diversity results to raster and return a multi-layer raster
  rl = purrr::map(metric, ~rast(x = traits[[1]],
                                vals = div[,.x]))
  out = do.call(c, rl)
  names(out) = metric
  return(out)
  
  ## NA value causes problems. may only be a problem if NAs for all traits
  tic()
  div <- dbFD(tr, abun, stand.x = T, 
              corr = "cailliez",
              calc.FRic = frich,
              m = pca_axes,
              calc.FGR = F,
              calc.CWM = F,
              calc.FDiv = F,
              messages = F)
  toc()
  
  div.full = div
  

  
  ## remap diversity results to raster and return a multi-layer raster
  rl = purrr::map(metric, ~rast(x = traits[[1]],
                                vals = div[[.x]]))
  out = do.call(c, rl)
  names(out) = metric
  return(out)
  
}
