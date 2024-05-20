## Purpose: Calculate focal (with moving window) landscape functional diversity
## Project: pyrodiversity

focal_fd = function(traits, #a list of rasters (SpatRaster) with the same extent and resolution
                    tr_wt = NULL, # relative weights for traits if using the FD package
                    # w, #odd number window size (rows/columns), passed to terra::focalValues
                    metric = "FDis", #character vector of FD metrics to return. Options: 'nbsp', 'FRic', 'FEve', 'FDis'
                    land_radius, # radius of the local landscape to consider for functional diversity, units are those of the length unit of the trait crs (usually meters)
                    frich = F, #logical, whether to also calculate function richness
                    pca_axes = "max", #number of PC dimensions to use when calculating FRic
                    mask_lyr = NULL, #optional mask layer path (e.g. remove non-flammable areas)
                    proot = getwd(), #optional project root needed to load point_fd and global_fd.R; temporary work around until this is a formal package
                    cores = 1, #number of cores to use for parallel processing
                    batch_size = 5000, #number of cells to process in each batch if cores > 1
                    # method = "fundiversity",
                    # na_replacement = T, #fundiversity will drop species that contain NAs; if TRUE na's will be replaced with the trait means
                    out_raster = NULL
                    )
{
  library(tidyverse)
  library(terra)
  library(furrr)
  source(paste0(proot, "/code/point_fd.R"))
  # library(FD)
  
  print("focal_fd is still under development, use with extreme caution! Multi-core processing is not yet implemented.")
  
  ## create points of each raster cell with an internal buffer equal to the land_radius
  ext.inner <- vect(ext(traits[[1]]), crs = crs(traits[[1]]))
  ext.inner <- buffer(ext.inner, -land_radius)
  r.inner <- crop(traits[[1]], ext.inner)
  
  pts <- as.points(crop(traits[[1]], r.inner), values = F, na.rm = F)
  
  ## check to make sure the number of points matches the number of pixels in the inner raster
  if(nrow(pts) != ncell(r.inner)) {
    stop("Number of points does not match number of pixels in inner raster")
  }
  
  ## if running in sequence just get to it
  if(cores == 1) {

    ## run through point_fd for each point
    pd <- point_fd(
      points = pts, # sf or vect object of points
      traits = traits, #vector of paths to trait rasters or rast objects
      tr_wt = tr_wt, # relative weights for traits
      land_radius = land_radius, # radius of the local landscape to consider for functional diversity, units are those of the length unit of the trait crs (usually meters)
      frich = frich, #logical, whether to also calculate function richness
      pca_axes = pca_axes, #number of PC dimensions to use when calculating FRic
      mask_lyr = mask_lyr, #optional mask layer path (e.g. remove non-flammable areas)
      proot = proot #optional project root needed to load global_fd.R; temporary work around until this is a formal package
    )
    
  }
  
  ## if running in parallel, first split up the landscape into batches/grid cells
  if(cores > 1) {
    
    ## can't export rasters so save a temp file to read in for each worker
    tmp_fl <- "traitrast.tif"
    writeRaster(traits, filename = tmp_fl, overwrite = T)

    ## split points into batches
    pts_split <- split(pts, ceiling(seq_along(pts)/batch_size)) %>% 
      ## can't pass vect objects through future_map either, use sf instead
      purrr::map(~st_as_sf(.x))
    
    ## run through point_fd for each batch using furrr
    ## set workers to number of cores selected or the number of batches, whichever is smaller
    workers <- ifelse(cores < length(pts_split), cores, length(pts_split))
    plan(multisession, workers = workers)
    
    pd.l <- future_map(pts_split, 
                       ~point_fd(
                         points = .x, # sf or vect object of points
                         traits = tmp_fl, #vector of paths to trait rasters or rast objects
                         tr_wt = tr_wt, # relative weights for traits
                         land_radius = land_radius, # radius of the local landscape to consider for functional diversity, units are those of the length unit of the trait crs (usually meters)
                         frich = frich, #logical, whether to also calculate function richness
                         pca_axes = pca_axes, #number of PC dimensions to use when calculating FRic
                         mask_lyr = mask_lyr, #optional mask layer path (e.g. remove non-flammable areas)
                         proot = proot #optional project root needed to load global_fd.R; temporary work around until this is a formal package
    ),
    options(future.globals.onReference = "error"))
    plan(sequential)
    
    ## delete temp directory
    unlink(tmp_fl)
    
    ## bind results
    pd <- bind_rows(pd.l)

  }
  
  ## map results to raster
  ## remap diversity results to raster and return a multi-layer raster
  ## if ntraits are all zero in pd, set vals to NA
  if(sum(pd[,"ntraits"]) == 0) {
    rl = purrr::map(metric, ~rast(x = r.inner,
                                  vals = NA))
  } else {
    rl = purrr::map(metric, ~rast(x = r.inner,
                                  vals = pd[,.x]))
  }
  
  out = do.call(c, rl)
  names(out) = metric

  if(!is.null(out_raster)) {
    writeRaster(out, filename = out_raster, overwrite = T)
  }

  return(out)
  
  # ## fundiversity functions cannot weight traits explicitly (test indirect method via abundance)
  # if(!is.null(tr_wt) & method == "fundiversity") stop("fundiversity functions cannot weight traits explicitly. Remove weighting ('tr_wt' argument) or use 'FD' method")
  # 
  # if(class(traits) != 'list') {
  #   traits = list(traits)
  # }
  # 
  # #### add option to pass either a list or multilayer SpatRaster file
  # 
  # ## add arbitrary names
  # tr_names = sapply(1:length(traits), function(x) paste0('trait', x))
  # 
  # ## slightly inelegant way to using multiple layers
  # ctab = tibble(name = tr_names, rast = traits) %>% 
  #   #focalValues returns a row for each cell ("community") w/ columns representing all neighbor values within window
  #   mutate(ctab = purrr::map2(rast, name, ~focalValues(.x, w = w) %>% 
  #                               as.data.frame() %>% 
  #                               rownames_to_column(var = "com") %>% 
  #                               mutate(com = as.integer(com)) %>% 
  #                               ## make long but hold onto com/cell values
  #                               pivot_longer(-com, values_to = "value") %>% 
  #                               mutate(trait = .y))) %>% 
  #   ## combine counts from each trait
  #   pull(ctab) %>% 
  #   bind_rows() %>% 
  #   ## somewhat slow step
  #   pivot_wider(id_cols = c(com, name), names_from = trait, values_from = value) %>% 
  #   ## remove the position ID
  #   dplyr::select(-name) %>%
  #   ## remove instances with only NA traits
  #   filter(if_any(2:last_col(), ~ !is.na(.))) %>%
  #   ## Count up species/unique fire histories in each community
  #   group_by_all() %>% 
  #   count() %>% 
  #   ungroup() %>% 
  #   ## name communities and species
  #   mutate(com = paste0('com', com)) 
  # 
  # ## Get species (unique combo of traits)
  # spp = dplyr::select(ctab, -com, -n) %>% 
  #   unique() %>% 
  #   rowid_to_column() %>% 
  #   mutate(spp = paste0("sp", rowid)) %>% 
  #   dplyr::select(-rowid)
  # 
  # ## add species number back to tibble
  # ctab2 = left_join(ctab, spp, by = tr_names)
  # 
  # ## Make the abundance matrix
  # abun = pivot_wider(ctab2, id_cols = com, names_from = spp, values_from = n) %>% 
  #   column_to_rownames("com") %>% 
  #   ## replace missing (NA) values with zero
  #   replace(is.na(.), 0)
  # 
  # ## make trait table
  # tr = dplyr::select(ctab2, -com, -n) %>% 
  #   unique() %>% 
  #   column_to_rownames("spp")
  # 
  # frich = ifelse('FRic' %in% c(metric), T, F)
  # 
  # ## If rasters are very large = we have many 'communities' to consider run subsets in parallel
  # ## not yet implimented
  # abun_split = mutate(abun,
  #                    gp = ceiling(row_number()/10000)) %>%
  #   rownames_to_column() %>% 
  #   group_split(gp, .keep = F) 
  # ## need rownames for fd_dis function
  # abun_split <- purrr::map(abun_split, ~column_to_rownames(.x))
  # 
  # # library(foreach)
  # # library(doParallel)
  # # 
  # # #### Some issues when splitting the data results in a species not occurring in any community
  # # #### How does this affect the outputs? Can we trick the function by adding in a dummy community then removing after?
  # # #### Does trimming the species list affect FDis and FRic? It would affect FEve.
  # # #### As we increase the number and resolution of traits this is likely to happen more and more
  # # 
  # # registerDoParallel(1)  # use multicore, set to the number of our cores
  # # tic()
  # # div_l = foreach (i=1:length(abun_split)) %dopar% {
  # #   
  # #   library(FD)
  # #   
  # #   abun_sub = abun_split[i][[1]]
  # #   
  # #   div_sub = dbFD(tr, abun_sub, stand.x = T, 
  # #                  corr = "cailliez",
  # #                  calc.FRic = frich,
  # #                  m = pca_axes,
  # #                  calc.FGR = F,
  # #                  calc.CWM = F,
  # #                  calc.FDiv = F,
  # #                  messages = F)
  # # }
  # # toc()
  # # 
  # # div = purrr::map(div_l, data.frame) %>% 
  # #   do.call(rbind, .)
  # # 
  # # 
  # # ## remap diversity results to raster and return a multi-layer raster
  # # rl = purrr::map(metric, ~rast(x = traits[[1]],
  # #                               vals = div[,.x]))
  # # out = do.call(c, rl)
  # # names(out) = metric
  # # return(out)
  # 
  # 
  # 
  # 
  # #### Testing out the Fundiversity version
  # if(method == "fundiversity")
  # {
  #   library(fundiversity)
  #   ## need to scale trait values first
  #   scale2 <- function(x, na.rm = TRUE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
  #   
  #   tr2 <- mutate_all(tr, scale2)
  #   
  #   ## What to do with NAs?
  #   if(na_replacement) {
  #     tr2 <- mutate_all(tr2, replace_na, replace = 0)
  #   }
  #   
  #   # future::plan(future::multisession, workers = 1)
  #   ## for large landscapes split apart so as not to max out your RAM
  #   if(length(abun_split) > 1) {
  #     dis <- purrr::map(abun_split, ~fd_fdis(tr2, as.matrix(.x))) %>% 
  #       bind_rows()
  #   } else {
  #     dis <- fd_fdis(tr2, as.matrix(abun))
  #   }
  # 
  #   ## Need to back-fill NAs before mapping to raster
  #   dis2 <- mutate(dis, id = as.integer(str_sub(site, 4)))
  #   dummy.fill <- data.frame(id = seq(1, ncell(traits[[1]])))
  #   dis2 <- full_join(dis2, dummy.fill) %>% 
  #     arrange(id)
  #   
  #   ## map values on
  #   rl = purrr::map(metric, ~rast(x = traits[[1]],
  #                                 vals = dis2[.x]))
  # } 
  # 
  # if(method %in% c("fd", "FD")) 
  #   {
  #   ## NA value causes problems. may only be a problem if NAs for all traits
  #   div <- dbFD(tr, abun, 
  #               w = tr_wt,
  #               stand.x = T, 
  #               corr = "cailliez",
  #               calc.FRic = frich,
  #               m = pca_axes,
  #               calc.FGR = F,
  #               calc.CWM = F,
  #               calc.FDiv = F,
  #               messages = F)
  #   
  #   ## remap diversity results to raster and return a multi-layer raster
  #   rl = purrr::map(metric, ~rast(x = traits[[1]],
  #                                 vals = div[[.x]]))
  # }
  # 
  # out = do.call(c, rl)
  # names(out) = metric
  # 
  # if(!is.null(out_raster)) {
  #   writeRaster(out, filename = out_raster, overwrite = T)
  # }
  # 
  # return(out)
  
  
  
}
