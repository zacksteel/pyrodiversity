## Purpose: Calculates functional diversity in a circular landscape around points
## Project: pyrodiversity

point_fd <- function(
  points, # sf or vect object of points
  traits, #vector of paths to trait rasters or rast objects
  tr_wt = NULL, # relative weights for traits
  land_radius, # radius of the local landscape to consider for functional diversity, units are those of the length unit of the trait crs
  frich = F, #logical, whether to also calculate function richness
  pca_axes = "max", #number of PC dimensions to use when calculating FRic
  mask = NULL #optional mask layer path (e.g. remove non-flammable areas)
  )
{
  ## read in necessary libraries and global_fd function
  library(tidyverse, quietly = T, warn.conflicts = F, verbose = F)
  library(sf, quietly = T, warn.conflicts = F, verbose = F)
  library(terra, quietly = T, warn.conflicts = F, verbose = F)
  
  #### sensitive to location of this function; need to fix if making into a package
  source(here("code/global_fd.R"))
  
  ## Read in rasters
  if(class(traits) != 'SpatRaster') {
    ts <- rast(traits)
  } else
  {
    ts <- traits
  }
  
  ## make sure points are in vect format and the same crs as traits
  if(!("SpatVector" %in% class(points))) {
    points <- vect(points)
  }
  
  if(crs(points) != crs(traits)) {
    points <- project(points, crs(traits)[1])
  }
  
  ## buffer the points per the user specified width and shape
  pts.buf <- buffer(points, width = land_radius)
  
  ## complete the rest on a point by point basis
  d <- data.frame()
  
  for(i in 1:length(pts.buf)) {
    pt <- pts.buf[i,]
    
    ts.sm <- crop(ts, pt) %>% 
      # touches = F means that a cell whose center is not within the circle is excluded
      mask(pt, touches = F)
    
    ## if no variation left in a trait drop it
    mm <- minmax(ts.sm)
    mm <- mm[2,] - mm[1,]
    keep <- ifelse(mm == 0, FALSE, TRUE)
    
    ## if none left move on
    if(TRUE %in% keep) {
      ts.sm <- ts.sm[[keep]]
      ## same for trait weights if not null
      if(!is.null(tr_wt)) {tr_wt <- tr_wt[keep]}
      
      nts <- nlyr(ts.sm)
      
      ## Run through global_fd
      d.sub <- global_fd(traits = ts.sm,
                         tr_wt = tr_wt,
                         frich = frich,
                         pca_axes = pca_axes,
                         mask = mask)
      d.sub$radius <- land_radius
      d.sub$ntraits <- nlyr(ts.sm)
      
      ## add to the running list
      d <- bind_rows(d, d.sub)
    } else {
      d.sub <- data.frame(radius = land_radius,
                          ntraits = 0)
      d <- bind_rows(d, d.sub)
    }
    
  }
  
  return(d)
}