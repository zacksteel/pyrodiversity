Fire Perimeter Preparation
================
Zack Steel

This document gives an example of how to prepare fire perimeters to be
used in [Google Earth Engine](https://tinyurl.com/CBImodel) with the
[Parks et al. 2019](https://www.mdpi.com/2072-4292/11/14/1735) model to
estimate fire severity. I will use California’s
[FRAP](https://frap.fire.ca.gov/mapping/gis-data/) interagency fire
perimeter database which includes fires down to 10 acres in size (4 ha)
back to 1950 (although we are only using fires back to 1985). I’ll focus
on just the Headwaters of the Merced River watershed of Yosemite.

``` r
library(tidyverse)
library(sf)

## Read in some perimeters - pre-filtered somewhat (central Sierra fires; >1984)
pers <- read_sf("data/spatial/firep.shp") %>% 
  mutate(Fire_Year = as.integer(YEAR_)) %>% 
  dplyr::select(FIRE_NAME, OBJECTID, Fire_Year)

## Get Merced HUC10
huc <- read_sf("data/spatial/yose_sheds.shp") %>% 
  filter(Name == "Headwaters Merced River") %>% 
  ## match coordinate system
  st_transform(crs = st_crs(pers)) %>% 
  dplyr::select(Name)

## Limit extent to the Upper Merced
keep <- st_intersects(pers, huc) %>% 
  apply(1, any)
pers1 <- pers[keep,]

## take a look
# mapview(huc, col.regions = "darkgreen") + mapview(pers1)
ggplot() +
  geom_sf(data = huc, fill = "darkgreen", color = NA, alpha = 0.3) +
  geom_sf(data = pers1, fill = "grey30", color = NA, alpha = 0.4) +
  theme_void()
```

![](PerimeterPrep_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Make the perimeter attributes match what the Park’s model expects. Make
sure Fire\_IDs are unique.

``` r
## Make names match Parks model
out <- mutate(pers1,
              state = "CA")
## Region-specific cut-offs for image search, see Parks et al. 2019 Table 2
dates <- data.frame(
  state = c("AZ", "CA", "CO", "ID", "KS", "MT", "ND", "NE", "NM", "NV", 
            "OK", "OR", "SD", "TX", "UT", "WA", "WY"),
  Start_Day = as.integer(c(91, 152, 152, 152, 121, 152, 152, 121, 91, 91,
                121, 152, 152, 121, 152, 152, 152)),
  End_Day = as.integer(c(181, 258, 258, 258, 212, 258, 258, 212, 181, 181,
              212, 218, 218, 212, 218, 218, 218)))
out2 <- merge(out, dates, by = "state") %>% 
  mutate(ID = 1:n()) %>% 
  dplyr::select(ID, Fire_Name = FIRE_NAME, Fire_ID = OBJECTID, Fire_Year, Start_Day, End_Day)
```

``` r
st_write(out2, "data/spatial/gee_pers.shp", append = F)
```

If modifying for your own fire perimeters, you’ll next need to upload
gee\_pers.shp (or whatever you name it) as an asset to your earth engine
account. Follow instructures for the [Parks CBI
model](https://tinyurl.com/CBImodel). You will likely need to modify
path to Asset file (line \~30), and path for outputs (line \~375). Also,
GEE doesn’t seem to understand subdirectories in Google Drive (or I
can’t figure them out), rather it looks for a unique base directory.
E.g. you just need to point to “dir”, regardless of whether it’s in “My
Drive/dir” or “My Drive/one/two/three/dir”. I’m not sure how it deals
with duplicate directory names.
