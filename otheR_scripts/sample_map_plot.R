#####==---- Plotting Figure 1 maps ----==##### -

## Loading world map below can take some time, so I saved it to my desktop before.
## Sub in the commented part for where I load "world_map.Rdata"

setwd("~/Desktop/shark_nf/")

## Define map areas ----

### Load population location data
loc <- read.csv("data/coordinates.csv")
names(loc) <- c("population", "x", "y")
#### Change population names on maps
levels(loc$population)[c(3, 4, 6, 8, 16)] <- c("Cocos (Keeling)", "Herald Cays", "Northern Lagoon",
                                            "Misool", "Southern Lagoon")
N <- c(22, 5, 23, 24, 24, 23, 19, 8, 21, 40, 54, 51, 46, 51, 39, 27, 36)
loc <- cbind(loc, N)

library(tidyverse)
loc$population <- gsub("_", " ", loc$population)
loc <- mutate(loc, isNC = ifelse(x < 152, FALSE, TRUE))
loc <- mutate(loc, isIndOc = ifelse(x < 100, TRUE, FALSE))

### New Caledonia map limits
NCxeast <- loc[loc$population == "Chesterfield", "x"]-1.14
NCxwest <- loc[loc$population == "Matthew", "x"]+1.14
NCynorth <- loc[loc$population == "Entrecasteaux", "y"]+1.56
NCysouth <- loc[loc$population == "Walpole", "y"]-1.56

### Torres Strait LGM map limits
Txwest <- loc[loc$population == "Misool", "x"]-25
Txeast <- loc[loc$population == "Misool", "x"]+25
Tynorth <- 0
Tysouth <- -25

### Load bathymetric data for...
library(marmap)

#### ...world map
load("~/Desktop/reference/world_map.Rdata") ## getNOAA.bathy(-180, 180, -90, 90, resolution = 4)

#### ...New Caledonia map
NC_bat <- getNOAA.bathy(NCxwest, NCxeast, NCysouth, NCynorth, resolution = 1)

#### ...Torres Strait LGM map
T_bat <- getNOAA.bathy(Txwest, Txeast, Tysouth, Tynorth, resolution = 4)

### Change projections...
library(sp) # for sample points
library(raster) # for bathymetric data

#### ...for sample map points...
ortho_project_points <- function(points_df, xmid, ymid) {
  projd_points <- points_df
  coordinates(projd_points) <- c("x", "y")
  proj4string(projd_points) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  projd_points <- spTransform(projd_points, CRS(paste0("+proj=ortho +lat_0=",
                                                       ymid, " +lon_0=", xmid,
                                                       " +x_0=0 +y_0=0")))
  return(projd_points)
}

##### ...for sampling locations in world map
projd_loc <- ortho_project_points(loc,
                                  mean(c(min(loc$x), max(loc$x))),
                                  mean(c(min(loc$y), max(loc$y))))

##### ...for sampling locations in NC
projd_loc <- ortho_project_points(loc,
                                  mean(c(min(loc$x), max(loc$x))),
                                  mean(c(min(loc$y), max(loc$y))))

##### ...for Torres Strait LGM sample map box
projd_T_box <- ortho_project_points(data.frame(x = c(Txwest, Txeast),
                                               y = c(Tysouth, Tynorth)),
                                    mean(c(min(loc$x), max(loc$x))),
                                    mean(c(min(loc$y), max(loc$y))))

##### ...for New Caledonia sample map box
projd_NC_box <- ortho_project_points(data.frame(x = c(NCxwest, NCxeast),
                                                y = c(NCysouth, NCynorth)),
                                     mean(c(min(loc$x), max(loc$x))),
                                     mean(c(min(loc$y), max(loc$y))))

#### ...for bathymetric data...
ortho_project_bathy <- function(bathy, xmid, ymid) {
  ras <- marmap::as.raster(bathy)
  ras_proj <- projectRaster(ras,
                            crs = paste0("+proj=ortho +lat_0=",
                                         ymid, " +lon_0=", xmid,
                                         " +x_0=0 +y_0=0"))
  bathy_proj <- as.bathy(ras_proj)
  return(bathy_proj)
}

##### ...world map
world_bat_projd <- ortho_project_bathy(world_map,
                                       mean(c(min(loc$x), max(loc$x))),
                                       mean(c(min(loc$y), max(loc$y))))

###### Load land shapefiles to cover weird world map bathy data
library(rnaturalearth)
land <- ne_countries(scale = "medium",
                     type = "map_units",
                     returnclass = "sf")

###### Load Antarctica shape to cover weird world map bathy data
antarctica <- ne_countries(scale = "medium",
                           type = "map_units",
                           country = "Antarctica",
                           returnclass = "sf")

###### Use some kind of magic to deal with potential infinity points
library(sf)
library(mapview)
land <- st_cast(land, 'MULTILINESTRING') %>%
  st_cast('LINESTRING', do_split = TRUE) %>%
  mutate(npts = npts(geometry, by_feature = TRUE)) %>%
  st_cast('POLYGON')
antarctica <- st_cast(antarctica, 'MULTILINESTRING') %>%
  st_cast('LINESTRING', do_split = TRUE) %>%
  mutate(npts = npts(geometry, by_feature = TRUE)) %>%
  st_cast('POLYGON')

##### ...for New Caledonia map
NC_bat_projd <- ortho_project_bathy(NC_bat,
                                    (NCxwest+NCxeast)/2,
                                    (NCysouth+NCynorth)/2)

## Plot sample map ----
proj_world_df <- fortify(world_bat_projd)
repel_df <- as.data.frame(projd_loc[which(loc$isNC == FALSE),])

library(scales)
library(ggnewscale)
library(ggrepel)

sample_map <- ggplot() +
  theme_void() +
  geom_raster(data = proj_world_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(colours = c("lightsteelblue4", "lightsteelblue3",
                                   "lightsteelblue2", "lightsteelblue2",
                                   grey(2/3), grey(2/3)),
                       values = rescale(c(min(proj_world_df$z, na.rm = TRUE),
                                          -200.0001, -200, 0,
                                          0.0001, max(proj_world_df$z, na.rm = TRUE))),
                       na.value = "transparent",
                       guide = "none") +
  geom_sf(data = land, color = "transparent", fill = grey(2/3)) +
  geom_sf(data = antarctica, color = "transparent", fill = "white") +
  geom_rect(data = as.data.frame(projd_T_box),
            aes(xmin = x[1], xmax = x[2],
                ymin = y[1], ymax = y[2]),
            color = "white", fill = "transparent",
            size = 0.4, linetype = "48") +
  geom_rect(data = as.data.frame(projd_NC_box),
            aes(xmin = x[1], xmax = x[2],
                ymin = y[1], ymax = y[2]),
            color = "white", fill = "transparent",
            size = 0.4, linetype = "48") +
  geom_point(data = repel_df, aes(x = x, y = y),
             color = "coral2", size = 2, shape = 1) +
  geom_text(size = 3.5, fontface = "bold", nudge_x = 3e5,
            aes(x = c(as.data.frame(projd_T_box)[2,1],
                      as.data.frame(projd_NC_box)[1,1]),
                y = c(as.data.frame(projd_T_box)[2,2],
                      as.data.frame(projd_NC_box)[2,2]),
                label = c("a", "c"))) +
  geom_text_repel(size = 2.8, point.padding = 1/2,
                  min.segment.length = 0, segment.size = 1/3,
                  nudge_x = ifelse(repel_df$population == "Ningaloo" |
                                     repel_df$population == "Rowley", -2.5e6,
                                   ifelse(repel_df$population == "Herald Cays" |
                                            repel_df$population == "South GBR", -1e6,
                                          ifelse(repel_df$population == "Scott", -2.5e5,
                                                 ifelse(repel_df$population == "Chagos", 1e4,
                                                        ifelse(repel_df$population == "Cocos (Keeling)" |
                                                                 repel_df$population == "North GBR", 5e5, 0))))),
                  nudge_y = ifelse(repel_df$population == "North GBR" |
                                     repel_df$population == "Cocos (Keeling)", 2e6,
                                   ifelse(repel_df$population == "Misool", 1.5e6,
                                          ifelse(repel_df$population == "Ningaloo" |
                                                   repel_df$population == "South GBR", -1e6, 0))),
                  aes(x = repel_df$x, y = repel_df$y,
                      label = paste0(repel_df$population, " (", repel_df$N, ")"))) +
  geom_text_repel(size = 2.8, #point.padding = 0,
                  min.segment.length = 0, segment.size = 1/3,
                  nudge_x = 3.75e5, nudge_y = 5e5,
                  aes(x = mean(as.data.frame(projd_NC_box)[,1]),
                      y = as.data.frame(projd_NC_box)[2,2],
                      label = paste0("New Caledonia (", sum(loc[which(loc$isNC==T), "N"]), ")"))) +
  coord_sf(crs = paste0("+proj=ortho +lat_0=", mean(c(min(loc$y), max(loc$y))),
                        " +lon_0=", mean(c(min(loc$x), max(loc$x))), " +x_0=0 +y_0=0"))
Torres_map <- ggplot() + 
  theme_void() +
  geom_raster(data = fortify(T_bat), aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(colours = c("lightsteelblue4", "lightsteelblue3",
                                   "lightsteelblue2", "lightsteelblue2",
                                   "#A4BCA4", "#A4BCA4", grey(2/3),
                                   grey(2/3)),
                       values = rescale(c(min(fortify(T_bat)$z, na.rm = TRUE),
                                          -320.0001, -320, -120.0001, -120,
                                          0, 0.0001, max(fortify(T_bat)$z, na.rm = TRUE))),
                       na.value = "transparent",
                       guide = "none") +
  geom_segment(aes(x = 142, y = -10, xend = 139, yend = -10.5), size = 3/8,
               arrow = arrow(length = unit(1.5, "mm"), ends = "first")) +
  geom_text(size = 2.9, aes(x = 137.5, y = -10.25, label = "Torres\nStrait")) +
  geom_point(data = loc[3:9,1:3], aes(x = x, y = y),
             color = "coral2", size = 2, shape = 1) +
  geom_text_repel(data = loc[3:9,1:3], point.padding = 0.25, size = 2.8,
                  nudge_x = ifelse(loc[3:9,1:3]$population == "Misool" |
                                     loc[3:9,1:3]$population == "North GBR" |
                                     loc[3:9,1:3]$population == "Ningaloo", 5, -5),
                  min.segment.length = 0, segment.size = 1/3,
                  aes(x = x, y = y, label = population))

NC_map <- ggplot() +
  theme_void() +
  geom_raster(data = fortify(NC_bat), aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(colours = c("lightsteelblue4", "lightsteelblue3",
                                   "lightsteelblue2", "lightsteelblue2",
                                   grey(2/3), grey(2/3)),
                       values = rescale(c(min(fortify(NC_bat)$z, na.rm = TRUE),
                                          -200.0001, -200, 0, 0.0001,
                                          max(fortify(NC_bat)$z, na.rm = TRUE))),
                       na.value = "transparent",
                       guide = "none") +
  geom_point(data = loc[which(loc$x > 152),], aes(x = x, y = y),
             color = "coral2", size = 2, shape = 1) +
  geom_text_repel(data = loc[which(loc$x > 152),],
                  size = 2.8, point.padding = 0.375, box.padding = 0.125,
                  min.segment.length = 0, segment.size = 1/3,
                  nudge_x = ifelse(loc[which(loc$x > 152),1] == "Entrecasteaux" |
                                     loc[which(loc$x > 152),1] == "Chesterfield", 1,
                                   ifelse(loc[which(loc$x > 152),1] == "Matthew", -1, 0)),
                  nudge_y = ifelse(loc[which(loc$x > 152),1] == "Entrecasteaux", 1, 0),
                  aes(x = x, y = y,
                      label = paste0(population, " (", N, ")")))

## Define layout
library(patchwork)
tgr_layout <- c(
  patchwork::area(l = 1, r = 12,
                  t = 1, b = 1),
  patchwork::area(l = 12, r = 23,
                  t = 1, b = 2),
  patchwork::area(l = 1, r = 12,
                  t = 2, b = 2)
)

## Organize and save together_plot ----
# ggplot() + ggplot() + ggplot() +
#   plot_layout(design = tgr_layout)

together_plot <- (Torres_map + theme(plot.tag.position = c(0.075,0.125))) +
  (sample_map + theme(plot.tag.position = c(0.875, 0.125))) +
  (NC_map + theme(plot.tag.position = c(0.075, 0.125))) +
  plot_layout(design = tgr_layout) + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 12.7, face = "bold"))

ggsave("results/sample_plot.pdf",
       plot = together_plot,
       device = "pdf",
       height = 4.8,
       width = 4.8*1.88162,
       units = "in")
