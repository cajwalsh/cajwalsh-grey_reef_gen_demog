#####==---- Wide-scale geographic analysis of genetic diversity ----==##### -

## Load and properly organize genetic and population location data ----
# setwd("~/Desktop/shark_nf/results/diversity")

### Diversity statistics
div_per_site <- read.table("div_per_site.txt", header = T)
colnames(div_per_site)[1] <- "population"
div_per_site$population <- c("Chagos", "Chesterfield", "Cocos", "Coral_Sea", "Entrecasteaux",
                             "GLN", "Grand_Astrolabe", "Indonesia", "Matthew", "Ningaloo",
                             "North_GBR", "Petit_Astrolabe", "Rowley", "Scott", "South_GBR",
                             "South_NC", "Walpole")

### Population coordinates
loc <- read.csv("../../data/coordinates.csv")
names(loc) <- c("population", "x", "y")

#### Edit strange thing that never popped up before saying the coords
#### for Indonesia population are above sea level but don't seem to be
#### irl (probably linked to some glitch in NOAA data, its new portal,
#### or corresponding marmap update)
loc[6,2] <- 130.5
loc[6,3] <- -2

diversity_pop_coords <- merge(div_per_site, loc, by = "population")

### Assign analysis region extremities
xwest <- 70
xeast <- 180
ynorth <- 26
ysouth <- -26

## Load and properly organize cartographic & oceanographic data ----
library(marmap)
bat <- getNOAA.bathy(xwest, xeast, ynorth, ysouth, resolution = 4)

### Create geographical point network and only keep points between -1m and -200m elevation
xcoords <- rep(seq(xwest+0.25,xeast-0.25, 0.5), length(seq(ynorth-0.25, ysouth+0.25, -0.5)))
xcoords <- xcoords[order(xcoords)]
ycoords <- rep(seq(ynorth-0.25, ysouth+0.25, -0.5), length(seq(xwest+0.25,xeast-0.25, 0.5)))
all_coords <- matrix(data = c(xcoords,ycoords), ncol = 2, dimnames = list(NULL, c("x", "y")))
filt_coords <- all_coords[get.depth(bat, locator = FALSE, all_coords[,1:2])[,3] <= -1 &
                        get.depth(bat, locator = FALSE, all_coords[,1:2])[,3] >= -200,]
rownames(filt_coords) <- 1:dim(filt_coords)[1]
filt_coords_list <- split(filt_coords, 1:nrow(filt_coords))

### Check map and bathymetry
# plot(bat, image = TRUE, drawlabels = FALSE, lwd = 0.1,
#      deep = -200, shallow = 0, step = 200, land = TRUE,
#      bpal = list(c(min(bat, na.rm = TRUE), -200, c("lightsteelblue4", "lightsteelblue2")),
#                  c(-200, 0, "lightsteelblue1"),
#                  c(0, max(bat, na.rm = TRUE), grey(0.65))))

### Check point filtering
# library(sp)
# sp_coords <- SpatialPoints(filt_coords, proj4string =
#                              CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# points(sp_coords, col = "deep pink 3", cex = 0.005, pch = 21)

### Make transition matrix
print("Making transition matrix")
mardist_matrix <- trans.mat(bat, min.depth = -1)

## Analyze correlation between nucleotide diversity and geographic distance ----
library(gdistance)
div_cor <- function(coord_vec) {
  cor(diversity_pop_coords$pi_per_site,
      matrix(costDistance(mardist_matrix,
                          matrix(coord_vec, ncol = 2),
                          as.matrix(diversity_pop_coords[,c('x', 'y')])), ncol = 1),
      method = "spearman")
}

start_time <- Sys.time()
print("Correlation analysis begun")

r_list <- lapply(filt_coords_list,
                 div_cor)

Sys.time() - start_time # took 3 hours 36 minutes

wide_r_div_coords <- cbind(filt_coords, unlist(r_list))

### Save results matrix as Rdata
save(wide_r_div_coords, file = "wide_r_div_coords.Rdata")

### Check plotting (load just results in to do it)
# load("wide_r_div_coords.Rdata")

### Reorganize bathymetric and results data for plotting
bat_df <- as.data.frame(as.SpatialGridDataFrame(bat))
colnames(bat_df) <- c("z", "x", "y")

r_df <- merge(all_coords, wide_r_div_coords,
              by = c("x", "y"), all = TRUE)
colnames(r_df) <- c("x", "y", "rho")

## Plot and save non-projected raster results map ## ----
library(ggplot2)
library(ggnewscale)
library(scales)
wide_div_map <- ggplot() +

  ### Create rho results raster and fill gradient
  geom_raster(data = r_df, aes(x = x, y = y, fill = rho)) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(7, "PuOr"),
                       na.value = "lightsteelblue2",
                       limits = c(-1, 1),
                       name = "rho",
                       breaks = c(round(min(r_df$rho, na.rm = TRUE), digits = 3), round(0),
                                  round(max(r_df$rho, na.rm = TRUE), digits = 3))) +

  ### Make room for another fill gradient on same graph for next geom
  new_scale("fill") +

  ### Create bathymetric map raster and fill gradient
  geom_raster(data = bat_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(colours = c("lightsteelblue4", "lightsteelblue2",
                                   "transparent", "transparent",
                                   grey(2/3), grey(2/3)),
                       values = rescale(c(min(bat, na.rm = TRUE),
                                          -200.0001, -200, 0, 0.0001,
                                          max(bat, na.rm = TRUE))),
                       guide = FALSE) +

  ### Add sampling point locations
  geom_point(data = loc, aes(x = x, y = y),
             size = 2, shape = 3, stroke = 3/4) +

  ### Add points for most negative location(s)
  geom_point(data = r_df[which(r_df$rho == min(r_df$rho, na.rm = T)),],
             aes(x = x, y = y),
             size = 2, stroke = 3/4,
             shape = 4, color = "white") +

  ### Set theme parameters
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 8),
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        legend.text = element_text(size = 8),
        legend.position = c(0.925, 0.75))

### Save plot as .pdf and .tif
ggsave("wide_div_map.pdf",
       plot = wide_div_map, device = "pdf",
       width = 6.5,
       height = 6.5/abs(diff(c(xwest, xeast))/diff(c(ysouth, ynorth))),
       units = "in")

ggsave("wide_div_map.tif",
       plot = wide_div_map, device = "tiff",
       dpi = 800,
       width = 6.5,
       height = 6.5/abs(diff(c(xwest, xeast))/diff(c(ysouth, ynorth))),
       units = "in")
