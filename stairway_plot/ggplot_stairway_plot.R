###==----- ggplot stairway plot results -----==###

## Set up workspace ----
library(tidyverse)
library(ggthemes)
library(ggrepel)

basedir <- "~/Desktop/shark_nf/"
pops <- dir(paste0(basedir, "results/stairway_plot/output"))[-3] # removing cocos which has very few samples

full_table <- data.frame()

## Create data collection function ----
for(pop in pops) {
  
  ### Assign populations to regions
  assign("region", ifelse(pop %in% c("indo", "scott", "rowley", "ningaloo", "chagos"), # list of pops in region 1
                          "West of Torres Strait", # name of region 1
                          ifelse(pop %in% c("north_gbr", "south_gbr", "coral_sea"), # list of pops in region 2
                                 "Great Barrier Reef", # name of region 2
                                 "New Caledonia"))) # name of region 3 which all other pops go into
  
  ### Assign data to populations
  assign("results", read.table(paste(basedir, "results/stairway_plot/output",
                                     pop, "stairways", paste0(pop, ".final.summary"),
                                     sep = "/"), header = TRUE)[,6:9])
  
  ### Make results table for graphing
  assign(paste0("results_table"),
         cbind(region = factor(region, levels = c("West of Torres Strait",
                                                  "Great Barrier Reef",
                                                  "New Caledonia")),
               population = as.character(pop),
               results))
  
  full_table <- rbind(full_table, results_table)
}

label_rows <- data.frame()
## Extremely messy just trying to get the job done fast
## Get row closest to x-value that should make good label point for each facet
for(a in as.character(unique(full_table$population))) {
  if(sum(+(full_table[full_table$population == a,1] == "Great Barrier Reef")) ==
     length(full_table[full_table$population == a,1])) {
    pop_subset <- full_table[full_table$population == a,]
    label_row <- pop_subset[which(abs(pop_subset[,3] - 7500) ==
                             min(abs(pop_subset[,3] - 7500)))[1],]
  }
  if(sum(+(full_table[full_table$population == a,1] == "West of Torres Strait")) ==
     length(full_table[full_table$population == a,1])) {
    pop_subset <- full_table[full_table$population == a,]
    label_row <- pop_subset[which(abs(pop_subset[,3] - 7300) ==
                                    min(abs(pop_subset[,3] - 7300)))[1],]
  }
  if(sum(+(full_table[full_table$population == a,1] == "New Caledonia")) ==
     length(full_table[full_table$population == a,1])) {
    pop_subset <- full_table[full_table$population == a,]
    label_row <- pop_subset[which(abs(pop_subset[,3] - 250) ==
                                    min(abs(pop_subset[,3] - 250)))[1],]
  }
  label_rows <- rbind(label_rows, label_row)
}

full_table <- full_table[-c(which(full_table$year < 140)),]
for(b in as.character(unique(full_table$population))) {
  pop_data <- full_table[which(full_table$population == b),]
  hundred_edge <- pop_data[1,]
  hundred_edge$year <- 140
  full_table <- rbind(full_table, hundred_edge)
}

final_labels <- label_rows[which(label_rows$population == "matthew" |
                                   label_rows$population == "walpole" |
                                   label_rows$population == "grand_astrolabe" |
                                   label_rows$population == "south_gbr" |
                                   label_rows$population == "ningaloo" |
                                   label_rows$population == "coral_sea" |
                                   label_rows$population == "north_gbr" |
                                   label_rows$population == "indo" |
                                   label_rows$population == "chagos"),]
# full_table <- full_table[order(full_table$year, decreasing = FALSE),]
# full_table <- full_table[order(full_table$population),]
# rownames(full_table) <- 1:dim(full_table)[1]

print(full_table[which(full_table$year == 100),c('population', 'Ne_median', 'Ne_2.5.', 'Ne_97.5.')])

stairway_plot <- ggplot() +
  geom_step(data = full_table, aes(x = year/1000/7,
                                   y = Ne_median/1000,
                                   group = population)) +
  geom_vline(xintercept = c(40, 90), linetype = "dotted") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(size = 8, vjust = -1),
        axis.text.y.left = element_text(size = 8, margin = margin(0,7,0,10)),
        axis.title = element_text(size = 9),
        legend.direction = "vertical",
        legend.key = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.box.margin = margin(c(-7.5,-2.5,-2.5,-2.5)),
        legend.position = c(0.15,0.355),
        legend.justification = c(0.5,0),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, face = 'italic'),
        legend.box.background = element_rect(colour = "black", size = 0.01),
        legend.background = element_blank(),
        strip.text.y = element_text(size = 8, colour = "black", angle = 270)) + 
  scale_x_continuous(trans = "log2",
                     limits = c(0.02,121), expand = c(0,0),
                     breaks = c(rep(seq(2,10,4),5)*10^sort(rep(seq(-2,2,1), 3))),
                     labels = c(rep(seq(2,10,4),5)*10^sort(rep(seq(-2,2,1), 3)))) +
  scale_y_continuous(breaks = seq(0, 200, 20)) +
  labs(x = "\nGenerations before present (thousands)",
       y = "Effective Population Size (thousands)") +
  scale_color_colorblind() +
  facet_grid(region~.) +
  geom_text(data = final_labels, size = 2.8, col = grey(0.5), show.legend = F,
            nudge_x = c(0, 0, 2.6, 0, 0, 0, 0, 0, 0),
            nudge_y = c(0, 0, 7.5, 0, -20, -20, 0, -20, 0),
            aes(x = year/1000/7,
                y = Ne_median/1000+10,
                label = c("Chagos", "Herald Cays",
                          "Grand Astrolabe",
                          "Misool", "Matthew",
                          "Ningaloo", "North GBR",
                          "South GBR", "Walpole")))
  # geom_text(data = final_labels, col = grey(0.5),
  #                 nudge_y = ifelse(final_labels$population == "indo", 1000, 0)
  #                                  # ifelse(final_labels$population == "walpole", 1000, 0)),
  #                 # nudge_x = ifelse(final_labels$population == "grand_astrolabe", 1000, 0),
  #                 aes(x = year/1000/7,
  #                     y = Ne_median/1000,
  #                     label = c("Chagos", "Herald Cays",
  #                               "Misool", "Matthew",
  #                               "Ningaloo", "North GBR",
  #                               "South GBR", "Walpole")),
  #                 size = 2.8, min.segment.length = 3/4,
  #                 point.padding = 1/3, box.padding = 1/4,
  #                 show.legend = FALSE) #+
  # geom_text(data = label_rows[which(label_rows$population == "grand_astrolabe"),],
  #           aes(x = year/1000/7, y = Ne_median/1000, label = "Grand Astrolabe"),
  #               size = 2.8, show.legend = FALSE, col = grey(0.5))

ggsave(paste0(basedir, 'results/stairway_plot/stairway_plot_all_pops.pdf'),
       plot = stairway_plot,
       width = 6.5,
       height = 6,
       unit = "in",
       device = "pdf")
