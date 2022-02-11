###==----- ggplot NeEstimator results -----==###

basedir <- "~/Desktop/shark_nf/"

## Population coordinates
loc <- read.csv(paste0(basedir, "data/coordinates.csv"))
names(loc) <- c("Population", "x", "y")
loc$Population <- sub("_", " ", loc$Population)

## Import NeEstimator results and arrange as needed
NeEstim <- read.csv(paste0(basedir, "results/NeEstimator/NeEstimator.csv"), header = T)
colnames(NeEstim)[1] <- "Population"
NeEstim <- merge(NeEstim, loc, by = "Population")
### Change population names on figure
levels(NeEstim$Population)[c(4, 6, 14)] <- c("Northern Lagoon", "Misool", "Southern Lagoon")
NeEstim <- NeEstim[order(NeEstim$x),]
NeEstim$Population <- factor(NeEstim$Population,
                                 levels = NeEstim$Population[order(NeEstim$x)])
NeEstim[which(is.na(NeEstim$ci_high_recalib)),'ci_high_recalib'] <- Inf

## Set up y-axis breaks and labels
yaxis <- c(20, 50, 100, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4)

## Plot NeEstimator results
library(ggplot2)
Ne_plot <- ggplot(data = NeEstim, aes(x = Population, y = ne_recalib)) +
  geom_point(show.legend = F, colour = "transparent") +
  geom_ribbon(aes(x = 1:15, ymin = ci_low_recalib,
                  ymax = ci_high_recalib, fill = "95% CI"),
              fill = grey(0.9), show.legend = F) +
  geom_point(show.legend = F) +
  geom_hline(yintercept = c(1e2, 1e3, 1e4), linetype = "48", alpha = 0.5) +
  scale_y_continuous(trans = "log2",
                     breaks = yaxis,
                     labels = yaxis,
                     name = expression(italic(N[e]))) +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 60,
                                   hjust = 1,
                                   vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.ticks.x = element_blank()) 

ggsave(paste0(basedir, "results/NeEstimator/Ne_plot.pdf"),
       plot = Ne_plot, device = "pdf",
       width = 4.33,
       height = 4.33/1.5,
       units = "in")
