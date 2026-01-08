

# Title: Creating additional model results figures for inclusion in
# the Supplemental Materials section of the publication

# Author: Nina Venuti

# Script inputs: a series of fitted generalized additive models (GAMs)
# of ABCO, PIPJ, and PILA cone counts one- and two-years after fire
# injury, stored in the Environment pane after running each, species-
# specific iteration of script #3 ("03_model_fitting_SPPCODE.R")

# Script outputs: two supplemental results figures - the first, 
# displaying the estimated marginal effects of tree size on predicted cone
# counts for each species-year base model (e.g., "ABCO_myr1", "PIPJ_myr2"),
# and the second, displaying the estimated marginal effects of mistletoe 
# infestation on predicted cone counts in PIPJ two years post-fire (the only
# species-year in which mistletoe was found to have a significant
# effect on total cone counts). All plots are exported as both .png and 
# .pdf files.

library("tidyverse")
library("gratia")
library("cowplot")


## compiling the six DBH partial effects plots (one from each species-year
## base model) into a single figure for publication...

# using the gratia package to recreate the partial effects plots
# for the smooth term "s(DBH_cm_scaled)" for each mgcv-fit, species-year 
# base model...

# starting with ABCO_myr1...
DBH_effects_plot_ABCOyr1 <- gratia::draw(ABCO_myr1, select = 1) +
  labs(title = NULL,
       caption = NULL,
       x = "DBH (scaled)") +
  scale_y_continuous(limits = c(-1, 2.2),
                     breaks = c(-1, 0, 1, 2)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
DBH_effects_plot_ABCOyr1

# moving on to ABCO_myr2...
DBH_effects_plot_ABCOyr2 <- gratia::draw(ABCO_myr2, select = 1) +
  labs(title = NULL,
       caption = NULL,
       x = "DBH (scaled)") +
  scale_y_continuous(limits = c(-1, 1.2),
                     breaks = c(-1, 0, 1)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
DBH_effects_plot_ABCOyr2

# then to PIPJ_myr1...
DBH_effects_plot_PIPJyr1 <- gratia::draw(PIPJ_myr1, select = 1) +
  labs(title = NULL,
       caption = NULL,
       x = "DBH (scaled)") +
  scale_y_continuous(limits = c(-2, 4),
                     breaks = c(-2, 0, 2, 4)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
DBH_effects_plot_PIPJyr1

# and PIPJ_myr2...
DBH_effects_plot_PIPJyr2 <- gratia::draw(PIPJ_myr2, select = 1) +
  labs(title = NULL,
       caption = NULL,
       x = "DBH (scaled)") +
  scale_y_continuous(limits = c(-2, 4.2),
                     breaks = c(-2, 0, 2, 4)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
DBH_effects_plot_PIPJyr2

# then to PILA_myr1...
DBH_effects_plot_PILAyr1 <- gratia::draw(PILA_myr1, select = 1) +
  labs(title = NULL,
       caption = NULL,
       x = "DBH (scaled)") +
  scale_y_continuous(limits = c(-3, 4),
                     breaks = c(-2, 0, 2, 4)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
DBH_effects_plot_PILAyr1

# and finally, PILA_myr2...
DBH_effects_plot_PILAyr2 <- gratia::draw(PILA_myr2, select = 1) +
  labs(title = NULL,
       caption = NULL,
       x = "DBH (scaled)") +
  scale_y_continuous(limits = c(-2, 4.2),
                     breaks = c(-2, 0, 2, 4)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
DBH_effects_plot_PILAyr2

# combining these six partial effects plots into a single panel for
# publication...
DBH_effects_plot_ALLSPPYRS <- plot_grid(DBH_effects_plot_ABCOyr1,
                                        DBH_effects_plot_ABCOyr2,
                                        DBH_effects_plot_PIPJyr1,
                                        DBH_effects_plot_PIPJyr2,
                                        DBH_effects_plot_PILAyr1,
                                        DBH_effects_plot_PILAyr2,
                                        nrow = 3,
                                        ncol = 2,
                                        labels = c("white fir year 1",
                                                   "white fir year 2",
                                                   "yellow pine year 1",
                                                   "yellow pine year 2",
                                                   "sugar pine year 1", 
                                                   "sugar pine year 2"),
                                        label_size = 14,
                                        hjust = c(-0.5, -0.5,
                                                  -0.41, -0.41,
                                                  -0.42, -0.42),
                                        vjust = 2.3,
                                        align = "hv")
DBH_effects_plot_ALLSPPYRS

# exporting this panel plot as both a .png and a .pdf file...
ggsave(plot = DBH_effects_plot_ALLSPPYRS,
       filename = "results/final/figures/DBH_partial_effects_plot_allspp.png",
       width = 18, 
       height = 22,
       units = "cm",
       dpi = 400)

ggsave(plot = DBH_effects_plot_ALLSPPYRS,
       filename = "results/final/figures/DBH_partial_effects_plot_allspp.pdf",
       width = 18, 
       height = 22,
       units = "cm",
       dpi = 400)



## extracting and saving the partial effects plot for the smooth
## term "s(MistletoeClumpsOrBranches_count_yr2_scaled)" from the relevant
## model of PIPJ fecundity two years post-fire (ie. PIPJ_m8)...

# recreating the Mistletoe_count partial effects plot using the gratia
# package...
mistletoe_effects_plot_PIPJyr2 <- gratia::draw(PIPJ_m8, select = 4, 
                                               rug = FALSE) +
  labs(title = NULL,
       caption = NULL,
       x = "Mistletoe infestation intensity (scaled)") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
mistletoe_effects_plot_PIPJyr2

# adding a rug plot by hand, to allow the tick marks to be jittered...
mistletoe_effects_plot_PIPJyr2_plusrug <- mistletoe_effects_plot_PIPJyr2 +
  geom_rug(data = PIPJ_data_yr2,
           aes(x = MistletoeClumpsOrBranches_count_yr2_scaled,
               y = 0),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           alpha = 0.5)
mistletoe_effects_plot_PIPJyr2_plusrug

# exporting the mistletoe partial effects plot (with rug) as both a 
# .png and a .pdf file...
ggsave(plot = mistletoe_effects_plot_PIPJyr2_plusrug,
       filename = "results/final/figures/Mistletoe_partial_effects_plot_PIPJyr2.png",
       width = 6, 
       height = 4,
       units = "in",
       dpi = 400)

ggsave(plot = mistletoe_effects_plot_PIPJyr2_plusrug,
       filename = "results/final/figures/Mistletoe_partial_effects_plot_PIPJyr2.pdf",
       width = 6, 
       height = 4,
       units = "in",
       dpi = 400)




