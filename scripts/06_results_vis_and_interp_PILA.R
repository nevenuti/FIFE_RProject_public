

# Title: Creating figures and tables of PILA (sugar pine) model results
# for publication

# Author: Nina Venuti

# Script inputs: 1) two fitted generalized additive models (GAMs) that 
# represent the "base models" for yr1 and yr2 PILA cone counts (ie.
# PILA_myr1 and PILA_myr2), and 2) the species and sample year-specific
# datasets used to fit those models (ie. PILA_data_yr1 and PILA_data_
# yr2), all of which should be stored in the Environment pane after
# running the script "03_model_fitting_PILA.R". Additionally, 3) an ordered
# list of results visualization and communication functions (Functions 
# G through L), stored in the Environment pane after running the script
# "04_results_functions.R".

# Script outputs: 1) two sets of final model results figures - the first,
# displaying predicted cone counts for small, medium, and big PILAs
# across the full range of crown volume damage, for each sample year
# separately; the second, displaying predicted cone counts for medium 
# PILAs (only) in sample years 1 and 2 on the same axes; all figures are
# exported as both .png and .pdf files. 2) two percent change tables,
# displaying predicted cone counts [±1SE] and the estimated percent
# change in those cone counts at key CrownDamage values for small,
# medium, and big PILAs in each sample year; percent change tables
# are exported as .csv files. And, 3) two contrasts summary tables, 
# displaying the results of pairwise comparisons of predicted cone 
# counts (ie. estimated mean differences in predicted counts [95% CI], 
# and associated p-values) at key CrownDamage values for small, medium, 
# and big PILAs in each sample year; contrasts summary tables are 
# exported as .csv files.

library("tidyverse")


## following the workflow outlined in script #4, we will use Functions F
## through H to generate the information needed to plot the results
## of each base model, plot those results, and then use Functions I
## through L to create the percent change tables, and contrasts summary
## tables, affiliated with each results figure...

## starting with the yr1 base model results...

# first, extracting the means and standard deviations for the predictors
# DBH_cm and CrownDamage_percentvolume (using the full, cross-species
# dataset) to enable backtransforming the scaled versions of these 
# variables during plotting and table creation below...
DBH_mean <- mean(model_data$DBH_cm)
DBH_sd <- sd(model_data$DBH_cm)
CrownDamage_mean <- mean(model_data$CrownDamage_percentvolume)
CrownDamage_sd <- sd(model_data$CrownDamage_percentvolume)

# running Function F to predict total conservative cone counts for
# three distinct tree sizes (the 10th, 50th, and 90th percentile-
# sized trees), along the full spectrum of CrownDamage, and at the
# median field site for yr1 PILA...

# double checking PILA_myr1's model summary table to verify the median
# fecundity site for yr1...
summary(PILA_myr1) # NORTH is the median fecundity site

# running Function F to generate a predictions data.frame...
PILA_myr1_pred_df <- F_generate_predictions_df(model = PILA_myr1,
                                               sppyr_data = PILA_data_yr1,
                                               median_site = "NORTH")

# running Function G2 to identify the threshold value (ie. the CrownDamage
# value at which predicted cone counts decline most quickly) along each
# distinct DBH curve...
PILA_myr1_thresholds <- PILA_myr1_pred_df %>% 
  group_by(DBH_cm_scaled) %>% 
  summarize(threshold_val = G2_identify_thresholds_minderiv2(
    x_values = CrownDamage_percentvolume_scaled,
    curve_values = fit,
    smooth = NULL))

# running Function H to predict total conservative cone counts at these
# identified threshold values...
PILA_myr1_pred_df_thresholds <- H_predict_cones_at_thresholds(
  model = PILA_myr1,
  median_site = "NORTH",
  threshold_table = PILA_myr1_thresholds)


# plotting PILA yr1 base model predictions; this plot includes points  
# that represent the predicted threshold values for each DBH curve, and 
# a rug plot that displays field-collected CrownDamage values (ie. the 
# range of raw data on CrownDamage) for yr1 PILA...
PILA_myr1_plot <- PILA_myr1_pred_df %>% 
  ggplot(aes(x = CrownDamage_percentvolume_unscaled,
             y = fit_backtransformed,
             color = DBH_cm_unscaled)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "YlGnBu",
                     labels = c("10th percentile (47.1)",
                                "50th percentile (68.1)",
                                "90th percentile (106.7)")) +
  geom_ribbon(aes(ymin = CI_lowerlimit,
                  ymax = CI_upperlimit,
                  fill = DBH_cm_unscaled),
              alpha = 0.2,
              linewidth = 0) +
  scale_fill_brewer(palette = "YlGnBu",
                    labels = c("10th percentile (47.1)",
                               "50th percentile (68.1)",
                               "90th percentile (106.7)")) +
  geom_point(data = PILA_myr1_pred_df_thresholds,
             mapping = aes(x = CrownDamage_percentvolume_unscaled,
                           y = fit_backtransformed,
                           color = DBH_cm_unscaled),
             size = 4) +
  geom_rug(data = PILA_data_yr1,
           mapping = aes(x = CrownDamage_percentvolume,
                         y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "grey") +
  coord_trans(y = "log") +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 40, 60, 80, 100)) +
  xlab("Crown volume damage (%)") +
  ylab("Number of cones per tree") +
  guides(color = guide_legend(title = "Tree size (cm DBH)", 
                              reverse = TRUE),
         fill = guide_legend(title = "Tree size (cm DBH)",
                             reverse = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))

PILA_myr1_plot


# running Function I to create a lookup table with key DBH values
# (affiliated with the 10th, 50th, and 90th percentile-sized trees),
# CrownDamage values (5%, 50%, and 95% crown volume damage, as well as 
# the threshold value, for each curve), and the median fecundity site 
# for yr1 PILA...
PILA_myr1_lookup_table <- I_create_lookup_table(
  sppyr_data = PILA_data_yr1,
  median_site = "NORTH",
  threshold_table = PILA_myr1_thresholds)

# running Function J to pivot this lookup table wider and add a column
# that holds the attributes of PILA_myr1...
PILA_myr1_lookup_table_wide <- J_pivot_lookup_table_wider(
  model = PILA_myr1,
  lookup_table = PILA_myr1_lookup_table)

# running Function K to create the PILA_myr1 percent change table...
PILA_myr1_percentchange_table <- K_create_percentchange_table(
  model = PILA_myr1,
  lookup_table = PILA_myr1_lookup_table,
  threshold_table = PILA_myr1_thresholds)

# rearranging the columns of the PILA_myr1 percent change table,
# to reflect the fact that two out of three of the DBH curves have 
# threshold values that are <50%, and to facilitate table interpretation...
PILA_myr1_percentchange_table <- PILA_myr1_percentchange_table %>% 
  select(c(1:3, 5, 4, 6, 8, 7, 10, 9, 11))

# running Functions L2 and L3 to create and then format the PILA_myr1
# contrasts summary table...
PILA_myr1_contrasts_summary_table <- L2_create_contrasts_summary_table(
  lookup_table_wide = PILA_myr1_lookup_table_wide)

PILA_myr1_contrasts_summary_table_wide <- L3_format_contrasts_summary_table(
  contrasts_summary_table = PILA_myr1_contrasts_summary_table,
  threshold_table = PILA_myr1_thresholds)

# again, rearranging the columns of the PILA_myr1 contrasts table, 
# to reflect the fact that two of the three threshold values are <50%, 
# and to facilitate table interpretation...
PILA_myr1_contrasts_summary_table_wide <-
  PILA_myr1_contrasts_summary_table_wide %>% 
  select(c(1:2, 5:6, 3:4, 9:10, 7:8, 11:12))



## moving on to the yr2 base model results...

# running Function F to predict total conservative cone counts for
# three distinct tree sizes (the 10th, 50th, and 90th percentile-
# sized trees), along the full spectrum of CrownDamage, and at the
# median field site for yr2 PILA...

# double checking PILA_myr2's model summary table to verify the median
# fecundity site for yr2...
summary(PILA_myr2) # CREEK is the median fecundity site

# running Function F to generate a predictions data.frame...
PILA_myr2_pred_df <- F_generate_predictions_df(model = PILA_myr2,
                                               sppyr_data = PILA_data_yr2,
                                               median_site = "CREEK")

# running Function G2 to identify the threshold value (ie. the CrownDamage
# value at which predicted cone counts decline most quickly) along each
# distinct DBH curve...
PILA_myr2_thresholds <- PILA_myr2_pred_df %>% 
  group_by(DBH_cm_scaled) %>% 
  summarize(threshold_val = G2_identify_thresholds_minderiv2(
    x_values = CrownDamage_percentvolume_scaled,
    curve_values = fit,
    smooth = NULL))

# running Function H to predict total conservative cone counts at these
# identified threshold values...
PILA_myr2_pred_df_thresholds <- H_predict_cones_at_thresholds(
  model = PILA_myr2,
  median_site = "CREEK",
  threshold_table = PILA_myr2_thresholds)


# plotting PILA yr2 base model predictions; this plot includes points  
# that represent the predicted threshold values for each DBH curve, and 
# a rug plot that displays field-collected CrownDamage values (ie. the 
# range of raw data on CrownDamage) for yr2 PILA...
PILA_myr2_plot <- PILA_myr2_pred_df %>% 
  ggplot(aes(x = CrownDamage_percentvolume_unscaled,
             y = fit_backtransformed,
             color = DBH_cm_unscaled)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "YlGnBu",
                     labels = c("10th percentile (48.6)",
                                "50th percentile (65.4)",
                                "90th percentile (106.7)")) +
  geom_ribbon(aes(ymin = CI_lowerlimit,
                  ymax = CI_upperlimit,
                  fill = DBH_cm_unscaled),
              alpha = 0.2,
              linewidth = 0) +
  scale_fill_brewer(palette = "YlGnBu",
                    labels = c("10th percentile (48.6)",
                               "50th percentile (65.4)",
                               "90th percentile (106.7)")) +
  geom_point(data = PILA_myr2_pred_df_thresholds,
             mapping = aes(x = CrownDamage_percentvolume_unscaled,
                           y = fit_backtransformed,
                           color = DBH_cm_unscaled),
             size = 4) +
  geom_rug(data = PILA_data_yr2,
           mapping = aes(x = CrownDamage_percentvolume,
                         y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "grey") +
  coord_trans(y = "log") +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 40, 60)) +
  xlim(0, 100) +
  xlab("Crown volume damage (%)") +
  ylab("Number of cones per tree") +
  guides(color = guide_legend(title = "Tree size (cm DBH)", 
                              reverse = TRUE),
         fill = guide_legend(title = "Tree size (cm DBH)",
                             reverse = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))

PILA_myr2_plot


# running Function I to create a lookup table with key DBH values
# (affiliated with the 10th, 50th, and 90th percentile-sized trees),
# CrownDamage values (5%, 50%, and 95% crown volume damage, as well as 
# the threshold value, for each curve), and the median fecundity site 
# for yr2 PILA...
PILA_myr2_lookup_table <- I_create_lookup_table(
  sppyr_data = PILA_data_yr2,
  median_site = "CREEK",
  threshold_table = PILA_myr2_thresholds)

# running Function J to pivot this lookup table wider and add a column
# that holds the attributes of PILA_myr2...
PILA_myr2_lookup_table_wide <- J_pivot_lookup_table_wider(
  model = PILA_myr2,
  lookup_table = PILA_myr2_lookup_table)

# running Function K to create the PILA_myr2 percent change table...
PILA_myr2_percentchange_table <- K_create_percentchange_table(
  model = PILA_myr2,
  lookup_table = PILA_myr2_lookup_table,
  threshold_table = PILA_myr2_thresholds)

# rearranging the columns of the PILA_myr2 percent change table,
# to reflect the fact that the threshold values for all three
# DBH curves are <50%, and to facilitate table interpretation...
PILA_myr2_percentchange_table <- PILA_myr2_percentchange_table %>% 
  select(c(1:3, 5, 4, 6, 8, 7, 10, 9, 11))

# running Functions L2 and L3 to create and then format the PILA_myr2
# contrasts summary table...
PILA_myr2_contrasts_summary_table <- L2_create_contrasts_summary_table(
  lookup_table_wide = PILA_myr2_lookup_table_wide)

PILA_myr2_contrasts_summary_table_wide <- L3_format_contrasts_summary_table(
  contrasts_summary_table = PILA_myr2_contrasts_summary_table,
  threshold_table = PILA_myr2_thresholds)

# again, rearranging the columns of the PILA_myr2 contrasts table, 
# to reflect the fact that threshold values are <50%, and to facilitate
# table interpretation...
PILA_myr2_contrasts_summary_table_wide <-
  PILA_myr2_contrasts_summary_table_wide %>% 
  select(c(1:2, 5:6, 3:4, 9:10, 7:8, 11:12))



## finally, plotting predicted cone counts for the median-sized
## PILA in sample years 1 and 2 on the same set of axes...

# filtering the yr1 and yr2 predictions data.frames to include only
# those rows pertaining to the median-sized tree for each year,
# and adding a column for SampleEffort_yr to label rows accordingly...
PILA_myr1_pred_df_mediantree <- PILA_myr1_pred_df %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr1"))

PILA_myr2_pred_df_mediantree <- PILA_myr2_pred_df %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr2"))

# joining these median tree predictions data.frames together...
PILA_mediantree_pred_df_BOTHYRS <- bind_rows(PILA_myr1_pred_df_mediantree,
                                             PILA_myr2_pred_df_mediantree)

# repeating these same steps for the yr1 and yr2 threshold predictions
# data.frames...
PILA_myr1_pred_thresholds_mediantree <- PILA_myr1_pred_df_thresholds %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr1"))

PILA_myr2_pred_thresholds_mediantree <- PILA_myr2_pred_df_thresholds %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr2"))

PILA_mediantree_pred_thresholds_BOTHYRS <- bind_rows(
  PILA_myr1_pred_thresholds_mediantree,
  PILA_myr2_pred_thresholds_mediantree)

# plotting the yr1 and yr2 median tree predictions on the same
# axes; note that this plot includes points that represent the predicted
# threshold values for each curve, but does not include a rug plot,
# due to the difficulty of legibly displaying the raw values of CrownDamage 
# for trees surveyed in both yr1 and yr2, all at once...
PILA_mediantree_plot <- PILA_mediantree_pred_df_BOTHYRS %>% 
  ggplot(aes(x = CrownDamage_percentvolume_unscaled,
             y = fit_backtransformed,
             color = SampleEffort_yr)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("yr1" = "paleturquoise4", 
                                "yr2" = "cadetblue3"),
                     labels = c("year 1", 
                                "year 2")) +
  geom_ribbon(aes(ymin = CI_lowerlimit,
                  ymax = CI_upperlimit,
                  fill = SampleEffort_yr),
              alpha = 0.2,
              linewidth = 0) +
  scale_fill_manual(values = c("yr1" = "paleturquoise4", 
                               "yr2" = "cadetblue3"),
                    labels = c("year 1",
                               "year 2")) +
  geom_point(data = PILA_mediantree_pred_thresholds_BOTHYRS,
             mapping = aes(x = CrownDamage_percentvolume_unscaled,
                           y = fit_backtransformed,
                           color = SampleEffort_yr),
             size = 4) +
  coord_trans(y = "log") +
  scale_y_continuous(limits = c(1.5, 40),
                     breaks = c(0, 5, 10, 20, 30, 40)) +
  xlab("Crown volume damage (%)") +
  ylab("Number of cones per tree") +
  theme_classic() +
  guides(color = guide_legend(title = "Sample year"),
         fill = guide_legend(title = "Sample year")) + 
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))

PILA_mediantree_plot

# adding a rug plot to the PILA_mediantree_plot, to test legibility...
PILA_mediantree_plot_wrug <- PILA_mediantree_plot + 
  geom_rug(data = PILA_data_yr1,
           aes(x = CrownDamage_percentvolume,
               y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "paleturquoise4",
           linewidth = 0.4) +
  geom_rug(data = PILA_data_yr2,
           aes(x = CrownDamage_percentvolume,
               y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "cadetblue3",
           linewidth = 0.4)

PILA_mediantree_plot_wrug

# adding a stacked, color-coded rug plot to the PILA_mediantree_plot, to 
# try to increase legibility of rug...
PILA_mediantree_plot_wrug_stacked <- PILA_mediantree_plot + 
  geom_point(data = PILA_data_yr1,
             aes(x = CrownDamage_percentvolume +
                   runif(nrow(PILA_data_yr1), min = -0.75, max = 0.75),
                 y = 1.62),
             shape = "|", 
             size = 3.5, 
             color = "paleturquoise4") +
  geom_point(data = PILA_data_yr2,
             aes(x = CrownDamage_percentvolume +
                   runif(nrow(PILA_data_yr2), min = -0.75, max = 0.75),
                 y = 1.5),
             shape = "|", 
             size = 3.5, 
             color = "cadetblue3")

PILA_mediantree_plot_wrug_stacked

# adding a grey-scale rug plot to the PILA_mediantree_plot, where 
# grey tick marks represent trees that died (or fell over, etc.)
# between years 1 and 2 (and were therefore included in the year 1 model,
# but excluded from the year 2 model), and black tick marks represent 
# trees that remained alive across both sample years (and were therefore
# included in both models)...
PILA_mediantree_plot_wrug_greyscale <- PILA_mediantree_plot + 
  geom_rug(data = subset(PILA_data_yr1, TreeStatus_yr2 != "A"),
           aes(x = CrownDamage_percentvolume,
               y = 5),
           color = "grey",
           sides = "b",
           position = position_jitter(width = 1),
           linewidth = 0.4) +
  geom_rug(data = subset(PILA_data_yr1, TreeStatus_yr2 == "A"),
           aes(x = CrownDamage_percentvolume,
               y = 5),
           color = "black", 
           sides = "b",
           position = position_jitter(width = 1),
           linewidth = 0.4)

PILA_mediantree_plot_wrug_greyscale



## exporting the final results figures, percent change tables, and
## contrasts summary tables for PILA yr1 and yr2...

# exporting the PILA_myr1 plot, the PILA_myr2 plot, and the cross-
# sample year median tree plot (with and without grey-scale and
# stacked rug plots) as both .pngs and .pdfs...
ggsave(plot = PILA_myr1_plot,
       filename = "results/final/figures/PILA_yr1_plot.png")
ggsave(plot = PILA_myr1_plot,
       filename = "results/final/figures/PILA_yr1_plot.pdf")

ggsave(plot = PILA_myr2_plot,
       filename = "results/final/figures/PILA_yr2_plot.png")
ggsave(plot = PILA_myr2_plot,
       filename = "results/final/figures/PILA_yr2_plot.pdf")

ggsave(plot = PILA_mediantree_plot,
       filename = "results/final/figures/PILA_mediantree_plot.png")
ggsave(plot = PILA_mediantree_plot,
       filename = "results/final/figures/PILA_mediantree_plot.pdf")

ggsave(plot = PILA_mediantree_plot_wrug_greyscale,
       filename = "results/final/figures/PILA_mediantree_plot_wrug.png")
ggsave(plot = PILA_mediantree_plot_wrug_greyscale,
       filename = "results/final/figures/PILA_mediantree_plot_wrug.pdf")

ggsave(plot = PILA_mediantree_plot_wrug_stacked,
       filename = "results/final/figures/PILA_mediantree_plot_wrug_stacked.png")
ggsave(plot = PILA_mediantree_plot_wrug_stacked,
       filename = "results/final/figures/PILA_mediantree_plot_wrug_stacked.pdf")

# exporting the percent change tables and the contrasts summary tables
# for both sample years as .csv files...
write.csv(PILA_myr1_percentchange_table,
          "results/final/tables/PILA_yr1_percentchange_table.csv",
          row.names = FALSE)
write.csv(PILA_myr1_contrasts_summary_table_wide,
          "results/final/tables/PILA_yr1_contrastssummary_table.csv",
          row.names = FALSE)

write.csv(PILA_myr2_percentchange_table,
          "results/final/tables/PILA_yr2_percentchange_table.csv",
          row.names = FALSE)
write.csv(PILA_myr2_contrasts_summary_table_wide,
          "results/final/tables/PILA_yr2_contrastssummary_table.csv",
          row.names = FALSE)





