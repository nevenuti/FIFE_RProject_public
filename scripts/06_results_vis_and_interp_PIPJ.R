

# Title: Creating figures and tables of PIPJ (yellow pine) model results
# for publication

# Author: Nina Venuti

# Script inputs: 1) two fitted generalized additive models (GAMs) that 
# represent the "base models" for yr1 and yr2 PIPJ cone counts (ie.
# PIPJ_myr1 and PIPJ_myr2), and 2) the species and sample year-specific
# datasets used to fit those models (ie. PIPJ_data_yr1_full and PIPJ_
# data_yr2), all of which should be stored in the Environment pane after
# running the script "03_model_fitting_PIPJ.R". Additionally, 3) an ordered
# list of results visualization and communication functions (Functions 
# G through L), stored in the Environment pane after running the script
# "04_results_functions.R".

# Script outputs: 1) two sets of final model results figures - the first,
# displaying predicted cone counts for small, medium, and big PIPJs
# across the full range of crown volume damage, for each sample year
# separately; the second, displaying predicted cone counts for medium 
# PIPJs (only) in sample years 1 and 2 on the same axes; all figures are
# exported as both .png and .pdf files. 2) two percent change tables,
# displaying predicted cone counts [±1SE] and the estimated percent
# change in those cone counts at key CrownDamage values for small,
# medium, and big PIPJs in each sample year; percent change tables
# are exported as .csv files. And, 3) two contrasts summary tables, 
# displaying the results of pairwise comparisons of predicted cone 
# counts (ie. estimated mean differences in predicted counts [95% CI], 
# and associated p-values) at key CrownDamage values for small, medium, 
# and big PIPJs in each sample year; contrasts summary tables are 
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
# median field site for yr1 PIPJ...

# double checking PIPJ_myr1's model summary table to verify the median
# fecundity site for yr1...
summary(PIPJ_myr1) # CREEK is the median fecundity site

# running Function F to generate a predictions data.frame...
PIPJ_myr1_pred_df <- F_generate_predictions_df(model = PIPJ_myr1,
                                               sppyr_data = PIPJ_data_yr1_full,
                                               median_site = "CREEK")

# running Function G2 to identify the threshold value (ie. the CrownDamage
# value at which predicted cone counts decline most quickly) along each
# distinct DBH curve...
PIPJ_myr1_thresholds <- PIPJ_myr1_pred_df %>% 
  group_by(DBH_cm_scaled) %>% 
  summarize(threshold_val = G2_identify_thresholds_minderiv2(
    x_values = CrownDamage_percentvolume_scaled,
    curve_values = fit,
    smooth = NULL))

# running Function H to predict total conservative cone counts at these
# identified threshold values...
PIPJ_myr1_pred_df_thresholds <- H_predict_cones_at_thresholds(
  model = PIPJ_myr1,
  median_site = "CREEK",
  threshold_table = PIPJ_myr1_thresholds)


# plotting PIPJ yr1 base model predictions; this plot includes points  
# that represent the predicted threshold values for each DBH curve, and 
# a rug plot that displays field-collected CrownDamage values (ie. the 
# range of raw data on CrownDamage) for yr1 PIPJ...
PIPJ_myr1_plot <- PIPJ_myr1_pred_df %>% 
  ggplot(aes(x = CrownDamage_percentvolume_unscaled,
             y = fit_backtransformed,
             color = DBH_cm_unscaled)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "YlGnBu",
                     labels = c("10th percentile (47.5)",
                                "50th percentile (68.9)",
                                "90th percentile (102.0)")) +
  geom_ribbon(aes(ymin = CI_lowerlimit,
                  ymax = CI_upperlimit,
                  fill = DBH_cm_unscaled),
              alpha = 0.2,
              linewidth = 0) +
  scale_fill_brewer(palette = "YlGnBu",
                    labels = c("10th percentile (47.5)",
                               "50th percentile (68.9)",
                               "90th percentile (102.0)")) +
  geom_point(data = PIPJ_myr1_pred_df_thresholds,
             mapping = aes(x = CrownDamage_percentvolume_unscaled,
                           y = fit_backtransformed,
                           color = DBH_cm_unscaled),
             size = 4) +
  geom_rug(data = PIPJ_data_yr1_full,
           mapping = aes(x = CrownDamage_percentvolume,
                         y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "grey") +
  coord_trans(y = "log") +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 30, 40, 50)) +
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

PIPJ_myr1_plot


# running Function I to create a lookup table with key DBH values
# (affiliated with the 10th, 50th, and 90th percentile-sized trees),
# CrownDamage values (5%, 50%, and 95% crown volume damage, as well as 
# the threshold value, for each curve), and the median fecundity site 
# for yr1 PIPJ...
PIPJ_myr1_lookup_table <- I_create_lookup_table(
  sppyr_data = PIPJ_data_yr1_full,
  median_site = "CREEK",
  threshold_table = PIPJ_myr1_thresholds)

# running Function J to pivot this lookup table wider and add a column
# that holds the attributes of PIPJ_myr1...
PIPJ_myr1_lookup_table_wide <- J_pivot_lookup_table_wider(
  model = PIPJ_myr1,
  lookup_table = PIPJ_myr1_lookup_table)

# running Function K to create the PIPJ_myr1 percent change table...
PIPJ_myr1_percentchange_table <- K_create_percentchange_table(
  model = PIPJ_myr1,
  lookup_table = PIPJ_myr1_lookup_table,
  threshold_table = PIPJ_myr1_thresholds)

# running Functions L2 and L3 to create and then format the PIPJ_myr1
# contrasts summary table...
PIPJ_myr1_contrasts_summary_table <- L2_create_contrasts_summary_table(
  lookup_table_wide = PIPJ_myr1_lookup_table_wide)

PIPJ_myr1_contrasts_summary_table_wide <- L3_format_contrasts_summary_table(
  contrasts_summary_table = PIPJ_myr1_contrasts_summary_table,
  threshold_table = PIPJ_myr1_thresholds)



## moving on to the yr2 base model results...

# running Function F to predict total conservative cone counts for
# three distinct tree sizes (the 10th, 50th, and 90th percentile-
# sized trees), along the full spectrum of CrownDamage, and at the
# median field site for yr2 PIPJ...

# double checking PIPJ_myr2's model summary table to verify the median
# fecundity site for yr2...
summary(PIPJ_myr2) # NORTH is the median fecundity site

# running Function F to generate a predictions data.frame...
PIPJ_myr2_pred_df <- F_generate_predictions_df(model = PIPJ_myr2,
                                               sppyr_data = PIPJ_data_yr2,
                                               median_site = "NORTH")

# running Function G2 to identify the threshold value (ie. the CrownDamage
# value at which predicted cone counts decline most quickly) along each
# distinct DBH curve...
PIPJ_myr2_thresholds <- PIPJ_myr2_pred_df %>% 
  group_by(DBH_cm_scaled) %>% 
  summarize(threshold_val = G2_identify_thresholds_minderiv2(
    x_values = CrownDamage_percentvolume_scaled,
    curve_values = fit,
    smooth = NULL))

# running Function H to predict total conservative cone counts at these
# identified threshold values...
PIPJ_myr2_pred_df_thresholds <- H_predict_cones_at_thresholds(
  model = PIPJ_myr2,
  median_site = "NORTH",
  threshold_table = PIPJ_myr2_thresholds)


# plotting PIPJ yr2 base model predictions; this plot includes points  
# that represent the predicted threshold values for each DBH curve, and 
# a rug plot that displays field-collected CrownDamage values (ie. the 
# range of raw data on CrownDamage) for yr2 PIPJ...
PIPJ_myr2_plot <- PIPJ_myr2_pred_df %>% 
  ggplot(aes(x = CrownDamage_percentvolume_unscaled,
             y = fit_backtransformed,
             color = DBH_cm_unscaled)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "YlGnBu",
                     labels = c("10th percentile (47.2)",
                                "50th percentile (66.7)",
                                "90th percentile (96.4)")) +
  geom_ribbon(aes(ymin = CI_lowerlimit,
                  ymax = CI_upperlimit,
                  fill = DBH_cm_unscaled),
              alpha = 0.2,
              linewidth = 0) +
  scale_fill_brewer(palette = "YlGnBu",
                    labels = c("10th percentile (47.2)",
                               "50th percentile (66.7)",
                               "90th percentile (96.4)")) +
  geom_point(data = PIPJ_myr2_pred_df_thresholds,
             mapping = aes(x = CrownDamage_percentvolume_unscaled,
                           y = fit_backtransformed,
                           color = DBH_cm_unscaled),
             size = 4) +
  geom_rug(data = PIPJ_data_yr2,
           mapping = aes(x = CrownDamage_percentvolume,
                         y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "grey") +
  coord_trans(y = "log") +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 30, 50, 70)) +
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

PIPJ_myr2_plot


# running Function I to create a lookup table with key DBH values
# (affiliated with the 10th, 50th, and 90th percentile-sized trees),
# CrownDamage values (5%, 50%, and 95% crown volume damage, as well as 
# the threshold value, for each curve), and the median fecundity site 
# for yr2 PIPJ...
PIPJ_myr2_lookup_table <- I_create_lookup_table(
  sppyr_data = PIPJ_data_yr2,
  median_site = "NORTH",
  threshold_table = PIPJ_myr2_thresholds)

# running Function J to pivot this lookup table wider and add a column
# that holds the attributes of PIPJ_myr2...
PIPJ_myr2_lookup_table_wide <- J_pivot_lookup_table_wider(
  model = PIPJ_myr2,
  lookup_table = PIPJ_myr2_lookup_table)

# running Function K to create the PIPJ_myr2 percent change table...
PIPJ_myr2_percentchange_table <- K_create_percentchange_table(
  model = PIPJ_myr2,
  lookup_table = PIPJ_myr2_lookup_table,
  threshold_table = PIPJ_myr2_thresholds)

# running Functions L2 and L3 to create and then format the PIPJ_myr2
# contrasts summary table...
PIPJ_myr2_contrasts_summary_table <- L2_create_contrasts_summary_table(
  lookup_table_wide = PIPJ_myr2_lookup_table_wide)

PIPJ_myr2_contrasts_summary_table_wide <- L3_format_contrasts_summary_table(
  contrasts_summary_table = PIPJ_myr2_contrasts_summary_table,
  threshold_table = PIPJ_myr2_thresholds)



## finally, plotting predicted cone counts for the median-sized
## PIPJ in sample years 1 and 2 on the same set of axes...

# filtering the yr1 and yr2 predictions data.frames to include only
# those rows pertaining to the median-sized tree for each year,
# and adding a column for SampleEffort_yr to label rows accordingly...
PIPJ_myr1_pred_df_mediantree <- PIPJ_myr1_pred_df %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr1"))

PIPJ_myr2_pred_df_mediantree <- PIPJ_myr2_pred_df %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr2"))

# joining these median tree predictions data.frames together...
PIPJ_mediantree_pred_df_BOTHYRS <- bind_rows(PIPJ_myr1_pred_df_mediantree,
                                             PIPJ_myr2_pred_df_mediantree)

# repeating these same steps for the yr1 and yr2 threshold predictions
# data.frames...
PIPJ_myr1_pred_thresholds_mediantree <- PIPJ_myr1_pred_df_thresholds %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr1"))

PIPJ_myr2_pred_thresholds_mediantree <- PIPJ_myr2_pred_df_thresholds %>% 
  filter(DBH_cm_scaled == median(DBH_cm_scaled)) %>% 
  mutate(SampleEffort_yr = as.factor("yr2"))

PIPJ_mediantree_pred_thresholds_BOTHYRS <- bind_rows(
  PIPJ_myr1_pred_thresholds_mediantree,
  PIPJ_myr2_pred_thresholds_mediantree)

# plotting the yr1 and yr2 median tree predictions on the same
# axes; note that this plot includes points that represent the predicted
# threshold values for each curve, but does not include a rug plot,
# due to the difficulty of legibly displaying the raw values of CrownDamage 
# for trees surveyed in both yr1 and yr2, all at once...
PIPJ_mediantree_plot <- PIPJ_mediantree_pred_df_BOTHYRS %>% 
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
  geom_point(data = PIPJ_mediantree_pred_thresholds_BOTHYRS,
             mapping = aes(x = CrownDamage_percentvolume_unscaled,
                           y = fit_backtransformed,
                           color = SampleEffort_yr),
             size = 4) +
  coord_trans(y = "log") +
  scale_y_continuous(limits = c(0.6, 40),
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

PIPJ_mediantree_plot

# adding a rug plot to the PIPJ_mediantree_plot, to test legibility...
PIPJ_mediantree_plot_wrug <- PIPJ_mediantree_plot + 
  geom_rug(data = PIPJ_data_yr1_full,
           aes(x = CrownDamage_percentvolume,
               y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "paleturquoise4",
           linewidth = 0.4) +
  geom_rug(data = PIPJ_data_yr2,
           aes(x = CrownDamage_percentvolume,
               y = 5),
           inherit.aes = FALSE,
           sides = "b",
           position = "jitter",
           color = "cadetblue3",
           linewidth = 0.4)

PIPJ_mediantree_plot_wrug

# adding a stacked, color-coded rug plot to the PIPJ_mediantree_plot, to 
# try to increase legibility of rug...
PIPJ_mediantree_plot_wrug_stacked <- PIPJ_mediantree_plot + 
  geom_point(data = PIPJ_data_yr1_full,
             aes(x = CrownDamage_percentvolume +
                   runif(nrow(PIPJ_data_yr1_full), min = -0.75, max = 0.75),
                 y = 0.67),
             shape = "|", 
             size = 3.5, 
             color = "paleturquoise4") +
  geom_point(data = PIPJ_data_yr2,
             aes(x = CrownDamage_percentvolume +
                   runif(nrow(PIPJ_data_yr2), min = -0.75, max = 0.75),
                 y = 0.6),
             shape = "|", 
             size = 3.5, 
             color = "cadetblue3")

PIPJ_mediantree_plot_wrug_stacked

# adding a grey-scale rug plot to the PIPJ_mediantree_plot, where 
# grey tick marks represent trees that died (or fell over, etc.)
# between years 1 and 2 (and were therefore included in the year 1 model,
# but excluded from the year 2 model), and black tick marks represent 
# trees that remained alive across both sample years (and were therefore
# included in both models)...
PIPJ_mediantree_plot_wrug_greyscale <- PIPJ_mediantree_plot + 
  geom_rug(data = subset(PIPJ_data_yr1_full, TreeStatus_yr2 != "A"),
           aes(x = CrownDamage_percentvolume,
               y = 5),
           color = "grey",
           sides = "b",
           position = position_jitter(width = 1),
           linewidth = 0.4) +
  geom_rug(data = subset(PIPJ_data_yr1_full, TreeStatus_yr2 == "A"),
           aes(x = CrownDamage_percentvolume,
               y = 5),
           color = "black", 
           sides = "b",
           position = position_jitter(width = 1),
           linewidth = 0.4)

PIPJ_mediantree_plot_wrug_greyscale



## exporting the final results figures, percent change tables, and
## contrasts summary tables for PIPJ yr1 and yr2...

# exporting the PIPJ_myr1 plot, the PIPJ_myr2 plot, and the cross-
# sample year median tree plot (with and without grey-scale and
# stacked rug plots) as both .pngs and .pdfs...
ggsave(plot = PIPJ_myr1_plot,
       filename = "results/final/figures/PIPJ_yr1_plot.png")
ggsave(plot = PIPJ_myr1_plot,
       filename = "results/final/figures/PIPJ_yr1_plot.pdf")

ggsave(plot = PIPJ_myr2_plot,
       filename = "results/final/figures/PIPJ_yr2_plot.png")
ggsave(plot = PIPJ_myr2_plot,
       filename = "results/final/figures/PIPJ_yr2_plot.pdf")

ggsave(plot = PIPJ_mediantree_plot,
       filename = "results/final/figures/PIPJ_mediantree_plot.png")
ggsave(plot = PIPJ_mediantree_plot,
       filename = "results/final/figures/PIPJ_mediantree_plot.pdf")

ggsave(plot = PIPJ_mediantree_plot_wrug_greyscale,
       filename = "results/final/figures/PIPJ_mediantree_plot_wrug.png")
ggsave(plot = PIPJ_mediantree_plot_wrug_greyscale,
       filename = "results/final/figures/PIPJ_mediantree_plot_wrug.pdf")

ggsave(plot = PIPJ_mediantree_plot_wrug_stacked,
       filename = "results/final/figures/PIPJ_mediantree_plot_wrug_stacked.png")
ggsave(plot = PIPJ_mediantree_plot_wrug_stacked,
       filename = "results/final/figures/PIPJ_mediantree_plot_wrug_stacked.pdf")

# exporting the percent change tables and the contrasts summary tables
# for both sample years as .csv files...
write.csv(PIPJ_myr1_percentchange_table,
          "results/final/tables/PIPJ_yr1_percentchange_table.csv",
          row.names = FALSE)
write.csv(PIPJ_myr1_contrasts_summary_table_wide,
          "results/final/tables/PIPJ_yr1_contrastssummary_table.csv",
          row.names = FALSE)

write.csv(PIPJ_myr2_percentchange_table,
          "results/final/tables/PIPJ_yr2_percentchange_table.csv",
          row.names = FALSE)
write.csv(PIPJ_myr2_contrasts_summary_table_wide,
          "results/final/tables/PIPJ_yr2_contrastssummary_table.csv",
          row.names = FALSE)





