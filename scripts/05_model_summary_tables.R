

# Title: Compiling ABCO, PIPJ, and PILA model results into summary tables
# for publication

# Author: Nina Venuti

# Script inputs: 1) a series of fitted generalized additive models (GAMs)
# of ABCO, PIPJ, and PILA cone counts one and two-year(s) after fire
# injury, stored in the Environment pane after running each, species-specific
# iteration of script #3 ("03_model_fitting_SPPCODE.R"), and 2) an 
# ordered list of model summary functions (Functions A through E),
# stored in the Environment pane after running script #4 ("04_results_
# functions.R")

# Script outputs: five model summary tables, exported as .csv files
# (see descriptions of the outputs of Functions A through E in script
# #4 for further details on the information contained in each summary
# table)

library("tidyverse")


## creating vectors of model names, to be used as inputs for Functions A
## through E below...

# creating a vector of base model names...
base_model_names <- c("ABCO_myr1", "ABCO_myr2", "PIPJ_myr1", "PIPJ_myr2",
                      "PILA_myr1", "PILA_myr2")

# creating a vector of BoleChar model names (ie. names of models that
# include the term BoleCharMeanHeight_m as a predictor); note that 
# such models are labeled as either SPPCODE_m2 (for yr1 models) or 
# SPPCODE_m6 (for yr2 models) in each iteration of script #3...
bolechar_model_names <- c("ABCO_m2", "ABCO_m6", "PIPJ_m2", "PIPJ_m6",
                          "PILA_m2", "PILA_m6")

# creating a vector of CanopyPosition model names (ie. names of models
# that include the term CanopyPosition as a predictor); note that 
# such models are labeled as either SPPCODE_m3 (for yr1 models) or 
# SPPCODE_m7 (for yr2 models) in each iteration of script #3...
canopyposition_model_names <- c("ABCO_m3", "ABCO_m7", "PIPJ_m3", "PIPJ_m7",
                                "PILA_m3", "PILA_m7")

# creating a vector of Mistletoe model names (ie. names of models that
# include the term MistletoeClumpsOrBranches_count as a predictor); note
# that such models are labeled as either SPPCODE_m4 (for yr1 models) or
# SPPCODE_m8 (for yr2 models) in each iteration of script #3, and that 
# Mistletoe_count was not included as a predictor in any PILA models
# (see the script "03_model_fitting_PILA.R" for further discussion of
# this decision)...
mistletoe_model_names <- c("ABCO_m4", "ABCO_m8", "PIPJ_m4", "PIPJ_m8")


## running Functions A through E (using relevant model name vectors)
## to generate model summary tables...

# Function A: create_basemodel_summary_table...
basemodel_summary_table <- A_create_basemodel_summary_table(base_model_names)

# Function B: compile_basemodel_site_effects...
basemodel_siteeffects_table <- B_compile_basemodel_site_effects(
  base_model_names)

# Function C: compile_bolechar_model_estimates...
bolechar_summary_table <- C_compile_bolechar_model_estimates(
  bolechar_model_names) 

# Function D: compile_canopyposition_model_estimates...
canopyposition_summary_table <- D_compile_canopyposition_model_estimates(
  canopyposition_model_names)

# Function E: compile_mistletoe_model_estimates...
mistletoe_summary_table <- E_compile_mistletoe_model_estimates(
  mistletoe_model_names)


## exporting each of these summary tables as .csv files...

# exporting the base model summary table...
write.csv(basemodel_summary_table,
          "results/final/tables/basemodel_summary_table.csv",
          row.names = FALSE)

# exporting the base model site effects table...
write.csv(basemodel_siteeffects_table,
          "results/final/tables/basemodel_siteeffects_table.csv",
          row.names = FALSE)

# exporting the BoleChar model summary table...
write.csv(bolechar_summary_table,
          "results/final/tables/bolecharmodel_summary_table.csv",
          row.names = FALSE)

# exporting the CanopyPosition model summary table...
write.csv(canopyposition_summary_table,
          "results/final/tables/canopypositionmodel_summary_table.csv",
          row.names = FALSE)

# exporting the Mistletoe_count model summary table...
write.csv(mistletoe_summary_table,
          "results/final/tables/mistletoemodel_summary_table.csv",
          row.names = FALSE)




