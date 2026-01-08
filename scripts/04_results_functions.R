

# Title: Defining functions that will be used to summarize and communicate
# model results for the fire injury-fecundity project

# Authors: Nina Venuti and Andrew Latimer

# Acknowledgements: sections of code within Functions F (generate_
# predictions_df), H (predict_cones_at_thresholds), and K (create_
# percentchange_table) were adapted from Claire Tortorelli

# Script inputs: NA

# Script outputs: an ordered list of functions that will be used to 
# summarize fire injury-fecundity model results in subsequent scripts

library("tidyverse")



## PART I: Model summary functions ##

## Function A: create_basemodel_summary_table

# this function extracts key information on model sample size, model
# deviance explained, estimated degrees of freedom (edf) for smooth
# terms, and p-values for both smooth and non-smooth (parametric) terms
# from the model summary tables for each spp-yr "base model" (ie.
# TotalConesConservative_count ~ s(DBH) + s(CrownDamage) + 
# ti(DBH:CrownDamage) + SiteName, for each species and year sampled),
# and compiles this information into a single data.frame for publication.

# function input:
# model_names (a character vector) = a list of base model names from 
#     which to extract key model summary metrics. Note that the three
#     smooth terms (DBH, CrownDamage, and DBH:CrownDamage) included 
#     in each base model must be listed in the same order (with DBH listed 
#     first, CrownDamage listed second, and DBH:CrownDamage listed third) 
#     in each model for the function to work properly.

# function output:
# base_model_summary_table_wstars (a data.frame) = a table containing 
#     information on sample size, deviance explained, smooth term edf 
#     values, and smooth and parametric term p-values (with asterisks 
#     assigned to p-values based on estimated model term significance 
#     at alpha = 0.05) for each model specified in the function input,
#     model_names.

A_create_basemodel_summary_table <- function(model_names) {
  
  # find the length of the vector "model_names" (ie. find the number of 
  # models to be summarized in the basemodel_summary_table) 
  n_values = length(model_names)
  
  # define individual numeric vectors of this length to hold the edf 
  # values of model smooth terms, the p-values of model smooth and 
  # parametric terms, the sample size of each model, and the deviance
  # explained value of each model
  DBH_edf = vector(mode = "numeric", length = n_values)
  DBH_pvalue = vector(mode = "numeric", length = n_values)
  CrownDamage_edf = vector(mode = "numeric", length = n_values)
  CrownDamage_pvalue = vector(mode = "numeric", length = n_values)
  Intx_edf = vector(mode = "numeric", length = n_values)
  Intx_pvalue = vector(mode = "numeric", length = n_values)
  SiteName_pvalue = vector(mode = "numeric", length = n_values)
  model_samplesize = vector(mode = "numeric", length = n_values)
  model_dev.expl = vector(mode = "numeric", length = n_values)
  
  # define "model_summary" as a list of individual model summary tables
  model_summary = list()
  
  # run a for loop that calls each model summary table, in turn,  
  # extracts the edf and p-values for all smooth terms (DBH, CrownDamage, 
  # and the interaction term, DBH:CrownDamage), the p-value for the
  # parametric term, SiteName, the model sample size, and the model deviance
  # explained value, and saves these values in their appropriate numeric
  # vectors
  for(i in 1:n_values) {
    model_summary[[i]] = summary(get(model_names[i]))
    DBH_edf[i] = model_summary[[i]]$edf[1]
    DBH_pvalue[i] = model_summary[[i]]$s.pv[1]
    CrownDamage_edf[i] = model_summary[[i]]$edf[2]
    CrownDamage_pvalue[i] = model_summary[[i]]$s.pv[2]
    Intx_edf[i] = model_summary[[i]]$edf[3]
    Intx_pvalue[i] = model_summary[[i]]$s.pv[3]
    SiteName_pvalue[i] = model_summary[[i]]$pTerms.pv
    model_samplesize[i] = model_summary[[i]]$n
    model_dev.expl[i] = model_summary[[i]]$dev.expl
  }
  
  # combine these vectors into a single data.frame
  basemodel_summary_table = data.frame(
    model_name = model_names,
    n = model_samplesize,
    DBH_edf = DBH_edf,
    DBH_pvalue = DBH_pvalue,
    CrownDamage_edf = CrownDamage_edf,
    CrownDamage_pvalue = CrownDamage_pvalue,
    Intx_edf = Intx_edf,
    Intx_pvalue = Intx_pvalue,
    SiteName_pvalue = SiteName_pvalue,
    deviance_explained = model_dev.expl)
  
  # round edf, deviance explained, and p-values, assign asterisks
  # to p-values to represent estimated model term significance (at 
  # alpha = 0.05), and rearrange data.frame columns, as appropriate
  basemodel_summary_table_wstars = basemodel_summary_table %>% 
    mutate(across(all_of(c(ends_with("_edf"), ends_with("_explained"))),
                  .fns = ~round(.x, digits = 2))) %>% 
    mutate(across(all_of(c(ends_with("_edf"), ends_with("_explained"))),
                  .fns = ~format(.x, nsmall = 2))) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "***",
                                    .x >= 0.001 & .x < 0.01 ~ "**",
                                    .x >= 0.01 & .x < 0.05 ~ "*",
                                    .x >= 0.05 & .x < 0.1 ~ "'",
                                    .default = ""),
                  .names = "{.col}_stars")) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "<0.001",
                                    .default = format(round(.x, 
                                                            digits = 3), 
                                                      nsmall = 3)))) %>% 
    mutate(DBH_pvalue = paste(DBH_pvalue,
                              DBH_pvalue_stars,
                              sep = "")) %>% 
    mutate(CrownDamage_pvalue = paste(CrownDamage_pvalue,
                                      CrownDamage_pvalue_stars,
                                      sep = "")) %>% 
    mutate(Intx_pvalue = paste(Intx_pvalue,
                               Intx_pvalue_stars,
                               sep = "")) %>% 
    mutate(SiteName_pvalue = paste(SiteName_pvalue,
                                   SiteName_pvalue_stars,
                                   sep = "")) %>% 
    select(-c(ends_with("_stars")))
  
  # return the final basemodel_summary_table, with asterisks representing
  # estimated model term significance
  return(basemodel_summary_table_wstars)
  
}



## Function B: compile_basemodel_site_effects

# this function extracts coefficients and p-values for each of the
# four levels of the variable SiteName (the only non-smooth, parametric
# term included in each spp-yr base model) from each base model's summary
# table, and compiles these values into a single data.frame for 
# publication. Each level of the variable SiteName represents a field
# site at which cone counts were conducted - Intercept (which includes
# the effects of the site AUGUST), CAPLES, CREEK, and NORTH.

# function input:
# model_names (a character vector) = a list of [base] model names from 
#     which to extract metrics on site-level effects.

# function output:
# basemodel_siteeffects_table_wstars (a data.frame) = a table containing 
#     information on the coefficients and p-values affiliated with each
#     level of the variable SiteName, for each model specified in the
#     function input, model_names. Asterisks are assigned to p-values
#     based on each site level's estimated significance at alpha = 0.05.

B_compile_basemodel_site_effects <- function(model_names) {
  
  # note, metrics on site-level effects will be extracted first for
  # ABCO and PIPJ models, each of which involve four levels of the 
  # variable SiteName (Intercept/AUGUST, CAPLES, CREEK, and NORTH),
  # and second for PILA models, each of which involve three levels of
  # the variable SiteName (Intercept/AUGUST, CREEK, and NORTH)
  
  # subset the vector "model_names" according to species code
  ABCO_model_names = str_subset(model_names, pattern = "ABCO")
  PIPJ_model_names = str_subset(model_names, pattern = "PIPJ")
  PILA_model_names = str_subset(model_names, pattern = "PILA")
  
  # combine the ABCO_ and PIPJ_model_names vectors into a single vector
  ABCOPIPJ_model_names = c(ABCO_model_names, PIPJ_model_names)
  
  # find the length of the vector "ABCOPIPJ_model_names" (ie. 
  # find the number of ABCO and PIPJ models to be summarized in the
  # final basemodel_siteeffects_table) 
  ABCOPIPJ_n_values = length(ABCOPIPJ_model_names) 

  # define individual numeric vectors of this length to hold the
  # coefficient values and p-values for each level of the 
  # parametric variable SiteName relevant to both ABCO and PIPJ
  # models (ie. Intercept, CAPLES, CREEK, and NORTH)
  Intercept_coeff_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  Intercept_pvalue_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  CAPLES_coeff_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  CAPLES_pvalue_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  CREEK_coeff_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  CREEK_pvalue_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  NORTH_coeff_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  NORTH_pvalue_ABPJ = vector(mode = "numeric", length = ABCOPIPJ_n_values)
  
  # define "model_summary_ABPJ" as a list of individual ABCO and
  # PIPJ model summary tables
  model_summary_ABPJ = list()
  
  # run a for loop that calls each ABCO and PIPJ model summary table, 
  # in turn, extracts the coefficients and p-values for each level of the 
  # parametric term, SiteName, relevant for those species (ie. Intercept,
  # CAPLES, CREEK, and NORTH), and saves these values in their appropriate
  # numeric vectors
  for(i in 1:ABCOPIPJ_n_values) {
    model_summary_ABPJ[[i]] = summary(get(ABCOPIPJ_model_names[i]))
    Intercept_coeff_ABPJ[i] = model_summary_ABPJ[[i]]$p.coeff[1]
    Intercept_pvalue_ABPJ[i] = model_summary_ABPJ[[i]]$p.pv[1]
    CAPLES_coeff_ABPJ[i] = model_summary_ABPJ[[i]]$p.coeff[2]
    CAPLES_pvalue_ABPJ[i] = model_summary_ABPJ[[i]]$p.pv[2]
    CREEK_coeff_ABPJ[i] = model_summary_ABPJ[[i]]$p.coeff[3]
    CREEK_pvalue_ABPJ[i] = model_summary_ABPJ[[i]]$p.pv[3]
    NORTH_coeff_ABPJ[i] = model_summary_ABPJ[[i]]$p.coeff[4]
    NORTH_pvalue_ABPJ[i] = model_summary_ABPJ[[i]]$p.pv[4]
  }
  
  # bind these vectors together into a single data.frame, summarizing
  # site-level effects in ABCO and PIPJ models
  ABCOPIPJ_siteeffects_table = data.frame(
    model_name = ABCOPIPJ_model_names,
    Intercept_coeff = Intercept_coeff_ABPJ,
    Intercept_pvalue = Intercept_pvalue_ABPJ,
    CAPLES_coeff = CAPLES_coeff_ABPJ,
    CAPLES_pvalue = CAPLES_pvalue_ABPJ,
    CREEK_coeff = CREEK_coeff_ABPJ,
    CREEK_pvalue = CREEK_pvalue_ABPJ,
    NORTH_coeff = NORTH_coeff_ABPJ,
    NORTH_pvalue = NORTH_pvalue_ABPJ)
  
  # repeat this process for PILA models
  # find the length of the vector "PILA_model_names" (ie. find the 
  # number of PILA models to be summarized in the final basemodel_
  # siteeffects_table) 
  PILA_n_values = length(PILA_model_names)
  
  # define individual numeric vectors of this length to hold the
  # coefficient values and p-values for each level of the 
  # parametric variable SiteName relevant to PILA models (ie. Intercept,
  # CREEK, and NORTH)
  Intercept_coeff_PILA = vector(mode = "numeric", length = PILA_n_values)
  Intercept_pvalue_PILA = vector(mode = "numeric", length = PILA_n_values)
  CREEK_coeff_PILA = vector(mode = "numeric", length = PILA_n_values)
  CREEK_pvalue_PILA = vector(mode = "numeric", length = PILA_n_values)
  NORTH_coeff_PILA = vector(mode = "numeric", length = PILA_n_values)
  NORTH_pvalue_PILA = vector(mode = "numeric", length = PILA_n_values)
  
  # define "model_summary_PILA" as a list of individual PILA model
  # summary tables
  model_summary_PILA = list()
  
  # run a for loop that calls each PILA model summary table, in turn,   
  # extracts the coefficients and p-values for each level of the parametric
  # term, SiteName, relevant for this species (ie. Intercept, CREEK, and
  # NORTH), and saves these values in their appropriate numeric vectors
  for(i in 1:PILA_n_values) {
    model_summary_PILA[[i]] = summary(get(PILA_model_names[i]))
    Intercept_coeff_PILA[i] = model_summary_PILA[[i]]$p.coeff[1]
    Intercept_pvalue_PILA[i] = model_summary_PILA[[i]]$p.pv[1]
    CREEK_coeff_PILA[i] = model_summary_PILA[[i]]$p.coeff[2]
    CREEK_pvalue_PILA[i] = model_summary_PILA[[i]]$p.pv[2]
    NORTH_coeff_PILA[i] = model_summary_PILA[[i]]$p.coeff[3]
    NORTH_pvalue_PILA[i] = model_summary_PILA[[i]]$p.pv[3]
  }
  
  # bind these vectors together into a single data.frame, summarizing
  # site-level effects in PILA models
  PILA_siteeffects_table = data.frame(
    model_name = PILA_model_names,
    Intercept_coeff = Intercept_coeff_PILA,
    Intercept_pvalue = Intercept_pvalue_PILA,
    CREEK_coeff = CREEK_coeff_PILA,
    CREEK_pvalue = CREEK_pvalue_PILA,
    NORTH_coeff = NORTH_coeff_PILA,
    NORTH_pvalue = NORTH_pvalue_PILA)
  
  # bind the ABCOPIPJ_siteeffects_table and the PILA_siteeffects_table
  # together into a single data.frame
  basemodel_siteeffects_table = bind_rows(ABCOPIPJ_siteeffects_table,
                                          PILA_siteeffects_table)
  
  # round all coefficient values and p-values, assign asterisks to 
  # p-values to represent estimated site level significance (at alpha =
  # 0.05), and rearrange data.frame columns, as appropriate
  basemodel_siteeffects_table_wstars = basemodel_siteeffects_table %>% 
    mutate(across(all_of(c(ends_with("_coeff"))),
                  .fns = ~round(.x, digits = 2))) %>% 
    mutate(across(all_of(c(ends_with("_coeff"))),
                  .fns = ~format(.x, nsmall = 2))) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "***",
                                    .x >= 0.001 & .x < 0.01 ~ "**",
                                    .x >= 0.01 & .x < 0.05 ~ "*",
                                    .x >= 0.05 & .x < 0.1 ~ "'",
                                    .default = ""),
                  .names = "{.col}_stars")) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "<0.001",
                                    .default = format(round(.x, 
                                                            digits = 3), 
                                                      nsmall = 3)))) %>% 
    mutate(Intercept_pvalue = paste(Intercept_pvalue,
                                    Intercept_pvalue_stars,
                                    sep = "")) %>% 
    mutate(CAPLES_pvalue = paste(CAPLES_pvalue,
                                 CAPLES_pvalue_stars,
                                 sep = "")) %>% 
    mutate(CREEK_pvalue = paste(CREEK_pvalue,
                                CREEK_pvalue_stars,
                                sep = "")) %>% 
    mutate(NORTH_pvalue = paste(NORTH_pvalue,
                                NORTH_pvalue_stars,
                                sep = "")) %>% 
    select(-c(ends_with("_stars")))
  
  # return the final basemodel_siteeffects_table, with asterisks 
  # representing estimated site level significance
  return(basemodel_siteeffects_table_wstars)
  
}



## Function C: compile_bolechar_model_estimates

# this function extracts key information on model sample size, model
# deviance explained, and the estimated degrees of freedom (edf) and
# p-values for the smooth term BoleCharMeanHeight_m (one of the three
# "auxiliary" variables added, sequentially, to each spp-yr base model) 
# from relevant model summary tables, and compiles this information  
# into a single data.frame for publication. 

# function input:
# bolechar_model_names (a character vector) = a list of model names from 
#     which to extract key model summary metrics; models must include 
#     BoleCharMeanHeight_m as a smooth term, and BoleCharMeanHeight_m
#     must be listed as the fourth smooth term in each model (ie. following
#     DBH, CrownDamage, and DBH:CrownDamage) for the function to work
#     properly.

# function output:
# bolechar_summary_table_wstars (a data.frame) = a table containing 
#     information on model sample size, model deviance explained, and
#     the edf and p-values affiliated with the smooth predictor 
#     BoleCharMeanHeight_m, for each model specified in the function
#     input, bolechar_model_names. Asterisks are assigned to p-values
#     based on estimated model term (BoleChar) significance at alpha = 
#     0.05.

C_compile_bolechar_model_estimates <- function(bolechar_model_names) {
  
  # find the length of the vector "bolechar_model_names" (ie. find
  # the number of models to be summarized in the bolechar_summary_table)
  n_values = length(bolechar_model_names)
  
  # define individual numeric vectors of this length to hold the 
  # edf and p-value of the BoleChar term, the sample size, and the 
  # deviance explained value for each model
  BoleChar_edf = vector(mode = "numeric", length = n_values)
  BoleChar_pvalue = vector(mode = "numeric", length = n_values)
  model_samplesize = vector(mode = "numeric", length = n_values)
  model_dev.expl = vector(mode = "numeric", length = n_values)
  
  # define "model_summary" as a list of individual model summary tables
  model_summary = list()
  
  # run a for loop that calls each model summary table, in turn, 
  # extracts the edf and p-value of the BoleChar term, the model
  # sample size, and the model deviance explained value, and saves
  # these values in their appropriate numeric vectors
  for(i in 1:n_values) {
    model_summary[[i]] = summary(get(bolechar_model_names[i]))
    BoleChar_edf[i] = model_summary[[i]]$edf[4]
    BoleChar_pvalue[i] = model_summary[[i]]$s.pv[4]
    model_samplesize[i] = model_summary[[i]]$n
    model_dev.expl[i] = model_summary[[i]]$dev.expl
  }
  
  # combine these vectors into a single data.frame
  bolechar_summary_table = data.frame(
    model_name = bolechar_model_names,
    n = model_samplesize,
    BoleChar_edf = BoleChar_edf,
    BoleChar_pvalue = BoleChar_pvalue,
    deviance_explained = model_dev.expl)
  
  # round edf, deviance explained, and p-values, assign asterisks
  # to p-values to represent estimated model term (BoleChar) significance
  # (at alpha = 0.05), and rearrange data.frame columns, as appropriate
  bolechar_summary_table_wstars = bolechar_summary_table %>% 
    mutate(across(all_of(c(ends_with("_edf"), ends_with("_explained"))),
           .fns = ~round(.x, digits = 2))) %>% 
    mutate(across(all_of(c(ends_with("_edf"), ends_with("_explained"))),
           .fns = ~format(.x, nsmall = 2))) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "***",
                                    .x >= 0.001 & .x < 0.01 ~ "**",
                                    .x >= 0.01 & .x < 0.05 ~ "*",
                                    .x >= 0.05 & .x < 0.1 ~ "'",
                                    .default = ""),
                  .names = "{.col}_stars")) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "<0.001",
                                    .default = format(round(.x, 
                                                            digits = 3), 
                                                      nsmall = 3)))) %>% 
    mutate(BoleChar_pvalue = paste(BoleChar_pvalue,
                                   BoleChar_pvalue_stars,
                                   sep = "")) %>% 
    select(-c(ends_with("_stars")))
  
  # return the bolechar_summary_table, with asterisks representing
  # estimated model term (BoleChar) significance
  return(bolechar_summary_table_wstars)
  
}



## Function D: compile_canopyposition_model_estimates

# this function extracts key information on model sample size, model
# deviance explained, and the coefficients and p-values for the non-
# smooth, parametric term CanopyPosition (one of the three "auxiliary" 
# variables added, sequentially, to each spp-yr base model) from relevant
# model summary tables, and compiles this information into a single 
# data.frame for publication. 

# function input:
# canopyposition_model_names (a character vector) = a list of model names 
#     from which to extract key model summary metrics; models must include 
#     CanopyPosition as a parametric term, and CanopyPosition must be 
#     listed as the first parametric term (ie. prior to SiteName) in each
#     model for the function to work properly.

# function output:
# canopyposition_summary_table_wstars (a data.frame) = a table containing 
#     information on model sample size, model deviance explained, and
#     the coefficients and p-values affiliated with the second level
#     of the parametric predictor CanopyPosition (ie. "CanopyPositionDO" - 
#     for dominant trees - representing the estimated effect of canopy
#     dominance, relative to canopy co-dominance, on cone production) for 
#     each model specified in the function input, canopyposition_model_names.
#     Asterisks are assigned to p-values based on estimated parametric 
#     term level (CanopyPositionDO) significance at alpha = 0.05.

D_compile_canopyposition_model_estimates <- function(
    canopyposition_model_names) {
  
  # find the length of the vector "canopyposition_model_names" (ie. find
  # the number of models to be summarized in the canopyposition_summary_
  # table)
  n_values = length(canopyposition_model_names)
  
  # define individual numeric vectors of this length to hold the 
  # coefficient value and p-value of CanopyPositionDO (the second level
  # of the parametric term CanopyPosition), the sample size, and the 
  # deviance explained value for each model
  CanopyPositionDO_coeff = vector(mode = "numeric", length = n_values)
  CanopyPositionDO_pvalue = vector(mode = "numeric", length = n_values)
  model_samplesize = vector(mode = "numeric", length = n_values)
  model_dev.expl = vector(mode = "numeric", length = n_values)
  
  # define "model_summary" as a list of individual model summary tables
  model_summary = list()
  
  # run a for loop that calls each model summary table, in turn, 
  # extracts the coefficient and p-value of CanopyPositionDO, the model
  # sample size, and the model deviance explained value, and saves
  # these values in their appropriate numeric vectors
  for(i in 1:n_values) {
    model_summary[[i]] = summary(get(canopyposition_model_names[i]))
    CanopyPositionDO_coeff[i] = model_summary[[i]]$p.coeff[2]
    CanopyPositionDO_pvalue[i] = model_summary[[i]]$p.pv[2]
    model_samplesize[i] = model_summary[[i]]$n
    model_dev.expl[i] = model_summary[[i]]$dev.expl
  }
  
  # combine these vectors into a single data.frame
  canopyposition_summary_table = data.frame(
    model_name = canopyposition_model_names,
    n = model_samplesize,
    CanopyPositionDO_coeff = CanopyPositionDO_coeff,
    CanopyPositionDO_pvalue = CanopyPositionDO_pvalue,
    deviance_explained = model_dev.expl)
  
  # round coefficients, deviance explained values, and p-values, assign
  # asterisks to p-values to represent estimated CanopyPosition level
  # significance (at alpha = 0.05), and rearrange data.frame columns,
  # as appropriate
  canopyposition_summary_table_wstars = canopyposition_summary_table %>% 
    mutate(across(all_of(c(ends_with("_coeff"), ends_with("_explained"))),
                  .fns = ~round(.x, digits = 2))) %>% 
    mutate(across(all_of(c(ends_with("_coeff"), ends_with("_explained"))),
                  .fns = ~format(.x, nsmall = 2))) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "***",
                                    .x >= 0.001 & .x < 0.01 ~ "**",
                                    .x >= 0.01 & .x < 0.05 ~ "*",
                                    .x >= 0.05 & .x < 0.1 ~ "'",
                                    .default = ""),
                  .names = "{.col}_stars")) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "<0.001",
                                    .default = format(round(.x, 
                                                            digits = 3), 
                                                      nsmall = 3)))) %>% 
    mutate(CanopyPositionDO_pvalue = paste(CanopyPositionDO_pvalue,
                                           CanopyPositionDO_pvalue_stars,
                                           sep = "")) %>% 
    select(-c(ends_with("_stars")))
  
  # return the canopyposition_summary_table, with asterisks representing
  # estimated CanopyPosition level significance
  return(canopyposition_summary_table_wstars)
  
}



## Function E: compile_mistletoe_model_estimates

# this function extracts key information on model sample size, model
# deviance explained, and the estimated degrees of freedom (edf) and 
# p-values for the sample year-specific (yr1 or yr2) smooth term
# MistletoeClumpsOrBranches_count (one of the three "auxiliary" variables 
# added, sequentially, to each spp-yr base model) from relevant model 
# summary tables, and compiles this information into a single data.frame 
# for publication.  

# function input:
# mistletoe_model_names (a character vector) = a list of model names from 
#     which to extract key model summary metrics; models must include 
#     MistletoeClumpsOrBranches_count_yr1 or _yr2 as a smooth term, 
#     and MistletoeClumpsOrBranches_count must be listed as the fourth 
#     smooth term in each model (ie. following DBH, CrownDamage, and 
#     DBH:CrownDamage) for the function to work properly. Note that this
#     list of model names should only include ABCO or PIPJ models, as the
#     incidence of mistletoe infestation in PILA was too rare to robustly
#     assess its impact on reproductive output in that species. 

# function output:
# mistletoe_summary_table_wstars (a data.frame) = a table containing 
#     information on model sample size, model deviance explained, and
#     the edf and p-values affiliated with the smooth predictor 
#     MistletoeClumpsOrBranches_count_yr1 or _yr2, for each model 
#     specified in the function input, mistletoe_model_names. Asterisks 
#     are assigned to p-values based on estimated model term (Mistletoe_
#     count) significance at alpha = 0.05.

E_compile_mistletoe_model_estimates <- function(mistletoe_model_names) {
  
  # find the length of the vector "mistletoe_model_names" (ie. find
  # the number of models to be summarized in the mistletoe_summary_table)
  n_values = length(mistletoe_model_names)
  
  # define individual numeric vectors of this length to hold the 
  # edf and p-value of the Mistletoe term, the sample size, and the 
  # deviance explained value for each model
  Mistletoe_edf = vector(mode = "numeric", length = n_values)
  Mistletoe_pvalue = vector(mode = "numeric", length = n_values)
  model_samplesize = vector(mode = "numeric", length = n_values)
  model_dev.expl = vector(mode = "numeric", length = n_values)
  
  # define "model_summary" as a list of individual model summary tables
  model_summary = list()
  
  # run a for loop that calls each model summary table, in turn, 
  # extracts the edf and p-value of the Mistletoe term, the model
  # sample size, and the model deviance explained value, and saves
  # these values in their appropriate numeric vectors
  for(i in 1:n_values) {
    model_summary[[i]] = summary(get(mistletoe_model_names[i]))
    Mistletoe_edf[i] = model_summary[[i]]$edf[4]
    Mistletoe_pvalue[i] = model_summary[[i]]$s.pv[4]
    model_samplesize[i] = model_summary[[i]]$n
    model_dev.expl[i] = model_summary[[i]]$dev.expl
  }
  
  # combine these vectors into a single data.frame
  mistletoe_summary_table = data.frame(
    model_name = mistletoe_model_names,
    n = model_samplesize,
    Mistletoe_edf = Mistletoe_edf,
    Mistletoe_pvalue = Mistletoe_pvalue,
    deviance_explained = model_dev.expl)
  
  # round edf, deviance explained, and p-values, assign asterisks
  # to p-values to represent estimated model term (Mistletoe_count)
  # significance (at alpha = 0.05), and rearrange data.frame columns,
  # as appropriate
  mistletoe_summary_table_wstars = mistletoe_summary_table %>%
    mutate(across(all_of(c(ends_with("_edf"), ends_with("_explained"))),
                  .fns = ~round(.x, digits = 2))) %>% 
    mutate(across(all_of(c(ends_with("_edf"), ends_with("_explained"))),
                  .fns = ~format(.x, nsmall = 2))) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "***",
                                    .x >= 0.001 & .x < 0.01 ~ "**",
                                    .x >= 0.01 & .x < 0.05 ~ "*",
                                    .x >= 0.05 & .x < 0.1 ~ "'",
                                    .default = ""),
                  .names = "{.col}_stars")) %>% 
    mutate(across(all_of(c(ends_with("_pvalue"))),
                  .fns = ~case_when(.x < 0.001 ~ "<0.001",
                                    .default = format(round(.x, 
                                                            digits = 3), 
                                                      nsmall = 3)))) %>% 
    mutate(Mistletoe_pvalue = paste(Mistletoe_pvalue,
                                    Mistletoe_pvalue_stars,
                                    sep = "")) %>% 
    select(-c(ends_with("_stars")))
    
  # return the mistletoe_summary_table, with asterisks representing
  # estimated model term (Mistletoe) significance
  return(mistletoe_summary_table_wstars)
  
}






## PART II: Results visualization and communication functions ## 

## Function F: generate_predictions_df

# this function uses a spp-yr "base model" (ie. TotalConesConservative_
# count ~ s(DBH) + s(CrownDamage) + ti(DBH:CrownDamage) + SiteName, for a
# specific species and sample year) to predict total conservative
# cone counts at three distinct DBH sizes (the 10th, 50th, and 90th
# percentile-sized trees), along the full spectrum of CrownDamage, and  
# at the median field site for the corresponding spp-yr dataset. It then
# combines model fit and predictor values into a single data.frame for 
# use in plotting and communicating model results.

# function inputs:
# model (a fitted GAM) = a spp-yr base model with which to predict total 
#     conservative cone counts. 
# sppyr_data (a data.frame) = a spp-yr-specific table of field-collected 
#     data on tree size, fire injury, and fecundity used to fit 
#     the spp-yr base model referenced above.
# median_site (a character string) = the name of the median cone production
#     site (ie. "AUGUST", "CREEK", "CAPLES", or "NORTH") for the relevant
#     spp-yr, as determined via visual inspection of the spp-yr base model's
#     summary table.

# function output:
# pred_df (a data.frame) = a table of predicted mean cone counts, predicted
#     confidence intervals around those mean cone counts, and the predictor
#     values used to generate such predictions. Predicted mean cone counts 
#     and confidence intervals are represented on both the log-transformed
#     and back-transformed (ie. response) scales. Predictor variables 
#     are represented as both scaled and unscaled (ie. original unit) values.

F_generate_predictions_df <- function(model, sppyr_data, median_site) {
  
  # create a data.frame with 3000 rows and 3 columns, populated with 
  # 1000 samples of the full range of crown volume damage (from min to 
  # max values of crown volume damage) contained in sppyr_data,
  # DBH values for the 10th, 50th, and 90th percentile-sized trees
  # contained in sppyr_data (ie. three distinct DBH values for each of
  # the 1000 samples of crown vol damage), and the median_site for the
  # predictor variables CrownDamage_percentvolume_scaled, DBH_cm_scaled, 
  # and SiteName, respectively
  pred_data = expand_grid(
    CrownDamage_percentvolume_scaled = seq(from = min(
      sppyr_data$CrownDamage_percentvolume_scaled), to = max(
        sppyr_data$CrownDamage_percentvolume_scaled),
      length.out = 1000),
    DBH_cm_scaled = quantile(sppyr_data$DBH_cm_scaled,
                             probs = c(0.1, 0.5, 0.9)),
    SiteName = as.factor(median_site))
  
  # predict total conservative cone counts using pred_data and the 
  # specified spp-yr model, and bind these predictions to pred_data. Then,
  # add columns for the unscaled predictor variables, as well as the
  # back-transformed model predictions (fit) and the upper and lower
  # limits of the predicted confidence intervals (ie. model fit +/- 1SE)
  # to enable plotting on more understandable axes. Lastly, make 
  # "DBH_cm_unscaled" a factor variable, to enable color coding by tree
  # size while plotting.
  pred_df = predict(object = model,
                    newdata = pred_data,
                    se.fit = TRUE) %>% 
    as_tibble() %>% 
    bind_cols(pred_data) %>% 
    mutate(DBH_cm_unscaled = (DBH_cm_scaled * DBH_sd) + DBH_mean) %>% 
    mutate(CrownDamage_percentvolume_unscaled = 
             (CrownDamage_percentvolume_scaled * CrownDamage_sd) +
             CrownDamage_mean) %>% 
    mutate(fit_backtransformed = exp(fit)) %>% 
    mutate(CI_upperlimit = exp(fit + se.fit)) %>% 
    mutate(CI_lowerlimit = exp(fit - se.fit)) %>% 
    mutate(DBH_cm_unscaled = as.factor(DBH_cm_unscaled))
  
  # return pred_df
  return(pred_df)
  
}



## Function G1: calculate_finite_second_derivative

# this function calculates the finite second derivative along a curve
# defined by specific x and y values, and produces a data.frame with three
# columns: one that contains x values, one that contains finite
# first derivative (ie. slope) values along the curve, and one that contains
# finite second derivative (ie. rate of change of the slope) values 
# along the curve.

# function inputs:
# x (a numeric vector) = a vector of x (ie. "run") values along which
#     a curve is plotted.
# y (a numeric vector) = a vector of y (ie. "rise") values that correspond
#     to the given x values. 
# smooth (a numeric value) = a number between 0 and 1, used to define
#     the width of a rolling window across which the mean finite second
#     derivative is calculated. While the default "smooth" value is set
#     to NULL, it can be helpful to specify a smoothing factor if 
#     small-scale fluctuations in a curve are producing misleading
#     finite second derivative values, and if averaging over such
#     small fluctuations would help with trend interpretation. Smaller 
#     smooth values (closer to 0) produce narrower rolling windows; larger
#     smooth values (closer to 1) produce wider rolling windows.

# function output:
# a data.frame with three columns, one that contains the x values 
#     for which finite second derivatives have been calculated, one
#     that contains the calculated finite first derivative (ie. slope) 
#     values along the curve, and one that contains the calculated finite 
#     second derivative (ie. rate of change of the slope) values along 
#     the curve.

G1_calculate_finite_second_derivative <- function(x, y, smooth = NULL) {
  
  # require the R package "zoo" to run this function; "zoo" contains
  # the function "rollmean," used below to find rolling window averages 
  # of deriv2
  require("zoo")
  
  # calculate the finite first derivative (ie. slope) between consecutive
  # points along the curve
  deriv1 = diff(y)/diff(x)
  
  # calculate the "rise" of the second derivative (ie. the difference
  # in consecutive values of deriv1)
  deriv2_rise = diff(deriv1)
  
  # calculate the "run" of the second derivative (ie. the difference 
  # in x values that are two indices/positions apart (indicated by
  # the "lag = 2" argument within diff()); these are the x values that
  # are relevant to the consecutive values of deriv1 used to calculate 
  # deriv2_rise)
  deriv2_run = diff(x, lag = 2)
  
  # calculate the finite second derivative
  deriv2 = deriv2_rise/deriv2_run 
  
  # find the x value that corresponds to each deriv2 value (defined here
  # as the midpoint of deriv2_run)
  deriv2_run_midpoint = x[1:(length(x) - 2)] + deriv2_run/2
  
  # if a smooth value is provided in the function inputs (ie. if
  # smooth is NOT NULL), calculate deriv2 as the mean finite second
  # derivative within a rolling window, whose width is determined
  # by length(deriv2)*smooth
  if(!is.null(smooth)) deriv2 = zoo::rollmean(deriv2, 
                                              k = length(deriv2)*smooth, 
                                              fill = NA)
  
  # return a data.frame containing three numeric columns: deriv2_run_
  # midpoint, deriv1, and deriv2; add NAs to the first, or first and second, 
  # rows in the data.frame to account for the fact that deriv1 and
  # deriv2 cannot be calculated for the first, and first and second,
  # values of x, respectively
  return(tibble(deriv2_run_midpoint = c(NA, NA, deriv2_run_midpoint),
                deriv1 = c(NA, deriv1),
                deriv2 = c(NA, NA, deriv2)))
  
}



## Function G2: identify_thresholds_minderiv2

# this function employs Function G1 (calculate_finite_second_derivative)
# to find the finite second derivative along a predicted cone
# count curve. It then identifies the minimum, negative value of the finite
# second derivative along that curve (ie. the point at which the slope of 
# the predicted cone count curve is decreasing most rapidly), and extracts 
# the x value (or the CrownDamage value) affiliated with this minimum
# point. This x value is defined here as the "threshold value"
# of the curve - ie. the level of crown volume damage at which cone 
# counts are decreasing most rapidly.

# function inputs:
# x_values (a numeric vector) = a vector of x (ie. "run") values along
#     which a curve is plotted; for our purposes, x_values will be
#     scaled values of the predictor variable CrownDamage_percentvolume.
# curve_values (a numeric vector) = a vector of y (ie. "rise") values
#     that correspond to the given x_values; for our purposes, curve_
#     values will be log-transformed values of predicted cone counts
#     (ie. model fits).
# smooth (a numeric value) = a number between 0 and 1, used to define
#     the width of a rolling window across which the mean finite second
#     derivative is calculated. See the description of the inputs for 
#     Function G1, above, for further details.

# function output: 
# threshold_val (a numeric value) = the x_value at which the minimum,
#     negative finite second derivative occurs; for our purposes, the 
#     threshold_val represents the level of crown volume damage at which 
#     cone counts begin to decrease rapidly.

G2_identify_thresholds_minderiv2 <- function(x_values, curve_values, 
                                             smooth = NULL) {
  
  # calculate the finite second derivative using specified x_values 
  # and curve_values
  deriv2_df = G1_calculate_finite_second_derivative(x = x_values,
                                                    y = curve_values,
                                                    smooth = smooth)
  
  # locate the minimum, negative finite second derivative value along
  # the curve, and extract the x_value affiliated with that minimum,
  # negative finite second derivative value (this is the threshold_val);
  # as we define the threshold_val as the point at which cone counts
  # are decreasing most rapidly, we restrict our search for the minimum
  # finite second derivative value to negative values only; if 
  # there are no negative finite second derivative values along the
  # curve, there is no threshold_val for that curve (and threshold_val
  # is returned as "NA" below)
  
  # if the minimum finite second derivative value calculated for the
  # curve is negative...
  if (min(deriv2_df$deriv2, na.rm = T) < 0) {
    
    # find the index (ie. the location) of the minimum second 
    # derivative value in deriv2_df
    min_index = which.min(deriv2_df$deriv2)
    
    # and find the x_value affiliated with that minimum second 
    # derivative value (this is the "threshold value")
    threshold_val = deriv2_df$deriv2_run_midpoint[min_index]
    
  }
  
  # if not, define the threshold value as "NA"
  else threshold_val = NA
  
  # return the identified threshold value
  return(threshold_val)
  
}



## Function H: predict_cones_at_thresholds

# this function predicts total conservative cone counts at identified
# threshold values (ie. CrownDamage_percentvolume_scaled values), their
# respective DBH values (ie. 10th, 50th, and 90th percentile-sized
# trees), and the median fecundity site for a given spp-yr base model
# (ie. TotalConesConservative_count ~ s(DBH) + s(CrownDamage) + 
# ti(DBH:CrownDamage) + SiteName, for a specific species and sample year).
# It returns a data.frame with model fit and predictor values that will
# be used during plotting.

# function inputs:
# model (a fitted GAM) = a spp-yr base model with which to predict total 
#     conservative cone counts.
# median_site (a character string) = the name of the median fecundity
#     site (ie. "AUGUST", "CREEK", "CAPLES", or "NORTH") for the relevant
#     spp-yr, as determined via visual inspection of the spp-yr base model's
#     summary table.
# threshold_table (a data.frame) = a table containing threshold values
#     identified using Functions G1 and G2 (described above) for 
#     predicted cone count curves (typically, for each of three distinct  
#     curves representing the 10th, 50th, and 90th percentile-sized trees) 
#     for the relevant spp-yr. 

# function output:
# pred_df_thresholds (a data.frame) = a table of predicted mean cone
#     counts (and corresponding confidence intervals) at identified 
#     threshold values. Predicted counts and confidence intervals
#     are represented on both the log-transformed and back-transformed
#     (ie. response) scales. Predictor variables are represented as 
#     both scaled and unscaled (ie. original unit) values.

H_predict_cones_at_thresholds <- function(model, median_site,
                                          threshold_table) {
  
  # create a data.frame that includes the identified threshold values,
  # their corresponding DBH values (ie. the 10th, 50th, and 90th percentile-
  # sized trees in the relevant spp-yr dataset), and the median_site
  pred_data_thresholds = threshold_table %>% 
    rename("CrownDamage_percentvolume_scaled" = threshold_val) %>% 
    mutate(SiteName = as.factor(median_site))
  
  # predict total conservative cone counts using pred_data_thresholds and
  # the specified spp-yr model, bind the predictions to
  # pred_data_thresholds, unscale the predictor variables, backtransform
  # the model fit and the predicted upper and lower bounds of the 
  # confidence intervals (ie. fit +/- 1SE), and transform "DBH_cm_unscaled"
  # into a factor variable to facilitate color coding during plotting
  pred_df_thresholds = predict(object = model,
                               newdata = pred_data_thresholds,
                               se.fit = TRUE) %>% 
    as_tibble() %>% 
    bind_cols(pred_data_thresholds) %>% 
    mutate(DBH_cm_unscaled = (DBH_cm_scaled * DBH_sd) + DBH_mean) %>% 
    mutate(CrownDamage_percentvolume_unscaled = 
             (CrownDamage_percentvolume_scaled * CrownDamage_sd) +
             CrownDamage_mean) %>% 
    mutate(fit_backtransformed = exp(fit)) %>% 
    mutate(CI_upperlimit = exp(fit + se.fit)) %>% 
    mutate(CI_lowerlimit = exp(fit - se.fit)) %>% 
    mutate(DBH_cm_unscaled = as.factor(DBH_cm_unscaled))
  
  # return pred_df_thresholds
  return(pred_df_thresholds)
  
}



## Function I: create_lookup_table

# this function creates a lookup table with key DBH, CrownDamage, and 
# SiteName values for a given spp-yr base model. DBH values represent 
# the 10th, 50th, and 90th percentile-sized trees in the relevant spp-yr 
# dataset; CrownDamage values are scaled values of 5%, 50%, and 95% 
# crown volume damage, as well as the threshold value (ie. the % crown
# volume damage at which the minimum, negative second derivative occurs) 
# for each tree size; SiteName is the median fecundity site for the
# relevant spp-yr. This lookup table will be used as an input for 
# subsequent functions (including Functions J and K) below.

# function inputs:
# sppyr_data (a data.frame) = a spp-yr-specific table of field-collected 
#     data on tree size, fire injury, and fecundity used to fit 
#     the relevant spp-yr base model.
# median_site (a character string) = the name of the median fecundity
#     site (ie. "AUGUST", "CREEK", "CAPLES", or "NORTH") for the relevant
#     spp-yr, as determined via visual inspection of the spp-yr base model's
#     summary table.
# threshold_table (a data.frame) = a table containing threshold values
#     for the predicted cone count curves of the 10th, 50th, and 90th
#     percentile-sized trees in the spp-yr. Threshold values are
#     identified using Functions G1 and G2 (described above). 

# function outputs: 
# lookup_table (a data.frame) = a table containing key DBH, CrownDamage,
#     and SiteName values of interest, for use in subsequent functions.
 
I_create_lookup_table <- function(sppyr_data, median_site, threshold_table) {
  
  # find scaled values of key crown volume damage levels of interest
  CrownDamage_05percent_scaled = (5 - CrownDamage_mean)/CrownDamage_sd
  CrownDamage_50percent_scaled = (50 - CrownDamage_mean)/CrownDamage_sd
  CrownDamage_95percent_scaled = (95 - CrownDamage_mean)/CrownDamage_sd
  
  # extract threshold values affiliated with each tree size (ie.
  # 10th, 50th, and 90th percentile-sized trees)
  smalltree_threshold_scaled = filter(threshold_table, DBH_cm_scaled == 
                                        min(DBH_cm_scaled))$threshold_val
  mediumtree_threshold_scaled = filter(threshold_table, DBH_cm_scaled ==
                                         median(DBH_cm_scaled))$threshold_val
  bigtree_threshold_scaled = filter(threshold_table, DBH_cm_scaled ==
                                      max(DBH_cm_scaled))$threshold_val
  
  # create a lookup table with the three DBH values of interest (ie. the 
  # 10th, 50th, and 90th percentile-sized trees in sppyr_data), the 
  # crown volume damage levels of interest (ie. the scaled values of 5%, 
  # 50%, and 95% crown volume damage, as well as the threshold value for 
  # each tree size), the median reproductive site for the relevant spp-yr
  # (median_site), and a column containing shorthand labels for each 
  # level of crown volume damage (ie. "5", "50", "threshold", and "95")
  lookup_table = tibble(
    DBH_cm_scaled = rep(quantile(sppyr_data$DBH_cm_scaled,
                                 probs = c(0.1, 0.5, 0.9)), each = 4),
    CrownDamage_label = rep(c("5", "50", "threshold", "95"), times = 3),
    CrownDamage_percentvolume_scaled = case_when(
      CrownDamage_label == "5" ~ CrownDamage_05percent_scaled,
      CrownDamage_label == "50" ~ CrownDamage_50percent_scaled,
      CrownDamage_label == "95" ~ CrownDamage_95percent_scaled,
      CrownDamage_label == "threshold" & DBH_cm_scaled == min(DBH_cm_scaled) ~
        smalltree_threshold_scaled,
      CrownDamage_label == "threshold" & DBH_cm_scaled == 
        median(DBH_cm_scaled) ~ mediumtree_threshold_scaled,
      CrownDamage_label == "threshold" & DBH_cm_scaled == max(DBH_cm_scaled) ~
        bigtree_threshold_scaled,
      .default = NA),
    SiteName = as.factor(median_site))
  
  # return the lookup table
  return(lookup_table)
  
}



## Function J: pivot_lookup_table_wider

# this function creates a wide version of the lookup table produced 
# by Function I (create_lookup_table), so that each DBH value of 
# interest occurs in its own row, and each CrownDamage level of
# of interest (corresponding to those DBH values) occurs in its own 
# column. This wide-version lookup table also includes a column
# containing the attributes of the specified spp-yr base model, 
# and will be used as an input for subsequent functions (e.g.,
# Function L2) below.

# function inputs:
# model (a fitted GAM) = a spp-yr base model.
# lookup_table (a data.frame) = a table containing key DBH, CrownDamage,
#     and SiteName values of interest for the specified spp-yr base model.
#     DBH values represent the 10th, 50th, and 90th percentile-sized trees
#     in the relevant spp-yr dataset; CrownDamage values are scaled values 
#     of 5%, 50%, and 95% crown volume damage, as well as the threshold value 
#     (ie. the % crown volume damage at which the minimum, negative second 
#     derivative occurs) for each tree size; SiteName is the median 
#     fecundity site for the relevant spp-yr.

# function output:
# lookup_table_wide (a data.frame) = a wide version of the lookup table
#    produced using Function I (described above), with DBH values represented
#    in rows, CrownDamage values represented in columns, and an additional
#    column containing base model attributes. This object will be used
#    in subsequent functions (e.g., Function L2), below.

J_pivot_lookup_table_wider <- function(model, lookup_table) {
  
  # pivot lookup table wider, so that each crown volume damage level
  # of interest (5%, 50%, 95%, and the threshold value) occurs in its own
  # column; add the column "model_to_test" to hold lists of spp-yr base
  # model attributes 
  lookup_table_wide = lookup_table %>% 
    group_by(DBH_cm_scaled) %>% 
    pivot_wider(id_cols = c("DBH_cm_scaled", "SiteName"),
                names_from = CrownDamage_label,
                values_from = CrownDamage_percentvolume_scaled,
                names_prefix = "CrownDamage_percentvolume_scaled_") %>% 
    mutate(model_to_test = list(model))
  
  # return the wide-version of the lookup table
  return(lookup_table_wide)
  
}



## Function K: create_percentchange_table

# this function predicts total conservative cone counts [±1SE] using
# key DBH, CrownDamage, and SiteName values of interest for a given
# spp-yr base model, calculates the percent change in cone counts across
# specific CrownDamage levels, and compiles predicted cone counts [±1SE]
# and percent change values into a table for publication.

# function inputs:
# model (a fitted GAM) = a spp-yr base model with which to predict total 
#     conservative cone counts.
# lookup_table (a data.frame) = a table containing key DBH, CrownDamage,
#     and SiteName values for the specified spp-yr base model. DBH
#     values represent the 10th, 50th, and 90th percentile-sized trees
#     in the relevant spp-yr dataset; CrownDamage values are scaled values 
#     of 5%, 50%, and 95% crown volume damage, as well as the threshold value 
#     (ie. the % crown volume damage at which the minimum, negative second 
#     derivative occurs) for each tree size; SiteName is the median 
#     fecundity site for the relevant spp-yr.
# threshold_table (a data.frame) = a table containing the threshold
#     values for the predicted cone count curves of the 10th, 50th, 
#     and 90th percentile-sized trees in the relevant spp-yr. 

# function output:
# percentchange_table (a data.frame) = a table displaying predicted
#     conservative cone counts [±1SE], and the estimated percent change
#     in predicted cone counts, at key CrownDamage values for the 10th, 
#     50th, and 90th percentile-sized trees in the spp-yr.

K_create_percentchange_table <- function(model, lookup_table,
                                         threshold_table) {
  
  # predict total conservative cone counts using the lookup_table and
  # the specified spp-yr model, bind the predictions to the lookup_table,
  # and add a column for the scaled threshold values. Then, unscale the
  # predictor variables, and backtransform the model fit and the predicted
  # upper and lower bounds of the confidence intervals (ie. fit +/- 1SE)
  pred_df = predict(object = model,
                    newdata = lookup_table,
                    se.fit = TRUE) %>% 
    as_tibble() %>% 
    bind_cols(lookup_table) %>%  
    left_join(., threshold_table, by = "DBH_cm_scaled") %>% 
    mutate(DBH_cm_unscaled = (DBH_cm_scaled * DBH_sd) + DBH_mean) %>% 
    mutate(CrownDamage_percentvolume_unscaled = 
             (CrownDamage_percentvolume_scaled * CrownDamage_sd) +
             CrownDamage_mean) %>% 
    mutate(threshold_val_unscaled = (threshold_val * CrownDamage_sd) +
             CrownDamage_mean) %>% 
    mutate(fit_backtransformed = exp(fit)) %>% 
    mutate(CI_upperlimit = exp(fit + se.fit)) %>% 
    mutate(CI_lowerlimit = exp(fit - se.fit))
  
  # pivot pred_df wider and calculate the percent change in predicted
  # cone counts at paired CrownDamage values of interest
  pred_df_wider = pred_df %>% 
    group_by(DBH_cm_unscaled) %>% 
    pivot_wider(id_cols = c("DBH_cm_unscaled", "threshold_val_unscaled"),
                names_from = CrownDamage_label,
                values_from = c(fit_backtransformed,
                                CI_lowerlimit,
                                CI_upperlimit)) %>% 
    mutate(percentchange_50v5 = ((fit_backtransformed_50 - 
           fit_backtransformed_5)/fit_backtransformed_5) * 100) %>% 
    mutate(percentchange_thresholdv5 = ((fit_backtransformed_threshold - 
           fit_backtransformed_5)/fit_backtransformed_5) * 100) %>% 
    mutate(percentchange_95v50 = ((fit_backtransformed_95 - 
           fit_backtransformed_50)/fit_backtransformed_50) * 100) %>% 
    mutate(percentchange_95vthreshold = ((fit_backtransformed_95 - 
           fit_backtransformed_threshold)/fit_backtransformed_threshold) 
           * 100) %>% 
    mutate(percentchange_95v5 = ((fit_backtransformed_95 - 
           fit_backtransformed_5)/fit_backtransformed_5) * 100)
  
  # round the model fits, the upper and lower limits of the fitted  
  # confidence intervals, the calculated percent change values, the DBH
  # values, and the threshold values; collapse the backtransformed model
  # fits and confidence intervals into individual columns (fit [±1SE]); 
  # label each DBH value according to its corresponding quantile; and 
  # rename and rearrange columns, as appropriate
  percentchange_table = pred_df_wider %>% 
    mutate(across(all_of(c(starts_with("fit"), starts_with("CI"))),
                  .fns = ~round(.x, digits = 0))) %>% 
    mutate(across(all_of(c(starts_with("percent"))),
                  .fns = ~round(.x, digits = 1))) %>% 
    mutate("Mean predicted cone count [±1SE] at 5% crown vol damage" = 
             paste(fit_backtransformed_5,
                   paste("[", paste(CI_lowerlimit_5, 
                                    CI_upperlimit_5, 
                                    sep = ", "), 
                         "]", sep = ""),
                   sep = " ")) %>% 
    mutate("Mean predicted cone count [±1SE] at 50% crown vol damage" = 
             paste(fit_backtransformed_50,
                   paste("[", paste(CI_lowerlimit_50, 
                                    CI_upperlimit_50, 
                                    sep = ", "), 
                         "]", sep = ""),
                   sep = " ")) %>% 
    mutate("Mean predicted cone count [±1SE] at threshold point" = 
             paste(fit_backtransformed_threshold,
                   paste("[", paste(CI_lowerlimit_threshold, 
                                    CI_upperlimit_threshold, 
                                    sep = ", "), 
                         "]", sep = ""),
                   sep = " ")) %>% 
    mutate("Mean predicted cone count [±1SE] at 95% crown vol damage" = 
             paste(fit_backtransformed_95,
                   paste("[", paste(CI_lowerlimit_95, 
                                    CI_upperlimit_95, 
                                    sep = ", "), 
                         "]", sep = ""),
                   sep = " ")) %>% 
    rename("Percent change in cone count (50% vs. 5% crown vol damage)" =
             percentchange_50v5,
           "Percent change in cone count (threshold point vs. 
           5% crown vol damage)" = percentchange_thresholdv5,
           "Percent change in cone count (95% vs. 50% crown vol damage)" = 
             percentchange_95v50,
           "Percent change in cone count (95% crown vol damage vs. 
           threshold point)" = percentchange_95vthreshold,
           "Percent change in cone count (95% vs. 5% crown vol damage)" = 
             percentchange_95v5) %>% 
    mutate(DBH_cm_unscaled = round(DBH_cm_unscaled, digits = 1)) %>% 
    mutate(DBH_cm_unscaled_chr = format(DBH_cm_unscaled, nsmall = 1)) %>% 
    ungroup() %>% 
    mutate("Tree size (cm DBH)" = case_when(
      DBH_cm_unscaled == min(DBH_cm_unscaled) ~ 
        paste("10th percentile (", 
              DBH_cm_unscaled_chr,
              ")", sep = ""),
      DBH_cm_unscaled == median(DBH_cm_unscaled) ~ 
        paste("50th percentile (",
              DBH_cm_unscaled_chr,
              ")", sep = ""),
      DBH_cm_unscaled == max(DBH_cm_unscaled) ~ 
        paste("90th percentile (",
              DBH_cm_unscaled_chr,
              ")", sep = ""),
      .default = NA)) %>% 
    mutate("Threshold point (% crown vol damage)" = 
             round(threshold_val_unscaled, digits = 1)) %>% 
    mutate(`Threshold point (% crown vol damage)` = format(
      `Threshold point (% crown vol damage)`, nsmall = 1)) %>% 
    select(c("Tree size (cm DBH)",
             "Threshold point (% crown vol damage)",
             starts_with("Mean"),
             starts_with("Percent")))
  
  # NOTE that the code above arranges mean predicted cone count columns
  # and percent change columns from left to right according to the 
  # level of crown volume damage with which they are associated (e.g., 
  # columns associated with 5% crown volume damage are situated furthest 
  # to the left, followed by those associated with 50% crown volume
  # damage, followed by those associated with threshold values, 
  # followed by those associated with 95% crown volume damage). For
  # some curves, threshold values occur BEFORE 50% crown volume damage
  # (ie. threshold_val < 50), and in such cases, it may make more sense 
  # to rearrange the columns of the percentchange_table to reflect this
  # (once the table is returned).
  
  # return the percentchange_table
  return(percentchange_table)
  
}



## Function L1: test_fecundity_differences

# this function sets up, runs, and summarizes custom pairwise comparisons
# of predicted cone counts at key levels of crown volume damage using 
# the function emmeans. Function L1 is used within Function L2 (described
# below) to run these custom contrasts on multiple rows of data within
# a data.frame (specifically, within the wide version of the relevant
# spp-yr lookup table, lookup_table_wide).

# function inputs:
# NOTE: all function inputs described below are column names from 
# the wide version of the spp-yr lookup table - ie. lookup_table_wide -
# generated via Function J.
# model_to_test (a list) = a list of attributes affiliated with a given
#     spp-yr base model, used to predict total conservative cone counts. 
# DBH_cm_scaled (a numeric value) = a scaled DBH value, typically 
#     associated with either the 10th, 50th, or 90th percentile-sized
#     tree in a given spp-yr dataset.
# CrownDamage_percentvolume_scaled_5 (a numeric value) = the scaled
#     value of 5% crown volume damage.
# CrownDamage_percentvolume_scaled_50 (a numeric value) = the scaled
#     value of 50% crown volume damage.
# CrownDamage_percentvolume_scaled_threshold (a numeric value) = the 
#     scaled value of the threshold (ie. the level of crown volume
#     damage at which the minimum, negative finite second derivative 
#     occurs along a curve).
# CrownDamage_percentvolume_scaled_95 (a numeric value) = the scaled
#     value of 95% crown volume damage.
# SiteName (a factor) = one of four field site names ("AUGUST", "CREEK",
#     "CAPLES", or "NORTH"); typically, the median cone production site
#     for the relevant spp-yr.

# function output: 
# contrasts_summary (a data.frame) = an emmeans-produced summary
#     table of contrasts results.

L1_test_fecundity_differences <- function(
    model_to_test, 
    DBH_cm_scaled, 
    CrownDamage_percentvolume_scaled_5, 
    CrownDamage_percentvolume_scaled_50,
    CrownDamage_percentvolume_scaled_threshold,
    CrownDamage_percentvolume_scaled_95,
    SiteName) {
  
  # require the R package "emmeans" to run this function; "emmeans"
  # contains the functions "emmeans" and "contrast," used below
  # to run custom pairwise comparisons of predicted cone counts
  # at key levels of crown volume damage
  require("emmeans")
  
  # combine crown damage levels of interest into a vector
  damage_levels = c(CrownDamage_percentvolume_scaled_5, 
                    CrownDamage_percentvolume_scaled_50,
                    CrownDamage_percentvolume_scaled_threshold,
                    CrownDamage_percentvolume_scaled_95)
  
  # if CrownDamage_percentvolume_scaled_threshold is "NA" - as is the
  # case when there is no minimum, negative finite second derivative
  # value for a given cone count curve (see Function G2 above for further
  # discussion) - replace this "NA" with an arbitrary number (e.g., 100). 
  # This allows emmeans to run below, calculating marginal means and 
  # contrasts for non-NA crown volume damage values, while also 
  # calculating marginal means and contrasts for this arbitrarily-assigned
  # threshold value. Mean differences and related metrics calculated for
  # this arbitrarily-assigned threshold value are replaced with "NA" in the 
  # contrasts_summary table below.
  if (is.na(CrownDamage_percentvolume_scaled_threshold)) {
    
    damage_levels[3] = 100
    
  }
  
  # find the estimated marginal means of predicted cone counts
  # at specified DBH, CrownDamage, and SiteName values
  means_grid = emmeans(model_to_test,
                       specs = ~ CrownDamage_percentvolume_scaled |
                         DBH_cm_scaled,
                       at = list(CrownDamage_percentvolume_scaled =
                                   damage_levels,
                                 DBH_cm_scaled = DBH_cm_scaled,
                                 SiteName = SiteName))
  
  # set up a custom contrasts list; for this to work, damage levels must
  # be in this exact order: damage_5, damage_50, damage_threshold, 
  # damage_95 in our damage_levels vector (ie. the marginal mean
  # affiliated with the threshold value must be listed third in 
  # the means_grid)
  custom_contrasts = list(
    "50% - 5%" = c(-1, 1, 0, 0),
    "threshold - 5%" = c(-1, 0, 1, 0),
    "95% - 50%" = c(0, -1, 0, 1),
    "95% - threshold" = c(0, 0, -1, 1),
    "95% - 5%" = c(-1, 0, 0, 1))
  
  # run the custom contrasts, comparing fecundity at specific crown
  # damage levels, using a bonferroni adjustment for multiple comparisons;
  # report the results on the response scale (ie. as ratios of total
  # cone counts, rather than differences in logged cone counts)
  contrasts_grid = contrast(means_grid, 
                            custom_contrasts, 
                            type = "response",
                            adjust = "bonferroni")
  
  # add confidence intervals to the contrasts summary table
  contrasts_summary = summary(contrasts_grid,
                              infer = TRUE)
  
  # if CrownDamage_percentvolume_scaled_threshold is "NA" (and was
  # therefore replaced with an arbitrary number (e.g., 100) above)...
  if (is.na(CrownDamage_percentvolume_scaled_threshold)) {
    
    # find the indices (ie. the locations) of any rows within the
    # contrasts_summary table that contain information related
    # to this threshold value
    NA_index = grep(pattern = "threshold",
                    x = contrasts_summary$contrast)
    
    # and replace all of the contrasts-related metrics contained 
    # within these rows with "NA"
    contrasts_summary[NA_index, 3:10] <- NA
    
  }
  
  # return the contrasts summary table, with confidence intervals
  # and NA's (where appropriate)
  return(contrasts_summary)
  
}



## Function L2: create_contrasts_summary_table

# this function applies Function L1 (test_fecundity_differences) to the
# wide-version of the spp-yr lookup table (lookup_table_wide, produced
# via Function J) to run custom pairwise comparisons of predicted
# cone counts at key crown volume damage levels for each specified
# tree size (ie. for the 10th, 50th, and 90th percentile-sized trees
# in the relevant spp-yr).

# function input:
# lookup_table_wide (a data.frame) = a wide version of the lookup table
#     produced using Function I (described above), with key DBH values
#     of interest represented in rows, key CrownDamage values of interest
#     represented in columns, and an additional column containing 
#     attributes of the relevant spp-yr base model.

# function output:
# contrasts_summary_table (a data.frame) = an emmeans-produced summary
#     table of contrasts results, grouped by tree size (ie. DBH value).

L2_create_contrasts_summary_table <- function(lookup_table_wide) {
  
  # apply Function L1 (test_fecundity_differences) to each row
  # of lookup_table_wide
  contrasts_summary_table = lookup_table_wide %>% 
    pmap_dfr(L1_test_fecundity_differences)
  
  # return the contrasts_summary_table
  return(contrasts_summary_table)
  
}



## Function L3: format_contrasts_summary_table

# this function formats the contrasts_summary_table produced by
# Function L2 (see above) for publication, pivoting the table
# wider (so that each row corresponds to a DBH value of interest,
# and each column holds important contrasts results), and reformatting,
# renaming, and rearranging columns, as appropriate.

# function inputs:
# contrasts_summary_table (a data.frame) = an emmeans-produced summary
#     table of contrasts results (ie. results from comparisons of predicted  
#     cone counts across key pairs of CrownDamage values), grouped by tree
#     size.
# threshold_table (a data.frame) = a table containing the threshold
#     values for the predicted cone count curves of the 10th, 50th, 
#     and 90th percentile-sized trees in the relevant spp-yr. 

# function output: 
# contrasts_summary_table_wide (a data.frame) = a wide version of the 
#     contrasts_summary_table produced by Function L2 (see above), with
#     estimated mean ratios of cone counts, 95% confidence intervals 
#     for those estimates, and p-values formatted for publication.

L3_format_contrasts_summary_table <- function(contrasts_summary_table,
                                              threshold_table) {
  
  # prepare threshold_table for merge by rounding DBH_cm_scaled values
  threshold_table_formerge = threshold_table %>% 
    mutate(DBH_cm_scaled = round(DBH_cm_scaled, digits = 6))
  
  # format contrasts_summary_table for publication by pivoting it
  # wider (so that each row corresponds to a DBH value of interest,
  # and each column holds important contrasts results); adding a
  # column for tree size-relevant threshold values; rounding 
  # estimated mean ratios, confidence intervals, and p-values;
  # collapsing mean ratios and associated confidence intervals
  # into individual columns (mean ratio [95% CI]); unscaling
  # DBH and threshold values; labeling DBH values according to their
  # corresponding quantile; assigning asterisks to p-values, 
  # according to estimated contrast significance at alpha = 0.05 (with 
  # bonferroni adjustment for each set of contrasts), and
  # renaming and rearranging columns, as appropriate
  contrasts_summary_table_wide = contrasts_summary_table %>% 
    select(-c("SE", "df", "null", "t.ratio")) %>% 
    mutate(DBH_cm_scaled = as.numeric(as.character(DBH_cm_scaled))) %>% 
    mutate(DBH_cm_scaled = round(DBH_cm_scaled, digits = 6)) %>% 
    left_join(., threshold_table_formerge, by = "DBH_cm_scaled") %>% 
    mutate(contrast = case_when(
      contrast == "50% / 5%" ~ "50div5",
      contrast == "threshold / 5%" ~ "thresholddiv5",
      contrast == "95% / 50%" ~ "95div50",
      contrast == "95% / threshold" ~ "95divthreshold",
      contrast == "95% / 5%" ~ "95div5",
      .default = NA)) %>% 
    pivot_wider(id_cols = c("DBH_cm_scaled", "threshold_val"),
                names_from = contrast,
                values_from = c(ratio,
                                lower.CL,
                                upper.CL,
                                p.value)) %>% 
    mutate(across(all_of(c(starts_with("ratio"), starts_with("lower"),
                           starts_with("upper"))),
                  .fns = ~format(round(.x, digits = 2), nsmall = 2))) %>% 
    mutate("Mean ratio [95% CI] in predicted cone counts 
           (50%/5% crown vol damage)" = paste(ratio_50div5,
                                              paste("[", 
                                                    paste(lower.CL_50div5,
                                                          upper.CL_50div5,
                                                          sep = ", "),
                                                    "]", sep = ""),
                                              sep = " ")) %>% 
    mutate("Mean ratio [95% CI] in predicted cone counts 
           (threshold/5% crown vol damage)" = paste(ratio_thresholddiv5,
                                                    paste("[", paste(
                                                      lower.CL_thresholddiv5,
                                                      upper.CL_thresholddiv5,
                                                      sep = ", "),
                                                      "]", sep = ""),
                                                    sep = " ")) %>% 
    mutate("Mean ratio [95% CI] in predicted cone counts 
           (95%/50% crown vol damage)" = paste(ratio_95div50,
                                               paste("[", 
                                                     paste(lower.CL_95div50,
                                                           upper.CL_95div50,
                                                           sep = ", "),
                                                     "]", sep = ""),
                                               sep = " ")) %>% 
    mutate("Mean ratio [95% CI] in predicted cone counts
           (95% crown vol damage/threshold)" = paste(ratio_95divthreshold,
                                                     paste("[", paste(
                                                       lower.CL_95divthreshold,
                                                       upper.CL_95divthreshold,
                                                       sep = ", "),
                                                       "]", sep = ""),
                                                     sep = " ")) %>% 
    mutate("Mean ratio [95% CI] in predicted cone counts
           (95%/5% crown vol damage)" = paste(ratio_95div5,
                                              paste("[", 
                                                    paste(lower.CL_95div5,
                                                          upper.CL_95div5,
                                                          sep = ", "),
                                                    "]", sep = ""),
                                              sep = " ")) %>% 
    mutate(across(all_of(c(starts_with("p.value"))),
                  .fns = ~case_when(.x < 0.001 ~ "***",
                                    .x >= 0.001 & .x < 0.01 ~ "**",
                                    .x >= 0.01 & .x < 0.05 ~ "*",
                                    .x >= 0.05 & .x < 0.1 ~ "'",
                                    .default = ""),
                  .names = "{.col}_stars")) %>% 
    mutate(across(c("p.value_50div5",
                    "p.value_thresholddiv5", 
                    "p.value_95div50", 
                    "p.value_95divthreshold",
                    "p.value_95div5"),
                  .fns = ~case_when(.x < 0.001 ~ "<0.001",
                                    .default = format(round(.x, 
                                                            digits = 3), 
                                                      nsmall = 3)))) %>% 
    mutate(p.value_50div5 = paste(p.value_50div5,
                                  p.value_50div5_stars,
                                  sep = "")) %>% 
    mutate(p.value_thresholddiv5 = paste(p.value_thresholddiv5,
                                         p.value_thresholddiv5_stars,
                                         sep = "")) %>% 
    mutate(p.value_95div50 = paste(p.value_95div50,
                                   p.value_95div50_stars,
                                   sep = "")) %>% 
    mutate(p.value_95divthreshold = paste(p.value_95divthreshold,
                                          p.value_95divthreshold_stars,
                                          sep = "")) %>% 
    mutate(p.value_95div5 = paste(p.value_95div5,
                                  p.value_95div5_stars,
                                  sep = "")) %>% 
    rename("p-value for contrast (50%/5% crown vol damage)" = 
             p.value_50div5,
           "p-value for contrast (threshold/5% crown vol damage)" = 
             p.value_thresholddiv5,
           "p-value for contrast (95%/50% crown vol damage)" = 
             p.value_95div50,
           "p-value for contrast (95% crown vol damage/threshold)" = 
             p.value_95divthreshold,
           "p-value for contrast (95%/5% crown vol damage)" = 
             p.value_95div5) %>% 
    mutate(DBH_cm_unscaled = (DBH_cm_scaled * DBH_sd) + DBH_mean) %>% 
    mutate(threshold_val_unscaled = (threshold_val * CrownDamage_sd) + 
             CrownDamage_mean) %>% 
    mutate(across(c("DBH_cm_unscaled", "threshold_val_unscaled"),
                  .fns = ~round(.x, digits = 1))) %>% 
    mutate(across(c("DBH_cm_unscaled", "threshold_val_unscaled"),
                  .fns = ~format(.x, nsmall = 1))) %>% 
    mutate("Tree size (cm DBH)" = case_when(
      DBH_cm_scaled == min(DBH_cm_scaled) ~ paste("10th percentile (",
                                                  DBH_cm_unscaled,
                                                  ")", sep = ""),
      DBH_cm_scaled == median(DBH_cm_scaled) ~ paste("50th percentile (",
                                                     DBH_cm_unscaled,
                                                     ")", sep = ""),
      DBH_cm_scaled == max(DBH_cm_scaled) ~ paste("90th percentile (",
                                                  DBH_cm_unscaled,
                                                  ")", sep = ""),
      .default = NA)) %>% 
    rename("Threshold point (% crown vol damage)" = threshold_val_unscaled) %>% 
    select(c(starts_with("Tree size"),
             starts_with("Threshold point"),
             starts_with("Mean ratio"),
             starts_with("p-value for"))) %>% 
    select(c(1:3, 8, 4, 9, 5, 10, 6, 11, 7, 12))
  
  # NOTE that the columns in this table are ordered as if the threshold 
  # values for each DBH class occur at points >50% crown volume damage;
  # thus, these columns may need to be reordered if threshold 
  # values are <50% crown volume damage, for ease of reading and
  # interpretation from left to right, once the table is returned. 
  
  # return the wide version of the contrasts_summary_table
  return(contrasts_summary_table_wide)
  
}




