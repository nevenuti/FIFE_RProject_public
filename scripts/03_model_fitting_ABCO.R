

# Title: Modeling the effects of fire injury on ABCO (white fir) 
# fecundity one- and two-year(s) post-fire

# Author: Nina Venuti

# Script inputs: a .csv file of fire injury-fecundity project data, 
# cleaned and prepared for modeling

# Script outputs: a series of fitted generalized additive models (GAMs)
# of ABCO cone counts one- and two-year(s) after fire injury


## loading relevant packages...
library("tidyverse")
library("corrplot")
library("mgcv")

## loading the final fire injury-fecundity model data...
model_data <- read_csv("data/FIFE_model_data_FINAL.csv")


## taking a few more steps to prepare the data for modeling by
## converting select character columns into factor columns, scaling 
## continuous predictor variables, and subsetting the data into species 
## and sample year-specific dataframes...

# first, converting select character columns into factor columns...
chr_predictors <- c("SiteName", "SppCode", "CanopyPosition", 
                    "TreeStatus_yr1", "TreeStatus_yr2")

model_data[chr_predictors] <- lapply(model_data[chr_predictors], as.factor)

# second, scaling all of the continuous predictor variables (using dataset-
# wide variable means and standard deviations), to facilitate model fitting
# and effects comparison...
continuous_predictors <- c("DBH_cm",
                           "BoleCharMeanHeight_m",
                           "CrownDamage_percentvolume",
                           "MistletoeClumpsOrBranches_count_yr1",
                           "MistletoeClumpsOrBranches_count_yr2")

model_data_scaled <- model_data %>% 
  mutate(across(all_of(continuous_predictors),
                .fns = ~c(scale(.x)),
                .names = "{.col}_scaled"))
# note that enclosing the function scale() inside a c() above converts 
# its output from a matrix to a vector and allows the output to be returned
# as columns within a dataframe (which is desirable for our purposes)

# finally, splitting the model data into smaller subsets that represent 
# live ABCO trees in the first and second year of sampling; modeling will
# be conducted separately for each sample year, using data collected on
# live trees only (ie. excluding trees that were dead, fallen, cut down, 
# broken-topped, or top-killed in either yr1 or yr2). While extracting
# and saving the yr1 dataset, remove any additional trees that lack data 
# for key predictors using the filter CanopyPosition != "blank"; this 
# excludes CAPLES trees that were windthrown between years 1 and 2, for which
# data on CanopyPosition and Mistletoe_count could not be collected during  
# yr2 and applied to yr1 (these two metrics were not integrated into the 
# field sampling protocol until 2021, ie. sample yr2 for CAPLES trees). 
# Trees with missing data are already removed from the yr2 dataset using 
# the filter TreeStatus_yr2 == "A".
ABCO_data_yr1 <- model_data_scaled %>% 
  filter(SppCode == "ABCO") %>% 
  filter(TreeStatus_yr1 == "A") %>% 
  filter(CanopyPosition != "blank")

ABCO_data_yr2 <- model_data_scaled %>% 
  filter(SppCode == "ABCO") %>% 
  filter(TreeStatus_yr2 == "A")


## calculating and inspecting correlation coefficients for pairs of predictor
## variables, to make sure predictors are not highly correlated with one
## another and may be included in the same models... 

# starting with the yr1 data...
# selecting key predictor variables for which to calculate correlation 
# coefficients...
ABCO_correlations_yr1 <- ABCO_data_yr1 %>% 
  select(c(6:9, 12))

# converting CanopyPosition (a factor variable) into a numeric variable
# temporarily, by assigning a score of either 2 (for "DO") or 
# 1 (for "CO") to each tree, to enable correlation coefficient
# calculations; afterwards, removing the original CanopyPosition column...
ABCO_correlations_yr1 <- ABCO_correlations_yr1 %>% 
  mutate(CanopyPositionScore = case_when(CanopyPosition == "DO" ~ 2,
                                         CanopyPosition == "CO" ~ 1)) %>% 
  select(-CanopyPosition)

# saving the yr1 correlations data.frame as a matrix...
ABCO_correlations_yr1_matrix <- as.matrix(ABCO_correlations_yr1)

# calculating correlation coefficients for all of the yr1 variable pairs...
ABCO_correlations_yr1_pairs <- cor(ABCO_correlations_yr1_matrix,
                                   use = "pairwise.complete.obs")

# making a correlation plot, to enable easy inspection of correlation
# coefficients...
corrplot(ABCO_correlations_yr1_pairs,
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black")
# --> this correlation plot demonstrates that yr1 predictor variables
# are only weakly or moderately correlated with one another (ie. all
# correlation coefficients are <0.6), and thus it should be ok to 
# include them in the same models. CrownDamage_percentvolume and
# BoleCharMeanHeight_m are moderately correlated with each other (coeff = 
# 0.44), as are DBH_cm and MistletoeClumpsOrBranches_count_yr1 (coeff =
# 0.42), and DBH_cm and CanopyPositionScore (coeff = 0.31). 

# repeating this process with the yr2 data...
# selecting key predictor variables for which to calculate correlation 
# coefficients...
ABCO_correlations_yr2 <- ABCO_data_yr2 %>% 
  select(c(6:9, 13))

# converting CanopyPosition (a factor variable) into a numeric variable
# temporarily, by assigning a score of either 2 (for "DO") or 
# 1 (for "CO") to each tree, to enable correlation coefficient
# calculations; afterwards, removing the original CanopyPosition column...
ABCO_correlations_yr2 <- ABCO_correlations_yr2 %>% 
  mutate(CanopyPositionScore = case_when(CanopyPosition == "DO" ~ 2,
                                         CanopyPosition == "CO" ~ 1)) %>% 
  select(-CanopyPosition)

# saving the yr2 correlations data.frame as a matrix...
ABCO_correlations_yr2_matrix <- as.matrix(ABCO_correlations_yr2)

# calculating correlation coefficients for all of the yr2 variable pairs...
ABCO_correlations_yr2_pairs <- cor(ABCO_correlations_yr2_matrix,
                                   use = "pairwise.complete.obs")

# making a correlation plot, to enable easy inspection of correlation
# coefficients...
corrplot(ABCO_correlations_yr2_pairs,
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black")
# --> this correlation plot demonstrates that yr2 predictor variables
# are only weakly or moderately correlated with one another (ie. all
# correlation coefficients are <0.6), and thus it should be ok to 
# include them in the same models. DBH_cm and BoleCharMeanHeight_m 
# are moderately correlated with each other (coeff = 0.37), as are
# DBH_cm and MistletoeClumpsOrBranches_count_yr2 (coeff = 0.44),
# DBH_cm and CanopyPositionScore (coeff = 0.34), and BoleCharMeanHeight_m
# and CrownDamage_percentvolume (coeff = 0.31).



## writing and evaluating models of ABCO cone counts...

## first, some notes about modeling protocols...

# for each focal species and sample year in this study, modeling
# begins with a "base model" which fits total conservative cone 
# counts (TotalConesConservative_count) as a function of tree size
# (DBH_cm), percent crown volume damage (CrownDamage_percentvolume),
# the interaction between these two terms (DBH:CrownDamage), and study
# site (SiteName), with the first three variables (DBH, CrownDamage, DBH:
# CrownDamage) included as smooth terms, and the last variable (SiteName)
# included as a parametric term (ie. TotalConesConservative_count ~ 
# s(DBH) + s(CrownDamage) + ti(DBH:CrownDamage) + SiteName). From there,
# each "auxiliary" variable of interest - BoleCharMeanHeight_m,
# CanopyPosition, and MistletoeClumpsOrBranches_count (_yr1 or
# _yr2, as appropriate) - is added to the base model, one by one,
# to test its effects on total cone counts. BoleChar and Mistletoe_
# count are added to the base model as smooth terms, while CanopyPosition
# is added as a parametric term. Each smooth term included in these
# models - except the interaction term, DBH:CrownDamage - is estimated
# using a maximum of five basis functions (ie. the parameter "k" is set to
# 5), in order to constrain the "wiggliness" of each smooth. The 
# interaction term, DBH:CrownDamage, is estimated using a maximum of 
# three basis functions (ie. k = 3), to reduce the estimated pairwise
# concurvity values between it and the other terms in each model.

# for every model written below, we examine the diagnostic plots
# (to evaluate if there are any problematic patterns in the model 
# residuals), the model summary table, the partial effects plots, and the
# estimated pairwise concurvity values for the model terms, systematically.
# We consider any estimated pairwise concurvity value >0.6 to be "high"
# and to trigger a reduction in the number of maximum basis functions
# allowed for the smooth term most recently added to the model (ie. k is 
# reduced first to 4, and, if necessary, to 3). We examine "estimated" 
# pairwise concurvity values - rather than "worst" or "observed" pairwise 
# concurvity values - to determine the degree of concurvity between model 
# terms because the "estimated" value is neither overly-pessimistic (as
# is the "worst" value) nor overly-optimistic (as is the "observed" value), 
# according to S. Wood, the author of the mgcv package.


## starting with yr1 ABCO models...

## ABCO_m1 (yr1 base model)...
ABCO_m1 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr1)

par(mfrow = c(2,2))
gam.check(ABCO_m1) 
# ABCO_m1 diagnostic plots look great; the histogram of residuals is 
# basically bell-shaped, and there are no strong patterns in the residuals 
# vs. linear predictor plot.

summary(ABCO_m1)
plot(ABCO_m1, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH is estimated to have a positive, linear relationship with 
# cone count (as DBH goes up, cone count goes up); CrownDamage is
# estimated to have a largely negative, nonlinear relationship
# with cone count (as CrownDamage goes up, cone count decreases
# nonlinearly - somewhat slowly at first, and then more quickly); 
# DBH:CrownDamage is estimated to have a linear effect on cone count. 
# DBH and CrownDamage are "significant" terms (at alpha = 0.05), 
# DBH:CrownDamage is not. "CREEK" is the median fecundity site for 
# yr1 ABCO.

ABCO_m1_concurvity <- concurvity(ABCO_m1, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## ABCO_m2 (yr1 base model + BoleCharMeanHeight_m)...
ABCO_m2 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(BoleCharMeanHeight_m_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr1)

gam.check(ABCO_m2) 
# ABCO_m2 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of ABCO_m1.

summary(ABCO_m2)
plot(ABCO_m2, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for ABCO_m1. BoleChar is estimated to have
# a linear and very slightly positive relationship with cone count
# (as BoleChar goes up, cone count goes up slightly), though the partial 
# effects plot displays very large confidence intervals around this 
# estimated effect, indicating high uncertainty. DBH and CrownDamage are 
# "significant" terms (at alpha = 0.05), DBH:CrownDamage and BoleChar 
# are not. "CREEK" is the median fecundity site.

ABCO_m2_concurvity <- concurvity(ABCO_m2, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## ABCO_m3 (yr1 base model + CanopyPosition)...
ABCO_m3 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 CanopyPosition +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr1)

gam.check(ABCO_m3) 
# ABCO_m3 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of ABCO_m1.

summary(ABCO_m3)
plot(ABCO_m3, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for ABCO_m1. DBH and CrownDamage are "significant"
# terms (at alpha = 0.05), DBH:CrownDamage is not. "DO" (dominant)
# trees are estimated to produce slightly more cones than "CO"
# (co-dominant) trees, but this effect is not significant (at alpha = 
# 0.05). "CREEK" is the median fecundity site.

ABCO_m3_concurvity <- concurvity(ABCO_m3, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## ABCO_m4 (yr1 base model + MistletoeClumpsOrBranches_count_yr1)...
ABCO_m4 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(MistletoeClumpsOrBranches_count_yr1_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr1)

gam.check(ABCO_m4) 
# ABCO_m4 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of ABCO_m1.

summary(ABCO_m4)
plot(ABCO_m4, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for ABCO_m1. Mistletoe_count is estimated to have
# a linear, negative relationship with cone count (as the level
# of mistletoe infestation goes up, cone count goes down), though the
# partial effects plot displays very large confidence intervals around
# this estimated effect, indicating high uncertainty. DBH and CrownDamage
# are "significant" terms (at alpha = 0.05), DBH:CrownDamage and 
# Mistletoe_count are not. "CREEK" is the median fecundity site.

ABCO_m4_concurvity <- concurvity(ABCO_m4, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).



## moving on to yr2 ABCO models...

## ABCO_m5 (yr2 base model)...
ABCO_m5 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr2)

gam.check(ABCO_m5)
# ABCO_m5 diagnostic plots look pretty good; the histogram of residuals
# is relatively bell-shaped, and there are no strong patterns in the
# residuals vs. linear predictor plot.

summary(ABCO_m5)
plot(ABCO_m5, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH is estimated to have a slightly positive, linear relationship 
# with cone count (as DBH goes up, cone count goes up slightly); CrownDamage
# is estimated to have a nonlinear, and overall negative relationship with 
# cone count (as CrownDamage goes up, cone count increases slightly/remains 
# stable, and then begins to decrease at low to moderate CrownDamage 
# values); DBH:CrownDamage is estimated to have a nonlinear effect on cone
# count. CrownDamage is a "significant" term (at alpha = 0.05), DBH 
# and DBH:CrownDamage are not. "CREEK" is the median fecundity site for
# yr2 ABCO, though it appears to be indistinguishable from "AUGUST" and 
# "NORTH".

ABCO_m5_concurvity <- concurvity(ABCO_m5, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## ABCO_m6 (yr2 base model + BoleCharMeanHeight_m)...
ABCO_m6 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(BoleCharMeanHeight_m_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr2)

gam.check(ABCO_m6) 
# ABCO_m6 diagnostic plots look pretty good; same basic patterns as those
# displayed in diagnostic plots of ABCO_m5.

summary(ABCO_m6)
plot(ABCO_m6, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for ABCO_m5. BoleChar is estimated to have a 
# negative, linear relationship with cone count (as BoleChar goes up,
# cone count goes down), though the partial effects plot displays
# large confidence intervals around this estimated effect, indicating
# high uncertainty. CrownDamage is "significant" (at alpha = 0.05),
# DBH, DBH:CrownDamage, and BoleChar are not. "CREEK" is the median
# fecundity site, though it appears to be indistinguishable from 
# "AUGUST" and "NORTH".

ABCO_m6_concurvity <- concurvity(ABCO_m6, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## ABCO_m7 (yr2 base model + CanopyPosition)
ABCO_m7 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 CanopyPosition +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr2)

gam.check(ABCO_m7)
# ABCO_m7 diagnostic plots look pretty good; same basic patterns as those
# displayed in diagnostic plots of ABCO_m5.

summary(ABCO_m7)
plot(ABCO_m7, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for ABCO_m5. CrownDamage is "significant" 
# (at alpha = 0.05), DBH and DBH:CrownDamage are not. "DO" (dominant)
# trees are estimated to produce slightly more cones than "CO" (co-
# dominant) trees, but this effect is not significant (at alpha =
# 0.05). "NORTH" is now the median fecundity site, though it appears 
# to be indistinguishable from "AUGUST" and "CREEK".

ABCO_m7_concurvity <- concurvity(ABCO_m7, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## ABCO_m8 (yr2 base model + MistletoeClumpsOrBranches_count_yr2)
ABCO_m8 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(MistletoeClumpsOrBranches_count_yr2_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr2)

gam.check(ABCO_m8) 
# ABCO_m8 diagnostic plots look pretty good; same basic patterns as those
# displayed in diagnostic plots of ABCO_m5.

summary(ABCO_m8)
plot(ABCO_m8, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for ABCO_m5. Mistletoe_count is estimated to
# have a nonlinear relationship with cone count (as the level of
# mistletoe infestation goes up, cone count increases, reaches a 
# maximum, and then decreases), though the partial effects plot
# displays very large confidence intervals around this estimated 
# effect, indicating high uncertainty. CrownDamage is "significant" 
# (at alpha = 0.05), DBH, DBH:CrownDamage, and Mistletoe_count are not.
# "NORTH" is now the median fecundity site, though it appears to be
# indistinguishable from "AUGUST" and "CREEK".

ABCO_m8_concurvity <- concurvity(ABCO_m8, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).

## --> none of the auxiliary variables were found to be significant
## (at alpha = 0.05) in any of the yr1 or yr2 ABCO models. Based on
## this lack of auxiliary variable significance in the yr1 models,
## any trees with missing auxiliary variable data, previously removed
## from the yr1 dataset (ie. CAPLES trees that were windthrown between
## sample years 1 and 2), can be added back into the dataset, and the 
## yr1 base model can be refit using this more complete, yr1 dataset
## (see below).


## refitting the yr1 base model, using the full yr1 dataset...

# first, saving a new yr1 dataset, inclusive of all previously-removed
# trees...
ABCO_data_yr1_full <- model_data_scaled %>% 
  filter(SppCode == "ABCO") %>% 
  filter(TreeStatus_yr1 == "A")

# then, refitting the yr1 base model using these data...
ABCO_myr1 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 SiteName,
               family = nb,
               method = "REML",
               data = ABCO_data_yr1_full)

# quickly checking the refit, yr1 base model's diagnostic plots,
# summary table, partial effects plots, and estimated pairwise 
# concurvity values...
gam.check(ABCO_myr1)
summary(ABCO_myr1)
plot(ABCO_myr1, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
ABCO_myr1_concurvity <- concurvity(ABCO_myr1, full = FALSE)$estimate
# all diagnostic plots look good; no issues with concurvity; estimated
# effects are stable/consistent across ABCO_m1 and ABCO_myr1.


## finally, renaming the yr2 base model following the same convention
## used above (ie. ABCO_myrX) for easier reference during model
## summary table creation (see script #5)...
ABCO_myr2 <- ABCO_m5



