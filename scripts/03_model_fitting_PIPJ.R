

# Title: Modeling the effects of fire injury on PIPJ (ponderosa and Jeffrey
# pine - ie. yellow pine) fecundity one- and two-year(s) post-fire

# Author: Nina Venuti

# Script inputs: a .csv file of fire injury-fecundity project data, 
# cleaned and prepared for modeling

# Script outputs: a series of fitted generalized additive models (GAMs)
# of PIPJ cone counts one- and two-year(s) after fire injury


## loading relevant packages...
library("tidyverse")
library("corrplot")
library("mgcv")

## loading the final fire injury-fecundity model data...
model_data <- read_csv("data/processed/FIFE_model_data_FINAL.csv")


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
# live PIPJ (yellow pine) trees in the first and second year of sampling; 
# modeling will be conducted for all yellow pines (PIPO and PIJE), together,
# but will be conducted separately for each sample year, using data collected
# on live trees only (ie. excluding trees that were dead, fallen, cut down, 
# broken-topped, or top-killed in either yr1 or yr2). While extracting
# and saving the yr1 dataset, remove any additional trees that lack data 
# for key predictors using the filter CanopyPosition != "blank"; this 
# excludes CAPLES trees that were windthrown between years 1 and 2, for which
# data on CanopyPosition and Mistletoe_count could not be collected during  
# yr2 and applied to yr1 (these two metrics were not integrated into the 
# field sampling protocol until 2021, ie. sample yr2 for CAPLES trees). 
# Trees with missing data are already removed from the yr2 dataset using 
# the filter TreeStatus_yr2 == "A".
PIPJ_data_yr1 <- model_data_scaled %>% 
  filter(SppCode %in% c("PIJE", "PIPO")) %>% 
  filter(TreeStatus_yr1 == "A") %>% 
  filter(CanopyPosition != "blank")

PIPJ_data_yr2 <- model_data_scaled %>% 
  filter(SppCode %in% c("PIJE", "PIPO")) %>% 
  filter(TreeStatus_yr2 == "A")


## calculating and inspecting correlation coefficients for pairs of predictor
## variables, to make sure predictors are not highly correlated with one
## another and may be included in the same models... 

# starting with the yr1 data...
# selecting key predictor variables for which to calculate correlation 
# coefficients...
PIPJ_correlations_yr1 <- PIPJ_data_yr1 %>% 
  select(c(6:9, 12))

# converting CanopyPosition (a factor variable) into a numeric variable
# temporarily, by assigning a score of either 2 (for "DO") or 
# 1 (for "CO") to each tree, to enable correlation coefficient
# calculations; afterwards, removing the original CanopyPosition column...
PIPJ_correlations_yr1 <- PIPJ_correlations_yr1 %>% 
  mutate(CanopyPositionScore = case_when(CanopyPosition == "DO" ~ 2,
                                         CanopyPosition == "CO" ~ 1)) %>% 
  select(-CanopyPosition)

# saving the yr1 correlations data.frame as a matrix...
PIPJ_correlations_yr1_matrix <- as.matrix(PIPJ_correlations_yr1)

# calculating correlation coefficients for all of the yr1 variable pairs...
PIPJ_correlations_yr1_pairs <- cor(PIPJ_correlations_yr1_matrix,
                                   use = "pairwise.complete.obs")

# making a correlation plot, to enable easy inspection of correlation
# coefficients...
corrplot(PIPJ_correlations_yr1_pairs,
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black")
# --> this correlation plot demonstrates that yr1 predictor variables
# are only weakly or moderately correlated with one another (ie. all
# correlation coefficients are <0.6), and thus it should be ok to 
# include them in the same models. CrownDamage_percentvolume and
# BoleCharMeanHeight_m are moderately correlated with each other (coeff = 
# 0.50), as are DBH_cm and BoleCharMeanHeight_m (coeff = 0.31). 

# repeating this process with the yr2 data...
# selecting key predictor variables for which to calculate correlation 
# coefficients...
PIPJ_correlations_yr2 <- PIPJ_data_yr2 %>% 
  select(c(6:9, 13))

# converting CanopyPosition (a factor variable) into a numeric variable
# temporarily, by assigning a score of either 2 (for "DO") or 
# 1 (for "CO") to each tree, to enable correlation coefficient
# calculations; afterwards, removing the original CanopyPosition column...
PIPJ_correlations_yr2 <- PIPJ_correlations_yr2 %>% 
  mutate(CanopyPositionScore = case_when(CanopyPosition == "DO" ~ 2,
                                         CanopyPosition == "CO" ~ 1)) %>% 
  select(-CanopyPosition)

# saving the yr2 correlations data.frame as a matrix...
PIPJ_correlations_yr2_matrix <- as.matrix(PIPJ_correlations_yr2)

# calculating correlation coefficients for all of the yr2 variable pairs...
PIPJ_correlations_yr2_pairs <- cor(PIPJ_correlations_yr2_matrix,
                                   use = "pairwise.complete.obs")

# making a correlation plot, to enable easy inspection of correlation
# coefficients...
corrplot(PIPJ_correlations_yr2_pairs,
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black")
# --> this correlation plot demonstrates that yr2 predictor variables
# are only weakly or moderately correlated with one another (ie. all
# correlation coefficients are <0.6), and thus it should be ok to 
# include them in the same models. CrownDamage_percentvolume and
# BoleCharMeanHeight_m are moderately correlated with each other (coeff = 
# 0.40), as are DBH_cm and BoleCharMeanHeight_m (coeff = 0.35). 



## writing and evaluating models of PIPJ cone counts...

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


## starting with yr1 PIPJ models...

## PIPJ_m1 (yr1 base model)...
PIPJ_m1 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr1)

par(mfrow = c(2,2))
gam.check(PIPJ_m1)
# PIPJ_m1 diagnostic plots look great; the histogram of residuals is 
# basically bell-shaped, and there are no strong patterns in the residuals 
# vs. linear predictor plot.

summary(PIPJ_m1)
plot(PIPJ_m1, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH is estimated to have a nonlinear, and mostly positive relationship 
# with cone count (as DBH increases, cone count increases, and then
# begins to plateau at high DBH values); CrownDamage is estimated
# to have a nonlinear, and overall negative relationship with cone
# count (as CrownDamage increases, cone count remains stable, and then
# begins to decrease at high CrownDamage values); DBH:CrownDamage
# is estimated to have a nonlinear effect on cone count. DBH and 
# CrownDamage are "significant" terms (at alpha = 0.05), DBH:CrownDamage
# is not. "CREEK" is the median fecundity site for yr1 PIPJ.

PIPJ_m1_concurvity <- concurvity(PIPJ_m1, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PIPJ_m2 (yr1 base model + BoleCharMeanHeight_m)...
PIPJ_m2 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(BoleCharMeanHeight_m_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr1)

gam.check(PIPJ_m2) 
# PIPJ_m2 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of PIPJ_m1.

summary(PIPJ_m2)
plot(PIPJ_m2, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PIPJ_m1. BoleChar is estimated to have
# a linear and very slightly negative relationship with cone count
# (as BoleChar increases, cone count decreases slightly), though the 
# estimated effect appears to be very close to a horizontal line.
# DBH and CrownDamage are "significant" terms (at alpha = 0.05), 
# DBH:CrownDamage and BoleChar are not. "CREEK" is the median fecundity 
# site.

PIPJ_m2_concurvity <- concurvity(PIPJ_m2, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PIPJ_m3 (yr1 base model + CanopyPosition)...
PIPJ_m3 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 CanopyPosition +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr1)

gam.check(PIPJ_m3) 
# PIPJ_m3 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of PIPJ_m1.

summary(PIPJ_m3)
plot(PIPJ_m3, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PIPJ_m1. DBH and CrownDamage are "significant"
# terms (at alpha = 0.05), DBH:CrownDamage is not. "DO" (dominant)
# trees are estimated to produce slightly more cones than "CO"
# (co-dominant) trees, but this effect is not significant (at alpha = 
# 0.05). "CREEK" is the median fecundity site.

PIPJ_m3_concurvity <- concurvity(PIPJ_m3, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PIPJ_m4 (yr1 base model + MistletoeClumpsOrBranches_count_yr1)...
PIPJ_m4 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(MistletoeClumpsOrBranches_count_yr1_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr1)

gam.check(PIPJ_m4)
# PIPJ_m4 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of PIPJ_m1.

summary(PIPJ_m4)
plot(PIPJ_m4, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PIPJ_m1. Mistletoe_count is estimated to have
# a nonlinear, almost parabolic relationship with cone count (as the level
# of mistletoe infestation increases, cone count decreases, hits a minimum,
# and then increases again), though the partial effects plot displays very
# large confidence intervals around this estimated effect, indicating high
# uncertainty. DBH and CrownDamage are "significant" terms (at alpha = 
# 0.05), DBH:CrownDamage and Mistletoe_count are not. "CREEK" is the
# median fecundity site.

PIPJ_m4_concurvity <- concurvity(PIPJ_m4, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).



## moving on to yr2 PIPJ models...

## PIPJ_m5 (yr2 base model)...
PIPJ_m5 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr2)

gam.check(PIPJ_m5)
# PIPJ_m5 diagnostic plots look decent; the histogram of residuals 
# is relatively bell-shaped (though a bit right-skewed), and there are 
# no strong patterns in the residuals vs. linear predictor plot 
# (though there may be a slight negative trend in mean residual values
# as fitted values increase).

summary(PIPJ_m5)
plot(PIPJ_m5, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH is estimated to have a positive, linear relationship with 
# cone count (as DBH increases, cone count increases). CrownDamage is
# estimated to have a nonlinear, and overall negative relationship with
# cone count (as CrownDamage increases, cone count remains steady/may
# even increase slightly, and then begins to decline at moderate 
# CrownDamage values). DBH:CrownDamage is estimated to have a nearly 
# linear effect on cone count. DBH and CrownDamage are "significant" 
# terms (at alpha = 0.05), DBH:CrownDamage is not. "NORTH" is the median
# fecundity site for yr2 PIPJ.

PIPJ_m5_concurvity <- concurvity(PIPJ_m5, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PIPJ_m6 (yr2 base model + BoleCharMeanHeight_m)...
PIPJ_m6 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(BoleCharMeanHeight_m_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr2)

gam.check(PIPJ_m6)
# PIPJ_m6 diagnostic plots look decent; same basic patterns as those
# displayed in diagnostic plots of PIPJ_m5.

summary(PIPJ_m6)
plot(PIPJ_m6, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PIPJ_m5. BoleChar is estimated to have a 
# nonlinear, nearly quadratic relationship with cone count (as
# BoleChar increases, cone count decreases slightly, reaches a 
# minimum, and then increases slightly again), though the 
# estimated effect appears to be very close to a horizontal line
# in the partial effects plot. DBH and CrownDamage are "significant"
# (at alpha = 0.05), DBH:CrownDamage is marginally so (p = 0.06),
# and BoleChar is not. "NORTH" is the median fecundity site.

PIPJ_m6_concurvity <- concurvity(PIPJ_m6, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PIPJ_m7 (yr2 base model + CanopyPosition)
PIPJ_m7 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 CanopyPosition +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr2)

gam.check(PIPJ_m7)
# PIPJ_m7 diagnostic plots look decent; same basic patterns as those
# displayed in diagnostic plots of PIPJ_m5.

summary(PIPJ_m7)
plot(PIPJ_m7, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PIPJ_m5. DBH and CrownDamage are "significant" 
# (at alpha = 0.05), DBH:CrownDamage is not. "DO" (dominant)
# trees are estimated to produce slightly fewer cones than "CO" (co-
# dominant) trees, and this effect is marginally significant at alpha =
# 0.05 (p = 0.08). "NORTH" is the median fecundity site.

PIPJ_m7_concurvity <- concurvity(PIPJ_m7, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PIPJ_m8 (yr2 base model + MistletoeClumpsOrBranches_count_yr2)
PIPJ_m8 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, CrownDamage_percentvolume_scaled, k = 3) +
                 s(MistletoeClumpsOrBranches_count_yr2_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PIPJ_data_yr2)

gam.check(PIPJ_m8)
# PIPJ_m8 diagnostic plots look decent; same basic patterns as those
# displayed in diagnostic plots of PIPJ_m5.

summary(PIPJ_m8)
plot(PIPJ_m8, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PIPJ_m5. Mistletoe_count is estimated to
# have a negative, linear relationship with cone count (as the level
# of mistletoe infestation increases, cone count decreases), though
# the partial effects plot displays large confidence intervals around
# this estimated effect, indicating some uncertainty in the trend. 
# DBH, CrownDamage, and Mistletoe_count are all "significant" terms
# at alpha = 0.05 (p = 0.049 for Mistletoe_count); DBH:CrownDamage is not.
# "NORTH" is the median fecundity site.

PIPJ_m8_concurvity <- concurvity(PIPJ_m8, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).

## --> none of the auxiliary variables were found to be significant
## (at alpha = 0.05) in any of the yr1 PIPJ models. In contrast,
## Mistletoe_count was found to be a significant predictor of yr2
## PIPJ cone counts (p = 0.049), and CanopyPosition was found to be a 
## marginally significant predictor of yr2 PIPJ cone counts (p = 0.08).

## --> based on the lack of auxiliary variable significance in the yr1
## models, any trees with missing auxiliary variable data, previously 
## removed from the yr1 dataset (ie. CAPLES trees that were windthrown 
## between sample years 1 and 2), can be added back into the dataset, and 
## the yr1 base model can be refit using this more complete, yr1 dataset
## (see below).


## refitting the yr1 base model, using the full yr1 dataset...

# first, saving a new yr1 dataset, inclusive of all previously-removed
# trees...
PIPJ_data_yr1_full <- model_data_scaled %>% 
  filter(SppCode %in% c("PIJE", "PIPO")) %>% 
  filter(TreeStatus_yr1 == "A")

# then, refitting the yr1 base model using these data...
PIPJ_myr1 <- gam(TotalConesConservative_count_yr1 ~
                   s(DBH_cm_scaled, k = 5) +
                   s(CrownDamage_percentvolume_scaled, k = 5) +
                   ti(DBH_cm_scaled, 
                      CrownDamage_percentvolume_scaled, k = 3) +
                   SiteName,
                 family = nb,
                 method = "REML",
                 data = PIPJ_data_yr1_full)

# quickly checking the refit, yr1 base model's diagnostic plots,
# summary table, partial effects plots, and estimated pairwise 
# concurvity values...
gam.check(PIPJ_myr1)
summary(PIPJ_myr1)
plot(PIPJ_myr1, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
PIPJ_myr1_concurvity <- concurvity(PIPJ_myr1, full = FALSE)$estimate
# all diagnostic plots look good; no issues with concurvity; estimated
# effects are stable/consistent across PIPJ_m1 and PIPJ_myr1.


## finally, renaming the yr2 base model following the same convention
## used above (ie. PIPJ_myrX) for easier reference during model
## summary table creation (see script #5)...
PIPJ_myr2 <- PIPJ_m5




