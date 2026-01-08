

# Title: Modeling the effects of fire injury on PILA (sugar pine) 
# fecundity one- and two-year(s) post-fire

# Author: Nina Venuti

# Script inputs: a .csv file of fire injury-fecundity project data, 
# cleaned and prepared for modeling

# Script outputs: a series of fitted generalized additive models (GAMs)
# of PILA cone counts one- and two-year(s) after fire injury


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
# live PILA trees in the first and second year of sampling; modeling will
# be conducted separately for each sample year, using data collected on
# live trees only (ie. excluding trees that were dead, fallen, cut down, 
# broken-topped, or top-killed in either yr1 or yr2). Note that no
# additional trees need to be removed from the yr1 or yr2 datasets
# at this point, as none of the PILA trees in these datasets are missing
# data for key predictor variables.
PILA_data_yr1 <- model_data_scaled %>% 
  filter(SppCode == "PILA") %>% 
  filter(TreeStatus_yr1 == "A")

PILA_data_yr2 <- model_data_scaled %>% 
  filter(SppCode == "PILA") %>% 
  filter(TreeStatus_yr2 == "A")


## calculating and inspecting correlation coefficients for pairs of predictor
## variables, to make sure predictors are not highly correlated with one
## another and may be included in the same models... 

# starting with the yr1 data...
# selecting key predictor variables for which to calculate correlation 
# coefficients...
PILA_correlations_yr1 <- PILA_data_yr1 %>% 
  select(c(6:9, 12))

# converting CanopyPosition (a factor variable) into a numeric variable
# temporarily, by assigning a score of either 2 (for "DO") or 
# 1 (for "CO") to each tree, to enable correlation coefficient
# calculations; afterwards, removing the original CanopyPosition column...
PILA_correlations_yr1 <- PILA_correlations_yr1 %>% 
  mutate(CanopyPositionScore = case_when(CanopyPosition == "DO" ~ 2,
                                         CanopyPosition == "CO" ~ 1)) %>% 
  select(-CanopyPosition)

# saving the yr1 correlations data.frame as a matrix...
PILA_correlations_yr1_matrix <- as.matrix(PILA_correlations_yr1)

# calculating correlation coefficients for all of the yr1 variable pairs...
PILA_correlations_yr1_pairs <- cor(PILA_correlations_yr1_matrix,
                                   use = "pairwise.complete.obs")

# making a correlation plot, to enable easy inspection of correlation
# coefficients...
corrplot(PILA_correlations_yr1_pairs,
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black")
# --> this correlation plot demonstrates that yr1 predictor variables
# are only weakly or moderately correlated with one another (ie. all
# correlation coefficients are <0.6), and thus it should be ok to 
# include them in the same models. DBH_cm and BoleCharMeanHeight_m 
# are moderately correlated with each other (coeff = 0.43), as are 
# CrownDamage_percentvolume and BoleCharMeanHeight_m (coeff = 0.42),
# and DBH_cm and CanopyPositionScore (coeff = 0.32).

# repeating this process with the yr2 data...
# selecting key predictor variables for which to calculate correlation 
# coefficients...
PILA_correlations_yr2 <- PILA_data_yr2 %>% 
  select(c(6:9, 13))

# converting CanopyPosition (a factor variable) into a numeric variable
# temporarily, by assigning a score of either 2 (for "DO") or 
# 1 (for "CO") to each tree, to enable correlation coefficient
# calculations; afterwards, removing the original CanopyPosition column...
PILA_correlations_yr2 <- PILA_correlations_yr2 %>% 
  mutate(CanopyPositionScore = case_when(CanopyPosition == "DO" ~ 2,
                                         CanopyPosition == "CO" ~ 1)) %>% 
  select(-CanopyPosition)

# saving the yr2 correlations data.frame as a matrix...
PILA_correlations_yr2_matrix <- as.matrix(PILA_correlations_yr2)

# calculating correlation coefficients for all of the yr2 variable pairs...
PILA_correlations_yr2_pairs <- cor(PILA_correlations_yr2_matrix,
                                   use = "pairwise.complete.obs")

# making a correlation plot, to enable easy inspection of correlation
# coefficients...
corrplot(PILA_correlations_yr2_pairs,
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black")
# --> this correlation plot demonstrates that yr2 predictor variables
# are only weakly or moderately correlated with one another (ie. all
# correlation coefficients are <0.6), and thus it should be ok to 
# include them in the same models. DBH_cm and BoleCharMeanHeight_m 
# are moderately correlated with each other (coeff = 0.49), as are
# DBH_cm and CanopyPositionScore (coeff = 0.32), and CrownDamage_
# percentvolume and BoleCharMeanHeight_m (coeff = 0.31).



## writing and evaluating models of PILA cone counts...

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
# to test its effects on total cone counts. NOTE, that for PILA, only
# the auxiliary variables BoleCharMeanHeight_m and CanopyPosition
# are tested, due to the very low incidence of mistletoe infestation 
# in PILA across both sample years (n = 2 and n = 6 trees with 
# mistletoe in their canopies in yr1 and yr2, respectively). Of
# the two auxiliary variables that are relevant for PILA models,
# one (BoleChar) is added to the base model as a smooth term, 
# while the other (CanopyPosition) is added as a parametric term. 
# Each smooth term included in these models - except the interaction 
# term, DBH:CrownDamage - is estimated using a maximum of five basis
# functions (ie. the parameter "k" is set to 5), in order to constrain
# the "wiggliness" of each smooth. The interaction term, DBH:CrownDamage, 
# is estimated using a maximum of three basis functions (ie. k = 3), 
# to reduce the estimated pairwise concurvity values between it and 
# the other terms in each model.

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
# is the "worst" value) nor overly-optimistic (as is the "observed"  
# value), according to S. Wood, the author of the mgcv package.


## starting with yr1 PILA models...

## PILA_m1 (yr1 base model)...
PILA_m1 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, 
                    CrownDamage_percentvolume_scaled, k = 3) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PILA_data_yr1)

par(mfrow = c(2,2))
gam.check(PILA_m1)
# PILA_m1 diagnostic plots look pretty good; the histogram of residuals  
# is basically bell-shaped, and there are no strong patterns in the 
# residuals vs. linear predictor plot.

summary(PILA_m1)
plot(PILA_m1, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH is estimated to have a largely positive, nonlinear relationship
# with cone count (as DBH increases, cone count increases, and then
# begins to plateau at high DBH values). CrownDamage is estimated 
# to have a slightly negative, linear relationship with cone count  
# (as CrownDamage increases, cone count decreases). DBH:CrownDamage
# is estimated to have a nonlinear effect on cone count. All three
# smooth terms - DBH, CrownDamage, and DBH:CrownDamage - are 
# "significant" at alpha = 0.05. "NORTH" is the median fecundity
# site for yr1 PILA.

PILA_m1_concurvity <- concurvity(PILA_m1, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PILA_m2 (yr1 base model + BoleCharMeanHeight_m)...
PILA_m2 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, 
                    CrownDamage_percentvolume_scaled, k = 3) +
                 s(BoleCharMeanHeight_m_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PILA_data_yr1)

gam.check(PILA_m2)
# PILA_m2 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of PILA_m1. Note that the basis
# dimension check produced by gam.check() indicates that the number
# of basis functions allowed for both DBH_cm and BoleCharMeanHeight_m
# may be too low to properly estimate the effects of these variables
# on cone count, but the edf values for these terms are not super
# close to k, so it should be ok to proceed (k is also purposely
# capped at 5, to constrain the "wiggliness" of each smooth, as 
# mentioned above in the notes on this project's modeling protocols). 

summary(PILA_m2)
plot(PILA_m2, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PILA_m1. BoleChar is estimated to have
# a nonlinear and overall negative relationship with cone count
# (as BoleChar increases, cone count largely remains steady, and
# then begins to decrease at very high BoleChar values). DBH and 
# CrownDamage are "significant" terms (at alpha = 0.05), 
# DBH:CrownDamage is marginally so (p = 0.059), and BoleChar is 
# not. "NORTH" is the median fecundity site.

PILA_m2_concurvity <- concurvity(PILA_m2, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PILA_m3 (yr1 base model + CanopyPosition)...
PILA_m3 <- gam(TotalConesConservative_count_yr1 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, 
                    CrownDamage_percentvolume_scaled, k = 3) +
                 CanopyPosition +
                 SiteName,
               family = nb,
               method = "REML",
               data = PILA_data_yr1)

gam.check(PILA_m3)
# PILA_m3 diagnostic plots look good; same basic patterns as those
# displayed in diagnostic plots of PILA_m1. Note, again, that the
# basis dimension check produced by gam.check() indicates that the
# number of basis functions allowed for DBH_cm may be too low to
# properly estimate the effects of DBH on cone count, but the 
# edf value for this term is not super close to k, so it should be
# ok to proceed (and, k is purposely capped at 5 to constrain 
# the "wiggliness" of each smooth in the model).

summary(PILA_m3)
plot(PILA_m3, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PILA_m1. All three smooth terms - DBH, 
# CrownDamage, and DBH:CrownDamage - are "significant" (at alpha =
# 0.05). "DO" (dominant) trees are estimated to produce slightly more
# cones than "CO" (co-dominant) trees, but this effect is not 
# significant (at alpha = 0.05). "NORTH" is the median fecundity site.

PILA_m3_concurvity <- concurvity(PILA_m3, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).

# **remember, the third auxiliary variable of interest, Mistletoe_count,
# is not added to PILA models, as the incidence of mistletoe in PILA
# canopies was very low in both sample years (n = 2 trees that had 
# mistletoe in their canopies in yr1, and n = 6 trees that had mistletoe 
# in their canopies in yr2). Thus, we can move on to modeling yr2
# cone counts at this point.**


## moving on to yr2 PILA models... 

## note that we've skipped the label "PILA_m4" (which would have been
## assigned to the yr1 base model + Mistletoe_count model), and started
## with the label "PILA_m5" below, for purposes of consistency, and ease
## of summary table creation (see script #5)...

## PILA_m5 (yr2 base model)...
PILA_m5 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, 
                    CrownDamage_percentvolume_scaled, k = 3) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PILA_data_yr2)

gam.check(PILA_m5)
# PILA_m5 diagnostic plots look decent; the histogram of residuals
# is relatively bell-shaped (though a bit right-skewed), and there 
# are no strong patterns in the residuals vs. linear predictor plot
# (though there may be a slight negative trend in mean residual
# values as fitted values increase).

summary(PILA_m5)
plot(PILA_m5, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH is estimated to have a largely positive, nonlinear relationship
# with cone count (as DBH increases, cone count increases, and then
# begins to plateau at high DBH values). CrownDamage is estimated to
# have a slightly negative, nonlinear relationship with cone count
# (as CrownDamage increases, cone count remains stable, and then begins
# to decrease slightly at low to moderate CrownDamage values).
# DBH:CrownDamage is estimated to have a linear effect on cone count.
# DBH is "significant" (at alpha = 0.05), CrownDamage and
# DBH:CrownDamage are not. "CREEK" is the median fecundity site
# for yr2 PILA, though it appears to be indistinguishable from 
# "NORTH".

PILA_m5_concurvity <- concurvity(PILA_m5, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PILA_m6 (yr2 base model + BoleCharMeanHeight_m)...
PILA_m6 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, 
                    CrownDamage_percentvolume_scaled, k = 3) +
                 s(BoleCharMeanHeight_m_scaled, k = 5) +
                 SiteName,
               family = nb,
               method = "REML",
               data = PILA_data_yr2)

gam.check(PILA_m6) 
# PILA_m6 diagnostic plots look pretty good; same basic patterns as those
# displayed in diagnostic plots of PILA_m5.

summary(PILA_m6)
plot(PILA_m6, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PILA_m5. BoleChar is estimated to have a 
# slightly positive, linear relationship with cone count (as 
# BoleChar increases, cone count increases slightly), though the
# estimated effect appears to be very close to a horizontal line
# in the partial effects plot. DBH is "significant" (at alpha = 0.05),
# CrownDamage, DBH:CrownDamage, and BoleChar are not. "CREEK" is
# the median fecundity site, though it appears to be indistinguishable 
# from "NORTH", or even "AUGUST", in this case.

PILA_m6_concurvity <- concurvity(PILA_m6, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).


## PILA_m7 (yr2 base model + CanopyPosition)
PILA_m7 <- gam(TotalConesConservative_count_yr2 ~
                 s(DBH_cm_scaled, k = 5) +
                 s(CrownDamage_percentvolume_scaled, k = 5) +
                 ti(DBH_cm_scaled, 
                    CrownDamage_percentvolume_scaled, k = 3) +
                 CanopyPosition +
                 SiteName,
               family = nb,
               method = "REML",
               data = PILA_data_yr2)

gam.check(PILA_m7)
# PILA_m7 diagnostic plots look pretty good; same basic patterns as those
# displayed in diagnostic plots of PILA_m5.

summary(PILA_m7)
plot(PILA_m7, all.terms = TRUE, shade = TRUE, scheme = 1, pages = 1)
# DBH, CrownDamage, and DBH:CrownDamage follow the same patterns
# as those described for PILA_m5. DBH is "significant" (at alpha = 
# 0.05), CrownDamage and DBH:CrownDamage are not. "DO" (dominant)
# trees are estimated to produce slightly fewer cones than "CO" (co-
# dominant) trees, but this effect is not significant (at alpha =
# 0.05). "NORTH" is now the median fecundity site, though it appears 
# to be indistinguishable from "CREEK".

PILA_m7_concurvity <- concurvity(PILA_m7, full = FALSE)$estimate
# no issues with concurvity here (no pairwise concurvity values
# are >0.6).

# **remember, the third auxiliary variable of interest, Mistletoe_count,
# is not added to PILA models, as the incidence of mistletoe in PILA
# canopies was very low in both sample years (n = 2 trees that had 
# mistletoe in their canopies in yr1, and n = 6 trees that had mistletoe 
# in their canopies in yr2). Thus, we can conclude the process of
# fitting PILA models with PILA_m7.**

## --> neither of the two auxiliary variables that were tested -
## BoleChar and CanopyPosition - were found to be significant (at
## alpha = 0.05) in any of the yr1 or yr2 PILA models.


## as a last step, renaming the yr1 and yr2 base models following the 
## same convention used for ABCO and PIPJ base models (ie. SPPCODE_myrX)
## for easier reference during model summary table creation (see
## script #5)...
PILA_myr1 <- PILA_m1
PILA_myr2 <- PILA_m5




