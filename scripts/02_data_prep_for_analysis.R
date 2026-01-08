

# Title: Preparing the fire injury-fecundity data for analysis and modeling

# Author: Nina Venuti

# Script inputs: a .csv file of cleaned fire injury-fecundity project data 
# from all three survey years (2020, 2021, and 2022)

# Script outputs: two .csv files, one of which contains the final data 
# used to fit and select models of the effects of fire injury on conifer
# reproductive capacity one- and two-year(s) post-fire, and one of which 
# includes extra variables (covariates and response variable components)
# that are not directly relevant to the model fitting process, but 
# may be relevant to future analyses. 


## loading relevant packages...
library("tidyverse")

## loading cleaned and compiled project data...
cone_data_allyrs <- read_csv("data/FIFE_data_allyrs_CLEAN.csv")

## simplifying the cleaned, compiled dataframe by removing columns 
## unnecessary for planned analyses...
cone_data_allyrs_simple <- cone_data_allyrs %>% 
  select(-c(SampleDate_mo, SampleDate_day, starts_with("Recorder"),
            starts_with("Scoper"), starts_with("Ranger"),
            starts_with("AbortedCones"), starts_with("ABCOCores"),
            FieldNotes, starts_with("Data")))


## adding columns to the dataframe for total conservative and generous 
## cone counts (model response variables)...

# first, accounting for the potential influence of nearby conspecifics on
# focal trees' ground cone counts -- reducing focal trees' ground counts 
# by the estimated proportion of canopy overlap between each focal tree
# and its nearest conspecific...

# converting ">20" values in the NearestConspp_m column to "20", and 
# converting the column into a numeric vector, to enable calculations below...
cone_data_allyrs_simple$NearestConspp_m[cone_data_allyrs_simple$NearestConspp_m == ">20"] <- 20.0
cone_data_allyrs_simple$NearestConspp_m <- as.numeric(cone_data_allyrs_simple$NearestConspp_m) # note, "blank" is coerced to "NA" here

# defining the GroundCount columns (the focus of the functions below)...
GroundCount_cols <- colnames(cone_data_allyrs_simple[21:24])

# defining the GroundCount_reduction function...

# this function calculates 1/2 of the central angle of a sector of a circle 
# (theta), and uses it to calculate the overlap in area between a symmetric
# lens (created by two overlapping circles of the same radius (r), whose
# centers are some distance (d) away from each other) and the focal circle
# (proportion_overlap). It then reduces the number of cones counted on the 
# ground within the bounds of the focal circle (ie. within the focal tree's
# dripline) by the calculated proportion_overlap.

# function inputs:
# d = the distance between the centers of two circles; for our purposes, 
#     d = NearestConspp_m (a numeric vector), the distance between the boles 
#     of each focal tree and its nearest, mature (≥40 cm) conspecific.
# r = the radius of the focal circle (whose nearest neighbor is assumed 
#     to have the same radius); for our purposes, r = DriplineRadius_m (a
#     numeric vector), the maximum dripline radius of each focal tree.
# groundcount = one of four GroundCount columns (all numeric vectors) = the 
#     number of cones tallied on the ground within each focal tree's 
#     maximum dripline radius.

# function output:
# groundcount_updated (a numeric vector) = new groundcount values (original
#     values, reduced by the estimated proportion_overlap of each focal tree
#     with its nearest conspecific).

GroundCount_reduction <- function(d, r, groundcount) {
  theta = acos(0.5*d/r)
  proportion_overlap = (2*r*theta - d*sin(theta))/(pi*r)
  groundcount_updated = groundcount/(1 + proportion_overlap)
  return(groundcount_updated)
}

# defining a second GroundCount function, to enable application of the 
# GroundCount_reduction function under specific conditions...

# this function applies the GroundCount_reduction function under specific
# conditions (when focal tree and nearest conspecific canopies overlap) and
# skips reduction calculations when such conditions are not met (when there 
# is no canopy overlap), returning the original, unchanged groundcount values
# in the groundcount_updated output. This function is purposely written to 
# loop through each row of the dataframe sequentially, to avoid NaN values
# erroneously produced when the function is vectorized.

# function inputs: see GroundCount_reduction function documentation above.

# function output: 
# groundcount_updated (a numeric vector) = either the original, unchanged
#     groundcount values (when conditions for the GroundCount_reduction
#     function are not met), or the new, reduced groundcount values (when
#     conditions for the GroundCount_reduction function are met).

GroundCount_reduction_w_conditions <- function(d, r, groundcount) {
  n_values = length(groundcount)
  groundcount_updated = vector(mode = "numeric", length = n_values)
  for(i in 1:n_values) {
    groundcount_updated[i] = ifelse(
      test = d[i] < 2*r[i],
      yes = GroundCount_reduction(d[i], r[i], groundcount[i]),
      no = groundcount[i]
    )
  }
  return(groundcount_updated)
}

# applying the GroundCount_reduction_w_conditions function to all four 
# GroundCount columns in the dataframe, and saving updated values in their
# own columns...
cone_data_allyrs_overlapcorrection <- cone_data_allyrs_simple %>% 
  mutate(across(all_of(GroundCount_cols),
                .fns = ~GroundCount_reduction_w_conditions(
                  d = NearestConspp_m,
                  r = DriplineRadius_m,
                  groundcount = .x),
                .names = "{.col}_UPDATED"))
  
# finally, computing total conservative and generous cone counts by summing
# across relevant count columns, and rounding totals to the nearest integer;
# note that we define the total conservative cone count as the sum of canopy-
# tallied healthy, ripening cones and ground-tallied (overlap-corrected)
# cut/stripped cones and freshly-opened cones; we define the total 
# generous cone count as the sum of these same categories, plus canopy-
# tallied unhealthy, ripening cones, canopy-tallied open cones of an 
# unknown age, and ground-tallied (overlap-corrected) open cones of an 
# unknown age... 
cone_data_allyrs_totalcounts <- cone_data_allyrs_overlapcorrection %>% 
  group_by(SampleDate_yr, TagNumber) %>% 
  mutate(TotalConesConservative_count = round(RipeningConesHealthy_count +
          GroundCount_StrippedOrMature_UPDATED + 
          GroundCount_OpenFresh_UPDATED)) %>% 
  mutate(TotalConesGenerous_count = round(RipeningConesHealthy_count +
          RipeningConesUnhealthy_count +
          OpenConesUnkAge_count +
          GroundCount_StrippedOrMature_UPDATED +
          GroundCount_OpenFresh_UPDATED +
          GroundCount_OpenUnkAge_UPDATED)) %>% 
  ungroup()
# note that TotalConesGenerous_count values are "NA" for rows pertaining 
# to the 2020 survey of CAPLES trees, because data on RipeningConesUnhealthy_
# counts were not collected until 2021.

    
## pivoting the dataframe wider, to allow each row in the dataframe to 
## represent a single, tagged tree, with repeat measures of canopy health and
## cone production stored in successive columns...

# assigning labels for survey-effort-year (survey yr1, yr2, or yr3 at each 
# field site) to each row in the dataframe, to facilitate pivot...
yr1 <- c("CAPLES:2020", "AUGUST:2021", "CREEK:2021", "NORTH:2021")
yr2 <- c("CAPLES:2021", "AUGUST:2022", "CREEK:2022", "NORTH:2022")
yr3 <- c("CAPLES:2022")

cone_data_allyrs_totalcounts_wlabels <- cone_data_allyrs_totalcounts %>% 
  mutate("SiteName:SampleDate_yr" = paste(SiteName, 
                                          SampleDate_yr, 
                                          sep = ":")) %>% 
  mutate("SurveyEffort_yr" = case_when(
    `SiteName:SampleDate_yr` %in% yr1 ~ "yr1",
    `SiteName:SampleDate_yr` %in% yr2 ~ "yr2",
    `SiteName:SampleDate_yr` %in% yr3 ~ "yr3"))

# extracting yr1 tree metrics (collected only once, during the initial survey 
# of trees at each field site), to be joined with the wider-format dataframe 
# below...
cone_data_yr1vars <- cone_data_allyrs_totalcounts_wlabels %>% 
  filter(SurveyEffort_yr == "yr1") %>% 
  select(SiteName, SppCode, TagNumber, TreeLatitude, TreeLongitude,
         DBH_cm, BoleCharMin_m, BoleCharMax_m, TotalHeight_m, CanopyPosition,
         PreFireNeedleHeight_m, CrownDamage_percent)

# pivoting the remaining, repeat tree metric columns wider...
cone_data_repeatvars <- cone_data_allyrs_totalcounts_wlabels %>% 
  select(SurveyEffort_yr, TagNumber, TreeStatus, GreenNeedleHeight_m,
         BrownNeedleVolumePostYr1_percent, starts_with("RipeningCones"),
         OpenConesUnkAge_count, DriplineRadius_m, starts_with("GroundCount"),
         NearestConspp_m, MistletoeClumpsOrBranches_count, 
         starts_with("TotalCones")) %>% 
  group_by(TagNumber) %>% 
  pivot_wider(names_from = SurveyEffort_yr,
              values_from = c(3:21))

# merging yr1 tree metrics with wider-format repeat metrics...
cone_data_allyrs_wide <- left_join(cone_data_yr1vars,
                                   cone_data_repeatvars,
                                   by = "TagNumber")


## adding columns to the dataframe for fire injury summary statistics (model
## predictor variables)...

# first, replacing "NA" values in the GreenNeedleHeight_m_yr2 and _yr3 columns
# with TotalHeight_m values for standing dead trees, to facilitate 
# CrownDead_percentlength and GreenCanopyLength_m calculations for those 
# sample years...
cone_data_allyrs_greencorrection <- cone_data_allyrs_wide %>% 
  group_by(TagNumber) %>% 
  mutate(GreenNeedleHeight_m_yr2 = case_when(
    TreeStatus_yr2 == "D" ~ TotalHeight_m,
    .default = GreenNeedleHeight_m_yr2)) %>% 
  mutate(GreenNeedleHeight_m_yr3 = case_when(
    TreeStatus_yr3 == "D" ~ TotalHeight_m,
    .default = GreenNeedleHeight_m_yr3))

# then, calculating fire injury summary variables (BoleCharMeanHeight_m, 
# CrownDamage_percentlength, CrownDead_percentlength_yr2 and _yr3, 
# and GreenCanopyLength_m_yr1, _yr2, and _yr3), and renaming, deleting, and
# reorganizing columns, as needed...
cone_data_allyrs_summarystats <- cone_data_allyrs_greencorrection %>% 
  mutate(BoleCharMeanHeight_m = (BoleCharMin_m + BoleCharMax_m)/2) %>% 
  mutate(CrownDamage_percentlength = (GreenNeedleHeight_m_yr1 -
    PreFireNeedleHeight_m)/(TotalHeight_m - PreFireNeedleHeight_m) * 100) %>% 
  mutate(CrownDead_percentlength_yr2 = (GreenNeedleHeight_m_yr2 -
    PreFireNeedleHeight_m)/(TotalHeight_m - PreFireNeedleHeight_m) * 100) %>% 
  mutate(CrownDead_percentlength_yr3 = (GreenNeedleHeight_m_yr3 -
    PreFireNeedleHeight_m)/(TotalHeight_m - PreFireNeedleHeight_m) * 100) %>% 
  mutate(GreenCanopyLength_m_yr1 = TotalHeight_m - GreenNeedleHeight_m_yr1) %>%
  mutate(GreenCanopyLength_m_yr2 = TotalHeight_m - GreenNeedleHeight_m_yr2) %>%
  mutate(GreenCanopyLength_m_yr3 = TotalHeight_m - GreenNeedleHeight_m_yr3) %>%
  rename(CrownDamage_percentvolume = CrownDamage_percent,
         BrownNeedle_percentvolume_yr2 = BrownNeedleVolumePostYr1_percent_yr2,
         BrownNeedle_percentvolume_yr3 = BrownNeedleVolumePostYr1_percent_yr3) %>% 
  select(-BrownNeedleVolumePostYr1_percent_yr1) %>% 
  select(c(1:8, 69, 10, 9, 11:18, 70:75, 19:68))


## making final adjustments to the model dataframe by collapsing the
## four CanopyPosition categories into two categories (to account for
## low sample sizes in certain categories)...

# converting the few (n = 3) "IN" (Intermediate) trees in the dataset into 
# "CO" (Co-dominant) trees, and converting the few (n = 4) "OP" (Open grown)
# trees in the dataset into "DO" (Dominant) trees; this lumps all trees
# in the dataset into two CanopyPosition categories - one that represents
# "lower light" conditions (encompassing IN and CO trees) and one that
# represents "higher light" conditions (encompassing OP and DO trees)...
cone_data_allyrs_summarystats <- cone_data_allyrs_summarystats %>% 
  mutate(CanopyPosition = case_when(CanopyPosition == "IN" ~ "CO",
                                    CanopyPosition == "OP" ~ "DO",
                                    .default = CanopyPosition))


## exporting this dataframe (the extended version of the model data, 
## inclusive of variables that are not used for final model fitting and
## selection) as a .csv...
# write.csv(cone_data_allyrs_summarystats,
#           "data/processed/FIFE_model_data_EXTENDED.csv", 
#           row.names = FALSE)



## creating a more compact version of the model dataframe, which includes
## only those variables relevant to the final model fitting and selection
## process...

# selecting only those columns that are relevant for final model
# fitting and selection...
cone_data_allyrs_summarystats_slim <- cone_data_allyrs_summarystats %>% 
  select(c(1:6, 9:10, 13:15, 67:68, 70:71, 73:74))

# exporting this dataframe (the final model data) as a .csv...
# write.csv(cone_data_allyrs_summarystats_slim,
#           "data/processed/FIFE_model_data_FINAL.csv", 
#           row.names = FALSE)


