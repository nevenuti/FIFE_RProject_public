

# Title: Cleaning and compiling data from the fire injury-fecundity project

# Author: Nina Venuti

# Script inputs: .csv files of raw data collected in 2020, 2021, and 2022 
# for the fire injury-fecundity project

# Script outputs: .csv files of cleaned (error-corrected) data from 2020, 
# 2021, and 2022, and a .csv file of cleaned data from all three survey
# years, combined


## loading relevant packages...
library("tidyverse")

## loading raw data files...
cone_data_2020 <- read_csv("data/raw/FIFE_data_2020_RAW.csv")
cone_data_2021 <- read_csv("data/raw/FIFE_data_2021_RAW.csv")
cone_data_2022 <- read_csv("data/raw/FIFE_data_2022_RAW.csv")


## preparing to fix erroneous values in the 2020, 2021, and 2022 raw
## data files...

# extracting original values for SppCode from the 2020 and 2021 dataframes
# to replace erroneous values recorded in the 2021 and 2022 dataframes...
Tree031_SppCode <- filter(cone_data_2020, TagNumber == "031")$SppCode
Tree054_SppCode <- filter(cone_data_2020, TagNumber == "054")$SppCode
Tree142_SppCode <- filter(cone_data_2020, TagNumber == "142")$SppCode
Tree168_SppCode <- filter(cone_data_2021, TagNumber == "168")$SppCode
Tree171_SppCode <- filter(cone_data_2021, TagNumber == "171")$SppCode
Tree177_SppCode <- filter(cone_data_2021, TagNumber == "177")$SppCode

# extracting corrected values for DBH_cm and NearestConspp_m from the 
# 2021 dataframe to replace erroneous values recorded in the 2020 and 
# 2022 dataframes...
ABCO147_DBH <- filter(cone_data_2021, TagNumber == "147")$DBH_cm
ABCO051_Conspp <- filter(cone_data_2021, TagNumber == "051")$NearestConspp_m
ABCO108_Conspp <- filter(cone_data_2021, TagNumber == "108")$NearestConspp_m
PILA195_Conspp <- filter(cone_data_2021, TagNumber == "195")$NearestConspp_m


## cleaning each of the raw data files, in turn, by removing rows with 
## NA's (purposely left blank during data entry, to aid in data 
## transcription and QAQC), and replacing erroneous values with 
## corrected values, as appropriate...

# cleaning the 2020 dataframe...
cone_data_2020_clean <- cone_data_2020 %>% 
  filter(!is.na(SampleStartDate_yr)) %>% 
  mutate(DBH_cm = case_when(TagNumber == "147" ~ as.numeric(ABCO147_DBH),
                            .default = DBH_cm)) %>% 
  mutate(NearestConspp_m = case_when(TagNumber == "051" ~ ABCO051_Conspp,
                                     TagNumber == "108" ~ ABCO108_Conspp,
                                     .default = NearestConspp_m))

# cleaning the 2021 dataframe...
cone_data_2021_clean <- cone_data_2021 %>% 
  filter(!is.na(SampleDate_yr)) %>% 
  mutate(SppCode = case_when(TagNumber == "142" ~ Tree142_SppCode,
                             .default = SppCode))

# cleaning the 2022 dataframe...
cone_data_2022_clean <- cone_data_2022 %>% 
  filter(!is.na(SampleDate_yr)) %>% 
  mutate(SppCode = case_when(TagNumber == "031" ~ Tree031_SppCode,
                             TagNumber == "054" ~ Tree054_SppCode,
                             TagNumber == "168" ~ Tree168_SppCode,
                             TagNumber == "171" ~ Tree171_SppCode,
                             TagNumber == "177" ~ Tree177_SppCode,
                             .default = SppCode)) %>% 
  mutate(NearestConspp_m = case_when(TagNumber == "195" ~ PILA195_Conspp,
                                     .default = NearestConspp_m))


## exporting each of the cleaned dataframes as .csv's...

# exporting the 2020 dataframe...
write.csv(cone_data_2020_clean,
          "data/cleaned/FIFE_data_2020_CLEAN.csv",
          row.names = FALSE)

# exporting the 2021 dataframe...
write.csv(cone_data_2021_clean,
          "data/cleaned/FIFE_data_2021_CLEAN.csv",
          row.names = FALSE)

# exporting the 2022 dataframe...
write.csv(cone_data_2022_clean,
          "data/cleaned/FIFE_data_2022_CLEAN.csv",
          row.names = FALSE)


## readying each of the cleaned dataframes for merge into a single
## dataframe, by ensuring column names and classes match across 
## all dataframes...

# preparing the 2020 dataframe for merge...
# first, by copying data on two key predictor variables (CanopyPosition,
# MistletoeClumpsOrBranches_count) from the 2021 dataframe and 
# appending these data to the 2020 dataframe (data for these two 
# variables were not collected until 2021 for CAPLES trees, but 
# values are thought to be relatively stable across years)...
caples_trees_missing_predictors <- cone_data_2021_clean %>% 
  filter(SiteName == "CAPLES") %>% 
  select(TagNumber, CanopyPosition, MistletoeClumpsOrBranches_count)

cone_data_2020_allpredictors <- left_join(cone_data_2020_clean,
                                          caples_trees_missing_predictors,
                                          by = "TagNumber")

# then, by removing unneeded date columns, and renaming and reclassifying 
# columns, as needed...
chr_cols_2020 <- colnames(cone_data_2020_allpredictors[c(33:35, 48)])

cone_data_2020_formerge <- cone_data_2020_allpredictors %>% 
  select(-c(SampleStartDate_yr, SampleStartDate_mo, SampleStartDate_day)) %>%
  rename(SampleDate_yr = SampleEndDate_yr,
         SampleDate_mo = SampleEndDate_mo,
         SampleDate_day = SampleEndDate_day,
         DataEnteredBy1 = DataEnteredBy) %>% 
  mutate_at(chr_cols_2020, as.numeric)
# note that the mutate_at statement above converts "blank" values in 
# select character columns (chr_cols_2020) to "NA" values, which is 
# ok for our purposes ("blank" is used in all of the raw data files for 
# this project to mean "no data available" or "column not applicable", 
# so "NA" is appropriate here). 

# preparing the 2021 dataframe for merge by reclassifying columns, 
# as needed...
chr_cols_2021 <- colnames(cone_data_2021_clean[c(19:24, 27:38, 40)])

cone_data_2021_formerge <- cone_data_2021_clean %>% 
  mutate_at(chr_cols_2021, as.numeric) # again, "blank" is coerced to "NA" 

# preparing the 2022 dataframe for merge by renaming and reclassifying
# columns, as needed...
chr_cols_2022 <- colnames(cone_data_2022_clean[c(23:35, 37)])
num_cols_2022 <- colnames(cone_data_2022_clean[3])

cone_data_2022_formerge <- cone_data_2022_clean %>% 
  mutate_at(chr_cols_2022, as.numeric) %>% # again, "blank" is coerced to "NA"
  mutate_at(num_cols_2022, as.character) %>% 
  rename(MistletoeClumpsOrBranches_count = MistletoeInfestedBranches_count,
         DataEnteredBy1 = DataEnteredBy)


## merging all three dataframes into one, reordering columns to keep like
## columns together, and replacing "NA" values with "blank" values for select
## character columns, for purposes of consistency (and in preparation for 
## dataset publication)...
cone_data_allyrs <- bind_rows(cone_data_2020_formerge,
                              cone_data_2021_formerge, 
                              cone_data_2022_formerge) %>% 
  select(c(1:8, 49, 9:12, 50, 13:16, 51, 17:18, 25, 19:24, 44, 26:28, 
           46, 29, 47, 30:32, 52:54, 33:39, 45, 40:42, 48, 43)) %>% 
  mutate(across(all_of(c("Recorder5", "Scoper5", "Ranger5")),
                .fns = ~case_when(SampleDate_yr %in% c(2020, 2021) ~ "blank",
                                  .default = .x))) %>% 
  mutate(DataEnteredBy2 = case_when(SampleDate_yr %in% c(2020, 2022) ~ "blank",
                                    .default = DataEnteredBy2))


## exporting this merged dataframe as a .csv...
write.csv(cone_data_allyrs,
          "data/cleaned/FIFE_data_allyrs_CLEAN.csv",
          row.names = FALSE)


