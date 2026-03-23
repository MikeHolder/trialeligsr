# Title: Systematic review evidence synthesis
# Author: Mike Holder
# Date: 2026-03-20
# Description: All code for evidence synthesis for systematic review of 
#              differences between trial-eligible and trial-ineligible groups
#              in real-world populations

# Load libraries ####

library(tidyverse)
library(janitor)
library(ggrepel)
library(readxl)
library(cowplot)
library(DT)
library(purrr)
library(tools)
library(kableExtra)
library(ggsci)
library(RColorBrewer)
library(gtsummary)
library(officer)
library(flextable)
library(ggforce)
library(patchwork)
library(ggh4x)
library(Hmisc)
library(here)

# Setup ####

# avoid scientific notation
options(scipen = 9999, digits = 4)

# set working directory
# setwd("path/to/your/folder") to add for your own folder

# Read in data ####

# Read in studies sheet
Studies <- read_excel("data.xlsx", sheet = "Studies") %>%
  # clean the data
  clean_names() %>%
  # simplify column names
  rename(
    "condition" = condition_s_of_interest,
    "category" = condition_category,
    "study_author" = first_author
  ) %>%
  # create column of broader group of specialties
  mutate(category_2 = factor(
    str_replace_all(
      category_2,
      c(
        "Ischaemic CVD" = "CVD",
        "HF" = "CVD",
        "AF" = "CVD",
        "ITU" = "CVD",
        "Renal" = "CVD" # This was a hypertension trial in patients with CKD
      )
    ),
    levels = c("CVD", "Cancer", "Resp", "MSK")
  ))

# Read in populations sheet
Populations <- read_excel("data.xlsx", sheet = "Population") %>%
  # clean names
  clean_names() %>%
  # arrange by author alphabetically
  arrange(study_author) %>%
  # Rename population name column
  rename("population" = name) %>%
  # remove NA values
  filter(!is.na(population))

# smaller dataset
Studiesbrief <- Studies %>% 
  select(study_author, condition, category_2)

# Read in trials dataset, including duplicates
Trials0 <- read_excel("data.xlsx", sheet = "Trials") %>% 
  # remove NAs
  filter(!is.na(`Trial Name`)) %>% 
  # clean column names
  clean_names() %>% 
  # better names
  rename("population" = population_name,
         "number" = trial_number,
         "trial" = trial_name,
         "trial_author" = first_author,
         "intervention" = intervention_s,
         "rob_representativeness" = rob_representativeness_to_target_population,
         "trial_rob" = trial_population_pair_risk_of_bias,
         "comorbs_excluded" = selected_co_morbidities_excluded,
         "participants" = number_of_participants_studied_randomised,
         "average_age" = trial_average_age,
         "prop_women" = trial_proportion_women_percent,
         "prop_ineligible" = proportion_of_cohort_i_neligible_percent,
         "eligible_average_age" = average_age_of_eligible,
         "eligible_prop_women" = proportion_of_women_among_eligible_percent,
         "ineligible_average_age" = average_age_of_ineligible,
         "ineligible_prop_women" = proportion_of_women_among_ineligible_percent,
         "trial_year" = year
  ) %>% 
  # arrange by prop ineligible (see next snippet for reason)
  arrange(prop_ineligible)
  
# deduplicate Trials data
Trials <- Trials0 %>% 
  # join to studies
  left_join(Studiesbrief, by = "study_author") %>% 
  distinct(trial, .keep_all = TRUE) %>% # previously sorted by prop ineligible
  # so first trial will be lowest number, therefore most conservative estimate
  
  # ensure numberic columns
  mutate(prop_ineligible = as.numeric(prop_ineligible),
         prop_women = as.numeric(prop_women),
         eligible_prop_women = as.numeric(eligible_prop_women),
         ineligible_prop_women = as.numeric(ineligible_prop_women),
         average_age = as.numeric(average_age),
         eligible_average_age = as.numeric(eligible_average_age),
         ineligible_average_age = as.numeric(ineligible_average_age)
  ) # will coerce some NAs where not present

# All duplicated trials
Trialsdup <- Trials0 %>% 
  group_by(trial) %>% 
  filter(n() > 1) %>% 
  summarise(trials = n())

# All duplicated populations
Popsdup <- Populations %>% 
  group_by(population) %>% 
  filter(n() > 1) %>% 
  summarise(populations = n())
  
# Read in comparisons sheet
Comparisons <- read_excel("data.xlsx", sheet = "Comparison") %>% 
  # arrange by study
  arrange(`Study Author`) %>% 
  # clean column names
  clean_names() %>%
  rename("population" = population_name,
         "trial" = trial_name)
  

# Big dataset of everything  
All <- Comparisons %>% 
  left_join(Trials0, by = c("study_author", "population", "trial")) %>% 
  full_join(Populations, by = c("study_author", "population")) %>% 
  full_join(Studies, by = "study_author")

# All unique trial-population pairs
TrialPop <- All %>% 
  select(trial, population) %>% 
  distinct(trial, population)


# Tidy Up ####

Tidy <- All %>%
  # make relevant columns numeric
  mutate(
    eligible_measurement_1 = as.numeric(eligible_measurement_1),
    eligible_measurement_2 = as.numeric(eligible_measurement_2),
    ineligible_measurement_1 = as.numeric(ineligible_measurement_1),
    ineligible_measurement_2 = as.numeric(ineligible_measurement_2),
    number_of_participants = as.numeric(number_of_participants),
    prop_ineligible = as.numeric(prop_ineligible)
  ) %>%
  # calculate n eligible and ineligible
  mutate(
    n_ineligible = round(number_of_participants * prop_ineligible / 100, 0),
    n_eligible = round(number_of_participants - n_ineligible),
    0
  ) %>%
  # standardise prop for women and men, and increased/decreased severity
  mutate(
    across(
      contains("_measurement_2"),
      ~ case_when(
        mapped %in% c("Men", "Decreased", "Lives with others")
        ~ 100 - .x,
        mapped == "Survival" &
          eligible_qualifier_2 == "%"
        ~ 100 - .x,
        TRUE ~ .x
      )
    ),
    across(
      contains("_measurement_1"),
      ~ case_when(
        mapped %in% c("Men", "Decreased", "Lives with others")
        ~ n_ineligible - .x,
        nature_of_comparison == "EQ-5D baseline"
        ~ 1 - .x,
        nature_of_comparison == "WIQ baseline"
        ~ 100 - .x,
        TRUE ~ .x
      )
    ),
    # calculate eGFR from creatinine (convert mg/dl to mmol/l)
    # algorithm based on CKD-EPI equation
    #
    # https://www.niddk.nih.gov/research-funding/research-programs/
    # kidney-clinical-research-epidemiology/laboratory/
    # glomerular-filtration-rate-equations/adults
    #
    # assuming female, aged 60 and creatinine measurements in mmol/l
    across(
      contains("_measurement_1"),
      ~ case_when(
        str_detect(nature_of_comparison, "creatinine (mg/dl)")
        ~ round(.x * 88.4, 1),
        .default = .x
      )
    ),
    across(
      contains("_measurement_1"),
      ~ case_when(
        str_detect(nature_of_comparison, "reat") &
          mapped == "Renal function" &
          .x > 30 ~
          round(142 * (.x / 61.88) ^ -1.2 * 0.696819, 1),
        .default = .x
      )
    ),
    # convert haemoglobin to UK units (g/l)
    across(
      contains("_measurement_1"),
      ~ case_when(nature_of_comparison == "Hb (g/dl)"
                  ~ .x * 10, .default = .x)
    ),
    # convert total cholesterol to UK units (mmol/l)
    across(
      contains("_measurement_1"),
      ~ case_when(
        nature_of_comparison == "Total chol (mg/dl)"
        ~ round(.x * 0.02586, 1),
        .default = .x
      )
    )
  ) %>%
  # Survival is measured in numbers and months so need to separate
  mutate(mapped = case_when(
    mapped == "Survival" & str_detect(nature_of_comparison, regex(
      "month", ignore_case = TRUE)) ~ "survival (months)",
    TRUE ~ mapped
  )) %>%
  # create column to standardise comparison names, as some are from subcategory
  # and some are from mapping
  mutate(
    comparison = case_when(
      subcategory %in% c(
        "Comorbidity",
        "Riskfactor",
        "Measurement",
        "Preventative",
        "Treatment",
        "Adverseevent",
        "Lengthofstay",
        "Mortality",
        "Treatmentresponse",
        "Changetotreatment",
        "Other"
      )
      ~ mapped,
      subcategory == "Age" &
        mapped == "Other"
      ~ "Proportion old",
      subcategory == "Living"
      ~ "Lives alone",
      str_detect(nature_of_comparison, "hite") &
        mapped == "Ethnicity"
      ~ "White race",
      str_detect(nature_of_comparison, "lack") &
        mapped == "Ethnicity"
      ~ "Black race",
      nature_of_comparison == "Not Hispanic" &
        mapped == "Ethnicity"
      ~ "For removal",
      str_detect(nature_of_comparison, "ispanic") &
        mapped == "Ethnicity"
      ~ "Hispanic",
      mapped == "Ethnicity"
      ~ "Other global majority",
      subcategory == "Sex" ~ "Women",
      .default = subcategory
    )
  ) %>%
  # calculate n where not given.
  # Calculate NA values (all NAs are where percentage given but no n)
  # done for completeness, shouldn't be necessary for any analysis
  mutate(
    eligible_measurement_1 = if_else(
      is.na(eligible_measurement_1),
      round(eligible_measurement_2 * n_eligible /
              100, 0),
      eligible_measurement_1
    ),
    eligible_qualifier_1 = if_else(eligible_qualifier_1 == "NA", "n", 
                                   eligible_qualifier_1),
    ineligible_measurement_1 = if_else(
      is.na(ineligible_measurement_1),
      round(ineligible_measurement_2 * n_ineligible /
              100, 0),
      ineligible_measurement_1
    ),
    ineligible_qualifier_1 = if_else(ineligible_qualifier_1 == "NA", "n", 
                                     ineligible_qualifier_1),
    measurement = case_when(
      eligible_qualifier_1 %in% c("mean", "Mean") ~ "mean",
      eligible_qualifier_1 %in% c("median", "Median") ~ "median",
      eligible_qualifier_2 == "%" ~ "prop",
      TRUE ~ NA
    ),
    # standardise means and proportions into same column
    elvalue = case_when(
      measurement %in% c("mean", "median")
      ~ eligible_measurement_1,
      measurement == "prop" ~ eligible_measurement_2
    ),
    inelvalue = case_when(
      measurement %in% c("mean", "median")
      ~ ineligible_measurement_1,
      measurement == "prop" ~ ineligible_measurement_2
    ),
  )


# For the synthesis, comparisons will be restricted to %, mean or median  
ComparisonTidy <- Tidy %>% 
  filter(!eligible_qualifier_1 %in% c("HR", "aHR", "IRR", "Rate", "rate",
                                      "cm", "%")) %>% 
  # filter out severity and treatmentresponse rows where direction is unclear
  filter(!(subcategory %in% c("Severity", "Treatmentresponse")
           & positive == "NA")) %>% 
  # remove absent measurements
  filter(!is.na(measurement)) %>% 
  # remove duplicates
  # group by population-trial pair, comparison and measurement (mean/median/%)
  group_by(population, trial, comparison, measurement) %>%
  filter(
    # Keep rows where positive is "No" or where measurements do not match a 
    # "No" row in the same group
    positive == "No" |
      !(eligible_measurement_1 %in% 
          eligible_measurement_1[positive == "No"] &
          eligible_measurement_2 %in% 
          eligible_measurement_2[positive == "No"] &
          ineligible_measurement_1 %in% 
          ineligible_measurement_1[positive == "No"] &
          ineligible_measurement_2 %in% 
          ineligible_measurement_2[positive == "No"])
  ) %>%
  ungroup() %>% 
  select(study_author, population, trial, category_2, treatment_category.y, 
         nature_of_comparison, category.x,subcategory, mapped, comparison, 
         measurement, positive, elvalue, inelvalue, 
         eligible_qualifier_1:ineligible_measurement_2, 
         number_of_participants, trial_rob)

# function for analysis
comparisonsfun <- function(data, condition, treatment) {
  # data to be used
  data %>%
    # group same comparisons together
    group_by(Comparison) %>%
    # summary stats
    summarise(
      nstuds = n_distinct(study_author),
      npops = n_distinct(population),
      ntrials = n_distinct(trial),
      ncomparisons = n(),
      ncondtreat = n_distinct(treatment_category.y),
      # median and IQR for eligible and ineligible
      Eligible = paste0(
        sprintf("%.1f", median(elvalue, na.rm = TRUE)),
        " (",
        sprintf("%.1f", quantile(elvalue, 0.25, na.rm = TRUE)),
        "-",
        sprintf("%.1f", quantile(elvalue, 0.75, na.rm = TRUE)),
        ")"
      ),
      Ineligible = paste0(
        sprintf("%.1f", median(inelvalue, na.rm = TRUE)),
        " (",
        sprintf("%.1f", quantile(inelvalue, 0.25, na.rm = TRUE)),
        "-",
        sprintf("%.1f", quantile(inelvalue, 0.75, na.rm = TRUE)),
        ")"
      ),
      # Wilcox test
      pvalue = wilcox.test(inelvalue, elvalue, 
                           paired = TRUE, # this bit makes it signed rank
                           exact = FALSE)$p.value,
      # take first value only (as all should be the same)
      Category = first(category.x),
      Subcategory = first(subcategory),
      nCVDstuds = n_distinct(study_author[category_2 == "CVD"]),
      nCVDtrial = n_distinct(trial[category_2 == "CVD"]),
      nCancerstuds = n_distinct(study_author[category_2 == "Cancer"]),
      nCancertrial = n_distinct(trial[category_2 == "Cancer"]),
      nRespstuds = n_distinct(study_author[category_2 == "Resp"]),
      nResptrial = n_distinct(trial[category_2 == "Resp"]),
      nMSKstuds = n_distinct(study_author[category_2 == "MSK"]),
      nMSKtrial = n_distinct(trial[category_2 == "MSK"]),
      .groups = "drop" # Ungroup after summarizing
    ) %>%
    mutate(Condition = condition, Treatment = treatment)
}

# Characteristics of studies, populations and trials ####

## Studies ####

StudiesInfo <- All %>% 
  # calculate n trials per study
  group_by(study_author) %>% 
  mutate(ntrials = length(unique(trial))) %>% 
  ungroup() %>% 
  # arrange by risk of bias (high first) so if any trial is high risk of bias,
  # this will show. The nremove duplicates
  arrange(trial_rob) %>% 
  distinct(study_author, .keep_all = TRUE) %>% 
  # select relevant categories
  select(category_2, year, trial_rob, ntrials) %>% 
  # group for decade, factorise risk-of-bias and group for trial numbers
  mutate(
    year = case_when(year < 2011 ~ "2000s",
                           year < 2021 ~ "2010s",
                           year < 2030 ~ "2020s"),
    trial_rob = factor(trial_rob, levels = c("High", "Low")),
    ntrials2 = factor(case_when(ntrials == 1 ~ "1",
                        ntrials == 2 ~ "2",
                        ntrials == 3 ~ "3",
                        ntrials <9 ~ "4-9",
                        .default = "10+"),
                     levels = c("1", "2", "3", "4-9", "10+"))
    )%>% 
  # Rename
  rename(
    "Condition" = category_2,
    "Study year" = year,
    "Risk of bias" = trial_rob,
    "Number of trials" = ntrials2,
    "Number of trials (range)" = ntrials)

# gt summary table
StudiesInfogt <- tbl_summary(
  StudiesInfo,
  digits = list(all_categorical() ~ c(0, 1)),
  statistic =
    list(all_continuous() ~
           "{median} ({min}, {max})")
)

StudiesInfogt

## Populations ####

PopsInfo <- All %>% 
  # remove duplicates
  distinct(population, .keep_all = TRUE) %>% 
  # keep nature and size of population
  select(nature_of_pop_research_clinical,
         number_of_participants) %>% 
  # factorise and group partially logarithmically
  mutate(
    popsize = factor(
      case_when(
        number_of_participants < 100 ~ "<100",
        number_of_participants < 251 ~ "100-250",
        number_of_participants < 1001 ~ "251-1000",
        number_of_participants < 2501 ~ "1001-2500",
        number_of_participants < 10001 ~ "2501-10000",
        .default = "10000+"),
        levels = c("<100", "100-250", "251-1000", "1001-2500", "2501-10000",
                   "10000+"))
  ) %>% 
  rename(
    "Nature" = nature_of_pop_research_clinical,
    "Population size (range)" = number_of_participants,
    "Population size" = popsize)

# gt summary table
PopsInfogt <- tbl_summary(PopsInfo, 
                             digits = list(all_categorical() ~ c(0, 1)),
                             statistic = 
                               list(all_continuous() ~ 
                                      "{median} ({min}, {p25}, {p75}, {max})"))

PopsInfogt


## Trials ####

TrialsInfo <- Trials %>% 
  select(category_2, funding, trial_year, drug_device, comparison_category, 
         location, participants) %>% 
  # group funding
  mutate(funding = case_when(funding == "Private" ~ "Industry",
                             funding == "Unknown" ~ "Unknown",
                             .default = "Public"),
         # year as decades
         trial_year = case_when(trial_year < 2000 ~ "1990s",
                                trial_year < 2010 ~ "2000s",
                                trial_year < 2020 ~ "2010s",
                                trial_year < 2030 ~ "2020s"),
         # rename drug/device variables
         drug_device = factor(if_else(drug_device == "Drug", 
                                         "Drug", "Device/Other"),
                                 levels = c("Drug", "Device/Other")),
         # regroup comparisons
         comparison_category = factor(case_when(comparison_category == "Placebo"
                                          ~ "Placebo",
                                          comparison_category == "NA"
                                          ~ "No control",
                                          .default = "Active"),
                                      levels = c("Active", "Placebo",
                                                 "No control")),
         # regroup location
         location = factor(
           case_when(
             location == "International" ~ "Multiple continents",
             location %in% c("Austrlia", "Japan", "New Zealand", 
                             "Oceania", "South Korea") ~ "Asia/Oceania",
             location %in% c("North America", "USA") ~ "North America",
             .default = "Europe"
           ), levels = c("Multiple continents", "Europe", "North America",
                         "Asia/Oceania")
         ))%>% 
  rename("Condition Group" = category_2,
         "Funding" = funding,
         "Trial decade" = trial_year,
         "Intervention" = drug_device,
         "Control" = comparison_category,
         "Location" = location,
         "Participants (median [IQR])" = participants)

# gtsummary table  
TrialsInfogt <- tbl_summary(TrialsInfo, 
                            digits = list(all_categorical() ~ c(0, 1)))

TrialsInfogt

# Combine comparisons ####

# list all capitalised abbreviations
Capitalise <- c(
  "bmi" = "BMI",
  "bnp" = "BNP",
  "ccb" = "CCB",
  "cvd" = "CVD",
  "dbp" = "DBP",
  "ef" = "EF",
  "hb" = "Hb",
  "hr" = "HR",
  "ph" = "pH",
  "sbp" = "SBP",
  "bblocker" = "β-blocker",
  "oac" = "OAC",
  "qol" = "QOL score",
  "raasb" = "RAAS inhibitors",
  "hf" = "HF",
  "mace" = "MACE",
  "sae" = "SAEs",
  "af" = "AF",
  "chd" = "CHD",
  "ckd" = "CKD",
  "copd" = "COPD",
  "htn" = "HTN",
  " mi" = " a previous MI",
  "pad" = "PAD",
  "stroke" = "a previous stroke",
  "bleed" = "a previous bleeding event"
)

# Retain comparisons that can be synthesises, including combining measurements
# within the same grouped comparison
CombinedComparisonsPrep <- ComparisonTidy %>%
  # regroup and improve readability
  mutate(
    # note capital C on this column
    Comparison = case_when(
      measurement == "mean"
      ~ paste0("Mean ", comparison),
      (measurement == "median")
      ~ paste0("Median ", comparison),
      mapped == "Coronary procedure"
      ~ "Proportion with a previous coronary procedure",
      subcategory == "Comorbidity"
      ~ paste0("Proportion with ", comparison),
      subcategory == "Preventative"
      ~ paste0("Proportion prescribed ", mapped),
      mapped == "Cardiac"
      ~ "Proportion requiring antiarrhythmia medications",
      mapped == "Vasoactive"
      ~ "Proportion requiring vasoactive agents",
      subcategory == "Treatment"
      ~ paste0("Proportion requiring ", mapped),
      subcategory == "Symptomatic"
      ~ paste0("Proportion prescribed ", mapped),
      mapped == "Death"
      ~ "Proportion who died",
      mapped == "Additional care"
      ~ "Proportion requiring additional care",
      subcategory == "Adverseevent"
      ~ paste0("Proportion who experienced ", mapped),
      subcategory == "Severity"
      ~ "Proportion with more severe disease pre-recruitment",
      mapped == "Smoker"
      ~ "Proportion who smoke",
      subcategory == "Treatmentresponse"
      ~ paste0("Proportion who experienced ", mapped),
      comparison == "Women" ~ "Proportion women",
      subcategory == "Changetotreatment"
      ~ paste0("Proportion who experienced ", mapped),
      mapped == "CVD Death"
      ~ "Proportion of CVD death",
      mapped == "Hospitalisation"
      ~ "Proportion re-hospitalised",
      comparison == "Lives alone"
      ~ "Proportion who live alone",
      mapped == "Ethnicity"
      ~ paste0("Proportion ", comparison),
      comparison == "White race"
      ~ comparison,
      comparison == "Obesity"
      ~ "Proportion with obesity",
      comparison == "Renal function"
      ~ "For removal",
      comparison == "Socioeconomics"
      ~ "Proportion who are less deprived or have a higher educational level",
      comparison == "Length of stay"
      ~ "Nights in hospital",
      comparison == "Survival"
      ~ "Proportion who died",
      TRUE ~ "For removal"
    )
  ) %>%
  # this amalgamates SAEs (MACE, Bleeds, etc)
  mutate(
    Comparison =
      if_else(
        mapped %in% c("Bleed", "MACE", "SAE"),
        "Proportion who experienced serious adverse events",
        Comparison
      )
  ) %>%
  # take out NAs (if there are any) and also filter out any specific ones that
  # aren't useable (manually checked)
  filter(
    !is.na(Comparison),!Comparison %in% c(
      "For removal",
      # Obvious
      "Proportion prescribed Other symptomatic",
      # Other
      "Mean Other haematology",
      # Other
      "Proportion requiring Other treatment",
      # Other
      "Mean Other measure",
      # Other
      "Median Other measure",
      # Other
      "Mean Other other",
      # Other
      "Proportion prescribed Other preventative",
      # Other
      "Proportion with Other comorbidity",
      # Other
      "Proportion with Other resp",
      # Other
      "Mean Coag",
      # confusing grouping
      "Mean Liver function",
      # inconsistent grouping and below power unpacked
      "Mean Treatment response",
      # confusing grouping
      "Mean QoL",
      # inconsistent grouping and below power unpacked
      "Mean Severity score",
      # inconsistent grouping and below power unpacked
      "Mean CVD risk score",
      # inconsistent grouping and below power unpacked
      "Mean Bleed risk score",
      # inconsistent grouping and below power unpacked
      "Proportion with Other CVD",
      # Other
      "Mean Cholesterol subset",
      # inconsistent grouping and below power unpacked
      "Median Cholesterol subset",
      # inconsistent grouping and below power unpacked
      "Proportion who experienced Change to medication" # confusing grouping
    )
  ) %>%
  # sort out grammar
  mutate(Comparison = str_to_sentence(Comparison)) %>%
  mutate(Comparison = str_replace_all(Comparison, Capitalise)) %>%
  # remaining capitalisation issues
  mutate(
    Comparison = str_replace(Comparison, "vasOACtive", "vasoactive"),
    Comparison = str_replace(
      Comparison,
      "prev a previous bleeding event",
      "a previous major bleed"
    )
  ) %>%
  # first address trial-clinical population pair duplication by finding the mean
  # for each mapped comparison (e.g. severity)
  group_by(trial, population, Comparison) %>%
  mutate(
    elvalue = mean(elvalue, na.rm = TRUE),
    inelvalue = mean(inelvalue, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # then remove duplicated values
  distinct(trial, Comparison, elvalue, .keep_all = TRUE)

# Combine further to keep average value for each trial-population pair, where
# trial is duplicated, and keep comparisons present in enough trials/studies
# for synthesis (>1 study, >2 trials)
CombinedComparisons <- CombinedComparisonsPrep %>% 
  # address trial duplication to find means of values where trial has
  # same comparisons more than once in different populations
  group_by(trial, Comparison) %>%
  mutate(
    elvalue = mean(elvalue, na.rm = TRUE),
    inelvalue = mean(inelvalue, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # and remove duplicate values again
  distinct(trial, Comparison, elvalue, .keep_all = TRUE) %>% 
  # now filter for measurements that are present in 2 or more studies and 3 or
  # more trials
  group_by(Comparison) %>%
  mutate(
    study_count = n_distinct(study_author),
    trial_count = n_distinct(trial)
  ) %>%
  ungroup() %>% 
  filter(study_count > 1, trial_count > 2)

# More trial and comparison stats #####

# create factor of trials in order of publication
trialtable <- read_excel("data.xlsx", sheet = "trialorder")

trialtableord <- factor(trialtable$Trial,
                        levels = trialtable$Trial,
                        ordered = TRUE)

# number of populations each trial is examined in
popcount <- Trials0 %>% 
  group_by(trial) %>% 
  summarise(npop = n()) %>%
  ungroup() %>% 
  arrange(match(trial, trialtableord))

# trial-population pair counts
trialpopcompcount <- ComparisonTidy %>% 
  group_by(trial, population) %>% 
  summarise(ntp = n(), .groups = "drop_last") %>% 
  ungroup()

# median number of comparisons (not necessarily synthesised) for all 
# trial-population pairs
median(trialpopcompcount$ntp)

# number of grouped comparisons synthesised for each trial
compcount <- CombinedComparisons %>% 
  group_by(trial) %>% 
  summarise(ncomp = sum(!is.na(study_author))) %>%
  ungroup() %>%
  arrange(match(trial, trialtableord))

allcount <- popcount %>% 
  full_join(compcount, by = "trial") %>% 
  replace_na(list(ncomp = 0))

# median number of comparisons per trial
median(allcount$ncomp)
# range
min(allcount$ncomp)
max(allcount$ncomp)

# Condition group denominators
TotalCVDtrialnum = n_distinct(Trials$trial[Trials$category_2 == "CVD"])
TotalCancertrialnum = n_distinct(Trials$trial[Trials$category_2 == "Cancer"])
TotalResptrialnum = n_distinct(Trials$trial[Trials$category_2 == "Resp"])
TotalMSKtrialnum = n_distinct(Trials$trial[Trials$category_2 == "MSK"])

# stats per study in each condition group
NumericalTrialStats <- All %>% 
  # group by study
  group_by(study_author) %>% 
  # number of trials per study
  summarise(ntrials = length(unique(trial)),
            category = first(category_2)) %>% 
  ungroup() %>% 
  # category being the condition group
  group_by(category) %>% 
  summarise(nstudies = n(),
            meantrials = mean(ntrials),
            mediantrials = median(ntrials),
            singletrials = sum(ntrials == 1),
            percsingle = round(singletrials/nstudies * 100, 1),
            range = paste0(min(ntrials), " - ", max(ntrials))
  )

# Synthesis ####

# compare all comparisons
AllStatsComparisons <- CombinedComparisons %>%
  # use bespoke function above to calculate p-values and convert to table
  comparisonsfun(condition = "All", treatment = "All") %>%
  # round p-values
  mutate(pvalue = round(pvalue, 3)) %>%
  arrange(Comparison) %>%
  rename(
    "Eligible median (IQR)" = Eligible,
    "Ineligible median (IQR)" = Ineligible,
    "p-value" = pvalue,
    "Studies (n)" = nstuds,
    "Populations (n)" = npops,
    "Trials (n)" = ntrials,
    "Comparisons (n)" = ncomparisons
  ) %>%
  mutate(Comparison = str_replace(Comparison, "renal function", "eGFR")) %>%
  # arrange so most comparisons are at the top
  arrange(Category, Subcategory, desc(`Comparisons (n)`)) %>%
  # adjust p-values
  mutate(
    `p-value` = str_replace(`p-value`, "^1$", ">0.999"),
    `p-value` = str_replace(`p-value`, "^0$", "<0.001"),
    # change some problematic wording
    Comparison = str_replace(Comparison, "median a previous bleeding event", 
                             "median bleeding"),
    Comparison = str_replace(
      Comparison,
      "experienced a previous bleeding event",
      "experienced a bleeding event"
    )
  ) %>%
  # total columns
  mutate(
    TotalAlltrial = nrow(Trials),
    TotalCVDtrial = TotalCVDtrialnum,
    TotalCancertrial = TotalCancertrialnum,
    TotalResptrial = TotalResptrialnum,
    TotalMSKtrial = TotalMSKtrialnum
  )


# function for creating better table
process_baseline <- function(data, specialties) {
  df <- data %>%
    # ensure categories kept together
    mutate(Category = factor(
      Category,
      levels = c(
        "Demography",
        "Clinicalcharacteristic",
        "Medication",
        "Outcome"
      )
    )) %>%
    # then order by most comparisons
    arrange(Category, desc(`Trials (n)`)) %>%
    select(
      Comparison,
      `Eligible median (IQR)`,
      `Ineligible median (IQR)`,
      `p-value`,
      `Trials (n)`,
      TotalAlltrial,
      nCVDtrial,
      TotalCVDtrial,
      nCancertrial,
      TotalCancertrial,
      nResptrial,
      TotalResptrial,
      nMSKtrial,
      TotalMSKtrial,
      Condition
    ) %>%
    # rename
    rename(nAlltrial = `Trials (n)`) %>%
    # create new percentage column for each condition group
    mutate(across(
      .cols = all_of(paste0("Total", specialties, "trial")),
      .fns = ~ formatC(get(paste0(
        "n", sub("Total", "", cur_column())
      )) / .x * 100, format = "f", digits = 1),
      .names = "{.col}_pct"
    ))

  # control over which specialties are included in subanalysis and calculate
  # number/pct of comparisons in that specialty
  for (s in specialties) {
    df[[paste0(s, "_summary")]] <- paste0(
      df[[paste0("n", s, "trial")]], " (", 
      df[[paste0("Total", s, "trial_pct")]], "%)")
  }
  
  # create the dataframe
  df %>%
    select(
      Comparison,
      Condition,
      `Eligible median (IQR)`,
      `Ineligible median (IQR)`,
      `p-value`,
      ends_with("_summary")
    )
}



# All trials results ####

# edit this for specialties of interest. For this analysis, using CVD 
# subanalysis only
specialties <- c("All", "CVD")

# create table using above function to calculate number of comaprisons in each
# specialty
ComparisonByTrial <- process_baseline(AllStatsComparisons,
                                      specialties = specialties) %>% 
  mutate(Comparison = str_remove_all(Comparison, 
                                 paste("Proportion with",
                                   "Proportion who experienced",
                                   "Proportion who experience",
                                   "Proportion requiring",
                                   "Proportion prescribed",
                                   "Proportion who are",
                                   "Proportion who",
                                   "Proportion prescribed",
                                   "Proportion requiring",
                                   "Proportion",
                                   "Mean or median",
                                   sep = "|")
                                 )) %>% 
  select(Comparison, All_summary, ends_with("summary"), `Eligible median (IQR)`,
         `Ineligible median (IQR)`, `p-value`)

# Convert data to a flextable
Baselineft <- flextable(ComparisonByTrial)

# Create Word document and add table
Baselinedoc <- read_docx() %>%
body_add_par("Title", style = "heading 1") %>%
body_add_flextable(Baselineft)

# Save Word file
# print(Baselinedoc, target = "BaselineComparisons.docx")

## Multimorbidity ####

Multimorb <- ComparisonTidy %>% 
  filter(mapped == "Multimorbidity")


# Granular comparisons ####

# Extracting all analysed comparisons down to finest granularity extracted

# Table to show all entries grouped together to form each comparison
Granularcomparisons1 <- CombinedComparisons %>%
  # make all lower case to remove duplication through different uses of cases
  mutate(
    nature_of_comparison = str_to_lower(nature_of_comparison),
    Comparison = str_to_lower(Comparison),
    Comparison = str_remove_all(
      Comparison,
      paste(
        "Proportion with ",
        "Proportion who experienced ",
        "Proportion who experience ",
        "Proportion requiring ",
        "Proportion prescribed ",
        "Proportion who are ",
        "Proportion who ",
        "Proportion prescribed ",
        "Proportion requiring ",
        "Proportion ",
        "Mean ",
        "Median ",
        "a previous ",
        sep = "|"
      )
    )
  ) %>%
  group_by(Comparison) %>%
  mutate(ntrials = n_distinct(trial), nstuds = n_distinct(study_author)) %>%
  # move rogue obesity comparison to join the others
  mutate(category.x = if_else(
    Comparison == "obesity", "Clinicalcharacteristic", category.x)) %>%
  ungroup() %>%
  # select granular and grouped characteristic columns as well as broad category
  select(nature_of_comparison,
         Comparison,
         category.x,
         ntrials,
         nstuds,
         trial) %>%
  # for each grouped comparison
  group_by(category.x, Comparison) %>%
  arrange(trial) %>%
    # collect which trials were examined for this comparison
  mutate(trial = paste(unique(trial), collapse = ", ")) %>%
    # keep all distinct granular comparisons
  distinct(nature_of_comparison, .keep_all = TRUE) %>%
    # list them all in one cell
  mutate(nature_of_comparison = paste(
    unique(nature_of_comparison), collapse = ", ")) %>%
  ungroup() %>%
    # deduplicate
  distinct(Comparison, category.x, .keep_all = TRUE) %>%
    # order by category of comparison
  mutate(category.x =
           factor(
             category.x,
             levels = c(
               "Demography",
               "Clinicalcharacteristic",
               "Medication",
               "Outcome"
             )
           )) %>%
  # most trials first
  arrange(category.x, desc(ntrials)) %>%
  # only the ones in the synthesis
  filter(ntrials > 2 & nstuds > 1)

# Convert data to a flextable
Granularft1 <- flextable(Granularcomparisons1)

# Create Word document and add table
Granulardoc1 <- read_docx() %>%
  body_add_par("Title", style = "heading 1") %>%
  body_add_flextable(Granularft1)

# Save Word file
# print(Granulardoc1, target = "GranularComparisons1.docx")


# Visualisation ####

# Proportion ineligible

# get median/IQR of prop ineligible
eligpercentiles <- Trials %>%
  drop_na(prop_ineligible) %>%
  summarise(
    n = n(),
    p25 = quantile(prop_ineligible, 0.25),
    median = quantile(prop_ineligible, 0.50),
    p75 = quantile(prop_ineligible, 0.75)
  )

# number of trials
elign <- eligpercentiles$n

eligplot <- Trials %>%
  # drop NAs and arrange in order
  drop_na(prop_ineligible) %>%
  arrange(prop_ineligible) %>%
  # create column to factorise order of trials
  mutate(trial = fct_inorder(trial)) %>%
  # make graph
  ggplot(aes(x = trial, y = prop_ineligible)) +
  geom_col(aes(fill = category_2), alpha = 0.5) +  # Adjust opacity
  # create lines for 25th, 50th and 75th centile
  # line starts at bottom of y axis and goes to trial number for each centile
  annotate(
    "segment",
    x = 0.5,
    xend = elign / 4,
    y = eligpercentiles$p25,
    yend = eligpercentiles$p25,
    linetype = "dotted",
    color = "black",
    linewidth = 1.2
  ) +
  annotate(
    "segment",
    x = 0.5,
    xend = elign / 2,
    y = eligpercentiles$median,
    yend = eligpercentiles$median,
    linetype = "dotted",
    color = "black",
    linewidth = 1.2
  ) +
  annotate(
    "segment",
    x = 0.5,
    xend = elign / 4 * 3,
    y = eligpercentiles$p75,
    yend = eligpercentiles$p75,
    linetype = "dotted",
    color = "black",
    linewidth = 1.2
  ) +
  # write labels for each
  annotate(
    "label",
    x = elign / 4 - 5,
    y = eligpercentiles$p25 + 11,
    label = paste0(
      "Three quarters of\ntrials would have\nexcluded at least\n",
      round(eligpercentiles$p25, 0),
      "% of patients"
    ),
    color = "black",
    fill = "white",
    size = 5,
    label.padding = unit(0.5, "lines"),
    # rounded corner
    label.r = unit(0.15, "lines")
  ) +
  annotate(
    "label",
    x = elign / 2 - 5,
    y = eligpercentiles$median + 11,
    label = paste0(
      "Half of trials\nwould have\nexcluded at least\n",
      round(eligpercentiles$median, 0),
      "% of patients"
    ),
    color = "black",
    fill = "white",
    size = 5,
    label.padding = unit(0.5, "lines"),
    label.r = unit(0.15, "lines")
  ) +
  annotate(
    "label",
    x = elign / 4 * 3 - 5,
    y = eligpercentiles$p75 + 11,
    label = paste0(
      "One quarter of\ntrials would have\nexcluded at least\n",
      round(eligpercentiles$p75, 0),
      "% of patients"
    ),
    color = "black",
    fill = "white",
    size = 5,
    label.padding = unit(0.5, "lines"),
    label.r = unit(0.15, "lines")
  ) +
  # labels for axes
  labs(x = "Trials, in order of proportion ineligible", y = "Percentage of comparison population ineligible", fill = "Category") +
  # Set x-axis breaks to intervals of 10
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    expand = c(0, 0)
  ) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "CVD" = "#D32F2F",
    # red
    "Cancer" = "#F48FB1",
    # pink
    "Resp" = "#1976D2",
    # darker blue (more distinct from purple)
    "MSK" = "#FFB74D"       # light pink (lighter to contrast with red)
  )) +
  theme(
    axis.text.x = element_blank(),
    # Hides x-axis text for cleaner appearance
    panel.grid.major = element_blank(),
    # Removes major grid lines
    panel.grid.minor = element_blank(),
    # Removes minor grid lines
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  )

eligplot

# ggsave("eligplot.svg", eligplot,
       # height = 250, width = 400, unit = "mm", dpi = 300)


# Violin plots for comparisons

# create long df: 2 rows for each
longCombinedComparisons <- CombinedComparisons %>%
  mutate(
    domain = case_when(
      subcategory %in% c("Comorbidity", "Measurement")
      ~ subcategory,
      subcategory == "Riskfactor"
      ~ "Demography",
      subcategory == "Severity"
      ~ "Comorbidity",
      .default = category.x
    )
  ) %>%
  select(Comparison,
         comparison,
         elvalue,
         inelvalue,
         domain,
         measurement,
         trial) %>%
  # 2 rows, one eligible, one ineligible
  pivot_longer(
    cols = c(elvalue, inelvalue),
    names_to = "Eligibility",
    values_to = "value"
  ) %>%
  # tinker
  mutate(
    # reword
    Eligibility = if_else(Eligibility == "elvalue", "eligible", "ineligible"),
    # group proportions and averages
    xmeas = if_else(measurement == "prop", "prop", "ave"),
    # make cleaner
    Comparison = str_remove_all(
      Comparison,
      paste(
        "Proportion with ",
        "Proportion who experienced ",
        "Proportion who experience ",
        "Proportion requiring ",
        "Proportion prescribed ",
        "Proportion who are ",
        "Proportion who ",
        "Proportion prescribed ",
        "Proportion requiring ",
        "Proportion ",
        "a previous ",
        " pre-recruitment",
        " or have a higher educational level",
        " race",
        sep = "|"
      )
    ),
    # If no capitals (e.g., abbrevs or Mean/Median) then change to sentence case
    Comparison = if_else(
      !str_detect(Comparison, "[A-Z]"),
      str_to_sentence(Comparison),
      Comparison
    ),
    # add units
    Comparison = case_when(
      xmeas == "prop" ~
        paste0(Comparison, " (%)"),
      comparison == "Age" ~
        paste0(Comparison, " (years)"),
      comparison == "BMI" ~
        paste0(Comparison, " (kg/m²)"),
      comparison == "Renal function" ~
        paste0(Comparison, " (ml/min/1.73m²"),
      comparison == "SBP" ~
        paste0(Comparison, " (mmHg)"),
      comparison == "DBP" ~
        paste0(Comparison, " (mmHg)"),
      comparison == "HR" ~
        paste0(Comparison, " (bpm)"),
      comparison == "EF" ~
        paste0(Comparison, " (%)"),
      comparison == "Hb" ~
        paste0(Comparison, " (g/l)"),
      comparison == "BNP" ~
        paste0(Comparison, " (pg/ml)"),
      comparison == "Total cholesterol" ~
        paste0(Comparison, " (mmol/l)"),
      comparison == "Lactate" ~
        paste0(Comparison, " (mmol/l)"),
      comparison == "Glucose" ~
        paste0(Comparison, " (mg/dl)"),
      comparison == "pH" ~
        Comparison,
      comparison == "pH" ~
        Comparison,
      Comparison == "Median survival (months)" ~
        Comparison
    ),
    Comparison = str_replace(Comparison, "renal function", "eGFR")
  ) %>%
  # group by lower case comparison to keep measures together (e.g., mean and
  # median age)
  group_by(comparison) %>%
  # add a column with group size to keep above together
  mutate(group_size = n()) %>%
  ungroup() %>%
  # put in the same order as the table (although not necessary with code below)
  mutate(domain = factor(
    domain,
    levels = c(
      "Demography",
      "Comorbidity",
      "Medication",
      "Measurement",
      "Outcome"
    )
  )) %>%
  # arrange by group size and factorise
  arrange(domain, desc(group_size)) %>%
  mutate(Comparison = factor(Comparison, levels = unique(Comparison)))

# list of categories of comparisons
domains <- unique(longCombinedComparisons$domain)

## All comparisons in violin plots ####
walk(domains, function(domain_name) {
  filtered_data <- longCombinedComparisons %>%
    filter(domain == domain_name) %>%
    mutate(Comparison = factor(Comparison)) # Ensure consistent facet order
  
  # Calculate p-values for each comparison
  p_valuesloop <- filtered_data %>%
    # calculate p-values for each comparison
    pivot_wider(names_from = Eligibility, values_from = value) %>% 
    group_by(Comparison) %>%
    summarise(
    pvalue = wilcox.test(eligible, ineligible, paired = TRUE, # this bit makes it signed rank
                         exact = FALSE)$p.value,
    .groups = 'drop'
  ) %>%
  # change notation for small p-values
  mutate(
    p_label = case_when(
      pvalue < 0.001 ~ "p < 0.001",
      TRUE ~ sprintf("p = %.3f", pvalue)
    )
  )
  
  
  filtered_data <- filtered_data %>%
    # join p-value labels onto data
    left_join(p_valuesloop, by = "Comparison") %>%
    # factor column to add p-value
    mutate(Comparison_with_p = fct_inorder(paste0(
      Comparison, " (", p_label, ")")))
  
  # ordered list of comparison labels for faceting
  comparisons <- levels(filtered_data$Comparison_with_p)
  n_facets <- length(comparisons)
  ncol <- 3
  nrow <- ceiling(n_facets / ncol)
  
  # output dimensions
  width <- ncol * 4   # 4 inches per column
  height <- nrow * 3  # 3 inches per row
  
  # Create list of y-axis scales per facet (prop and average)
  y_limits_list2 <- map(comparisons, function(comp) {
    xmeas_type <- filtered_data %>%
      filter(Comparison_with_p == comp) %>%
      distinct(xmeas) %>%
      pull()
    
    if (xmeas_type == "prop") {
      as.formula(paste0("Comparison_with_p == '", comp, 
                        "' ~ scale_y_continuous(limits = c(0, 100))"))
    } else {
      as.formula(paste0("Comparison_with_p == '", comp, 
                        "' ~ scale_y_continuous()"))
    }
  })
 
  # Build violin and boxplot faceted by comparison group 
  p <- ggplot(filtered_data, aes(
    x = Eligibility, y = value, fill = Eligibility)) +
    # violin
    geom_violin(width = 0.8, alpha = 0.6, trim = FALSE) +
    # box
    geom_boxplot(width = 0.1, show.legend = FALSE, 
                 fill = "white", colour = "black") +
    # facet by comparison with free y-axis
    ggh4x::facet_wrap2(~ Comparison_with_p, ncol = ncol, scales = "free_y") +
    # apply limits specified above (props 0-100)
    ggh4x::facetted_pos_scales(y = y_limits_list2) +
    labs(title = domain_name, x = "", y = "") +
    theme_minimal() +
    scale_fill_brewer(palette = "Dark2") +
    theme(
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      strip.text = element_text(size = 12)
    )
  
  p
  
  # save each individually
  svg_filename <- paste0("plot_", domain_name, ".svg")
  svg(svg_filename, width = width, height = height)
  print(p)
  dev.off()
})

## Filtered violin plots #####

shortCombinedComparisons <- longCombinedComparisons %>% 
  # create column for number of comparisons
  # tolower() as had an issue with bBlocker/Bblockers
  group_by(Comparison) %>% # more granular (e.g. mean and median separate)
  mutate(n = n()) %>% 
  ungroup %>% 
  group_by(tolower(comparison)) %>% # combined, to get medians where mean is in
  mutate(max_n = max(n)) %>% 
  ungroup() %>% 
  filter(max_n >= 40, domain != "Outcome") %>% 
  mutate(domain = factor(domain, levels = c(
    "Demography", "Comorbidity", "Measurement", "Medication"
  ))) %>% 
  arrange(domain, desc(max_n)) %>% 
  mutate(Comparison = factor(Comparison), # Ensure consistent facet order
         Comparison = str_replace(
           Comparison, "age", "age (years)"),
         Comparison = str_replace(
           Comparison, "BMI", "BMI (kg/m²)"),
         Comparison = str_replace(
           Comparison, "eGFR", "eGFR (ml/min/1.73m²)")
  )

# Calculate p-values for each comparison
p_values <- shortCombinedComparisons %>%
  # calculate p-values for each comparison
  pivot_wider(names_from = Eligibility, values_from = value) %>% 
  group_by(Comparison) %>%
  summarise(
    pvalue = wilcox.test(eligible, ineligible, paired = TRUE, # this bit makes 
                         # it signed rank
                         exact = FALSE)$p.value,
    .groups = 'drop'
  ) %>%
  # change notation for small p-values
  mutate(
    p_label = case_when(
      pvalue < 0.001 ~ "p < 0.001",
      TRUE ~ sprintf("p = %.3f", pvalue)
    )
  )

# make the label include the p-value
shortCombinedComparisons <- shortCombinedComparisons %>%
  left_join(p_values, by = "Comparison") %>%
  mutate(Comparison_with_p = fct_inorder(paste0(Comparison, 
                                                " (", p_label, ")")))

comparisons <- unique(shortCombinedComparisons$Comparison_with_p)

# Create list of y-axis scales per facet (different for averages and props)
y_limits_list2 <- map(comparisons, function(comp) {
  xmeas_type <- shortCombinedComparisons %>%
    filter(Comparison_with_p == comp) %>%
    distinct(xmeas) %>%
    pull()
  
  if (xmeas_type == "prop") {
    as.formula(paste0("Comparison_with_p == '", comp, 
                      "' ~ scale_y_continuous(limits = c(0, 100))"))
  } else {
    as.formula(paste0("Comparison_with_p == '", comp, 
                      "' ~ scale_y_continuous()"))
  }
})

# plot for comparisons in >=20 trials
p <- ggplot(shortCombinedComparisons, aes(
  x = Eligibility, y = value, fill = Eligibility)) +
  geom_violin(width = 0.8, alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, show.legend = FALSE, 
               fill = "white", colour = "black") +
  ggh4x::facet_wrap2(~ Comparison_with_p, ncol = 4, scales = "free_y") +
  # use list above to customise y-axis scales
  ggh4x::facetted_pos_scales(y = y_limits_list2) +
  labs(x = "", y = "") +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2") +
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text = element_text(size = 12,)
  )
p

# ggsave("violinplot-mostcommon.svg", plot = p,
      # height = 400, width = 400, unit = "mm")
