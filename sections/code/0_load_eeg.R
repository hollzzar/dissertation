####################
## Set up session ##
####################

# Load packages and variables
source("0_global.R")

# Save file path
eeg_path <- "/Users/hollyzaharchuk/Mirror/dissertation/13_eeg"

######################
## Set up contrasts ##
######################

# Create contrast codes: Match type
contrast_match <- matrix(c(1/2, -1/2, 0,
                           -1/3, -1/3, 2/3), 
                         ncol = 2,
                         dimnames = list(c("control", "competitor", "identity"), 
                                         c("cont_comp", "mis_match")))

# Create contrast codes: Target type
contrast_target <- matrix(c(1/2, -1/2), 
                          ncol = 1,
                          dimnames = list(c("test", "filler"), 
                                          c("test_fill")))

#######################
## Load offline data ##
#######################

# Load language background (keep separate)
eeg_lhq <- read.csv(paste(eeg_path, "output/clean_lhq.csv", sep = "/")) %>%
  replace_na(list(lang_eng_dialect = "Not provided"))

# Get contrast code for dialect experience
lhq_dat <- eeg_lhq %>%
  mutate(eng_cc = case_when(
    care_eng_l1 == "Yes" ~ -1/2,
    care_eng_l1 == "No, neither one" ~ 1/2,
    TRUE ~ 0)) %>%
  select(id, eng_cc)

# Load offline measures
behave_dat <- read.csv(paste(eeg_path, "output/clean_axcpt.csv", sep = "/")) %>%
  left_join(read.csv(paste(eeg_path, "output/clean_offline.csv", sep = "/")), by = "id") %>%
  left_join(lhq_dat, by = "id") %>%
  replace_na(list(spk_fluency_2 = 0))

##########################
## Load behavioral data ##
##########################

# Load table with trigger and condition information
code_tab <- read.csv(paste(eeg_path, "input/code_tab.csv", sep = "/")) %>%
  mutate(prime_cond = factor(prime_cond, 
                             levels = c("identity", "competitor", "control"),
                             labels = c("identity", "competitor", "unrelated")))

# Load stimulus info
word_info <- read.csv("~/Mirror/dissertation/4_stims/output/selections/all_behave.csv") 

# Load in behavioral data
all_trial_dat <- read.csv(paste(eeg_path, "output/clean_eprime.csv", sep = "/")) %>% 
  left_join(code_tab) %>%
  mutate(vot = T2OnsetTime - T1OnsetTime,
         target = gsub("\\.wav", "", Target)) %>%
  rename(prime = Prime) %>%
  left_join(word_info, by = c("target" = "stim")) %>%
  select(trial, id, prime, target, true_rt, acc, prime_cond, target_cond, vot, target, FreqZipfUS) %>%
  group_by(id) %>% 
  mutate(exp_amount = cumsum(target_cond == "exposure")) %>%
  ungroup()

# Add contrast codes
contrasts(all_trial_dat$prime_cond) <- contrast_match

##########################
## Behavioral: accuracy ##
##########################

# Pull test and filler trials
behave_trials_all <- all_trial_dat %>%
  dplyr::filter(target_cond == "test" | target_cond == "filler")

# Clean
behave_trials_clean <- behave_trials_all %>%
  dplyr::filter(acc == 1 | true_rt > 250) %>%
  mutate(across(.cols = c(id, target, target_cond, prime),
                as.factor),
         exp_cent = exp_amount - mean(exp_amount),
         vot_cent = vot - mean(vot),
         freq_cent = FreqZipfUS - mean(FreqZipfUS))

# Set contrasts
contrasts(behave_trials_clean$target_cond) <- contrast_target

######################
## Behavioral: test ##
######################

# Pull test trials
test_trial_dat <- all_trial_dat %>%
  dplyr::filter(target_cond == "test")

# Prep data
test_trial_dat_clean <- test_trial_dat %>%
  dplyr::filter(acc == 1 | true_rt > 250) %>%
  mutate(across(.cols = c(id, target, prime),
                as.factor),
         trial_cent = trial - mean(trial),
         vot_cent = vot - mean(vot),
         freq_cent = FreqZipfUS - mean(FreqZipfUS),
         exp_cent = exp_amount - mean(exp_amount))

###################
## Load ERP data ##
###################

# Load in N1
erp_n1_all <- read.csv(paste(eeg_path, "output/erp_n1.csv", sep = "/")) %>% 
  rename(prime = Prime,
         target = Target) %>%
  mutate(target = gsub("\\.wav", "", target)) %>%
  left_join(all_trial_dat, by = join_by(id, trial, prime, target))

# Load in P2
erp_p2_all <- read.csv(paste(eeg_path, "output/erp_p2.csv", sep = "/")) %>% 
  rename(prime = Prime,
         target = Target) %>%
  mutate(target = gsub("\\.wav", "", target)) %>%
  left_join(all_trial_dat, by = join_by(id, trial, prime, target))

# Load in N400
erp_n4_all <- read.csv(paste(eeg_path, "output/erp_n4.csv", sep = "/")) %>% 
  rename(prime = Prime,
         target = Target) %>%
  mutate(target = gsub("\\.wav", "", target)) %>%
  left_join(all_trial_dat, by = join_by(id, trial, prime, target))

##################
## Get channels ##
##################

# Pivot
erp_n1_long <- erp_n1_all %>%
  select(-T1Label) %>%
  pivot_longer(F7:P8, names_to = "chan", values_to = "volt")

# N1
chans_n1 <- erp_n1_long %>%
  group_by(chan) %>%
  summarise(mean_eff = mean(volt)) %>%
  dplyr::filter(mean_eff < -1) %>%
  pull(chan)

# Pivot
erp_p2_long <- erp_p2_all %>%
  select(-T1Label) %>%
  pivot_longer(F7:P8, names_to = "chan", values_to = "volt")

# P2
chans_p2 <- erp_p2_long %>%
  group_by(chan) %>%
  summarise(mean_eff = mean(volt)) %>%
  dplyr::filter(mean_eff > 1) %>%
  pull(chan)

# Pivot
erp_n4_long <- erp_n4_all %>%
  select(-T1Label) %>%
  pivot_longer(F7:P8, names_to = "chan", values_to = "volt")

# N400
chans_n4 <- erp_n4_long %>%
  group_by(chan) %>%
  summarise(mean_eff = mean(volt)) %>%
  dplyr::filter(mean_eff < -1) %>%
  pull(chan)

###############
## Prep data ##
###############

# N1
erp_n1 <- erp_n1_long %>%
  dplyr::filter(chan %in% chans_n1) %>%
  mutate(across(.cols = c(chan, id, target, prime),
                as.factor),
         trial_cent = trial - mean(trial),
         vot_cent = vot - mean(vot),
         freq_cent = FreqZipfUS - mean(FreqZipfUS),
         exp_cent = exp_amount - mean(exp_amount))

# P2
erp_p2 <- erp_p2_long %>%
  dplyr::filter(chan %in% chans_p2) %>%
  mutate(across(.cols = c(chan, id, target, prime),
                as.factor),
         trial_cent = trial - mean(trial),
         vot_cent = vot - mean(vot),
         freq_cent = FreqZipfUS - mean(FreqZipfUS),
         exp_cent = exp_amount - mean(exp_amount))

# N400
erp_n4 <- erp_n4_long %>%
  dplyr::filter(chan %in% chans_n4) %>%
  mutate(across(.cols = c(chan, id, target, prime),
                as.factor),
         trial_cent = trial - mean(trial),
         vot_cent = vot - mean(vot),
         freq_cent = FreqZipfUS - mean(FreqZipfUS),
         exp_cent = exp_amount - mean(exp_amount))

##################
## Correlations ##
##################

# Mean N1
mean_n1 <- erp_n1 %>%
  group_by(id, prime_cond) %>%
  summarise(mean = mean(volt), .groups = "keep") %>%
  pivot_wider(id_cols = id, names_from = prime_cond, values_from = mean) %>%
  mutate(mismatch = mean(competitor, unrelated),
         mean_n1 = abs(identity - mismatch))

# Mean P2
mean_p2 <- erp_p2 %>%
  group_by(id, prime_cond) %>%
  summarise(mean = mean(volt), .groups = "keep") %>%
  pivot_wider(id_cols = id, names_from = prime_cond, values_from = mean) %>%
  mutate(mismatch = mean(competitor, unrelated),
         mean_p2 = abs(identity - mismatch))

# Mean N400
mean_n4 <- erp_n4 %>%
  group_by(id, prime_cond) %>%
  summarise(mean = mean(volt), .groups = "keep") %>%
  pivot_wider(id_cols = id, names_from = prime_cond, values_from = mean) %>%
  mutate(mismatch = mean(competitor, unrelated),
         mean_n4 = abs(identity - mismatch))

# Accuracy
mean_acc <- test_trial_dat_clean %>%
  group_by(id, prime_cond) %>%
  summarise(mean = mean(acc), .groups = "keep") %>%
  pivot_wider(id_cols = id, names_from = prime_cond, values_from = mean) %>%
  mutate(mismatch = mean(competitor, unrelated),
         mean_acc = identity - mismatch)

# Combine
corr_dat <- left_join(behave_dat %>% mutate(id = as.character(id)), 
                      mean_n1 %>% select(id, mean_n1), by = "id") %>%
  left_join(mean_p2 %>% select(id, mean_p2), by = "id") %>%
  left_join(mean_n4 %>% select(id, mean_n4), by = "id") %>%
  left_join(mean_acc %>% select(id, mean_acc), by = "id") %>%
  select(-c(acc_l1_rom, acc_region, id))



