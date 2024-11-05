library(tidyverse)
library(stringr)
library(reticulate)
library(viridis)
library(epiR)
library(progress)
np = import('numpy')
pd = import('pandas')


########################################
##  READ IN RESULTS FROM SIMULATIONS  ##
########################################


# Read in all parameter sets
data =
  read.csv('simulations/parameters.csv') %>%
  as_tibble()
#data_set2 = data %>% filter(symptoms_ural_male > 0.50 & symptoms_ural_male > symptoms_ural_female & symptoms_ural_female < 0.5)
indices = which(data$symptoms_ural_male>0.5 & data$symptoms_ural_male>data$symptoms_ural_female & data$symptoms_ural_female<0.5)
data_set2 <- data[indices, ]

# Read in the target prevalence
target =
  read.csv('data/calibration_prevalence.csv') %>%
  as_tibble()


# Read in the target testing rates
target_test =
  read.csv('data/testing_rates.csv') %>%
  as_tibble() %>%
  mutate(prob = 100 * prob,
         lab = case_when(gender == 0 ~ str_c('f', age_group, sep = ''),
                         gender == 1 ~ str_c('m', age_group, sep = ''))) %>%
  select(-age_group, -gender) %>%
  spread(lab, prob) %>%
  mutate(tot = 19,
         m = 14,
         f = 23)


# Initilise columns for computing prevalence
data_set2 =
  data_set2 %>%
  mutate(set = indices,
         prev_tot = NA,
         prev_m = NA,
         prev_f = NA,
         prev_m0 = NA,
         prev_m1 = NA,
         prev_m2 = NA,
         prev_m3 = NA,
         prev_f0 = NA,
         prev_f1 = NA,
         prev_f2 = NA,
         prev_f3 = NA,
         test_tot = NA,
         test_m = NA,
         test_f = NA,
         not_equlib = TRUE)


# Read in data and determine the equilibrium point
scenario = 3
pb = progress_bar$new(total = nrow(data_set2))
for (i in indices-1){
  pb$tick()


  # Load simulation data
  file_name = str_c('simulations/calibration/scenario_', scenario, '/simulation_', i, '_output_prevalence.npy')
  file_name_2 = str_c('simulations/calibration/scenario_', scenario, '/simulation_', i, '_output_environment.pkl')
  if ( file.exists(file_name) ){


    # Load prevalence
    prev = np$load(file_name) %>% as_tibble()
    names(prev) = c('tot', 'm', 'f', 'm0', 'm1', 'm2', 'm3', 'm4', 'f0', 'f1', 'f2', 'f3', 'f4')
    prev = 100 * prev


    # Load testing rates
    test = 100 * pd$read_pickle(file_name_2)$tstt


    # Look for the maximum t where the prevalence seems to have converged to equilibrium
    prev = prev[seq(round(nrow(prev)/2), nrow(prev), 50),]
    test = test[seq(round(nrow(test)/2), nrow(test), 50),]
    n = nrow(prev)
    prev$t = 1:n
    for ( t0 in 1:(n-10) ){


      # Fit a linear regression
      input_data = prev %>% mutate(t = t - t0)
      fit = lm('tot ~ 1 + t', data = input_data[t0:n, ])


      # Check to see if 0 is in the CI for the slope
      ci = confint(fit, 't', level = 0.95)
      if ( ci[1] < 0 ){


        # Compute the mean prevalence
        prev = prev[t0:n, ]
        test = test[t0:n, c(1, 5, 9)]
        data_set2$prev_tot[i+1] = mean(prev$tot)
        data_set2$prev_m[i+1] = mean(prev$m)
        data_set2$prev_f[i+1] = mean(prev$f)
        data_set2$prev_m0[i+1] = mean(prev$m0)
        data_set2$prev_m1[i+1] = mean(prev$m1)
        data_set2$prev_m2[i+1] = mean(prev$m2)
        data_set2$prev_m3[i+1] = mean(prev$m3)
        data_set2$prev_f0[i+1] = mean(prev$f0)
        data_set2$prev_f1[i+1] = mean(prev$f1)
        data_set2$prev_f2[i+1] = mean(prev$f2)
        data_set2$prev_f3[i+1] = mean(prev$f3)
        data_set2$test_tot[i+1] = mean(test[,1])
        data_set2$test_m[i+1] = mean(test[,2])
        data_set2$test_f[i+1] = mean(test[,3])
        break

      }
    }
  }
}


# Drop incomplete runs
cat(str_c(sum(is.na(data_set2$prev_f)), ' (', 100*sum(is.na(data_set2$prev_f))/nrow(data_set2), '%) runs not finished'))
data_set2 =
  data_set2 %>%
  filter(!is.na(prev_tot))


# Summary of overall prevalence
data_set2 %>%
  select(starts_with('prev') | starts_with('test')) %>%
  gather(series, prev) %>%
  ggplot(aes(x=prev)) +
  geom_histogram(bins = 20) +
  facet_wrap(~series) +
  labs(x = 'Prevalence',
       y = 'Number of Simualtions',
       title = 'Distribution of Prevalence Amongst Sub-populations')


############################
##  SENSITIVITY ANALYSIS  ##
############################


# Compute correlation with different variables
x = data_set2 %>% filter(!is.na(prev_tot)) %>% select(-starts_with('prev'), -starts_with('test'), -starts_with('mean'), -set)
y = data_set2 %>% filter(!is.na(prev_tot)) %>% select(starts_with('prev'), starts_with('test'))
out = tibble(var = names(x),
             cor_tot = NA,
             cor_m = NA,
             cor_f = NA,
             cor_m0 = NA,
             cor_m1 = NA,
             cor_m2 = NA,
             cor_m3 = NA,
             cor_f0 = NA,
             cor_f1 = NA,
             cor_f2 = NA,
             cor_f3 = NA,
             cor_ttot = NA,
             cor_tm = NA,
             cor_tf = NA)
compute_cor = function(x, y){
  out = epi.prcc(tibble(x, y))
  return(out$est)
}
for ( i in 1:ncol(x) ){
  if( !all(is.na(x[,i])) ){
    #x[,i] = x[,i] - mean(data.matrix(x[,i]), na.rm = T)
    #x[,i] = x[,i]/sqrt(var(data.matrix(x[,i]), na.rm = T))
    out$cor_tot[i] = compute_cor(x[,i], y[,1])
    out$cor_m[i] = compute_cor(x[,i], y[,2])
    out$cor_f[i] = compute_cor(x[,i], y[,3])
    out$cor_m0[i] = compute_cor(x[,i], y[,4])
    out$cor_m1[i] = compute_cor(x[,i], y[,5])
    out$cor_m2[i] = compute_cor(x[,i], y[,6])
    out$cor_m3[i] = compute_cor(x[,i], y[,7])
    out$cor_f0[i] = compute_cor(x[,i], y[,8])
    out$cor_f1[i] = compute_cor(x[,i], y[,9])
    out$cor_f2[i] = compute_cor(x[,i], y[,10])
    out$cor_f3[i] = compute_cor(x[,i], y[,11])
    out$cor_ttot[i] = compute_cor(x[,i], y[,12])
    out$cor_tm[i] = compute_cor(x[,i], y[,13])
    out$cor_tf[i] = compute_cor(x[,i], y[,14])
  }
}


# Plot the correlation for each prevalence type
out %>%
  mutate(var = case_when(var == 'pup' ~ 'trans_prob_ural_phar',
                         var == 'ppp' ~ 'trans_prob_phar_phar',
                         var == 'pru' ~ 'trans_prob_rect_ural',
                         var == 'ppu' ~ 'trans_prob_phar_ural',
                         var == 'pur' ~ 'trans_prob_ural_rect',
                         var == 'prp' ~ 'trans_prob_rect_phar',
                         var == 'ppr' ~ 'trans_prob_phar_rect',
                         var == 'puu' ~ 'trans_prob_ural_ural',
                         var == 'p_anal_MM' ~ 'prob_anal_M2M',
                         var == 'p_sex_MM' ~ 'prob_sex_M2M',
                         var == 'p_rim' ~ 'prob_rim',
                         var == 'p_kiss' ~ 'prob_kiss',
                         var == 'p_oral_MM' ~ 'prob_oral_M2M',
                         var == 'p_oral_MF' ~ 'prob_oral_M2F',
                         var == 'p_oral_FM' ~ 'prob_oral_F2M',
                         var == 'p_oral_FF' ~ 'prob_oral_F2F',
                         var == 'p_sex_MF' ~ 'prob_sex_M2F',
                         var == 'p_anal_MF' ~ 'prob_anal_M2F',
                         var == 'p_sex_FF' ~ 'prob_sex_F2F',
                         T ~ var),
         var = reorder(var, cor_tot)) %>%
  gather(type, val, -var) %>%
  mutate(type = ordered(type,
                        levels=c('cor_tot', 'cor_m', 'cor_f', 'cor_m0', 'cor_m1', 'cor_m2', 'cor_m3', 'cor_f0', 'cor_f1', 'cor_f2', 'cor_f3', 'cor_ttot', 'cor_tm', 'cor_tf'),
                        labels=c('Overall', 'Males', 'Females', 'Male 16-19', 'Male 20-24', 'Male 25-29', 'Male 30-35', 'Female 16-19', 'Female 20-24', 'Female 25-29', 'Female 30-35', 'Test Overall', 'Test Males', 'Test Females')),
         # val = case_when(val > 0.2 ~ 0.2,
         #                 val > 0.1 ~ 0.1,
         #                 val > 0 ~ 0,
         #                 val > -0.1 ~ -0.1,
         #                 val > -0.2 ~ -0.2,
         #                 val > -0.3 ~ -0.3,
         #                 val > -0.4 ~ -0.4,
         #                 val > -0.5 ~ -0.5,
         #                 val > -0.6 ~ -0.6,
         #                 val > -0.7 ~ -0.7)
         ) %>%
  filter(!is.na(val)) %>%
  ggplot(aes(x=type, y=var, fill=val)) +
  geom_tile() +
  # scale_fill_viridis() +
  scale_fill_distiller(palette = 'Spectral') + #, limits = c(-.9, .9)) +
  labs(y = 'Parameters being calibrated',
       x = 'Prevalence Category',
       title = 'Partial-rank Correlation Coefficient Between Simulation Results and Input Parameters',
       fill = 'Correlation') +
  theme(axis.text.x = element_text(angle = -30, hjust = 0.1))


################################
##  FIND THE BEST PARAMETERS  ##
################################


# Compute how good each parameter set is
overall_weight = 1/11
gender_weight = 1/11
age_gender_weight_m = 1/11
age_gender_weight_f = 1/11
data_set2 =
  data_set2 %>%
  mutate(prev_ss = (overall_weight * (prev_tot - target$tot)^2 +
                    gender_weight * (prev_m - target$m)^2 +
                    gender_weight * (prev_f - target$f)^2 +
                    age_gender_weight_m * (prev_m0 - target$m0)^2 +
                    age_gender_weight_m * (prev_m1 - target$m1)^2 +
                    age_gender_weight_m * (prev_m2 - target$m2)^2 +
                    age_gender_weight_m * (prev_m3 - target$m3)^2 +
                    age_gender_weight_f * (prev_f0 - target$f0)^2 +
                    age_gender_weight_f * (prev_f1 - target$f1)^2 +
                    age_gender_weight_f * (prev_f2 - target$f2)^2 +
                    age_gender_weight_f * (prev_f3 - target$f3)^2),
         test_ss = (1/6 * (prev_tot - target$tot)^2 +
                    1/6 * (prev_m - target$m)^2 +
                    1/6 * (prev_f - target$f)^2 +
                    1/6 * (test_tot - target_test$tot)^2 +
                    1/6 * (test_m - target_test$m)^2 +
                    1/6 * (test_f - target_test$f)^2))


# Make a graph of all RSS results
data_set2 %>%
  filter(!is.na(test_ss)) %>%
  ggplot(aes(x=test_ss)) +
  geom_histogram(bins = 40) +
  labs(x = 'Mean Residual Sum of Squares',
       y = 'Count',
       title = 'Distribution of the Goodness of Fit for all Simulations')


# Pull out the best 50
calibrated =
  data %>%
  arrange(test_ss) %>%
  slice_head(n = 50)
write.csv(calibrated, str_c('simulations/calibrated_scenario_', scenario, '.csv'), row.names = F)


###################################################
##  DISTRIBUTION OF THE BEST PREVALENCE RESULTS  ##
###################################################
# Look at the distribution of prevalence for the top 50


# Preallocate
n = nrow(prev)
n = 19
prev_overall = array(0, c(n, 50))
prev_m = array(0, c(n, 50))
prev_f = array(0, c(n, 50))


# Read in data and determine the equilibrium point
# scenario = 1
pb = progress_bar$new(total = nrow(calibrated))
for (j in 1:50){
  pb$tick()
  i = calibrated$set[j]


  # Load prevalence data
  file_name = str_c('simulations/calibration/scenario_', scenario, '/simulation_', i, '_output_prevalence.npy')
  if ( file.exists(file_name) ){
    prev = np$load(file_name) %>% as_tibble()
    names(prev) = c('tot', 'm', 'f', 'm0', 'm1', 'm2', 'm3', 'm4', 'f0', 'f1', 'f2', 'f3', 'f4')
    prev = prev[round(nrow(prev)/2):nrow(prev),]
    prev = prev[seq(2, nrow(prev), 50),]
    prev = 100 * prev


    # Compute the average prevalence
    prev_overall[,j] = prev$tot
    prev_m[,j] = prev$m
    prev_f[,j] = prev$f
  }
}


# Plot Overall
row_mean = rowMeans(prev_overall)
prev_overall %>%
  as_tibble() %>%
  mutate(t = 1:nrow(prev_overall),
         row_mean = row_mean) %>%
  gather(sim, val, -t, -row_mean) %>%
  mutate(sim = as.factor(str_extract(sim, '[0-9]+'))) %>%
  ggplot(aes(x=t)) +
  geom_line(aes(y=val, group=sim), alpha = 0.1) +
  geom_line(aes(y=row_mean), colour = 'red') +
  geom_hline(aes(yintercept = target$tot), lty = 'dashed') +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = 'Time step (every 50)',
       y = 'Prevalence',
       title = 'Comparison Rolling Mean of Overall Simulated Prevalence to Overall STRIVE Prevalence')


# Plot Male
row_mean = rowMeans(prev_m)
prev_m %>%
  as_tibble() %>%
  mutate(t = 1:nrow(prev_m),
         row_mean = row_mean) %>%
  gather(sim, val, -t, -row_mean) %>%
  mutate(sim = as.factor(str_extract(sim, '[0-9]+'))) %>%
  ggplot(aes(x=t)) +
  geom_line(aes(y=val, group=sim), alpha = 0.1) +
  geom_line(aes(y=row_mean), colour = 'red') +
  geom_hline(aes(yintercept = target$m), lty = 'dashed') +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = 'Days from 2000th time step',
       y = 'Prevalence',
       title = 'Comparison Rolling Mean of Male Simulated Prevalence to Male STRIVE Prevalence')


# Plot Female
row_mean = rowMeans(prev_f)
col_mean = mean(colMeans(prev_f))
prev_f %>%
  as_tibble() %>%
  mutate(t = 1:nrow(prev_f),
         row_mean = row_mean) %>%
  gather(sim, val, -t, -row_mean) %>%
  mutate(sim = as.factor(str_extract(sim, '[0-9]+'))) %>%
  ggplot(aes(x=t)) +
  geom_line(aes(y=val, group=sim), alpha = 0.1) +
  geom_line(aes(y=row_mean), colour = 'red') +
  geom_hline(aes(yintercept = target$f), lty = 'dashed') +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = 'Days from 2000th time step',
       y = 'Prevalence',
       title = 'Comparison Rolling Mean of Female Simulated Prevalence to Female STRIVE Prevalence')


##############################################
##  PLOT THE DISTRIBUTION OF THE BEST FITS  ##
##############################################


# Compute how good each parameter set is
data_set2 %>%
  filter(!is.na(prev_tot)) %>%
  mutate(prev_tot = prev_tot - target$tot,
         prev_m = prev_m - target$m,
         prev_f = prev_f - target$f,
         prev_m0 = prev_m0 - target$m0,
         prev_m1 = prev_m1 - target$m1,
         prev_m2 = prev_m2 - target$m2,
         prev_m3 = prev_m3 - target$m3,
         prev_f0 = prev_f0 - target$f0,
         prev_f1 = prev_f1 - target$f1,
         prev_f2 = prev_f2 - target$f2,
         prev_f3 = prev_f3 - target$f3,
         test_tot = test_tot - target_test$tot,
         test_m = test_m - target_test$m,
         test_f = test_f - target_test$f,
         chosen = case_when(set %in% calibrated$set ~ T, T ~ F)) %>%
  select(chosen, starts_with('prev'), starts_with('test'), -ends_with('ss')) %>%
  gather(case, val, -chosen) %>%
  mutate(cases = as.factor(case),
         case = ordered(case,
                        levels=c('prev_tot', 'prev_m', 'prev_f', 'prev_m0', 'prev_m1', 'prev_m2', 'prev_m3', 'prev_f0', 'prev_f1', 'prev_f2', 'prev_f3', 'test_tot', 'test_m', 'test_f'),
                        labels=c('Overall', 'Males', 'Females', 'Male 16-19', 'Male 20-24', 'Male 25-29', 'Male 30-35', 'Female 16-19', 'Female 20-24', 'Female 25-29', 'Female 30-35', 'Test Overall', 'Test Male', 'Test Female'))) %>%
  ggplot(aes(fill = chosen, x = val)) +
  geom_vline(aes(xintercept = 0), lty = 'dashed') +
  geom_density(alpha = 0.6) +
  facet_wrap(~case) +
  labs(x = 'Distribution of Difference Between Simulated and STRIVE Prevalence',
       y = 'Density',
       title = 'Analysis of Systematic Bias Within Prevalence Categories',
       fill = 'Is in top 50: ') +
  theme(legend.position = 'bottom') +
  scale_x_continuous(limits = c(-10, 10))


# Change target for plotting
target_long =
  target %>%
  gather(case, prev) %>%
  mutate(case = ordered(case,
                        levels=c('tot', 'm', 'f', 'm0', 'm1', 'm2', 'm3', 'f0', 'f1', 'f2', 'f3'),
                        labels=c('Overall', 'Males', 'Females', 'Male 16-19', 'Male 20-24', 'Male 25-29', 'Male 30-35', 'Female 16-19', 'Female 20-24', 'Female 25-29', 'Female 30-35')))


# Compute how good each parameter set is
data_set2 %>%
  filter(!is.na(prev_tot)) %>%
  mutate(chosen = case_when(set %in% calibrated$set ~ T, T ~ F)) %>%
  select(chosen, starts_with('prev'), starts_with('test'), -ends_with('ss')) %>%
  gather(case, val, -chosen) %>%
  mutate(cases = as.factor(case),
         case = ordered(case,
                        levels=c('prev_tot', 'prev_m', 'prev_f', 'prev_m0', 'prev_m1', 'prev_m2', 'prev_m3', 'prev_f0', 'prev_f1', 'prev_f2', 'prev_f3', 'test_tot', 'test_m', 'test_f'),
                        labels=c('Overall', 'Males', 'Females', 'Male 16-19', 'Male 20-24', 'Male 25-29', 'Male 30-35', 'Female 16-19', 'Female 20-24', 'Female 25-29', 'Female 30-35', 'Test Overall', 'test Male', 'Test Female'))) %>%
  ggplot(aes(fill = chosen, x = val)) +
  geom_vline(aes(xintercept = prev), target_long, lty = 'dashed') +
  geom_density(alpha = 0.6) +
  facet_wrap(~case) +
  labs(x = 'Distribution of Difference Between Simulated and STRIVE Prevalence',
       y = 'Density',
       title = 'Analysis of Systematic Bias Within Prevalence Categories',
       fill = 'Is in top 50: ') +
  theme(legend.position = 'bottom')


#######################################
##  EXAMINE PARAMETER DISTRIBUTIONS  ##
#######################################


data_set2 %>%
  mutate(chosen = case_when(set %in% calibrated$set ~ T, T ~ F)) %>%
  select(-starts_with('prev'), -starts_with('test'), -set) %>%
  gather(param, val, -chosen) %>%
  mutate(param = case_when(param == 'pup' ~ 'Urethra => Pharynx',
                           param == 'ppp' ~ 'Pharynx => Pharynx',
                           param == 'pru' ~ 'Rectum => Urethra',
                           param == 'ppu' ~ 'Pharynx => Urethra',
                           param == 'pur' ~ 'Urethra => Rectum',
                           param == 'prp' ~ 'Rectum => Pharynx',
                           param == 'ppr' ~ 'Pharynx => Rectum',
                           param == 'puu' ~ 'Urethra => Urethra',
                           param == 'symptoms_rectal' ~ 'Testing rate rectal',
                           param == 'symptoms_ural_male' ~ 'Testing rate urethral male',
                           param == 'symptoms_ural_female' ~ 'Testing rate urethral female',
                           T ~ param)) %>%
  ggplot(aes(x = val, fill = chosen)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~param) +
  theme(legend.position = 'bottom') +
  labs(x = 'Paramter Value',
       y = 'Density',
       fill = 'Is in top 50: ',
       title = 'Distribution of Sample Parameters Compared to Selected Parameters')


###########################################################################
##  MAKE A COPY OF THE CALIBRATION SIMS TO HELP WITH STORAGE MANAGEMENT  ##
###########################################################################


# Make a list of all the output file names
fnames = list.files(str_c('simulations/calibration/scenario_', scenario, '/'))


# Iterate over the chosen calibration sims
pb = progress_bar$new(total = nrow(calibrated))
for ( i in 1:nrow(calibrated) ){
  pb$tick()


  # Find all files to copy
  from = str_c('simulations/calibration/scenario_', scenario, '/')
  to = str_c('simulations/calibration_start_files/scenario_', scenario, '/')
  to_copy = fnames[str_detect(fnames, str_c('simulation_', calibrated$set[i], '_'))]


  # Copy over
  for ( j in 1:length(fnames) ){
    file.copy(str_c(from, to_copy[j]),
              str_c(to, to_copy[j]),
              overwrite = T, copy.mode = T, copy.date = T)
  }
}










