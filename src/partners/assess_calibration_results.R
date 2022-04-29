source('D:/Dropbox/Projects/setup.R')


# Read in target partner rates
partners_target = 
  read.csv('data/calibration_partnership_rates.csv',
         col.names = c('0', '1', '2', '3', '4')) %>%
  as_tibble() %>%
  mutate(num_partners = ordered(c('0', '1', '2-4', '5+'))) %>%
  gather('age_group', 'rate_target', -num_partners) %>%
  mutate(age_group = ordered(substr(age_group, 2, 2), labels = c('16-19', '20-24', '25-29', '30-34', '35'))) %>%
  select(age_group, num_partners, rate_target)


# Read in parameters
param =
  read.csv('simulations/parameters_partnerships.csv') %>%
  as_tibble()


# Find all output files
files =
  list.files("D:\\OneDrive - UNSW\\ng_remote\\simulations\\partnerships_calibration\\scenario_1") %>%
  enframe(name = NULL, value = 'file_name') %>%
  mutate(par_no = as.numeric(str_extract(file_name, '[0-9]+')),
         file_name = str_c('simulations/partnerships_calibration/scenario_1/', file_name)) %>%
  arrange(par_no)


# Read in files and assess how well they did
var_names = c('sim', 
              'A0P0', 'A0P1', 'A0P2', 'A0P3', 
              'A1P0', 'A1P1', 'A1P2', 'A1P3',
              'A2P0', 'A2P1', 'A2P2', 'A2P3',
              'A3P0', 'A3P1', 'A3P2', 'A3P3',
              'A4P0', 'A4P1', 'A4P2', 'A4P3',
              'runtime')
for ( i in 1:nrow(files) ){
  cat(red(i, 'of', nrow(files), '\n'))
  
  
  # Read in result
  temp = 
    read.csv(files$file_name[i]) %>%
    as_tibble()
  
  
  # Clean it up
  names(temp) = var_names
  temp = 
    temp %>%
    gather('group', 'rate', -sim, -runtime) %>%
    mutate(age_group = ordered(substr(group, 2, 2), levels = c('0', '1', '2', '3', '4'), labels = c('16-19', '20-24', '25-29', '30-34', '35')),
           num_partners = ordered(substr(group, 4, 4), levels = c('0', '1', '2', '3'), labels = c('0', '1', '2-4', '5+')),
           param = i) %>%
    select(param, sim, runtime, age_group, num_partners, rate) %>%
    arrange(sim, age_group, num_partners) %>%
    group_by(param, age_group, num_partners) %>%
    summarise(runtime = mean(runtime), rate = mean(rate)) %>%
    ungroup()

  
  # Store for later
  if ( i == 1 ){
    out = temp
  }else{
    out = rbind(out, temp)
  }
}


# Bring in target rates
out = 
  out %>%
  left_join(partners_target, by = c('age_group', 'num_partners'))


# Some cuts of the data
out %>%
  ggplot(aes(x = num_partners, y = (rate - rate_target)^2)) +
  geom_boxplot() +
  facet_wrap(~age_group) +
  labs(x = 'Number of Partners',
       y = 'Relative Error',
       title = 'Comparison of Simulated Partnership rates to Target Rates')


################################
##  FIND THE BEST PARAMETERS  ##
################################


# Compute the overall rss
wt = c(10, 10, 10, 10)
wt = wt/sum(wt)
rss = 
  out %>%
  filter(age_group != '35') %>%
  mutate(dist = (rate - rate_target)/rate_target) %>%
  group_by(param, age_group) %>%
  summarise(dist = 
              (num_partners=='0') * dist * wt[1] +
              (num_partners=='1') * dist * wt[2] +
              (num_partners=='2-4') * dist * wt[3] +
              (num_partners=='5+') * dist * wt[4]) %>%
  ungroup(age_group) %>%
  summarise(dist = sum(dist)) %>%
  arrange(abs(dist))


# Show the distribution
rss %>%
  ggplot(aes(x = dist)) +
  geom_histogram(bins = 50) +
  labs(x = 'Residual Sum of Squares',
       y = 'Count',
       title = 'Distribution of the RSS')


# Take a list of 'best' ones
best =
  rss %>%
  slice(1:10)


# Plot
out %>%
  filter(param %in% best$param) %>%
  ggplot(aes(x = num_partners)) +
  geom_col(aes(y = rate_target/nrow(best))) +
  geom_jitter(aes(y = rate, colour = as.factor(param)), width = 0.1) +
  facet_wrap(~age_group) +
  labs(x = 'Number of Partners',
       y = 'Partnership Rates',
       title = 'Comparison of Target Rates to Calibration Rates') +
  theme(legend.position = 'none')


##############################################
##  LOOK AT THE SELECTED PARAMETER CHOICES  ##
##############################################


# Which were the best
param = 
  read.csv('simulations/parameters_partnerships.csv') %>%
  mutate(param = 1:20000) %>%
  filter(param %in% out$param) %>%
  mutate(best = case_when(param %in% best$param ~ T, T ~ F)) %>%
  select(-l4, -h4) %>%
  gather(var, val, -param, -best) %>%
  mutate(risk = case_when(str_detect(var, 'h') ~ 'High-risk', T ~ 'Low-risk'),
         age = case_when(substr(var, 2, 2) == 0 ~ '16-19',
                         substr(var, 2, 2) == 1 ~ '20-24',
                         substr(var, 2, 2) == 2 ~ '25-29',
                         substr(var, 2, 2) == 3 ~ '30-34',
                         T ~ '35'))

# Make graphs
param %>%
  ggplot(aes(x = val)) +
  geom_density() +
  geom_point(aes(y = 0), data = filter(param, best)) +
  facet_grid(rows = vars(risk), cols = vars(age), scales = 'free') +
  labs(x = 'Parameter',
       y = 'Calibration Set',
       title = 'Comparison of Optimal Parameters to Full Calibration Set')


# Save a list of the optimal params
best %>%
  select(prt_param_no = param) %>%
  arrange(prt_param_no) %>%
  write.csv('simulations/parameters_partnerships_calibrated.csv', row.names = F)





