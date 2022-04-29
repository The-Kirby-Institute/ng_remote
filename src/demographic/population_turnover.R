source('D:/Dropbox/Projects/setup.R')


# Load in 2011 LGAs
lgas = 
  st_read('D:/Data/abs/spatial/LGA_2011_AUST.gpkg', quiet = T) %>%
  filter(str_detect(STE_NAME11, 'Northern'))


# Read in turnover rates
turnover =
  read_excel('data/turnover/pop_turnover.xls') %>%
  mutate(LGA_CODE11 = str_extract(lga, '[0-9]+'),
         turnover = 100 * turnover/(1000 * 5),
         arrivals = 100 * (arrivals/pop_2006_corrected)/5,
         departures = 100 * (departures/pop_2006_corrected)/5) %>%
  select(LGA_CODE11, turnover, arrivals, departures)


# Convert to spatial
turnover = 
  lgas %>%
  left_join(turnover, by = 'LGA_CODE11')


# Read in communites
communities = 
  st_read('data/2016_ILOC_shape/ILOC_2016_AUST.shp', quiet = T) %>%
  filter(str_detect(STATE_NAME, 'Northern')) %>%
  select(ILOC_CODE) %>%
  left_join(readRDS('data/communities_for_analysis.Rds'), by = 'ILOC_CODE') %>%
  filter(!is.na(group)) %>%
  st_centroid()


# Make map
turnover %>%
  tm_shape() +
  tm_polygons('turnover', 
              style = 'cont', 
              title = 'Mean Annual Turnover (%)') +
  tm_shape(communities) +
  tm_dots(title = 'Simulation Scenario',
          col = 'group',
          size = 1,
          palette = 'Dark2') +
  tm_layout(legend.outside = T)


# Compute mean turnover rates for each scenario
turnover = 
  communities %>%
  st_centroid() %>%
  st_intersection(turnover)
st_geometry(turnover) = NULL  
turnover %>%
  as_tibble() %>%
  group_by(group) %>%
  summarise(t = mean(turnover),
            t_ci = qnorm(0.975)*sqrt(var(turnover))/sqrt(n()),
            a = mean(arrivals),
            a_ci = qnorm(0.975)*sqrt(var(arrivals))/sqrt(n()),
            d = mean(departures),
            d_ci = qnorm(0.975)*sqrt(var(departures))/sqrt(n())) %>%
  gather(est, val, -group) %>%
  mutate(qty = case_when(str_detect(est, 'ci') ~ 'ci', T ~ 'm'),
         est = case_when(str_detect(est, 't') ~ 'Turnover',
                         str_detect(est, 'a') ~ 'Arrivals',
                         T ~ 'Departures')) %>%
  spread(qty, val) %>%
  ggplot(aes(x = est, y = m)) +
  geom_errorbar(aes(ymin = m-ci, ymax = m+ci), width = 0.2) +
  facet_wrap(~group) +
  geom_point(colour = 'red') +
  labs(x = 'Component of Population Turnover',
       y = 'Annual Turnover Rate (Estimate +/- 95% CI)',
       title = 'Estimated Population Turnover Rate (2006-11 Census Data)')


turnover %>%
  as_tibble() %>%
  group_by(group) %>%
  summarise(t = mean(turnover),
            t_ci = qnorm(0.975)*sqrt(var(turnover))/sqrt(n()),
            a = mean(arrivals),
            a_ci = qnorm(0.975)*sqrt(var(arrivals))/sqrt(n()),
            d = mean(departures),
            d_ci = qnorm(0.975)*sqrt(var(departures))/sqrt(n())) %>%
  ggplot(aes(x = group, y = t)) +
  geom_errorbar(aes(ymin = t-t_ci, ymax = t+t_ci), width = 0.2) +
  geom_point(colour = 'red') +
  labs(x = 'Population Scenario',
       y = 'Annual Turnover Rate (Estimate +/- 95% CI)',
       title = 'Estimated Population Turnover Rate (2006-11 Census Data)')


turnover %>%
  as_tibble() %>%
  group_by(group) %>%
  summarise(t = mean(turnover),
            t_ci = qnorm(0.975)*sqrt(var(turnover))/sqrt(n()),
            a = mean(arrivals),
            a_ci = qnorm(0.975)*sqrt(var(arrivals))/sqrt(n()),
            d = mean(departures),
            d_ci = qnorm(0.975)*sqrt(var(departures))/sqrt(n())) %>%
  ggplot(aes(x = group, y = t)) +
  geom_errorbar(aes(ymin = t-t_ci, ymax = t+t_ci), width = 0.2) +
  geom_point(colour = 'red') +
  labs(x = 'Population Scenario',
       y = 'Annual Turnover Rate (Estimate +/- 95% CI)',
       title = 'Estimated Population Turnover Rate (2006-11 Census Data)')

