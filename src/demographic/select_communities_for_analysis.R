library(tidyverse)
library(sf)



################################
##  PROCESS DEMOGRAPHIC DATA  ##
################################


# Read in stats by age and sex
data = 
  read.csv("data/2016_ATSIP_ILOC_for_AUS_short-header/2016 Census ATSIP Indigenous Locations for AUST/2016Census_I01A_AUS_ILOC.csv") %>%
  as_tibble() %>%
  select(ILOC_CODE_2016, Tot_p_TotP, Tot_p_Indig_P, Tot_p_NonInd_P, Tot_p_IndigStat_nsP)


# Read in spatial boundaries
iloc = 
  st_read("data/2016_ILOC_shape/ILOC_2016_AUST.shp", quiet = T) %>%
  filter(STATE_NAME == "Northern Territory") %>%
  filter(AREA_SQKM > 0)


# Filter data to just the NT
data =
  data %>%
  filter(ILOC_CODE_2016 %in% iloc$ILOC_CODE)


# Reshape so that Indig status is a factor
data = 
  data %>%
  rename(totPop = Tot_p_TotP,
         indig = Tot_p_Indig_P,
         nonInd = Tot_p_NonInd_P,
         noStat = Tot_p_IndigStat_nsP) %>%
  gather(response, count, -ILOC_CODE_2016, -totPop) %>%
  mutate(proportion = 100 * count/totPop)


# Look at the ratio of Indig to non-Indig
data %>%
  ggplot(aes(x = fct_reorder(as.factor(ILOC_CODE_2016), desc(totPop)))) +
  geom_col(aes(y = count, fill = response))


#############################################
##  INFORM DEMOGRAPHIC BY REMOTENESS AREA  ##
#############################################


# Looks good, now I'll quantify remoteness
ra = 
  st_read("data/2016_RA_shape/RA_2016_AUST.shp", quiet = T) %>%
  select(RA_NAME16)


# Tag ILOC by RA
temp =
  iloc %>%
  select(ILOC_CODE) %>%
  st_point_on_surface() %>%
  st_intersection(ra)
st_geometry(temp) = NULL
iloc = left_join(iloc, temp, by = "ILOC_CODE")
data = left_join(data, temp, by = c("ILOC_CODE_2016" = "ILOC_CODE"))


# Make a quick map to check the result
# tm_shape(ra) +
#   tm_borders() +
#   tm_shape(st_point_on_surface(iloc)) +
#   tm_dots("RA_NAME16")


# Plot population by remoteness area
data %>%
  ggplot(aes(x = fct_reorder(as.factor(ILOC_CODE_2016), desc(totPop)))) +
  geom_col(aes(y = count, fill = response)) +
  facet_wrap(~RA_NAME16, scales = "free")


#######################################################
##  INFORM DEMOGRAPHIC BY DISTANCE FROM URBAN AREAS  ##
#######################################################


# Read in Urban Centers and Localities
ucl = 
  st_read("data/2016_UCL_shape/UCL_2016_AUST.shp", quiet = T) %>%
  filter(!str_detect(UCL_NAME16, "Remainder"))


# Read in significant urban areas
sua = 
  st_read("data/2016_SUA_shape/SUA_2016_AUST.shp", quiet = T) %>%
  filter(!str_detect(SUA_NAME16, "Not in any"))


# Setup storage of all the UCL and SUA distances
temp = 
  iloc %>%
  select(ILOC_CODE) %>%
  mutate(nearest_ucl = NA, nearest_sua = NA) %>%
  st_centroid()


# Now compute the distance of all the ILOCs to their nearest one
pb = txtProgressBar(1, nrow(temp), style = 3)
for( i in 1:nrow(temp) ){
  setTxtProgressBar(pb, i)
  
  suppressMessages({
    # Distance to nearest UCL
    j = st_nearest_feature(temp[i,], ucl)
    temp$nearest_ucl[i] = st_distance(temp[i,], ucl[j,])
    
    # Distance to nearest SUA
    j = st_nearest_feature(temp[i,], sua)
    temp$nearest_sua[i] = st_distance(temp[i,], sua[j,])
  })
}


# Join data back up with spatial data
st_geometry(temp) = NULL
iloc = left_join(iloc, temp, by = "ILOC_CODE")


# Join remoteness data up with demographic data
temp = 
  iloc %>%
  select(ILOC_CODE, AREA_SQKM, nearest_ucl, nearest_sua)
st_geometry(temp) = NULL
data = left_join(data, temp, by = c("ILOC_CODE_2016" = "ILOC_CODE"))


# Plot results
data %>%
  ggplot(aes(x = fct_reorder(as.factor(ILOC_CODE_2016), desc(totPop)))) +
  geom_col(aes(y = count, fill = response)) +
  geom_point(aes(y = (15000/log(1 + max(nearest_sua, na.rm = T))) * log(1 + nearest_sua), colour = RA_NAME16))


# Another plot
data %>%
  filter(response == "indig") %>%
  ggplot(aes(x = proportion, 
             y = log(AREA_SQKM),
             size = nearest_sua,
             colour = RA_NAME16)) +
  geom_abline(aes(intercept = log(100), slope = 0), lty = "dashed") +
  geom_vline(aes(xintercept = 50), lty = "dashed") +
  geom_point(alpha = 0.6) +
  labs(x = "Proportion of Population Reported as Aboriginal",
       y = expression(paste("Log of ILOC Area (k", m^{2}, ")")),
       colour = "Remoteness") +
  theme(legend.position = "bottom") +
  scale_size_continuous(guide = F)
ggsave("graphs/ilocs_by_proportion_indegenous_and_remoteness.jpeg", height = 6, width = 6)


# Another plot
data.st = 
  iloc %>%
  select(-RA_NAME16, -AREA_SQKM, -starts_with("nearest")) %>%
  left_join(data, by = c("ILOC_CODE" = "ILOC_CODE_2016")) %>%
  st_as_sf()


# Data
data = data.st
st_geometry(data) = NULL  
data %>% 
  filter(response == "indig",
         proportion > 50,
         RA_NAME16 == "Very Remote Australia") %>%
  ggplot(aes(x = proportion, y = AREA_SQKM)) + geom_point()


################################
##  DECIDE WHICH ONES TO USE  ##
################################


# Use this to make sure youre happy with the map
data.st %>%
  filter(response == "indig",
         proportion > 50,
         AREA_SQKM < 100,
         RA_NAME16 == "Very Remote Australia") %>%
  tm_shape() +
  tm_borders() +
  tm_fill("totPop",
          alpha = 0.6)


# Extract the corresponding IDs
temp = 
  data.st %>%
  filter(response == "indig",
         proportion > 50,
         AREA_SQKM < 100,
         RA_NAME16 == "Very Remote Australia") %>%
  select(ILOC_CODE)
st_geometry(temp) = NULL


# Extract the data
data = 
  data %>%
  as_tibble() %>%
  filter(ILOC_CODE %in% temp$ILOC_CODE) %>%
  select(-proportion) %>%
  spread(response, count)


#################################
##  DECIDE ON COMMUNITY SIZES  ##
#################################


# Aggregate up to find town size
towns =
  data %>%
  group_by(ILOC_CODE) %>%
  summarise(pop = sum(indig + nonInd + noStat))


# Work out how to group the towns
clusters = kmeans(towns$pop, centers = 3)
clusters$centers


# Save cluster info for later
scenario = 
  clusters$centers %>%
  as_tibble() %>%
  mutate(cluster = 1:3) %>%
  rename(mean_size = V1) %>%
  arrange((mean_size)) %>%
  mutate(scenario_num = 1:3,
         scenario_num = ordered(scenario_num, 1:4),
         scenario_use = c(200, 750, 2000))


# Save for later
write.csv(scenario, "data/scenarios.csv", row.names = F)


# Put groupings in with the data
towns = 
  towns %>%
  mutate(cluster = clusters$cluster) %>%
  left_join(scenario, by = "cluster") %>%
  select(ILOC_CODE, pop, group = scenario_num)


# Make a plot of the groupings
towns %>%
  ggplot(aes(x=pop, fill = as.factor(group))) + 
  geom_histogram(binwidth = 50, palette = 1) +
  labs(x = "Population Size",
       y = "Number of ILOCs With This Population Size",
       fill = "Range") +
  scale_fill_discrete(labels = c("0-449", "450-1,999", "2,000 or more")) +
  theme(legend.position = "bottom")
ggsave("graphs/iloc_population_distrbution.jpeg", width = 6, height = 6)


# Put in with the full dataset
towns = select(towns, -pop)
data =
  data %>%
  left_join(towns, by = "ILOC_CODE")


# Save for later
saveRDS(data, "data/communities_for_analysis.Rds")


##########################################
##  JOIN UP SCENARIOS WITH CENSUS DATA  ##
##########################################


# Read in data on which towns are being analysed
data = readRDS("data/communities_for_analysis.Rds")


# Extract corresponding demographic data
census = 
  read.csv("data/2016_ATSIP_ILOC_for_AUS_short-header/2016 Census ATSIP Indigenous Locations for AUST/2016Census_I03A_AUS_ILOC.csv") %>%
  left_join(read.csv("data/2016_ATSIP_ILOC_for_AUS_short-header/2016 Census ATSIP Indigenous Locations for AUST/2016Census_I03B_AUS_ILOC.csv"), by = "ILOC_CODE_2016") %>%
  left_join(read.csv("data/2016_ATSIP_ILOC_for_AUS_short-header/2016 Census ATSIP Indigenous Locations for AUST/2016Census_I03C_AUS_ILOC.csv"), by = "ILOC_CODE_2016") %>%
  as_tibble() %>%
  filter(ILOC_CODE_2016 %in% data$ILOC_CODE) %>%
  select(ILOC_CODE_2016, matches("Age_[0-9]+_[0-9]+_TotM") | matches("Age_[0-9]+_[0-9]+_TotF"))
# select(ILOC_CODE_2016, starts_with("Age") & (ends_with("TotM") | ends_with("TotF")))


# Work out the distributions by age group
census = 
  census %>%
  gather(group, count, -ILOC_CODE_2016) %>%
  mutate(sex = str_extract(group, "[MF]"),
         age_lower = str_extract(group, "Age_[0-9]+"),
         age_lower = as.numeric(str_extract(age_lower, "[0-9]+")),
         age_upper = str_extract(group, "[0-9]+_Tot"),
         age_upper = as.numeric(str_extract(age_upper, "[0-9]+")),
         age_upper = case_when(age_upper == age_lower ~ age_upper + 1, T ~ age_upper)) %>%
  select(-group)


# Join all the data up
data = 
  data %>%
  select(ILOC_CODE, totPop, scenario = group) %>%
  left_join(census, by = c("ILOC_CODE" = "ILOC_CODE_2016")) 


######################################
##  IDENTIFY THE AGE DISTRIBUTIONS  ##
######################################


# Work out the totals by town
data = 
  data %>%
  spread(sex, count) %>%
  group_by(ILOC_CODE) %>%
  summarise(totM = sum(M),
            totF = sum(F),
            tot = sum(M + F)) %>%
  left_join(data, by = "ILOC_CODE") %>%
  select(-totPop) %>%
  mutate(weight = case_when(sex == "M" ~ count/totM,
                            sex == "F" ~ count/totF))


# Compute the distributions
means = 
  data %>%
  group_by(scenario, sex, age_lower, age_upper) %>%
  summarise(count = median(count)) %>%
  ungroup(age_lower) %>%
  mutate(ave_smooth = rollapply(count, 4, mean, fill=NA, align='center'),
         ave = count / sum(count)) %>%
  ungroup()


# Check these distributions
data %>%
  ggplot(aes(x=age_lower, y=count)) +
  geom_line(aes(group=ILOC_CODE)) +
  stat_summary(aes(y = count), fun=mean, colour="red", geom="line", size=2) +
  facet_grid(rows = vars(sex),
             cols = vars(scenario),
             scales = 'free')


# Make a different plot
means %>%
  ggplot(aes(x=age_lower, colour=sex)) +
  geom_line(aes(y=ave)) +
  facet_wrap(~scenario)


# Save the smoothed version
means = 
  means %>%
  select(-ave, -count) %>%
  rename(ave = ave_smooth) %>%
  filter(!is.na(ave)) %>%
  group_by(scenario, sex) %>%
  mutate(ave = ave / sum(ave)) %>%
  ungroup()


# Save for later
saveRDS(means, "data/age_distributions.Rds")
write.csv(means, "data/age_distributions.csv", row.names = F)


############################################
##  WORK OUT THE DISTRIBUTIONS BY GENDER  ##
############################################


# Compute the sex distribution for each scenario
sex = 
  readRDS("data/communities_for_analysis.Rds") %>%
  select(ILOC_CODE_2016 = ILOC_CODE,
         scenario = group) %>%
  left_join(census, by = "ILOC_CODE_2016") %>%
  spread(sex, count) %>%
  group_by(ILOC_CODE_2016, scenario) %>%
  summarise(M = sum(M), F = sum(F)) %>%
  ungroup() %>%
  mutate(pMale = M / (M + F)) %>%
  group_by(scenario) %>%
  summarise(M = sum(M),
            F = sum(F),
            pMale = mean(pMale))


# Save for later
write.csv(sex, "data/sex_distributions.csv", row.names = F)
saveRDS(sex, "data/sex_distributions.Rds")


########################################
##  PLOTS FOR THE TECHNICAL ANALYSIS  ##
########################################


## PLOT OF THE INCLUDED ILOCS
library(tmaptools)
tmap_mode("plot")


# Form a boundary for the remoteness areas
nt = 
  ra[41:43,] %>%
  mutate(Remoteness = ordered(RA_NAME16, 
                              c("Outer Regional Australia", "Remote Australia", "Very Remote Australia"),
                              c("Outer Regional", "Remote", "Very Remote")))


# Convert some regions to points
plot_dat =
  iloc %>%
  left_join(towns, by = "ILOC_CODE") %>%
  filter(!is.na(group)) %>%
  mutate(group = ordered(group, 1:3, c("0-449", "450-1,999", "Over 2000")))
plot_dat = st_centroid(plot_dat)


# Make the plot
p = 
  tm_shape(nt) +
  tm_borders() +
  tm_fill("Remoteness",
          palette = "Blues") +
  tm_shape(plot_dat) +
  tm_dots("group",
          title = "Population",
          palette = "YlOrRd",
          size = 0.5) +
  tm_legend("",
            scale = 1,
            legend.outside = T,
            legend.outside.position = c("bottom"),
            legend.stack = "horizontal",
            legend.bg.alpha = 1,
            legend.bg.color = "white") +
  tm_layout(frame = F)
p
tmap_save(p, "graphs/ilocs.jpeg", 6, 6)


# Excluded ILOCS
tmap_mode("view")
iloc %>%
  # filter(!(ILOC_CODE %in% towns$ILOC_CODE)) %>%
  left_join(read.csv("data/2016_ATSIP_ILOC_for_AUS_short-header/2016 Census ATSIP Indigenous Locations for AUST/2016Census_I01A_AUS_ILOC.csv"), 
            by = c("ILOC_CODE" = "ILOC_CODE_2016")) %>%
  mutate(decision = case_when(RA_NAME16 != "Very Remote Australia" ~ "Not 'Very Remote'",
                              AREA_SQKM > 100 ~ "Larger than 15,000km2 (multiple communities)",
                              Tot_p_Indig_P/Tot_p_TotP < 0.5 ~ "less than 50% Aboriginal",
                              T ~ "Used for analysis")) %>%
  tm_shape() +
  tm_fill("decision", alpha = 0.6) +
  tm_borders() +
  tm_basemap("OpenStreetMap")


##  PLOT OF DEMOGRAPHIC DISTRIBUTIONS
p = 
  data %>%
  ungroup() %>%
  mutate(scenario = ordered(scenario, 1:3, str_c("Scenario ", c("1: Population size 0 to 449", "2: Population size 450-1,999", "3: Population size over 2000")))) %>%
  ggplot(aes(x = age_lower + 2, y = ave, fill = sex)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~scenario, nrow = 3, ncol = 1) +
  labs(x = "Age Group",
       y = "Proportion",
       fill = "Sex") +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0, 65, by = 5))
p
ggsave("graphs/demographic_distribution.jpeg", p,  width = 6, height = 7.5)  








