# -*- coding: utf-8 -*-
"""
CODE FOR DOING ALL THE DEMOGRAPHIC STUFF

Created on Mon Mar 15 11:04:41 2021

@author: nicol
"""


#%%  Setup Libraries


# Load libraries
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Read in parameters
sim_parameters = pd.read_csv("data/param.csv")
scenario_global = sim_parameters.scenario[0]


#%% FUN setup_data()
def setup_data(scenario = scenario_global, run_mode = 'serial'):


    # What doing
    ( print('Demographic attributes: scenario ' + str(scenario)) if run_mode == 'serial' else [] )


    # Load distribution of community size
    size = pd.read_csv("data/scenarios.csv")
    size = size.loc[size["scenario_num"] == scenario,]
    n = size.scenario_use.iloc[0]


    # Load distribution of sex
    sex_dist = pd.read_csv("data/sex_distributions.csv")
    sex_dist = sex_dist.loc[sex_dist["scenario"] == scenario,]


    # Load distribution of age
    age_dist = pd.read_csv("data/age_distributions.csv")
    age_dist = age_dist.loc[age_dist["scenario"] == scenario,]


    # Convert sex from M/F to 1/0
    age_dist.loc[age_dist.sex == 'M', 'sex'] = 1
    age_dist.loc[age_dist.sex == 'F', 'sex'] = 0


    # Modify age distribution to agree with the age bounds 16 and 35
    age_dist = age_dist.loc[age_dist.age_lower >= 15,:]
    age_dist = age_dist.loc[age_dist.age_lower <= 35,:]
    age_dist.loc[age_dist.age_lower == 15, 'ave'] = (4/5) * age_dist.loc[age_dist.age_lower == 15, 'ave']
    age_dist.loc[age_dist.age_lower == 35, 'ave'] = (1/5) * age_dist.loc[age_dist.age_lower == 35, 'ave']


    # Renormalise age distribution
    age_dist.loc[age_dist.sex == 1, 'ave'] = (1/np.sum(age_dist.loc[age_dist.sex == 1, 'ave'])) * age_dist.loc[age_dist.sex == 1, 'ave']
    age_dist.loc[age_dist.sex == 0, 'ave'] = (1/np.sum(age_dist.loc[age_dist.sex == 0, 'ave'])) * age_dist.loc[age_dist.sex == 0, 'ave']


    # Modify the quoted upper and lower bounds
    age_dist.loc[age_dist.age_lower == 15, 'age_lower'] = 16
    age_dist.loc[:, 'age_upper'] = age_dist.loc[:, 'age_upper'] + 1
    age_dist.loc[age_dist.age_upper == 40, 'age_upper'] = 36


    # Compute the CDF for the age distribution
    age_dist.loc[:, 'cdf'] = np.cumsum(age_dist.ave)
    age_dist.loc[age_dist.sex == 1, 'cdf'] = age_dist.loc[age_dist.sex == 1, 'cdf'] - 1
    
    
    # Load population turnover rate
    pop_turnover = pd.read_csv('data/population_dynamics.csv')
    

    # Load distribution of sexual orientation
    orientation_dist = pd.read_csv("data/orientation_distribution.csv")


    # Load distribution of number of sexual partners
    # partners_dist = pd.read_csv("data/partners_distribution_DELETE_THIS.csv")
    partners_dist = pd.read_csv("data/calibration_partnership_rates.csv")
    
    
    # Load distribution of number of sexual partners
    # partners_dist = pd.read_csv("data/partners_distribution_DELETE_THIS.csv")
    sexual_acts_dist = pd.read_csv("data/calibration_sexual_acts.csv")


    # Compute target distributions
    targetF = age_dist.ave[age_dist.sex == 0].reset_index(drop=True) * (1-sex_dist['pMale'].iloc[0]) * n
    targetM = age_dist.ave[age_dist.sex == 1].reset_index(drop=True) * sex_dist['pMale'].iloc[0] * n
    target = {'0': targetF, '1': targetM}
    
    
    # Load test and treat dynamics
    testing_rates = pd.read_csv("data/testing_rates.csv")
    testing_param = pd.read_csv("data/testing_dynamics.csv")


    # Store all population data in a big dictionary
    pop_parameters = {'scenario': scenario,
                      'size': size,
                      'n': n,
                      'sex_dist': sex_dist,
                      'age_dist': age_dist,
                      'orientation_dist': orientation_dist,
                      'partners_dist': partners_dist,
                      'sexual_acts_dist': sexual_acts_dist,
                      'target': target,
                      'pop_turnover': pop_turnover,
                      'testing_param': testing_param,
                      'testing_rates': testing_rates}


    # Return the data
    return pop_parameters


#%% FUN generate_population()
def generate_population(pop_parameters,
                        n_generate = 'initilise',
                        prop_infected = sim_parameters.init_prob_exposed[0],
                        t = sim_parameters.partner_burn_in[0]):


    # Set the number of individuals in the population
    n = pop_parameters['n'] if n_generate == 'initilise' else int(n_generate)


    ## INITILISE DATA.FRAME
    meta_init = initilise_meta(n)


    ## SET GENDER
    meta_init.loc[meta_init.gender < pop_parameters['sex_dist'].pMale.iloc[0], 'gender'] = 1
    meta_init.loc[meta_init.gender < 1, 'gender'] = 0


    ## SET AGE
    meta_init = initilise_age(meta_init, pop_parameters['age_dist'], n_generate)


    ## SET SEXUAL ORIENTATION
    meta_init = initilise_orientation(meta_init, pop_parameters['orientation_dist'])


    ## SET INFECTION RISK
    meta_init = initilise_risk(meta_init, pop_parameters['partners_dist'])
    
    
    ## SET SEXUAL ACTS
    meta_init = initilise_sexual_acts(meta_init, pop_parameters['sexual_acts_dist'])


    ## SET INFECTION STATUS
    meta_init = initilise_infections(meta_init, prop_infected, t)


    ## CLEANUP VARIABLE TYPES
    meta_init = cleanup_dtypes(meta_init)


    # End it
    return meta_init



#%% FUN initilise_partner_matrix()
def initilise_partner_matrix(pop_parameters):
    n = int(1.5 * pop_parameters['n'])
    out = np.zeros((n, n))
    return out


#%% FUN initilise_partner_durations()
def initilise_partner_duration(pop_parameters):
    n = int(1.5 * pop_parameters['n'])
    out =  float("inf") * np.ones((n, n))
    return out


#%%  HELPER initilise_meta()
#
#
# Function which makes a new version of meta_init
#
#
def initilise_meta(n):


    # Initilise data.frame containing attributes of the population
    meta_init = pd.DataFrame(columns = ["gender",               # Binary 0/1 are they male
                                   "age",                  # Age of the individual
                                   "age_group",            # Age group (16-19, 20-24, 25-29, 30-36)
                                   "orientation",          # Sexual orientation
                                   "risk",                 # Infection risk level
                                   "partner",              # The indivdual's long-term partner
                                   "counter",              # Counter of how many partners they've had
                                   "has_oral",             # Indicator variable of whether they have oral
                                   "has_anal",             # Indicator variable of whether they have anal
                                   "has_sex",              # Indicator variable of whether they have sex
                                   "state",                # Their current infection state (S,I,R)
                                   "site0",                # Binary 0/1 are they infected at anatomical site 0
                                   "site1",                # Binary 0/1 are they infected at anatomical site 1
                                   "site2",                # Binary 0/1 are they infected at anatomical site 2
                                   "site0_t0",             # The simulation time that they became infected at site 0
                                   "site1_t0",             # The simulation time that they became infected at site 1
                                   "site2_t0",             # The simulation time that they became infected at site 0
                                   "site0_t1",             # The simulation time that they recover at site 0
                                   "site1_t1",             # The simulation time that they recover at site 1
                                   "site2_t1",             # The simulation time that they recover at site 2
                                   'site0_symptoms',
                                   'site1_symptoms',
                                   'site2_symptoms',
                                   "treatment_threshold",  # Threshold in [0,1] indicating when they'll get treatment
                                   "recovery_time",        # Simulation time that they get treatment
                                   'import_time',           # Simulation time that they were imported
                                   'test_time',             # Simulation time of next test
                                   'test_time_last',        # Time of last test
                                   'test_reason_last',        # Reason for last test
                                   'treatment_time',        # Simulation time of next treatment
                                   'vaccinated',             # Booian are they currently considered vaccinated
                                   'vaccination_t0',        # Time they were vaccinated
                                   'vaccination_t1',        # Time their vaccine wears off
                                   'booster_t0',            # Time they get a booster
                                   ])         
    
    
    # Set variable types
    meta_init.gender.astype("int64")
    meta_init.age.astype("float64")
    meta_init.age_group.astype("int64")
    meta_init.orientation.astype("int64")
    meta_init.risk.astype("int64")
    meta_init.partner.astype("int64")
    meta_init.counter.astype("int64")
    meta_init.has_oral.astype('bool')
    meta_init.has_anal.astype('bool')
    meta_init.has_sex.astype('bool')
    meta_init.state.astype("category")
    meta_init.site0.astype("int64")
    meta_init.site1.astype("int64")
    meta_init.site2.astype("int64")
    meta_init.site0_t0.astype("float64")
    meta_init.site1_t0.astype("float64")
    meta_init.site2_t0.astype("float64")
    meta_init.site0_t1.astype("float64")
    meta_init.site1_t1.astype("float64")
    meta_init.site2_t1.astype("float64")
    meta_init.site0_symptoms.astype('bool')
    meta_init.site1_symptoms.astype('bool')
    meta_init.site2_symptoms.astype('bool')
    meta_init.treatment_threshold.astype("float64")
    meta_init.recovery_time.astype("float64")
    meta_init.import_time.astype("float64")
    meta_init.test_time.astype("float64")
    meta_init.test_time_last.astype("float64")
    meta_init.test_reason_last.astype("int64")
    meta_init.treatment_time.astype("float64")
    meta_init.vaccinated.astype("bool")
    meta_init.vaccination_t0.astype("float64")
    meta_init.vaccination_t1.astype("float64")
    meta_init.booster_t0.astype("float64")


    # Set default values
    meta_init.loc[:, 'gender'] = np.random.random(n)
    meta_init.loc[:, 'age'] = np.random.random(n)
    meta_init.loc[:, 'age_group'] = -1
    meta_init.loc[:, 'orientation'] = np.random.random(n)
    meta_init.loc[:, 'risk'] = np.random.random(n)
    meta_init.loc[:, 'partner'] = -1
    meta_init.loc[:, 'counter'] = 0
    meta_init.loc[:, 'has_oral'] = False
    meta_init.loc[:, 'has_anal'] = False
    meta_init.loc[:, 'has_sex'] = False
    meta_init.loc[:, 'state'] = 'S'
    meta_init.loc[:, 'site0'] = 0
    meta_init.loc[:, 'site1'] = 0
    meta_init.loc[:, 'site2'] = 0
    meta_init.loc[:, 'site0_t0'] = float("inf")
    meta_init.loc[:, 'site1_t0'] = float("inf")
    meta_init.loc[:, 'site2_t0'] = float("inf")
    meta_init.loc[:, 'site0_t1'] = float("inf")
    meta_init.loc[:, 'site1_t1'] = float("inf")
    meta_init.loc[:, 'site2_t1'] = float("inf")
    meta_init.loc[:, 'site0_symptoms'] = False
    meta_init.loc[:, 'site1_symptoms'] = False
    meta_init.loc[:, 'site2_symptoms'] = False
    meta_init.loc[:, 'treatment_threshold'] = np.random.random(n)
    meta_init.loc[:, 'recovery_time'] = float("inf")
    meta_init.loc[:, 'import_time'] = 0
    meta_init.loc[:, 'test_time'] = float('Inf')
    meta_init.loc[:, 'test_time_last'] = -float('Inf')
    meta_init.loc[:, 'test_reason_last'] = 0
    meta_init.loc[:, 'treatment_time'] = float('Inf')
    meta_init.loc[:, 'vaccinated'] = False
    meta_init.loc[:, 'vaccination_t0'] = float("inf")
    meta_init.loc[:, 'vaccination_t1'] = float("inf")
    meta_init.loc[:, 'booster_t0'] = float("inf")


    return meta_init


#%% HELPER initilise_age()
#
#
# Sets the ages
#
#
def initilise_age(meta_init, age_dist, n_generate):


    # Female 16-19
    who = ((meta_init.gender == 0) & (meta_init.age < age_dist.cdf.iloc[0]))
    meta_init.loc[who, 'age'] = 16 + 4 * np.random.random(sum(who))


    # Female 20-24
    who = ((meta_init.gender == 0) & (meta_init.age < age_dist.cdf.iloc[1]))
    meta_init.loc[who, 'age'] = 20 + 5 * np.random.random(sum(who))


    # Female 25-29
    who = ((meta_init.gender == 0) & (meta_init.age < age_dist.cdf.iloc[2]))
    meta_init.loc[who, 'age'] = 25 + 5 * np.random.random(sum(who))


    # Female 30-34
    who = ((meta_init.gender == 0) & (meta_init.age < age_dist.cdf.iloc[3]))
    meta_init.loc[who, 'age'] = 30 + 5 * np.random.random(sum(who))


    # Female 35
    who = ((meta_init.gender == 0) & (meta_init.age < age_dist.cdf.iloc[4]))
    meta_init.loc[who, 'age'] = 35 + np.random.random(sum(who))


    # Male 16-19
    who = ((meta_init.gender == 1) & (meta_init.age < age_dist.cdf.iloc[5]))
    meta_init.loc[who, 'age'] = 16 + 4 * np.random.random(sum(who))


    # Male 20-24
    who = ((meta_init.gender == 1) & (meta_init.age < age_dist.cdf.iloc[6]))
    meta_init.loc[who, 'age'] = 20 + 5 * np.random.random(sum(who))


    # Male 25-29
    who = ((meta_init.gender == 1) & (meta_init.age < age_dist.cdf.iloc[7]))
    meta_init.loc[who, 'age'] = 25 + 5 * np.random.random(sum(who))


    # Male 30-34
    who = ((meta_init.gender == 1) & (meta_init.age < age_dist.cdf.iloc[8]))
    meta_init.loc[who, 'age'] = 30 + 5 * np.random.random(sum(who))


    # Male 35
    who = ((meta_init.gender == 1) & (meta_init.age < age_dist.cdf.iloc[9]))
    meta_init.loc[who, 'age'] = 35 + np.random.random(sum(who))


    # Age group
    meta_init.age_group = np.floor((meta_init.age - 15)/5)


    # Upon initilisation, make sure there's somebody in each age group
    if n_generate == 'initilise':
        n_oldest = sum(meta_init.age_group == 4)
        if n_oldest == 0:


            # If not, pick somebody at random and put them into that age group
            i = round(float(n_generate) * np.random.random(1)[0])
            meta_init.loc[i, 'age'] = 35.5
            meta_init.loc[i, 'age_group'] = 4


    # Lump the 4th age group (35-39) in with the 3rd age group (30-34) anyway
    meta_init.loc[meta_init.age_group == 4, 'age_group'] = 3


    return meta_init


#%% HELPER initilise_orientation()
#
#
# Sets the orientation in meta_init
#
#
def initilise_orientation(meta_init, orientation_dist):


    # Individuals ages 16-19
    meta_init.loc[(meta_init.age_group == 0) & (meta_init.orientation <= orientation_dist.hetero.iloc[0]), 'orientation'] = 0
    meta_init.loc[(meta_init.age_group == 0) & (meta_init.orientation <= orientation_dist.homo.iloc[0]) & (meta_init.orientation != 0), 'orientation'] = 1
    meta_init.loc[(meta_init.age_group == 0) & (meta_init.orientation <= orientation_dist.bi.iloc[0]) & (meta_init.orientation != 0) & (meta_init.orientation != 1), 'orientation'] = 2


    # Individuals ages 20-24
    meta_init.loc[(meta_init.age_group == 1) & (meta_init.orientation <= orientation_dist.hetero.iloc[1]), 'orientation'] = 0
    meta_init.loc[(meta_init.age_group == 1) & (meta_init.orientation <= orientation_dist.homo.iloc[1]) & (meta_init.orientation != 0), 'orientation'] = 1
    meta_init.loc[(meta_init.age_group == 1) & (meta_init.orientation <= orientation_dist.bi.iloc[1]) & (meta_init.orientation != 0) & (meta_init.orientation != 1), 'orientation'] = 2


    # Individuals aged 25 and older
    meta_init.loc[(meta_init.age_group > 1) & (meta_init.orientation <= orientation_dist.hetero.iloc[2]), 'orientation'] = 0
    meta_init.loc[(meta_init.age_group > 1) & (meta_init.orientation <= orientation_dist.homo.iloc[2]) & (meta_init.orientation != 0), 'orientation'] = 1
    meta_init.loc[(meta_init.age_group > 1) & (meta_init.orientation <= orientation_dist.bi.iloc[2]) & (meta_init.orientation != 0) & (meta_init.orientation != 1), 'orientation'] = 2


    return meta_init


#%% HELPER initilise_risk()
#
#
# Initilises risk
#
#
def initilise_risk(meta_init, partners_dist):


    ## HETEROSEXUALS
    hetero = meta_init.orientation == 0


    # Individuals aged 16-19
    meta_init.loc[hetero & (meta_init.age_group == 0) & (meta_init.risk > partners_dist['16-19'].iloc[3]), 'risk'] = 0
    meta_init.loc[hetero & (meta_init.age_group == 0) & (meta_init.risk > 0), 'risk'] = 1


    # Individuals aged 20-24
    meta_init.loc[hetero & (meta_init.age_group == 1) & (meta_init.risk > partners_dist['20-24'].iloc[3]), 'risk'] = 0
    meta_init.loc[hetero & (meta_init.age_group == 1) & (meta_init.risk > 0), 'risk'] = 1


    # Individuals aged over 25
    meta_init.loc[hetero & (meta_init.age_group > 1) & (meta_init.risk > partners_dist['25-29'].iloc[3]), 'risk'] = 0
    meta_init.loc[hetero & (meta_init.age_group > 1) & (meta_init.risk > 0), 'risk'] = 1


    ## BI and HOMOSEXUALS
    non_hetero = hetero == False


    # Individuals aged 16-19
    meta_init.loc[non_hetero & (meta_init.age_group == 0) & (meta_init.risk > 4*partners_dist['16-19'].iloc[3]), 'risk'] = 0
    meta_init.loc[non_hetero & (meta_init.age_group == 0) & (meta_init.risk > 0), 'risk'] = 1


    # Individuals aged 20-24
    meta_init.loc[non_hetero & (meta_init.age_group == 1) & (meta_init.risk > 4*partners_dist['20-24'].iloc[3]), 'risk'] = 0
    meta_init.loc[non_hetero & (meta_init.age_group == 1) & (meta_init.risk > 0), 'risk'] = 1


    # Individuals aged over 25
    meta_init.loc[non_hetero & (meta_init.age_group > 1) & (meta_init.risk > 4*partners_dist['25-29'].iloc[3]), 'risk'] = 0
    meta_init.loc[non_hetero & (meta_init.age_group > 1) & (meta_init.risk > 0), 'risk'] = 1


    return meta_init


#%% HELPER initilise_sexual_acts()
#
#
# Initilises who does which sexual acts
#
#
def initilise_sexual_acts(meta_init, acts_dist):


    ## AGE GROUP 16-19
    who = meta_init.age_group == 0
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[0][1]), 'has_oral'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[1][1]), 'has_sex'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[2][1]), 'has_anal'] = True
    
    
    ## AGE GROUP 20-24
    who = meta_init.age_group == 1
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[0][2]), 'has_oral'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[1][2]), 'has_sex'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[2][2]), 'has_anal'] = True
    
    
    ## AGE GROUP 25-29
    who = meta_init.age_group == 2
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[0][3]), 'has_oral'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[1][3]), 'has_sex'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[2][3]), 'has_anal'] = True
    
    
    ## AGE GROUP 30-34
    who = meta_init.age_group == 3
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[0][4]), 'has_oral'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[1][4]), 'has_sex'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[2][4]), 'has_anal'] = True
    
    
    ## AGE GROUP 35
    who = meta_init.age_group == 4
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[0][5]), 'has_oral'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[1][5]), 'has_sex'] = True
    meta_init.loc[who & (np.random.random(len(who)) < acts_dist.values[2][5]), 'has_anal'] = True


    return meta_init


#%% HELPER initilise_infections()
#
#
#
#
#
def initilise_infections(meta_init, prop_infected, t):


    # Choose people at random to infect
    infected = meta_init.treatment_threshold < prop_infected
    meta_init.loc[infected, 'state'] = 'E'


    # Choose one site of infection for these individuals
    site = np.random.choice([0, 1, 2], len(meta_init))
    site0 = infected & (site == 0)
    site1 = infected & (site == 1)
    site2 = infected & (site == 2)


    # Set infection status
    meta_init.loc[site0, 'site0'] = 1
    meta_init.loc[site1, 'site1'] = 1
    meta_init.loc[site2, 'site2'] = 1


    # Set the time of infection to the end of the burn in period
    meta_init.loc[site0, 'site0_t0'] = t
    meta_init.loc[site1, 'site1_t0'] = t
    meta_init.loc[site2, 'site2_t0'] = t


    # Set the duration of infection
    latent = sim_parameters.init_duration_exposed[0] * np.random.random(len(meta_init))
    meta_init.loc[site0, 'site0_t1'] = t + latent[site0]
    meta_init.loc[site1, 'site1_t1'] = t + latent[site1]
    meta_init.loc[site2, 'site2_t1'] = t + latent[site2]


    # Set symptomatic status
    meta_init.loc[site0 & (meta_init.treatment_threshold[site0] < sim_parameters.init_prob_exposed[0]), 'site0_symptoms'] = True
    meta_init.loc[site1 & (meta_init.treatment_threshold[site1] < sim_parameters.init_prob_exposed[0]), 'site1_symptoms'] = True
    meta_init.loc[site2 & (meta_init.treatment_threshold[site2] < sim_parameters.init_prob_exposed[0]), 'site2_symptoms'] = True


    return meta_init


#%% HELPER cleanup_dypes()
#
#
# Some dypes get confused during setup and need to be remade
#
#
def cleanup_dtypes(meta_init):


    # Set variable types
    meta_init.loc[:,'gender'] = meta_init.gender.astype("int64")
    meta_init.loc[:,'age'] = meta_init.age.astype("float64")
    meta_init.loc[:,'age_group'] = meta_init.age_group.astype("int64")
    meta_init.loc[:,'orientation'] = meta_init.orientation.astype("int64")
    meta_init.loc[:,'risk'] = meta_init.risk.astype("int64")
    meta_init.loc[:,'partner'] = meta_init.partner.astype("int64")
    meta_init.loc[:,'counter'] = meta_init.counter.astype("int64")
    meta_init.loc[:,'has_oral'] = meta_init.has_oral.astype("bool")
    meta_init.loc[:,'has_anal'] = meta_init.has_anal.astype("bool")
    meta_init.loc[:,'has_sex'] = meta_init.has_sex.astype("bool")
    meta_init.loc[:,'state'] = meta_init.state.astype(CategoricalDtype(categories=['S', 'E', 'I', 'R', 'V', 'T'], ordered=True))
    meta_init.loc[:,'site0'] = meta_init.site0.astype("int64")
    meta_init.loc[:,'site1'] = meta_init.site1.astype("int64")
    meta_init.loc[:,'site2'] = meta_init.site2.astype("int64")
    meta_init.loc[:,'site0_t0'] = meta_init.site0_t0.astype("float64")
    meta_init.loc[:,'site1_t0'] = meta_init.site1_t0.astype("float64")
    meta_init.loc[:,'site2_t0'] = meta_init.site2_t0.astype("float64")
    meta_init.loc[:,'site0_t1'] = meta_init.site0_t1.astype("float64")
    meta_init.loc[:,'site1_t1'] = meta_init.site1_t1.astype("float64")
    meta_init.loc[:,'site2_t1'] = meta_init.site2_t1.astype("float64")
    meta_init.loc[:,'site0_symptoms'] = meta_init.site0_symptoms.astype('bool')
    meta_init.loc[:,'site1_symptoms'] = meta_init.site1_symptoms.astype('bool')
    meta_init.loc[:,'site2_symptoms'] = meta_init.site2_symptoms.astype('bool')
    meta_init.loc[:,'treatment_threshold'] = meta_init.treatment_threshold.astype("float64")
    meta_init.loc[:,'recovery_time'] = meta_init.recovery_time.astype("float64")
    meta_init.loc[:,'import_time'] = meta_init.import_time.astype("float64")


    return meta_init


#%% GRAPH graph_population()
#
#
#
#
#
def graph_population(pop_parameters, meta_init, save_dir):


    # Setup demographic data
    n = pop_parameters['n']
    sex_dist = pop_parameters['sex_dist']
    age_dist = pop_parameters['age_dist']
    orientation_dist = pop_parameters['orientation_dist']
    partners_dist = pop_parameters['partners_dist']


    #% Sex


    # Sex distribution
    simulated = [np.sum( (meta_init.gender == 1) ), np.sum( (meta_init.gender == 0) )]
    target = [n * sex_dist.pMale.iloc[0], n * (1 - sex_dist.pMale.iloc[0])]
    labels = ['Males', 'Females']
    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots()
    ax.bar(x-width/2, simulated, width, label = "Simulated")
    ax.bar(x+width/2, target, width, label = "Target")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    ax.set_ylabel('Number of Individuals')
    ax.set_title('Sex Distribution: Simulated vs. Target')
    ax.legend()
    plt.savefig(save_dir + '_sex.png', bbox_inches='tight')
    plt.close()
    #plt.show()


    #% Sexual orientation by age group


    # Setup the simulated data
    meta_init.loc[meta_init.age > 35, 'age_group'] = 4
    counted = meta_init.pivot_table(index='age_group', columns='orientation', fill_value=0, aggfunc='count')['age'].unstack()
    hetero = counted[0]
    homo = counted[1]
    bi = counted[2]
    labels = ['[16, 20)', '[20, 25)', '[25, 30)', '[30, 35)', '[35, 36)']
    x = np.arange(len(labels))

    # Setup the target data
    target_dist = orientation_dist.append(orientation_dist.loc[3,:], ignore_index=True)
    target_dist = target_dist.drop_duplicates()
    target_dist.age_group.iloc[4] = 4
    target_dist.bi = target_dist.bi - target_dist.homo
    target_dist.homo = target_dist.homo - target_dist.hetero
    target_dist.loc[0,:] = n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[5] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[0]) * target_dist.loc[0,:]
    target_dist.loc[1,:] = n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[6] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[1]) * target_dist.loc[1,:]
    target_dist.loc[2,:] = n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[7] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[2]) * target_dist.loc[2,:]
    target_dist.loc[3,:] = n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[8] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[3]) * target_dist.loc[3,:]
    target_dist.loc[4,:] = n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[9] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[4]) * target_dist.loc[4,:]


    # Make figure
    fig, ax = plt.subplots()
    ax.bar(x-width/2, bi, width, label = 'Simulated', color = 'tab:green')
    ax.bar(x-width/2, homo, width, bottom = bi, label = 'Simulated', color = 'tab:orange')
    ax.bar(x-width/2, hetero, width, bottom = bi + homo, label = 'Simulated', color = 'tab:blue')
    ax.bar(x+width/2, target_dist.bi, width, label = 'Target', color = 'tab:green', alpha = 0.6)
    ax.bar(x+width/2, target_dist.homo, width, bottom = target_dist.bi, label = 'Target', color = 'tab:orange', alpha = 0.6)
    ax.bar(x+width/2, target_dist.hetero, width, bottom = target_dist.bi + target_dist.homo, label = 'Target', color = 'tab:blue', alpha = 0.6)

    ax.set_xlabel('Age Group')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of Individuals')
    # ax.set_title('Distribution of Age by Sexual Orientation: Simulated (left) vs Target (right)')
    ax.legend(handles = [mpatches.Patch(color='tab:blue', label = "Heterosexual"),
                         mpatches.Patch(color='tab:orange', label = "Homosexual"),
                         mpatches.Patch(color='tab:green', label = "Bisexual")])
    plt.savefig(save_dir + '_age_and_orientation.png', bbox_inches='tight')
    plt.close()
    #plt.show()
    
    
    #% Sexual acts
    oral_y = [sum((meta_init.age_group == 0) & (meta_init.has_oral)),
              sum((meta_init.age_group == 1) & (meta_init.has_oral)),
              sum((meta_init.age_group == 2) & (meta_init.has_oral)),
              sum((meta_init.age_group == 3) & (meta_init.has_oral)),
              sum((meta_init.age_group == 4) & (meta_init.has_oral))]
    anal_y = [sum((meta_init.age_group == 0) & (meta_init.has_anal)),
              sum((meta_init.age_group == 1) & (meta_init.has_anal)),
              sum((meta_init.age_group == 2) & (meta_init.has_anal)),
              sum((meta_init.age_group == 3) & (meta_init.has_anal)),
              sum((meta_init.age_group == 4) & (meta_init.has_anal))]
    sex_y = [sum((meta_init.age_group == 0) & (meta_init.has_sex)),
             sum((meta_init.age_group == 1) & (meta_init.has_sex)),
             sum((meta_init.age_group == 2) & (meta_init.has_sex)),
             sum((meta_init.age_group == 3) & (meta_init.has_sex)),
             sum((meta_init.age_group == 4) & (meta_init.has_sex))]
    tot = [sum((meta_init.age_group == 0)),
           sum((meta_init.age_group == 1)),
           sum((meta_init.age_group == 2)),
           sum((meta_init.age_group == 3)),
           sum((meta_init.age_group == 4))]
    
    # Make figure
    fig, ax = plt.subplots(3)
    ax[0].bar(x-width/2, oral_y, width, label = 'Simulated')
    ax[1].bar(x-width/2, anal_y, width, label = 'Simulated')
    ax[2].bar(x-width/2, sex_y, width, label = 'Simulated')
    ax[0].bar(x+width/2, np.multiply(pop_parameters['sexual_acts_dist'].values[0][1:], np.array(tot)), width, label = 'Target')
    ax[1].bar(x+width/2, np.multiply(pop_parameters['sexual_acts_dist'].values[2][1:], np.array(tot)), width, label = 'Target')
    ax[2].bar(x+width/2, np.multiply(pop_parameters['sexual_acts_dist'].values[1][1:], np.array(tot)), width, label = 'Target')
    
    ax[2].set_xlabel('Age Group')
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(labels)
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    
    ax[0].set_ylabel('Number of Individuals')
    ax[1].set_ylabel('Number of Individuals')
    ax[2].set_ylabel('Number of Individuals')
    
    ax[0].set_title('People Who Engage in Oral')
    ax[1].set_title('People Who Engage in Anal')
    ax[2].set_title('People Who Engage in Sex')
    # ax.set_title('Distribution of Age by Sexual Orientation: Simulated (left) vs Target (right)')
    ax[2].legend()
    plt.savefig(save_dir + '_sexual_acts.png', bbox_inches='tight')
    plt.close()
    #plt.show()
    


    #% Risk level by age group

    counted = meta_init.pivot_table(index='age_group', columns='risk', fill_value=0, aggfunc='count')['age'].unstack()
    low = counted[0]
    high = counted[1]


    # Setup the simulated data
    #low = meta_init.loc[meta_init.risk == 0, :]. groupby('age_group')['age_group'].count()
    #high = meta_init.loc[meta_init.risk == 1, :]. groupby('age_group')['age_group'].count()
    labels = ['[16, 20)', '[20, 25)', '[25, 30)', '[30, 35)', '[35, 36)']
    x = np.arange(len(labels))


    # Setup for working out the partners_dist target
    risk = partners_dist.iloc[3, :]


    # Setup the target data
    target_dist = partners_dist.append(partners_dist.loc[3,:], ignore_index=True)
    target_dist = target_dist.drop_duplicates()
    target_high = [n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[5] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[0]) * risk[0],
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[6] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[1]) * risk[1],
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[7] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[2]) * risk[2],
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[8] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[3]) * risk[3],
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[9] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[4]) * risk[3]]
    target_low =  [n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[5] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[0]) * (1-risk[0]),
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[6] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[1]) * (1-risk[1]),
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[7] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[2]) * (1-risk[2]),
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[8] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[3]) * (1-risk[3]),
                   n * (sex_dist.pMale.iloc[0] * age_dist.ave.iloc[9] + (1-sex_dist.pMale.iloc[0]) * age_dist.ave.iloc[4]) * (1-risk[3])]


    # Make figure
    fig, ax = plt.subplots()
    ax.bar(x-width/2, high, width, label = 'Simulated', color = 'tab:orange')
    ax.bar(x-width/2, low, width, bottom = high, label = 'Simulated', color = 'tab:blue')
    ax.bar(x+width/2, target_high, width, label = 'Target', color = 'tab:orange', alpha = 0.6)
    ax.bar(x+width/2, target_low, width, bottom = target_high, label = 'Target', color = 'tab:blue', alpha = 0.6)

    ax.set_xlabel('Age Group')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of Individuals')
    # ax.set_title('Distribution of Age by Sexual Orientation: Simulated (left) vs Target (right)')
    ax.legend(handles = [mpatches.Patch(color='tab:blue', label = "Low Infection-risk"),
                         mpatches.Patch(color='tab:orange', label = "High Infection-risk")])
    plt.savefig(save_dir + '_age_and_risk.png', bbox_inches='tight')
    plt.close()
    #plt.show()