# -*- coding: utf-8 -*-
"""
Created on Wed May  5 12:59:05 2021

@author: Nicolas Rebuli

This script will generate a number of synthetic populations and then run
the partnership dynamics for a burn in period before saving the output.

The parameters are as follows:
    param.n_populations = the number of populations to simulate
    param.partner_burn_in = the number of days to run the partnership dynamics for
"""


#%% SETUP Load Modules


# Load libraries
import numpy as np
import pandas as pd
# import feather
# import matplotlib.pyplot as plt
import tqdm as tqdm
import matplotlib.pyplot as plt


# For clearing out directories
import os
import glob


# Parallel processing
from multiprocessing import Pool


# Load modules for simulation script
import src.demographic.generate_population as pop
import src.partners.partners as prt
import src.demographic.population_dynamics as demo
import src.calibration.setup as setup
# import src.infections.ng as ng
# import src.treatment.simple as trt


# Read in simulation parameters
param = pd.read_csv("data/param.csv")


# How many cores to run on
n_cores = 25


#%% SETUP Script parameters


# Which things should be generated here
generate_populations = True
burn_in_partnerships = True
track_partnership_rates = True
generate_parameters = True

# Set whether or not you want to overwrite the existing data
recovery_mode = True


#%% RUN Generate Populations


# Continue if populations have been asked to be regenrated
if generate_populations:


    # Print an update
    print('Generating ' + str(param.n_populations[0]) + ' synthetic populations for each scenario\n', flush = True)


    # Purge the directory if not in recovery mode
    if (__name__ == '__main__') & (recovery_mode == False):


        # Identify all existing population data
        files = \
            glob.glob('simulations/populations/scenario_1/*') + \
            glob.glob('simulations/populations/scenario_2/*') + \
            glob.glob('simulations/populations/scenario_3/*')


        # Delete them all
        for f in files:
            os.remove(f)


    # Define function for generating the requested population data
    def generate_population_fun(i):
        print('Generating population ' + str(i), flush=True)


        # Iterate over the scenario number
        for scenario in [1, 2, 3]:


            # Parse demographic parameters for population
            pop_parameters = pop.setup_data(scenario, 'parallel')


            # Check if the file is there already or not
            file_name = 'simulations/populations/scenario_' + str(scenario) + '/population_' + str(i)
            if os.path.exists(file_name + '.ftr') == False:


                # Generate population
                meta = pop.generate_population(pop_parameters)


                # Store population
                meta.to_feather(file_name + '.ftr')


                # Make graphs of the population
                pop.graph_population(pop_parameters, meta, file_name)


    # Define function for handling the parallel pool
    def pool_handler():
        p = Pool(n_cores)
        p.map(generate_population_fun, range(0, param.n_populations[0]))


    # Run generate_population() function in parallel
    if __name__ == '__main__':
        pool_handler()


#%% RUN Burn in Partnership Networks


# Continue if this has been asked for
if burn_in_partnerships:


    # Purge the directory if not in recovery mode
    if (__name__ == '__main__') & (recovery_mode == False):


        # Identify all existing partnership data
        files = \
            glob.glob('simulations/partnerships/scenario_1/*') + \
            glob.glob('simulations/partnerships/scenario_2/*') + \
            glob.glob('simulations/partnerships/scenario_3/*')


        # Delete them all
        for f in files:
            os.remove(f)


    # Define function for burning in the population
    def burn_in_partnerships_fun(i):


        # Iterate over the scenario number
        for scenario in [1, 2, 3]:


            # Parse parameters for population
            pop_parameters = pop.setup_data(scenario, 'parallel')
            prt_parameters = prt.setup_data(pop_parameters, 'parallel')
            inf_parameters = setup.setup_transmission_parameters(set = 'default', run_mode = 'parallel')


            # Iterate over the simulation number
            print('Burn in for population ' + str(i) + ' of scenario ' + str(scenario))


            # Skip file if it is already there
            save_dir = 'simulations/partnerships/scenario_' + str(scenario) + '/population_' + str(i)
            if os.path.exists(save_dir + '_matrix.npy') == False:


                # Read in population data
                n_days = param.partner_burn_in[0]
                partner_matrix = pop.initilise_partner_matrix(pop_parameters)
                partner_expire = pop.initilise_partner_duration(pop_parameters)
                file_name_pop = 'simulations/populations/scenario_' + str(scenario) + '/population_' + str(i) + '.ftr'
                meta = pd.read_feather(file_name_pop)


                # Setup population mobility dynamics
                pop_parameters = demo.initilise_demographic_dynamics(pop_parameters, inf_parameters, meta, partner_matrix, 0)
                

                # Initlise variables for making graphs
                if track_partnership_rates:
                    partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.initilise_partnership_trackers(n_days)
                    compartments, import_count, export_count, demographic, import_time, infections, active_age = demo.initilise_demographic_trackers(n_days)


                # Run Partnership Dynamics
                for t in range(0, n_days):
                    meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, inf_parameters, meta, partner_matrix, partner_expire)
                    meta, partner_matrix, partner_expire = prt.update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire)
                    if track_partnership_rates:
                        compartments, import_count, demographic, import_time, infections, export_count, active_age = demo.update_demographic_tracker(t, pop_parameters, meta, compartments, import_count, demographic, import_time, infections, export_count, active_age)
                        partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.update_partnership_trackers(t, pop_parameters, meta, partner_matrix, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot)



                # Store data for later
                meta.to_feather(save_dir + '_meta.ftr')
                np.save(save_dir + '_matrix.npy', partner_matrix)
                np.save(save_dir + '_expire.npy', partner_expire)
                np.save(save_dir + '_partners.npy', prt_parameters['prt_calib_set'])


                # Graph Partnership dynamics
                if track_partnership_rates:
                    tt = list(range(0, n_days))
                    demo.make_demographic_graphs(tt, pop_parameters, compartments, import_count, demographic, import_time, infections, export_count, active_age, save_dir)
                    prt.make_partnership_graphs(tt, pop_parameters, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot, save_dir)


    # Define function for handling the parallel pool
    def pool_handler():
        p = Pool(n_cores)
        p.map(burn_in_partnerships_fun, range(0, param.n_populations[0]))


    # Run generate_population() function in parallel
    if __name__ == '__main__':
        pool_handler()


#%% RUN Generate parameters


# Only run if asked to
if (__name__ == '__main__') & (generate_parameters == True):


    # How many parameter sets to simulate
    n_sim = param['n_parameter_sets'][0]


    # Probability that an imported individual is infectious
    # import_prob = np.random.beta(10, 40, n_sim)
    # import_prob = 0.5
    # plt.hist(import_prob)


    # Latent period
    # latent_period = np.round(14 * np.random.beta(2, 5, n_sim))
    # latent_period = 4
    # plt.hist(latent_period)


    # Duration of natural infection
    #https://pubmed.ncbi.nlm.nih.gov/26886136/#:~:text=Pharyngeal%20gonorrhoea%20(114%2D138%20days,days)%20compared%20with%20gonorrhoea%20infection.
    # mean_rectal = 360
    # var_rectal = 35
    # symptoms_rectal = 0.15
    # mean_ural_male = 185
    # var_ural_male = 35
    # symptoms_ural_male = 0.9
    # mean_ural_female = 185
    # var_ural_female = 35
    # symptoms_ural_female = 0.11
    # mean_phar = 84
    # var_phar = 21
    # symptoms_phar = 0
    

    # Probability of seeking treatment given a rectal infection
    # Assumed to be the same distribution for a female urogenital infection
    # m = 0.05
    # v = 0.001
    # alpha = ((m*(1-m)/v - 1)) * m
    # beta = ((m*(1-m)/v - 1)) * (1 - m)
    # symptoms_rectal = np.random.beta(alpha, beta, n_sim)
    # symptoms_ural_female = np.random.beta(alpha, beta, n_sim)
    symptoms_rectal = np.random.random(n_sim)
    symptoms_ural_female = np.random.random(n_sim)
    
    
    # Probability of a male seeking treatment given a urogenetal infection
    # m = 0.09
    # v = 0.002
    # alpha = ((m*(1-m)/v - 1)) * m
    # beta = ((m*(1-m)/v - 1)) * (1 - m)
    # symptoms_ural_male = np.random.beta(alpha, beta, n_sim)
    symptoms_ural_male = np.random.random(n_sim)
    

    # Probabilities of transmission given a sexual event
    pru = np.random.random(n_sim)
    prp = np.random.random(n_sim)
    pur = np.random.random(n_sim)
    puu = np.random.random(n_sim)
    pup = np.random.random(n_sim)
    ppr = np.random.random(n_sim)
    ppu = np.random.random(n_sim)
    ppp = np.random.random(n_sim)


    # Probabilities of engaging in different acts
    # p_kiss = np.random.random(n_sim)
    # p_oral_MF = np.random.random(n_sim)
    # p_oral_MM = np.random.random(n_sim)
    # p_oral_FF = np.random.random(n_sim)
    # p_oral_FM = np.random.random(n_sim)
    # p_sex_MF = np.random.random(n_sim)
    # p_sex_MM = np.random.random(n_sim)
    # p_sex_FF = np.random.random(n_sim)
    # p_anal_MF = np.random.random(n_sim)
    # p_anal_MM = np.random.random(n_sim)
    # p_rim = np.random.random(n_sim)


    # Parameters regarding the likelihood of seeking treatment
    # treat_mean = (21 * np.random.beta(10, 5, n_sim))
    # plt.hist(treat_mean)
    # treat_var = 2.2


    # Parameters around immunity from treatment
    # immune_mean	= 7 * np.random.beta(10, 5, n_sim)
    # plt.hist(immune_mean)
    # immune_var = 2


    # Construct parameter set
    sim_parameters = pd.DataFrame({ # 'prob_import_infectious': import_prob,
                                    # 'latent_period': latent_period,
                                    # 'mean_rectal': mean_rectal,
                                    # 'var_rectal': var_rectal,
                                    'symptoms_rectal': symptoms_rectal,
                                    # 'mean_ural_male': mean_ural_male,
                                    # 'var_ural_male': var_ural_male,
                                    'symptoms_ural_male': symptoms_ural_male,
                                    # 'mean_ural_female': mean_ural_female,
                                    # 'var_ural_female': var_ural_female,
                                    'symptoms_ural_female': symptoms_ural_female,
                                    # 'mean_phar': mean_phar,
                                    # 'var_phar': var_phar,
                                    # 'symptoms_phar': symptoms_phar,
                                    'pru': pru,
                                    'prp': prp,
                                    "pur": pur,
                                    'puu': puu,
                                    'pup': pup,
                                    'ppr': ppr,
                                    'ppu': ppu,
                                    'ppp': ppp,
                                    # 'p_kiss': p_kiss,
                                    # 'p_oral_MF': p_oral_MF,
                                    # 'p_oral_MM': p_oral_MM,
                                    # 'p_oral_FF': p_oral_FF,
                                    # 'p_oral_FM': p_oral_FM,
                                    # 'p_sex_MF': p_sex_MF,
                                    # 'p_sex_MM': p_sex_MM,
                                    # 'p_sex_FF': p_sex_FF,
                                    # 'p_anal_MF': p_anal_MF,
                                    # 'p_anal_MM': p_anal_MM,
                                    # 'p_rim': p_rim,
                                    # 'treat_mean': treat_mean,
                                    # 'treat_var': treat_var,
                                    # 'immune_mean': immune_mean,
                                    # 'immune_var': immune_var
                                    })


    # Check to see if there is an old version lying around
    current_version = 'simulations/parameters.csv'
    if os.path.exists(current_version) == True:


        # If so, rename the file
        for v in range(0, 100):
            old_version = 'simulations/parameters-old_version_' + str(v) + '.csv'
            if os.path.exists(old_version) == False:
                os.rename(current_version, old_version)
                break


    # Save the new version of the parameter set
    sim_parameters.to_csv(current_version, index = False)











#%%  TEST DISTRIBUTIONS


# # Distribution for probability of testing rectal infection
# m = 0.05
# v = 0.001

# alpha = ((m*(1-m)/v - 1)) * m
# beta = ((m*(1-m)/v - 1)) * (1 - m)

# x = np.random.beta(alpha, beta, 10000)
# plt.hist(x)


# # Distribution for probability of testing rectal infection
# m = 0.09
# v = 0.002

# alpha = ((m*(1-m)/v - 1)) * m
# beta = ((m*(1-m)/v - 1)) * (1 - m)

# x = np.random.beta(alpha, beta, 10000)
# plt.hist(x)










