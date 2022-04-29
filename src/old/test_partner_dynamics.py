# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 09:57:13 2021

@author: nicol
"""


#%% SETUP Modules


# Standard modules
import numpy as np
import pandas as pd
import tqdm as tqdm
import matplotlib.pyplot as plt
import time


# My modules
import src.demographic.generate_population as pop
import src.demographic.population_dynamics as demo
import src.partners.partners as prt
import src.calibration.setup as setup


run_mode = 'serial'
scenario = 1
population_no = 0
param_set = 'calibration'
t0_sim = 0


#%% SETUP Parse simulation data
print('PARSING PARAMETERS\n------------------')


# Parse population parameters
pop_parameters = pop.setup_data(scenario, run_mode)


# Read in population data
meta, partner_expire, partner_matrix = setup.parse_population_data(scenario, population_no, run_mode)


# Parse partnership parameters
prt_parameters = prt.setup_data(pop_parameters)


# Parse infection parameters
# Pass this function which numbered parameter set you want to use
# Will default back to the baseline parameter set
inf_parameters = setup.parse_parameters(param_set, scenario, run_mode = run_mode)


# Preallocate importations
pop_parameters = demo.initilise_demographic_dynamics(pop_parameters, inf_parameters, meta, partner_matrix, t0_sim)


#%% RUN Simulation


# Track partnership numbers
n_years = 10
n_steps = n_years*365 + 1
tt = range(t0_sim, n_steps)


# Initlise variables for making graphs
compartments, import_count, export_count, demographic, import_time, infections = demo.initilise_demographic_trackers(n_steps)
partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.initilise_partnership_trackers(n_steps)


# Time code
t0 = time.time()


# Check to see if this dataset has been run to completion
( print('Running simulation...\n', flush = True) if run_mode == 'serial' else [] )
for t in tqdm.tqdm(tt):


    # Update demographic
    meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, inf_parameters, meta, partner_matrix, partner_expire)


    # Update partnerships
    meta, partner_matrix, partner_expire = prt.update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire)


    # Update trackers for graphs
    compartments, import_count, demographic, import_time, infections, export_count = demo.update_demographic_tracker(t, pop_parameters, meta, compartments, import_count, demographic, import_time, infections, export_count)
    partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.update_partnership_trackers(t, pop_parameters, meta, partner_matrix, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot)


# Print update
runtime = (time.time() - t0)/60
rt_min = int(np.floor(runtime))
rt_sec = int(60 * (runtime - rt_min))
print('\nRuntime: ' + str(rt_min) + ' min ' + str(rt_sec) + ' seconds')


#%% GRAPHS


# Demographic graphs
demo.make_demographic_graphs(tt, pop_parameters, compartments, import_count, demographic, import_time, infections, export_count)


# Partnership graphs
prt.make_partnership_graphs(tt, pop_parameters, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot, show_graphs = True)












