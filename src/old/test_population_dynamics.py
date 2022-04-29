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
import src.infections.ng as ng
import src.calibration.setup as setup


run_mode = 'serial'
scenario = 3
population_no = 0
t0_sim = 0


#%% SETUP Parse simulation data
print('PARSING PARAMETERS\n------------------')


# Read in simulation parameters
sim_parameters = pd.read_csv('data/param.csv')
    

# Parse demographic parameters
pop_parameters = pop.setup_data(scenario)


# Parse partnership parameters
prt_parameters = prt.setup_data(pop_parameters)


# Parse infection parameters
inf_parameters = setup.parse_default_transmission_parameters()


# Read in a population file
population_no = np.random.randint(0, sim_parameters.n_populations[0]-1)
meta, partner_expire, partner_matrix = setup.parse_population_data(scenario, population_no, run_mode)


# Preallocate importations
t0_sim = sim_parameters.partner_burn_in[0]
pop_parameters = demo.initilise_demographic_dynamics(pop_parameters, inf_parameters, meta, partner_matrix, t0_sim)


#%% RUN Simulation


# Track partnership numbers
n_years = 20
n_steps = n_years*365
tt = range(t0_sim, n_steps)


# Initlise variables for making graphs
compartments, import_count, export_count, demographic, import_time, active_age = demo.initilise_demographic_trackers(n_steps, t0_sim = t0_sim)


# Time code
t0 = time.time()


# Check to see if this dataset has been run to completion
( print('Running simulation...\n', flush = True) if run_mode == 'serial' else [] )
# for t in tqdm.tqdm(range(sim_parameters.partner_burn_in[0], sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0])):
# for t in range(sim_parameters.partner_burn_in[0], sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0]):
for t in tqdm.tqdm(tt):


    # Update population
    meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, prt_parameters, inf_parameters, meta, partner_matrix, partner_expire)


    # Update cohort trackers
    compartments, import_count, demographic, import_time, export_count, active_age = demo.update_demographic_tracker(t - t0_sim, pop_parameters, meta, compartments, import_count, demographic, import_time, export_count, active_age, t0_sim = t0_sim)


# Print update
runtime = (time.time() - t0)/60
rt_min = int(np.floor(runtime))
rt_sec = int(60 * (runtime - rt_min))
print('\nRuntime: ' + str(rt_min) + ' min ' + str(rt_sec) + ' seconds')


#%% GRAPHS


# Make graphs
demo.make_demographic_graphs(tt, pop_parameters, compartments, import_count, demographic, import_time, export_count, active_age)


# Any other graphs?












