# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 11:49:42 2021

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
scenario = 1


#%% SETUP Parse simulation data
print('PARSING PARAMETERS\n------------------')


# Read in simulation parameters
sim_parameters = pd.read_csv('data/param.csv')
    

# Parse demographic parameters
pop_parameters = pop.setup_data(scenario)


# Parse partnership parameters
prt_parameters = prt.setup_data(pop_parameters)


# Parse infection parameters
inf_parameters = setup.parse_parameters(scenario = scenario, run_mode = 'serial')


# Read in a population file
population_no = np.random.randint(0, sim_parameters.n_populations[0]-1)
meta, partner_expire, partner_matrix = setup.parse_population_data(scenario, population_no, run_mode)


# Preallocate importations
t0_sim = sim_parameters.partner_burn_in[0]
pop_parameters = demo.initilise_demographic_dynamics(pop_parameters, inf_parameters, meta, partner_matrix, t0_sim)


#%% RUN Simulation


# Setup time variables
n_steps = sim_parameters.simulation_length[0]
n_steps = 5*365
tt = range(t0_sim, t0_sim + n_steps)


# Initlise variables for making graphs
compartments, import_count, export_count, demographic, import_time, infections = demo.initilise_demographic_trackers(n_steps)
partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.initilise_partnership_trackers(n_steps)
yt, i0lt, i0ht, i1lt, i1ht, i2lt, i2ht, i3lt, i3ht = ng.initilise_infection_trackers(n_steps)


# Time code
t0 = time.time()


# Check to see if this dataset has been run to completion
( print('Running simulation...\n', flush = True) if run_mode == 'serial' else [] )
for t in tqdm.tqdm(tt):


    # Update demographic
    meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, inf_parameters, meta, partner_matrix, partner_expire)


    # Update partnerships
    meta, partner_matrix, partner_expire = prt.update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire)
    
    
    # Update infections
    meta = ng.update_infections(inf_parameters, meta, partner_matrix, t)
    

    # Update trackers for graphs
    compartments, import_count, demographic, import_time, infections, export_count = demo.update_demographic_tracker(t - t0_sim, pop_parameters, meta, compartments, import_count, demographic, import_time, infections, export_count)
    partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.update_partnership_trackers(t - t0_sim, pop_parameters, meta, partner_matrix, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot)
    yt, i0lt, i0ht, i1lt, i1ht, i2lt, i2ht, i3lt, i3ht = ng.update_infection_trackers(meta, t - t0_sim, yt, i0lt, i0ht, i1lt, i1ht, i2lt, i2ht, i3lt, i3ht)
    

# Print update
runtime = (time.time() - t0)/60
rt_min = int(np.floor(runtime))
rt_sec = int(60 * (runtime - rt_min))
print('\nRuntime: ' + str(rt_min) + ' min ' + str(rt_sec) + ' seconds')


#%% GRAPHS


# Demographic graphs
demo.make_demographic_graphs(tt, pop_parameters, compartments, import_count, demographic, import_time, infections, export_count)


# Partnership graphs
prt.make_partnership_graphs(tt, pop_parameters, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot)


# Infection graphs
ng.make_infection_graphs(yt, i0lt, i0ht, i1lt, i1ht, i2lt, i2ht, i3lt, i3ht, tt)










