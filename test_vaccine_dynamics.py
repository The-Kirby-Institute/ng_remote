# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:21:09 2022

@author: Nicolas Rebuli


Test script for the development of the vaccine dynamics implemented
inplimented in the model.
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
import src.vaccinations.vaccination as vax


make_graphs = True
run_mode = 'serial'
scenario = 1


#%% SETUP Parse simulation data
print('PARSING PARAMETERS\n------------------')


# Setup simulation data
sim_parameters, pop_parameters, prt_parameters, inf_parameters, meta, partner_expire, partner_matrix, population_no = \
    setup.setup_data(scenario = scenario,
                     run_mode = run_mode,
                     inf_param_set = 'calibration')


#%% SETUP Seed infections to ensure transmission


# Set star time for simulation
t0_sim = sim_parameters.partner_burn_in[0]


# Just seed a tonne of infections
site0 = np.random.random(len(meta)) < 0.25
site1 = np.random.random(len(meta)) < 0.25
site2 = np.random.random(len(meta)) < 0.25
# infected = np.random.random(len(meta)) <= sim_parameters['init_prob_exposed'][0]
infected = np.random.random(len(meta)) <= 1
infected = (meta.risk == 1) & infected & (site0 | site1 | site2)
meta.loc[infected, 'state'] = 'I'
meta.loc[infected & site0, 'site0'] = int(1)
meta.loc[infected & site1, 'site1'] = int(1)
meta.loc[infected & site2, 'site2'] = int(1)
meta.loc[(infected) & (site0), 'site0_t0'] = t0_sim
meta.loc[(infected) & (site1), 'site1_t0'] = t0_sim
meta.loc[(infected) & (site2), 'site2_t0'] = t0_sim
meta.loc[(infected) & (site0), 'site0_t1'] = t0_sim + ng.duration_rectal(inf_parameters, meta.loc[(infected) & (site0), ], np.sum((infected) & (site0)))
meta.loc[(infected) & (site1), 'site1_t1'] = t0_sim + ng.duration_urethral(inf_parameters, meta.loc[(infected) & (site1), ], np.sum((infected) & (site1)))
meta.loc[(infected) & (site2), 'site2_t1'] = t0_sim + ng.duration_pharyngeal(inf_parameters, meta.loc[(infected) & (site2), ], np.sum((infected) & (site2)))



#%% SETUP Initilise parameters for vaccination rollout


# Deployment modes:
#  0 vaccinations at the same time as treatments
#  1 everyone vaccinated at the age of 16
#  2 a proportion of each age group is vaccinated


# Effectiveness
#  0 complete immunity for the duration of the vaccine
#  1 transmission probability is reduced
#  2 the probability of being asymptomatic is increased
#  3 the duration of infection is decreased
#  4 a combination of all of the above


# Setup parameters for vaccination
vax_parameters = {'duration_mean': 200*365,
                  'duration_var': 1,
                  'prop_effective': 0,      # Probability of vaccination working
                  'effect': 0,
                  'site0_trans_reduce': 0,   # Scaling of the transmission probability
                  'site1_trans_reduce': 0,   # Scaling of the transmission probability
                  'site2_trans_reduce': 0,   # Scaling of the transmission probability
                  'site0_symp_reduce': 0.5,    # Scaling of the probability of symptoms
                  'site1_symp_reduce': 0.5,    # Scaling of the probability of symptoms
                  'site2_symp_reduce': 0.5,    # Scaling of the probability of symptoms
                  'site0_duration_reduce': 0.5, # Scaling of the duration of inection
                  'site1_duration_reduce': 0.5, # Scaling of the duration of inection
                  'site2_duration_reduce': 0.5, # Scaling of the duration of inection
                  'deployment': 0,
                  'prop_vax_0': 0.5,  # Proportion of 16-19 year olds vaccinated
                  'prop_vax_1': 0.4,  # Proportion of 20-24 year olds vaccinated
                  'prop_vax_2': 0.1,  # Proportion of 25-29 year olds vaccinated
                  'prop_vax_3': 0.0}  # Proportion of over 30 year olds vaccinated


# Setup method for updating infections
update_infections = vax.set_function_for_updating_infections(vax_parameters)



#%% RUN Simulation


# Setup time variables
# n_steps = sim_parameters.simulation_length[0]
n_steps = int(30*365)
tt = range(t0_sim, t0_sim + n_steps)


track_vax = np.zeros((n_steps, 6*2))


# Initlise variables for making graphs
if make_graphs:
    partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.initilise_partnership_trackers(n_steps)
    compartments_demo, import_count, export_count, demographic, import_time, active_age = demo.initilise_demographic_trackers(n_steps, t0_sim = t0_sim)
    inf_tracker = ng.initilise_infection_trackers(n_steps)


# Record code runtime
t0 = time.time()


# Check to see if this dataset has been run to completion
( print('Running simulation...\n', flush = True) if run_mode == 'serial' else [] )
for t in tqdm.tqdm(tt):


    # Update demographic
    meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, prt_parameters, inf_parameters, meta, partner_matrix, partner_expire)
    
    
    # Update partnerships
    meta, partner_matrix, partner_expire = prt.update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire)
    
    
    # Update infections
    meta = update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix)
    
    
    # Update vaccine tracker
    ii = 0
    for state in ['S', 'E', 'I', 'R', 'T', 'V']:
        track_vax[t-t0_sim, ii] = np.sum((meta.vaccinated) & (meta.state == state))
        ii = ii + 1
        
    for state in ['S', 'E', 'I', 'R', 'T', 'V']:
        track_vax[t-t0_sim, ii] = np.sum((~meta.vaccinated) & (meta.state == state))
        ii = ii + 1
    
    
    # Update trackers for graphs
    if make_graphs:
        pop_parameters = demo.update_attribute_tracker(pop_parameters, meta, partner_matrix)
        compartments_demo, import_count, demographic, import_time, export_count, active_age = demo.update_demographic_tracker(t - t0_sim, pop_parameters, meta, compartments_demo, import_count, demographic, import_time, export_count, active_age, t0_sim = t0_sim)
        partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.update_partnership_trackers(t - t0_sim, pop_parameters, meta, partner_matrix, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot)
        inf_tracker = ng.update_infection_trackers(t - t0_sim, pop_parameters, meta, inf_tracker)
    
    

        

# Print update
runtime = (time.time() - t0)/60
rt_min = int(np.floor(runtime))
rt_sec = int(60 * (runtime - rt_min))
print('\nRuntime: ' + str(rt_min) + ' min ' + str(rt_sec) + ' seconds')


#%% GRAPHS


if make_graphs:
    demo.make_demographic_graphs(tt, pop_parameters, compartments_demo, import_count, demographic, import_time, export_count, active_age)
    prt.make_partnership_graphs(tt, pop_parameters, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot)
    ng.make_infection_graphs(tt, inf_tracker, pop_parameters)


# # Setup graph
# cmap = plt.get_cmap("tab10")
# fig, ax = plt.subplots(1)

# # Put on high-level state numbers
# states = ['S', 'E', 'I', 'R', 'T', 'V']
# for ii in range(1, 6):
#     ax.plot(tt, track_vax[:, ii+6], label = states[ii] + '-NoVax', color = cmap(ii))
#     ax.plot(tt, track_vax[:, ii], label = states[ii] + '-Vax', color = cmap(ii), linestyle = '--')


# # Labels
# ax.set_title('Number of People in Each Infection State\nBy Whether They are Vaccinated')
# ax.legend(ncol = 5)
# plt.xlabel('Day')
# plt.ylabel('Number of People')

















