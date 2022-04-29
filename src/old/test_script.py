# -*- coding: utf-8 -*-
"""
Created on Fri May 21 16:27:20 2021

@author: nicol
"""


#%% SETUP Modules


# Standard modules
import numpy as np
import pandas as pd
import tqdm as tqdm
import matplotlib.pyplot as plt
import random
import pickle
from pathlib import Path
import shutil
import os
import time


# My modules
import src.demographic.generate_population as pop
import src.demographic.population_dynamics as demo
import src.partners.partners as prt
import src.calibration.setup as setup
import src.infections.ng as ng
import src.vaccinations.deployment_0 as vax


run_mode = 'serial'
scenario = 1
calibrated_no = 0


#%% SETUP Parse simulation data


# Parse infection parameters
# Pass this function which numbered parameter set you want to use
# Will default back to the baseline parameter set
inf_parameters = setup.parse_parameters('calibrated', scenario, calibrated_no)


# Read in data from the end of the simulation
output_dat = open('simulations/calibration/scenario_' + str(scenario) + '/simulation_' + str(inf_parameters['calibration_set']) + '/output.pkl', 'rb')
output = pickle.load(output_dat)
output_dat.close()


# Parse model and simulation parameters
parameter_no = output['parameter_no']
population_no = output['population_no']
inf_parameters = output['inf_parameters']
pop_parameters = output['pop_parameters']
sim_parameters = output['sim_parameters']


# Parse data from simulation
meta = output['meta']
partner_expire = output['partner_expire']
partner_matrix = output['partner_matrix']
sim_t0 = output['t']


# Print
if run_mode == 'serial':
    print('Parsing scenario ' + str(scenario) + ' with calibrated parameter set ' + str(calibrated_no) + '\n')
    print('Restoring previous simulation parameters\n----------------------------------------')
    print('Scenario: ' + str(scenario))
    print('Calibrated set: ' + str(calibrated_no))
    print('Parameter set: ' + str(parameter_no))
    print('Population set: ' + str(population_no))
    print('\n')


# Clean up workspace
del output, output_dat


#%% SETUP Modify meta to contain data on vaccination rollout


# Insert vaccination status
meta.insert(0, 'vaccinated', len(meta) * [False])


# Insert time of first vaccination
meta.insert(1, 'vaccination_t0', len(meta) * [float('inf')])


# Insert duration of vaccination
meta.insert(2, 'vaccination_t1', len(meta) * [float('inf')])


# Insert time of booster
meta.insert(3, 'booster_t0', len(meta) * [float('inf')])


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
vax_parameters = {'duration_mean': 2*365,
                  'duration_var': 1,
                  'prop_effective': 0.3,      # Probability of vaccination working
                  'effect': 1,
                  'site0_trans_reduce': 0.1,   # Scaling of the transmission probability
                  'site1_trans_reduce': 0.1,   # Scaling of the transmission probability
                  'site2_trans_reduce': 0.1,   # Scaling of the transmission probability
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


# Track partnership numbers
n_steps = 1000
tt = range(0, n_steps)
partners = np.zeros((n_steps, 4))
compartments = np.zeros((n_steps, 12))


# Time code
t0 = time.time()


# Check to see if this dataset has been run to completion
out_dir = 'simulations/calibration/scenario_' + str(scenario) +'/test'
( print('Running simulation...\n', flush = True) if run_mode == 'serial' else [] )
# for t in tqdm.tqdm(range(sim_parameters.partner_burn_in[0], sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0])):
# for t in range(sim_parameters.partner_burn_in[0], sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0]):
for t in tqdm.tqdm(tt):


    # Update population
    meta, partner_matrix, partner_expire = demo.update_population(pop_parameters, inf_parameters, meta, partner_matrix, partner_expire, sim_t0 + t)


    # Update partnerships
    meta, partner_matrix, partner_expire = prt.update_partnerships(meta, partner_matrix, partner_expire, sim_t0 + t)


    # Update infections
    meta = update_infections(inf_parameters, vax_parameters, meta, partner_matrix, sim_t0 + t)


    # Dump simulation output
    # meta.to_feather(out_dir + '/timestep' + str(t) + '.ftr')


    # How many partners do they have?
    partners[t, 0] = sum(meta.partner == -1)
    partners[t, 1] = sum(meta.partner > -1)
    partners[t, 2] = len(meta)
    partners[t, 3] = sum(sum(partner_matrix>0)) - sum(meta.partner > -1)


    # How many are in each compartment?
    compartments[t, 0] = sum(meta.state == 'S')
    compartments[t, 1] = sum(meta.state == 'E')
    compartments[t, 2] = sum(meta.state == 'I')
    compartments[t, 3] = sum(meta.state == 'R')
    compartments[t, 4] = sum(meta.state == 'V')
    compartments[t, 5] = sum(meta.state == 'T')
    compartments[t, 6] = sum((meta.site0 == 1) & (meta.site0_symptoms == True))
    compartments[t, 7] = sum((meta.site0 == 1) & (meta.site0_symptoms == False))
    compartments[t, 8] = sum((meta.site1 == 1) & (meta.site1_symptoms == True))
    compartments[t, 9] = sum((meta.site1 == 1) & (meta.site1_symptoms == False))
    compartments[t, 10] = sum((meta.site2 == 1) & (meta.site2_symptoms == True))
    compartments[t, 11] = sum((meta.site2 == 1) & (meta.site2_symptoms == False))


# Print update
runtime = (time.time() - t0)/60
print('\nRuntime: ' + str(runtime) + ' min')


#%% GRAPHS


# Check on aggregate partner numbers
plt.plot(tt, partners[:,0], label = 'single')
plt.plot(tt, partners[:,1], label = 'long-term')
plt.plot(tt, partners[:,3], label = 'short-term')
plt.plot(tt, partners[:,2], label = 'total')
plt.plot(tt, partners[:,0] + partners[:,1], label = 'total check')
plt.title('Partnership Dynamics')
plt.legend()
plt.show()


# Check on aggregate compartment numbers
# plt.plot(tt, compartments[:,0], label = 'S')
plt.plot(tt, compartments[:,1], label = 'E')
# plt.plot(tt, compartments[:,2], label = 'I')
plt.plot(tt, compartments[:,7] + compartments[:,9] + compartments[:,11], label = 'Ia')
plt.plot(tt, compartments[:,6] + compartments[:,8] + compartments[:,10], label = 'Is')
# plt.plot(tt, compartments[:,3], label = 'R')
plt.plot(tt, compartments[:,4], label = 'V')
plt.plot(tt, compartments[:,5], label = 'T')
# plt.plot(tt, partners[:,2], label = 'total')
# plt.plot(tt, np.sum(compartments, axis = 1), label = 'check')
plt.xlabel('Day')
plt.ylabel('Number of Individuals')
plt.title('One Simulation Run of the Proposed Vaccination Scenario')
plt.legend()
plt.show()