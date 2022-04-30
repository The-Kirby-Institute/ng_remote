# -*- coding: utf-8 -*-
"""
Created on Fri May 21 09:19:55 2021

@author: nicol
"""

#%% SETUP Modules


# Standard modules
from pathlib import Path
import pandas as pd
import numpy as np
import random
import shelve
import shutil
import time
import tqdm
import os



# My modules
import src.demographic.generate_population as pop
import src.demographic.population_dynamics as demo
import src.partners.partners as prt
import src.calibration.setup as setup
import src.infections.ng as ng



run_mode = 'serial'
scenario = 1
parameter_no = 0
population_no = 0
# population_no = random.randint(0, sim_parameters.n_populations[0]-1)


#%% SETUP Import data


# Read in simulation parameters
( print('Parsing global simulation settings.\n') if run_mode == 'serial' else [] )
sim_parameters = pd.read_csv('data/param.csv')


# Parse infection parameters
# Pass this function which numbered parameter set you want to use
# Will default back to the baseline parameter set
inf_parameters = setup.parse_parameters('calibrated', scenario, parameter_no)


# Parse demographic parameters for population
# Pass this function the scenario number
# Will default back to scenario 0
pop_parameters = pop.setup_data(scenario, run_mode)


# Parse simulated population
# Pass this function the scenario number and which simulated population to use
# Will default back to scenario 1 and population 0 if none selected
meta, partner_expire, partner_matrix = setup.parse_population_data(scenario, population_no)


#%% RUN Simulation


# Check to see if this dataset has been run to completion
out_dir = 'simulations/vaccination/scenario_' + str(scenario) +'/simulation_' + str(parameter_no)
last_file = out_dir + '/timestep' + str(sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0] - 1) + '.ftr'
if os.path.exists(last_file) == False:


    #% SETUP Output Directory


    # Clean out directory from any previous results
    ( print('Purging simulation output directory.\n') if run_mode == 'serial' else [] )
    dirpath = Path(out_dir)
    if dirpath.exists() and dirpath.is_dir():
        shutil.rmtree(dirpath)


    # Remake the directory
    os.mkdir(dirpath)


    #% RUN


    # Time code
    t0 = time.time()


    ( print('Running simulation...\n', flush = True) if run_mode == 'serial' else [] )
    for t in tqdm.tqdm(range(sim_parameters.partner_burn_in[0], sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0])):
    # for t in range(sim_parameters.partner_burn_in[0], sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0]):


        # Update population
        meta, partner_matrix, partner_expire = demo.update_population(pop_parameters, inf_parameters, meta, partner_matrix, partner_expire, t)


        # Update partnerships
        meta, partner_matrix, partner_expire = prt.update_partnerships(meta, partner_matrix, partner_expire, t)


        # Update infections
        meta = ng.update_infections(inf_parameters, meta, partner_matrix, t)


        # Dump simulation output
        meta.to_feather(out_dir + '/timestep' + str(t) + '.ftr')


    ( print('\n\nSimulation complete!\n') if run_mode == 'serial' else [] )


    #%% SAVE OUTPUT


    # Shelve variables
    ( print('Shelving environment variables.\n') if run_mode == 'serial' else [] )
    my_shelf = shelve.open(out_dir + '_workspace.pkl', 'n')
    for k in dir():
        print(k)
        try:
            my_shelf[k] = globals()[k]
        except:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
            ( print('ERROR shelving: {0}'.format(key)) if run_mode == 'serial' else [] )
    my_shelf.close()
    ( print('\nComplete!\n\n\n\n') if run_mode == 'serial' else [] )


# Print update
# runtime = (time.time() - t0)/60
# ( print('Running parameter set: ' + str(parameter_no) + ' on scenario: ' + str(scenario) + ' with population: ' + str(population_no) + ' - runtime: ' + str(runtime) + ' min') if run_mode == 'parallel' else [] )
