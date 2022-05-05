# -*- coding: utf-8 -*-
"""
Created on Thu May 13 16:32:39 2021

@author: nicol
"""
#%% SETUP Modules


# Standard modules
import numpy as np
import pandas as pd
import tqdm as tqdm
import random
import shelve
from pathlib import Path
import shutil
import os


# Parallel processing
from multiprocessing import Pool


# My modules
import src.demographic.generate_population as pop
import src.demographic.population_dynamics as demo
import src.partners.partners as prt
import src.calibration.setup as calibrate
import src.infections.ng as ng


# Parameters
sim_parameters = pd.read_csv('data/param.csv')


#%% Run a test


# Just run one scenario
# calibrate.run_one_simulation(scenario = 1, parameter_no = 10)


# for i in range(0, sim_parameters['n_parameter_sets'][0]):
#     calibrate.run_one_simulation(scenario = 1, parameter_no = i)




#%% Run the script


# Identify the total number of parameter sets to run
sim_parameters = list(range(0, sim_parameters['n_parameter_sets'][0]))
# sim_parameters = [0]


# Define function for running the code
def worker(it):
    # calibrate.run_one_simulation(scenario = 1, parameter_no = it, run_mode = 'parallel')
    # calibrate.run_one_simulation(scenario = 2, parameter_no = it, run_mode = 'parallel')
    calibrate.run_one_simulation(scenario = 3, parameter_no = it, run_mode = 'parallel')


# Define function for handling the parallel pool
def pool_handler():
    p = Pool(16)
    p.map(worker, sim_parameters)


if __name__ == '__main__':
    pool_handler()


