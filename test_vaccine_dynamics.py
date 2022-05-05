# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:21:09 2022

@author: Nicolas Rebuli


Test script for the development of the vaccine dynamics implemented
inplimented in the model.
"""



#%% SETUP Modules


# Modules
from multiprocessing import Pool
import src.vaccinations.vaccination as vax


# Import vaccination parameters
# from simulations.vaccination.vaccination_scenarios.test_parameters import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_1_vax_distributed_no_effect import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_2_vax_means_no_infections import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_3_vax_means_infection_half_as_likely import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_4_vax_means_no_transmission import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_5_vax_means_half_transmission import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_6_vax_means_no_detection import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_7_vax_means_half_detection import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_8_vax_means_no_infectious_period import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_9_vax_means_half_infectious_period import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_10_vax_means_all_effects import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_11_vax_means_no_effects import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_12_vax_means_no_rectal import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_13_vax_means_no_pharyngeal import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_14_vax_means_no_urethral import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_15_vax_means_all_effects_at_debut import vax_parameters
# from simulations.vaccination.vaccination_scenarios.test_16_vax_means_all_effects_all_cohorts import vax_parameters
from simulations.vaccination.vaccination_scenarios.test_17_vax_everyone import vax_parameters






#%% RUN SIMULATIONS


# Set parameters
run_mode = 'serial'
n_cores = max(16, vax_parameters['n_sims'])
sim_run = list(range(0, vax_parameters['n_sims']))




# If running in serial
if run_mode == 'serial':
    
    
    # Run all the sims
    for sim_no in sim_run:
        vax.run_vaccine_scenario(vax_parameters, sim_no, run_mode = 'serial')




# If running in parallel
if run_mode == 'parallel':
    
    
    # Define worker for running the simulation
    def worker(it):
        vax.run_vaccine_scenario(vax_parameters, it, run_mode = 'parallel')
    
    
    # Define pool handler
    def pool_handler():
        p = Pool(n_cores)
        p.map(worker, sim_run)
    
    
    # Batch scripts
    if __name__ == '__main__':
        pool_handler()




















