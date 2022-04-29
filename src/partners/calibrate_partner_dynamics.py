#%% SETUP LIBRARIES


# Load standard modules
import os
import time
import numpy as np
import pandas as pd
import multiprocessing as mp


# Load modules for simulation script
import src.demographic.generate_population as pop
import src.partners.partners as prt
import src.demographic.population_dynamics as demo
import src.calibration.setup as setup


#%% SETUP WORKER FUNCTION FOR RUNNING IN PARALLEL


# Base simulation variables
scenario = 1
n_years = 3
n_days = 365 * n_years
n_cores = 25
sim_parameters = pd.read_csv("data/param.csv")
    

# Function
def worker(sim_no, calib_set, return_dict):
    print('\tsimulation number ' + str(sim_no))
    
    
    # Base simulation variables
    scenario = 1
    n_years = 3
    n_days = 365 * n_years
    n_cores = 25
    
    
    # Read in simulation parameters
    sim_parameters = pd.read_csv("data/param.csv")
    pop_parameters = pop.setup_data(scenario, 'parallel')
    prt_parameters = prt.setup_data(pop_parameters)
    inf_parameters = setup.parse_parameters('default', scenario)
    
    
    # Parse parameters for population
    pop_parameters = pop.setup_data(scenario, 'parallel')
    prt_parameters = prt.setup_data(pop_parameters, prt_calib_set = calib_set)
    inf_parameters = setup.parse_parameters('default', scenario)
    

    # Setup population data
    pop_no = np.random.choice(sim_parameters['n_populations'][0])
    partner_matrix = pop.initilise_partner_matrix(pop_parameters)
    partner_expire = pop.initilise_partner_duration(pop_parameters)
    file_name_pop = 'simulations/populations/scenario_' + str(scenario) + '/population_' + str(pop_no) + '.ftr'
    meta = pd.read_feather(file_name_pop)
    
    
    # Initilise import/export data
    pop_parameters = demo.initilise_demographic_dynamics(pop_parameters, inf_parameters, meta, partner_matrix, 0)
    partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.initilise_partnership_trackers(n_days)
    
    
    # Run Partnership Dynamics
    t0 = time.time()
    for t in range(0, n_days):
        
        
        # Update partnership dynamics
        meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, inf_parameters, meta, partner_matrix, partner_expire)
        meta, partner_matrix, partner_expire = prt.update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire)
        
        
        # Update partnership trackers
        partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.update_partnership_trackers(t, pop_parameters, meta, partner_matrix, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot)
    
    
    # Construct output
    t1 = time.time()
    out = partners_cum_tot[n_days-1,:]
    out = np.hstack((out, [t1-t0]))
    return_dict[sim_no] = out
    
    
    # Done
    return return_dict


#%% RUN WORKERS IN PARALLEL


if __name__ == '__main__':
    
    
    # Setup some parallel stuff
    manager = mp.Manager()
    return_dict = manager.dict()
    jobs = []
    
    
    # Iterate over the calibration sets
    for calib_set in range(0, sim_parameters['n_prt_param_sets'][0]):
        print('Running calibration set: ' + str(calib_set))
        
        
        # Test to see if output is there
        output_name = 'simulations/partnerships_calibration/scenario_' + str(scenario) + '/partners_cum_' + str(calib_set) + '.csv'
        if os.path.isfile(output_name) == False:
            
            
            # Run the thing
            for sim_no in range(n_cores):
                p = mp.Process(target = worker, args = (sim_no, calib_set, return_dict))
                jobs.append(p)
                p.start()
    
    
            # Not sure
            for proc in jobs:
                proc.join()
        
    
            # Save output
            out = pd.DataFrame.from_dict(data = return_dict, orient = 'index')
            out.to_csv(output_name)
            # print(out)

