# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:22:21 2021

@author: Nicolas Rebuli
"""

#%% SETUP Modules


# Standard modules
import numpy as np
import pandas as pd
import tqdm as tqdm
import random
import pickle
from pathlib import Path
import shutil
import os
import time
import matplotlib.pyplot as plt


# My modules
import src.demographic.generate_population as pop
import src.demographic.population_dynamics as demo
import src.partners.partners as prt
import src.calibration.setup as setup
import src.infections.ng as ng


# Set parameters
run_mode = 'serial'




#%% FUN run_one_simulation()
#
#
# Function which does a simulation run
#
#
def run_one_simulation(scenario = 1, parameter_no = 0, run_mode = run_mode):
    ( print('BATCHING: scenario: ' + str(scenario) + ' - parameter set: ' + str(parameter_no)) if run_mode != 'serial' else [] )
    ( print('Setting up simulation...') if run_mode == 'serial' else [] )
    

    # Time code
    t0 = time.time()


    # Setup simulation data
    sim_parameters, pop_parameters, prt_parameters, inf_parameters, meta, partner_expire, partner_matrix, population_no = \
        setup_data(scenario, parameter_no, run_mode, 'calibration')
    
    
    # Check to see if this dataset has been run to completion
    out_dir = 'simulations/calibration/scenario_' + str(scenario) +'/simulation_' + str(parameter_no)
    last_file = out_dir + '_prevalence.png'
    if os.path.exists(last_file) == False:
        
        
        # Initilise arrays for tracking transmission
        t0_sim = sim_parameters.partner_burn_in[0]
        tt = range(t0_sim, t0_sim + sim_parameters.simulation_length[0])
        yt = np.zeros((len(tt), 5))
        # popt = np.zeros((len(tt), 5))
        # part = np.zeros((len(tt), 2))
        inft = np.zeros((len(tt), 8))
        prvt = np.zeros((len(tt), 13))
        tstt = np.zeros((len(tt), 9))
        
        
        #% Run the simulation
        ( print('\nRunning simulation...\n') if run_mode == 'serial' else [] )
        for t in tqdm.tqdm(tt):
        # for t in tt:


            # Update population
            meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, prt_parameters, inf_parameters, meta, partner_matrix, partner_expire)


            # Update partnerships
            meta, partner_matrix, partner_expire = prt.update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire)


            # Update infections
            meta = ng.update_infections(pop_parameters, inf_parameters, meta, partner_matrix, t)


            # Update infection trackers
            yt, inft, prvt, tstt = update_tracking_data(meta, t - t0_sim, yt, inft, prvt, tstt)
            
    
        # Progress
        ( print('Tidying up...') if run_mode == 'serial' else [] )
        
        
        # Make ouput graphs
        ( print('Making summary graphs') if run_mode == 'serial' else [] )
        make_tracking_graphs(tt, parameter_no, yt, inft, prvt, tstt, out_dir)
        
        
        # Saving all environment variables
        ( print('Saving all environment variables') if run_mode == 'serial' else [] )
        del pop_parameters['lookup'], pop_parameters['imports'], pop_parameters['leavers']
        output = {'inf_parameters': inf_parameters,
                  'last_file': last_file,
                  # 'meta': meta,
                  'out_dir': out_dir,
                  'parameter_no': parameter_no,
                  # 'partner_expire': partner_expire,
                  # 'partner_matrix': partner_matrix,
                  'pop_parameters': pop_parameters,
                  'population_no': population_no,
                  'scenario': scenario,
                  'sim_parameters': sim_parameters,
                  't': t,
                  'yt': yt,
                  # 'popt': popt,
                  'inft': inft,
                  # 'prvt': prvt
                  'tstt': tstt
                  }
        out_file = open(out_dir + '_output_environment.pkl', 'wb')
        pickle.dump(output, out_file)
        out_file.close()
        
        
        # Saving some things in a more convenient format
        np.save(out_dir + '_output_prevalence.npy', prvt)
        np.save(out_dir + '_output_partner_matrix.npy', partner_matrix)
        np.save(out_dir + '_output_partner_expire.npy', partner_expire)
        meta.to_feather(out_dir + '_output_meta.ftr')
        


        # Print update
        runtime = (time.time() - t0)/60
        if run_mode != 'serial':
            print('COMPLETE: scenario: ' + str(scenario) + ' - parameter set: ' + str(parameter_no) + ' - population: ' + str(population_no) + ' - runtime: ' + str(runtime) + ' min')
        else:
            print('Run completed in ' + str(runtime) + ' min\n')
    
    
    # Done!
    return None



#%% FUN setup_data()
#
#
# Gets all the data ready in order to run the simulation script
#
#
def setup_data(scenario = 1, parameter_no = 0, run_mode = 'serial', inf_param_set = 'default'):
    
    
    # Read in global simulation parameters
    sim_parameters = pd.read_csv('data/param.csv')
    sim_parameters.loc[:, 'simulation_length'] = sim_parameters.loc[0, 'simulation_length_' + str(scenario)]
    
    
    # Read in demographic parameters
    pop_parameters = pop.setup_data(scenario = scenario,
                                    run_mode = run_mode)
    
    
    # Setup partnership parameters
    prt_parameters = prt.setup_data(pop_parameters,
                                    run_mode = run_mode)
    
    
    # Setup transmission parameters
    inf_parameters = setup.setup_transmission_parameters(set = inf_param_set,
                                                         scenario = scenario, 
                                                         parameter_no = parameter_no,
                                                         run_mode = run_mode)
    
    
    # Setup simulation population
    population_no = random.randint(0, sim_parameters.n_populations[0]-1)
    meta, partner_expire, partner_matrix = setup.parse_population_data(scenario = scenario, 
                                                                       population_no = population_no,
                                                                       run_mode = run_mode)
    
    
    # Preallocate importations
    pop_parameters = demo.initilise_demographic_dynamics(pop_parameters, 
                                                         inf_parameters, 
                                                         meta, 
                                                         partner_matrix, 
                                                         sim_parameters.partner_burn_in[0])
    
        
    # Done
    return sim_parameters, pop_parameters, prt_parameters, inf_parameters, meta, partner_expire, partner_matrix, population_no




#%% FUN setup_transmission_parameters()
#
#
# A function to load parameters.
#
# Set:
#      default = the baseline set of parameters used during initial testing
#
#      calibrated =
#
#
def setup_transmission_parameters(set = 'default', scenario = 3, parameter_no = 0, run_mode = run_mode):


    # This will just read in the default parameters stored in '~/data/'
    if set == 'default':


        # Read in the default set
        ( print('Infection parameters: default') if run_mode == 'serial' else [] )
        parameters = parse_default_transmission_parameters()

    
    # This will read in a specified set of calibration parameters
    elif set == 'calibration':


        # Read in the calibration parameters
        calib_param = pd.read_csv('simulations/parameters.csv')
        calib_param = calib_param.iloc[parameter_no,:]


        # Format the parameters
        ( print('Infection parameters: calibration set ' + str(parameter_no)) if run_mode == 'serial' else [] )
        parameters = parse_calibration_transmission_parameters(calib_param, parameter_no, set)

    
    # This will read in a specific set of calibration parameters
    elif set == 'calibrated':
        
        
        # Read in the calibrated parameters
        # calib_param = pd.read_csv('simulations/calibrated_scenario_' + str(scenario) + '.csv')
        calib_param = pd.read_csv('simulations/calibrated_scenario_3.csv')
        parameter_no = np.random.choice(len(calib_param), 1)[0]
        calib_param = calib_param.iloc[parameter_no,:]


        # Format the parameters
        ( print('Infection parameters: calibrated set ' + str(parameter_no)) if run_mode == 'serial' else [] )
        parameters = parse_calibration_transmission_parameters(calib_param, parameter_no, set)


    return parameters




#%% FUN parse_default_transmission_parameters()
#
#
# A function to parse the standard model parameters stored in the data folder.
#
#
def parse_default_transmission_parameters():


    # Parse parameters controlling infection dynamics and regular clinical treatment
    infection = pd.read_csv('data/ng_anatomical_site_infectious_period.csv')


    # Parse parameters controlling the probability of engaging in different
    # sexual acts. Read matrices as the probability of somebody of sex i
    # engaging in the given sex act with somebody of sex j. Where the order
    # of i and j are female then male.
    p_anal = pd.read_csv('data/sexual_act_probability_anal.csv').to_numpy()
    p_oral = pd.read_csv('data/sexual_act_probability_oral.csv').to_numpy()
    p_kiss = pd.read_csv('data/sexual_act_probability_kiss.csv').to_numpy()
    p_rim = pd.read_csv('data/sexual_act_probability_rim.csv').to_numpy()
    p_sex = pd.read_csv('data/sexual_act_probability_sex.csv').to_numpy()
    # p_event = pd.read_csv('data/sexual_act_rates.csv').to_numpy()[0]


    # Parse parameters controlling the probability of anatomical site specific
    # transmission as a result of engaging in particular sex acts


    # Form the site-to-site transmission probability matrix for convenience.
    # Read as the probability of transmission from site i to site j
    # Sites are in the order: rectal, urethral, pharangeal
    s2s = pd.read_csv('data/ng_anatomical_site_to_site_transmission_probabilities.csv').to_numpy()


    # Compose matrix with transmission probabilities as a result of anal sex
    trans_anal = np.array([[0,          s2s[0,1],   0],
                           [s2s[1,0],   0,          0],
                           [0,          0,          0]])


    # Compose matrix with transmission probabilities as a result of oral sex
    trans_oral = np.array([[0,          0,          0],
                           [0,          0,          s2s[1,2]],
                           [0,          s2s[2,1],   0]])


    # Compose matrix with transmission probabilities as a result of kissing
    trans_kiss = np.array([[0,          0,          0],
                           [0,          0,          0],
                           [0,          0,          s2s[2,2]]])


    # Compose matrix with transmission probabilities as a result of rimming
    trans_rim = np.array([[0,           0,          s2s[0,2]],
                          [0,           0,          0],
                          [s2s[2,0],    0,          0]])


    # Compose matrix with transmission probabilities as a result of sex
    trans_sex = np.array([[0,           0,          0],
                          [0,           s2s[1,1],   0],
                          [0,           0,          0]])


    # Combine everything into a dictionary
    parameters = {'set': 'default',
                  'infection': infection,
                  'p_anal': p_anal,
                  'p_oral': p_oral,
                  'p_kiss': p_kiss,
                  'p_sex': p_sex,
                  'p_rim': p_rim,
                  # 'p_event': p_event,
                  'trans_anal': trans_anal,
                  'trans_oral': trans_oral,
                  'trans_kiss': trans_kiss,
                  'trans_sex': trans_sex,
                  'trans_rim': trans_rim}


    # Return the parameters in dictionary form
    return parameters




#%% FUN parse_calibration_transmission_parameters()
#
#
# A function to parse the model parameters stored in a csv to the
# format that the code requires.
#
#
def parse_calibration_transmission_parameters(calib_param, parameter_no, set):
    
    
    # Read in the default transmission parameters
    parameters = parse_default_transmission_parameters()
    
    
    # Update the parameter set
    parameters['set'] = [set + ' set ' + str(parameter_no)]
    
    
    # Form the site-to-site transmission probability matrix for convenience.
    # Read as the probability of transmission from site i to site j
    # Sites are in the order: rectal, urethral, pharangeal
    s2s = np.array([[0,               calib_param.pru, calib_param.prp],
                    [calib_param.pur, calib_param.puu, calib_param.pup],
                    [calib_param.ppr, calib_param.ppu, calib_param.ppp]])


    # Compose matrix with transmission probabilities as a result of anal sex
    trans_anal = np.array([[0,          s2s[0,1],   0],
                           [s2s[1,0],   0,          0],
                           [0,          0,          0]])


    # Compose matrix with transmission probabilities as a result of oral sex
    trans_oral = np.array([[0,          0,          0],
                           [0,          0,          s2s[1,2]],
                           [0,          s2s[2,1],   0]])


    # Compose matrix with transmission probabilities as a result of kissing
    trans_kiss = np.array([[0,          0,          0],
                           [0,          0,          0],
                           [0,          0,          s2s[2,2]]])


    # Compose matrix with transmission probabilities as a result of rimming
    trans_rim = np.array([[0,           0,          s2s[0,2]],
                          [0,           0,          0],
                          [s2s[2,0],    0,          0]])


    # Compose matrix with transmission probabilities as a result of sex
    trans_sex = np.array([[0,           0,          0],
                          [0,           s2s[1,1],   0],
                          [0,           0,          0]])
    
    
    # Update all the calibration parameters
    parameters['trans_anal'] = trans_anal
    parameters['trans_oral'] = trans_oral
    parameters['trans_kiss'] = trans_kiss
    parameters['trans_sex'] = trans_sex
    parameters['trans_rim'] = trans_rim
    parameters['infection'].loc[:, 'symptoms_rectal'] = calib_param['symptoms_rectal']
    parameters['infection'].loc[0, 'symptoms_urethral'] = calib_param['symptoms_ural_male']
    parameters['infection'].loc[1, 'symptoms_urethral'] = calib_param['symptoms_ural_female']


    # Return the parameters in dictionary form
    return parameters




#%% FUN parse_population_data()
#
#
# A function which loads in a set of population data
#
#
def parse_population_data(scenario = 1, population_no = 0, run_mode = run_mode):


    # Read in a specific set
    ( print('Population set: ' + str(population_no), flush = True) if run_mode == 'serial' else [] )


    # Read in the population metadata dataframe
    meta = pd.read_feather('simulations/partnerships/scenario_' + str(scenario) + '/population_' + str(population_no) + '_meta.ftr')


    # Check that the duration exposed is correct in meta
    sim_parameters = pd.read_csv('data/param.csv')
    meta.loc[meta.site0 == 1, 'site0_t0'] = sim_parameters.partner_burn_in[0] + sim_parameters.init_duration_exposed[0] * np.random.random(sum(meta.site0 == 1))
    meta.loc[meta.site1 == 1, 'site1_t0'] = sim_parameters.partner_burn_in[0] + sim_parameters.init_duration_exposed[0] * np.random.random(sum(meta.site1 == 1))
    meta.loc[meta.site2 == 1, 'site2_t0'] = sim_parameters.partner_burn_in[0] + sim_parameters.init_duration_exposed[0] * np.random.random(sum(meta.site2 == 1))
    
    
    # Check for some new variables
    if 'last_test_time' not in meta.columns:
        meta.loc[:, 'test_reason_last'] = 0
        meta.loc[:, 'test_time_last'] = -float('Inf')
        meta.loc[:, 'test_time'] = float('Inf')
        meta.loc[:, 'treatment_time'] = float('Inf')
        meta.loc[:, 'vaccinated'] = False
        meta.loc[:, 'vaccination_t0'] = float("inf")
        meta.loc[:, 'vaccination_t1'] = float("inf")
        meta.loc[:, 'booster_t0'] = float("inf")
    

    # Read in the simulated partnership data arrays
    partner_expire = np.load('simulations/partnerships/scenario_' + str(scenario) + '/population_' + str(population_no) + '_expire.npy')
    partner_matrix = np.load('simulations/partnerships/scenario_' + str(scenario) + '/population_' + str(population_no) + '_matrix.npy')


    # Return data
    return meta, partner_expire, partner_matrix




#%% FUN update_tracking_data()
#
#
# A function which updates the arrays used to track transmission
#
#
def update_tracking_data(meta, t, yt, inft, prvt, tstt):
    
    
    # Compute the current population size
    n = len(meta)
    
    
    # Compute proportion of population in each category
    yt[t,:] = [sum(meta.state == 'S')/n,
               sum(meta.state == 'E')/n,
               sum(meta.state == 'I')/n,
               sum(meta.state == 'R')/n,
               sum(meta.state == 'T')/n]
    
    
    # # Compute number of individuals in each age group
    # popt[t,:] = [sum(meta.age_group == 0),
    #             sum(meta.age_group == 1),
    #             sum(meta.age_group == 2),
    #             sum(meta.age_group == 3),
    #             n]
    
    
    # # Compute the total number of long-term relationships
    # part[t,:] = [sum(meta.partner == -1),
    #             sum(meta.partner > -1)]
    
    
    # Compute the proportion of the population with an infection by anatomical site
    site0 = meta.site0 == 1
    site1 = meta.site1 == 1
    site2 = meta.site2 == 1
    inft[t,:] = [sum( ~site0 & ~site1 & ~site2 )/n,
                sum( site0 & ~site1 & ~site2 )/n,
                sum( ~site0 & site1 & ~site2 )/n,
                sum( ~site0 & ~site1 & site2 )/n,
                sum( site0 & site1 & ~site2 )/n,
                sum( ~site0 & site1 & site2 )/n,
                sum( site0 & ~site1 & site2 )/n,
                sum( site0 & site1 & site2 )/n]
    
    
    # Compute comparisons to STRIVE data - prevalence of urethral infections by sub-populations
    I = meta.state == 'I'
    males = meta.gender == 1
    prvt[t,:] = [sum(I)/n,
                 sum(I & (males))/sum(males),
                 sum(I & (~males))/sum(~males),
                 sum(I & (males) & (meta.age_group == 0))/max(1, sum( (meta.gender == 1) & (meta.age_group == 0) ) ),
                 sum(I & (males) & (meta.age_group == 1))/max(1, sum( (meta.gender == 1) & (meta.age_group == 1) ) ),
                 sum(I & (males) & (meta.age_group == 2))/max(1, sum( (meta.gender == 1) & (meta.age_group == 2) ) ),
                 sum(I & (males) & (meta.age_group == 3) & (meta.age < 35))/max(1, sum( (males) & (meta.age_group == 3) & (meta.age < 35) )),
                 sum(I & (males) & (meta.age_group == 3) & (meta.age > 35))/max(1, sum( (males) & (meta.age_group == 3) & (meta.age > 35) )),
                 sum(I & (~males) & (meta.age_group == 0))/max(1, sum( (~males) & (meta.age_group == 0) ) ),
                 sum(I & (~males) & (meta.age_group == 1))/max(1, sum( (~males) & (meta.age_group == 1) ) ),
                 sum(I & (~males) & (meta.age_group == 2))/max(1, sum( (~males) & (meta.age_group == 2) ) ),
                 sum(I & (~males) & (meta.age_group == 3) & (meta.age < 35))/max(1, sum( (~males) & (meta.age_group == 3) & (meta.age < 35) )),
                 sum(I & (~males) & (meta.age_group == 3) & (meta.age > 35))/max(1, sum( (~males) & (meta.age_group == 3) & (meta.age > 35) ))]
    
    
    # Compute the proportion of the population who have tested in the last 365 days
    tested = meta.test_time_last >= (t-365)
    tstt[t,:] = [np.sum(tested) / n,
                 np.sum(males & tested & (meta.test_reason_last == int(1))) / np.sum(males),
                 np.sum(males & tested & (meta.test_reason_last == int(2))) / np.sum(males),
                 np.sum(males & tested & (meta.test_reason_last == int(3))) / np.sum(males),
                 np.sum(males & tested) / np.sum(males),
                 np.sum(~males & tested & (meta.test_reason_last == int(1))) / np.sum(~males),
                 np.sum(~males & tested & (meta.test_reason_last == int(2))) / np.sum(~males),
                 np.sum(~males & tested & (meta.test_reason_last == int(3))) / np.sum(~males),
                 np.sum(~males & tested) / np.sum(~males),
                 ]
    
    
    # Done
    return yt, inft, prvt, tstt




#%% make_tracking_graphs()
#
#
# Makes graphs
#
#
def make_tracking_graphs(tt, sim, yt, inft, prvt, tstt, out_dir):


    # Plot aggregate infection levels
    new_fig = plt.figure(figsize = [12, 8])
    yt = 100 * yt
    plt.plot(tt, yt[:,0], label = 'S')
    plt.plot(tt, yt[:,1], label = 'E')
    plt.plot(tt, yt[:,2], label = 'I')
    plt.plot(tt, yt[:,3], label = 'R')
    plt.plot(tt, yt[:,4], label = 'T')
    plt.title('Proportion of the Population in Each Model Compartment (Parameter Set ' + str(sim) + ')')
    plt.legend()
    plt.xlabel('Global Simulation Timestep')
    plt.ylabel('Proportion of Population (%)')
    plt.savefig(out_dir + '_graph_aggregate_infections.png')
    plt.close(new_fig)


    # # Plot the number in each age group
    # new_fig = plt.figure(figsize = [12, 8])
    # plt.plot(tt, popt[:,0], label = '16-19')
    # plt.plot(tt, popt[:,1], label = '20-24')
    # plt.plot(tt, popt[:,2], label = '25-30')
    # plt.plot(tt, popt[:,3], label = 'Over 30')
    # plt.plot(tt, np.sum(popt, axis = 1), label = 'Total')
    # plt.title('Number of People in Each Age Group (Parameter Set ' + str(sim) +')')
    # plt.xlabel('Global Simulation Timestep')
    # plt.ylabel('Number of Individuals')
    # plt.legend()
    # plt.savefig(out_dir + '_graph_age_group.png')
    # plt.close()


    # # Plot the number in a long term relationship
    # plt.plot(tt, part[:,0], label = 'Single')
    # plt.plot(tt, popt[:,1], label = 'Long-term relationship')
    # plt.title('Long-term Relationships - Parameter Set ' + str(sim))
    # plt.legend()
    # plt.savefig(file_name + '_long_term_relationships.png')
    # plt.close()


    # Graph of infections by anatomical site
    # Setup graph
    fig, axs = plt.subplots(3, 1, figsize = [12, 8])
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Proportion of Population with Infection by Anatomical-Site (Parameter Set ' + str(sim) + ')')
    
    
    # make first pannel
    #inft = 100 * inft
    axs[0].plot(tt, inft[:,1] + inft[:,4] + inft[:,6] + inft[:,7], label = 'Rectal')
    axs[0].plot(tt, inft[:,2] + inft[:,4] + inft[:,5] + inft[:,7], label = 'Urethral')
    axs[0].plot(tt, inft[:,3] + inft[:,5] + inft[:,6] + inft[:,7], label = 'Pharyngeal')
    axs[0].legend()
    axs[0].set_title('Infections Including Each Site')
    
    
    # make second pannel
    axs[1].plot(tt, inft[:,1], label = 'Rectal')
    axs[1].plot(tt, inft[:,2], label = 'Urethral')
    axs[1].plot(tt, inft[:,3], label = 'Pharyngeal')
    axs[1].legend()
    axs[1].set_title('Infections at Just One Site')
    
    
    # Make third pannel
    axs[2].plot(tt, inft[:,4], label = 'Sites Rec and Ure')
    axs[2].plot(tt, inft[:,5], label = 'Sites Ure and Pha')
    axs[2].plot(tt, inft[:,6], label = 'Sites Rec and Pha')
    axs[2].plot(tt, inft[:,7], label = 'All sites')
    axs[2].set_title('Infections at Multiple Sites')
    
    
    # Label output
    axs[2].legend()
    axs[2].set_xlabel('Global Simulation Timestep')
    axs[1].set_ylabel('Proportion of the Population (%)')
    
    
    # Save output
    plt.savefig(out_dir + '_graph_anatomical_site.png')
    plt.close()


    # PLOTS OF PREVALENCE


    # Overall prevalence
    prvt = prvt * 100
    fig, axs = plt.subplots(3, 1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Aggregate Prevalence of Urogenetal Infections Compared to Prevalence in STRIVE (Parameter Set ' + str(sim) + ')')
    
    
    # Overall
    axs[0].plot(tt, prvt[:,0], label = 'Simulated')
    axs[0].plot(tt, len(tt) * [9.5], color = 'black', linestyle = '--', label = 'STRIVE')
    axs[0].set_title('Overall Prevalence')
    axs[0].legend()
    
    
    # Males
    axs[1].plot(tt, prvt[:,1])
    axs[1].plot(tt, len(tt) * [10.4], color = 'black', linestyle = '--')
    axs[1].set_title('Prevalence Amongst Males')
    
    
    # Females
    axs[2].plot(tt, prvt[:,2])
    axs[2].plot(tt, len(tt) * [8.9], color = 'black', linestyle = '--')
    axs[2].set_title('Prevalence Amongst Females')
    
    
    # Labels
    axs[2].set_xlabel('Global Simulation Timestep')
    axs[1].set_ylabel('Proportion of the Population (%)')
    
    
    # Save output
    plt.savefig(out_dir + '_graph_prevalence_aggregate.png')
    plt.close()


    # Prevalence by age group and gender
    fig, axs = plt.subplots(4, 2)
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Simulated Prevalence Compared to Prevalence in STRIVE (Parameter Set ' + str(sim) + ')')
    
    
    # Males
    axs[0, 0].set_title('Males')
    axs[0, 0].plot(tt, prvt[:,3], label = 'Simulated')
    axs[0, 0].plot(tt, len(tt) * [21.6], color = 'black', linestyle = '--', label = 'STRIVE')
    axs[1, 0].plot(tt, prvt[:,4])
    axs[1, 0].plot(tt, len(tt) * [17.4], color = 'black', linestyle = '--')
    axs[2, 0].plot(tt, prvt[:,5])
    axs[2, 0].plot(tt, len(tt) * [11.6], color = 'black', linestyle = '--')
    axs[3, 0].plot(tt, prvt[:,6])
    axs[3, 0].plot(tt, len(tt) * [8.1], color = 'black', linestyle = '--')
    
    
    # Females
    axs[0, 1].set_title('Females')
    axs[0, 1].plot(tt, prvt[:,8])
    axs[0, 1].plot(tt, len(tt) * [20.1], color = 'black', linestyle = '--')
    axs[1, 1].plot(tt, prvt[:,9])
    axs[1, 1].plot(tt, len(tt) * [15.4], color = 'black', linestyle = '--')
    axs[2, 1].plot(tt, prvt[:,10])
    axs[2, 1].plot(tt, len(tt) * [7.3], color = 'black', linestyle = '--')
    axs[3, 1].plot(tt, prvt[:,11])
    axs[3, 1].plot(tt, len(tt) * [7], color = 'black', linestyle = '--')
    
    
    # Axis labels
    axs[0, 0].legend
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    fig.text(0.5, 0.05, 'Global Simulation Timestep', ha='center')
    fig.text(0.05, 0.5, 'Proportion of Sub-Populations (%)', va='center', rotation='vertical')
    
    
    # Save output
    plt.savefig(out_dir + '_graph_prevalence_age_group.png')
    plt.close()
    
    
    ## PLOT OF TESTING
    
    
    # Overall testing rates
    tstt = tstt * 100
    fig, axs = plt.subplots(3, 1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Overall Testing Rates Compared to Mean Testing Rates in STRIVE (Parameter Set ' + str(sim) + ')')
    
    
    # Overall
    axs[0].plot(tt, tstt[:,0], label = 'Simulated')
    axs[0].plot(tt, len(tt) * [17.6], color = 'black', linestyle = '--', label = 'STRIVE')
    axs[0].set_title('Overall Prevalence')
    axs[0].legend()
    
    
    # Males
    axs[1].plot(tt, tstt[:,4], label = 'Overall')
    axs[1].plot(tt, tstt[:,1], label = 'Symptoms')
    axs[1].plot(tt, tstt[:,2], label = 'Background')
    axs[1].plot(tt, tstt[:,3], label = 'Contact Tracing')
    axs[1].plot(tt, len(tt) * [13.4], color = 'black', linestyle = '--')
    axs[1].set_title('Testing Rates Amongst Males by Reason for Last Test')
    
    
    # Females
    axs[2].plot(tt, tstt[:,8], label = 'Overall')
    axs[2].plot(tt, tstt[:,5], label = 'Symptoms')
    axs[2].plot(tt, tstt[:,6], label = 'Background')
    axs[2].plot(tt, tstt[:,7], label = 'Contact Tracing')
    axs[2].plot(tt, len(tt) * [21.8], color = 'black', linestyle = '--')
    axs[2].set_title('Testing Rates Amongst Females by Reason for Last Test')
    
    
    # Labels
    axs[2].set_xlabel('Global Simulation Timestep')
    axs[1].set_ylabel('Proportion of Sub-Population (%)')
    axs[1].legend()
    
    
    # Save output
    plt.savefig(out_dir + '_graph_testing_rates.png')
    plt.close()


    # Done
    return None
    
    





































