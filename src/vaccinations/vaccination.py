# -*- coding: utf-8 -*-
"""
Created on Mon May 24 09:10:03 2021

@author: nicol
"""


#%% SETUP Modules


# Standard libaries
import numpy as np
import pandas as pd
import scipy.stats as sp
import tqdm as tqdm
import matplotlib.pyplot as plt
import time
import os
import pickle


# My modules
import src.demographic.generate_population as pop
import src.demographic.population_dynamics as demo
import src.partners.partners as prt
import src.infections.ng as ng
import src.calibration.setup as setup




#%% NG run_vaccine_scenario()
# FUNCTION TO RUN A VACINE SCENARIO
#
#
# Top-level function which actually runs the simulations.
#
#
def run_vaccine_scenario(vax_parameters, sim_no, run_mode = 'serial'):
    ( print('BATCHING: scenario: ' + vax_parameters['vax_scenario'] + ' - simulation number: ' + str(sim_no)) if run_mode != 'serial' else [] )
    ( print('Setting up simulation...') if run_mode == 'serial' else [] )
    

    # Time code
    t0 = time.time()


    # Setup transmission functions
    ( print('Vaccination scenario: ' + vax_parameters['vax_scenario']) if run_mode == 'serial' else [] )
    update_infections = set_function_for_updating_infections(vax_parameters)


    # Set top-level simulation parameters
    scenario = vax_parameters['pop_scenario']
    inf_param_set = 'calibrated'


    # Setup simulation data
    sim_parameters, pop_parameters, prt_parameters, inf_parameters, meta, partner_expire, partner_matrix, population_no = \
        setup.setup_data(scenario = scenario,
                        run_mode = run_mode,
                        inf_param_set = inf_param_set)


    # Setup simulation runtime
    n_steps = 365*vax_parameters['n_years'] if 'n_days' in vax_parameters else  sim_parameters.simulation_length[0]
    t0_sim = sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0]
    tt = range(t0_sim, t0_sim + n_steps)
    
    
    # Check to see if the output directory is there
    out_dir = 'simulations/vaccination/' + vax_parameters['vax_scenario']
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
    out_dir = out_dir + '/simulation_' + str(sim_no) + '_'
    
    
    # Initilise arrays for tracking transmission
    inf_tracker = initilise_tracking_data(len(tt))
    
    
    #% Run the simulation
    ( print('\nRunning simulation...\n') if run_mode == 'serial' else [] )
    # for t in tqdm.tqdm(tt):
    for t in tt:


        # Update population
        meta, partner_matrix, partner_expire = demo.update_population(t, pop_parameters, prt_parameters, inf_parameters, meta, partner_matrix, partner_expire)


        # Update partnerships
        meta, partner_matrix, partner_expire = prt.update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire)


        # Update infections
        meta = update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix)


        # Update infection trackers
        inf_tracker = update_tracking_data(meta, t0_sim, t - t0_sim, inf_tracker)
    
    
    # Progress
    ( print('Tidying up...') if run_mode == 'serial' else [] )
    
    
    # Make ouput graphs
    ( print('Making summary graphs') if run_mode == 'serial' else [] )
    make_tracking_graphs(tt, sim_no, inf_tracker, out_dir)
    
    
    # Saving all environment variables
    ( print('Saving all environment variables') if run_mode == 'serial' else [] )
    del pop_parameters['lookup'], pop_parameters['imports'], pop_parameters['leavers']
    output = {  'scenario': scenario,
                'population_no': population_no,
                'sim_no': sim_no,
                'sim_parameters': sim_parameters,
                'pop_parameters': pop_parameters,
                'prt_parameters': prt_parameters,
                'inf_parameters': inf_parameters,
                'vax_parameters': vax_parameters,
                'tt': tt,
                'inf_tracker': inf_tracker,
                }
    out_file = open(out_dir + 'output_environment.pkl', 'wb')
    pickle.dump(output, out_file)
    out_file.close()
    
    
    # Saving some things in a more convenient format
    np.save(out_dir + 'output_prevalence.npy', inf_tracker['prvt'])
    np.save(out_dir + 'output_partner_matrix.npy', partner_matrix)
    np.save(out_dir + 'output_partner_expire.npy', partner_expire)
    meta.to_feather(out_dir + 'output_meta.ftr')
    

    # Print update
    runtime = (time.time() - t0)
    h = int(np.floor(runtime/(60*60)))
    m = int(np.floor(runtime/60 - 60*h))
    s = int(np.floor(runtime - 60*60*h - 60*m))
    if run_mode != 'serial':
        print('COMPLETE: scenario: ' + vax_parameters['vax_scenario'] + ' - simulation number: ' + str(sim_no) + ' - runtime: ' + str(h) + 'h ' + str(m) + 'm ' + str(s) + 's')
    else:
        print('Run completed in ' + str(h) + 'h ' + str(m) + 'm ' + str(s) + 's\n')


    # Done!
    return None




#%% FUN set_function_for_updating_infections()
#
#
# Top-level function which handles how the vaccine parameters
# are interpreted by the simulation code
#
#
def set_function_for_updating_infections(vax_parameters):

    
    ##########################################
    ##  ALL EFFECTS WITH DEPLOYMENT MODE 0  ##
    ##########################################
    
    
    # Where distribution occurrs during treatment
    if vax_parameters['deployment'] == 0:


        # Where the effect is total immunity
        if vax_parameters['effect'] == 0:


            # Situation where vaccinations given during treatment and makes people immune
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_infection_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce transmission
        elif vax_parameters['effect'] == 1:


            # Situation where vaccines given during treatment and result is to
            # decrease transmission probability
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce symptoms
        elif vax_parameters['effect'] == 2:


            # Situation where vaccines given during treatment and result in
            # people being more likely to be asymptomatic
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         symptoms_prob = symptoms_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is the reduce the duration of infection
        elif vax_parameters['effect'] == 3:


            # Situation where vaccines given during treatment and result in
            # people having a shorter duration of infection
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta
            
            
        # Where the effect includes all of the above
        elif vax_parameters['effect'] == 4:
            
            
            # Situation where vaccines have a combination of all of the
            # above effects
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_and_inf_prob, 
                                         symptoms_prob = symptoms_vax,
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta
            
            
    ##########################################
    ##  ALL EFFECTS WITH DEPLOYMENT MODE 1  ##
    ##########################################
    
    
    # Where distribution occurs at the age of 16
    elif vax_parameters['deployment'] == 1:


        # Where the effect is total immunity
        if vax_parameters['effect'] == 0:


            # Situation where vaccinations given during treatment and makes people immune
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_infection_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce transmission
        elif vax_parameters['effect'] == 1:


            # Situation where vaccines given during treatment and result is to
            # decrease transmission probability
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce symptoms
        elif vax_parameters['effect'] == 2:


            # Situation where vaccines given during treatment and result in
            # people being more likely to be asymptomatic
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         symptoms_prob = symptoms_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is the reduce the duration of infection
        elif vax_parameters['effect'] == 3:


            # Situation where vaccines given during treatment and result in
            # people having a shorter duration of infection
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta
            
            
        # Where the effect includes all of the above
        elif vax_parameters['effect'] == 4:
            
            
            # Situation where vaccines have a combination of all of the
            # above effects
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_and_inf_prob, 
                                         symptoms_prob = symptoms_vax,
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta
            
            
    ##########################################
    ##  ALL EFFECTS WITH DEPLOYMENT MODE 2  ##
    ##########################################
    
    
    # Where distribution occurs at the age of 16
    elif vax_parameters['deployment'] == 2:


        # Where the effect is total immunity
        if vax_parameters['effect'] == 0:


            # Situation where vaccinations given during treatment and makes people immune
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_infection_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce transmission
        elif vax_parameters['effect'] == 1:


            # Situation where vaccines given during treatment and result is to
            # decrease transmission probability
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce symptoms
        elif vax_parameters['effect'] == 2:


            # Situation where vaccines given during treatment and result in
            # people being more likely to be asymptomatic
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         symptoms_prob = symptoms_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is the reduce the duration of infection
        elif vax_parameters['effect'] == 3:


            # Situation where vaccines given during treatment and result in
            # people having a shorter duration of infection
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta
            
            
        # Where the effect includes all of the above
        elif vax_parameters['effect'] == 4:
            
            
            # Situation where vaccines have a combination of all of the
            # above effects
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_and_inf_prob, 
                                         symptoms_prob = symptoms_vax,
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta
            
            
    ##########################################
    ##  ALL EFFECTS WITH DEPLOYMENT MODE 3  ##
    ##########################################
    
    
    # Where distribution occurs at the age of 16
    elif vax_parameters['deployment'] == 3:


        # Where the effect is total immunity
        if vax_parameters['effect'] == 0:


            # Situation where vaccinations given during treatment and makes people immune
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_infection_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce transmission
        elif vax_parameters['effect'] == 1:


            # Situation where vaccines given during treatment and result is to
            # decrease transmission probability
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_prob, 
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce symptoms
        elif vax_parameters['effect'] == 2:


            # Situation where vaccines given during treatment and result in
            # people being more likely to be asymptomatic
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         symptoms_prob = symptoms_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is the reduce the duration of infection
        elif vax_parameters['effect'] == 3:


            # Situation where vaccines given during treatment and result in
            # people having a shorter duration of infection
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta
            
            
        # Where the effect includes all of the above
        elif vax_parameters['effect'] == 4:
            
            
            # Situation where vaccines have a combination of all of the
            # above effects
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = vax_at_debut(vax_parameters, meta, t)
                meta = vax_by_targets(vax_parameters, meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t, 
                                         trans_prob_fun = vaccine_reduced_trans_and_inf_prob, 
                                         symptoms_prob = symptoms_vax,
                                         duration_mod = duration_vax,
                                         vax_parameters = vax_parameters)
                meta = ng.progress_state_of_infection(meta, t)
                meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta

    
    # Done
    return update_infections




#%% FUN progress_state_of_vaccination()
#
#
#
#
#
def progress_state_of_vaccination(meta, t):


    # Check for new people in the population
    new = meta.vaccinated.isin([True, False]) == False
    meta.loc[new, 'vaccinated'] = False
    meta.loc[new, 'boosted'] = False
    meta.loc[new, 'vaccination_t0'] = float('inf')
    meta.loc[new, 'vaccination_t1'] = float('inf')
    meta.loc[new, 'booster_t0'] = float('inf')
    meta.loc[new, 'booster_t1'] = float('inf')


    # Remove vaccine effect
    waned = meta.loc[:, 'vaccination_t1'] < t
    meta.loc[waned, 'vaccinated'] = False
    meta.loc[waned, 'vaccination_t1'] = float('inf')
    
    
    # Check for boosters
    boost = (meta.loc[:, 'booster_t0'] < t) & (meta.loc[:, 'boosted'] == False)
    meta.loc[boost, 'vaccinated'] = True
    meta.loc[boost, 'boosted'] = True
    
    
    # Remove booster effect
    waned = (meta.loc[:, 'booster_t1'] < t) & (meta.loc[:, 'vaccinated'])
    meta.loc[waned, 'vaccinated'] = False  
    
    
    # Done
    return meta




#%% FUN seek_treatment()
# FUNCTION FOR GIVING VACCINE AT TIME OF TREATMENT
#
#
# Just calles ng.seek_treatment() but then does give_vaccine() on those
# who get treated.
# 
#
def seek_treatment(pop_parameters, parameters, vax_parameters, meta, partner_matrix, t):

    
    # Run the seek treatment function from ng but ask for the index values of those
    # who are treated in response
    meta, treat = ng.seek_treatment(pop_parameters, parameters, meta, partner_matrix, t, return_treated_ids=True)
    
    
    # People getting treated who haven't been vaccinated yet get a vaccine
    vax = treat[np.isinf(meta.vaccination_t0[treat])]
    

    # Implement vaccinations
    if len(vax) > 0:
        meta = give_vaccine(vax_parameters, meta, t, vax)


    # Done
    return meta




#%% FUN give_vaccine()
#
#
# What happens when somebody is vaccinated 
#
#
def give_vaccine(vax_parameters, meta, t, treat):


    # Decide which vaccinations are effective
    uu = np.random.random(len(treat))
    uu = uu < vax_parameters['prop_effective']


    # Update their vaccination data
    meta.loc[treat[uu], 'vaccinated'] = True
    meta.loc[treat[uu], 'vaccination_t0'] = t
    meta.loc[treat[uu], 'vaccination_t1'] = t + vaccination_duration(vax_parameters, len(treat[uu]))


    # Set booster time
    meta.loc[treat[uu], 'booster_t0'] = t + vax_parameters['booster_delay']
    meta.loc[treat[uu], 'booster_t1'] = t + vax_parameters['booster_delay'] + vaccination_duration(vax_parameters, len(treat[uu]))


    # Pass back meta
    return meta




#%% FUN vaccination_duration()
#
#
# Duration of time the vaccine is effective
# Drawn from Gamma distribution.
#
#
def vaccination_duration(vax_parameters, n):


    # Sample duration
    duration = np.random.gamma(vax_parameters['duration_mean'] / \
                               vax_parameters['duration_var'], \
                               vax_parameters['duration_var'], n)

    
    # Done
    return duration




#%% FUN vaccine_reduced_infection_prob()
#
#
# Compute the site-specific transmission probabilities
#
#
def vaccine_reduced_infection_prob(inf_parameters, vax_parameters, meta, i, j):


    # Compute the usual transmission probability
    trans_prob = ng.transmission_probability(inf_parameters, vax_parameters, meta, i, j)

    
    # Check if infectee is vaccinated
    vaccinated = 1.0 * meta.vaccinated.iloc[j]
    
    
    # Multiply transmission probability
    vax_reduction = np.array([[vax_parameters['site0_inf_mult'] * vaccinated + (1-vaccinated)], \
                              [vax_parameters['site1_inf_mult'] * vaccinated + (1-vaccinated)], \
                              [vax_parameters['site2_inf_mult'] * vaccinated + (1-vaccinated)]])
                
    
    # Multiply the transmission probability by the vaccine reduction
    out = vax_reduction * trans_prob

    
    # Done
    return out




#%% FUN vaccine_reduced_trans_prob()
#
#
# Compute the site-specific transmission probabilities
#
#
def vaccine_reduced_trans_prob(inf_parameters, vax_parameters, meta, i, j):


    # Compute the usual transmission probability
    trans_prob = ng.transmission_probability(inf_parameters, vax_parameters, meta, i, j)
    
    
    # Identify if the infector is vaccinated
    vaccinated = 1.0 * meta.vaccinated.iloc[i]
    
    
    # Construct the reduction vector
    vax_reduction = np.array([[vax_parameters['site0_trans_mult'] * vaccinated + (1-vaccinated)], \
                              [vax_parameters['site1_trans_mult'] * vaccinated + (1-vaccinated)], \
                              [vax_parameters['site2_trans_mult'] * vaccinated + (1-vaccinated)]])
                
    
    # Multiply the transmission probability by the vaccine reduction
    out = vax_reduction * trans_prob

    
    # Done
    return out




#%% FUN vaccine_reduced_trans_and_inf_prob()
#
#
# Compute the site-specific transmission probabilities
#
#
def vaccine_reduced_trans_and_inf_prob(inf_parameters, vax_parameters, meta, i, j):


    # Compute the usual transmission probability
    trans_prob = ng.transmission_probability(inf_parameters, vax_parameters, meta, i, j)
    
    
    # Transmission reduction
    infector = 1.0 * meta.vaccinated.iloc[i]
    vax_reduction = np.array([[vax_parameters['site0_trans_mult'] * infector + (1-infector)], \
                              [vax_parameters['site1_trans_mult'] * infector + (1-infector)], \
                              [vax_parameters['site2_trans_mult'] * infector + (1-infector)]])
        
        
    # Infection reduction
    infectee = 1.0 * meta.vaccinated.iloc[j]
    vax_reduction[0] = vax_reduction[0] * (vax_parameters['site0_inf_mult'] * infectee + (1-infectee))
    vax_reduction[1] = vax_reduction[1] * (vax_parameters['site1_inf_mult'] * infectee + (1-infectee))
    vax_reduction[2] = vax_reduction[2] * (vax_parameters['site2_inf_mult'] * infectee + (1-infectee))
                
    
    # Multiply the transmission probability by the vaccine reduction
    out = vax_reduction * trans_prob

    
    # Done
    return out



#%% FUN symptoms_rectal
#
#
# Reduces the probability of symptoms
#
#
def symptoms_rectal(inf_parameters, vax_parameters, meta, i, j):

    
    # Compute the baseline porbability
    p_baseline = inf_parameters['infection'].symptoms_rectal[meta.gender.iloc[j]]


    # Modify the probability using the vaxination parameters
    modifier = meta.vaccinated.iloc[j] * vax_parameters['site0_symp_mult'] + (1-meta.vaccinated.iloc[j])
    prob = modifier * p_baseline


    # Decide whether or not they will be symptomatic
    symptoms = np.random.random() < prob

    
    # Done
    return symptoms




#%% FUN symptoms_pharynx
#
#
# Reduces the probability of symptoms
#
#
def symptoms_pharynx(inf_parameters, vax_parameters, meta, i, j):

    
    # Compute the baseline porbability
    p_baseline = inf_parameters['infection'].symptoms_pharyngeal[meta.gender.iloc[j]]


    # Modify the probability using the vaxination parameters
    modifier = meta.vaccinated.iloc[j] * vax_parameters['site1_symp_mult'] + (1-meta.vaccinated.iloc[j])
    prob = modifier * p_baseline


    # Decide whether or not they will be symptomatic
    symptoms = np.random.random() < prob


    # Done
    return symptoms




#%% FUN symptoms_urethra
#
#
# Reduces the probability of symptoms
#
#
def symptoms_urethra(inf_parameters, vax_parameters, meta, i, j):

    
    # Compute the baseline porbability
    p_baseline = inf_parameters['infection'].symptoms_urethral[meta.gender.iloc[j]]


    # Modify the probability using the vaxination parameters
    modifier = meta.vaccinated.iloc[j] * vax_parameters['site2_symp_mult'] + (1-meta.vaccinated.iloc[j])
    prob = modifier * p_baseline


    # Decide whether or not they will be symptomatic
    symptoms = np.random.random() < prob


    # Done
    return symptoms




#%% VAR symptoms_vax
#
#
# Combines the functions for reducing the prob of symptoms
#
#
symptoms_vax = {'site0': symptoms_rectal,
                'site1': symptoms_pharynx,
                'site2': symptoms_urethra}




#%% VAR duration_vax
#
#
# Defines collection fo functions for modifying the duration of infection
#
#


# Duration for rectal
def duration_vax_site0(vax_parameters, meta, j):
    dur_reduce = meta.vaccinated.iloc[j] * vax_parameters['site0_duration_mult'] + (1-meta.vaccinated.iloc[j])
    return dur_reduce


# Duration for pharyngeal
def duration_vax_site1(vax_parameters, meta, j):
    dur_reduce = meta.vaccinated.iloc[j] * vax_parameters['site1_duration_mult'] + (1-meta.vaccinated.iloc[j])
    return dur_reduce


# Duration for urogenital
def duration_vax_site2(vax_parameters, meta, j):
    dur_reduce = meta.vaccinated.iloc[j] * vax_parameters['site2_duration_mult'] + (1-meta.vaccinated.iloc[j])
    return dur_reduce


# Combine
duration_vax = {'site0': duration_vax_site0,
                'site1': duration_vax_site1,
                'site2': duration_vax_site2}




#%% FUN vax_at_debut()
#
#
# vaccinate everyone at the age of 16
#
#
def vax_at_debut(vax_parameters, meta, t):
    
    
    # Find everyone who is 16
    vax = meta.index[meta.age <= (16 + 1/365)]
    
    
    # Vaccinate them
    meta = give_vaccine(vax_parameters, meta, t, vax)
    
    
    # Done
    return meta




#%% FUN vax_by_targets()
#
#
# vaccinate everyone according to age group targets
#
#
def vax_by_targets(vax_parameters, meta, t):
    
    
    # Iterate over age groups
    for age in [0, 1, 2, 3, 4]:
        
        
        # Pull out that age group
        who = meta.loc[meta.age_group == age, :]
        
        
        # Compute proportion vaccinated
        unvax = np.isinf(who.vaccination_t0)
        prop = np.sum(unvax)/len(who)
        
        
        # Target to vaccinate
        target = vax_parameters['prop_vax_' + str(age)] - prop
        target = max(0, target/14)
        
        
        # Sample people to vaccinate
        temp = unvax.index[np.random.random(len(unvax)) < target]
        if age == 0:
            vax = temp
        else:
            vax = np.append(vax, temp)
    
    
    # Give the vaccines
    meta = give_vaccine(vax_parameters, meta, t, vax)
    
    
    # Done
    return meta



























#%% FUN GRAPHS initilise_tracking_data()
#
#
def initilise_tracking_data(n):
    yt = np.zeros((n, 7))
    # popt = np.zeros((n, 5))
    # part = np.zeros((n, 2))
    inft = np.zeros((n, 8))
    prvt = np.zeros((n, 13))
    tstt = np.zeros((n, 9))
    # Initlise variables for making graphs
    # if make_graphs:
    #     partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.initilise_partnership_trackers(n_steps)
    #     compartments_demo, import_count, export_count, demographic, import_time, active_age = demo.initilise_demographic_trackers(n_steps, t0_sim = t0_sim)
    #     inf_tracker = ng.initilise_infection_trackers(n_steps)
    inf_tracker = {'yt': yt, 'inft': inft, 'prvt': prvt, 'tstt': tstt}
    return inf_tracker
    
    

#%% FUN GRAPHS update_tracking_data()
#
#
# A function which updates the arrays used to track transmission
#
#
def update_tracking_data(meta, t0_sim, t, inf_tracker):
    
    
    # Extract trackers
    yt = inf_tracker['yt']
    inft = inf_tracker['inft']
    prvt = inf_tracker['prvt']
    tstt = inf_tracker['tstt']
    
    
    # Compute the current population size
    n = len(meta)
    
    
    # Compute proportion of population in each category
    yt[t,:] = [sum(meta.state == 'S')/n,
               sum(meta.state == 'E')/n,
               sum(meta.state == 'I')/n,
               sum(meta.state == 'R')/n,
               sum(meta.state == 'T')/n,
               sum(meta.vaccinated)/n,
               sum(meta.boosted)/n]
    
    
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
    tested = meta.test_time_last >= (t+t0_sim-365)
    tstt[t,:] = [np.sum(tested) / n,
                 np.sum((males & tested) & (meta.test_reason_last == int(1))) / np.sum(males),
                 np.sum((males & tested) & (meta.test_reason_last == int(2))) / np.sum(males),
                 np.sum((males & tested) & (meta.test_reason_last == int(3))) / np.sum(males),
                 np.sum(males & tested) / np.sum(males),
                 np.sum((~males & tested) & (meta.test_reason_last == int(1))) / np.sum(~males),
                 np.sum((~males & tested) & (meta.test_reason_last == int(2))) / np.sum(~males),
                 np.sum((~males & tested) & (meta.test_reason_last == int(3))) / np.sum(~males),
                 np.sum(~males & tested) / np.sum(~males),
                 ]
    
    
    # Rebundle trackers
    inf_tracker.update({'yt': yt})
    inf_tracker.update({'inft': inft})
    inf_tracker.update({'prvt': prvt})
    inf_tracker.update({'tstt': tstt})
    
    
    # Done
    return inf_tracker




#%% FUN GRAPHS make_tracking_graphs()
#
#
# Makes graphs
#
#
def make_tracking_graphs(tt, sim, inf_tracker, out_dir):
    
    
    # Extract trackers
    yt = inf_tracker['yt']
    inft = inf_tracker['inft']
    prvt = inf_tracker['prvt']
    tstt = inf_tracker['tstt']
    

    # Plot aggregate infection levels
    new_fig = plt.figure(figsize = [12, 8])
    yt = 100 * yt
    plt.plot(tt, yt[:,0], label = 'S')
    plt.plot(tt, yt[:,1], label = 'E')
    plt.plot(tt, yt[:,2], label = 'I')
    plt.plot(tt, yt[:,3], label = 'R')
    plt.plot(tt, yt[:,4], label = 'T')
    plt.plot(tt, yt[:,5], label = 'V')
    plt.plot(tt, yt[:,6], label = 'B')
    plt.title('Proportion of the Population in Each Model Compartment (Parameter Set ' + str(sim) + ')')
    plt.legend()
    plt.xlabel('Global Simulation Timestep')
    plt.ylabel('Proportion of Population (%)')
    plt.savefig(out_dir + 'graph_aggregate_infections.png')
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
    inft = 100 * inft
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
    plt.savefig(out_dir + 'graph_anatomical_site.png')
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
    plt.savefig(out_dir + 'graph_prevalence_aggregate.png')
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
    plt.savefig(out_dir + 'graph_prevalence_age_group.png')
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
    plt.savefig(out_dir + 'graph_testing_rates.png')
    plt.close()


    # Done
    return None




# #%% NG MOD progress_state_of_infection()
# # PROGRESS THE STATE OF INFECTION
# # This is almost identical to the version ng.progress_state_of_infection()
# # except that it runs vac.progress_state_of_vaccination() at the end.
# # I've removed this from the code but not deleted it yet.
# #
# #
# # Simply checks the duration of each compartment and progresses the
# # individual if the duration is up.
# #
# #
# # INPUT
# #   meta, t
# #
# # OUTPUT
# #   meta
# #
# #
# def progress_state_of_infection(meta, t):


#     # Test for new rectal infections
#     infectious = (meta["site0_t0"] < t) & (meta["site0"] == 0)
#     meta.loc[infectious, "state"] = "I"
#     meta.loc[infectious, "site0"] = 1


#     # Test for new urethral infections
#     infectious = (meta["site1_t0"] < t) & (meta["site1"] == 0)
#     meta.loc[infectious, "state"] = "I"
#     meta.loc[infectious, "site1"] = 1


#     # Test for new pharyngeal infections
#     infectious = (meta["site2_t0"] < t) & (meta["site2"] == 0)
#     meta.loc[infectious, "state"] = "I"
#     meta.loc[infectious, "site2"] = 1


#     # Remove expired rectal infections
#     recovery0 = meta["site0_t1"] <= t
#     meta.at[recovery0, "site0"] = 0
#     meta.at[recovery0, "site0_t0"] = float("inf")
#     meta.at[recovery0, "site0_t1"] = float("inf")
#     meta.at[recovery0, 'site0_symptoms'] = False


#     # Remove expired urethral infections
#     recovery1 = meta["site1_t1"] <= t
#     meta.at[recovery1, "site1"] = 0
#     meta.at[recovery1, "site1_t0"] = float("inf")
#     meta.at[recovery1, "site1_t1"] = float("inf")
#     meta.at[recovery1, 'site1_symptoms'] = False


#     # Remove expired pharengynal infections
#     recovery2 = meta["site2_t1"] <= t
#     meta.at[recovery2, "site2"] = 0
#     meta.at[recovery2, "site2_t0"] = float("inf")
#     meta.at[recovery2, "site2_t1"] = float("inf")
#     meta.at[recovery2, 'site2_symptoms'] = False


#     # Check if anybody has achieved a natural recovery
#     natural_recovery = ((recovery0 | recovery1) | recovery2) & (meta.site0==0) & (meta.site1==0) & (meta.site2==0)
#     meta.at[natural_recovery, 'state'] = 'S'


#     # Identify those where treatment has worn off
#     waned_immunity = meta["recovery_time"] < t
#     meta.loc[waned_immunity, "state"] = "S"
#     meta.loc[waned_immunity, "recovery_time"] = float("inf")


#     # Progress the state of vaccinated people
#     meta = progress_state_of_vaccination(meta, t)


#     return meta