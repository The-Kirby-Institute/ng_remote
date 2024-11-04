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


    # Overwrite testing rates
    t_max = np.max(meta.test_time_last)
    for i in range(0, len(pop_parameters['testing_rates'])):
        age = pop_parameters['testing_rates'].age_group.iloc[i]
        gender = pop_parameters['testing_rates'].gender.iloc[i]
        who = meta.loc[(meta.age_group == age) & (meta.gender == gender), :]
        who_test_background = who.loc[(who.test_reason_last == int(2)) & (t_max - who.test_time_last <= 365), :]
        pop_parameters['testing_rates'].loc[i, 'prob'] = len(who_test_background)/len(who)


    # Setup simulation runtime
    n_steps = int(365*vax_parameters['n_years']) if 'n_years' in vax_parameters else  sim_parameters.simulation_length[0]
    t0_sim = sim_parameters.partner_burn_in[0] + sim_parameters.simulation_length[0]
    tt = range(t0_sim, t0_sim + n_steps)


    # Check to see if the output directory is there
    out_dir = 'simulations/vaccination/' + vax_parameters['vax_scenario']
    if os.path.isdir(out_dir) == False:
        os.makedirs(out_dir)
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
        inf_tracker = update_tracking_data(meta, partner_matrix, t0_sim, t - t0_sim, inf_tracker)


    # Progress
    ( print('Tidying up...') if run_mode == 'serial' else [] )


    # Make ouput graphs
   # ( print('Making summary graphs') if run_mode == 'serial' else [] )
    #make_tracking_graphs(pop_parameters, tt, sim_no, inf_tracker, out_dir)


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
    #np.save(out_dir + 'output_prevalence.npy', inf_tracker['prvt'])
    #np.save(out_dir + 'output_partner_matrix.npy', partner_matrix)
    #np.save(out_dir + 'output_partner_expire.npy', partner_expire)
    #np.save(out_dir + 'output_sites.npy', inf_tracker['sites'])
    meta.to_feather(out_dir + 'output_meta.ftr')

    meta.orientation.to_csv(out_dir + 'output_meta_orientation.csv')

    pv = pd.DataFrame(inf_tracker['prvt'])
    pv.to_csv(out_dir + "output_prevalence.csv")

    s = pd.DataFrame(inf_tracker['sites'])
    s.to_csv(out_dir + "output_prevalence_by_sites.csv")

    prop_pop = pd.DataFrame(inf_tracker['pop_size'])
    prop_pop.to_csv(out_dir + "output_pop_size.csv")

    pv_uro = pd.DataFrame(inf_tracker['inft'])
    pv_uro.to_csv(out_dir + "output_prevalence_by_sites_aggregated.csv")

    inf_states = pd.DataFrame(inf_tracker['states'])
    inf_states.to_csv(out_dir + "output_infection_vax_states.csv")#by stratification

    vax_protected = pd.DataFrame(inf_tracker['vax'])
    vax_protected.to_csv(out_dir + "output_vaccinated.csv")

    Infc_states = pd.DataFrame(inf_tracker['yt'])
    Infc_states.to_csv(out_dir + "Infection_states.csv")#overall pop

    testing_info = pd.DataFrame(inf_tracker['testing'])
    testing_info.to_csv(out_dir + "Testing.csv")

    testing_aggregated = pd.DataFrame(inf_tracker['tstt'])
    testing_aggregated.to_csv(out_dir + " Testing_rates_aggregated_and_by_pop.csv")

    partner_info = pd.DataFrame(inf_tracker['partners'])
    partner_info.to_csv(out_dir + " Partner_numbers.csv")


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
# Function was retired because the who things became unnecessary when I
# finished witht he testing phase.
#
#
# Top-level function which handles how the vaccine parameters
# are interpreted by the simulation code
#
#
def set_function_for_updating_infections(vax_parameters):


    ##########################################
    ##  ALL EFFECTS WITH DEPLOYMENT MODE 3  ##
    ##########################################


    # Assume all vaccine scenarios run from the multi-effect channels
    if (vax_parameters['deployment'] == 3) &  (vax_parameters['effect'] == 4):


            # This is the sequencing of code which will result in the desired
            # effects being implemented
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


    # Put in an error check just for fun
    else:
        print('ERROR: Please implement vaccine deployment and effects through options 4 and 3, respectively.')
        update_infections = []


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
    # who are tested in response
    meta, vax = ng.seek_treatment(pop_parameters, parameters, meta, partner_matrix, t, return_tested_ids=True, project=True)


    # But not all of them
    u = np.random.random(len(vax)) <= vax_parameters['p_test_to_vax']


    # Implement vaccinations
    if len(vax) > 0:
        meta = give_vaccine(vax_parameters, meta, t, vax[u])
        meta.loc[vax[u], 'vaccine_source'] = int(1)


    # Done
    return meta




#%% FUN give_vaccine()
#
#
# What happens when somebody is vaccinated
#
#
def give_vaccine(vax_parameters, meta, t, vax):


    # Must be unvaccinated
    vax = vax[np.isinf(meta.vaccination_t0[vax])]


    # Decide which vaccinations are effective
    uu = np.random.random(len(vax))
    uu = uu < vax_parameters['prop_effective']


    # Update their vaccination data
    meta.loc[vax[uu], 'vaccinated'] = True
    meta.loc[vax[uu], 'vaccination_t0'] = t
    meta.loc[vax[uu], 'vaccination_t1'] = t + vaccination_duration(vax_parameters, len(vax[uu]))


    # Set booster time
    meta.loc[vax[uu], 'booster_t0'] = t + vax_parameters['booster_delay']
    meta.loc[vax[uu], 'booster_t1'] = t + vax_parameters['booster_delay'] + vaccination_duration(vax_parameters, len(vax[uu]))


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


    # But only vaccinate a certain proportion
    u = np.random.random(len(vax)) <= vax_parameters['p_vax_16']


    # Vaccinate them
    meta = give_vaccine(vax_parameters, meta, t, vax[u])
    meta.loc[vax[u], 'vaccine_source'] = int(2)


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
        #prop = 1 - (np.sum(unvax)/len(who))

        # Target to vaccinate
        target = vax_parameters['prop_vax_' + str(age)]/365
        #target = vax_parameters['prop_vax_' + str(age)] - prop


        # Sample people to vaccinate
        temp = unvax.index[np.random.random(len(unvax)) < target]
        if age == 0:
            vax = temp
        else:
            vax = np.append(vax, temp)


    # Give the vaccines
    if t > (15*365) + (365*20):
        meta = give_vaccine(vax_parameters, meta, t, vax)
        meta.loc[vax, 'vaccine_source'] = int(3)
    # meta = give_vaccine(vax_parameters, meta, t, vax)
    # meta.loc[vax, 'vaccine_source'] = int(3)

    # Done
    return meta




#%% FUN GRAPHS initilise_tracking_data()
#
#
def initilise_tracking_data(n):
    yt = np.zeros((n, 59))
    # popt = np.zeros((n, 5))
    # part = np.zeros((n, 2))
    inft = np.zeros((n, 8))
    prvt = np.zeros((n, 30))
    tstt = np.zeros((n, 9))

    pop_size = np.zeros((n, 10*3))
    partners = np.zeros((n, 10*4))
    states = np.zeros((n, 10*6))
    sites = np.zeros((n, 10*8))
    testing = np.zeros((n, 10*4))
    vax = np.zeros((n, 10*4))
    # Initlise variables for making graphs
    # if make_graphs:
    #     partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot = prt.initilise_partnership_trackers(n_steps)
    #     compartments_demo, import_count, export_count, demographic, import_time, active_age = demo.initilise_demographic_trackers(n_steps, t0_sim = t0_sim)
    #     inf_tracker = ng.initilise_infection_trackers(n_steps)
    inf_tracker = {'yt': yt, 'inft': inft, 'prvt': prvt, 'tstt': tstt, 'pop_size': pop_size, 'partners': partners, 'states': states, \
                   'sites': sites, 'testing': testing, 'vax': vax}
    return inf_tracker



#%% FUN GRAPHS update_tracking_data()
#
#
# A function which updates the arrays used to track transmission
#
#
def update_tracking_data(meta, partner_matrix, t0_sim, t, inf_tracker):


    # Extract trackers
    yt = inf_tracker['yt']
    inft = inf_tracker['inft']
    prvt = inf_tracker['prvt']
    tstt = inf_tracker['tstt']
    pop_size = inf_tracker['pop_size']
    partners = inf_tracker['partners']
    states = inf_tracker['states']
    sites = inf_tracker['sites']
    testing = inf_tracker['testing']
    vax = inf_tracker['vax']


    ##  UPDATE AGGREGATE TRACKERS


    # Compute the current population size
    n = len(meta)
    males = meta.gender == 1
    S = meta.state == 'S'
    E = meta.state == 'E'
    I = meta.state == 'I'
    T = meta.state == 'T'
    V = meta.vaccinated

    # Aggregate infectious state tracker
    yt[t,:] = [sum(meta.state == 'S')/n,
               sum(meta.state == 'E')/n,
               sum(meta.state == 'I')/n,
               sum(meta.state == 'R')/n,
               sum(meta.state == 'T')/n,
               sum(meta.vaccinated)/n,
               sum(meta.boosted)/n,
               sum(~(meta.vaccinated))/n,
               sum((S | E | I |T) & ~V)/n,
               n,
               sum(males),
               sum(~males),
               sum(I & V),
               sum(I & ~V),
               sum(V),
               sum(~V),
               sum(meta.age_group == 0),
               sum(meta.age_group == 1),
               sum(meta.age_group == 2),
               sum(meta.age_group == 3),
               sum(meta.age_group == 4),
               sum(V & (meta.age_group == 0)),
               sum(V & (meta.age_group == 1)),
               sum(V & (meta.age_group == 2)),
               sum(V & (meta.age_group == 3)),
               sum(V & (meta.age_group == 4)),
               sum(~V & (meta.age_group == 0)),
               sum(~V & (meta.age_group == 1)),
               sum(~V & (meta.age_group == 2)),
               sum(~V & (meta.age_group == 3)),
               sum(~V & (meta.age_group == 4)),
               sum((meta.age_group == 0) & males),
               sum((meta.age_group == 1) & males),
               sum((meta.age_group == 2) & males),
               sum((meta.age_group == 3) & males),
               sum((meta.age_group == 4) & males),
               sum((meta.age_group == 0) & ~males),
               sum((meta.age_group == 1) & ~males),
               sum((meta.age_group == 2) & ~males),
               sum((meta.age_group == 3) & ~males),
               sum((meta.age_group == 4) & ~males),
               sum(V & I & (meta.age_group == 0)),
               sum(V & I & (meta.age_group == 1)),
               sum(V & I & (meta.age_group == 2)),
               sum(V & I & (meta.age_group == 3)),
               sum(V & I & (meta.age_group == 4)),
               sum(~V & I & (meta.age_group == 0)),
               sum(~V & I & (meta.age_group == 1)),
               sum(~V & I & (meta.age_group == 2)),
               sum(~V & I & (meta.age_group == 3)),
               sum(~V & I & (meta.age_group == 4)),
               sum(V & I & (males)),
               sum(V & males),
               sum(~V & I & (males)),
               sum(~V & males),
               sum(V & I & (~males)),
               sum(V & (~males)),
               sum(~V & I & (~males)),
               sum(~V & (~males))]

    #print(sum((meta.age >= 16) & (meta.age <17)))

    # Aggregate anatomical-site specific infection tracker
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


    # Prevalence aggregate and by sub-population
    I = meta.state == 'I'
    V = meta.vaccinated

    pp = (sum(meta.age_group == 1) + sum(meta.age_group == 2) + sum(meta.age_group == 3) + sum(meta.age_group == 4))
    prvt[t,:] = [sum(I & V)/(max(1,sum(V))),
                 sum(I & ~V)/(sum(~V)),
                 sum(I)/n,
                 sum(I & (males))/sum(males),
                 sum(I & (~males))/sum(~males),
                 sum(I & (males) & (meta.age_group == 0))/max(1, sum((males) & (meta.age_group == 0))),
                 sum(I & (males) & (meta.age_group == 1))/max(1, sum((males) & (meta.age_group == 1))),
                 sum(I & (males) & (meta.age_group == 2))/max(1, sum((males) & (meta.age_group == 2))),
                 sum(I & (males) & (meta.age_group == 3))/max(1, sum((males) & (meta.age_group == 3))),
                 sum(I & (males) & (meta.age_group == 4))/max(1, sum((males) & (meta.age_group == 4))),
                 sum(I & (~males) & (meta.age_group == 0))/max(1, sum((~males) & (meta.age_group == 0))),
                 sum(I & (~males) & (meta.age_group == 1))/max(1, sum((~males) & (meta.age_group == 1))),
                 sum(I & (~males) & (meta.age_group == 2))/max(1, sum((~males) & (meta.age_group == 2))),
                 sum(I & (~males) & (meta.age_group == 3))/max(1, sum((~males) & (meta.age_group == 3))),
                 sum(I & (~males) & (meta.age_group == 4))/max(1, sum((~males) & (meta.age_group == 4))),
                 sum(V & I & (meta.age_group == 0))/max(1, sum((meta.age_group == 0) & V)),
                 sum(V & I & (meta.age_group == 1))/max(1, sum((meta.age_group == 1) & V)),
                 sum(V & I & (meta.age_group == 2))/max(1, sum((meta.age_group == 2) & V)),
                 sum(V & I & (meta.age_group == 3))/max(1, sum((meta.age_group == 3) & V)),
                 sum(V & I & (meta.age_group == 4))/max(1, sum((meta.age_group == 4) & V)),
                 sum(~V & I & (meta.age_group == 0))/max(1, sum((meta.age_group == 0) & ~V)),
                 sum(~V & I & (meta.age_group == 1))/max(1, sum((meta.age_group == 1) & ~V)),
                 sum(~V & I & (meta.age_group == 2))/max(1, sum((meta.age_group == 2) & ~V)),
                 sum(~V & I & (meta.age_group == 3))/max(1, sum((meta.age_group == 3) & ~V)),
                 sum(~V & I & (meta.age_group == 4))/max(1, sum((meta.age_group == 4) & ~V)),
                 sum(I & (meta.age_group == 0))/sum(meta.age_group == 0),
                 (sum(I & (meta.age_group == 1)) + sum(I & (meta.age_group == 2)) + sum(I & (meta.age_group == 3))+ sum(I & (meta.age_group == 4)))/pp,
                 sum(meta.orientation == 0),
                 sum(meta.orientation == 1),
                 sum(meta.orientation == 2)]


    # Testing rates aggregate and by sub-population
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


    ##  UPDATE SUB-COHORT TRACKERS


    # Initilise indexes
    ii_pop_size = 0
    ii_partners = 0
    ii_states = 0
    ii_sites = 0
    ii_test = 0
    ii_vax = 0


    # Population size
    for a in [0, 1, 2, 3, 4]:
        for g in [0, 1]:


            # Extract relevent population
            who = meta.loc[(meta.age_group == a) & (meta.gender == g), :]
            n = len(who)


            # Compute population size
            pop_size[t, ii_pop_size] = n
            pop_size[t, ii_pop_size + 1] = np.sum((who.import_time <= t0_sim))
            pop_size[t, ii_pop_size + 2] = np.sum((who.import_time > t0_sim))
            ii_pop_size = ii_pop_size + 3


            # Compute the number of partners
            pmatrow = np.sum(partner_matrix[who.index, 0:len(meta)], axis = 1)
            partners[t, ii_partners] = np.sum(pmatrow == 0)/n
            partners[t, ii_partners + 1] = np.sum((pmatrow == 1) & (who.partner != -1))/n
            partners[t, ii_partners + 2] = np.sum((pmatrow == 1) & (who.partner == -1))/n
            partners[t, ii_partners + 3] = np.sum(pmatrow > 1)/n
            ii_partners = ii_partners + 4


            # Compute infectious state
            states[t, ii_states] = np.sum(who.state == 'S')/n
            states[t, ii_states + 1] = np.sum(who.state == 'E')/n
            states[t, ii_states + 2] = np.sum(who.state == 'I')/n
            states[t, ii_states + 3] = np.sum(who.state == 'T')/n
            states[t, ii_states + 4] = np.sum(who.vaccinated)/n
            states[t, ii_states + 5] = np.sum(who.boosted)/n
            ii_states = ii_states + 6


            # Compute infection prevalence by anatomical-site
            sites[t, ii_sites] = np.sum(who.site0 & ~who.site1 & ~who.site2)/n
            sites[t, ii_sites + 1] = np.sum(~who.site0 & who.site1 & ~who.site2)/n
            sites[t, ii_sites + 2] = np.sum(~who.site0 & ~who.site1 & who.site2)/n
            sites[t, ii_sites + 3] = np.sum(who.site0 & who.site1 & ~who.site2)/n
            sites[t, ii_sites + 4] = np.sum(who.site0 & ~who.site1 & who.site2)/n
            sites[t, ii_sites + 5] = np.sum(~who.site0 & who.site1 & who.site2)/n
            sites[t, ii_sites + 6] = np.sum(who.site0 & who.site1 & who.site2)/n
            ii_sites = ii_sites + 7


            # Compute testing rates
            last_365 = who.loc[(t + t0_sim - who.test_time_last) <= 365, ]
            testing[t, ii_test] = (len(last_365))/n
            testing[t, ii_test + 1] = np.sum(last_365.test_reason_last == int(1))/n
            testing[t, ii_test + 2] = np.sum(last_365.test_reason_last == int(2))/n
            testing[t, ii_test + 3] = np.sum(last_365.test_reason_last == int(3))/n
            ii_test = ii_test + 4


            # Compute vaccination rates by source
            vax[t, ii_vax] = np.sum((who.vaccination_t0 <= (t0_sim + t)) | (who.booster_t0 <= (t0_sim + t)))/n
            vax[t, ii_vax + 1] = np.sum(((who.vaccination_t0 <= (t0_sim + t)) | (who.booster_t0 <= (t0_sim + t))) & (who.vaccine_source == int(1)))/n
            vax[t, ii_vax + 2] = np.sum(((who.vaccination_t0 <= (t0_sim + t)) | (who.booster_t0 <= (t0_sim + t))) & (who.vaccine_source == int(2)))/n
            vax[t, ii_vax + 3] = np.sum(((who.vaccination_t0 <= (t0_sim + t)) | (who.booster_t0 <= (t0_sim + t))) & (who.vaccine_source == int(3)))/n
            ii_vax = ii_vax + 4


    # Rebundle trackers
    inf_tracker.update({'yt': yt})
    inf_tracker.update({'inft': inft})
    inf_tracker.update({'prvt': prvt})
    inf_tracker.update({'tstt': tstt})
    inf_tracker.update({'pop_size': pop_size})
    inf_tracker.update({'partners': partners})
    inf_tracker.update({'states': states})
    inf_tracker.update({'sites': sites})
    inf_tracker.update({'testing': testing})
    inf_tracker.update({'vax': vax})


    # Done
    return inf_tracker




## %% FUN GRAPHS make_tracking_graphs()-----REMOVED


## Makes graphs



#



#%% REMOVED set_function_for_updating_infections()
# Function was retired because the who things became unnecessary when I
# finished witht he testing phase.
#
#
# Top-level function which handles how the vaccine parameters
# are interpreted by the simulation code
#
#
# def set_function_for_updating_infections(vax_parameters):


#     ##########################################
#     ##  ALL EFFECTS WITH DEPLOYMENT MODE 0  ##
#     ##########################################


#     # Where distribution occurrs during treatment
#     if vax_parameters['deployment'] == 0:


#         # Where the effect is total immunity
#         if vax_parameters['effect'] == 0:


#             # Situation where vaccinations given during treatment and makes people immune
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_infection_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce transmission
#         elif vax_parameters['effect'] == 1:


#             # Situation where vaccines given during treatment and result is to
#             # decrease transmission probability
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce symptoms
#         elif vax_parameters['effect'] == 2:


#             # Situation where vaccines given during treatment and result in
#             # people being more likely to be asymptomatic
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          symptoms_prob = symptoms_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is the reduce the duration of infection
#         elif vax_parameters['effect'] == 3:


#             # Situation where vaccines given during treatment and result in
#             # people having a shorter duration of infection
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect includes all of the above
#         elif vax_parameters['effect'] == 4:


#             # Situation where vaccines have a combination of all of the
#             # above effects
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_and_inf_prob,
#                                          symptoms_prob = symptoms_vax,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#     ##########################################
#     ##  ALL EFFECTS WITH DEPLOYMENT MODE 1  ##
#     ##########################################


#     # Where distribution occurs at the age of 16
#     elif vax_parameters['deployment'] == 1:


#         # Where the effect is total immunity
#         if vax_parameters['effect'] == 0:


#             # Situation where vaccinations given during treatment and makes people immune
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_infection_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce transmission
#         elif vax_parameters['effect'] == 1:


#             # Situation where vaccines given during treatment and result is to
#             # decrease transmission probability
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce symptoms
#         elif vax_parameters['effect'] == 2:


#             # Situation where vaccines given during treatment and result in
#             # people being more likely to be asymptomatic
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          symptoms_prob = symptoms_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is the reduce the duration of infection
#         elif vax_parameters['effect'] == 3:


#             # Situation where vaccines given during treatment and result in
#             # people having a shorter duration of infection
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect includes all of the above
#         elif vax_parameters['effect'] == 4:


#             # Situation where vaccines have a combination of all of the
#             # above effects
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_and_inf_prob,
#                                          symptoms_prob = symptoms_vax,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#     ##########################################
#     ##  ALL EFFECTS WITH DEPLOYMENT MODE 2  ##
#     ##########################################


#     # Where distribution occurs at the age of 16
#     elif vax_parameters['deployment'] == 2:


#         # Where the effect is total immunity
#         if vax_parameters['effect'] == 0:


#             # Situation where vaccinations given during treatment and makes people immune
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_infection_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce transmission
#         elif vax_parameters['effect'] == 1:


#             # Situation where vaccines given during treatment and result is to
#             # decrease transmission probability
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce symptoms
#         elif vax_parameters['effect'] == 2:


#             # Situation where vaccines given during treatment and result in
#             # people being more likely to be asymptomatic
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          symptoms_prob = symptoms_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is the reduce the duration of infection
#         elif vax_parameters['effect'] == 3:


#             # Situation where vaccines given during treatment and result in
#             # people having a shorter duration of infection
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect includes all of the above
#         elif vax_parameters['effect'] == 4:


#             # Situation where vaccines have a combination of all of the
#             # above effects
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_and_inf_prob,
#                                          symptoms_prob = symptoms_vax,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
#                 return meta


#     ##########################################
#     ##  ALL EFFECTS WITH DEPLOYMENT MODE 3  ##
#     ##########################################


#     # Where distribution occurs at the age of 16
#     elif vax_parameters['deployment'] == 3:


#         # Where the effect is total immunity
#         if vax_parameters['effect'] == 0:


#             # Situation where vaccinations given during treatment and makes people immune
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_infection_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce transmission
#         elif vax_parameters['effect'] == 1:


#             # Situation where vaccines given during treatment and result is to
#             # decrease transmission probability
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_prob,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is to reduce symptoms
#         elif vax_parameters['effect'] == 2:


#             # Situation where vaccines given during treatment and result in
#             # people being more likely to be asymptomatic
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          symptoms_prob = symptoms_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect is the reduce the duration of infection
#         elif vax_parameters['effect'] == 3:


#             # Situation where vaccines given during treatment and result in
#             # people having a shorter duration of infection
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#         # Where the effect includes all of the above
#         elif vax_parameters['effect'] == 4:


#             # Situation where vaccines have a combination of all of the
#             # above effects
#             def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
#                 meta = vax_at_debut(vax_parameters, meta, t)
#                 meta = vax_by_targets(vax_parameters, meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
#                                          trans_prob_fun = vaccine_reduced_trans_and_inf_prob,
#                                          symptoms_prob = symptoms_vax,
#                                          duration_mod = duration_vax,
#                                          vax_parameters = vax_parameters)
#                 meta = ng.progress_state_of_infection(meta, t)
#                 meta = progress_state_of_vaccination(meta, t)
#                 meta = seek_treatment(pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix, t)
#                 return meta


#     # Done
#     return update_infections


#%% REMOVED progress_state_of_infection()
# Function was retired by some nifty coding, but it's here in case I still need it.
#
# PROGRESS THE STATE OF INFECTION
# This is almost identical to the version ng.progress_state_of_infection()
# except that it runs vac.progress_state_of_vaccination() at the end.
# I've removed this from the code but not deleted it yet.
#
#
# Simply checks the duration of each compartment and progresses the
# individual if the duration is up.
#
#
# INPUT
#   meta, t
#
# OUTPUT
#   meta
#
#
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