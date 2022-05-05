# -*- coding: utf-8 -*-
"""
CODE FOR DOING ALL OF THE INFECTION STUFF

Created on Thu Nov 12 13:00:07 2020

@author: nicol
"""


# %% Setup Libraries
import numpy as np
import pandas as pd
import scipy.stats as sp
from scipy.stats import gamma
import matplotlib.pyplot as plt


# %% FUN update_infections()
#
#
# Function to update the state of infections
#
#
def update_infections(pop_parameters, inf_parameters, meta, partner_matrix, t):

    
    # Implement a transmission event
    meta = new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t)


    # Update infectious states
    meta = progress_state_of_infection(meta, t)


    # Implement treatment
    meta = seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)


    return meta


# %% FUN duration_rectal()
# FUNCTION FOR DURATION OF RECTAL INFECTIONS
#
#
# Latent period
#    constant
#
# Duration of rectal infection
#    Infections can be symptomatic or asymptomatic.
#    Duration is drawn from a Gamma distribution with specified mean and var
#
# Duration of urethral infection
#    Infections can be symptomatic or asymptomatic
#    The probability of a symptomatic infection depends on gender
#    Duration is drawn from a Gamma distribution with specified mean and var
#    Mean and var is the same for both genders
#
# Duration of pharyngeal infection
#    All infections assumed to be asymptomatic
#
#
# INPUT
#    infectee = meta.loc[infectee] - used in some functions for convenience
#
#    mean_rectal
#    var_rectal
#    symptomatic_rectal
#
#    mean_urethral
#    var_urethral
#    symptomatic_urethral - 1x2 list, columns correspond to gender
#
#    mean_phar
#    var_rectal
#
#
# OUTPUT
#    duration
#
#
def duration_rectal(inf_parameters, infectee, n=1):

    # Simulate duration of natural infection
    duration = np.random.gamma(inf_parameters['infection'].mean_rectal[infectee.gender] /
                               inf_parameters['infection'].var_rectal[infectee.gender],
                               inf_parameters['infection'].var_rectal[infectee.gender],
                               n)

    # Return duration
    return duration


# %% FUN duration_urethral()
def duration_urethral(inf_parameters, infectee, n=1):

    # Simulate duration of natural infection
    duration = np.random.gamma(inf_parameters['infection'].mean_urethral[infectee.gender] /
                               inf_parameters['infection'].var_urethral[infectee.gender],
                               inf_parameters['infection'].var_urethral[infectee.gender],
                               n)

    # Return duration
    return duration


# %% FUN duration_pharyngeal()
def duration_pharyngeal(inf_parameters, infectee, n=1):

    # Simulate duration of natural infection
    duration = np.random.gamma(inf_parameters['infection'].mean_pharyngeal[infectee.gender] /
                               inf_parameters['infection'].var_pharyngeal[infectee.gender],
                               inf_parameters['infection'].var_pharyngeal[infectee.gender],
                               n)

    # Return duration
    return duration


# %% FUN duration_treatment_immunity()
def duration_treatment_immunity(parameters, treat):

    # All treatment results results in a Gamma distributed immune period
    duration = np.random.gamma(parameters['infection'].immunity_mean[0] /
                               parameters['infection'].immunity_var[0],
                               parameters['infection'].immunity_var[0], len(treat))

    # Return duration
    return duration


# %% FUN latent_period()
def latent_period(inf_parameters, infectee):

    # Just using a constant latent period for now
    out = inf_parameters['infection'].latent_period[infectee.gender]
    out = np.repeat(out, 1)

    return out


# %% FUN delay_test_to_treatment()
def delay_test_to_treatment(pop_parameters, n_sim):

    # Simulate delay
    # u = np.random.random(n_sim)
    # treat_delay = \
    #     pop_parameters['testing_param'].delay_x1_baseline[0] + \
    #     np.log(u/(1-u))/pop_parameters['testing_param'].delay_x0_baseline[0]
    # treat_delay[treat_delay < 0] = 0

    # David just wants POC testing / treatment
    treat_delay = np.array(n_sim*[0])

    # Done
    return treat_delay


# %% FUN transmission_probability()
#
#
# Compute the site-specific transmission probabilities
#
#
def transmission_probability(inf_parameters, vax_parameters, meta, i, j):


    # # Determine whether any sexual event occurs
    # is_long = int(meta.at[i, 'partner'] == j)
    # p_event = inf_parameters['p_event'][is_long]


    # Pull out some indicies for convenience
    g0 = int(meta.at[i, "gender"])
    g1 = int(meta.at[j, "gender"])


    # Identify which acts they can engage in
    can_oral = meta.at[i, 'has_oral'] & meta.at[j, 'has_oral']
    can_sex = meta.at[i, 'has_sex'] & meta.at[j, 'has_sex']
    can_anal = meta.at[i, 'has_anal'] & meta.at[j, 'has_anal']


    # Decide which sexual acts they engage in
    u = np.random.random(5)
    acts = np.array([[1.0 * (can_anal & (u[0] < inf_parameters['p_anal'][g0, g1])),
                      1.0 * (can_oral & (u[1] < inf_parameters['p_oral'][g0, g1])),
                      1.0 * (u[2] < inf_parameters['p_kiss'][g0, g1]),
                      1.0 * (can_anal & (u[3] < inf_parameters['p_rim'][g0, g1])),
                      1.0 * (can_sex & (u[4] < inf_parameters['p_sex'][g0, g1]))]])


    # Compose the transition probability matrix
    trans_prob = \
        acts[0, 0] * inf_parameters['trans_anal'] + \
        acts[0, 1] * inf_parameters['trans_oral'] + \
        acts[0, 2] * inf_parameters['trans_kiss'] + \
        acts[0, 3] * inf_parameters['trans_rim'] + \
        acts[0, 4] * inf_parameters['trans_sex']


    return trans_prob


# %% FUN symptoms_rectal()
#
#
def symptoms_rectal(inf_parameters, vax_parameters, meta, i, j):
    prob = np.random.random() < inf_parameters['infection'].symptoms_rectal[meta.gender.iloc[j]]
    return prob


# %% FUN symptoms_pharynx()
#
#
def symptoms_pharynx(inf_parameters, vax_parameters, meta, i, j):
    prob = np.random.random() < inf_parameters['infection'].symptoms_pharyngeal[meta.gender.iloc[j]]
    return prob


# %% FUN symptoms_urethra()
#
#
def symptoms_urethra(inf_parameters, vax_parameters, meta, i, j):
    prob = np.random.random() < inf_parameters['infection'].symptoms_urethral[meta.gender.iloc[j]]
    return prob


# %% VAR symptoms_baseline
#
#
symptoms_baseline = {'site0': symptoms_rectal,
                     'site1': symptoms_pharynx,
                     'site2': symptoms_urethra}


# %% VAR duration_baseline
#
#
def duration_baseline_site0(vax_parameters, meta, j): return 1
def duration_baseline_site1(vax_parameters, meta, j): return 1
def duration_baseline_site2(vax_parameters, meta, j): return 1


duration_baseline = {'site0': duration_baseline_site0,
                     'site1': duration_baseline_site1,
                     'site2': duration_baseline_site2}


# %% FUN new_infections()
# LOOK AT ALL PARTNERSHIPS AND SEED NEW INFECTIONS
#
#
# Uses matrices for describing the probability of engaging in certain acts as
# a function of gender. In particular, p_anal, p_oral, p_kiss, p_rim, and p_sex.
# In all cases, read entry (i, j) as the probability that an individual of
# sex i (0=F, 1=M) engaging in that particular act with an individual of sex j.
#
#

#
#
# Decision tree:
#
#    1. Find all individuals in the population who have at least 1 partner
#        and are infections.
#
#    2. For a given person i, look at all their partners j and check that:
#        a. Individual j is not immune
#        b. Individual j is not infected at all 3 sites
#
#    3. Given an infectious individual i and a potentially infectable partner j
#        use the act-specific probabilities to determine whether or not a
#        particular sexual act takes place.
#
#    4. Given a range of sexual acts, use the site-specific probabilities to
#        determine which sites are to be infected.
#
#    5. Seed those sites with an infection. They will not become infectious
#        until the latent period is over.
#
#
# INPUT
#   meta, partner_matrix, t
#   p_anal, The probability of a given sexual act (Act-specific probabilities)
#   p_oral,
#   p_kiss,
#   p_rim,
#   p_sex,
#   trans_anal, The probability of a site-to-site transmission (Site-specific probabilities)
#   trans_oral,
#   trans_kiss,
#   trans_rim,
#   trans_sex
#
# OUTPUT
#   meta
#
#
def new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t,
                   trans_prob_fun=transmission_probability,
                   vax_parameters=[],
                   symptoms_prob=symptoms_baseline,
                   duration_mod=duration_baseline):

    
    # Determine whether or not any new infections could possibly occur
    infectors = meta[(meta["state"] == "I") & (np.sum(partner_matrix[0:len(meta), 0:len(meta)], axis=1) > 0)]


    # Iterate over all infectors
    for i in infectors.index:

        
        # Extract all partners of the infector
        partners = np.asarray(np.where(partner_matrix[i, :] == 1))


        # Remove partners who are either immune or already infected at all sites
        partner_susceptible = meta.loc[partners[0, :], 'state'].isin(['S', 'E', 'I'])
        partner_site = (meta.loc[partners[0, :], "site0"] == 1) & \
                       (meta.loc[partners[0, :], "site1"] == 1) & \
                       (meta.loc[partners[0, :], "site2"] == 1)
        partners = partners[0, partner_susceptible & (~partner_site)]


        # Iterate over all partnerships with potential for transmission
        for j in partners:

            
            # Pull out some indices for convenience
            trans_prob = trans_prob_fun(inf_parameters, vax_parameters, meta, i, j)


            # Determine if any transmissions have occured this iteration
            sites = np.array(meta.loc[i, ["site0", "site1", "site2"]])
            M = sites * np.transpose(trans_prob)
            U = np.random.random((3, 3))
            N = np.sum(U < M, axis=1) > 0
            # print(M, U, N)
            new_inf = np.flatnonzero(N)


            # Make sure any new infections don't overlap with current infections
            infectee = meta.loc[j, :]
            sites = infectee[["site0", "site1", "site2"]] == 1
            exposures = infectee[["site0_t0", "site1_t0", "site2_t0"]] < float("inf")
            old_inf = np.where(sites.to_numpy() | exposures.to_numpy())
            new_inf = new_inf[~np.isin(new_inf, old_inf)]


            # Seed new infections
            if len(new_inf) > 0:

                
                # Update state of infectee
                meta.at[j, 'state'] = ('E' if infectee.state == 'S' else infectee.state)


                # Set infection parameters for the new cases
                for k in new_inf:

                    
                    # Set duration of latent period (assumed to be the same for all sites)
                    end_latent = t + latent_period(inf_parameters, infectee)
                    meta.at[j, "site" + str(int(k)) + "_t0"] = end_latent


                    # Set site-specific infection parameters
                    if k == 0:

                        
                        # Infection at rectum
                        meta.at[j, "site0_t1"] = end_latent + duration_mod['site0'](vax_parameters, meta, j) * duration_rectal(inf_parameters, infectee)
                        meta.at[j, 'site0_symptoms'] = symptoms_prob['site0'](inf_parameters, vax_parameters, meta, i, j)


                    elif k == 1:

                        
                        # Infection at pharynx
                        meta.at[j, "site1_t1"] = end_latent + duration_mod['site1'](vax_parameters, meta, j) * duration_urethral(inf_parameters, infectee)
                        meta.at[j, 'site1_symptoms'] = symptoms_prob['site1'](inf_parameters, vax_parameters, meta, i, j)


                    else:

                        
                        # Infection at urethra
                        meta.at[j, "site2_t1"] = end_latent + duration_mod['site2'](vax_parameters, meta, j) * duration_pharyngeal(inf_parameters, infectee)
                        meta.at[j, 'site2_symptoms'] = symptoms_prob['site2'](inf_parameters, vax_parameters, meta, i, j)


                    # Set delay from symptoms to treatment
                    symptoms = meta.at[j, 'site0_symptoms'] | meta.at[j, 'site1_symptoms'] | meta.at[j, 'site2_symptoms']
                    if symptoms:

                        
                        # Pull the time of test from Gamma CDF
                        test_delay = gamma.ppf(meta.at[j, 'treatment_threshold'],
                                               inf_parameters['infection'].treatment_mean[0] /
                                               inf_parameters['infection'].treatment_var[0],
                                               inf_parameters['infection'].treatment_var[0])
                        meta.loc[j, 'test_time'] = end_latent + test_delay


                        # Decide whether they'll get treated
                        p_treat = pop_parameters['testing_param'].test_sensitivity[0] * pop_parameters['testing_param'].p_treat_baseline[0]
                        if np.random.random(1)[0] < p_treat:

                            
                            # Sample a treatment delay
                            u = np.random.random(1)[0]
                            treat_delay = delay_test_to_treatment(pop_parameters, 1)
                            meta.loc[j, 'treatment_time'] = meta.loc[j, 'test_time'] + treat_delay

    # Return meta
    return meta


# %% FUN progress_state_of_infection()
# PROGRESS THE STATE OF INFECTION
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
def progress_state_of_infection(meta, t):

    # Test for new rectal infections
    infectious = (meta["site0_t0"] < t) & (meta["site0"] == 0)
    meta.loc[infectious, "state"] = "I"
    meta.loc[infectious, "site0"] = 1

    # Test for new urethral infections
    infectious = (meta["site1_t0"] < t) & (meta["site1"] == 0)
    meta.loc[infectious, "state"] = "I"
    meta.loc[infectious, "site1"] = 1

    # Test for new pharyngeal infections
    infectious = (meta["site2_t0"] < t) & (meta["site2"] == 0)
    meta.loc[infectious, "state"] = "I"
    meta.loc[infectious, "site2"] = 1

    # Remove expired rectal infections
    recovery0 = meta["site0_t1"] <= t
    meta.at[recovery0, "site0"] = 0
    meta.at[recovery0, "site0_t0"] = float("inf")
    meta.at[recovery0, "site0_t1"] = float("inf")
    meta.at[recovery0, 'site0_symptoms'] = False

    # Remove expired urethral infections
    recovery1 = meta["site1_t1"] <= t
    meta.at[recovery1, "site1"] = 0
    meta.at[recovery1, "site1_t0"] = float("inf")
    meta.at[recovery1, "site1_t1"] = float("inf")
    meta.at[recovery1, 'site1_symptoms'] = False

    # Remove expired pharengynal infections
    recovery2 = meta["site2_t1"] <= t
    meta.at[recovery2, "site2"] = 0
    meta.at[recovery2, "site2_t0"] = float("inf")
    meta.at[recovery2, "site2_t1"] = float("inf")
    meta.at[recovery2, 'site2_symptoms'] = False

    # Check if anybody has achieved a natural recovery
    natural_recovery = ((recovery0 | recovery1) | recovery2) & (meta.site0 == 0) & (meta.site1 == 0) & (meta.site2 == 0)
    meta.at[natural_recovery, 'state'] = 'S'

    # Remove treatment-conferred immunity
    waned_immunity = meta["recovery_time"] < t
    meta.loc[waned_immunity, "state"] = "S"
    meta.loc[waned_immunity, "recovery_time"] = float("inf")

    return meta


# %% FUN seek_treatment()
# FUNCTION FOR TIME UNTIL TREATMENT
#
#
# Time until treatment
#    Gamma distributed with specified mean and variance.
#    Upon treatment, the indivdual's parter will also be treated.
#    Immunity is conferred for a specied period.
#
#
# INPUT
#    meta, t
#
#    Parameters of the distribution of time until an individual seeks treatment
#    treat_mean
#    treat_var
#
#    Parameters of the distribution of time for which an individual is immune
#    immune_mean
#    immune_var
#
#
# OUTPUT
#    meta
#
#
def seek_treatment(pop_parameters, parameters, meta, partner_matrix, t, return_treated_ids = False):


    ###############
    ##  TESTING  ##
    ###############


    # SYMPTOM BASED TESTING


    # Determine if anybody is to get tested today due to their symptoms
    test_symptoms = meta.treatment_time <= t
    test_symptoms = meta.index[test_symptoms]


    # Update testing data
    meta.loc[test_symptoms, 'test_time'] = float('Inf')
    meta.loc[test_symptoms, 'test_time_last'] = t
    meta.loc[test_symptoms, 'test_reason_last'] = int(1)
    
    
    # BACKGROUND TESTING


    # Determine if anyone is to get background tested today
    test_background = np.array(list(), dtype='int64')
    rates = pop_parameters['testing_rates']
    for age in [0, 1, 2, 3, 4]:
        for gender in [0, 1]:

            
            # Count up the proportion of the population that has been tested in the last 12 months
            who = meta.loc[(meta.age_group == age) & (meta.gender == gender), :]
            last_365 = (t - who.test_time_last) <= 365
            # test_reason = who.test_reason_last != int(2)
            # prop = (np.sum(test_reason & test_time))/max(len(who), 1)
            prop = (np.sum(last_365))/max(len(who), 1)


            # Calculate the current shortfall proportion
            background = rates.prob[(rates.age_group == age) & (rates.gender == gender)].values[0] - prop
            background = max(0, background/(7))


            # Sample new people to test
            who_test = who.index[last_365 == False]
            who_test = who_test[np.random.random(len(who_test)) <= background]
            test_background = np.append(test_background, who_test)


    # Decide if any will be treated
    p_treat = pop_parameters['testing_param'].test_sensitivity[0] * pop_parameters['testing_param'].p_treat_baseline[0]
    treat_background = (meta.state[test_background] == 'I') & (np.random.random(len(test_background)) <= p_treat)
    treat_background = test_background[treat_background]


    # Update data for background testing and treatment
    meta.loc[test_background, 'test_time'] = float('Inf')
    meta.loc[test_background, 'test_time_last'] = t
    meta.loc[test_background, 'test_reason_last'] = int(2)
    meta.loc[treat_background, 'treatment_time'] = t + delay_test_to_treatment(pop_parameters, len(treat_background))


    # CONTACT TRACING


    # Pull out the partners of everybody tested today
    test = np.append(test_symptoms, test_background)
    part = meta.partner.iloc[test]
    part = part[part > -1]
    part = part.values[meta.state.iloc[part] == 'I']
    p_treat = pop_parameters['testing_param'].test_sensitivity[0] * pop_parameters['testing_param'].p_treat_partner[0]
    treat_part = part[np.random.random(len(part)) < p_treat]


    # Assume they just present for treatment one day
    treat_delay = delay_test_to_treatment(pop_parameters, len(treat_part))
    meta.loc[treat_part, 'test_reason_last'] = int(3)
    meta.loc[treat_part, 'test_time_last'] = t + treat_delay
    meta.loc[treat_part, 'treatment_time'] = t + treat_delay


    #################
    ##  TREATMENT  ##
    #################


    # Pull out everyone with treatment today
    treat = meta.index[meta.treatment_time <= t]


    # Give them treatment
    meta.loc[treat, "state"] = "T"
    meta.loc[treat, "recovery_time"] = t + duration_treatment_immunity(parameters, treat)
    meta.loc[treat, "site0"] = 0
    meta.loc[treat, "site1"] = 0
    meta.loc[treat, "site2"] = 0
    meta.loc[treat, "site0_t0"] = float("inf")
    meta.loc[treat, "site1_t0"] = float("inf")
    meta.loc[treat, "site2_t0"] = float("inf")
    meta.loc[treat, "site0_t1"] = float("inf")
    meta.loc[treat, "site1_t1"] = float("inf")
    meta.loc[treat, "site2_t1"] = float("inf")
    meta.loc[treat, 'site0_symptoms'] = False
    meta.loc[treat, 'site1_symptoms'] = False
    meta.loc[treat, 'site2_symptoms'] = False
    meta.loc[treat, 'treatment_time'] = float('Inf')


    # Added this in so that the function will return the variable treat
    # if it's needed by the vaccine code
    if return_treated_ids:
        return meta, treat
    else:
        return meta


# %% GRAPH initilise_infection_trackers
#
#
#
def initilise_infection_trackers(n_days):

    # Tracking the infectious states
    yt = pd.DataFrame(columns=["S", 'E', "I", "R", 'T', 'V',
                      "site0s", "site1s", "site2s", "site0a", "site1a", "site2a"])

    # Tracking headline infection process by age group and risk group
    it = np.zeros((n_days, 7 * 2 * 5))

    # Tracking headline infections by anatomical site
    infections = np.zeros((n_days, 7*10))

    # Prevalence trackers
    inf_state_age_risk = np.zeros((n_days, 7*10))
    inf_site_age_gender = np.zeros((n_days, 7*10))
    pop_size_age_risk = np.zeros((n_days, 5*2))
    pop_size_age_gender = np.zeros((n_days, 5*2))
    infections_site_age_gender_symptom = np.zeros((n_days, 8*10))

    # Tracking testing rates
    track_treat = np.zeros((n_days, 5*2*4))

    # Put variables into dict for tracking
    infection_tracker = {'aggregate_infections': yt,
                         'infections_state_age_gender': it,
                         'infections_state_age_risk': inf_state_age_risk,
                         'infections_site_age_gender': inf_site_age_gender,
                         'infections_site_age_risk': infections,
                         'pop_size_age_risk': pop_size_age_risk,
                         'pop_size_age_gender': pop_size_age_gender,
                         'infections_site_age_gender_symptom': infections_site_age_gender_symptom,
                         'track_treat': track_treat
                         }

    return infection_tracker


# %% GRAPH update_infection_trackers()
#
#
#
def update_infection_trackers(t, pop_parameters, meta, inf_tracker, t0_sim):

    # Update population size trackers
    temp1 = inf_tracker['pop_size_age_risk']
    temp2 = inf_tracker['pop_size_age_gender']
    ii = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            temp1[t, ii] = np.sum(pop_parameters['lookup']
                                  ['RA' + str(l2) + str(age)])
            temp2[t, ii] = np.sum(pop_parameters['lookup']
                                  ['GA' + str(l2) + str(age)])
            ii = ii + 1
    inf_tracker.update({'pop_size_age_risk': temp1})
    inf_tracker.update({'pop_size_age_gender': temp2})

    # Update aggregate infection graph
    yt = inf_tracker['aggregate_infections']
    yt.at[t, "S"] = np.sum(pop_parameters['lookup']['SS'])
    yt.at[t, "E"] = np.sum(pop_parameters['lookup']['SE'])
    yt.at[t, "I"] = np.sum(pop_parameters['lookup']['SI'])
    yt.at[t, "R"] = np.sum(pop_parameters['lookup']['SR'])
    yt.at[t, "T"] = np.sum(pop_parameters['lookup']['ST'])
    yt.at[t, "V"] = np.sum(pop_parameters['lookup']['SV'])
    yt.at[t, "site0s"] = np.sum(
        pop_parameters['lookup']['I01'] & meta.site0_symptoms)
    yt.at[t, "site1s"] = np.sum(
        pop_parameters['lookup']['I11'] & meta.site1_symptoms)
    yt.at[t, "site2s"] = np.sum(
        pop_parameters['lookup']['I21'] & meta.site2_symptoms)
    yt.at[t, "site0a"] = np.sum(
        pop_parameters['lookup']['I01'] & ~meta.site0_symptoms)
    yt.at[t, "site1a"] = np.sum(
        pop_parameters['lookup']['I11'] & ~meta.site1_symptoms)
    yt.at[t, "site2a"] = np.sum(
        pop_parameters['lookup']['I21'] & ~meta.site2_symptoms)
    inf_tracker.update({'yt': yt})

    # Update infectious state by age and gender
    temp = inf_tracker['infections_state_age_gender']
    ii = 0
    for age in [0, 1, 2, 3, 4]:
        for gender in [0, 1]:
            for state in ['S', 'E', 'I', 'R', 'T', 'V']:
                temp[t, ii] = np.sum(
                    pop_parameters['lookup']['S' + state + 'GA' + str(gender) + str(age)])
                ii = ii + 1
    inf_tracker.update({'infections_state_age_gender': temp})

    # Update infectious state by age and risk
    temp = inf_tracker['infections_state_age_risk']
    ii = 0
    for age in [0, 1, 2, 3, 4]:
        for risk in [0, 1]:
            for state in ['S', 'E', 'I', 'R', 'T', 'V']:
                temp[t, ii] = np.sum(
                    pop_parameters['lookup']['S' + state + 'RA' + str(risk) + str(age)])
                ii = ii + 1
    inf_tracker.update({'infections_state_age_risk': temp})

    # Update infections by site for age and gender
    temp = inf_tracker['infections_site_age_gender']
    ii = 0
    for age in [0, 1, 2, 3, 4]:
        for gender in [0, 1]:
            site0 = pop_parameters['lookup']['I0GA1' + str(gender) + str(age)]
            site1 = pop_parameters['lookup']['I1GA1' + str(gender) + str(age)]
            site2 = pop_parameters['lookup']['I2GA1' + str(gender) + str(age)]
            temp[t, ii + 0] = np.sum(site0 & ~site1 & ~site2)
            temp[t, ii + 1] = np.sum(~site0 & site1 & ~site2)
            temp[t, ii + 2] = np.sum(~site0 & ~site1 & site2)
            temp[t, ii + 3] = np.sum(site0 & site1 & ~site2)
            temp[t, ii + 4] = np.sum(site0 & ~site1 & site2)
            temp[t, ii + 5] = np.sum(~site0 & site1 & site2)
            temp[t, ii + 6] = np.sum(site0 & site1 & site2)
            ii = ii + 7
    inf_tracker.update({'infections_site_age_gender': temp})

    # Update infections by site for age and risk
    temp = inf_tracker['infections_site_age_risk']
    ii = 0
    for age in [0, 1, 2, 3, 4]:
        for risk in [0, 1]:
            site0 = pop_parameters['lookup']['I0RA1' + str(risk) + str(age)]
            site1 = pop_parameters['lookup']['I1RA1' + str(risk) + str(age)]
            site2 = pop_parameters['lookup']['I2RA1' + str(risk) + str(age)]
            temp[t, ii + 0] = np.sum(site0 & ~site1 & ~site2)
            temp[t, ii + 1] = np.sum(~site0 & site1 & ~site2)
            temp[t, ii + 2] = np.sum(~site0 & ~site1 & site2)
            temp[t, ii + 3] = np.sum(site0 & site1 & ~site2)
            temp[t, ii + 4] = np.sum(site0 & ~site1 & site2)
            temp[t, ii + 5] = np.sum(~site0 & site1 & site2)
            temp[t, ii + 6] = np.sum(site0 & site1 & site2)
            ii = ii + 7
    inf_tracker.update({'infections_site_age_risk': temp})

    # Update infections by site and symptoms age and gender
    temp = inf_tracker['infections_site_age_gender_symptom']
    ii = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            site0 = pop_parameters['lookup']['I0GA1' + str(l2) + str(age)]
            site1 = pop_parameters['lookup']['I1GA1' + str(l2) + str(age)]
            site2 = pop_parameters['lookup']['I2GA1' + str(l2) + str(age)]
            co_inf = np.sum([site0, site1, site2], axis=0) > 1
            site0s = meta.site0_symptoms
            site1s = meta.site1_symptoms
            site2s = meta.site2_symptoms
            any_symp = site0s | site1s | site2s
            temp[t, ii + 0] = np.sum(site0 & site0s)
            temp[t, ii + 1] = np.sum(site0 & ~site0s)
            temp[t, ii + 2] = np.sum(site1 & site1s)
            temp[t, ii + 3] = np.sum(site1 & ~site1s)
            temp[t, ii + 4] = np.sum(site2 & site2s)
            temp[t, ii + 5] = np.sum(site2 & ~site2s)
            temp[t, ii + 6] = np.sum(co_inf & any_symp)
            temp[t, ii + 7] = np.sum(co_inf & ~any_symp)
            ii = ii + 8
    inf_tracker.update({'infections_site_age_gender_symptom': temp})

    # Update testing tracker
    temp = inf_tracker['track_treat']
    ii = 0
    for a in [0, 1, 2, 3, 4]:
        for g in [0, 1]:

            # Track people currently in the population
            who = meta.loc[(meta.age_group == a) & (meta.gender == g), :]
            last_365 = (t+t0_sim - who.test_time_last) <= 365
            in_count1 = np.sum(last_365 & (who.test_reason_last == int(1)))
            in_count2 = np.sum(last_365 & (who.test_reason_last == int(2)))
            in_count3 = np.sum(last_365 & (who.test_reason_last == int(3)))
            in_count_sum = np.sum(last_365)
            in_tot = len(who)

            # Update tracker
            n = max(1, in_tot)
            temp[t, ii] = (in_count_sum)/n
            temp[t, ii+1] = (in_count1)/n
            temp[t, ii+2] = (in_count2)/n
            temp[t, ii+3] = (in_count3)/n
            ii = ii + 4
    inf_tracker.update({'track_treat': temp})

    return inf_tracker


# %% GRAPH make_infection_graphs
def make_infection_graphs(tt, inf_tracker, pop_parameters, save_loc='graphs/output/'):

    # GRAPH SHOWING AGGREGATE INFECTIONS
    cmap = plt.get_cmap("tab10")
    nt = np.sum(inf_tracker['pop_size_age_risk'], axis=1)/100
    yt = inf_tracker['aggregate_infections']
    fig, ax = plt.subplots(2)
    # ax[0].plot(tt, yt["S"]/nt, label = "S")
    ax[0].plot(tt, yt["E"]/nt, label="E")
    ax[0].plot(tt, yt["I"]/nt, label="I")
    ax[0].plot(tt, yt["R"]/nt, label="R")
    ax[0].plot(tt, yt["T"]/nt, label="T")
    ax[0].plot(tt, yt["V"]/nt, label="V")
    ax[0].set_title("Proportion of Population in Each State")
    ax[0].legend()

    ax[1].plot(tt, yt["site0s"]/nt, label="Rectal")
    ax[1].plot(tt, yt["site1s"]/nt, label="Urogenital")
    ax[1].plot(tt, yt["site2s"]/nt, label="Pharyngeal")
    ax[1].plot(tt, yt["site0a"]/nt, linestyle='--',
               label="Rectal - asymptomatic", color=cmap(0))
    ax[1].plot(tt, yt["site1a"]/nt, linestyle='--',
               label="Urogenital - asymptomatic", color=cmap(1))
    ax[1].plot(tt, yt["site2a"]/nt, linestyle='--',
               label="Pharyngeal - asymptomatic", color=cmap(2))
    ax[1].set_title("Prevalence by Anatomical Site")
    ax[1].legend(ncol=2)

    # Save the graph
    fig.savefig(save_loc + 'infections_aggregate.png', dpi=200)
    plt.close(fig)

    # GRAPH SHOWING TOTAL INFECTIONS BY STATE AGE AND GENDER
    fig, ax = plt.subplots(5, 2)
    it = inf_tracker['infections_state_age_gender']
    ii = 0
    jj = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            denom = inf_tracker['pop_size_age_gender'][:, jj]/100
            # ax[age, l2].plot(tt, it[:, ii + 0]/denom, label = 'S')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 1], denom), label='E')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 2], denom), label='I')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 3], denom), label='R')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 4], denom), label='T')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 5], denom), label='V')
            ii = ii + 6
            jj = jj + 1

    # Labels
    ax[0, 0].set_ylabel('Age 16-19')
    ax[1, 0].set_ylabel('Age 20-24')
    ax[2, 0].set_ylabel('Age 25-29')
    ax[3, 0].set_ylabel('Age 30-34')
    ax[4, 0].set_ylabel('Age 35')

    fig.suptitle('Prevalence of Each Infectious State')
    ax[4, 0]. set_xlabel('Day')
    ax[4, 1]. set_xlabel('Day')
    ax[0, 0].set_title('Females')
    ax[0, 1].set_title('Males')
    ax[4, 1].legend(loc='upper center', ncol=6,
                    bbox_to_anchor=(-0.1, -0.2), fancybox=True)

    # Save graph
    fig.savefig(save_loc + 'infections_prev_by_age_and_gender.png', dpi=200)
    plt.close(fig)

    # GRAPH SHOWING TOTAL INFECTIONS BY STATE AGE AND RISK
    fig, ax = plt.subplots(5, 2)
    it = inf_tracker['infections_state_age_risk']
    ii = 0
    jj = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            denom = inf_tracker['pop_size_age_risk'][:, jj]/100
            # ax[age, l2].plot(tt, it[:, ii + 0]/denom, label = 'S')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 1], denom), label='E')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 2], denom), label='I')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 3], denom), label='R')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 4], denom), label='T')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 5], denom), label='V')
            ii = ii + 6
            jj = jj + 1

    # Labels
    ax[0, 0].set_ylabel('Age 16-19')
    ax[1, 0].set_ylabel('Age 20-24')
    ax[2, 0].set_ylabel('Age 25-29')
    ax[3, 0].set_ylabel('Age 30-34')
    ax[4, 0].set_ylabel('Age 35')

    fig.suptitle('Prevalence of Each Infectious State')
    ax[4, 0]. set_xlabel('Day')
    ax[4, 1]. set_xlabel('Day')
    ax[0, 0].set_title('Low-Risk')
    ax[0, 1].set_title('High-Risk')
    ax[4, 1].legend(loc='upper center', ncol=6,
                    bbox_to_anchor=(-0.1, -0.2), fancybox=True)

    # Save graph
    fig.savefig(save_loc + 'infections_prev_by_age_and_risk.png', dpi=200)
    plt.close(fig)

    # GRAPH SHOWING TOTAL INFECTIONS BY SITE AGE AND GENDER
    fig, ax = plt.subplots(5, 2)
    it = inf_tracker['infections_site_age_gender']
    ii = 0
    jj = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            denom = inf_tracker['pop_size_age_gender'][:, jj]/100
            # n_g = np.sum(it[:, ii:(ii+5)], axis = 1)
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 0], denom), label='Rectum')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 3], denom),
                             label='Rectum and Urethra')
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 6], denom), label='All Sites')
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 1], denom), label='Urethra')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 4], denom),
                             label='Rectum and Pharynx')
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 2], denom), label='Pharynx')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 5], denom),
                             label='Urethra and Pharynx')
            ii = ii + 7
            jj = jj + 1

    # Labels
    ax[0, 0].set_ylabel('Age 16-19')
    ax[1, 0].set_ylabel('Age 20-24')
    ax[2, 0].set_ylabel('Age 25-29')
    ax[3, 0].set_ylabel('Age 30-34')
    ax[4, 0].set_ylabel('Age 35')

    fig.suptitle('Proportion of People Infected at Each Anatomical Site')
    ax[4, 0]. set_xlabel('Day')
    ax[4, 1]. set_xlabel('Day')
    ax[0, 0].set_title('Female')
    ax[0, 1].set_title('Male')
    ax[4, 1].legend(loc='upper center', ncol=3,
                    bbox_to_anchor=(-0.1, -0.2), fancybox=True)

    # Save graph
    fig.savefig(save_loc + 'infections_prev_site_by_age_and_gender.png', dpi=200)
    plt.close(fig)

    # GRAPH SHOWING TOTAL INFECTIONS BY SITE AGE AND RISK
    fig, ax = plt.subplots(5, 2)
    it = inf_tracker['infections_site_age_risk']
    ii = 0
    jj = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            denom = inf_tracker['pop_size_age_risk'][:, jj]/100
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 0], denom), label='Rectum')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 3], denom),
                             label='Rectum and Urethra')
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 6], denom), label='All Sites')
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 1], denom), label='Urethra')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 4], denom),
                             label='Rectum and Pharynx')
            ax[age, l2].plot(tt, np.divide(
                it[:, ii + 2], denom), label='Pharynx')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 5], denom),
                             label='Urethra and Pharynx')
            ii = ii + 7
            jj = jj + 1

    # Labels
    ax[0, 0].set_ylabel('Age 16-19')
    ax[1, 0].set_ylabel('Age 20-24')
    ax[2, 0].set_ylabel('Age 25-29')
    ax[3, 0].set_ylabel('Age 30-34')
    ax[4, 0].set_ylabel('Age 35')

    fig.suptitle('Proportion of People Infected at Each Anatomical Site')
    ax[4, 0]. set_xlabel('Day')
    ax[4, 1]. set_xlabel('Day')
    ax[0, 0].set_title('Low-Risk')
    ax[0, 1].set_title('High-Risk')
    ax[4, 1].legend(loc='upper center', ncol=3,
                    bbox_to_anchor=(-0.1, -0.2), fancybox=True)

    # Save graph
    fig.savefig(save_loc + 'infections_prev_site_by_age_and_risk.png', dpi=200)
    plt.close(fig)

    # GRAPH SHOWING TOTAL INFECTIONS BY SITE AGE AND GENDER AND SYMPTOMS
    fig, ax = plt.subplots(5, 2)
    it = inf_tracker['infections_site_age_gender_symptom']
    ii = 0
    jj = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            denom = inf_tracker['pop_size_age_gender'][:, jj]/100
            ax[age, l2].plot(tt, np.divide(it[:, ii + 0], denom),
                             label='Rectum', color=cmap(0))
            ax[age, l2].plot(tt, np.divide(it[:, ii + 1], denom),
                             label='Rectum - Asymptomatic', color=cmap(0), linestyle='--')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 2], denom),
                             label='Pharynx', color=cmap(1))
            ax[age, l2].plot(tt, np.divide(it[:, ii + 3], denom),
                             label='Pharynx - Asymptomatic', color=cmap(1), linestyle='--')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 4], denom),
                             label='Urethra', color=cmap(2))
            ax[age, l2].plot(tt, np.divide(it[:, ii + 5], denom),
                             label='Urethra - Asymptomatic', color=cmap(2), linestyle='--')
            ax[age, l2].plot(tt, np.divide(it[:, ii + 6], denom),
                             label='Co-infection', color=cmap(3))
            ax[age, l2].plot(tt, np.divide(it[:, ii + 7], denom),
                             label='Co-infection - Asymptomatic', color=cmap(3), linestyle='--')
            ii = ii + 8
            jj = jj + 1

    # Labels
    ax[0, 0].set_ylabel('Age 16-19')
    ax[1, 0].set_ylabel('Age 20-24')
    ax[2, 0].set_ylabel('Age 25-29')
    ax[3, 0].set_ylabel('Age 30-34')
    ax[4, 0].set_ylabel('Age 35')

    fig.suptitle('Proportion of People Infected at Each Anatomical Site')
    ax[4, 0]. set_xlabel('Day')
    ax[4, 1]. set_xlabel('Day')
    ax[0, 0].set_title('Female')
    ax[0, 1].set_title('Male')
    ax[4, 1].legend(loc='upper center', ncol=4,
                    bbox_to_anchor=(-0.1, -0.2), fancybox=True)

    # Save graph
    fig.savefig(
        save_loc + 'infections_prev_site_by_age_and_risk_and_symptoms.png', dpi=200)
    plt.close(fig)

    # GRAPH SHOWING TESTING PREVALENCE BY AGE AND GENDER
    fig, ax = plt.subplots(5, 2)
    it = inf_tracker['track_treat']
    rates = pop_parameters['testing_rates']
    cmap = plt.get_cmap("tab10")
    ii = 0
    for age in [0, 1, 2, 3, 4]:
        for l2 in [0, 1]:
            target = rates.prob[(rates.age_group == age) & (rates.gender == l2)].values[0]
            ax[age, l2].plot(tt, len(tt) * [100 * target], linestyle='--', label='STRIVE Baseline Proportion', color=cmap(0))
            ax[age, l2].plot(tt, 100 * it[:, ii], label='Overall Testing', color=cmap(0))
            ax[age, l2].plot(tt, 100 * it[:, ii+1], label='Symptoms', color=cmap(1))
            ax[age, l2].plot(tt, 100 * it[:, ii+2], label='Background', color=cmap(2))
            ax[age, l2].plot(tt, 100 * it[:, ii+3], label='Contact Tracing', color=cmap(3))
            ii = ii + 4

    # Labels
    ax[0, 0].set_ylabel('Age 16-19')
    ax[1, 0].set_ylabel('Age 20-24')
    ax[2, 0].set_ylabel('Age 25-29')
    ax[3, 0].set_ylabel('Age 30-34')
    ax[4, 0].set_ylabel('Age 35')

    fig.suptitle(
        'Proportion of Sub-Populations Presenting for an STI Test\nBy Reason for Test')
    ax[4, 0]. set_xlabel('Day')
    ax[4, 1]. set_xlabel('Day')
    ax[0, 0].set_title('Female')
    ax[0, 1].set_title('Male')
    ax[4, 1].legend(loc='upper center', ncol=5,
                    bbox_to_anchor=(-0.1, -0.2), fancybox=True)

    # Save graph
    fig.savefig(save_loc + 'infections_testing_rates_by_age_gender_and_reason.png', dpi=200)
    plt.close(fig)
