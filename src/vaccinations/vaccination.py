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
# import matplotlib.pyplot as plt


# Import the baseline NG library for doing all the stuff that hasn't been changed
import src.infections.ng as ng


#%% FUN vaccine_reduced_trans_prob()
#
#
# Compute the site-specific transmission probabilities
#
#
def vaccine_reduced_trans_prob(inf_parameters, vax_parameters, meta, i, j):


    # Compute the usual transmission probability
    trans_prob = ng.transmission_probability(inf_parameters, vax_parameters, meta, i, j)


    # Construct the reduction vector
    vaccinated = 1.0 * (meta.vaccinated.iloc[i] | meta.vaccinated.iloc[j])
    vax_reduction = np.array([[vax_parameters['site0_trans_reduce'] * vaccinated + (1-vaccinated)], \
                              [vax_parameters['site1_trans_reduce'] * vaccinated + (1-vaccinated)], \
                              [vax_parameters['site2_trans_reduce'] * vaccinated + (1-vaccinated)]])
                
    
    # Multiply the transmission probability by the vaccine reduction
    out = vax_reduction * trans_prob


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
    modifyer = meta.vaccinated.iloc[j] * vax_parameters['site0_symp_reduce'] + (1-meta.vaccinated.iloc[j])
    prob = modifyer * p_baseline

    # Decide whether or not they will be symptomatic
    symptoms = np.random.random() < prob

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
    modifyer = meta.vaccinated.iloc[j] * vax_parameters['site1_symp_reduce'] + (1-meta.vaccinated.iloc[j])
    prob = modifyer * p_baseline

    # Decide whether or not they will be symptomatic
    symptoms = np.random.random() < prob

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
    modifyer = meta.vaccinated.iloc[j] * vax_parameters['site2_symp_reduce'] + (1-meta.vaccinated.iloc[j])
    prob = modifyer * p_baseline

    # Decide whether or not they will be symptomatic
    symptoms = np.random.random() < prob

    return symptoms


#%% VAR symptoms_vax
#
#
symptoms_vax = {'site0': symptoms_rectal,
                'site1': symptoms_pharynx,
                'site2': symptoms_urethra}


#%% VAR duration_vax
#
#
def duration_vax_site0(vax_parameters, meta, j):
    dur_reduce = meta.vaccinated.iloc[j] * vax_parameters['site0_duration_reduce'] + (1-meta.vaccinated.iloc[j])
    return dur_reduce

def duration_vax_site1(vax_parameters, meta, j):
    dur_reduce = meta.vaccinated.iloc[j] * vax_parameters['site1_duration_reduce'] + (1-meta.vaccinated.iloc[j])
    return dur_reduce

def duration_vax_site2(vax_parameters, meta, j):
    dur_reduce = meta.vaccinated.iloc[j] * vax_parameters['site2_duration_reduce'] + (1-meta.vaccinated.iloc[j])
    return dur_reduce

duration_vax = {'site0': duration_vax_site0,
                'site1': duration_vax_site1,
                'site2': duration_vax_site2}


#%% FUN set_function_for_updating_infections()
#
#
# Function to update the state of infections
#
#
def set_function_for_updating_infections(vax_parameters):


    # Where distribution occurrs during treatment
    if vax_parameters['deployment'] == 0:


        # Where the effect is total immunity
        if vax_parameters['effect'] == 0:


            # Situation where vaccinations given during treatment and makes people immune
            def update_infections(t, pop_parameters, inf_parameters, vax_parameters, meta, partner_matrix):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(pop_parameters, inf_parameters, meta, partner_matrix, t)
                meta = progress_state_of_infection(meta, t)
                # meta = progress_state_of_vaccination(meta, t)
                meta = ng.seek_treatment(pop_parameters, inf_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce transmission
        elif vax_parameters['effect'] == 1:


            # Situation where vaccines given during treatment and result is to
            # decrease transmission probability
            def update_infections(inf_parameters, vax_parameters, meta, partner_matrix, t):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(inf_parameters, meta, partner_matrix, t, ['S', 'E', 'I', 'V'], vaccine_reduced_trans_prob, vax_parameters)
                meta = progress_state_of_infection(meta, t)
                # meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is to reduce symptoms
        elif vax_parameters['effect'] == 2:


            # Situation where vaccines given during treatment and result in
            # people being more likely to be asymptomatic
            def update_infections(inf_parameters, vax_parameters, meta, partner_matrix, t):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(inf_parameters, meta, partner_matrix, t, ['S', 'E', 'I', 'V'], vaccine_reduced_trans_prob, vax_parameters, symptoms_vax)
                meta = progress_state_of_infection(meta, t)
                # meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


        # Where the effect is the reduce the duration of infection
        elif vax_parameters['effect'] == 3:


            # Situation where vaccines given during treatment and result in
            # people having a shorter duration of infection
            def update_infections(inf_parameters, vax_parameters, meta, partner_matrix, t):
                meta = progress_state_of_vaccination(meta, t)
                meta = ng.new_infections(inf_parameters, meta, partner_matrix, t, ['S', 'E', 'I', 'V'], vaccine_reduced_trans_prob, duration_mod = duration_vax)
                meta = progress_state_of_infection(meta, t)
                # meta = progress_state_of_vaccination(meta, t)
                meta = seek_treatment(inf_parameters, vax_parameters, meta, partner_matrix, t)
                return meta


    return update_infections



#%% FUN vaccination_duration()
#
#
# Drawn from Gamma distribution
#
#
def vaccination_duration(vax_parameters, n):


    # Sample duration
    duration = np.random.gamma(vax_parameters['duration_mean'] / \
                               vax_parameters['duration_var'], \
                               vax_parameters['duration_var'], n)


    return duration


#%% FUN give_vaccine()
#
#
# What happens when somebody is vaccinated while being treated for an
# active NG infection
#
#
def give_vaccine(inf_parameters, vax_parameters, meta, t, treat):


    # Decide which vaccinations are effective
    uu = np.random.random(len(treat))
    uu = uu < vax_parameters['prop_effective']


    # Update their vaccination data
    # meta.loc[treat[uu], 'state'] = 'V'
    meta.loc[treat[uu], 'vaccinated'] = True
    meta.loc[treat[uu], 'vaccination_t0'] = t


    # Set the time until the vaccine wears off
    meta.loc[treat[uu], 'vaccination_t1'] = t + vaccination_duration(vax_parameters, len(treat[uu]))


    # Pass back meta
    return meta


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
    meta.loc[new, 'vaccination_t0'] = float('inf')
    meta.loc[new, 'vaccination_t1'] = float('inf')
    meta.loc[new, 'booster_t0'] = float('inf')


    # Remove vaccine-conferred immunity
    waned_immunity = meta["vaccination_t1"] < t
    meta.loc[waned_immunity, "state"] = "S"
    meta.loc[waned_immunity, "vaccinated"] = False


    return meta


#%% NG MOD seek_treatment()
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
def seek_treatment(parameters, vax_parameters, meta, partner_matrix, t):


    # Work out who is symptomatic
    symp0 = meta.site0_symptoms == True
    symp1 = meta.site1_symptoms == True
    symp2 = meta.site2_symptoms == True


    # Work out their duration of infectiousness
    dur0 = t - meta.loc[symp0, "site0_t0"]
    dur1 = t - meta.loc[symp1, "site1_t0"]
    dur2 = t - meta.loc[symp2, "site2_t0"]


    # Work out whose infectious period has exceeded their tolerance
    treat0 = sp.gamma.cdf(dur0, parameters['infection'].treatment_mean[0]/parameters['infection'].treatment_var[0], parameters['infection'].treatment_var[0]) >= meta.loc[symp0, "treatment_threshold"]
    treat1 = sp.gamma.cdf(dur1, parameters['infection'].treatment_mean[0]/parameters['infection'].treatment_var[0], parameters['infection'].treatment_var[0]) >= meta.loc[symp1, "treatment_threshold"]
    treat2 = sp.gamma.cdf(dur2, parameters['infection'].treatment_mean[0]/parameters['infection'].treatment_var[0], parameters['infection'].treatment_var[0]) >= meta.loc[symp2, "treatment_threshold"]


    # Pull out their identifiers
    treat0 = treat0.index[treat0]
    treat1 = treat1.index[treat1]
    treat2 = treat2.index[treat2]


    # Combine these into one big list
    treat = np.append(treat0, treat1)
    treat = np.append(treat, treat2)
    treat = np.unique(treat)


    # Have their current long-term partners get treated as well in X% of partnerships
    # part = np.asarray(np.where(partner_matrix[treat,:] == 1)) <- this was commented out already
    part = meta.partner.iloc[treat]
    part = part[part > -1]
    part = part[np.random.random(len(part)) < parameters['infection'].prob_partner_treatment[0]]
    treat = np.append(treat, part)
    

    # Make amendments to meta
    meta.loc[treat, "state"] = "T"
    meta.loc[treat, "recovery_time"] = t + ng.duration_treatment_immunity(parameters, treat)
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


    # Implement vaccinations
    if len(treat) > 0:
        meta = give_vaccine(parameters, vax_parameters, meta, t, treat)


    # Return duration
    return meta


#%% NG MOD progress_state_of_infection()
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
    natural_recovery = ((recovery0 | recovery1) | recovery2) & (meta.site0==0) & (meta.site1==0) & (meta.site2==0)
    meta.at[natural_recovery, 'state'] = 'S'


    # Identify those where treatment has worn off
    waned_immunity = meta["recovery_time"] < t
    meta.loc[waned_immunity, "recovery_time"] = float("inf")


    # Work out if they should progress to vacinated or susceptible
    meta.loc[waned_immunity & (meta.vaccinated == False), "state"] = 'S'
    meta.loc[waned_immunity & (meta.vaccinated == True), "state"] = 'V'


    # Progress the state of vaccinated people
    meta = progress_state_of_vaccination(meta, t)

    return meta



