# -*- coding: utf-8 -*-
"""
Created on Mon May 2
@author: Nicolas Rebuli

This is a top-level script containing parameters for vaccination simulations. 
These are provided in two forms: vaccine distribution, and vaccine effect. 
Both are described further below.

Vaccine effect:
0 infection probability is reduced
1 transmission probability is reduced
2 the probability of being detected is reduced (asymptomatic cases increased)
3 the duration of infection is decreased
4 a combination of the above

Vaccine distribution:
0 vaccines are given at the time of treatment.
1 everybody is vaccinated at the age of 16.
2 a proportion of each age group and gender is vaccinated.
3 a combination of the above

Further parameters contributing to each of these is contained in the dictionary variable.

Notes:
* The vaccine is assumed to be 'effective' for a period of time T, where T comes from a Gamma distribution with mean and variance `duration_mean` and `duration_var`.


THIS ONE OF THE SCENARIO SCRIPTS
vaccine reduces the probability of infection and transmission by 100%


"""


# Setup parameters for vaccination
vax_parameters = {
    
    # TOP-LEVEL SIMULATION PARAMETERS
    'vax_scenario': 'vax_scenario_0_baseline',  # Give the parameter set a name
    'pop_scenario': int(2),             # Which population size to use? 1 (N=200), 2 (N=750), 3 (N=2200)
    'n_sims': 1,                       # Set the desired number of simulations to gather
    'n_years': 1,                      # Set the number of years you want the simulation to run for
    
    # CONTROLS THE DURATION THE VACCINE IS EFFECTIVE
    'duration_mean': 2*365,           # Mean of distribution
    'duration_var': 30,                  # Variance of distribution
    
    # VACCINE EFFECT PARAMETERS
    'effect': 4,                        # LEAVE THIS
    'site0_inf_mult': 1,                # Scaling of the susceptibility (rectal)
    'site1_inf_mult': 1,                # Scaling of the susceptibility (urethral)
    'site2_inf_mult': 1,                # Scaling of the susceptibility (pharyngeal)
    'site0_trans_mult': 1,              # Scaling of the transmission probability (rectal)
    'site1_trans_mult': 1,              # Scaling of the transmission probability (urethral)
    'site2_trans_mult': 1,              # Scaling of the transmission probability (pharyngeal)
    'site0_symp_mult': 1,             # Scaling of the probability of testing
    'site1_symp_mult': 1,             # Scaling of the probability of testing
    'site2_symp_mult': 1,             # Scaling of the probability of testing
    'site0_duration_mult': 1,         # Scaling of the duration of inection
    'site1_duration_mult': 1,         # Scaling of the duration of inection
    'site2_duration_mult': 1,         # Scaling of the duration of inection
    
    # VACCINE DEPLOYMENT PARAMETERS
    'deployment': 3,                    # LEAVE THIS
    'booster_delay': float('inf'),      # Baseline assumes no boosted vaccination
    'p_treat_to_vax': 0,                # Proportion of people treated who get vaccinated
    'p_vax_16': 0.75,                   # Proportion of people vaccinated at sexual debut
    'prop_vax_0': 0.5,                  # Proportion of 16-19 year olds vaccinated
    'prop_vax_1': 0.5,                  # Proportion of 20-24 year olds vaccinated
    'prop_vax_2': 0.5,                  # Proportion of 25-29 year olds vaccinated
    'prop_vax_3': 0.5,                  # Proportion of over 30 year olds vaccinated
    'prop_vax_4': 0.5,                  # Proportion of over 30 year olds vaccinated
    
    # Old redundant parameters
    'prop_effective': 1,                # Changes the proportion of vaccines which work
                                        # This is useful from a testing perspective sometimes
                                        # so I have left it in.
    }




# check the tesing rates by age group worked - done!
# make better graphs and outputs - put in trackers just need to make the graphs
# formalise run scripts for scenarios