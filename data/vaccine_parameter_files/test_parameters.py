# -*- coding: utf-8 -*-
"""
Created on Mon May 2
@author: Nicolas Rebuli

This is a top-level script containing parameters for vaccination simulations. These are provided in two forms: vaccine distribution, and vaccine effect. Both are described further below.

Vaccine distribution:
0 vaccines are given at the time of treatment.
1 everybody is vaccinated at the age of 16.
2 a proportion of each age group and gender is vaccinated.

Vaccine effect:
0 complete immunity for the duration of the vaccine
1 transmission probability is reduced
2 the probability of being asymptomatic is increased
3 the duration of infection is decreased
4 a combination of the above

Further parameters contributing to each of these is contained in the dictionary variable.

Notes:
* The vaccine is assumed to be 'effective' for a period of time T, where T comes from a Gamma distribution with mean and variance `duration_mean` and `duration_var`.
"""


# Setup parameters for vaccination
vax_parameters = {
    'duration_mean': 200*365,   # Mean of distribution of T
    'duration_var': 1,          # Variance of distribution of T
    'prop_effective': 0,        # Probability of vaccination working
    'effect': 0,
    'site0_trans_reduce': 0,   # Scaling of the transmission probability
    'site1_trans_reduce': 0,   # Scaling of the transmission probability
    'site2_trans_reduce': 0,   # Scaling of the transmission probability
    'site0_symp_reduce': 0.5,    # Scaling of the probability of symptoms
    'site1_symp_reduce': 0.5,    # Scaling of the probability of symptoms
    'site2_symp_reduce': 0.5,    # Scaling of the probability of symptoms
    'site0_duration_reduce': 0.5, # Scaling of the duration of inection
    'site1_duration_reduce': 0.5, # Scaling of the duration of inection
    'site2_duration_reduce': 0.5, # Scaling of the duration of inection
    'deployment': 0,
    'prop_vax_0': 0.5,  # Proportion of 16-19 year olds vaccinated
    'prop_vax_1': 0.4,  # Proportion of 20-24 year olds vaccinated
    'prop_vax_2': 0.1,  # Proportion of 25-29 year olds vaccinated
    'prop_vax_3': 0.0}  # Proportion of over 30 year olds vaccinated