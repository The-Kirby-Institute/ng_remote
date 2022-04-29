# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:38:04 2021

@author: nicol

Treatment occurs roughly 7 days after symptoms
Both partners get treated
"""


#%% IMPORT REQUIRED MODULES
import numpy as np
import scipy.stats as sp


#%% FUN seek_treatment()
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


#######################
##  GETTING TREATED  ##
#######################
def seek_treatment(parameters, meta, partner_matrix, t):


    # Work out who is symptomatic
    symp0 = np.isfinite(meta["site0_t0"]) & np.isfinite(meta["site0_t1"])
    symp1 = np.isfinite(meta["site1_t0"]) & np.isfinite(meta["site1_t1"])
    symp2 = np.isfinite(meta["site2_t0"]) & np.isfinite(meta["site2_t1"])


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


    # Look up their current partners
    # part = np.asarray(np.where(partner_matrix[treat,:] == 1))
    # treat = np.append(treat, part[1,:]).tolist()


    # Make amendments to meta
    meta.loc[treat, "state"] = "R"
    meta.loc[treat, "recovery_time"] = t + np.random.gamma(parameters['infection'].immunity_mean[0]/parameters['infection'].immunity_var[0], parameters['infection'].immunity_var[0], len(treat))
    meta.loc[treat, "site0"] = 0
    meta.loc[treat, "site1"] = 0
    meta.loc[treat, "site2"] = 0
    meta.loc[treat, "site0_t0"] = float("inf")
    meta.loc[treat, "site1_t0"] = float("inf")
    meta.loc[treat, "site2_t0"] = float("inf")
    meta.loc[treat, "site0_t1"] = float("inf")
    meta.loc[treat, "site1_t1"] = float("inf")
    meta.loc[treat, "site2_t1"] = float("inf")


    # Return duration
    return meta









