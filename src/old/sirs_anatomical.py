# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 09:51:25 2020

@author: nicol
"""

import math
import random
import matplotlib.pyplot as plt
from operator import attrgetter
import time
import pandas as pd
import numpy as np

plt.rcParams['figure.figsize'] = (12.0, 7.5)

#%%###########################################################################
##  SETUP THE PERSON CLASS
##############################################################################


# MAKE THE PERSON CLASS
class person:
    
    # Define the probability of the person causing a new infection
    def transmission_probability(self):
        return 0.8
    
    # Define method for calculating the duration of infection
    def duration_infectious(self):
        return random.expovariate(1/7)
    
    # Define methods for calculating the duration of a partnership
    def duration_partnership(self):
        return random.expovariate(1/7)
    
    # Define probability of getting partnered up
    def partnership_probability(self):
        return 0.1
    
    # Define partnership biasing
    def partnership_bias(self, available):
        return random.choice(available.index)


# Function for setting when somebody decides to get treatment
def seek_treatment(meta, t):
    meta["treatment"] = \
        meta.apply(lambda row: row.site0/(1 + np.exp(0.5*(14-(t-row.site0_t0))))    \
                   + row.site1/(1 + np.exp(1.5*(4-(t-row.site1_t0))))               \
                   + row.site2/(1 + 10*np.exp(0.05*(60-(t-row.site2_t0))))          \
                   > row.treatment_threshold,
                   axis=1)
    return meta

# Define method for calculating the duration of immunity
def duration_removed(meta):
    return random.expovariate(1/14)

# Define transmission 
trans_sites = [[1, 2, 3, 1],
               [1, 2, 3, 1],
               [1, 2, 3, 1]]


#%%###########################################################################
##  INITILISE A LIST OF PEOPLE
##############################################################################

# Initilise population metadata
start_sim = 100
n_people = 100
people = []
meta = pd.DataFrame(index = range(0, n_people), 
                    columns = ["age", 
                               "sex", 
                               "location",
                               "state",
                               "site0",
                               "site1",
                               "site2",
                               "site0_t0",
                               "site1_t0",
                               "site2_t0",
                               "treatment_threshold",
                               "recovery_time",
                               "partner",
                               "partner_t1",
                               "methods"])


# Fill in metadata matrix
for i in range(0, n_people):
    
    # Basic demography
    meta.at[i, "age"] = random.uniform(0, 100)
    meta.at[i, "sex"] = int(random.random() > 0.5)
    meta.at[i, "location"] = 1
    meta.at[i, "methods"] = person()
    
    # Infection status
    if random.random() > 0.1:
        # Set as infected
        meta.at[i, "state"] = "I"
        
        # Set treatment threshold
        meta.at[i, "treatment_threshold"] = random.random()
        
        # Set duration of infection
        
        # Decide on one anatomical site to start infection
        u = np.cumsum(np.random.random(3))
        u = np.min(np.where(u/u[2] > np.random.random()))
        site = [0, 0, 0]
        site[u] = 1
        meta.at[i, "site0"] = site[0]
        meta.at[i, "site1"] = site[1]
        meta.at[i, "site2"] = site[2]
        meta.at[i, "site0_t0"] = start_sim
        meta.at[i, "site1_t0"] = start_sim
        meta.at[i, "site2_t0"] = start_sim
    else:
        meta.at[i, "state"] = "S"
        meta.at[i, "site0"] = 0
        meta.at[i, "site1"] = 0
        meta.at[i, "site2"] = 0
        meta.at[i, "site0_t0"] = 0
        meta.at[i, "site1_t0"] = 0
        meta.at[i, "site2_t0"] = 0
        
    # Partnerships
    if random.random() < meta.at[i, "methods"].partnership_probability():
        
        # Work out who is available
        single = meta[pd.isna(meta["partner"])]
        if len(single) > 0:
            
            # Decide on a partner
            ii = meta.at[i, "methods"].partnership_bias(single)
            meta.at[i, "partner"] = ii
            meta.at[ii, "partner"] = i
            
            # Decide on the duration of their partnership
            dur = meta.at[i, "methods"].duration_partnership()
            meta.at[i, "partner_t1"] = dur
            meta.at[ii, "partner_t1"] = dur
            
# Check on metadata
# plt.hist(meta.age)
# plt.hist(meta.sex)
# plt.hist(meta.location)
# plt.hist(meta.state)
# plt.hist(meta.state_t1)  
# plt.hist(meta.partner)      


#%%###########################################################################
##  SIMULATE AN INFECTION PROCESS
##############################################################################

# For plotting
T = range(0, 1000, 1)
xt = pd.DataFrame(index = range(0, len(T)), 
                  columns = ["S", "I", "R", "single", "partnered", "site0", "site1", "site2"])

# Iterate over each time point
for t in T:
    
    # Check for expired partnerships
    breakups = meta[meta["partner_t1"] < t]
    if len(breakups) > 0:
        i = breakups.index
        meta.loc[i, "partner"] = pd.NA
        meta.loc[i, "partner_t1"] = pd.NA
    
    # Make new partnerships
    single = meta[pd.isna(meta["partner"])]
    for i in single.index:
        if random.random() < meta.at[i, "methods"].partnership_probability():
        
            # Decide on a partner
            ii = meta.at[i, "methods"].partnership_bias(single)
            meta.at[i, "partner"] = ii
            meta.at[ii, "partner"] = i
            
            # Decide on the duration of their partnership
            dur = meta.at[i, "methods"].duration_partnership()
            meta.at[i, "partner_t1"] = t + dur
            meta.at[ii, "partner_t1"] = t + dur
            
            # Update the singles list
            single = meta[pd.isna(meta["partner"])]
            if len(single) < 2:
                break
    
    
    # Simulate infections after burn in
    if t > start_sim:
        
        # Remove infections
        meta = seek_treatment(meta, t)
        recoveries = meta[meta["treatment"] == True]
        if len(recoveries) > 1:
            for ii in recoveries.index:
                meta.at[ii, "state"] = "R"
                meta.at[ii, "recovery_time"] = t + duration_removed(meta)
                meta.at[ii, "site0"] = 0
                meta.at[ii, "site1"] = 0
                meta.at[ii, "site2"] = 0
                meta.at[ii, "site0_t0"] = 0
                meta.at[ii, "site1_t0"] = 0
                meta.at[ii, "site2_t0"] = 0
        
        # Remove immunity
        waning = meta[(meta["recovery_time"] < t) & (meta["state"] == "R")]
        if len(waning) > 1:
            for ii in waning.index:
                meta.at[ii, "state"] = "S"
                meta.at[ii, "recovery_time"] = pd.NA
            
        # Cause new infection
        possible_infection = meta[(meta["state"] == "I") & (pd.isna(meta["partner"]) == False)]
        for i in possible_infection.index:
            
            # Set treatment threshold
            meta.at[i, "treatment_threshold"] = random.random()
            
            # Choose a new site for infection
            dist = np.matmul([meta.at[i, "site0"], meta.at[i, "site1"], meta.at[i, "site2"]], trans_sites).cumsum()
            dist = dist/dist[3]
            u = random.random()
            site = np.min(np.where(dist >= u))
            
            
            # Do the infection
            if site < 3:
                ii = meta.at[i, "partner"]
                meta.at[ii, "state"] = "I"
                if meta.at[ii, "site" + str(site)] == 0:
                    meta.at[ii, "site" + str(site)] = 1
                    meta.at[ii, "site" + str(site) + "_t0"] = t
        
                
    #how long go for
    #prop sympo
    #treatment seeking
    
    
    # Updat plot vector
    xt.loc[t] = [sum(meta["state"]=="S"), 
                 sum(meta["state"]=="I"), 
                 sum(meta["state"]=="R"),
                 sum(pd.isna(meta["partner"])), 
                 n_people - sum(pd.isna(meta["partner"])),
                 sum(meta["site0"] == 1),
                 sum(meta["site1"] == 1),
                 sum(meta["site2"] == 1)]
    
    
# Make a plot
fig, ax = plt.subplots(3)
ax[0].plot(T, xt["S"], label = "S")
ax[0].plot(T, xt["I"], label = "I")
ax[0].plot(T, xt["R"], label = "R")
ax[0].set_title("Compartment Status")
ax[0].legend()

ax[1].plot(T, xt["single"], label = "Single")
ax[1].plot(T, xt["partnered"], label = "Partnered")
ax[1].set_title("Partnership Status")
ax[1].legend()

ax[2].plot(T, xt["site0"], label = "Rectal")
ax[2].plot(T, xt["site1"], label = "Urogenital")
ax[2].plot(T, xt["site2"], label = "Pharyngeal")
ax[2].set_title("Prevalence by Anatomical Site")
ax[2].legend()
                
            






#%%









