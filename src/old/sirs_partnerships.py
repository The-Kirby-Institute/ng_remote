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

    # Define attributes of the person
    #def __init__(self):
    
    # Set what to print when somebody asks about this object
    #def __str__(self):
        #return f"id: {self.id}\n age: {self.age}\n sex: {self.sex}\n location: {self.location}\n compartment:   {self.compartment}\n"

    # Instance method (a function which can be run)
    #def description(self):
        #return f"Person {self.id} is a {self.age} year old {self.sex} located in {self.location}"
    
    # Define the probability of the person causing a new infection
    def transmission_probability(self):
        return 0.8
    
    # Define method for calculating the duration of infection
    def duration_infectious(self):
        return random.expovariate(1/7)
    
    # Define method for calculating the duration of immunity
    def duration_removed(self):
        return random.expovariate(1/14)
    
    # Define methods for calculating the duration of a partnership
    def duration_partnership(self):
        return random.expovariate(1/7)
    
    # Define probability of getting partnered up
    def partnership_probability(self):
        return 0.1
    
    # Define partnership biasing
    def partnership_bias(self, available):
        return random.choice(available.index)
        


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
                               "state_t0",
                               "state_dt",
                               "state_t1",
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
        meta.at[i, "state"] = "S"
    else:
        meta.at[i, "state"] = "I"
        meta.at[i, "state_dt"] = meta.at[i, "methods"].duration_infectious()
        meta.at[i, "state_t0"] = start_sim
        meta.at[i, "state_t1"] = meta.at[i, "state_t0"] + meta.at[i, "state_dt"]
        
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
                  columns = ["S", "I", "R", "single", "partnered"])

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
        recoveries = meta[(meta["state_t1"] < t) & (meta["state"] == "I")]
        if len(recoveries) > 1:
            for ii in recoveries.index:
                meta.loc[ii, "state"] = "R"
                meta.loc[ii, "state_t1"] = t + meta.at[ii, "methods"].duration_removed()
        
        # Remove immunity
        waning = meta[(meta["state_t1"] < t) & (meta["state"] == "R")]
        if len(waning) > 1:
            for ii in waning.index:
                meta.loc[ii, "state"] = "S"
                meta.loc[ii, "state_t1"] = pd.NA
            
        # Cause new infection
        possible_infection = meta[(meta["state"] == "I") & (pd.isna(meta["partner"]) == False)]
        for i in possible_infection.index:
            if random.random() < meta.at[i, "methods"].transmission_probability():
                ii = meta.at[i, "partner"]
                meta.at[ii, "state"] = "I"
                meta.at[ii, "state_t1"] = t + meta.at[ii, "methods"].duration_infectious()
        
                
    #how long go for
    #prop sympo
    #treatment seeking
    
    
    # Updat plot vector
    xt.loc[t] = [sum(meta["state"]=="S"), sum(meta["state"]=="I"), sum(meta["state"]=="R"),
                 sum(pd.isna(meta["partner"])), n_people - sum(pd.isna(meta["partner"]))]
    
    
# Make a plot
fig, ax = plt.subplots(2)
ax[0].plot(T, xt["S"], label = "S")
ax[0].plot(T, xt["I"], label = "I")
ax[0].plot(T, xt["R"], label = "R")
ax[0].set_title("Compartment Status")
ax[0].legend()

ax[1].plot(T, xt["single"], label = "Single")
ax[1].plot(T, xt["partnered"], label = "Partnered")
ax[1].set_title("Partnership Status")
ax[1].legend()
                
            






#%%









