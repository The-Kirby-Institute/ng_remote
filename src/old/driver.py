# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:19:43 2021

@author: nicol
"""



#%% SETUP Load Modules


# Load libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
from tqdm import tqdm


# Load modules for simulation script
import src.demographic.generate_population as pop
import src.partners.partners as prt
import src.infections.ng as ng
import src.treatment.simple as trt


# Load in graphing modules
import src.partners.summary_stats as pstat
import src.infections.summary_stats as istat


#%% SETUP Load Data

# Read in simulation parameters
param = pd.read_csv("data/param.csv")


# Generate population
meta = pop.generate_population()


# Initilise partnership network
partner_matrix = pop.initilise_partner_matrix()
partner_expire = pop.initilise_partner_duration()


#%% SETUP Initilise Tracking Arrays


# Indicators of what to track
track_partnership_stats = False
track_partnership_rates = True
track_infection_rates = True


# Set time to start counting from
n_days = param.partner_burn_in[0] + param.simulation_length[0]


# Tracking the number of people in each partnership type
p0ht = np.zeros((n_days, 5))
p0lt = np.zeros((n_days, 5))
p1ht = np.zeros((n_days, 5))
p1lt = np.zeros((n_days, 5))
p2ht = np.zeros((n_days, 5))
p2lt = np.zeros((n_days, 5))
p3ht = np.zeros((n_days, 5))
p3lt = np.zeros((n_days, 5))


# Tracking the cumulative number of partners
xt = np.zeros((n_days, 4))
g0t = np.zeros((n_days, 4))
g1t = np.zeros((n_days, 4))
g2t = np.zeros((n_days, 4))
g3t = np.zeros((n_days, 4))


# Tracking partnership durations
d0t = []
d1t = []
d2t = []
d3t = []


# Tracking the infection process
yt = pd.DataFrame(columns = ["S", "I", "R", "single", "partnered", "site0", "site1", "site2", "E"])


# Tracking the infection process by age group and risk group
i0lt = np.zeros((n_days, 6))
i0ht = np.zeros((n_days, 6))
i1lt = np.zeros((n_days, 6))
i1ht = np.zeros((n_days, 6))
i2lt = np.zeros((n_days, 6))
i2ht = np.zeros((n_days, 6))
i3lt = np.zeros((n_days, 6))
i3ht = np.zeros((n_days, 6))


#%% RUN Simulation


# Iterate over all time points
for t in tqdm(range(0, n_days)):


    # Create a new relationships
    meta, partner_matrix, partner_expire, d0ti, d1ti, d2ti, d3ti = prt.new_partnership(meta, partner_matrix, partner_expire, t)


    # Update partnership duration tracking
    if track_partnership_stats:
        d0t.append(d0ti)
        d1t.append(d1ti)
        d2t.append(d2ti)
        d3t.append(d3ti)


    # If tracking the infection process
    if track_infection_rates:


        # Update infection status after
        if t > param.partner_burn_in[0]:


            # Seed new infections
            meta = ng.new_infections(meta, partner_matrix, t)


            # Progress the state of infections
            meta = ng.progress_state_of_infection(meta, t)


            # Allow people to seek treatment
            meta = trt.seek_treatment(meta, partner_matrix, t)


        # Update summary stats on prevalence
        if track_infection_rates:
            yt = istat.update_prevalence(meta, t, yt)
            i0lt, i0ht, i1lt, i1ht, i2lt, i2ht, i3lt, i3ht = istat.update_infections_by_group(meta, t, i0lt, i0ht, i1lt, i1ht, i2lt, i2ht, i3lt, i3ht)


    # Identify and remove expired relationships
    meta, partner_matrix, partner_expire = prt.old_partnerships(meta, partner_matrix, partner_expire, t)


    # Update summary stats on the number of annual partners
    if track_partnership_rates:
        xt, g0t, g1t, g2t, g3t = pstat.update_cumulative_partners(meta, partner_matrix, partner_expire, t, pop.n, xt, g0t, g1t, g2t, g3t)
        p0ht, p0lt, p1ht, p1lt, p2ht, p2lt, p3ht, p3lt = pstat.update_partnership_types(meta, partner_matrix, t, p0ht, p0lt, p1ht, p1lt, p2ht, p2lt, p3ht, p3lt)
        if t == param.partner_burn_in[0]: meta.counter = np.sum(partner_matrix, axis=1)


#%% GRAPH Prevalence Age Group


# If tracking prevalence
if track_infection_rates:


    # Setup graph
    t = range(0, n_days)
    fig, ax = plt.subplots(4, 2)


    # Graph [0,0]
    n_g = sum((meta.age_group == 0) & (meta.risk == 1))
    ax[0, 0].plot(t, i0ht[:,0]/n_g)
    ax[0, 0].plot(t, i0ht[:,1]/n_g)
    ax[0, 0].plot(t, i0ht[:,2]/n_g)
    ax[0, 0].plot(t, i0ht[:,3]/n_g)
    ax[0, 0].plot(t, i0ht[:,4]/n_g)
    ax[0, 0].plot(t, i0ht[:,5]/n_g)
    #ax[0, 0].plot(t, i0ht[:,6]/n_g)
    ax[0, 0].set_title("High-risk")
    ax[0, 0].set_ylabel("Age 16-19")
    ax[0, 0].set_xlim([param.partner_burn_in[0], n_days])
    ax[0, 0].set_ylim([0, 1])


    # Graph [0,1]
    n_g = sum((meta.age_group == 0) & (meta.risk == 0))
    ax[0, 1].plot(t, i0lt[:,0]/n_g)
    ax[0, 1].plot(t, i0lt[:,1]/n_g)
    ax[0, 1].plot(t, i0lt[:,2]/n_g)
    ax[0, 1].plot(t, i0lt[:,3]/n_g)
    ax[0, 1].plot(t, i0lt[:,4]/n_g)
    ax[0, 1].plot(t, i0lt[:,5]/n_g)
    #ax[0, 1].plot(t, i0lt[:,6]/n_g)
    ax[0, 1].set_title("Low-risk")
    ax[0, 1].set_xlim([param.partner_burn_in[0], n_days])
    ax[0, 1].set_ylim([0, 1])


    # Graph [1,0]
    n_g = sum((meta.age_group == 1) & (meta.risk == 1))
    ax[1, 0].plot(t, i1ht[:,0]/n_g)
    ax[1, 0].plot(t, i1ht[:,1]/n_g)
    ax[1, 0].plot(t, i1ht[:,2]/n_g)
    ax[1, 0].plot(t, i1ht[:,3]/n_g)
    ax[1, 0].plot(t, i1ht[:,4]/n_g)
    ax[1, 0].plot(t, i1ht[:,5]/n_g)
    #ax[1, 0].plot(t, i1ht[:,6]/n_g)
    ax[1, 0].set_ylabel("Age 20-24")
    ax[1, 0].set_xlim([param.partner_burn_in[0], n_days])
    ax[1, 0].set_ylim([0, 1])


    # Graph [1,1]
    n_g = sum((meta.age_group == 1) & (meta.risk == 0))
    ax[1, 1].plot(t, i1lt[:,0]/n_g)
    ax[1, 1].plot(t, i1lt[:,1]/n_g)
    ax[1, 1].plot(t, i1lt[:,2]/n_g)
    ax[1, 1].plot(t, i1lt[:,3]/n_g)
    ax[1, 1].plot(t, i1lt[:,4]/n_g)
    ax[1, 1].plot(t, i1lt[:,5]/n_g)
    #ax[1, 1].plot(t, i1lt[:,6]/n_g)
    ax[1, 1].set_xlim([param.partner_burn_in[0], n_days])
    ax[1, 1].set_ylim([0, 1])


    # Graph [2,0]
    n_g = sum((meta.age_group == 2) & (meta.risk == 1))
    ax[2, 0].plot(t, i2ht[:,0]/n_g)
    ax[2, 0].plot(t, i2ht[:,1]/n_g)
    ax[2, 0].plot(t, i2ht[:,2]/n_g)
    ax[2, 0].plot(t, i2ht[:,3]/n_g)
    ax[2, 0].plot(t, i2ht[:,4]/n_g)
    ax[2, 0].plot(t, i2ht[:,5]/n_g)
    #ax[2, 0].plot(t, i2ht[:,6]/n_g)
    ax[2, 0].set_ylabel("Age 24-29")
    ax[2, 0].set_xlim([param.partner_burn_in[0], n_days])
    ax[2, 0].set_ylim([0, 1])


    # Graph [2,1]
    n_g = sum((meta.age_group == 2) & (meta.risk == 0))
    ax[2, 1].plot(t, i2lt[:,0]/n_g)
    ax[2, 1].plot(t, i2lt[:,1]/n_g)
    ax[2, 1].plot(t, i2lt[:,2]/n_g)
    ax[2, 1].plot(t, i2lt[:,3]/n_g)
    ax[2, 1].plot(t, i2lt[:,4]/n_g)
    ax[2, 1].plot(t, i2lt[:,5]/n_g)
    #ax[2, 1].plot(t, i2lt[:,6]/n_g)
    ax[2, 1].set_xlim([param.partner_burn_in[0], n_days])
    ax[2, 1].set_ylim([0, 1])


    # Graph [3,0]
    n_g = sum((meta.age_group == 3) & (meta.risk == 1))
    ax[3, 0].plot(t, i3ht[:,0]/n_g, label = "Rectal")
    ax[3, 0].plot(t, i3ht[:,1]/n_g, label = "Urethral")
    ax[3, 0].plot(t, i3ht[:,2]/n_g, label = "Pharyngeal")
    ax[3, 0].plot(t, i3ht[:,3]/n_g, label = "Removed")
    ax[3, 0].plot(t, i3ht[:,4]/n_g, label = "Susceptible")
    ax[3, 0].plot(t, i3ht[:,5]/n_g, label = "Exposed")
    #ax[3, 0].plot(t, i3ht[:,6]/n_g, label = "Immune")
    ax[3, 0].legend(loc = "lower left", bbox_to_anchor=(0.18,-0.05), bbox_transform=fig.transFigure, ncol=3)
    ax[3, 0].set_xlim([param.partner_burn_in[0], n_days])
    ax[3, 0].set_ylabel("Age 29 or Over")
    ax[3, 0].set_xlabel("Time (days)")
    ax[3, 0].set_ylim([0, 1])


    # Graph [3,1]
    n_g = sum((meta.age_group == 3) & (meta.risk == 0))
    ax[3, 1].plot(t, i3lt[:,0]/n_g)
    ax[3, 1].plot(t, i3lt[:,1]/n_g)
    ax[3, 1].plot(t, i3lt[:,2]/n_g)
    ax[3, 1].plot(t, i3lt[:,3]/n_g)
    ax[3, 1].plot(t, i3lt[:,4]/n_g)
    ax[3, 1].plot(t, i3lt[:,5]/n_g)
    #ax[3, 1].plot(t, i3lt[:,6]/n_g)
    ax[3, 1].set_xlim([param.partner_burn_in[0], n_days])
    ax[3, 1].set_ylim([0, 1])
    ax[3, 1].set_xlabel("Time (days)")



#%% GRAPH Prevalence Overall


# If tracking prevalence
if track_infection_rates:


    # Setup graph
    fig, ax = plt.subplots(3)
    t = range(0, n_days)


    # Plot of S, I, R occupancy
    ax[0].plot(t, yt["S"], label = "S")
    ax[0].plot(t, yt["I"], label = "I")
    ax[0].plot(t, yt["R"], label = "R")
    ax[0].set_title("Compartment Status")
    ax[0].legend()


    # Graph of the number partnered
    ax[1].plot(t, yt["single"], label = "Single")
    ax[1].plot(t, yt["partnered"], label = "Partnered")
    ax[1].set_title("Partnership Status")
    ax[1].legend()


    # Graph of prevalence by anatomical site of infection
    ax[2].plot(t, yt["site0"], label = "Rectal")
    ax[2].plot(t, yt["site1"], label = "Urogenital")
    ax[2].plot(t, yt["site2"], label = "Pharyngeal")
    ax[2].set_title("Prevalence by Anatomical Site")
    ax[2].legend()


#%% GRAPH Calibration Timeseries


# If tracking partnership rates
if track_partnership_rates:


    # Read in data to calibrate to
    partner_prob_raw = pd.read_csv('data/calibration_partnership_rates.csv')


    # Setup graph
    fig, ax = plt.subplots(4)


    # Plot long vs short term partnerships over time
    # ax[0].plot(range(0, n_days), xt[:,0], label = "single")
    # ax[0].plot(range(0, n_days), xt[:,1], label = "long term")
    # ax[0].plot(range(0, n_days), xt[:,2], label = "short term")
    # ax[0].plot(range(0, n_days), xt[:,3], label = "time left")
    # ax[0].set_ylim([0, n])
    # ax[0].legend(loc = "upper right")


    # Plot the number of people in each partnership group over time
    ax[0].plot(range(0, n_days), g0t[:,0], label = "0", color = "blue")
    ax[0].plot(range(0, n_days), g0t[:,1], label = "1", color = "red")
    ax[0].plot(range(0, n_days), g0t[:,2], label = "2-4", color = "green")
    ax[0].plot(range(0, n_days), g0t[:,3], label = "5 more", color = "purple")
    weight = sum(g0t[0,:])
    ax[0].axhline(weight*partner_prob_raw.iloc[0,0], color = "blue", linestyle = "--")
    ax[0].axhline(weight*partner_prob_raw.iloc[1,0], color = "red", linestyle = "--")
    ax[0].axhline(weight*partner_prob_raw.iloc[2,0], color= "green", linestyle = "--")
    ax[0].axhline(weight*partner_prob_raw.iloc[3,0], color = "purple", linestyle = "--")
    ax[0].set_xlim([param.partner_burn_in[0], n_days])
    ax[0].legend(loc = "upper right")
    ax[0].set_title("Comparison of 12-month partner count to GOANNA data")


    # Plot the number of people in each partnership group over time
    ax[1].plot(range(0, n_days), g1t[:,0], label = "0", color = "blue")
    ax[1].plot(range(0, n_days), g1t[:,1], label = "1", color = "red")
    ax[1].plot(range(0, n_days), g1t[:,2], label = "2-4", color = "green")
    ax[1].plot(range(0, n_days), g1t[:,3], label = "5 more", color = "purple")
    weight = sum(g1t[0,:])
    ax[1].axhline(weight*partner_prob_raw.iloc[0,1], color = "blue", linestyle = "--")
    ax[1].axhline(weight*partner_prob_raw.iloc[1,1], color = "red", linestyle = "--")
    ax[1].axhline(weight*partner_prob_raw.iloc[2,1], color= "green", linestyle = "--")
    ax[1].axhline(weight*partner_prob_raw.iloc[3,1], color = "purple", linestyle = "--")
    ax[1].set_xlim([param.partner_burn_in[0], n_days])


    # Plot the number of people in each partnership group over time
    ax[2].plot(range(0, n_days), g2t[:,0], label = "0", color = "blue")
    ax[2].plot(range(0, n_days), g2t[:,1], label = "1", color = "red")
    ax[2].plot(range(0, n_days), g2t[:,2], label = "2-4", color = "green")
    ax[2].plot(range(0, n_days), g2t[:,3], label = "5 more", color = "purple")
    weight = sum(g2t[0,:])
    ax[2].axhline(weight*partner_prob_raw.iloc[0,2], color = "blue", linestyle = "--")
    ax[2].axhline(weight*partner_prob_raw.iloc[1,2], color = "red", linestyle = "--")
    ax[2].axhline(weight*partner_prob_raw.iloc[2,2], color= "green", linestyle = "--")
    ax[2].axhline(weight*partner_prob_raw.iloc[3,2], color = "purple", linestyle = "--")
    ax[2].set_xlim([param.partner_burn_in[0], n_days])


    # Plot the number of people in each partnership group over time
    ax[3].plot(range(0, n_days), g3t[:,0], label = "0", color = "blue")
    ax[3].plot(range(0, n_days), g3t[:,1], label = "1", color = "red")
    ax[3].plot(range(0, n_days), g3t[:,2], label = "2-4", color = "green")
    ax[3].plot(range(0, n_days), g3t[:,3], label = "5 more", color = "purple")
    weight = sum(g3t[0,:])
    ax[3].axhline(weight*partner_prob_raw.iloc[0,3], color = "blue", linestyle = "--")
    ax[3].axhline(weight*partner_prob_raw.iloc[1,3], color = "red", linestyle = "--")
    ax[3].axhline(weight*partner_prob_raw.iloc[2,3], color= "green", linestyle = "--")
    ax[3].axhline(weight*partner_prob_raw.iloc[3,3], color = "purple", linestyle = "--")
    ax[3].set_xlim([param.partner_burn_in[0], n_days])


#%% GRAPH Calibration Bar Graph


# If tracking partnership rates
if track_partnership_rates:


    # Setup graph
    fig, ax = plt.subplots(4)
    fig.suptitle("Comparison of Simulated Partnership Rates to Observed Data")
    labs = ["0", "1", "2-4", "5 or more"]
    sim = g0t[len(g0t)-1,:]/sum(g0t[0,:])
    goanna = partner_prob_raw.iloc[0,:]
    x = np.arange(len(labs))
    width = 0.35


    # Graph of partnership rates for 16-20 age group
    ax[0].bar(x-width/2, g0t[len(g0t)-1,:]/sum(g0t[0,:]), width, label = "Simulated")
    ax[0].bar(x+width/2, partner_prob_raw.iloc[:,0], width, label = "GOANNA")
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(labs)
    ax[0].set_title("Age 16-20")
    ax[0].set_ylabel("Proportion")
    ax[0].legend()


    # Graph of partnership rates for 21-24 age group
    ax[1].bar(x-width/2, g1t[len(g1t)-1,:]/sum(g1t[0,:]), width, label = "Simulated")
    ax[1].bar(x+width/2, partner_prob_raw.iloc[:,1], width, label = "GOANNA")
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(labs)
    ax[1].set_title("Age 21-24")
    ax[1].set_ylabel("Proportion")


    # Graph of partnership rates for 25-29 age group
    ax[2].bar(x-width/2, g2t[len(g2t)-1,:]/sum(g2t[0,:]), width, label = "Simulated")
    ax[2].bar(x+width/2, partner_prob_raw.iloc[:,2], width, label = "GOANNA")
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(labs)
    ax[2].set_title("Age 25-29")
    ax[2].set_ylabel("Proportion")


    # Graph of partnership rates for over 30 age group
    ax[3].bar(x-width/2, g3t[len(g3t)-1,:]/sum(g3t[0,:]), width, label = "Simulated")
    ax[3].bar(x+width/2, partner_prob_raw.iloc[:,3], width, label = "GOANNA")
    ax[3].set_xticks(x)
    ax[3].set_xticklabels(labs)
    ax[3].set_xlabel("Number of sexual partners in last 12 months")
    ax[3].set_title("Age 30 and Over")
    ax[3].set_ylabel("Proportion")


#%% GRAPH Partnership Statistics


# If tracking partnership stats
if track_partnership_stats:


    # Unnest the duration lists
    d0t = [item for sublist in d0t for item in sublist]
    d1t = [item for sublist in d1t for item in sublist]
    d2t = [item for sublist in d2t for item in sublist]
    d3t = [item for sublist in d3t for item in sublist]
    max_val = 1000
    d0t_plot = [val if val < max_val else max_val for val in d0t]
    d1t_plot = [val if val < max_val else max_val for val in d1t]
    d2t_plot = [val if val < max_val else max_val for val in d2t]
    d3t_plot = [val if val < max_val else max_val for val in d3t]


    # Setup graph
    fig, ax = plt.subplots(4, 2)
    fig.suptitle("Summary Distributions of Partnership Dynamics")


    # Duration for 16-20 age group
    ax[0, 0].hist(d0t_plot, bins = 10, range = (0, max_val))
    ax[0, 0].set_ylabel("Age 16-20")
    ax[0, 0].set_title("Partnership Duration")


    # Duration for 21-24 age group
    ax[1, 0].hist(d1t_plot, bins = 10, range = (0, max_val))
    ax[1, 0].set_ylabel("Age 21-24")


    # Duration for 25-29 age group
    ax[2, 0].hist(d2t_plot, bins = 10, range = (0, max_val))
    ax[2, 0].set_ylabel("Age 25-29")


    # Duration for 30 or older
    ax[3, 0].hist(d3t_plot, bins = 10, range = (0, max_val))
    ax[3, 0].set_ylabel("Age 30 and Over")
    ax[3, 0].set_xlabel("Duration (days)")


    # Partnership numbers for all age groups
    ax[0, 1].hist(meta.loc[meta["age_group"] == 0, "counter"], range = (0, max(meta.counter)))
    ax[0, 1].set_title("Number of Partners")
    ax[1, 1].hist(meta.loc[meta["age_group"] == 1, "counter"], range = (0, max(meta.counter)))
    ax[2, 1].hist(meta.loc[meta["age_group"] == 2, "counter"], range = (0, max(meta.counter)))
    ax[3, 1].hist(meta.loc[meta["age_group"] == 3, "counter"], range = (0, max(meta.counter)))
    ax[3, 1].set_xlabel("Number")


#%% GRAPH Relationship Occupancy


# If tracking partnership rates
if track_partnership_rates:


    # Setup graph
    t = range(0, n_days)
    fig, ax = plt.subplots(4, 2)


    # Graph of partnershipd in high risk 16-20
    n_g = sum((meta.age_group == 0) & (meta.risk == 1))
    ax[0, 0].bar(t, (p0ht[:,0] + p0ht[:,1] + p0ht[:,2] + p0ht[:,3] + p0ht[:,4])/n_g, color = "navy", width = 1)
    ax[0, 0].bar(t, (p0ht[:,1] + p0ht[:,2] + p0ht[:,3] + p0ht[:,4])/n_g, color = "darkred", width = 1)
    ax[0, 0].bar(t, (p0ht[:,1] + p0ht[:,2] + p0ht[:,3])/n_g, color = "red", width = 1)
    ax[0, 0].bar(t, (p0ht[:,1] + p0ht[:,2])/n_g, color = "darkgreen", width = 1)
    ax[0, 0].bar(t, p0ht[:,1]/n_g, color = "green", width = 1)
    ax[0, 0].set_title("High-risk")
    ax[0, 0].set_ylabel("Age 16-19")
    ax[0, 0].set_xlim([0, n_days])
    ax[0, 0].set_ylim([0, 1])


    # Graph of partnershipd in low risk 16-20
    n_g = sum((meta.age_group == 0) & (meta.risk == 0))
    ax[0, 1].bar(t, (p0lt[:,0] + p0lt[:,1] + p0lt[:,2] + p0lt[:,3] + p0lt[:,4])/n_g, color = "navy", width = 1)
    ax[0, 1].bar(t, (p0lt[:,1] + p0lt[:,2] + p0lt[:,3] + p0lt[:,4])/n_g, color = "darkred", width = 1)
    ax[0, 1].bar(t, (p0lt[:,1] + p0lt[:,2] + p0lt[:,3])/n_g, color = "red", width = 1)
    ax[0, 1].bar(t, (p0lt[:,1] + p0lt[:,2])/n_g, color = "darkgreen", width = 1)
    ax[0, 1].bar(t, p0lt[:,1]/n_g, color = "green", width = 1)
    ax[0, 1].set_title("Low-risk")
    ax[0, 1].set_xlim([0, n_days])
    ax[0, 1].set_ylim([0, 1])


    # Graph of partnershipd in high risk 21-24
    n_g = sum((meta.age_group == 1) & (meta.risk == 1))
    ax[1, 0].bar(t, (p1ht[:,0] + p1ht[:,1] + p1ht[:,2] + p1ht[:,3] + p1ht[:,4])/n_g, color = "navy", width = 1)
    ax[1, 0].bar(t, (p1ht[:,1] + p1ht[:,2] + p1ht[:,3] + p1ht[:,4])/n_g, color = "darkred", width = 1)
    ax[1, 0].bar(t, (p1ht[:,1] + p1ht[:,2] + p1ht[:,3])/n_g, color = "red", width = 1)
    ax[1, 0].bar(t, (p1ht[:,1] + p1ht[:,2])/n_g, color = "darkgreen", width = 1)
    ax[1, 0].bar(t, p1ht[:,1]/n_g, color = "green", width = 1)
    ax[1, 0].set_ylabel("Age 20-24")
    ax[1, 0].set_xlim([0, n_days])
    ax[1, 0].set_ylim([0, 1])


    # Graph of partnershipd in low risk 21-24
    n_g = sum((meta.age_group == 1) & (meta.risk == 0))
    ax[1, 1].bar(t, (p1lt[:,0] + p1lt[:,1] + p1lt[:,2] + p1lt[:,3] + p1lt[:,4])/n_g, color = "navy", width = 1)
    ax[1, 1].bar(t, (p1lt[:,1] + p1lt[:,2] + p1lt[:,3] + p1lt[:,4])/n_g, color = "darkred", width = 1)
    ax[1, 1].bar(t, (p1lt[:,1] + p1lt[:,2] + p1lt[:,3])/n_g, color = "red", width = 1)
    ax[1, 1].bar(t, (p1lt[:,1] + p1lt[:,2])/n_g, color = "darkgreen", width = 1)
    ax[1, 1].bar(t, p1lt[:,1]/n_g, color = "green", width = 1)
    ax[1, 1].set_xlim([0, n_days])
    ax[1, 1].set_ylim([0, 1])


    # Graph of partnershipd in high risk 25-29
    n_g = sum((meta.age_group == 2) & (meta.risk == 1))
    ax[2, 0].bar(t, (p2ht[:,0] + p2ht[:,1] + p2ht[:,2] + p2ht[:,3] + p2ht[:,4])/n_g, label = "Single", color = "navy", width = 1)
    ax[2, 0].bar(t, (p2ht[:,1] + p2ht[:,2] + p2ht[:,3] + p2ht[:,4])/n_g, label = "Short-term concurrent", color = "darkred", width = 1)
    ax[2, 0].bar(t, (p2ht[:,1] + p2ht[:,2] + p2ht[:,3])/n_g, label = "Short-term", color = "red", width = 1)
    ax[2, 0].bar(t, (p2ht[:,1] + p2ht[:,2])/n_g, label = "Long-term concurrent", color = "darkgreen", width = 1)
    ax[2, 0].bar(t, p2ht[:,1]/n_g, label = "Long-term", color = "green", width = 1)
    ax[2, 0].legend(loc = "lower left", bbox_to_anchor=(0.18,-0.05), bbox_transform=fig.transFigure, ncol=3)
    ax[2, 0].set_xlim([0, n_days])
    ax[2, 0].set_ylabel("Age 25-29")
    ax[2, 0].set_xlabel("Time (days)")
    ax[2, 0].set_ylim([0, 1])


    # Graph of partnershipd in low risk 25-29
    n_g = sum((meta.age_group == 2) & (meta.risk == 0))
    ax[2, 1].bar(t, (p2lt[:,0] + p2lt[:,1] + p2lt[:,2] + p2lt[:,3] + p2lt[:,4])/n_g, color = "navy", width = 1)
    ax[2, 1].bar(t, (p2lt[:,1] + p2lt[:,2] + p2lt[:,3] + p2lt[:,4])/n_g, color = "darkred", width = 1)
    ax[2, 1].bar(t, (p2lt[:,1] + p2lt[:,2] + p2lt[:,3])/n_g, color = "red", width = 1)
    ax[2, 1].bar(t, (p2lt[:,1] + p2lt[:,2])/n_g, color = "darkgreen", width = 1)
    ax[2, 1].bar(t, p2lt[:,1]/n_g, color = "green", width = 1)
    ax[2, 1].set_xlim([0, n_days])
    ax[2, 1].set_ylim([0, 1])
    ax[2, 1].set_xlabel("Time (days)")


    # Graph of partnershipd in high risk 30 or older
    n_g = sum((meta.age_group == 3) & (meta.risk == 1))
    ax[3, 0].bar(t, (p3ht[:,0] + p3ht[:,1] + p3ht[:,2] + p3ht[:,3] + p3ht[:,4])/n_g, label = "Single", color = "navy", width = 1)
    ax[3, 0].bar(t, (p3ht[:,1] + p3ht[:,2] + p3ht[:,3] + p3ht[:,4])/n_g, label = "Short-term concurrent", color = "darkred", width = 1)
    ax[3, 0].bar(t, (p3ht[:,1] + p3ht[:,2] + p3ht[:,3])/n_g, label = "Short-term", color = "red", width = 1)
    ax[3, 0].bar(t, (p3ht[:,1] + p3ht[:,2])/n_g, label = "Long-term concurrent", color = "darkgreen", width = 1)
    ax[3, 0].bar(t, p3ht[:,1]/n_g, label = "Long-term", color = "green", width = 1)
    ax[3, 0].legend(loc = "lower left", bbox_to_anchor=(0.18,-0.05), bbox_transform=fig.transFigure, ncol=3)
    ax[3, 0].set_xlim([0, n_days])
    ax[3, 0].set_ylabel("Age Over 29")
    ax[3, 0].set_xlabel("Time (days)")
    ax[3, 0].set_ylim([0, 1])


    # Graph of partnershipd in low risk 30 or older
    n_g = sum((meta.age_group == 3) & (meta.risk == 0))
    ax[3, 1].bar(t, (p3lt[:,0] + p3lt[:,1] + p3lt[:,2] + p3lt[:,3] + p3lt[:,4])/n_g, color = "navy", width = 1)
    ax[3, 1].bar(t, (p3lt[:,1] + p3lt[:,2] + p3lt[:,3] + p3lt[:,4])/n_g, color = "darkred", width = 1)
    ax[3, 1].bar(t, (p3lt[:,1] + p3lt[:,2] + p3lt[:,3])/n_g, color = "red", width = 1)
    ax[3, 1].bar(t, (p3lt[:,1] + p3lt[:,2])/n_g, color = "darkgreen", width = 1)
    ax[3, 1].bar(t, p3lt[:,1]/n_g, color = "green", width = 1)
    ax[3, 1].set_xlim([0, n_days])
    ax[3, 1].set_ylim([0, 1])
    ax[3, 1].set_xlabel("Time (days)")


    # Save graph
    fig.savefig("graphs/relationship_distribution.pdf")











