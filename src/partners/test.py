# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:27:30 2020

@author: nicol

Script for setting up how partnerships are made.
"""

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import src.partners.partners as ng
import src.partners.summary_stats as pstat
import src.demographic.goanna as pop


#%% Set up meta population

n = 2000
meta = pop.goanna(n)


#%% Test partnership biasing function


# Make a big array for setting partnership conncetions
partner_matrix = np.zeros((n,n))


# Run algorithm for a bunch of people
if False:
    for k in range(0, 200):
        for i in range(0, n):


            # Find a partner
            j = ng.find_partner(meta, partner_matrix, i)


            if j != -1:


                # Decide on their relationship
                is_short = ng.choose_relationship(meta, i, j)


                # Update population array
                partner_matrix[i,j] = 1
                partner_matrix[j,i] = 1


                # Update partners array
                if is_short == 0:
                    meta.at[i, "partner"] = j
                    meta.at[j, "partner"] = i

    # plt.imshow(partner_matrix)


#%% Test Partnership aquisition rates


# Default partnership durations to Inf
partner_expire = float("inf") * np.ones((n, n))


# Set simulation length
t_reset = 300


# Set time to start counting from
n_days = t_reset + 365


# Preallocate for summary statistics
# Number of partners last year
xt = np.zeros((n_days, 4))
g0t = np.zeros((n_days, 4))
g1t = np.zeros((n_days, 4))
g2t = np.zeros((n_days, 4))


# Duration of partnerships
d0t = []
d1t = []
d2t = []


# In each partnership
p0ht = np.zeros((n_days, 5))
p0lt = np.zeros((n_days, 5))
p1ht = np.zeros((n_days, 5))
p1lt = np.zeros((n_days, 5))
p2ht = np.zeros((n_days, 5))
p2lt = np.zeros((n_days, 5))


# Iterate over all time points
n_people = min(n, 2000)
for t in tqdm(range(0, n_days)):


    # Create a new relationship
    meta, partner_matrix, partner_expire, d0ti, d1ti, d2ti = ng.new_partnership(meta, partner_matrix, partner_expire, i, t)

    d0t.append(d0ti)
    d1t.append(d1ti)
    d2t.append(d2ti)

    # Identify and remove expired relationships
    meta, partner_matrix, partner_expire = ng.old_partnerships(meta, partner_matrix, partner_expire, t)


    # Update summary stats on the number of annual partners
    xt, g0t, g1t, g2t = pstat.update_cumulative_partners(meta, partner_matrix, partner_expire, t, n_people, xt, g0t, g1t, g2t)
    p0ht, p0lt, p1ht, p1lt, p2ht, p2lt = pstat.update_partnership_types(meta, partner_matrix, t, p0ht, p0lt, p1ht, p1lt, p2ht, p2lt)
    if t == t_reset: meta.counter = np.sum(partner_matrix, axis=1)


# Plot of the partnership matrix
plt.imshow(partner_matrix)






#%%
# Plot long vs short term partnerships over time
fig, ax = plt.subplots(4)
ax[0].plot(range(0, n_days), xt[:,0], label = "single")
ax[0].plot(range(0, n_days), xt[:,1], label = "long term")
ax[0].plot(range(0, n_days), xt[:,2], label = "short term")
ax[0].plot(range(0, n_days), xt[:,3], label = "time left")
ax[0].set_ylim([0, n])
ax[0].legend(loc = "upper right")


# Plot the number of people in each partnership group over time
ax[1].plot(range(0, n_days), g0t[:,0], label = "0", color = "blue")
ax[1].plot(range(0, n_days), g0t[:,1], label = "1", color = "red")
ax[1].plot(range(0, n_days), g0t[:,2], label = "2-4", color = "green")
ax[1].plot(range(0, n_days), g0t[:,3], label = "5 more", color = "purple")
weight = sum(g0t[0,:])
ax[1].axhline(weight*pop.partner_prob_raw[0,0], color = "blue", linestyle = "--")
ax[1].axhline(weight*pop.partner_prob_raw[0,1], color = "red", linestyle = "--")
ax[1].axhline(weight*pop.partner_prob_raw[0,2], color= "green", linestyle = "--")
ax[1].axhline(weight*pop.partner_prob_raw[0,3], color = "purple", linestyle = "--")
ax[1].set_xlim([t_reset, n_days])
ax[1].legend(loc = "upper right")
ax[0].set_title("Comparison of 12-month partner count to GOANNA data")


# Plot the number of people in each partnership group over time
ax[2].plot(range(0, n_days), g1t[:,0], label = "0", color = "blue")
ax[2].plot(range(0, n_days), g1t[:,1], label = "1", color = "red")
ax[2].plot(range(0, n_days), g1t[:,2], label = "2-4", color = "green")
ax[2].plot(range(0, n_days), g1t[:,3], label = "5 more", color = "purple")
weight = sum(g1t[0,:])
ax[2].axhline(weight*pop.partner_prob_raw[1,0], color = "blue", linestyle = "--")
ax[2].axhline(weight*pop.partner_prob_raw[1,1], color = "red", linestyle = "--")
ax[2].axhline(weight*pop.partner_prob_raw[1,2], color= "green", linestyle = "--")
ax[2].axhline(weight*pop.partner_prob_raw[1,3], color = "purple", linestyle = "--")
ax[2].set_xlim([t_reset, n_days])


# Plot the number of people in each partnership group over time
ax[3].plot(range(0, n_days), g2t[:,0], label = "0", color = "blue")
ax[3].plot(range(0, n_days), g2t[:,1], label = "1", color = "red")
ax[3].plot(range(0, n_days), g2t[:,2], label = "2-4", color = "green")
ax[3].plot(range(0, n_days), g2t[:,3], label = "5 more", color = "purple")
weight = sum(g2t[0,:])
ax[3].axhline(weight*pop.partner_prob_raw[2,0], color = "blue", linestyle = "--")
ax[3].axhline(weight*pop.partner_prob_raw[2,1], color = "red", linestyle = "--")
ax[3].axhline(weight*pop.partner_prob_raw[2,2], color= "green", linestyle = "--")
ax[3].axhline(weight*pop.partner_prob_raw[2,3], color = "purple", linestyle = "--")
ax[3].set_xlim([t_reset, n_days])



#%%
# Partner durations
# fig, ax = plt.subplots(3)

# ax[0].hist(d0lt, alpha = 0.6, bins = int(n_days/5), range = (0, n_days), label = "low-risk", density = True)
# ax[0].hist(d0ht, alpha = 0.6, bins = int(n_days/5), range = (0, n_days), label = "high-risk", density = True)
# ax[0].legend()

# ax[1].hist(d1lt, alpha = 0.6, bins = int(n_days/5), range = (0, n_days), label = "low-risk", density = True)
# ax[1].hist(d1ht, alpha = 0.6, bins = int(n_days/5), range = (0, n_days), label = "high-risk", density = True)

# ax[2].hist(d2lt, alpha = 0.6, bins = int(n_days/5), range = (0, n_days), label = "low-risk", density = True)
# ax[2].hist(d2ht, alpha = 0.6, bins = int(n_days/5), range = (0, n_days), label = "high-risk", density = True)



#%% Number of people in each relationship type
t = range(0, n_days)
fig, ax = plt.subplots(3, 2)
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


n_g = sum((meta.age_group == 0) & (meta.risk == 0))
ax[0, 1].bar(t, (p0lt[:,0] + p0lt[:,1] + p0lt[:,2] + p0lt[:,3] + p0lt[:,4])/n_g, color = "navy", width = 1)
ax[0, 1].bar(t, (p0lt[:,1] + p0lt[:,2] + p0lt[:,3] + p0lt[:,4])/n_g, color = "darkred", width = 1)
ax[0, 1].bar(t, (p0lt[:,1] + p0lt[:,2] + p0lt[:,3])/n_g, color = "red", width = 1)
ax[0, 1].bar(t, (p0lt[:,1] + p0lt[:,2])/n_g, color = "darkgreen", width = 1)
ax[0, 1].bar(t, p0lt[:,1]/n_g, color = "green", width = 1)
ax[0, 1].set_title("Low-risk")
ax[0, 1].set_xlim([0, n_days])
ax[0, 1].set_ylim([0, 1])


n_g = sum((meta.age_group == 1) & (meta.risk == 1))
ax[1, 0].bar(t, (p1ht[:,0] + p1ht[:,1] + p1ht[:,2] + p1ht[:,3] + p1ht[:,4])/n_g, color = "navy", width = 1)
ax[1, 0].bar(t, (p1ht[:,1] + p1ht[:,2] + p1ht[:,3] + p1ht[:,4])/n_g, color = "darkred", width = 1)
ax[1, 0].bar(t, (p1ht[:,1] + p1ht[:,2] + p1ht[:,3])/n_g, color = "red", width = 1)
ax[1, 0].bar(t, (p1ht[:,1] + p1ht[:,2])/n_g, color = "darkgreen", width = 1)
ax[1, 0].bar(t, p1ht[:,1]/n_g, color = "green", width = 1)
ax[1, 0].set_ylabel("Age 20-24")
ax[1, 0].set_xlim([0, n_days])
ax[1, 0].set_ylim([0, 1])

n_g = sum((meta.age_group == 1) & (meta.risk == 0))
ax[1, 1].bar(t, (p1lt[:,0] + p1lt[:,1] + p1lt[:,2] + p1lt[:,3] + p1lt[:,4])/n_g, color = "navy", width = 1)
ax[1, 1].bar(t, (p1lt[:,1] + p1lt[:,2] + p1lt[:,3] + p1lt[:,4])/n_g, color = "darkred", width = 1)
ax[1, 1].bar(t, (p1lt[:,1] + p1lt[:,2] + p1lt[:,3])/n_g, color = "red", width = 1)
ax[1, 1].bar(t, (p1lt[:,1] + p1lt[:,2])/n_g, color = "darkgreen", width = 1)
ax[1, 1].bar(t, p1lt[:,1]/n_g, color = "green", width = 1)
ax[1, 1].set_xlim([0, n_days])
ax[1, 1].set_ylim([0, 1])


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

n_g = sum((meta.age_group == 2) & (meta.risk == 0))
ax[2, 1].bar(t, (p2lt[:,0] + p2lt[:,1] + p2lt[:,2] + p2lt[:,3] + p2lt[:,4])/n_g, color = "navy", width = 1)
ax[2, 1].bar(t, (p2lt[:,1] + p2lt[:,2] + p2lt[:,3] + p2lt[:,4])/n_g, color = "darkred", width = 1)
ax[2, 1].bar(t, (p2lt[:,1] + p2lt[:,2] + p2lt[:,3])/n_g, color = "red", width = 1)
ax[2, 1].bar(t, (p2lt[:,1] + p2lt[:,2])/n_g, color = "darkgreen", width = 1)
ax[2, 1].bar(t, p2lt[:,1]/n_g, color = "green", width = 1)
ax[2, 1].set_xlim([0, n_days])
ax[2, 1].set_ylim([0, 1])
ax[2, 1].set_xlabel("Time (days)")

fig.savefig("graphs/relationship_distribution.pdf")






