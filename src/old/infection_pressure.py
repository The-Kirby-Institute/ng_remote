# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 15:45:42 2020

@author: nicol
"""

import random
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# Make some infection pressure functions that look roughly correct
t0 = 0
site0 = 0
site1 = 1
site2 = 0
site0_t0 = 0
site1_t0 = 1
site2_t0 = 2

t = np.linspace(0, 30)

pressure0 = 1/(1 + np.exp(0.5*(14-t)))
pressure1 = 1/(1 + np.exp(1.5*(4-t)))
pressure2 = 1/(1 + 10*np.exp(0.05*(60-t)))

fig, ax = plt.subplots(2)
ax[0].plot(t, pressure0, label="Anal")
ax[0].plot(t, pressure1, label="Genital")
ax[0].plot(t, pressure2, label="Pharynx")
ax[0].set_title("Infection Pressure/Discomfort")
ax[0].legend()


# Try sampling from the infection pressure distribution
df = pd.DataFrame(data = {"col0": pressure0,
                          "col1": pressure1,
                          "col2": pressure2})
df["pressure"] = df.apply(lambda row: np.max([row.col0, row.col1]), axis=1)

u = pd.DataFrame(data = {"rand": sp.random.random(10000)})
u["treatment"] = u.apply(lambda row: np.min(np.where(row.rand < df["pressure"])), axis=1)

ax[1].hist(u["treatment"], bins = 30, range=[0, 30])
ax[1].set_title("Treatment Time Distribution")


# Test the actual function for the simulation code
pressure = site0/(1 + np.exp(0.5*(14-site0_t0)))\
           + site1/(1 + np.exp(0.7*(7-site1_t0)))\
           + site2/(1 + 10*np.exp(0.05*(60-site2_t0)))

# Draw a random number from this distribution
u = random.random()
i = np.searchsorted(pressure, u)















