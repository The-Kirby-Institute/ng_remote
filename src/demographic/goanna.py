# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 10:51:38 2020

@author: nicol

Module to create a meta population in line with the GOANNA survey
"""

#%% Packages
import pandas as pd
import numpy as np
# import plotnine as pn


#%% Set constants


# Probability matrix of sexual orientation by gender and age group
# Data from Table 3-3 of the GOANNA study - demographic characteristics by gender and age group
# Note that this data is only given by either gender or age group
# rows = age group and columns = orientation
# Orientation: 0=heterosexual, 1=homosexual, 2=bisexual
# orientation_by_age = np.array([[1156/1231, 30/1231, 45/1231],
#                                [792/881, 55/881, 34/881],
#                                [641/714, 42/714, 31/714]])


# Update for second GOANNA survey
orientation_by_age = np.array([[418/469, 20/469, 31/469],
                               [334/374, 18/374, 22/374],
                               [239/283, 22/283, 22/283],
                               [239/283, 22/283, 22/283]])
orientation_by_age = np.cumsum(orientation_by_age, axis = 1)


# Probability distribution of sexual experience
# Data from Table 5-1
# Read as: age group (row) and number of sexual partners (column)
# Number of partners state space: {0 (none), 1, (one), 2 (two-four), 3 (>=five)}
# partner_prob_raw = np.array([[78, 346, 367, 76],
#                              [55, 366, 292, 64],
#                              [58, 363, 210, 27]])


# Updated for second GOANNA survey
partner_prob_raw = np.array([[13, 106, 112, 59],
                             [22, 154, 114, 70],
                             [17, 131, 71, 44],
                             [17, 131, 71, 44]])
partner_prob_raw = partner_prob_raw/partner_prob_raw.sum(axis=1, keepdims = True)




#%%###########################################################################
##                                  GOANNA                                  ##
##############################################################################
# SET UP META POPULATION IN LINE WITH GOANNA SURVEY
#
# gender and age from a uniform distribution.
#
# orientation in line with GOANNA population.
#
# risk in line with GOANNA population. In particular, the ones who said
# they had 5 or more partners.
#
#
# INPUT
#   n = population size
#
# OUTPUT
#   meta
#
def goanna(n, burn_in):

    #%% Initilise meta
    #
    # gender: 0=male, 1=female
    # age
    # age_group: 0=(16-19), 1=(20-24), 2=(25-29), 3=(>29)
    # orientation: 0=hetero, 1=homo, 2=bi
    # risk: 0=low, 1=high
    # partner: -1=no partner, otherwise the index of their partner
    # counter: (summary statistic - number of partners)
    #


    # Population size
    # n = 2000


    # Initilise meta population matrix
    meta = pd.DataFrame(columns = ["gender", "age", "age_group", "orientation", "risk", "partner", "counter",
                                   "state", "site0", "site1", "site2", "site0_t0", "site1_t0", "site2_t0",
                                   "site0_t1", "site1_t1", "site2_t1", "treatment_threshold", "recovery_time"])


    # Set variable types
    meta.gender.astype("int64")
    meta.age.astype("float64")
    meta.age_group.astype("int64")
    meta.orientation.astype("int64")
    meta.risk.astype("int64")
    meta.partner.astype("int64")
    meta.counter.astype("int64")
    meta.state.astype("category")
    meta.site0.astype("int64")
    meta.site1.astype("int64")
    meta.site2.astype("int64")
    meta.site0_t0.astype("float64")
    meta.site1_t0.astype("float64")
    meta.site2_t0.astype("float64")
    meta.site0_t1.astype("float64")
    meta.site1_t1.astype("float64")
    meta.site2_t1.astype("float64")
    meta.treatment_threshold.astype("float64")
    meta.recovery_time.astype("float64")


    # Initilise columns to be filled in during the simulation
    meta.partner = n * [-1]
    meta.counter = n * [0]


    #%% Set age, age_group and gender


    # Gender
    meta.gender = np.round(np.random.random(n))


    # Age
    meta.age = 16 + (35-16)*np.random.random(n)


    # Age-group
    for i in range(0, n):
        if meta.at[i, "age"] <= 19:
            meta.at[i, "age_group"] = int(0)
        elif meta.at[i, "age"] <= 24:
            meta.at[i, "age_group"] = int(1)
        elif meta.at[i, "age"] <= 29:
            meta.at[i, "age_group"] = int(2)
        else:
            meta.at[i, "age_group"] = int(3)


    #%% Set Orientation


    # Work out sexual orientations
    meta.orientation = np.random.random(n)
    for i in range(0, 3):


        # Extract orientation data
        x = meta.loc[meta["age_group"] == i, "orientation"]


        # Determine where each random number falls in the distribution
        x.loc[x <= orientation_by_age[i, 0]] = 0
        x.loc[x >= orientation_by_age[i, 1]] = 2
        x.loc[~np.isin(x, [0, 2])] = 1


        # Put back into the data frame
        meta.loc[x.index, "orientation"] = x


    #%% Set risk


    # Risk
    meta.risk = np.random.random(n)
    meta.loc[(meta["age_group"] == 0) & (meta["risk"] < partner_prob_raw[0,3]), "risk"] = 1
    meta.loc[(meta["age_group"] == 1) & (meta["risk"] < partner_prob_raw[1,3]), "risk"] = 1
    meta.loc[(meta["age_group"] == 2) & (meta["risk"] < partner_prob_raw[2,3]), "risk"] = 1
    meta.loc[(meta["age_group"] == 3) & (meta["risk"] < partner_prob_raw[3,3]), "risk"] = 1
    meta.loc[meta["risk"] != 1, "risk"] = 0


    #%% Set infection status


    # Set default values
    meta.state = n * ["S"]
    meta.site0 = n * [0]
    meta.site1 = n * [0]
    meta.site2 = n * [0]
    meta.site0_t0 = n * [float("inf")]
    meta.site1_t0 = n * [float("inf")]
    meta.site2_t0 = n * [float("inf")]
    meta.site0_t1 = n * [float("inf")]
    meta.site1_t1 = n * [float("inf")]
    meta.site2_t1 = n * [float("inf")]
    meta.treatment_threshold = np.random.random(n)
    meta.recovery_time = n * [float("inf")]


    # Seed some infections
    for i in range(0, n):


        # Set infection rate to 90% of the population
        if np.random.random() < 0.9:


            # Set as infected
            meta.at[i, "state"] = "E"


            # Choose a site at random
            u = np.cumsum(np.random.random(3))
            u = np.min(np.where(u/u[2] > np.random.random()))
            meta.at[i, "site" + str(u) + "_t0"] = burn_in + 100*np.random.random(1)


    return meta


#%%###########################################################################
##                  MAKE SUMMARY GRAPHS TO TEST META                        ##
##############################################################################


# def test_goanna(meta):
    # MAKE SUMMARY GRAPHS OF THE GOANNA META POPULATION


    # Check distribution
    # groupings = meta.groupby(["gender", "age_group"])["age"].count().reset_index(name="count")
    # fig = (\
    #     pn.ggplot(groupings, pn.aes(x="age_group", y="count", group="gender", fill="factor(gender)")) +\
    #     pn.geom_col(position = "dodge") \
    #     )
    # print(fig)


    # # Check on age distributions
    # groupings = meta.groupby(["gender", "age_group", "orientation"])["age"].count().reset_index(name="count")
    # fig = (\
    #     pn.ggplot(groupings, pn.aes(x="age_group", y="count", fill="factor(orientation)")) +\
    #     pn.geom_col(position = "stack") +\
    #     pn.facet_wrap("gender")
    #     )
    # print(fig)


    # # Check risk distributions
    # groupings = meta.groupby(["age_group", "risk", "gender"])["age"].count().reset_index(name="count")
    # fig = (\
    #     pn.ggplot(groupings, pn.aes(x="age_group", y="count", fill="factor(risk)")) +\
    #     pn.geom_col(position = "stack") +\
    #     pn.facet_wrap("gender")
    #     )
    # print(fig)


    # # Check infection numbers
    # groupings = meta.groupby(["age_group", "risk", "state"])["age"].count().reset_index(name="count")
    # fig = (\
    #     pn.ggplot(groupings, pn.aes(x="age_group", y="count", fill="factor(state)")) +\
    #     pn.geom_col(position = "stack") +\
    #     pn.facet_wrap("risk")
    #     )
    # print(fig)

