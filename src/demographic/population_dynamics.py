# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 09:57:13 2021

@author: nicol
@profile
"""


#%% Setup Libraries


# Load libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy


# My modules
import src.demographic.generate_population as pop
import src.partners.partners as prt
import src.infections.ng as ng


# Parse general simulation parameters
sim_parameters = pd.read_csv('data/param.csv')




#%% FUN update_population()
#
#
# Function to implement all of the population dynamics
#
#
def update_population(t, pop_parameters, prt_parameters, inf_parameters, meta, partner_matrix, partner_expire):


    # Mobility dynamics
    meta, partner_matrix, partner_expire, pop_parameters = mobility(t, pop_parameters, prt_parameters, inf_parameters, meta, partner_matrix, partner_expire)
    
    
    # Change in sexual dynamics
    meta = update_sexual_acts(pop_parameters, meta)
    

    # Make people older by a day
    meta.loc[:, 'age'] = meta.loc[:, 'age'] + (1/365)
    meta.loc[:, 'age_group'] = np.floor((meta.age - 15)/5)


    # Update demographic indexing
    # pop_parameters = update_attribute_tracker(pop_parameters, meta, partner_matrix)


    return meta, partner_matrix, partner_expire




#%% FUN mobility()
#
#
# Implement an importation event.
#
# The imported case is infectious at one site and immediately joins a short
# term relationship at random.
#
#
def mobility(t, pop_parameters, prt_parameters, inf_parameters, meta, partner_matrix, partner_expire):


    # Make things a little easier
    imports = pop_parameters['imports']
    # lookup = pop_parameters['lookup']


    # Iterate over all age and gender combinations
    array_in = []
    list_out = []
    for gender in [0, 1]:


        # Age groups
        for age_group in [0, 1, 2, 3, 4]:


            # Set target
            target = pop_parameters['target'][str(gender)]


            # Sample the number of new people to come in
            rate = target[age_group] * pop_parameters['pop_turnover']['annual_arrivals'][0] / 365
            n_new = np.random.poisson(rate)


            # If this is the 16 y/o group, add in some more in
            if age_group == 0:
                n_debut = np.random.poisson(target[age_group]/(4*365))
            else:
                n_debut = 0


            # sample the people to add to the population
            temp = sample_from_imports(imports, age_group, gender, n_new + n_debut)


            # Change the age to 16 for some individuals
            if n_debut > 0:
                temp.loc[0:n_debut, 'age'] = float(16)


            # Concatinate with any other new people to come in
            if (len(temp) > 0) & (len(array_in) == 0):
                array_in = temp
            elif len(temp) > 0:
                array_in = array_in.append(temp).reset_index(drop = True)


            # Sample the number of people to leave
            rate = target[age_group] * pop_parameters['pop_turnover']['annual_departures'][0] / 365
            n_out = np.random.poisson(rate)
            
            
            # Sample the people to leave
            # cohort = meta.index[lookup['GA' + str(gender) + str(age_group)]]
            cohort = meta.index[(meta.gender == gender) & (meta.age_group == age_group)]
            list_out = list_out + ([] if (n_out == 0) | (len(cohort) == 0) else np.random.choice(cohort, n_out).tolist())


    # Identify anybody over 36 to also leave
    list_out = list_out + list(meta.index[meta.age > 36])


    # Change some attributes of the imports
    n_new = len(array_in)
    if n_new > 0:

        # Change infectious people to already be infectious upon entering the population
        infected = array_in.state == 'E'
        site0 = array_in.site0 == int(1)
        site1 = array_in.site1 == int(1)
        site2 = array_in.site2 == int(1)
        array_in.loc[infected, 'state'] = 'I'
        array_in.loc[(infected) & (site0), 'site0_t0'] = t
        array_in.loc[(infected) & (site1), 'site1_t0'] = t
        array_in.loc[(infected) & (site2), 'site2_t0'] = t
        array_in.loc[(infected) & (site0), 'site0_t1'] = t + ng.duration_rectal(inf_parameters, array_in.loc[(infected) & (site0), ], np.sum((infected) & (site0)))
        array_in.loc[(infected) & (site1), 'site1_t1'] = t + ng.duration_urethral(inf_parameters, array_in.loc[(infected) & (site1), ], np.sum((infected) & (site1)))
        array_in.loc[(infected) & (site2), 'site2_t1'] = t + ng.duration_pharyngeal(inf_parameters, array_in.loc[(infected) & (site2), ], np.sum((infected) & (site2)))
        array_in.loc[:, 'import_time'] = t
        
        
        # Set the sexual activity status of the new people
        array_in.loc[:, 'has_oral'] = np.random.random(n_new) < pop_parameters['sexual_acts_dist'].loc[0].values[1:6][list(array_in.age_group.astype('int32'))]
        array_in.loc[:, 'has_sex'] = np.random.random(n_new) < pop_parameters['sexual_acts_dist'].loc[1].values[1:6][list(array_in.age_group.astype('int32'))]
        array_in.loc[:, 'has_anal'] = np.random.random(n_new) < pop_parameters['sexual_acts_dist'].loc[2].values[1:6][list(array_in.age_group.astype('int32'))]
        
        
        # Default 16 year olds to no activity
        array_in.loc[0:n_debut, 'has_oral'] = False
        array_in.loc[0:n_debut, 'has_sex'] = False
        array_in.loc[0:n_debut, 'has_anal'] = False


        # Update the meta-population data
        meta, partner_matrix, partner_expire = add_to_meta(pop_parameters, prt_parameters, meta, partner_matrix, partner_expire, array_in, t)


    # Take some people out if needed
    if len(list_out) > 0:


        # Store data on people leaving
        list_out = list(np.unique(list_out))
        pop_parameters = record_data_of_leavers(pop_parameters, meta, list_out, t)


        # Update partner indicies in meta
        meta, partner_matrix, partner_expire = remove_from_meta(meta, partner_matrix, partner_expire, list_out)


    return meta, partner_matrix, partner_expire, pop_parameters



#%% FUN initilise_demographic_dynamics()
#
#
# Pre-compute a bunch of people to use as imports later.
#
#
def initilise_demographic_dynamics(pop_parameters, inf_parameters, meta, partner_matrix, t0):


    # Preallocate
    n_gen = 10000
    t0 = -1
    imports = dict()


    # Generate list of people for future importations
    for i in [0, 1]:
        for j in [0, 1, 2, 3, 4]:
            imports.update({str(i) + str(j) : import_particular_age_group_and_gender(pop_parameters, inf_parameters, j, i, n_gen, t0)})


    # Save in pop_parameters
    pop_parameters.update({'imports': imports})


    # Preallocate for storing people who left the population
    leavers = copy.deepcopy(meta.iloc[[0]])
    leavers['export_time'] = t0
    pop_parameters.update({'leavers': leavers})


    # Run the demographic indexer
    pop_parameters = update_attribute_tracker(pop_parameters, meta, partner_matrix)


    return pop_parameters





#%% FUN update_sexual_acts()
#
#
# Function which updates whether or not people engage in different sexual acts
#
#
def update_sexual_acts(pop_parameters, meta):
    
    
    # Extract the importation rate
    beta = pop_parameters['pop_turnover'].annual_arrivals[0]/365
    
    
    #######################
    ##  16-19 AGE GROUP  ##
    #######################
    
    
    # Extract who engages in what
    sub_group = meta.loc[(meta.age_group == 0), :]
    oral = sub_group.index[sub_group.has_oral == False]
    anal = sub_group.index[sub_group.has_anal == False]
    sex = sub_group.index[sub_group.has_sex == False]
    
    
    # Sample people to switch
    alpha = 1/(4*365)
    oral_new = np.random.random(len(oral)) < sexual_acts_turn_on_0(alpha, beta, pop_parameters['sexual_acts_dist'].loc[0].values[1])
    sex_new = np.random.random(len(sex)) < sexual_acts_turn_on_0(alpha, beta, pop_parameters['sexual_acts_dist'].loc[1].values[1])
    anal_new = np.random.random(len(anal)) < sexual_acts_turn_on_0(alpha, beta, pop_parameters['sexual_acts_dist'].loc[2].values[1])
    
    
    # Update meta
    meta.loc[oral[oral_new], 'has_oral'] = True
    meta.loc[anal[anal_new], 'has_anal'] = True
    meta.loc[sex[sex_new], 'has_sex'] = True
    
    
    #######################
    ##  20-24 AGE GROUP  ##
    #######################
    
    
    # Update whether people engage in different acts
    sub_group = meta.loc[(meta.age_group == 1), :]
    oral = sub_group.index[sub_group.has_oral == False]
    anal = sub_group.index[sub_group.has_anal == False]
    sex = sub_group.index[sub_group.has_sex == False]
    
    
    # Sample people to switch
    alpha = 1/(5*365)
    oral_new = np.random.random(len(oral)) < sexual_acts_turn_on_1(alpha, beta, pop_parameters['sexual_acts_dist'].loc[0].values[2], pop_parameters['sexual_acts_dist'].loc[0].values[1])
    sex_new = np.random.random(len(sex)) < sexual_acts_turn_on_1(alpha, beta, pop_parameters['sexual_acts_dist'].loc[1].values[2], pop_parameters['sexual_acts_dist'].loc[1].values[1])
    anal_new = np.random.random(len(anal)) < sexual_acts_turn_on_1(alpha, beta, pop_parameters['sexual_acts_dist'].loc[2].values[2], pop_parameters['sexual_acts_dist'].loc[2].values[1])
    
    
    # Update meta
    meta.loc[oral[oral_new], 'has_oral'] = True
    meta.loc[anal[anal_new], 'has_anal'] = True
    meta.loc[sex[sex_new], 'has_sex'] = True
    
    
    #######################
    ##  25-29 AGE GROUP  ##
    #######################
    # No change to aggregate oral and anal proportions in this cohort.
    
    
    # Update whether people engage in different acts
    sub_group = meta.loc[(meta.age_group == 2), :]
    # oral = sub_group.index[sub_group.has_oral == False]
    # anal = sub_group.index[sub_group.has_anal == False]
    sex = sub_group.index[sub_group.has_sex == False]
    
    
    # Sample people to switch
    # oral_new = np.random.random(len(oral)) < compute_prob(alpha, beta, pop_parameters['sexual_acts_dist'].loc[0].values[3] - pop_parameters['sexual_acts_dist'].loc[0].values[2])
    # anal_new = np.random.random(len(anal)) < compute_prob(alpha, beta, pop_parameters['sexual_acts_dist'].loc[2].values[3] - pop_parameters['sexual_acts_dist'].loc[2].values[2])
    sex_new = np.random.random(len(sex)) < sexual_acts_turn_on_2(alpha, beta, \
                                                                 pop_parameters['sexual_acts_dist'].loc[1].values[3], \
                                                                 pop_parameters['sexual_acts_dist'].loc[1].values[2], \
                                                                 pop_parameters['sexual_acts_dist'].loc[1].values[1])
        
    
    # Update meta
    # meta.loc[oral[oral_new], 'has_oral'] = True
    # meta.loc[anal[anal_new], 'has_anal'] = True
    meta.loc[sex[sex_new], 'has_sex'] = True
    
    
    # Done!
    return meta




#%% HELPER sexual_acts_turn_on_0()
#
#
# Computes the turn on rate for people in age group 0
#
#
def sexual_acts_turn_on_0(alpha, beta, p0):
    beta = (1-p0) * beta
    gamma0 = ( p0 * (alpha + beta))/((1 - p0) * (1 - alpha - beta))
    return gamma0




#%% HELPER sexual_acts_turn_on_1()
#
#
# Computes the turn on rate for people in age group 1
#
#
def sexual_acts_turn_on_1(alpha, beta, p1, p0):
    gamma0 = sexual_acts_turn_on_0(alpha, beta, p0)
    alpha1 = (1-gamma0) * alpha
    beta = (1-p1) * beta
    gamma1 = (alpha1 * (p1 - p0) + p1 * beta) / ((1 - p1) * (1 - beta - alpha1))
    return gamma1




#%% HELPER sexual_acts_turn_on_2()
#
#
# Computes the turn on rate for people in age group 2
#
#
def sexual_acts_turn_on_2(alpha, beta, p2, p1, p0):
    gamma0 = sexual_acts_turn_on_0(alpha, beta, p0)
    gamma1 = sexual_acts_turn_on_1(alpha, beta, p1, p0)
    alpha2 = (1-gamma0) * (1-gamma1) * alpha
    beta = (1-p2) * beta
    gamma2 = (alpha2 * (p2 - p1) + p2 * beta) / ((1 - p2) * (1 - beta - alpha2))
    return gamma2




#%% HELPER update_attribute_tracker()
#
#
#
#
#
def update_attribute_tracker(pop_parameters, meta, partner_matrix):
    # Experimental: some way of cutting out the number of logicals being used


    # Compute row sums for the partner matrix
    pmatrow = np.sum(partner_matrix[0:len(meta), 0:len(meta)], axis = 0)


    # Identifier for all the base line age and sex cohorts
    lookup = {
            'G0': meta.gender == 0,
            'G1': meta.gender == 1,
            'A0': meta.age_group == 0,
            'A1': meta.age_group == 1,
            'A2': meta.age_group == 2,
            'A3': (meta.age_group == 3) & (meta.age < 35),
            'A4': meta.age >= 35,
            'R0': meta.risk == 0,
            'R1': meta.risk == 1,
            'O0': meta.orientation == 0,
            'O1': meta.orientation == 1,
            'O2': meta.orientation == 2,
            'PL0': meta.partner == -1,
            'PL1': meta.partner > -1,
            'PS0': (pmatrow == 0) | ((pmatrow == 1) & (meta.partner > -1)),
            'PS1': (pmatrow >= 2) | ((pmatrow == 1) & (meta.partner == -1)),
            'I00': meta.site0 == 0,
            'I01': meta.site0 == 1,
            'I10': meta.site1 == 0,
            'I11': meta.site1 == 1,
            'I20': meta.site2 == 0,
            'I21': meta.site2 == 1,
            'SS': meta.state == 'S',
            'SE': meta.state == 'E',
            'SI': meta.state == 'I',
            'SR': meta.state == 'R',
            'ST': meta.state == 'T',
            'SV': meta.state == 'V',
            }


    # Identifiers 2-way interactions
    lookup.update({
                # All gender - age combinations
                'GA00': lookup['G0'] & lookup['A0'],
                'GA01': lookup['G0'] & lookup['A1'],
                'GA02': lookup['G0'] & lookup['A2'],
                'GA03': lookup['G0'] & lookup['A3'],
                'GA04': lookup['G0'] & lookup['A4'],
                'GA10': lookup['G1'] & lookup['A0'],
                'GA11': lookup['G1'] & lookup['A1'],
                'GA12': lookup['G1'] & lookup['A2'],
                'GA13': lookup['G1'] & lookup['A3'],
                'GA14': lookup['G1'] & lookup['A4'],
                
                # All risk - age combinations
                'RA00': lookup['R0'] & lookup['A0'],
                'RA01': lookup['R0'] & lookup['A1'],
                'RA02': lookup['R0'] & lookup['A2'],
                'RA03': lookup['R0'] & lookup['A3'],
                'RA04': lookup['R0'] & lookup['A4'],
                'RA10': lookup['R1'] & lookup['A0'],
                'RA11': lookup['R1'] & lookup['A1'],
                'RA12': lookup['R1'] & lookup['A2'],
                'RA13': lookup['R1'] & lookup['A3'],
                'RA14': lookup['R1'] & lookup['A4'],

                # # All orientation - age combinations
                # 'OA00': lookup['O0'] & lookup['A0'],
                # 'OA01': lookup['O0'] & lookup['A1'],
                # 'OA02': lookup['O0'] & lookup['A2'],
                # 'OA03': lookup['O0'] & lookup['A3'],
                # 'OA04': lookup['O0'] & lookup['A4'],
                # 'OA10': lookup['O1'] & lookup['A0'],
                # 'OA11': lookup['O1'] & lookup['A1'],
                # 'OA12': lookup['O1'] & lookup['A2'],
                # 'OA13': lookup['O1'] & lookup['A3'],
                # 'OA14': lookup['O1'] & lookup['A4'],
                # 'OA20': lookup['O2'] & lookup['A0'],
                # 'OA21': lookup['O2'] & lookup['A1'],
                # 'OA22': lookup['O2'] & lookup['A2'],
                # 'OA23': lookup['O2'] & lookup['A3'],
                # 'OA24': lookup['O2'] & lookup['A4'],

                # All long-term partner - age combinations
                'PLA00': lookup['PL0'] & lookup['A0'],
                'PLA01': lookup['PL0'] & lookup['A1'],
                'PLA02': lookup['PL0'] & lookup['A2'],
                'PLA03': lookup['PL0'] & lookup['A3'],
                'PLA04': lookup['PL0'] & lookup['A4'],
                'PLA10': lookup['PL1'] & lookup['A0'],
                'PLA11': lookup['PL1'] & lookup['A1'],
                'PLA12': lookup['PL1'] & lookup['A2'],
                'PLA13': lookup['PL1'] & lookup['A3'],
                'PLA14': lookup['PL1'] & lookup['A4'],

                # All short-term partner - age combinations
                'PSA00': lookup['PS0'] & lookup['A0'],
                'PSA01': lookup['PS0'] & lookup['A1'],
                'PSA02': lookup['PS0'] & lookup['A2'],
                'PSA03': lookup['PS0'] & lookup['A3'],
                'PSA04': lookup['PS0'] & lookup['A4'],
                'PSA10': lookup['PS1'] & lookup['A0'],
                'PSA11': lookup['PS1'] & lookup['A1'],
                'PSA12': lookup['PS1'] & lookup['A2'],
                'PSA13': lookup['PS1'] & lookup['A3'],
                'PSA14': lookup['PS1'] & lookup['A4'],

                # Other partner combinations
                'P0': lookup['PL0'] & lookup['PS0'],
                'P1': lookup['PL1'] | lookup['PS1'],
                'PC0': lookup['PL0'] | (lookup['PL1'] & lookup['PS0']),
                'PC1': lookup['PL1'] & lookup['PS1'],
                
                # All infection status - age combinations
                'I0A00': lookup['I00'] & lookup['A0'],
                'I0A01': lookup['I00'] & lookup['A1'],
                'I0A02': lookup['I00'] & lookup['A2'],
                'I0A03': lookup['I00'] & lookup['A3'],
                'I0A04': lookup['I00'] & lookup['A4'],
                'I0A10': lookup['I01'] & lookup['A0'],
                'I0A11': lookup['I01'] & lookup['A1'],
                'I0A12': lookup['I01'] & lookup['A2'],
                'I0A13': lookup['I01'] & lookup['A3'],
                'I0A14': lookup['I01'] & lookup['A4'],
                
                'I1A00': lookup['I10'] & lookup['A0'],
                'I1A01': lookup['I10'] & lookup['A1'],
                'I1A02': lookup['I10'] & lookup['A2'],
                'I1A03': lookup['I10'] & lookup['A3'],
                'I1A04': lookup['I10'] & lookup['A4'],
                'I1A10': lookup['I11'] & lookup['A0'],
                'I1A11': lookup['I11'] & lookup['A1'],
                'I1A12': lookup['I11'] & lookup['A2'],
                'I1A13': lookup['I11'] & lookup['A3'],
                'I1A14': lookup['I11'] & lookup['A4'],
                
                'I2A00': lookup['I20'] & lookup['A0'],
                'I2A01': lookup['I20'] & lookup['A1'],
                'I2A02': lookup['I20'] & lookup['A2'],
                'I2A03': lookup['I20'] & lookup['A3'],
                'I2A04': lookup['I20'] & lookup['A4'],
                'I2A10': lookup['I21'] & lookup['A0'],
                'I2A11': lookup['I21'] & lookup['A1'],
                'I2A12': lookup['I21'] & lookup['A2'],
                'I2A13': lookup['I21'] & lookup['A3'],
                'I2A14': lookup['I21'] & lookup['A4'],
                })


    # Identifiers of 3-way interactions
    lookup.update({
                # All gender - age - risk combinations
                'GAR000': lookup['GA00'] & lookup['R0'],
                'GAR010': lookup['GA01'] & lookup['R0'],
                'GAR020': lookup['GA02'] & lookup['R0'],
                'GAR030': lookup['GA03'] & lookup['R0'],
                'GAR040': lookup['GA04'] & lookup['R0'],
                'GAR100': lookup['GA10'] & lookup['R0'],
                'GAR110': lookup['GA11'] & lookup['R0'],
                'GAR120': lookup['GA12'] & lookup['R0'],
                'GAR130': lookup['GA13'] & lookup['R0'],
                'GAR140': lookup['GA14'] & lookup['R0'],
                
                'GAR001': lookup['GA00'] & lookup['R1'],
                'GAR011': lookup['GA01'] & lookup['R1'],
                'GAR021': lookup['GA02'] & lookup['R1'],
                'GAR031': lookup['GA03'] & lookup['R1'],
                'GAR041': lookup['GA04'] & lookup['R1'],
                'GAR101': lookup['GA10'] & lookup['R1'],
                'GAR111': lookup['GA11'] & lookup['R1'],
                'GAR121': lookup['GA12'] & lookup['R1'],
                'GAR131': lookup['GA13'] & lookup['R1'],
                'GAR141': lookup['GA14'] & lookup['R1'],
                
                # All gender - age - ORIENTATION combinations
                'GAO000': lookup['GA00'] & lookup['O0'],
                'GAO010': lookup['GA01'] & lookup['O0'],
                'GAO020': lookup['GA02'] & lookup['O0'],
                'GAO030': lookup['GA03'] & lookup['O0'],
                'GAO040': lookup['GA04'] & lookup['O0'],
                'GAO100': lookup['GA10'] & lookup['O0'],
                'GAO110': lookup['GA11'] & lookup['O0'],
                'GAO120': lookup['GA12'] & lookup['O0'],
                'GAO130': lookup['GA13'] & lookup['O0'],
                'GAO140': lookup['GA14'] & lookup['O0'],
                   
                'GAO001': lookup['GA00'] & lookup['O1'],
                'GAO011': lookup['GA01'] & lookup['O1'],
                'GAO021': lookup['GA02'] & lookup['O1'],
                'GAO031': lookup['GA03'] & lookup['O1'],
                'GAO041': lookup['GA04'] & lookup['O1'],
                'GAO101': lookup['GA10'] & lookup['O1'],
                'GAO111': lookup['GA11'] & lookup['O1'],
                'GAO121': lookup['GA12'] & lookup['O1'],
                'GAO131': lookup['GA13'] & lookup['O1'],
                'GAO141': lookup['GA14'] & lookup['O1'],
                   
                'GAO002': lookup['GA00'] & lookup['O2'],
                'GAO012': lookup['GA01'] & lookup['O2'],
                'GAO022': lookup['GA02'] & lookup['O2'],
                'GAO032': lookup['GA03'] & lookup['O2'],
                'GAO042': lookup['GA04'] & lookup['O2'],
                'GAO102': lookup['GA10'] & lookup['O2'],
                'GAO112': lookup['GA11'] & lookup['O2'],
                'GAO122': lookup['GA12'] & lookup['O2'],
                'GAO132': lookup['GA13'] & lookup['O2'],
                'GAO142': lookup['GA14'] & lookup['O2'],
                   
                # All gender - age - infection combinations
                'I0GA000': lookup['I0A00'] & lookup['G0'],
                'I0GA001': lookup['I0A01'] & lookup['G0'],
                'I0GA002': lookup['I0A02'] & lookup['G0'],
                'I0GA003': lookup['I0A03'] & lookup['G0'],
                'I0GA004': lookup['I0A04'] & lookup['G0'],
                'I0GA010': lookup['I0A00'] & lookup['G1'],
                'I0GA011': lookup['I0A01'] & lookup['G1'],
                'I0GA012': lookup['I0A02'] & lookup['G1'],
                'I0GA013': lookup['I0A03'] & lookup['G1'],
                'I0GA014': lookup['I0A04'] & lookup['G1'],
                   
                'I0GA100': lookup['I0A10'] & lookup['G0'],
                'I0GA101': lookup['I0A11'] & lookup['G0'],
                'I0GA102': lookup['I0A12'] & lookup['G0'],
                'I0GA103': lookup['I0A13'] & lookup['G0'],
                'I0GA104': lookup['I0A14'] & lookup['G0'],
                'I0GA110': lookup['I0A10'] & lookup['G1'],
                'I0GA111': lookup['I0A11'] & lookup['G1'],
                'I0GA112': lookup['I0A12'] & lookup['G1'],
                'I0GA113': lookup['I0A13'] & lookup['G1'],
                'I0GA114': lookup['I0A14'] & lookup['G1'],
                   
                'I1GA000': lookup['I1A00'] & lookup['G0'],
                'I1GA001': lookup['I1A01'] & lookup['G0'],
                'I1GA002': lookup['I1A02'] & lookup['G0'],
                'I1GA003': lookup['I1A03'] & lookup['G0'],
                'I1GA004': lookup['I1A04'] & lookup['G0'],
                'I1GA010': lookup['I1A00'] & lookup['G1'],
                'I1GA011': lookup['I1A01'] & lookup['G1'],
                'I1GA012': lookup['I1A02'] & lookup['G1'],
                'I1GA013': lookup['I1A03'] & lookup['G1'],
                'I1GA014': lookup['I1A04'] & lookup['G1'],
                   
                'I1GA100': lookup['I1A10'] & lookup['G0'],
                'I1GA101': lookup['I1A11'] & lookup['G0'],
                'I1GA102': lookup['I1A12'] & lookup['G0'],
                'I1GA103': lookup['I1A13'] & lookup['G0'],
                'I1GA104': lookup['I1A14'] & lookup['G0'],
                'I1GA110': lookup['I1A10'] & lookup['G1'],
                'I1GA111': lookup['I1A11'] & lookup['G1'],
                'I1GA112': lookup['I1A12'] & lookup['G1'],
                'I1GA113': lookup['I1A13'] & lookup['G1'],
                'I1GA114': lookup['I1A14'] & lookup['G1'],
                   
                'I2GA000': lookup['I0A00'] & lookup['G0'],
                'I2GA001': lookup['I0A01'] & lookup['G0'],
                'I2GA002': lookup['I0A02'] & lookup['G0'],
                'I2GA003': lookup['I0A03'] & lookup['G0'],
                'I2GA004': lookup['I0A04'] & lookup['G0'],
                'I2GA010': lookup['I0A00'] & lookup['G1'],
                'I2GA011': lookup['I0A01'] & lookup['G1'],
                'I2GA012': lookup['I0A02'] & lookup['G1'],
                'I2GA013': lookup['I0A03'] & lookup['G1'],
                'I2GA014': lookup['I0A04'] & lookup['G1'],
                   
                'I2GA100': lookup['I2A10'] & lookup['G0'],
                'I2GA101': lookup['I2A11'] & lookup['G0'],
                'I2GA102': lookup['I2A12'] & lookup['G0'],
                'I2GA103': lookup['I2A13'] & lookup['G0'],
                'I2GA104': lookup['I2A14'] & lookup['G0'],
                'I2GA110': lookup['I2A10'] & lookup['G1'],
                'I2GA111': lookup['I2A11'] & lookup['G1'],
                'I2GA112': lookup['I2A12'] & lookup['G1'],
                'I2GA113': lookup['I2A13'] & lookup['G1'],
                'I2GA114': lookup['I2A14'] & lookup['G1'],
                
                # All risk - age - infection combinations
                'I0RA000': lookup['I0A00'] & lookup['R0'],
                'I0RA001': lookup['I0A01'] & lookup['R0'],
                'I0RA002': lookup['I0A02'] & lookup['R0'],
                'I0RA003': lookup['I0A03'] & lookup['R0'],
                'I0RA004': lookup['I0A04'] & lookup['R0'],
                'I0RA010': lookup['I0A00'] & lookup['R1'],
                'I0RA011': lookup['I0A01'] & lookup['R1'],
                'I0RA012': lookup['I0A02'] & lookup['R1'],
                'I0RA013': lookup['I0A03'] & lookup['R1'],
                'I0RA014': lookup['I0A04'] & lookup['R1'],
                   
                'I0RA100': lookup['I0A10'] & lookup['R0'],
                'I0RA101': lookup['I0A11'] & lookup['R0'],
                'I0RA102': lookup['I0A12'] & lookup['R0'],
                'I0RA103': lookup['I0A13'] & lookup['R0'],
                'I0RA104': lookup['I0A14'] & lookup['R0'],
                'I0RA110': lookup['I0A10'] & lookup['R1'],
                'I0RA111': lookup['I0A11'] & lookup['R1'],
                'I0RA112': lookup['I0A12'] & lookup['R1'],
                'I0RA113': lookup['I0A13'] & lookup['R1'],
                'I0RA114': lookup['I0A14'] & lookup['R1'],
                   
                'I1RA000': lookup['I1A00'] & lookup['R0'],
                'I1RA001': lookup['I1A01'] & lookup['R0'],
                'I1RA002': lookup['I1A02'] & lookup['R0'],
                'I1RA003': lookup['I1A03'] & lookup['R0'],
                'I1RA004': lookup['I1A04'] & lookup['R0'],
                'I1RA010': lookup['I1A00'] & lookup['R1'],
                'I1RA011': lookup['I1A01'] & lookup['R1'],
                'I1RA012': lookup['I1A02'] & lookup['R1'],
                'I1RA013': lookup['I1A03'] & lookup['R1'],
                'I1RA014': lookup['I1A04'] & lookup['R1'],
                   
                'I1RA100': lookup['I1A10'] & lookup['R0'],
                'I1RA101': lookup['I1A11'] & lookup['R0'],
                'I1RA102': lookup['I1A12'] & lookup['R0'],
                'I1RA103': lookup['I1A13'] & lookup['R0'],
                'I1RA104': lookup['I1A14'] & lookup['R0'],
                'I1RA110': lookup['I1A10'] & lookup['R1'],
                'I1RA111': lookup['I1A11'] & lookup['R1'],
                'I1RA112': lookup['I1A12'] & lookup['R1'],
                'I1RA113': lookup['I1A13'] & lookup['R1'],
                'I1RA114': lookup['I1A14'] & lookup['R1'],
                   
                'I2RA000': lookup['I0A00'] & lookup['R0'],
                'I2RA001': lookup['I0A01'] & lookup['R0'],
                'I2RA002': lookup['I0A02'] & lookup['R0'],
                'I2RA003': lookup['I0A03'] & lookup['R0'],
                'I2RA004': lookup['I0A04'] & lookup['R0'],
                'I2RA010': lookup['I0A00'] & lookup['R1'],
                'I2RA011': lookup['I0A01'] & lookup['R1'],
                'I2RA012': lookup['I0A02'] & lookup['R1'],
                'I2RA013': lookup['I0A03'] & lookup['R1'],
                'I2RA014': lookup['I0A04'] & lookup['R1'],
                   
                'I2RA100': lookup['I2A10'] & lookup['R0'],
                'I2RA101': lookup['I2A11'] & lookup['R0'],
                'I2RA102': lookup['I2A12'] & lookup['R0'],
                'I2RA103': lookup['I2A13'] & lookup['R0'],
                'I2RA104': lookup['I2A14'] & lookup['R0'],
                'I2RA110': lookup['I2A10'] & lookup['R1'],
                'I2RA111': lookup['I2A11'] & lookup['R1'],
                'I2RA112': lookup['I2A12'] & lookup['R1'],
                'I2RA113': lookup['I2A13'] & lookup['R1'],
                'I2RA114': lookup['I2A14'] & lookup['R1'],
                   
                # All gender - age - partner combinations
                'PGA000': lookup['P0'] & lookup['GA00'],
                'PGA001': lookup['P0'] & lookup['GA01'],
                'PGA002': lookup['P0'] & lookup['GA02'],
                'PGA003': lookup['P0'] & lookup['GA03'],
                'PGA004': lookup['P0'] & lookup['GA04'],
                'PGA010': lookup['P0'] & lookup['GA10'],
                'PGA011': lookup['P0'] & lookup['GA11'],
                'PGA012': lookup['P0'] & lookup['GA12'],
                'PGA013': lookup['P0'] & lookup['GA13'],
                'PGA014': lookup['P0'] & lookup['GA14'],
                   
                'PLGA100': lookup['PLA10'] & lookup['G0'] & lookup['PC0'],
                'PLGA101': lookup['PLA11'] & lookup['G0'] & lookup['PC0'],
                'PLGA102': lookup['PLA12'] & lookup['G0'] & lookup['PC0'],
                'PLGA103': lookup['PLA13'] & lookup['G0'] & lookup['PC0'],
                'PLGA104': lookup['PLA14'] & lookup['G0'] & lookup['PC0'],
                'PLGA110': lookup['PLA10'] & lookup['G1'] & lookup['PC0'],
                'PLGA111': lookup['PLA11'] & lookup['G1'] & lookup['PC0'],
                'PLGA112': lookup['PLA12'] & lookup['G1'] & lookup['PC0'],
                'PLGA113': lookup['PLA13'] & lookup['G1'] & lookup['PC0'],
                'PLGA114': lookup['PLA14'] & lookup['G1'] & lookup['PC0'],
                
                'PSGA100': lookup['PSA10'] & lookup['G0'] & lookup['PC0'],
                'PSGA101': lookup['PSA11'] & lookup['G0'] & lookup['PC0'],
                'PSGA102': lookup['PSA12'] & lookup['G0'] & lookup['PC0'],
                'PSGA103': lookup['PSA13'] & lookup['G0'] & lookup['PC0'],
                'PSGA104': lookup['PSA14'] & lookup['G0'] & lookup['PC0'],
                'PSGA110': lookup['PSA10'] & lookup['G1'] & lookup['PC0'],
                'PSGA111': lookup['PSA11'] & lookup['G1'] & lookup['PC0'],
                'PSGA112': lookup['PSA12'] & lookup['G1'] & lookup['PC0'],
                'PSGA113': lookup['PSA13'] & lookup['G1'] & lookup['PC0'],
                'PSGA114': lookup['PSA14'] & lookup['G1'] & lookup['PC0'],
                
                'PCGA100': lookup['PC1'] & lookup['GA00'],
                'PCGA101': lookup['PC1'] & lookup['GA01'],
                'PCGA102': lookup['PC1'] & lookup['GA02'],
                'PCGA103': lookup['PC1'] & lookup['GA03'],
                'PCGA104': lookup['PC1'] & lookup['GA04'],
                'PCGA110': lookup['PC1'] & lookup['GA10'],
                'PCGA111': lookup['PC1'] & lookup['GA11'],
                'PCGA112': lookup['PC1'] & lookup['GA12'],
                'PCGA113': lookup['PC1'] & lookup['GA13'],
                'PCGA114': lookup['PC1'] & lookup['GA14'],
                
                # All risk - age - partner combinations
                'PRA000': lookup['P0'] & lookup['RA00'],
                'PRA001': lookup['P0'] & lookup['RA01'],
                'PRA002': lookup['P0'] & lookup['RA02'],
                'PRA003': lookup['P0'] & lookup['RA03'],
                'PRA004': lookup['P0'] & lookup['RA04'],
                'PRA010': lookup['P0'] & lookup['RA10'],
                'PRA011': lookup['P0'] & lookup['RA11'],
                'PRA012': lookup['P0'] & lookup['RA12'],
                'PRA013': lookup['P0'] & lookup['RA13'],
                'PRA014': lookup['P0'] & lookup['RA14'],
                
                'PLRA100': lookup['PLA10'] & lookup['R0'] & lookup['PC0'],
                'PLRA101': lookup['PLA11'] & lookup['R0'] & lookup['PC0'],
                'PLRA102': lookup['PLA12'] & lookup['R0'] & lookup['PC0'],
                'PLRA103': lookup['PLA13'] & lookup['R0'] & lookup['PC0'],
                'PLRA104': lookup['PLA14'] & lookup['R0'] & lookup['PC0'],
                'PLRA110': lookup['PLA10'] & lookup['R1'] & lookup['PC0'],
                'PLRA111': lookup['PLA11'] & lookup['R1'] & lookup['PC0'],
                'PLRA112': lookup['PLA12'] & lookup['R1'] & lookup['PC0'],
                'PLRA113': lookup['PLA13'] & lookup['R1'] & lookup['PC0'],
                'PLRA114': lookup['PLA14'] & lookup['R1'] & lookup['PC0'],
                
                'PSRA100': lookup['PSA10'] & lookup['R0'] & lookup['PC0'],
                'PSRA101': lookup['PSA11'] & lookup['R0'] & lookup['PC0'],
                'PSRA102': lookup['PSA12'] & lookup['R0'] & lookup['PC0'],
                'PSRA103': lookup['PSA13'] & lookup['R0'] & lookup['PC0'],
                'PSRA104': lookup['PSA14'] & lookup['R0'] & lookup['PC0'],
                'PSRA110': lookup['PSA10'] & lookup['R1'] & lookup['PC0'],
                'PSRA111': lookup['PSA11'] & lookup['R1'] & lookup['PC0'],
                'PSRA112': lookup['PSA12'] & lookup['R1'] & lookup['PC0'],
                'PSRA113': lookup['PSA13'] & lookup['R1'] & lookup['PC0'],
                'PSRA114': lookup['PSA14'] & lookup['R1'] & lookup['PC0'],
                
                'PCRA100': lookup['PC1'] & lookup['RA00'],
                'PCRA101': lookup['PC1'] & lookup['RA01'],
                'PCRA102': lookup['PC1'] & lookup['RA02'],
                'PCRA103': lookup['PC1'] & lookup['RA03'],
                'PCRA104': lookup['PC1'] & lookup['RA04'],
                'PCRA110': lookup['PC1'] & lookup['RA10'],
                'PCRA111': lookup['PC1'] & lookup['RA11'],
                'PCRA112': lookup['PC1'] & lookup['RA12'],
                'PCRA113': lookup['PC1'] & lookup['RA13'],
                'PCRA114': lookup['PC1'] & lookup['RA14'],

                # All age - risk - state combinations
                'SSRA00': lookup['SS'] & lookup['RA00'],
                'SSRA01': lookup['SS'] & lookup['RA01'],
                'SSRA02': lookup['SS'] & lookup['RA02'],
                'SSRA03': lookup['SS'] & lookup['RA03'],
                'SSRA04': lookup['SS'] & lookup['RA04'],
                'SSRA10': lookup['SS'] & lookup['RA10'],
                'SSRA11': lookup['SS'] & lookup['RA11'],
                'SSRA12': lookup['SS'] & lookup['RA12'],
                'SSRA13': lookup['SS'] & lookup['RA13'],
                'SSRA14': lookup['SS'] & lookup['RA14'],

                'SERA00': lookup['SE'] & lookup['RA00'],
                'SERA01': lookup['SE'] & lookup['RA01'],
                'SERA02': lookup['SE'] & lookup['RA02'],
                'SERA03': lookup['SE'] & lookup['RA03'],
                'SERA04': lookup['SE'] & lookup['RA04'],
                'SERA10': lookup['SE'] & lookup['RA10'],
                'SERA11': lookup['SE'] & lookup['RA11'],
                'SERA12': lookup['SE'] & lookup['RA12'],
                'SERA13': lookup['SE'] & lookup['RA13'],
                'SERA14': lookup['SE'] & lookup['RA14'],

                'SIRA00': lookup['SI'] & lookup['RA00'],
                'SIRA01': lookup['SI'] & lookup['RA01'],
                'SIRA02': lookup['SI'] & lookup['RA02'],
                'SIRA03': lookup['SI'] & lookup['RA03'],
                'SIRA04': lookup['SI'] & lookup['RA04'],
                'SIRA10': lookup['SI'] & lookup['RA10'],
                'SIRA11': lookup['SI'] & lookup['RA11'],
                'SIRA12': lookup['SI'] & lookup['RA12'],
                'SIRA13': lookup['SI'] & lookup['RA13'],
                'SIRA14': lookup['SI'] & lookup['RA14'],

                'STRA00': lookup['ST'] & lookup['RA00'],
                'STRA01': lookup['ST'] & lookup['RA01'],
                'STRA02': lookup['ST'] & lookup['RA02'],
                'STRA03': lookup['ST'] & lookup['RA03'],
                'STRA04': lookup['ST'] & lookup['RA04'],
                'STRA10': lookup['ST'] & lookup['RA10'],
                'STRA11': lookup['ST'] & lookup['RA11'],
                'STRA12': lookup['ST'] & lookup['RA12'],
                'STRA13': lookup['ST'] & lookup['RA13'],
                'STRA14': lookup['ST'] & lookup['RA14'],

                'SRRA00': lookup['SR'] & lookup['RA00'],
                'SRRA01': lookup['SR'] & lookup['RA01'],
                'SRRA02': lookup['SR'] & lookup['RA02'],
                'SRRA03': lookup['SR'] & lookup['RA03'],
                'SRRA04': lookup['SR'] & lookup['RA04'],
                'SRRA10': lookup['SR'] & lookup['RA10'],
                'SRRA11': lookup['SR'] & lookup['RA11'],
                'SRRA12': lookup['SR'] & lookup['RA12'],
                'SRRA13': lookup['SR'] & lookup['RA13'],
                'SRRA14': lookup['SR'] & lookup['RA14'],
                
                'SVRA00': lookup['SV'] & lookup['RA00'],
                'SVRA01': lookup['SV'] & lookup['RA01'],
                'SVRA02': lookup['SV'] & lookup['RA02'],
                'SVRA03': lookup['SV'] & lookup['RA03'],
                'SVRA04': lookup['SV'] & lookup['RA04'],
                'SVRA10': lookup['SV'] & lookup['RA10'],
                'SVRA11': lookup['SV'] & lookup['RA11'],
                'SVRA12': lookup['SV'] & lookup['RA12'],
                'SVRA13': lookup['SV'] & lookup['RA13'],
                'SVRA14': lookup['SV'] & lookup['RA14'],
                
                # All age - gender - state combinations
                'SSGA00': lookup['SS'] & lookup['GA00'],
                'SSGA01': lookup['SS'] & lookup['GA01'],
                'SSGA02': lookup['SS'] & lookup['GA02'],
                'SSGA03': lookup['SS'] & lookup['GA03'],
                'SSGA04': lookup['SS'] & lookup['GA04'],
                'SSGA10': lookup['SS'] & lookup['GA10'],
                'SSGA11': lookup['SS'] & lookup['GA11'],
                'SSGA12': lookup['SS'] & lookup['GA12'],
                'SSGA13': lookup['SS'] & lookup['GA13'],
                'SSGA14': lookup['SS'] & lookup['GA14'],

                'SEGA00': lookup['SE'] & lookup['GA00'],
                'SEGA01': lookup['SE'] & lookup['GA01'],
                'SEGA02': lookup['SE'] & lookup['GA02'],
                'SEGA03': lookup['SE'] & lookup['GA03'],
                'SEGA04': lookup['SE'] & lookup['GA04'],
                'SEGA10': lookup['SE'] & lookup['GA10'],
                'SEGA11': lookup['SE'] & lookup['GA11'],
                'SEGA12': lookup['SE'] & lookup['GA12'],
                'SEGA13': lookup['SE'] & lookup['GA13'],
                'SEGA14': lookup['SE'] & lookup['GA14'],

                'SIGA00': lookup['SI'] & lookup['GA00'],
                'SIGA01': lookup['SI'] & lookup['GA01'],
                'SIGA02': lookup['SI'] & lookup['GA02'],
                'SIGA03': lookup['SI'] & lookup['GA03'],
                'SIGA04': lookup['SI'] & lookup['GA04'],
                'SIGA10': lookup['SI'] & lookup['GA10'],
                'SIGA11': lookup['SI'] & lookup['GA11'],
                'SIGA12': lookup['SI'] & lookup['GA12'],
                'SIGA13': lookup['SI'] & lookup['GA13'],
                'SIGA14': lookup['SI'] & lookup['GA14'],

                'STGA00': lookup['ST'] & lookup['GA00'],
                'STGA01': lookup['ST'] & lookup['GA01'],
                'STGA02': lookup['ST'] & lookup['GA02'],
                'STGA03': lookup['ST'] & lookup['GA03'],
                'STGA04': lookup['ST'] & lookup['GA04'],
                'STGA10': lookup['ST'] & lookup['GA10'],
                'STGA11': lookup['ST'] & lookup['GA11'],
                'STGA12': lookup['ST'] & lookup['GA12'],
                'STGA13': lookup['ST'] & lookup['GA13'],
                'STGA14': lookup['ST'] & lookup['GA14'],

                'SRGA00': lookup['SR'] & lookup['GA00'],
                'SRGA01': lookup['SR'] & lookup['GA01'],
                'SRGA02': lookup['SR'] & lookup['GA02'],
                'SRGA03': lookup['SR'] & lookup['GA03'],
                'SRGA04': lookup['SR'] & lookup['GA04'],
                'SRGA10': lookup['SR'] & lookup['GA10'],
                'SRGA11': lookup['SR'] & lookup['GA11'],
                'SRGA12': lookup['SR'] & lookup['GA12'],
                'SRGA13': lookup['SR'] & lookup['GA13'],
                'SRGA14': lookup['SR'] & lookup['GA14'],
                
                'SVGA00': lookup['SV'] & lookup['GA00'],
                'SVGA01': lookup['SV'] & lookup['GA01'],
                'SVGA02': lookup['SV'] & lookup['GA02'],
                'SVGA03': lookup['SV'] & lookup['GA03'],
                'SVGA04': lookup['SV'] & lookup['GA04'],
                'SVGA10': lookup['SV'] & lookup['GA10'],
                'SVGA11': lookup['SV'] & lookup['GA11'],
                'SVGA12': lookup['SV'] & lookup['GA12'],
                'SVGA13': lookup['SV'] & lookup['GA13'],
                'SVGA14': lookup['SV'] & lookup['GA14'],
                })

    # Store in pop parameters
    pop_parameters.update({'lookup': lookup})


    return pop_parameters


#%% HELPER import_particular_age_group_and_gender()
#
#
# Helper for just importing people from one age group
#
#
def import_particular_age_group_and_gender(pop_parameters, inf_parameters, age_group, gender, n, t):


    # Set the desired age distribution
    temp = copy.deepcopy(pop_parameters)
    temp['age_dist'].loc[:, 'cdf'] = 0
    temp['age_dist'].loc[temp['age_dist'].age_upper == min(20 + 5*age_group, 36), 'cdf'] = 1


    # Set the desired gender distribution
    temp['sex_dist'].loc[:, 'pMale'] = gender


    # Generate new people for the population
    new = pop.generate_population(temp, n, inf_parameters['infection'].prob_import_infectious[0], 0)
    new.loc[:, 'import_time'] = t


    return new


#%% HELPER sample_from_imports()
#
#
# Select n individuals from imports dataset
#
#
def sample_from_imports(imports, age, gender, n_sample):


    # Select a few rows from the imports sample cohort at random
    cohort = str(gender) + str(age)
    temp = np.random.choice(imports[cohort].index, n_sample)
    temp = imports[cohort].loc[temp, :]
    temp = temp.reset_index(drop = True)


    return temp


#%% HELPER add_to_meta()
#
#
# Function which adds a new person into the meta dataframe and the
# associated matrices
#
#
def add_to_meta(pop_parameters, prt_parameters, meta, partner_matrix, partner_expire, new_person, t):


    # Put the new person into the meta-population
    meta = meta.append(new_person).reset_index(drop = True)
    pop_tot = len(meta)
    n_new = len(new_person)


    # Put the new person into the partner matrix
    partner_matrix = add_to_array(partner_matrix, pop_tot, 0)


    # Put the new person into the partner duration matrix
    partner_expire = add_to_array(partner_expire, pop_tot, float('inf'))


    # Test to see if they want a partner
    new_cases = new_person.index[(new_person.state == 'I') & (new_person.risk == 0)]
    if len(new_cases) > 0:
        
        
        # Work out which cases have a partner
        p_no_partner = pop_parameters['partners_dist'].to_numpy()[0,:]
        p_no_partner = p_no_partner[new_person.age_group.astype('int').values[new_cases]]
        has_partner = np.random.random(len(new_cases)) > p_no_partner
        has_partner = new_cases[has_partner]
        
        
        # Find the partners
        for i in has_partner:
            ii = pop_tot - (n_new-(i+1)) - 1
            jj = prt.find_partner(prt_parameters, meta, partner_matrix, ii)
    
    
            # Check somebody was found
            if jj != -1:
    
    
                # Sample a duration
                duration = prt.relationship_duration(prt_parameters, meta, ii, jj, 0)
    
    
                # Update partnership array
                partner_matrix[ii, jj] = 1
                partner_matrix[jj, ii] = 1
    
    
                # Update partnership duration matrix
                partner_expire[ii, jj] = t + duration
                partner_expire[jj, ii] = t + duration


    return meta, partner_matrix, partner_expire


#%% HELPER record_data_of_leavers()
#
#
#
#
#
def record_data_of_leavers(pop_parameters, meta, leaving, t):


    # Preallocate for storing people who left the population
    temp = copy.deepcopy(meta.iloc[leaving])
    temp['export_time'] = t
    pop_parameters.update({'leavers': pop_parameters['leavers'].append(temp).reset_index(drop = True)})


    return pop_parameters


#%% HELPER remove_from_meta()
def remove_from_meta(meta, partner_matrix, partner_expire, leave):


    # Take them out of the population
    # partner_matrix = np.delete(partner_matrix, leave, 0)
    # partner_matrix = np.delete(partner_matrix, leave, 1)
    # partner_expire = np.delete(partner_expire, leave, 0)
    # partner_expire = np.delete(partner_expire, leave, 1)
    partner_matrix = delete_from_array(partner_matrix, leave, len(meta), 0)
    partner_expire = delete_from_array(partner_expire, leave, len(meta), float('inf'))


    # Terminate long term relationships with these guys
    meta.loc[meta.partner.isin(leave), 'partner'] = -1


    # Take them out of meta
    meta = meta.drop(index = leave, axis = 0).reset_index()


    # Isolate how the indicies have changed
    old_index = meta['index']
    new_index = meta.index


    # Adjust the partnership indicies in meta
    old_partners = meta.partner[meta.partner > -1]
    for ii in old_partners.index:
        meta.loc[ii, 'partner'] = new_index[old_index == old_partners[ii]]
    meta = meta.drop('index', axis = 1)


    return meta, partner_matrix, partner_expire


#%% HELPER delete_from_array()
#
#
# Helper function to delete a row and column from an array
#
#
def delete_from_array(a, leave, n, values):


    # Drop those entries
    # a = a[keep, :][:, keep]


    # Remove values
    a[leave, :] = values
    a[:, leave] = values


    # Shuffle unwanted rows and columns to the end of the array but don't delete them
    keep = list(set(range(0, n)) - set(leave))
    permutation = np.array(keep + leave)
    a = a[permutation, :][:, permutation]


    return a


def add_to_array(a, tot_needed, values):


    # Check first that there isn't already an empty row and column lying around
    to_add = tot_needed - len(a)


    # If more rows are needed then do a vstack/hstack
    if to_add > 0:


        # Put the new person into the partner matrix
        new_rows = values * np.ones((to_add, len(a)))
        new_cols = values * np.ones((tot_needed, to_add))

        a = np.vstack((a, new_rows))
        a = np.hstack((a, new_cols))


    return a


#%% GRAPHS initilise_demographic_trackers()
#
#
#
#
#
def initilise_demographic_trackers(n_steps, t0_sim = 0):
    n_years = int(np.ceil(n_steps/365)) + int(t0_sim/365)
    compartments = np.zeros((n_steps, 10))
    import_count = np.zeros((n_steps, 10))
    export_count = np.zeros((n_steps, 10))
    demographic = np.zeros((n_steps, 5*10))
    import_time = np.zeros((n_steps, (n_years + 2)*10))
    active_age = np.zeros((n_steps, 5*3))
    return compartments, import_count, export_count, demographic, import_time, active_age


#%% GRAPHS update_demographic_tracker()
#
#
#
#
#
def update_demographic_tracker(t, pop_parameters, meta, compartments, import_count, demographic, import_time, export_count, active_age, t0_sim = 0):


    # How many are in each compartment?
    compartments[t, 0] = sum(pop_parameters['lookup']['GA00'])
    compartments[t, 1] = sum(pop_parameters['lookup']['GA01'])
    compartments[t, 2] = sum(pop_parameters['lookup']['GA02'])
    compartments[t, 3] = sum(pop_parameters['lookup']['GA03'])
    compartments[t, 4] = sum(pop_parameters['lookup']['GA04'])
    compartments[t, 5] = sum(pop_parameters['lookup']['GA10'])
    compartments[t, 6] = sum(pop_parameters['lookup']['GA11'])
    compartments[t, 7] = sum(pop_parameters['lookup']['GA12'])
    compartments[t, 8] = sum(pop_parameters['lookup']['GA13'])
    compartments[t, 9] = sum(pop_parameters['lookup']['GA14'])


    # Specifically, how many have been imported?
    new_import = meta.import_time == t + t0_sim
    import_count[t, 0] = sum((pop_parameters['lookup']['GA00']) & new_import)
    import_count[t, 1] = sum((pop_parameters['lookup']['GA01']) & new_import)
    import_count[t, 2] = sum((pop_parameters['lookup']['GA02']) & new_import)
    import_count[t, 3] = sum((pop_parameters['lookup']['GA03']) & new_import)
    import_count[t, 4] = sum((pop_parameters['lookup']['GA04']) & new_import)
    import_count[t, 5] = sum((pop_parameters['lookup']['GA10']) & new_import)
    import_count[t, 6] = sum((pop_parameters['lookup']['GA11']) & new_import)
    import_count[t, 7] = sum((pop_parameters['lookup']['GA12']) & new_import)
    import_count[t, 8] = sum((pop_parameters['lookup']['GA13']) & new_import)
    import_count[t, 9] = sum((pop_parameters['lookup']['GA14']) & new_import)


    # How many are being exported
    new_export = pop_parameters['leavers'].export_time == t + t0_sim
    export_count[t, 0] = sum(new_export & (pop_parameters['leavers'].gender == 0) & (pop_parameters['leavers'].age_group == 0))
    export_count[t, 1] = sum(new_export & (pop_parameters['leavers'].gender == 0) & (pop_parameters['leavers'].age_group == 1))
    export_count[t, 2] = sum(new_export & (pop_parameters['leavers'].gender == 0) & (pop_parameters['leavers'].age_group == 2))
    export_count[t, 3] = sum(new_export & (pop_parameters['leavers'].gender == 0) & (pop_parameters['leavers'].age_group == 3))
    export_count[t, 4] = sum(new_export & (pop_parameters['leavers'].gender == 0) & (pop_parameters['leavers'].age_group == 4))
    export_count[t, 5] = sum(new_export & (pop_parameters['leavers'].gender == 1) & (pop_parameters['leavers'].age_group == 0))
    export_count[t, 6] = sum(new_export & (pop_parameters['leavers'].gender == 1) & (pop_parameters['leavers'].age_group == 1))
    export_count[t, 7] = sum(new_export & (pop_parameters['leavers'].gender == 1) & (pop_parameters['leavers'].age_group == 2))
    export_count[t, 8] = sum(new_export & (pop_parameters['leavers'].gender == 1) & (pop_parameters['leavers'].age_group == 3))
    export_count[t, 9] = sum(new_export & (pop_parameters['leavers'].gender == 1) & (pop_parameters['leavers'].age_group == 4))


    # Update the orientation and risk by cohort trackers
    ii = 0
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            denom = sum(pop_parameters['lookup']['GA' + str(gender) + str(age)])
            if denom > 0:
                demographic[t, ii + 0] = sum(pop_parameters['lookup']['GAR' + str(gender) + str(age) + '0'])/denom
                demographic[t, ii + 1] = sum(pop_parameters['lookup']['GAR' + str(gender) + str(age) + '1'])/denom
                demographic[t, ii + 2] = sum(pop_parameters['lookup']['GAO' + str(gender) + str(age) + '0'])/denom
                demographic[t, ii + 3] = sum(pop_parameters['lookup']['GAO' + str(gender) + str(age) + '1'])/denom
                demographic[t, ii + 4] = sum(pop_parameters['lookup']['GAO' + str(gender) + str(age) + '2'])/denom
            ii = ii + 5


    # Update the proportion of individuals by the year that they entered the simulation
    ii = 0
    _, ncol = import_time.shape
    n_years = int(ncol/10) - 2
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            denom = sum(pop_parameters['lookup']['GA' + str(gender) + str(age)])
            if denom > 0:
                import_time[t, ii + 0] = sum((meta.import_time == 0) & pop_parameters['lookup']['GA' + str(gender) + str(age)]) / denom
                import_time[t, ii + 1] = sum((np.floor(meta.import_time/365) == 0) & (meta.import_time != 0) & pop_parameters['lookup']['GA' + str(gender) + str(age)]) / denom
                for year in range(1, n_years):
                    import_time[t, ii + 1 + year] = sum((np.floor(meta.import_time/365) == year) & pop_parameters['lookup']['GA' + str(gender) + str(age)]) / denom
            ii = ii + n_years + 2
                    
                    
    # Track the number of sexually active poeple
    ii = 0
    for a in [0, 1, 2, 3, 4]:
        who = meta.loc[(meta.age_group == a) & (meta.gender == 1), :]
        n = len(who)
        if n > 0:
            active_age[t, ii + 0] = np.sum(who.has_oral)/n
            active_age[t, ii + 1] = np.sum(who.has_sex)/n
            active_age[t, ii + 2] = np.sum(who.has_anal)/n
        ii = ii + 3


    return compartments, import_count, demographic, import_time, export_count, active_age


#%% GRAPHS make_demographic_graphs()
#
#
#
#
#
def make_demographic_graphs(tt, pop_parameters, compartments, import_count, demographic, import_time, export_count, active_age, save_loc = 'graphs/output/'):


    ## GRAPH OF AGE AND GENDER COHORT SIZES OVER TIME

    # Initilise figure
    fig, axs = plt.subplots(5, 2)


    # Number in the compartment
    axs[0, 0].plot(tt, compartments[:,0])
    axs[1, 0].plot(tt, compartments[:,1])
    axs[2, 0].plot(tt, compartments[:,2])
    axs[3, 0].plot(tt, compartments[:,3])
    axs[4, 0].plot(tt, compartments[:,4])

    axs[0, 1].plot(tt, compartments[:,5])
    axs[1, 1].plot(tt, compartments[:,6])
    axs[2, 1].plot(tt, compartments[:,7])
    axs[3, 1].plot(tt, compartments[:,8])
    axs[4, 1].plot(tt, compartments[:,9])


    # Target number for compartments
    targetF = pop_parameters['age_dist'].ave[pop_parameters['age_dist'].sex == 0].reset_index(drop=True) * (1-pop_parameters['sex_dist']['pMale'].iloc[0]) * pop_parameters['n']
    targetM = pop_parameters['age_dist'].ave[pop_parameters['age_dist'].sex == 1].reset_index(drop=True) * pop_parameters['sex_dist']['pMale'].iloc[0] * pop_parameters['n']
    axs[0, 0].plot(tt, ([targetF[0]] * len(tt)), linestyle = '--')
    axs[1, 0].plot(tt, ([targetF[1]] * len(tt)), linestyle = '--')
    axs[2, 0].plot(tt, ([targetF[2]] * len(tt)), linestyle = '--')
    axs[3, 0].plot(tt, ([targetF[3]] * len(tt)), linestyle = '--')
    axs[4, 0].plot(tt, ([targetF[4]] * len(tt)), linestyle = '--')

    axs[0, 1].plot(tt, ([targetM[0]] * len(tt)), linestyle = '--')
    axs[1, 1].plot(tt, ([targetM[1]] * len(tt)), linestyle = '--')
    axs[2, 1].plot(tt, ([targetM[2]] * len(tt)), linestyle = '--')
    axs[3, 1].plot(tt, ([targetM[3]] * len(tt)), linestyle = '--')
    axs[4, 1].plot(tt, ([targetM[4]] * len(tt)), linestyle = '--')


    # Axis labels
    plt.suptitle('Demographic Sizes Over Time')
    axs[0, 0].set_title('Females')
    axs[0, 1].set_title('Males')
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')

    
    # Save the graph
    fig.savefig(save_loc + 'demographic_by_age.png', dpi = 200)
    plt.close(fig)


    ## GRAPH OF THE NUMBER OF INDIVIDUALS IMPORTED AT EACH TIME STEP


    # Initilise figure
    fig, axs = plt.subplots(5, 2)


    # Number of imports
    import_count[0,:] = 0
    axs[0, 0].plot(tt, import_count[:,0])
    axs[1, 0].plot(tt, import_count[:,1])
    axs[2, 0].plot(tt, import_count[:,2])
    axs[3, 0].plot(tt, import_count[:,3])
    axs[4, 0].plot(tt, import_count[:,4])

    axs[0, 1].plot(tt, import_count[:,5])
    axs[1, 1].plot(tt, import_count[:,6])
    axs[2, 1].plot(tt, import_count[:,7])
    axs[3, 1].plot(tt, import_count[:,8])
    axs[4, 1].plot(tt, import_count[:,9])


    # Mean import_count
    import_count[0,:] = 0
    axs[0, 0].plot(tt, len(tt) * [np.mean(import_count[:,0])], linestyle = '--')
    axs[1, 0].plot(tt, len(tt) * [np.mean(import_count[:,1])], linestyle = '--')
    axs[2, 0].plot(tt, len(tt) * [np.mean(import_count[:,2])], linestyle = '--')
    axs[3, 0].plot(tt, len(tt) * [np.mean(import_count[:,3])], linestyle = '--')
    axs[4, 0].plot(tt, len(tt) * [np.mean(import_count[:,4])], linestyle = '--')

    axs[0, 1].plot(tt, len(tt) * [np.mean(import_count[:,5])], linestyle = '--')
    axs[1, 1].plot(tt, len(tt) * [np.mean(import_count[:,6])], linestyle = '--')
    axs[2, 1].plot(tt, len(tt) * [np.mean(import_count[:,7])], linestyle = '--')
    axs[3, 1].plot(tt, len(tt) * [np.mean(import_count[:,8])], linestyle = '--')
    axs[4, 1].plot(tt, len(tt) * [np.mean(import_count[:,9])], linestyle = '--')


    # Number of exports
    export_count[0,:] = 0
    axs[0, 0].plot(tt, -export_count[:,0])
    axs[1, 0].plot(tt, -export_count[:,1])
    axs[2, 0].plot(tt, -export_count[:,2])
    axs[3, 0].plot(tt, -export_count[:,3])
    axs[4, 0].plot(tt, -export_count[:,4])

    axs[0, 1].plot(tt, -export_count[:,5])
    axs[1, 1].plot(tt, -export_count[:,6])
    axs[2, 1].plot(tt, -export_count[:,7])
    axs[3, 1].plot(tt, -export_count[:,8])
    axs[4, 1].plot(tt, -export_count[:,9])


    # Mean export_count
    export_count[0,:] = 0
    axs[0, 0].plot(tt, len(tt) * [np.mean(-export_count[:,0])], linestyle = '--')
    axs[1, 0].plot(tt, len(tt) * [np.mean(-export_count[:,1])], linestyle = '--')
    axs[2, 0].plot(tt, len(tt) * [np.mean(-export_count[:,2])], linestyle = '--')
    axs[3, 0].plot(tt, len(tt) * [np.mean(-export_count[:,3])], linestyle = '--')
    axs[4, 0].plot(tt, len(tt) * [np.mean(-export_count[:,4])], linestyle = '--')

    axs[0, 1].plot(tt, len(tt) * [np.mean(-export_count[:,5])], linestyle = '--')
    axs[1, 1].plot(tt, len(tt) * [np.mean(-export_count[:,6])], linestyle = '--')
    axs[2, 1].plot(tt, len(tt) * [np.mean(-export_count[:,7])], linestyle = '--')
    axs[3, 1].plot(tt, len(tt) * [np.mean(-export_count[:,8])], linestyle = '--')
    axs[4, 1].plot(tt, len(tt) * [np.mean(-export_count[:,9])], linestyle = '--')


    # Axis labels
    plt.suptitle('Number of Imports and Exports')
    axs[0, 0].set_title('Females')
    axs[0, 1].set_title('Males')
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')
    
    
    # Save the graph
    fig.savefig(save_loc + 'demographic_import_export_age_gender.png', dpi = 200)
    plt.close(fig)


    ## GRAPH OF THE NET CHANGE IN COHORT SIZES AFTER TAKING OUT IMPORTS


    # Initilise figure
    fig, axs = plt.subplots(5, 2)


    # Net change in population size
    axs[0, 0].plot(tt[1:len(tt)], (compartments[1:len(tt),0] - compartments[0:(len(tt)-1),0] - import_count[1:len(tt),0] + export_count[1:len(tt),0]))
    axs[1, 0].plot(tt[1:len(tt)], (compartments[1:len(tt),1] - compartments[0:(len(tt)-1),1] - import_count[1:len(tt),1] + export_count[1:len(tt),1]))
    axs[2, 0].plot(tt[1:len(tt)], (compartments[1:len(tt),2] - compartments[0:(len(tt)-1),2] - import_count[1:len(tt),2] + export_count[1:len(tt),2]))
    axs[3, 0].plot(tt[1:len(tt)], (compartments[1:len(tt),3] - compartments[0:(len(tt)-1),3] - import_count[1:len(tt),3] + export_count[1:len(tt),3]))
    axs[4, 0].plot(tt[1:len(tt)], (compartments[1:len(tt),4] - compartments[0:(len(tt)-1),4] - import_count[1:len(tt),4] + export_count[1:len(tt),4]))

    axs[0, 1].plot(tt[1:len(tt)], (compartments[1:len(tt),5] - compartments[0:(len(tt)-1),5] - import_count[1:len(tt),5] + export_count[1:len(tt),5]))
    axs[1, 1].plot(tt[1:len(tt)], (compartments[1:len(tt),6] - compartments[0:(len(tt)-1),6] - import_count[1:len(tt),6] + export_count[1:len(tt),6]))
    axs[2, 1].plot(tt[1:len(tt)], (compartments[1:len(tt),7] - compartments[0:(len(tt)-1),7] - import_count[1:len(tt),7] + export_count[1:len(tt),7]))
    axs[3, 1].plot(tt[1:len(tt)], (compartments[1:len(tt),8] - compartments[0:(len(tt)-1),8] - import_count[1:len(tt),8] + export_count[1:len(tt),8]))
    axs[4, 1].plot(tt[1:len(tt)], (compartments[1:len(tt),9] - compartments[0:(len(tt)-1),9] - import_count[1:len(tt),9] + export_count[1:len(tt),9]))


    # Mean net change
    axs[0, 0].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),0] - compartments[0:(len(tt)-1),0] - import_count[1:len(tt),0] + export_count[1:len(tt),0]))], linestyle = '--')
    axs[1, 0].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),1] - compartments[0:(len(tt)-1),1] - import_count[1:len(tt),1] + export_count[1:len(tt),1]))], linestyle = '--')
    axs[2, 0].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),2] - compartments[0:(len(tt)-1),2] - import_count[1:len(tt),2] + export_count[1:len(tt),2]))], linestyle = '--')
    axs[3, 0].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),3] - compartments[0:(len(tt)-1),3] - import_count[1:len(tt),3] + export_count[1:len(tt),3]))], linestyle = '--')
    axs[4, 0].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),4] - compartments[0:(len(tt)-1),4] - import_count[1:len(tt),4] + export_count[1:len(tt),4]))], linestyle = '--')

    axs[0, 1].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),5] - compartments[0:(len(tt)-1),5] - import_count[1:len(tt),5] + export_count[1:len(tt),5]))], linestyle = '--')
    axs[1, 1].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),6] - compartments[0:(len(tt)-1),6] - import_count[1:len(tt),6] + export_count[1:len(tt),6]))], linestyle = '--')
    axs[2, 1].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),7] - compartments[0:(len(tt)-1),7] - import_count[1:len(tt),7] + export_count[1:len(tt),7]))], linestyle = '--')
    axs[3, 1].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),8] - compartments[0:(len(tt)-1),8] - import_count[1:len(tt),8] + export_count[1:len(tt),8]))], linestyle = '--')
    axs[4, 1].plot(tt, len(tt) * [np.mean((compartments[1:len(tt),9] - compartments[0:(len(tt)-1),9] - import_count[1:len(tt),9] + export_count[1:len(tt),9]))], linestyle = '--')


    # Axis labels
    plt.suptitle('Change in Age Groups (Age Dynamics)')
    axs[0, 0].set_title('Females')
    axs[0, 1].set_title('Males')
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')

    
    # Save the graph
    fig.savefig(save_loc + 'demographic_ageing_age_gender.png', dpi = 200)
    plt.close(fig)


    ## GRAPH OF THE BASIS POINT DIFFERENCE BETWEEN THE TARGET AND ACTUAL
    ## ORIENTATION AND RISK COHORT SIZES


    # Initilise figure
    fig, axs = plt.subplots(5, 2)


    # Net change in population size
    ii = 0
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            homo = pop_parameters['orientation_dist'].homo - pop_parameters['orientation_dist'].hetero
            bi = pop_parameters['orientation_dist'].bi - pop_parameters['orientation_dist'].homo
            axs[age, gender].plot(tt, pop_parameters['orientation_dist'].hetero[age] - demographic[:, ii + 2], label = 'Hetero-sexual')
            axs[age, gender].plot(tt, (1-pop_parameters['partners_dist'].iloc[3][age]) - demographic[:, ii + 0], label = 'Low-Risk', linestyle = '--')
            axs[age, gender].plot(tt, homo[age] - demographic[:, ii + 3], label = 'Homo-sexual')
            axs[age, gender].plot(tt, pop_parameters['partners_dist'].iloc[3][age] - demographic[:, ii + 1], label = 'High-Risk', linestyle = '--')
            axs[age, gender].plot(tt, bi[age] - demographic[:, ii + 4], label = 'Bi-sexual')
            axs[age, gender].plot([min(tt), max(tt)], [0, 0], color = 'black', linestyle = ':', label = '')
            ii = ii + 5


    # Axis labels
    plt.suptitle('Difference Between Target Proportions and Actual Proportions\nFor Orientation and Risk Cohorts')
    axs[0, 0].set_title('Females')
    axs[0, 1].set_title('Males')
    axs[4, 0].legend(loc = 'upper center', ncol = 3, bbox_to_anchor = (1.1, -0.05), fancybox = True)
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')
    
    
    # Save the graph
    fig.savefig(save_loc + 'demographic_deviation_from_target_orientations_age_gender.png', dpi = 200)
    plt.close(fig)


    ## GRAPH SHOWING THE PROPORTION OF INDIVIDUALS IN THE SIMULATION BY
    ## THE YEAR THAT THEY ENTERED THE SIMULATION

    # Initilise figure
    fig, axs = plt.subplots(5, 2)


    # Net change in population size
    _, ncol = import_time.shape
    n_years = int(ncol/10) - 2
    ii = 0
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            for year in range(0, n_years + 1):
                axs[age, gender].plot(tt, import_time[:, ii + year])
            ii = ii + n_years + 2


    # Axis labels
    plt.suptitle('Proportion of Individuals by the Year that they Entered the Population\nFor Each Gender and Age Group Cohort')
    axs[0, 0].set_title('Females')
    axs[0, 1].set_title('Males')
    # axs[4, 0].legend(loc = 'upper center', ncol = 3, bbox_to_anchor = (1.1, -0.05), fancybox = True)
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')
    
    
    # Save the graph
    fig.savefig(save_loc + 'demographic_population_by_year_of_arrival_age_gender.png', dpi = 200)
    plt.close(fig)
    
    
    ##  GRAPH SHOWING POPULATION TURNOVER
    
    
    # Some function for computing a rolling sum
    def rolling_sum(a, n=4) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret
    
    
    # Initilise figure
    fig, axs = plt.subplots(5, 2)


    # Net change in population size
    ii = 0
    col = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            
            # Plot imports
            if age > 0:
                target = 100 * pop_parameters['pop_turnover']['annual_arrivals'][0]
                imports = 100 * rolling_sum(import_count[:, ii], 365) / pop_parameters['target'][str(gender)][age]
                axs[age, gender].plot(tt, imports, color = col[0])
                axs[age, gender].plot(tt, [target] * len(tt), color = col[0], linestyle = '--')
            
            
            # Plot exports
            if age < 4:
                target = 100 * pop_parameters['pop_turnover']['annual_departures'][0]
                exports = 100 * rolling_sum(export_count[:, ii], 365) / pop_parameters['target'][str(gender)][age]
                axs[age, gender].plot(tt, -exports, color = col[1])
                axs[age, gender].plot(tt, [-target] * len(tt), color = col[1], linestyle = '--')
            
            
            ii = ii + 1
        
    
    # Axis labels
    plt.suptitle('Population Turnover\nRolling 1-year Turnover Rate Compared to Target Rate')
    axs[0, 0].set_title('Females')
    axs[0, 1].set_title('Males')
    # axs[4, 0].legend(loc = 'upper center', ncol = 3, bbox_to_anchor = (1.1, -0.05), fancybox = True)
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')
    
    
    # Save the graph
    fig.savefig(save_loc + 'demographic_population_turnover.png', dpi = 200)
    plt.close(fig)
    
    
    
    
    ##  SEXUALLY ACTIVE BY AGE - COMPARISON TO GOANNA
    
    
    # Set colour cycle
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    
    
    # Make plot trackng sexual acts by age
    fig, ax = plt.subplots(5, 1)
    teval = np.array(tt) - tt[0]
    ii = 0
    for a in [0, 1, 2, 3, 4]:
        ax[a].plot(teval, active_age[:, ii + 0] - pop_parameters['sexual_acts_dist'].loc[0][1+np.min([2, a])], label = 'Oral')
        ax[a].plot(teval, active_age[:, ii + 1] - pop_parameters['sexual_acts_dist'].loc[1][1+np.min([2, a])], label = 'Vaginal')
        ax[a].plot(teval, active_age[:, ii + 2] - pop_parameters['sexual_acts_dist'].loc[2][1+np.min([2, a])], label = 'Anal')
        ax[a].plot(teval, len(tt) * [0], linestyle = '--', color = 'black')    
        ii = ii + 3
    
    
    # Label the plot
    ax[4].set_xlabel('Day')
    ax[0].set_ylabel('16-19')
    ax[1].set_ylabel('20-24')
    ax[2].set_ylabel('25-29')
    ax[3].set_ylabel('30-34')
    ax[4].set_ylabel('35')
    ax[0].set_title('Proportion of the Population Engaging in Each Sexual Act - GOANNA Proportions')
    # plt.suptitle('Comparison of Modelled Sexual Act Rates to GOANNA Rates')
    ax[0].legend()
    
    
    # Save the graph
    fig.savefig(save_loc + 'demographic_engages_in_sexual_act.png', dpi = 200)
    plt.close(fig)
    
    
    
    
    
    
    


