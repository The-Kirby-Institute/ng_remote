#  -*- coding: utf-8 -*-
"""
Created on Mon May 17 09:33:12 2021

@author: nicol
"""


# Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
from glob import glob
import os
import re
from multiprocessing import Pool


# Model parameters
param = pd.read_csv('data/param.csv')
tt = range(0, param['simulation_length'][0])


#%% FUN extract_prevalence()
#
#
# Just runs through the prevalence data and extracts that for further analysis
#
#
def extract_prevalence(sim):


    # Try all three scenarios
    for scenario in [1, 2, 3]:


        # print update
        print('PREVALENCE - Scenario ' + str(scenario) + ' and parameter set ' + str(sim), flush = True)


        # Check that the requested simulation has been run
        if check_for_output_data(scenario, sim) == True:


            # Check to see that the prevalence hasn't already been calculated
            output_name = 'simulations/prevalence/scenario_' + str(scenario) + '/data' + str(sim) + '.npy'
            if os.path.exists(output_name) == False:


                # Preallocate
                prvt = np.zeros((param['simulation_length'][0], 13))


                # Iterate over all time steps
                for t in range(0, param['simulation_length'][0], 50):


                    # Read in meta
                    meta = pd.read_feather('simulations/calibration/scenario_' + str(scenario) + '/simulation_' + str(sim) + '/timestep' + str(param['partner_burn_in'][0] + t) + '.ftr')


                    # Compute prevalence
                    # I = meta.state == 'I'
                    I = (meta.site1 == 1)
                    prvt[t,:] = [100 * sum(I)/len(meta),
                                100 * sum(I & (meta.gender == 1))/sum(meta.gender == 1),
                                100 * sum(I & (meta.gender == 0))/sum(meta.gender == 0),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 0))/max(1, sum( (meta.gender == 1) & (meta.age_group == 0) ) ),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 1))/max(1, sum( (meta.gender == 1) & (meta.age_group == 1) ) ),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 2))/max(1, sum( (meta.gender == 1) & (meta.age_group == 2) ) ),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 3) & (meta.age < 35))/max(1, sum( (meta.gender == 1) & (meta.age_group == 3) & (meta.age < 35) )),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 3) & (meta.age > 35))/max(1, sum( (meta.gender == 1) & (meta.age_group == 3) & (meta.age > 35) )),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 0))/max(1, sum( (meta.gender == 0) & (meta.age_group == 0) ) ),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 1))/max(1, sum( (meta.gender == 0) & (meta.age_group == 1) ) ),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 2))/max(1, sum( (meta.gender == 0) & (meta.age_group == 2) ) ),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 3) & (meta.age < 35))/max(1, sum( (meta.gender == 0) & (meta.age_group == 3) & (meta.age < 35) )),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 3) & (meta.age > 35))/max(1, sum( (meta.gender == 0) & (meta.age_group == 3) & (meta.age > 35) ))]


                # Save data
                np.save(output_name, prvt)



#%% FUN run_standard_checks()
#
#
# Makes a bunch of graphs to check that the simulation ran as expected
#
#
def run_standard_checks(sim):


    # Try all three scenarios
    for scenario in [1, 2, 3]:


        # print update
        print('CHECKS - Scenario ' + str(scenario) + ' and parameter set ' + str(sim), flush = True)


        # Check that this one hasn't already been done
        if os.path.exists('simulations/checks/scenario_' + str(scenario) + '/simulation_' + str(sim) + '_prevalence_age_group.png') == False:


            # Check that the requested set has been run
            if check_for_output_data(scenario, sim) == True:


                # Preallocate for output
                file_name = 'simulations/calibration/scenario_' + str(scenario) + '/simulation_'
                yt = np.zeros((param['simulation_length'][0], 5))
                popt = np.zeros((param['simulation_length'][0], 4))
                part = np.zeros((param['simulation_length'][0], 2))
                inft = np.zeros((param['simulation_length'][0], 8))
                prvt = np.zeros((param['simulation_length'][0], 13))


                # Iterate over each time step
                for t in tt:


                    # Read in file
                    meta = pd.read_feather(file_name + str(sim) + '/timestep' + str(param['partner_burn_in'][0] + t) + '.ftr')


                    # Compute the number in each infectious state
                    yt[t,:] = [sum(meta.state == 'S'),
                               sum(meta.state == 'E'),
                               sum(meta.state == 'I'),
                               sum(meta.state == 'R'),
                               sum(meta.state == 'T')]


                    # Compute the number in each age group
                    popt[t,:] = [sum(meta.age_group == 0),
                                sum(meta.age_group == 1),
                                sum(meta.age_group == 2),
                                sum(meta.age_group == 3)]


                    # Compute the total number of long-term relationships
                    part[t,:] = [sum(meta.partner == -1),
                                sum(meta.partner > -1)]


                    # Compute the total nuber of infections by site
                    site0 = meta.site0 == 1
                    site1 = meta.site1 == 1
                    site2 = meta.site2 == 1
                    inft[t,:] = [sum( ~site0 & ~site1 & ~site2 ),
                                sum( site0 & ~site1 & ~site2 ),
                                sum( ~site0 & site1 & ~site2 ),
                                sum( ~site0 & ~site1 & site2 ),
                                sum( site0 & site1 & ~site2 ),
                                sum( ~site0 & site1 & site2 ),
                                sum( site0 & ~site1 & site2 ),
                                sum( site0 & site1 & site2 )]


                    # Compute prevalence by sex and age group
                    I = meta.state == 'I'
                    prvt[t,:] = [100 * sum(I)/len(meta),
                                100 * sum(I & (meta.gender == 1))/sum(meta.gender == 1),
                                100 * sum(I & (meta.gender == 0))/sum(meta.gender == 0),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 0))/max(1, sum( (meta.gender == 1) & (meta.age_group == 0) ) ),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 1))/max(1, sum( (meta.gender == 1) & (meta.age_group == 1) ) ),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 2))/max(1, sum( (meta.gender == 1) & (meta.age_group == 2) ) ),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 3) & (meta.age < 35))/max(1, sum( (meta.gender == 1) & (meta.age_group == 3) & (meta.age < 35) )),
                                100 * sum(I & (meta.gender == 1) & (meta.age_group == 3) & (meta.age > 35))/max(1, sum( (meta.gender == 1) & (meta.age_group == 3) & (meta.age > 35) )),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 0))/max(1, sum( (meta.gender == 0) & (meta.age_group == 0) ) ),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 1))/max(1, sum( (meta.gender == 0) & (meta.age_group == 1) ) ),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 2))/max(1, sum( (meta.gender == 0) & (meta.age_group == 2) ) ),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 3) & (meta.age < 35))/max(1, sum( (meta.gender == 0) & (meta.age_group == 3) & (meta.age < 35) )),
                                100 * sum(I & (meta.gender == 0) & (meta.age_group == 3) & (meta.age > 35))/max(1, sum( (meta.gender == 0) & (meta.age_group == 3) & (meta.age > 35) ))]


                # Make graphs
                make_graphs(scenario, sim, yt, popt, part, inft, prvt)


#%% FUN make_graphs()
def make_graphs(scenario, sim, yt, popt, part, inft, prvt):


    # Where to save output data
    file_name = 'simulations/checks/scenario_' + str(scenario) + '/simulation_' + str(sim)


    # Plot aggregate infection levels
    plt.plot(tt, yt[:,0], label = 'S')
    plt.plot(tt, yt[:,1], label = 'E')
    plt.plot(tt, yt[:,2], label = 'I')
    plt.plot(tt, yt[:,3], label = 'R')
    plt.plot(tt, yt[:,4], label = 'T')
    plt.title('Aggregate Infectious State - Parameter Set ' + str(sim))
    plt.legend()
    plt.savefig(file_name + '_aggregate_infections.png')
    plt.close()


    # Plot the number in each age group
    plt.plot(tt, popt[:,0], label = '16-19')
    plt.plot(tt, popt[:,1], label = '20-24')
    plt.plot(tt, popt[:,2], label = '25-30')
    plt.plot(tt, popt[:,3], label = 'Over 30')
    plt.plot(tt, np.sum(popt, axis = 1), label = 'Total')
    plt.title('Age Group Occupancy - Parameter Set ' + str(sim))
    plt.legend()
    plt.savefig(file_name + '_age_group.png')
    plt.close()


    # Plot the number in a long term relationship
    plt.plot(tt, part[:,0], label = 'Single')
    plt.plot(tt, popt[:,1], label = 'Long-term relationship')
    plt.title('Long-term Relationships - Parameter Set ' + str(sim))
    plt.legend()
    plt.savefig(file_name + '_long_term_relationships.png')
    plt.close()


    # Graph of infections by anatomical site
    # Setup graph
    fig, axs = plt.subplots(3, 1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Aggregate Infections by Anatomical-site - Parameter Set ' + str(sim))


    # make first pannel
    axs[0].plot(tt, inft[:,1] + inft[:,4] + inft[:,6] + inft[:,7], label = 'Rectal')
    axs[0].plot(tt, inft[:,2] + inft[:,4] + inft[:,5] + inft[:,7], label = 'Urethral')
    axs[0].plot(tt, inft[:,3] + inft[:,5] + inft[:,6] + inft[:,7], label = 'Pharyngeal')
    axs[0].legend()
    axs[0].set_title('Infections Including Each Site')


    # make second pannel
    axs[1].plot(tt, inft[:,1], label = 'Rectal')
    axs[1].plot(tt, inft[:,2], label = 'Urethral')
    axs[1].plot(tt, inft[:,3], label = 'Pharyngeal')
    axs[1].legend()
    axs[1].set_title('Infections at Just One Site')


    # Make third pannel
    axs[2].plot(tt, inft[:,4], label = 'Sites Rec and Ure')
    axs[2].plot(tt, inft[:,5], label = 'Sites Ure and Pha')
    axs[2].plot(tt, inft[:,6], label = 'Sites Rec and Pha')
    axs[2].plot(tt, inft[:,7], label = 'All sites')
    axs[2].legend()
    axs[2].set_title('Infections at Multiple Sites')


    # Save output
    plt.savefig(file_name + '_anatomical_site.png')
    plt.close()


    # PLOTS OF PREVALENCE


    # Overall prevalence
    fig, axs = plt.subplots(3, 1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Simulated Prevalence Compared to Prevalence in STRIVE - Parameter Set ' + str(sim))


    # Overall
    axs[0].plot(tt, prvt[:,0])
    axs[0].plot(tt, len(tt) * [9.5], color = 'black', linestyle = '--')
    axs[0].set_title('Overall Prevalence')


    # Males
    axs[1].plot(tt, prvt[:,1])
    axs[1].plot(tt, len(tt) * [10.4], color = 'black', linestyle = '--')
    axs[1].set_title('Prevalence Amongst Males')


    # Females
    axs[2].plot(tt, prvt[:,2])
    axs[2].plot(tt, len(tt) * [8.9], color = 'black', linestyle = '--')
    axs[2].set_title('Prevalence Amongst Females')


    # Save output
    plt.savefig(file_name + '_prevalence_aggregate.png')
    plt.close()


    # Prevalence by age group and gender
    fig, axs = plt.subplots(4, 2)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Simulated Prevalence Compared to Prevalence in STRIVE - Parameter Set ' + str(sim))


    # Males
    axs[0, 0].set_title('Males')
    axs[0, 0].plot(tt, prvt[:,3])
    axs[0, 0].plot(tt, len(tt) * [21.6], color = 'black', linestyle = '--')
    axs[1, 0].plot(tt, prvt[:,4])
    axs[1, 0].plot(tt, len(tt) * [17.4], color = 'black', linestyle = '--')
    axs[2, 0].plot(tt, prvt[:,5])
    axs[2, 0].plot(tt, len(tt) * [11.6], color = 'black', linestyle = '--')
    axs[3, 0].plot(tt, prvt[:,6])
    axs[3, 0].plot(tt, len(tt) * [8.1], color = 'black', linestyle = '--')


    # Females
    axs[0, 1].set_title('Females')
    axs[0, 1].plot(tt, prvt[:,8])
    axs[0, 1].plot(tt, len(tt) * [20.1], color = 'black', linestyle = '--')
    axs[1, 1].plot(tt, prvt[:,9])
    axs[1, 1].plot(tt, len(tt) * [15.4], color = 'black', linestyle = '--')
    axs[2, 1].plot(tt, prvt[:,10])
    axs[2, 1].plot(tt, len(tt) * [7.3], color = 'black', linestyle = '--')
    axs[3, 1].plot(tt, prvt[:,11])
    axs[3, 1].plot(tt, len(tt) * [7], color = 'black', linestyle = '--')


    # Axis labels
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')


    # Save output
    plt.savefig(file_name + '_prevalence_age_group.png')
    plt.close()


#%% FUN check_for_output_data()
#
#
# Simply looks to see if the last output file is there for the specified
# scenario and parameter set combination.
#
#
def check_for_output_data(scenario, sim):
    return os.path.exists('simulations/calibration/scenario_' + str(scenario) + '/simulation_' + str(sim) + '/timestep' + str(param['partner_burn_in'][0] + param['simulation_length'][0] - 1) + '.ftr')


#%% FUN check_completion()
def check_completion():
    for scenario in [1, 2, 3]:
        done = 0
        for sim in range(0, 1000):
            if check_for_output_data(scenario, sim):
                done = done + 1
        print('Scenario ' + str(scenario) + ': ' + str(100*done/1000) + '% complete\n' )
