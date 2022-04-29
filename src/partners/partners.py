# -*- coding: utf-8 -*-
"""
         Collection of functions for the remote communities model

INDEX
    find_partner: decides who a given person will partner with
    choose_relationship: decides if a given relationship will be long or short


"""

#%% SETUP Load Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Parse general simulation parameters
sim_parameters = pd.read_csv('data/param.csv')


#%% FUN setup_data()
#
#
# Setup partnership data
#
#
def setup_data(pop_parameters, run_mode = 'serial'):
    # Probability matrix of age biasing of last sexual encounter
    # Data from Table 5-3 of the GOANNA study - characteristics of last sexual encounter
    # Read as: from (row) to (column)
    # State space: {16-19, 20-24, 25-29, >29}
    # Note that the <16 demographics are excluded as they aren't in the model
    # The rows only correspond to the 16-19, 20-24, 25-29 age groups
    bias_age = pd.read_csv('data/age_partnership_distribution.csv').to_numpy()


    # Probability distribution of taking a high-risk partner by risk-group
    # Probabilties are made up
    # Columns: 0 (low-risk), 1(high-risk)
    # Row: the probability of taking a high-risk partner
    p_risky = pd.read_csv('data/probability_high_risk_partner.csv').to_numpy()


    # Probability distribution of cheating by risk group
    # Probabilities are made up
    # Columns: 0 (low-risk), 1 (high risk)
    # Row: the probability of cheating on a long-term partner
    p_cheat = pd.read_csv('data/probability_cheat.csv').to_numpy()


    # Probability distribution of nature of sexual partnerships
    # Data from Table 5-3 as above
    # Read as: age group (row) and relationship (column)
    # Relationship state space: {0 (long term), 1 (short term)}
    # Updated for new GOANNA survey (Table 5.5)
    bias_relationship = pd.read_csv('data/probability_relationship.csv').to_numpy()


    # The probability of a long-term relationship by relationship risk-group
    # The relationship risk-group defined as follows
    #    0 = 0 x high-risk | 2 x low-risk
    #    1 = 1 x high-risk | 1 x low-risk
    #    2 = 2 x high-risk | 0 x low-risk
    # Column: relationship risk-group
    aversion = pd.read_csv('data/scaling_long_term_by_risk_group.csv').to_numpy()


    # Probability of forming a new partnership
    # Adjusted to agree with the GOANNA survey data
    # Row: risk level, column: age-group
    
    # Reads in the default set which used to be the standard
    (print('Partner parameters: calibrated') if run_mode == 'serial' else [])
    p_new_partner = pd.read_csv('data/partnership_rates_scaling.csv')
    p_new_partner = (1/365) * p_new_partner
    # p_new_partner = (1/365) * pop_parameters['n'] * p_new_partner
    # p_new_partner.loc[0, :] = partner_rates.low[0] * p_new_partner.loc[0, :]
    # p_new_partner.loc[1, :] = partner_rates.high[0] * p_new_partner.loc[1, :]
    p_new_partner = p_new_partner.to_numpy()
    
        
    # Partnership duration parameters
    #
    #
    # Parameters of sample distribution for long-term relationships
    # Note that this is a Gamma distribution with seperate parameters for
    # each relationship risk-group
    # Row: parameter, column: relationship risk-group (see above)
    # duration_params = {"long": np.array([[365/100, 2*30/10, 2*30/10, 2*30/10],
    #                                      [100, 10, 10, 10]]),
    #                    "short": 14}
    #
    #
    partner_durations = pd.read_csv('data/partnership_durations.csv')
    duration_params = {'long': np.array([4 * [partner_durations.long_mean[0]/partner_durations.long_var[0]],
                                         4 * [partner_durations.long_var[0]]]),
                       'short': partner_durations.short[0]}


    # Combine parameters
    prt_parameters = {'bias_age': bias_age,
                      'p_risky': p_risky,
                      'p_cheat': p_cheat,
                      'bias_relationship': bias_relationship,
                      'aversion': aversion,
                      'p_new_partner': p_new_partner,
                      'duration_params': duration_params}


    return prt_parameters


#%% FUN update_partnerships()
#
#
# Top-level function
# This drives all the partnership updating
# It's basic function is to just call functions to make new partnerships
# and then to call funcitons to delete old partnerships.
#
#
def update_partnerships(t, prt_parameters, meta, partner_matrix, partner_expire):


    # Cleanup potential changes in the dtypes
    meta = meta.astype({'gender': 'int64',
                        'age_group': 'int64',
                        'orientation': 'int64',
                        'risk': 'int'})


    # Update partnership network
    meta, partner_matrix, partner_expire = new_partnership(prt_parameters, meta, partner_matrix, partner_expire, t)


    # Remove expired partnerships
    meta, partner_matrix, partner_expire = old_partnerships(meta, partner_matrix, partner_expire, t)


    return meta, partner_matrix, partner_expire


#%% FUN new_partnership()
#
# FUNCTION FOR MAKING NEW PARTNERSHIPS
#
# CREATE A NEW RELATIONSHIP FOR A GIVEN PERSON
#
#
# 1. Sample prob_partnership() to decide whether or not to make a new relationship
#
# 2. Run find_partner() to decide who the partner will be
#
# 3. Run choose_relationship() to decide if it will be a long or short-term relationship
#
# 4. Sample a duration from relationship_duration()
#
# 5. Update meta, partner_matrix, partner_expire accordingly.
#      Note that entering a long-term relationship causes an end to all
#      current short-term relationships
#
#
# INPUT
#   meta, partner_matrix, partner_expire, i, t
#
# OUTPUT
#   meta, partner_matrix, partner_expire
def new_partnership(prt_parameters, meta, partner_matrix, partner_expire, t):


    # Sample indicies of all who look for a new partner this iteration
    seek_partner = prob_partnership(prt_parameters, meta)
    seek_partner = np.random.random(len(meta)) < seek_partner
    seek_partner = np.where(seek_partner)[0]


    # Iterate over all people
    for i in range(0, len(seek_partner)):
        ii = seek_partner[i]


        # Find a new partner
        j = find_partner(prt_parameters, meta, partner_matrix, ii)


        # Check that a partner was indeed found
        if j != -1:


            # Decide on their relationship type
            is_short = choose_relationship(prt_parameters, meta, ii, j)


            # Sample a duration
            duration = relationship_duration(prt_parameters, meta, ii, j, is_short)


            # Updates for long-term relationships
            if is_short == 0:


                # Update partnership status
                meta.at[ii, "partner"] = j
                meta.at[j, "partner"] = ii


                # End all other relationships
                partner_matrix[ii, :] = 0
                partner_matrix[j, :] = 0
                partner_expire[ii, :] = float("inf")
                partner_expire[j, :] = float("inf")


            # Put partnership into matrix
            partner_matrix[ii, j] = 1
            partner_matrix[j, ii] = 1


            # Update partner counter
            meta.at[i, "counter"] = meta.at[i, "counter"] + 1
            meta.at[j, "counter"] = meta.at[j, "counter"] + 1


            # Update partnership duration matrix
            # print(i, j, is_short, duration)
            partner_expire[ii, j] = t + duration
            partner_expire[j, ii] = t + duration


    # Results
    return meta, partner_matrix, partner_expire


#%% FUN prob_partnership()
#
# FUNCTION FOR DECIDING THE PROBABILITY OF FORMING A RELATIONSHIP
#
# THE PER-DAY PROBABILITY OF AN INDIVIDUAL FORMING ANY KIND OF RELATIONSHIP
#
#
# These probabilities are based on a comparison between the resulting partnership
# aquisition rates and those sampled in the GOANNA survey.
#
# Some notes:
#
#   1. The probability of a partnership is given by
#
#          (expected number of partnerships per year) x scaling(age, risk) x cheat(risk)
#
#       In this equation, the expected number of partners is different for
#       the high-risk and the low-risk groups. The the scaling for the low-risk
#       group has been based on the GOANNA numbers and the number for the
#       high-risk group has been tuned.
#
#       The scaling function adjusts the partnership aquisition rates for
#       discrepencies between age groups. For example, it was determined
#       that for the middle age group low-risk partnerships are relatively
#       infrequent while high-risk partnerships are more frequent.
#
#       The cheating function adjusts the partnership aquisition rate to
#       account for cheating behaviour. In particular, how likely is each
#       risk group to cheat on a long-term partner. Note that concurrent
#       short-term partnerships is not considered cheating.
#
#
# INPUT
#   meta = the population array
#   i = the index of the person under consideration
#
# OUTPUT
#   the probability of this individual entering a new partnership
def prob_partnership(prt_parameters, meta):


    # Vectorised version
    p_partner_it = ( 1.0 * (meta.partner == -1) + \
                   ( 1.0 * (meta.partner != -1) * prt_parameters['p_cheat'][0, meta.risk] ) ) * \
                   prt_parameters['p_new_partner'][meta.risk, meta.age_group]


    # Return the partnership formation probability
    return p_partner_it


#%% FUN find_partner()
#
# FUNCTION FOR MAKING PARTNERSHIPS
#
# FIND A SEXUAL PARTNER FOR PERSON i FROM THE POOL OF ELIGABLE SINGLES IN meta
#
#
# Decision tree is as follows:
#
#   1. Use the sexual orientation of the bachelor to narrow the population
#       down into the individuals who would they would be interested in
#       partnering with.
#
#   2. The population is put into three groups based on age:
#       16-19, 20-24, 25-29 and > 29.
#
#   3. Using data from the GOANNA study, the age-preferences of the bachelor
#       are identified based on their age group. This distribution is then
#       multiplied by the age group preferences of their gender.
#       ie.
#       P(partner age|bachelor age, bachelor gender) =
#           P(partner age|bachelor age) X P(partner age|bachelor gender)
#
#   4. A selection of partners is extracted from the pool of possible
#       partners based on the age group distribution above.
#
#   5. A partner is then selected at random from that pool.
#
#
# INPUT
#   meta = the population array
#   i = the index of the bachelor
#
# OUTPUT
#   j = the index of their new partner.
#
#
def find_partner(prt_parameters, meta, partner_matrix, bachelor_index):


    ######################################################
    ##  NARROW DOWN POPULATION TO AVAILABLE CANDIDATES  ##
    ######################################################
    # Pick out the bachelor
    bachelor = meta.loc[bachelor_index, :]


    # Pick out available partners
    if bachelor.orientation == 0:
        # Pick out partners for heterosexuals
        # Could be a heterosexual or a bisexual but must be of the opposite gender
        partners = meta[(meta.orientation != 1) &\
                        (meta.gender != bachelor.gender) &\
                        (partner_matrix[bachelor_index, 0:len(meta)]==0) &\
                        (meta.index != bachelor_index)]


    elif bachelor.orientation == 1:
        # Pick out partners for homosexuals
        # Could be a homosexual or a bisexual of the same gender
        partners = meta[(meta.orientation != 0) &\
                        (meta.gender == bachelor.gender) &\
                        (partner_matrix[bachelor_index,0:len(meta)]==0) &\
                        (meta.index != bachelor_index)]


    else:
        # Pick out partners for bisexuals
        # Could be anyone but a hetero of the same sex or a homo of the other sex
        partners = meta[~((meta.orientation == 0) & (meta.gender == bachelor.gender)) &\
                        ~((meta.orientation == 1) & (meta.gender != bachelor.gender)) &\
                        (partner_matrix[bachelor_index,0:len(meta)]==0) &\
                        (meta.index != bachelor_index)]


    #############################################
    ##  DECIDE WHICH AGE GROUP TO PARTER WITH  ##
    #############################################
    if len(partners) > 0:


        # Pull out age biasing distribution
        partner_dist = prt_parameters['bias_age'][bachelor.age_group, :] #* bias_sex[int(bachelor["gender"])]


        # Check that there's a parter from each age group available
        to_ignore = partners.groupby("age_group")["age_group"].count().reset_index(name = "count")
        to_ignore = to_ignore[to_ignore["count"] == 0]


        # Make sure there are still some partners available
        if len(to_ignore) < 5:


            # Some age groups may still be empty - set frequency of such
            # groups to zero so it's impossible to end up choosing them
            partner_dist[to_ignore.age_group.values] = 0


            # Calculate CDF of gender distribution
            partner_dist = np.cumsum(partner_dist)
            partner_dist = partner_dist/max(partner_dist)


            # Decide which age-group to partner with
            partner_age_group = np.searchsorted(partner_dist, np.random.random(1))
            partners = partners[partners.age_group == partner_age_group[0]]


            ###############################################
            ##  DECIDE WHICH RISK GROUP TO PARTNER WITH  ##
            ###############################################


            # Test that there are some
            if len(np.unique(partners.risk)) > 1:


                # Decide which risk-group to partner with
                if np.random.random(1) < prt_parameters['p_risky'][0, bachelor.risk]:

                    # High-risk
                    partners = partners[partners["risk"] == 1]

                else:

                    # Low-risk
                    partners = partners[partners["risk"] == 0]


            ###################################################
            ##  DECIDE WHICH PARTNER STATUS TO PARTNER WITH  ##
            ###################################################


            # Test that there are some of each group
            if len(np.unique(partners.partner)) > 1:


                # Decide whether or not to cheat
                if np.random.random(1) < prt_parameters['p_cheat'][0, bachelor.risk]:

                    # Cheating
                    partners = partners[partners["partner"] != -1]

                else:
                    # Not cheating
                    partners = partners[partners["partner"] == -1]


            ######################
            ##  FINAL DECISION  ##
            ######################


            # Now just choose one at random
            if len(partners) > 0:
                partner = np.random.choice(partners.index, 1)
                partner = int(partner)
            else:
                partner = -1
        else:
            partner = -1
    else:
        partner = -1


    # Return the meta array but with the updated partner status
    return partner


#%% FUN choose_relationship()
#
# FUNCTION FOR DECIDING ON THE TYPE OF RELATIONSHIP
#
# DECIDE ON THE RELATIONSHIP (CASUAL vs LONG TERM) BETWEEN TWO INDIVIDUALS
#
#
# Decision tree is relatively simple:
#
#   1. Look to see if either i or j are in a long term relationship.
#
#   2.a If at least one of i and j are in a long term relationship, then
#         this is a short term relationship.
#
#   2.b Otherwise make a decision using data from the GOANNA study.
#
#
# INPUT
#   meta = the population array
#   i = the index of the person finding a partner
#   j = the index of their partner
#
# OUTPUT
#   relationship = {0 (long term), 1 (short term)}
#
def choose_relationship(prt_parameters, meta, i, j):


    ##################################################
    ##  DECIDE ON THE NATURE OF THEIR RELATIONSHIP  ##
    ##################################################


    # Check that neither of the individuals are in a long term relationship
    if meta.at[i, "partner"] + meta.at[j, 'partner'] > -2:


        # If already in a long term relationship - set this to be a short term
        relationship = 1


    # Otherwise make a decision at random based on the age group of i
    # and their combined risk grouping
    else:


        # Decide if relationship is long term or short term
        risk_group = meta.at[i, "risk"] + meta.at[j, "risk"]
        p_relationship = prt_parameters['bias_relationship'][meta.age_group[i]]
        p_relationship = np.cumsum(p_relationship)
        p_relationship[0] = prt_parameters['aversion'][0, int(risk_group)] * p_relationship[0]


        # Now decide on a relationship type at random
        relationship = int(np.searchsorted(p_relationship, np.random.random(1)))


    # print(i, j, relationship)
    # Return which type of relationship has been chosen
    return relationship


#%% FUN relationship_duration()
#
# FUNCTION FOR SAMPLING RELATIONSHIP DURATION
#
# SAMPLE THE DURATION OF A PARTNERSHIP
#
#
# Short term relationships are sampled from an exponential distribution with
# a mean of 5 days.
#
# Long term relationships are sampled from a Gamma distribution with a seperate
# set of parameters for each relationhip risk-group.
#
#
# INPUT
#   meta, i, j
#   is_short = {0=long term relationship, 1=short term relationship}
#
# OUTPUT
#   the duration of the relationship
def relationship_duration(prt_parameters, meta, i, j, is_short):


    # Sample a duration
    if is_short == 0:
        risk_group = meta.at[i, "risk"] + meta.at[j, "risk"]
        duration = np.random.gamma(prt_parameters['duration_params']['long'][0, int(risk_group)],
                                   prt_parameters['duration_params']['long'][1, int(risk_group)])
        #duration = 1000
    else:
        duration = np.random.exponential(prt_parameters['duration_params']['short'])
        #duration = 1


    # Return duration
    return duration


#%% FUN old_partnerships()
#
# FUNCTION FOR REMOVING OLD PARTNERSHIPS
#
# REMOVE EXPIRED RELATIONSHIPS FROM THE POPULATION
#
#
# Check the partner_expire array to see if any relationships have expired.
# Remove all expired relationships from the array and the partner_matrix array.
# Check meta to see if any of the expired relationships are long-term
# and update meta accordingly.
#
# INPUT
#   meta, partner_matrix, partner_expire, i
#
# OUTPUT
#   meta, partner_matrix, partner_expire
def old_partnerships(meta, partner_matrix, partner_expire, t):


    # Identify partnerships to terminate
    [ii, jj] = np.where(partner_expire < t)


    # Update meta to remove long term partnerships
    ll = meta.partner[ii].values
    ll = ll[ll == jj]
    meta.loc[ll, 'partner'] = -1


    # Update matrices to remove short term partnerships
    partner_matrix[ii, jj] = 0
    partner_expire[ii, jj] = float("inf")


    # Return output
    return meta, partner_matrix, partner_expire


#%% GRAPH initilise_partnership_trackers()
#
#
#
#
#
def initilise_partnership_trackers(n_steps):
    partners_cohort = np.zeros((n_steps, 4*10))
    partners_risk = np.zeros((n_steps, 4*10))
    partners_cum_risk = np.zeros((n_steps, 4*10))
    partners_cum_age = np.zeros((n_steps, 4*10))
    partners_cum_tot = np.zeros((n_steps, 4*5))
    return partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot


#%% GRAPH update_partnership_trackers()
#
#
#
#
#
def update_partnership_trackers(t, pop_parameters, meta, partner_matrix, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot):

    
    # Identifier for people who have been in the population for over 365 days
    age_okay = meta.import_time < sim_parameters['partner_burn_in'][0] + t - 365
    

    # Run partnership counters by cohorts of age/risk and gender
    ii = 0
    kk = 0
    ll = 0
    for gORr in [0, 1]:
        for a in [0, 1, 2, 3, 4]:


            # Cumulative number of partners by age group
            if gORr == 0:
                who = pop_parameters['lookup']['A' + str(a)] & age_okay
                n_cohort = sum(who)
                if n_cohort > 0:
                    partners_cum_tot[t, ll] = sum(who & (meta.counter==0))/n_cohort
                    partners_cum_tot[t, ll + 1] = sum(who & (meta.counter==1))/n_cohort
                    partners_cum_tot[t, ll + 2] = sum(who & (meta.counter.isin([2, 3, 4])))/n_cohort
                    partners_cum_tot[t, ll + 3] = sum(who & (meta.counter>=5))/n_cohort
                ll = ll + 4


            # Cumulative number of partners by risk group
            who = pop_parameters['lookup']['RA' + str(gORr) + str(a)] & age_okay
            n_cohort = sum(who)
            if n_cohort > 0:
                partners_cum_risk[t, kk] = sum(who & (meta.counter == 0))/n_cohort
                partners_cum_risk[t, kk + 1] = sum(who & (meta.counter == 1))/n_cohort
                partners_cum_risk[t, kk + 2] = sum(who & (meta.counter.isin([2, 3, 4])))/n_cohort
                partners_cum_risk[t, kk + 3] = sum(who & (meta.counter >= 5))/n_cohort


            # Cumulative number of partners by age group
            who = pop_parameters['lookup']['GA' + str(gORr) + str(a)] & age_okay
            n_cohort = sum(who)
            if n_cohort > 0:
                partners_cum_age[t, kk] = sum(who & (meta.counter == 0))/n_cohort
                partners_cum_age[t, kk + 1] = sum(who & (meta.counter == 1))/n_cohort
                partners_cum_age[t, kk + 2] = sum(who & (meta.counter.isin([2, 3, 4])))/n_cohort
                partners_cum_age[t, kk + 3] = sum(who & (meta.counter >= 5))/n_cohort
            kk = kk + 4


            # Partners by relationship type
            for p in ['P', 'PL', 'PS', 'PC']:
                jj = 0 if p == 'P' else 1
                denom = sum(pop_parameters['lookup']['GA' + str(gORr) + str(a)])
                if denom > 0:
                    partners_cohort[t, ii] = sum(pop_parameters['lookup'][p + 'GA' + str(jj) + str(gORr) + str(a)]) / denom
                denom = sum(pop_parameters['lookup']['RA' + str(gORr) + str(a)])
                if denom > 0:
                    partners_risk[t, ii] = sum(pop_parameters['lookup'][p + 'RA' + str(jj) + str(gORr) + str(a)]) / denom
                ii = ii + 1


    # Reset the partner counter after 365 iterations
    if np.mod(t, 365) == 0:
        meta.loc[:, 'counter'] = np.sum(partner_matrix, axis = 0)[0:len(meta)]


    return partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot


#%% GRAPH make_partnership_graphs()
#
#
#
#
#
def make_partnership_graphs(tt, pop_parameters, partners_cohort, partners_risk, partners_cum_risk, partners_cum_age, partners_cum_tot, save_loc = 'graphs/output/', show_graphs = False):


    ## GRAPH SHOWING THE PROPORTION OF PEOPLE IN EACH PARTNERSHIP TYPE BY AGE AND GENDER COHORT


    # Initilise figure
    fig, axs = plt.subplots(5, 2)
    
    
    # Net change in population size
    ii = 0
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            axs[age, gender].plot(tt, partners_cohort[:, ii + 0], label = 'No Partner')
            axs[age, gender].plot(tt, partners_cohort[:, ii + 1], label = 'Long-Term')
            axs[age, gender].plot(tt, partners_cohort[:, ii + 2], label = 'Short-Term')
            axs[age, gender].plot(tt, partners_cohort[:, ii + 3], label = 'Cheating')
            ii = ii + 4
            
    
    # Axis labels
    plt.suptitle('Proportion of Popultion with in Each Partnership Configuration\nBy Age and Gender Cohort')
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
    if show_graphs:
        plt.show()
    else:
        fig.savefig(save_loc + 'partnerships_by_age_and_gender.png', dpi = 200)
        plt.close(fig)
    

    ## GRAPH SHOWING THE PROPORTION OF PEOPLE IN EACH PARTNERSHIP TYPE BY AGE AND RISK COHORT


    # Initilise figure
    fig, axs = plt.subplots(5, 2)
    
    
    # Net change in population size
    ii = 0
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            axs[age, gender].plot(tt, partners_risk[:, ii + 0], label = 'No Partner')
            axs[age, gender].plot(tt, partners_risk[:, ii + 1], label = 'Long-Term')
            axs[age, gender].plot(tt, partners_risk[:, ii + 2], label = 'Short-Term')
            axs[age, gender].plot(tt, partners_risk[:, ii + 3], label = 'Cheating')
            ii = ii + 4
            
    
    # Axis labels
    plt.suptitle('Proportion of Popultion with in Each Partnership Configuration\nBy Age and Risk Cohort')
    axs[0, 0].set_title('Low-Risk')
    axs[0, 1].set_title('High-Risk')
    axs[4, 0].legend(loc = 'upper center', ncol = 3, bbox_to_anchor = (1.1, -0.05), fancybox = True)
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')
    
    
    # Save the graph
    if show_graphs:
        plt.show()
    else:
        fig.savefig(save_loc + 'partnerships_by_age_and_risk.png', dpi = 200)
        plt.close(fig)


    ## GRAPH SHOWING THE NUMBER OF PEOPLE WITH 0, 1, 2-4 and >5 PARTNERS
    ## IN THE LAST 12 MONTHS BY AGE AND RISK COHORT


    # Initilise figure
    fig, axs = plt.subplots(5, 2)
    
    
    # Net change in population size
    ii = 0
    for risk in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            axs[age, risk].plot(tt, partners_cum_risk[:, ii + 0], label = 'No Partners')
            axs[age, risk].plot(tt, partners_cum_risk[:, ii + 1], label = '1 partner')
            axs[age, risk].plot(tt, partners_cum_risk[:, ii + 2], label = '2-4 partners')
            axs[age, risk].plot(tt, partners_cum_risk[:, ii + 3], label = '5 or more partners')
            ii = ii + 4
            
            
    # Add lines to show
    goanna = pd.read_csv('data/calibration_partnership_rates.csv')
    col = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for age in [0, 1, 2, 3, 4]:
        for risk in [0, 1]:
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[0, 0], len(tt)), color = col[0], linestyle = '--')
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[1, 0], len(tt)), color = col[1], linestyle = '--')
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[2, 0], len(tt)), color = col[2], linestyle = '--')
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[3, 0], len(tt)), color = col[3], linestyle = '--')
    
    
    # Axis labels
    plt.suptitle('Comparison of the Proportion of Popultion with in Each Partnership Configuration to GOANNA\nBy Age and Risk Cohort')
    axs[0, 0].set_title('Low-Risk')
    axs[0, 1].set_title('High-Risk')
    axs[4, 0].legend(loc = 'upper center', ncol = 3, bbox_to_anchor = (1.1, -0.05), fancybox = True)
    axs[0, 0].set_ylabel('16-19')
    axs[1, 0].set_ylabel('20-24')
    axs[2, 0].set_ylabel('25-29')
    axs[3, 0].set_ylabel('30-34')
    axs[4, 0].set_ylabel('35')
    axs[4, 0].set_xlabel('Day')
    axs[4, 1].set_xlabel('Day')
    
    
    # Save the graph
    if show_graphs:
        plt.show()
    else:
        fig.savefig(save_loc + 'partnerships_cumulative_vs_goanna_by_age_risk.png', dpi = 200)
        plt.close(fig)


    ## GRAPH SHOWING THE NUMBER OF PEOPLE WITH 0, 1, 2-4 and >5 PARTNERS
    ## IN THE LAST 12 MONTHS BY AGE AND GENDER COHORT


    # Initilise figure
    fig, axs = plt.subplots(5, 2)
    
    
    # Net change in population size
    ii = 0
    for gender in [0, 1]:
        for age in [0, 1, 2, 3, 4]:
            
            # Graphs for individual age and gender groups
            axs[age, gender].plot(tt, partners_cum_age[:, ii + 0], label = 'No Partners')
            axs[age, gender].plot(tt, partners_cum_age[:, ii + 1], label = '1 partner')
            axs[age, gender].plot(tt, partners_cum_age[:, ii + 2], label = '2-4 partners')
            axs[age, gender].plot(tt, partners_cum_age[:, ii + 3], label = '5 or more partners')
            
            # Graphs of
            
            ii = ii + 4
            
            
    # Add lines to show
    # goanna = pd.read_csv('data/calibration_partnership_rates.csv')
    col = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for age in [0, 1, 2, 3, 4]:
        for risk in [0, 1]:
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[0, 0], len(tt)), color = col[0], linestyle = '--')
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[1, 0], len(tt)), color = col[1], linestyle = '--')
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[2, 0], len(tt)), color = col[2], linestyle = '--')
            axs[age, risk].plot(tt, np.repeat(goanna.iloc[3, 0], len(tt)), color = col[3], linestyle = '--')
            
            
    # Axis labels
    plt.suptitle('Comparison of the Proportion of Popultion with in Each Partnership Configuration to GOANNA\nBy Age and Gender Cohort')
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
    if show_graphs:
        plt.show()
    else:
        fig.savefig(save_loc + 'partnerships_cumulative_vs_goanna_by_age_gender.png', dpi = 200)
        plt.close(fig)


    ## GRAPH SHOWING THE NUMBER OF PEOPLE WITH 0, 1, 2-4 and >5 PARTNERS
    ## IN THE LAST 12 MONTHS BY AGE COHORT


    # Initilise figure
    fig, axs = plt.subplots(5, 1)
    
    
    # Net change in population size
    labs = ['No partners', '1 partner', '2-4 partners', '5 or more partners']
    for age in [0, 1, 2, 3, 4]:
        for j in [0, 1, 2, 3]:
            
            col = next(axs[age]._get_lines.prop_cycler)['color']
            y = partners_cum_tot[:, int(4*age + j)]
            z = len(tt) * [pop_parameters['partners_dist'].iloc[int(j), age]]
            axs[age].plot(tt, y, label = labs[int(j)], color = col)
            axs[age].plot(tt, z, linestyle='--', color = col)
    
    
    # Axis labels
    plt.suptitle('Proportion of Popultion with in Each Partnership Configuration by Age Cohort')
    axs[4].legend(loc = 'upper center', ncol = 2, bbox_to_anchor = (0.5, -0.05), fancybox = True)
    axs[0].set_ylabel('16-19')
    axs[1].set_ylabel('20-24')
    axs[2].set_ylabel('25-29')
    axs[3].set_ylabel('30-34')
    axs[4].set_ylabel('35')
    axs[4].set_xlabel('Day')
    
    
    # Save the graph
    if show_graphs:
        plt.show()
    else:
        fig.savefig(save_loc + 'partnerships_cumulative_vs_goanna_by_age.png', dpi = 200)
        plt.close(fig)


    ## BAR GRAPH SHOWING THE NUMBER OF PEOPLE WITH 0, 1, 2-4 and >5 PARTNERS
    ## IN THE LAST 12 MONTHS BY AGE COHORT


    # Initilise figure
    fig, axs = plt.subplots(5, 1)


    # Net change in population size
    width = 0.35
    x = np.arange(4)
    n = len(partners_cum_tot)-1
    labs = ['No partners', '1 partner', '2-4 partners', '5 or more partners']
    for age in [0, 1, 2, 3, 4]:

        # Get data for graph
        y = [partners_cum_tot[n, 4*age],
             partners_cum_tot[n, 4*age + 1],
             partners_cum_tot[n, 4*age + 2],
             partners_cum_tot[n, 4*age + 3]]

        # Make graph
        axs[age].bar(x - width/2, y, width, label = 'Simulated')
        axs[age].bar(x + width/2, pop_parameters['partners_dist'].iloc[:, int(age)], width, label = 'Target')
        axs[age].set_xticks(x)
        axs[age].set_ylim(0, 0.6)
        if age < 4:
            axs[age].set_xticklabels(4*[''])
        else:
            axs[age].set_xticklabels(labs)


    # Axis labels
    plt.suptitle('Calibration of Partnership Rates')
    axs[0].set_title('The Proportion of the Population in Each Partnership Rate Category')
    axs[0].set_ylabel('16-19')
    axs[1].set_ylabel('20-24')
    axs[2].set_ylabel('25-29')
    axs[3].set_ylabel('30-34')
    axs[4].set_ylabel('35')
    axs[4].set_xlabel('Day')
    axs[4].legend()
    fig.tight_layout()
    
    
    # Save the graph
    if show_graphs:
        plt.show()
    else:
        fig.savefig(save_loc + 'partnerships_cumulative_vs_goanna_by_age_bar.png', dpi = 200)
        plt.close(fig)




#%% CALIBRATION generate_parameters()
#
#
# Makes a parameter set to test with calibration
#
#
def generate_calibration_parameters():
    
    
    # Grid for low-risk people
    prt_lower_low = 0.5
    prt_upper_low = 3
    n_low = (prt_upper_low - prt_lower_low)/(2/3)
    grid_low = np.linspace(prt_lower_low, prt_upper_low, int(n_low) + 1)
    
    
    # Grid for high-risk people
    prt_lower_high = 5
    prt_upper_high = 20
    n_high = (prt_upper_high - prt_lower_high)/5
    grid_high = np.linspace(prt_lower_high, prt_upper_high, int(n_high) + 1)
    
    
    # Compute the full mesh-grid
    # l0, l1, l2, l3, l4, h0, h1, h2, h3, h4 = \
    #     np.meshgrid(grid_low, grid_low, grid_low, grid_low, [0], grid_high, grid_high, grid_high, grid_high, [0])
    
    
    # Too slow. Use sampling instead.
    n_sam = sim_parameters['n_prt_param_sets'][0]
    l0 = prt_lower_low + (prt_upper_low - prt_lower_low) * np.random.random(n_sam)
    l1 = prt_lower_low + (prt_upper_low - prt_lower_low) * np.random.random(n_sam)
    l2 = prt_lower_low + (prt_upper_low - prt_lower_low) * np.random.random(n_sam)
    l3 = prt_lower_low + (prt_upper_low - prt_lower_low) * np.random.random(n_sam)
    h0 = prt_lower_high + (prt_upper_high - prt_lower_high) * np.random.random(n_sam)
    h1 = prt_lower_high + (prt_upper_high - prt_lower_high) * np.random.random(n_sam)
    h2 = prt_lower_high + (prt_upper_high - prt_lower_high) * np.random.random(n_sam)
    h3 = prt_lower_high + (prt_upper_high - prt_lower_high) * np.random.random(n_sam)
    z = [0]*n_sam
    
    # Reshape into a csv for the model input
    param = pd.DataFrame({'l0': l0, \
                          'l1': l1, \
                          'l2': l2, \
                          'l3': l3, \
                          'l4': z, \
                          'h0': h0, \
                          'h1': h1, \
                          'h2': h2, \
                          'h3': h3, \
                          'h4': z})
    
        
    # Store parameters
    param.to_csv('simulations/parameters_partnerships.csv', index = False)
    
    
#%% CALIBRATION parse_calibration_parameters()
#
#
# Unpacks calibration parameters
#
#
def parse_calibration_parameters(par_set):
    

    # Read in the specific parameter set
    p = pd.read_csv('simulations/parameters_partnerships.csv')
    p = p.iloc[par_set, :]
    
    
    # Convert to the correct format
    p_new_partner = np.array([[p.iloc[0], p.iloc[1], p.iloc[2], p.iloc[3], p.iloc[4]], \
                              [p.iloc[5], p.iloc[6], p.iloc[7], p.iloc[8], p.iloc[9]]])
    p_new_partner = (1/365) * p_new_partner
    
    
    # Done
    return p_new_partner















