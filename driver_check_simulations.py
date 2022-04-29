#  -*- coding: utf-8 -*-
"""
Created on Mon May 17 09:33:12 2021

@author: nicol
"""


# Libraries
from multiprocessing import Pool
import src.calibration.output as out

out.check_completion()


#%% RUN


# Set the list of simulations to consider
sim = list(range(0, 10000))


# Define function for handling the parallel pool
def pool_handler_prevalence():
    p = Pool(30)
    p.map(out.extract_prevalence, sim)


def pool_handler_checks():
    p = Pool(30)
    p.map(out.run_standard_checks, sim)


if __name__ == '__main__':
    pool_handler_prevalence()
    pool_handler_checks()







