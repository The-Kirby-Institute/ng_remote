# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:30:19 2022

@author: nicol
"""


# Import the module for running calibration
import numpy as np
import src.calibration.setup as cal


# Run the function which does on calibration run
random_no = np.random.choice(range(0, 10000))
cal.run_one_simulation(scenario = 1, 
                       parameter_no = random_no, 
                       run_mode = 'serial')
