# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:30:19 2022

@author: nicol
"""


# Import the module for running calibration
import src.calibration.setup as cal


# Run the function which does on calibration run
cal.run_one_simulation(scenario = 3, 
                       parameter_no = 0, 
                       run_mode = 'serial')
