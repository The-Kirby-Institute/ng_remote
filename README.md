# Neisseria gonorrhoeae Transmission Model

This repository contains an individual-based mathematical model of *Neisseria gonorrhoeae* transmission among Aboriginal and Torres Strait Islander peoples living in remote areas of Australia. The model is used to assess the impact of vaccination in controlling *N. gonorrhoeae* prevalence among this population.

## Key Features

### Programming Languages
- **Python** (currently version 3.11.5): Main model implementation.
- **R**: Model fitting.

### Main Model Components
- **N. gonorrhoeae infection dynamics**: `ng.py`
- **Vaccination dynamics**: `vaccination.py`
- **Demographics**: `generate_population.py` and `population_dynamics.py`
- **Sexual partnership dynamics**: `partners.py`
  - Sexual behaviour data is based on the [GOANNA Study](https://youngdeadlyfree.org.au/wp-content/uploads/2021/01/GoannaSurvey2-FINAL.pdf).

## Usage Instructions

### Model Fitting
Use the R script `driver_calibration.R` to fit the model to the age and sex-specific prevalence reported in the [STRIVE Study](http://doi.org/10.1136/sextrans-2014-051617).

### Running Vaccination Scenarios
1. Execute `driver_vaccine_scenario.py` to run vaccination scenarios.
2. Define vaccination scenarios in `vax_example.py`.

### Data
The model parameterisation data is located in the `data` folder.

### Contributors to this Modelling Project:
The code for this project was originally developed by **Nicolas Rebuli** (School of Population Health/Kirby Institute, UNSW Sydney) and later taken over and modified by **Thilini Padeniya** (Institute for Biomedicine and Glycomics, Griffith University).

### Other Contributors:
- **Ben Hui** (Kirby Institute, UNSW Sydney)
- **David Regan** (Kirby Institute, UNSW Sydney)
- **James Wood** (School of Population Health, UNSW Sydney)
- **Kate Seib** (Institute for Biomedicine and Glycomics, Griffith University)



