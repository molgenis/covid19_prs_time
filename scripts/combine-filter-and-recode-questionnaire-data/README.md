# Fit the baseline models and perform the sensitivity analysis method

The scripts located in this directory were used to combine the questionnaire data of the LifeLines corona research project and filter and recode the questionnaire data

# Prerequisites

This program is developed in Python 3.7

The program requires the following packages to be installed:

* pandas ([v0.25.1](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))
* numpy ([v1.17.2](https://github.com/numpy/numpy/releases); [BSD License](https://www.numpy.org/license.html))* numpy ([v1.17.2](https://github.com/numpy/numpy/releases); [BSD License](https://www.numpy.org/license.html))
* statsmodels ([v0.12.2](https://www.statsmodels.org/stable/index.html); [BSD 3-Clause License](https://github.com/statsmodels/statsmodels/blob/main/LICENSE.txt))
* matplotlib ([v3.1.0](https://github.com/matplotlib/matplotlib/releases); [PSF License](https://matplotlib.org/3.1.0/users/license.html))
* seaborn ([v0.9.0](https://github.com/mwaskom/seaborn); [BSD 3-Clause License](https://github.com/mwaskom/seaborn/blob/master/LICENSE))


# Scripts

The scripts use commandline arguments for the input files. Description of all the necessary input files are available in the first lines of the main method of the scripts.

### Combine questionnaire data
This script combine the questionnaire data
```
python3 combine_questionairs_data.py
```

### filter and recode the questionnaire data
Script to filter and questionnaire data and recode the answers options to fit the models
```
python3 questionnaire_filter_participants_and_recode.py
```

### Combine PGS data
Combine the polygenic risk score files
```
python3 combine_pgs_scores.py
```