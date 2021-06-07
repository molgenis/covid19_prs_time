# Fit the baseline models and perform the sensitivity analysis method

The scripts located in this this directory were used to calculate the baseline moddels between questionnaires and Polygenic risk scores.

The pipeline contains four steps
1) Fit the properate models between the questionnaire items and the PGS
2) Perform the meta analysis between the results of the two genetic chips
3) Identify the time relationship of the correlation coefficients over time to perform the sensitivity analysis of the longitudinal results
4) Plot the outcomes of the baseline associations

The first two steps were used for both identify the baseline associations and the sensitivity analysis. The thirth step were only used in the validation method. The fourth step were used to plot the baseline outcome associations.

# Prerequisites

This program is developed in Python 3.7

The program requires the following packages to be installed:

* pandas ([v0.25.1](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))
* numpy ([v1.17.2](https://github.com/numpy/numpy/releases); [BSD License](https://www.numpy.org/license.html))* numpy ([v1.17.2](https://github.com/numpy/numpy/releases); [BSD License](https://www.numpy.org/license.html))
* statsmodels ([v0.12.2](https://www.statsmodels.org/stable/index.html); [BSD 3-Clause License](https://github.com/statsmodels/statsmodels/blob/main/LICENSE.txt))
* matplotlib ([v3.1.0](https://github.com/matplotlib/matplotlib/releases); [PSF License](https://matplotlib.org/3.1.0/users/license.html))
* seaborn ([v0.9.0](https://github.com/mwaskom/seaborn); [BSD 3-Clause License](https://github.com/mwaskom/seaborn/blob/master/LICENSE))

# Scripts

Every step as described earlier contains it's own python script. 
The scripts use commandline arguments for the input files. Description of all the necessary input files are availbale in the first lines of the main method of the scripts.

### The first step
Fit the properate models between the questionnaire items and the PGS
```
python3 fiteBaselineModelsBetweenQuestionsAndPGS.py
```

### The second step
Perform the meta analysis between the results of the two genetic chips
```
python3 PerformMetaAnalysis.py
```

### The third step
Identify the time relationship of the correlation coefficients over time to validate the longitudinal results
```
python3 performValidationAnalysisOverBaselineModels.py
```
### The fourth step
plot the baseline outcomes between the questionnaire items and the PGS scores
```
python3 plots_pgs_baseline_outcome_items.py
```





