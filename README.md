# covid19_prs_time (later we can put the title here)

## Quality control

script: 

## PGS calculation

Polygenic scores (PGSs) are calculated using the scripts in `./scripts/pgs-calculations`.
The scripts are described [here](scripts/pgs-calculations).

## Relation PGSs and traits

### Longitudinal models

We selected ## questions for longitudinal models. The data was converted to long format and other preprocessing needed for the mixed models.
After fitting the models the plots are made using this script: ####

### Validation models

## Correlation with nation-wide statistics

A number of publicly available data sources reflecting the development
of the pandemic were selected for correlation with subjective quality of life.
After preprocessing of the PGSs and traits, `./scripts/nationWideCorrelations.R` 
is used for creating figures and correlations.

## Correlating outcome variables and PGSs in baseline sample set
We correlated 5 outcome variables that were assumed to tag the same genetics * time effect.
Additionally, we correlated polygenic scores for all PGSs in the baseline sample set.
After preprocessing of the PGSs and traits, `./scripts/correlationsForBaselineSamples.R`
is used for creating figures and calculating these correlations.

## Comparison of included samples and invited samples for polygenic scores

A comparison was done betwen samples that were invited in the studies and samples
that were included for the polygenic scores.
After preprocessing of the PGSs and traits, `./scripts/pgsParticipationComparison.R`
is used for creating figures and performing significance tests.

## Visualizing distributions of questions

After preprocesssing of the PGSs and traits, `./scripts/plotOutcomeDistributions.R`
can be used to plot the distributions of each of the outcome variables for which
we attempted to fit a time interaction.
