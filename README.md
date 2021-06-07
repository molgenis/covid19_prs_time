# covid19_prs_time (later we can put the title here)

## Quality control

Sample quality control is done in: `./scripts/sampleQc.R`

## PGS calculation

Polygenic scores (PGSs) are calculated using the scripts in `./scripts/pgs-calculations`.
The scripts are described [here](scripts/pgs-calculations).

## Relation PGSs and traits

### Longitudinal models

The data was converted to long format and other preprocessing needed for the mixed models are done here: `/scripts/prepareForLongitudinalModels.R`

The longitudinal models are fitted by: `/scripts/longitudinalModels.R` 

The plots of the longitudinal models are made using: `/scripts/plotLongitudinalModels.R`

### Sensitivity analysis

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

## Manuscript figures

* Fig2: `./scripts/plotHeatmapBaselineModels.R`
* Fig3: `./scripts/qolPlot.R`
* Fig4: `./scripts/c19Plot.R`
