# PGS-calculations

## Description
Scripts that are used to calculate polygenic scores

## Files

`sumChromosomeScores.sh`:
  Shell script for merging PGSs over chromomsomes


`sumChromosomeScores.R`: 
  Corresponding R script that performs the actual merging of PGSs over chromosomes
  
  
`calculatePolygenicScores.sh`:
  Overarching script that submits `calculatePolygenicScoresForChromosome.sh` for every GWAS, chromosome
  
  
`calculatePolygenicScoresForChromosome.sh`
  Script that calculates PGSs for a particular GWAS and chromosome
