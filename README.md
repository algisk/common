# Perform multiple mutations using FoldX

## Description
This R script is used to replace specified multiple amino acids to 
a selected amino acid. It is useful for exploratory analysis and 
selecting potential proteins to be used in future experiments.

Initial purpose of this script was to select a target protein to be used
as a template for a Poly-Phe protein in 2016 Vilnius iGEM project.

---

## Dependencies
To use this R script, you need:
+ R ("data.table" and "foreach" packages)
+ FoldX
+ rotabase.txt
+ target .pdb file

All of these need to be in the same folder

---

## Usage
Simply type `Rscript --vanilla FoldX_template.R target.pdb number_of_runs omits.txt`. The number of runs is an integer and `omits.txt` is an optional
file containing a list of positions to ignore.
