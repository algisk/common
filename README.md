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
+ [FoldX](http://foldxsuite.crg.eu/)
+ rotabase.txt (comes with FoldX)
+ target.pdb file

All of these need to be in the same folder

---

## Usage
Simply type `Rscript --vanilla FoldX_template.R target.pdb number_of_runs omits.txt`. The number of runs is an integer and `omits.txt` is an optional
file containing a list of positions to ignore.

---

## Output
The script creates a folder called `out_[PDB_name]` which includes:
+ A list of performed mutations with free energy change (cumulative).
+ A plot showing energy change after each mutation compared to wild-type.
+ Another plot showing local minimums.
+ A folder containing FASTA files for each of these minimums.

---

## Notes
Some pdb files simply do not work because of major flaws in FoldX itself.
Currently, only the A chain of target protein will be mutated.
The target protein sequence must start at position 1.