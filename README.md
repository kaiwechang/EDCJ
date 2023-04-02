## DCJ_Scaffolder

A DCJ-based genome scaffolder with maximum matching/intermediate/exemplar models.

### Requirements
1. `Gurobi`, the ILP Optimizer
2. `fmt` library, the modern cpp format library

### Usage
```bash
	>>> ./DCJ_Scaffolder -r <ref genome> -t <tar genome> [optional options]

	required options:

	 -r <ref genome> : reference genome path (".all" files from Sibelia)

	 -t <tar genome> : target genome path (".all" files from Sibelia)

	optional options:

	 -o <output dir> : output directory (default: output)

	 -m / -i / -e    : maximum matching / intermediate / exemplar model 
					   (default: exemplar, prioity: e > i > m)
	 -x              : extended version of speedup 3 (default: off)

	 -g <gap>        : gap tolerance for Gurobi (default: 0.0001)

	 -p <# of procs> : number of processors (default: all processors)

	 -l <time limit> : time limit for Gurobi (default: 1800)

	 -h              : usage instructions
```
