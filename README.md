## EDCJ-Scaffolder

A scaffolding tool based on exemplar double-cut-and-join distance

### Dependency

- required
	- [`Gurobi`](https://www.gurobi.com/downloads/), the ILP Optimizer
	- [`fmt`](https://github.com/fmtlib/fmt) library, the modern cpp format library
- optional
	- [`Sibelia`](https://github.com/bioinf/Sibelia), the synteny blocks (markers) finding tools
	- [`MUMmer`](https://github.com/mummer4/mummer), the genome alignment tool

### Usage

```bash
	$ ./EDCJ_Scaffolder -r <ref genome> -t <tar genome> [optional options]
```

option details:
```
	required options:
	
	 -r <ref genome> : reference genome path (".all" files from Sibelia)
	
	 -t <tar genome> : target genome path (".all" files from Sibelia)
	
	optional options:
	
	 -o <output dir> : output directory (default: output)
	
	 -x              : extended version of speedup 3 (default: off)
	
	 -g <gap>        : gap tolerance for Gurobi (default: 0.0001)
	
	 -p <# of procs> : number of processors (default: all processors)
	
	 -l <time limit> : time limit for Gurobi (default: 1800)
	
	 -h              : usage instructions
```

### Output Format

`scaffolds.txt` in the output directory
```
	> Scaffold_1
	<contig_name_a> <orientation_a>
	<contig_name_b> <orientation_b>
	<contig_name_c> <orientation_c>
	...
	
	> Scaffold_2
	<contig_name_d> <orientation_d>
	<contig_name_e> <orientation_e>
	<contig_name_f> <orientation_f>
	...
```

### Evaluation

compare scaffolds output from EDCJ_Scaffolder with ground truth answer
```bash
	$ tools/misJoin_eval <ground truth> <scaffolds from EDCJ_Scaffolder>
```

which outputs to `stdout` with the following format:
```
	Scaffolds:        %d
	Multi-contig:     %d
	MisJoins:         %d
	Joins:            %d
	Sensitivity:    %.3f
	Precision:      %.3f 
	Fscore:         %.3f 
```

### Demo with Example testcase

```bash
	$ make
	$ ./EDCJ_Scaffolder -r example/reference.all -t example/target.all -o output
	$ tools/misJoin_eval example/answerToAll output/scaffolds.txt
```

or just
```bash
	$ make run_example
```

then, we should see the following output
```bash
	Scaffolds:         11
	Multi-contig:       6
	MisJoins:           0
	Joins:             47
	Sensitivity:    0.825
	Precision:      1.000
	Fscore:         0.904
```

