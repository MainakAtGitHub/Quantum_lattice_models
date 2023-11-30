#! /bin/bash
#SBATCH -p housewives
#SBATCH -n 1
#SBATCH --mem-per-cpu=2800
  parameter1=input_template.txt
  parameter2=
  parameter3=
subdir=out
  run_BdG_impurity_v3.sh /home/software/matlabR2012a-64/ $parameter1 $parameter2 $parameter3 > ./${subdir}/$parameter1$parameter2$parameter3.out 2>&1
