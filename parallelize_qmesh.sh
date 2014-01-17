#/bin/bash

# bash script to parallize (trivially) the k-mesh calculation of impurity_dos
# Usage parallize_qmesh.sh parameterfile number_of_tasks submit_script_template

replace_input() {
#/bin/bash
if [ $# -eq 3 ]
then
	# replace all lines containing $2 with a line that sets $2 = $3
	string='/'$2'/c\'$2' = '$3
else
	string=' / '$2'/c\  '$2'='$3
fi
#echo $string
sed -i "$string" $1
}


# check whether all inputs are set correctly and necessary files exist

if [ $# -lt 1 ]
	then
		echo "Too few arguments"
	echo "Usage parallelize_qmesh.sh parameterfile number_of_tasks submit_script_template"
	exit
fi
if [ $# -lt 2 ]
then
	echo "Setting number_of_tasks to 2"
	num_tasks=2
else
	if [ $2 -lt 2 ]
	then      
		echo "number_of_tasks has to be greater than one"
		echo "Using number_of_tasks=2 now."
		num_tasks=2
	else
		num_tasks=$2
	fi
fi

if [ ! -e "$1" ]
then
	echo "Error: parameterfile $1 does not exist"
	exit
fi
inputfile=$1

template=$3

batchcommand='sbatch'
holdflag='--hold'
subdir='out'
if [ ! -e "$3" ]
then
	echo "Error submit_script_template $3 does not exist, using standart ones:"
	while true;
       	do
		read -p "[I]TP (Slurm), [H]PC-script " ih
		case $ih in
			[iI]* ) 
				batchfile='script.sh'
				echo "#! /bin/bash" > ${batchfile}
				# check for available queue for current user
				numdfg=`sinfo | grep dfg | wc -l`
				if [ "$numdfg" -eq "0" ]
					echo "#SBATCH -p housewives" >> ${batchfile}
				else
					echo "#SBATCH -p dfg" >> ${batchfile}
				fi
				echo "#SBATCH -n 1" >> ${batchfile}
				echo "#SBATCH --mem-per-cpu=3800" >> ${batchfile}
				echo "  parameter1=standart_input_imp_dos_M_20.mat" >> ${batchfile}
				echo "  parameter2=" >> ${batchfile}
				echo "  parameter3=" >> ${batchfile}
				echo "subdir=${subdir}" >> ${batchfile}
				echo '  ./run_impurity_dos.sh /home/software/matlabR2012a-64/ $parameter1 $parameter2 $parameter3 > ./${subdir}/$parameter1$parameter2$parameter3.out 2>&1' >> ${batchfile}
				releasecommand='scontrol release <job_id>'
				template=${batchfile}
				break;;
			[Hh]* ) 
				batchcommand='qsub'
				batchfile='script.pbs'
				holdflag='-h'
				echo "#! /bin/bash" > ${batchfile}
				echo "#PBS -N impurity_dos" >> ${batchfile}
				echo "#PBS -o nout.out" >> ${batchfile}
				echo "#PBS -e error.err" >> ${batchfile}
				#echo "#PBS -M kreisel@phys.ufl.edu" >> ${batchfile}
				echo "#PBS -r n" >> ${batchfile}
				echo "#PBS -l walltime=12:00:00" >> ${batchfile}
				echo "#PBS -l nodes=1:ppn=1" >> ${batchfile}
				echo "#PBS -l pmem=3500mb" >> ${batchfile}
				echo 'cd $PBS_O_WORKDIR' >> ${batchfile}
				echo "  parameter1=input_10Band_fese_000GP_exp_ce.dat" >> ${batchfile}
				echo "  parameter2=" >> ${batchfile}
				echo "  parameter3=" >> ${batchfile}
				echo "subdir=${subdir}" >> ${batchfile}
				echo "module load matlab/2013a" >> ${batchfile}
				echo './run_impurity_dos.sh ${MATLAB} $parameter1 $parameter2 $parameter3 > ./${subdir}/$parameter1$parameter2$parameter3.out 2>&1' >> ${batchfile}
				template=${batchfile}
				releasecommand='qrls  <job_identifier>'
				break;;
				* ) echo "Please answer correctly.";;
			esac
		done
	fi

#executable=run_impurity_dos_p.sh
# create subdirectory for the output files
#case0=`./get_input.sh $1 "Case"`
#sub_dir='kmeshparallel_'$2
if [ ! -d "$sub_dir" ]
then                    
	echo "Creating subdirectory $subdir."
	mkdir $subdir
fi                                                                                                      
if [ ! -d "$subdir" ]
then                    
	echo "Error creating subdirectory."
	exit
fi

# to do: copy input files into subdirectory



# create number_of_tasks inputfiles in subdirectory
for (( task=1; task<=num_tasks+1; task++ ))
do
	#submit_script=${template}'_'${inputfile}'_task_'${task}
	# do not store all submit scripts (only last one that is set on hold)
	submit_script=${template}'_'${inputfile}
	cp $template $submit_script
	replace_input $submit_script parameter1 ${inputfile} space
	replace_input $submit_script parameter3 ${task} space
	replace_input $submit_script parameter2 ${num_tasks} space
	## check whether chi0 is already calculated if not submit a job
	if [ "$task" -le "$num_tasks" ]
	 then
		 echo "Submitting job..."
		 sleep 1
	 	 ${batchcommand} ${submit_script}
	 else 
		 #echo "Please submit one job ${batchcommand} ${submit_script} when precalculation has been finished."
		 ${batchcommand} ${holdflag} ${submit_script}
		 echo "The last job is on hold, please release when precalculation has been finished."
		 echo ${releasecommand}
		 # do the actual submitting with holding
 	 fi
done
