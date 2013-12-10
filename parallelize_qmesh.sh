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
	echo "Usage parallize_qmesh.sh parameterfile number_of_tasks submit_script_template"
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

if [ ! -e "$3" ]
then
	echo "Error submit_script_template $3 does not exist, using standart ones:"
	while true;
       	do
		read -p "[I]TP (Slurm), [H]PC-script " ih
		case $ih in
			[iI]* ) 
				echo "#! /bin/bash" > chi0_itp.sh
				echo "#SBATCH -p dfg" >> chi0_itp.sh
				echo "#SBATCH -n 1" >> chi0_itp.sh
				echo "#SBATCH --mem-per-cpu=3800" >> chi0_itp.sh
				echo "  parameter1=standart_input_imp_dos_M_20.mat" >> chi0_itp.sh
				echo "  parameter2=" >> chi0_itp.sh
				echo "  parameter3=" >> chi0_itp.sh
				echo '  run_impurity_dos.sh /home/software/matlabR2012a-64/ $parameter1 $parameter2 $parameter3 > $parameter1$parameter2$parameter3.out 2>&1' >> chi0_itp.sh
				template=chi0_itp.sh
				break;;
			[Hh]* ) 
				echo "not implemented yet"
				exit;;
			* ) echo "Please answer correctly.";;
		esac
	done
fi

batchcommand=sbatch
executable=run_impurity_dos_p.sh
# create subdirectory number_of_tasks_Casename
#case0=`./get_input.sh $1 "Case"`
sub_dir='kmeshparallel_'$2
if [ ! -d "$sub_dir" ]        
then                    
	echo "Creating subdirectory $sub_dir."
	mkdir $sub_dir
fi                                                                                                      
if [ ! -d "$sub_dir" ]
then                    
	echo "Error"
	exit
fi

# to do: copy input files into subdirectory



# create number_of_tasks inputfiles in subdirectory
for (( task=1; task<=num_tasks+1; task++ ))
do
	submit_script=${template}'_task_'${task}
	cp $template $submit_script
	replace_input $submit_script parameter1 ${inputfile} space
	replace_input $submit_script parameter3 ${task} space
	replace_input $submit_script parameter2 ${num_tasks} space
	## check whether chi0 is already calculated if not submit a job
	if [ "$task" -le "$num_tasks" ]
	 then
		 echo "Submitting job..."
		 sleep 2
	 	 ${batchcommand} ${submit_script}
	 else 
		 echo "Please submit one job ${batchcommand} ${submit_script} when precalculation has been finished."
 	 fi
done


