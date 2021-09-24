#!/bin/sh

jobName=""
scriptPrefix=""
mem=1.8
priority=0
queue=""
arguments=""

while [ "$1" != "" ]; 
do
	case $1 in
		--work_dir )				shift
									workDir=$1
									;;
		--job_name )				shift
									jobName=$1
									;;
		--script_prefix )			shift
									scriptPrefix=$1
									;;
		--mem )						shift
									mem=$1
									;;
		--priority )				shift
									priority=$1
									;;
		--cpu )	        			shift
									cpu=$1
									;;
		--queue )					shift
									queue="-q "$1
									;;
		--program )					shift
									program=$1
									;;
		* )							arguments=$arguments" "$1
									;;
	esac
	shift
done

fileScript=$scriptPrefix".qsub.sh"
fileScriptLog=$scriptPrefix".qsub.log"
fileScriptError=$scriptPrefix".qsub.error"

echo \# path list > $fileScript
echo \#!/bin/bash >> $fileScript
echo \# >> $fileScript
echo \# resource list  >> $fileScript
echo \#PBS -l walltime=600:00:00  >> $fileScript
echo \# >> $fileScript
echo \# job name  >> $fileScript
echo \#PBS -N $jobName >> $fileScript
echo \# >> $fileScript
echo >> $fileScript
echo cd $workDir >> $fileScript
echo $program $arguments >> $fileScript

chmod 775 $fileScript

qsub -l mem=${mem}gb,ncpus=${cpu} -p $priority $queue -o $fileScriptLog -e $fileScriptError $fileScript
