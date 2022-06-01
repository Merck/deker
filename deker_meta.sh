#!/bin/bash
set -e
set -u
set -o pipefail
############################################################################
usage() {
    echo "Usage: ./$(basename $0) [ options ] [ control_file ]"
}
help() {
    echo "
Summary: Script for deker network inference on SGE cluster "
    usage
    
    echo "
Options:
  -h                      Help.
  -v                      Verbose; writes a log file 

Arguments:

  control_file         DATA_FILE       Control file to use; plain text file with at minimum a line 'data_file [path/to/file]'. 

When using, please cite the following contributors:
Sean Hayes
"
}
############################################################################
############################################################################
# Safety - create a lock directory to prevent multiple instances
############################################################################
############################################################################
lock_dir=deker_meta_lock
if ! mkdir $lock_dir; then
  echo "script locked, can't run"
  exit
fi 
############################################################################
############################################################################
# Read in optional arguments
############################################################################
############################################################################
#Defaults
VERBOSE=false
############################################################################
while getopts 'hvcr:' OPTION; do
  case "$OPTION" in
    h)
      help
      exit
      ;;
    v)
      VERBOSE=true
      ;;
    ?)
      usage
      exit 1
      ;;
  esac
done
echo "Starting..."
############################################################################
############################################################################
# Read control file arguments
############################################################################
############################################################################
shift $(($OPTIND - 1))
CONTROL_FILE=$1
if [ ! -f "$CONTROL_FILE" ]; then
  echo "Control file not found"
  exit
fi
############################################################################
############################################################################
# Final pre-run setup
############################################################################
############################################################################
dt_id=$(date '+%S%M%H%d%m%Y')
digest_file="deker.$dt_id/digest"
edgelist_file="edge_out_$dt_id.csv"
n_responses=$(./deker -s $CONTROL_FILE $digest_file)
echo "Running feature selection for $n_responses features in $1; key $dt_id"
#if verbose, set a log_file
if [ "$VERBOSE" = true ]; then
  log_file="$dt_id.control_log"
fi
#make intermediate directories
mkdir -p logs
############################################################################
############################################################################
# Initial Parameter ID Submit Loop
############################################################################
############################################################################
for ((i=0;i<n_responses;i++)); do
  n_iters[$i]=0
  skip_list[$i]=false
  job_slots[$i]=-1 
done
############################################################################
############################################################################
# Submit & Monitoring Loop
############################################################################
############################################################################
max_pending_jobs=2000
check_pending_jobs=1000
####
rm -f $dt_id.ind_out/*
exit_while=false
while [ "$exit_while" = false ]; do
  ###############################################################################
  if [ ! -z ${log_file+x} ]; then
	  echo "Checking for completed jobs" >> $log_file
	fi
	#get the job ids of running jobs in qstat
	set +e
	jobs_on_grid=$(qstat -s rp | grep -o '^.*deker' | grep -o '^[[:space:]]*[[:digit:]]*' | grep -o '[[:digit:]]*' | uniq | paste -s -d " ")
	n_pending_jobs=$(qstat -s p | grep -c '^.*deker')
	set -e
	if [ ! -z ${log_file+x} ]; then
	  echo $jobs_on_grid >> $log_file
	fi
	if [ $n_pending_jobs -lt $check_pending_jobs ]; then
  	exit_while=true
  	#############################################################################
    for ((i=0;i<n_responses;i++)); do
  	  if [ "${skip_list[$i]}" = false ]; then
  	    exit_while=false
  	    #############################################################################
  	    if [ ! -z ${log_file+x} ]; then
  	      echo "Checking feature $i" >> $log_file
  	      echo "Slot ${job_slots[$i]}" >> $log_file
  	    fi
  		#############################################################################
  		#check if saved job for given feature is found in qstat 
  		  job_complete=true
  		  for j in $jobs_on_grid; do 
  			  if [ "${job_slots[$i]}" = "$j" ]
  			  then
  			    job_complete=false
  			    exit_while=false
  			    break
  			  fi
  		  done
  	    #############################################################################
  		#submit next set of needed iterations if previous set is done
  	    if [ "$job_complete" = true ]; then
      	##########################################################################################
          check_status=$(./deker -c $i $digest_file $edgelist_file)
          if [ $check_status -eq "1" ]; then
	          skip_list[$i]=true
	          continue
	      fi
          ####
    	    set +e
    	    n_pending_jobs=$(qstat -s p | grep -c '^.*deker')
    	    set -e
    	    new_pending_jobs=$(($n_pending_jobs+1))
	        if [ $new_pending_jobs -gt $max_pending_jobs ]; then
	          break
	        fi
          ####
      	  ((++n_iters[$i]))
      	  job_slots[$i]=$(qsub -terse ./deker_qsub.sh $i $digest_file | grep -o "^[[:digit:]]*")
      	  if [ ! -z ${log_file+x} ]; then
      	    echo "Submitting $i" >> $log_file
      	  fi
      	  sleep .05 #throttle submission a little so we don't load in too many pending before they're able to start running
    	  ##########################################################################################
  	    fi #End if job complete
  	  fi #End if not skipped
    done #End feature loop
  fi #End if pending check
  #Wait before checking for more finished jobs
  sleep 30
done #End submit loop
rm -rf $lock_dir
exit