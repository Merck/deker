#!/bin/bash
#
#$ -cwd -V -b y
#$ -P long
#$ -N deker
#$ -l mem_reserve=10G -l h_vmem=40G -l fourperhost=1
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.o -e logs/$JOB_NAME.$JOB_ID.$TASK_ID.e

# $1 : which_response
echo "which_response $1"
# $2 : digest_dir
echo "digest_dir $2"

./deker -f $1 $2