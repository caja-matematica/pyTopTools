#!/bin/sh

#cd /sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/rbc_analysis

var=`cat cellnames2`

for line in $var; 
do
    echo $line;
    qsub -v cell=$line -l nodes=1:x5672:ppn=1 -l walltime=4:00:00 run_rbc_persistence.sh;
done

