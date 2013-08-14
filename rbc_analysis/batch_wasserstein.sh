#!/bin/sh

#cd /sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/rbc_analysis

var=`cat cellnames`

for line in $var; 
do
    echo $line;
    qsub -v cell=$line,lag=10 -l nodes=1:x5672:ppn=1 -l walltime=01:00:00 run_wasserstein.sh;
done

