#!/bin/csh 

cd /sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/catastrophes

# run 50 realizations;  trial start=1 and end=50
foreach trial (`seq 1 30`) 
    qsub -v run=$trial -l nodes=1:x5672:ppn=1 -l walltime=20:00:00 saddle_node_run.sh;
end








