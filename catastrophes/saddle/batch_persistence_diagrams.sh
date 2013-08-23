#!/bin/csh 

cd /sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/catastrophes/saddle

# run 50 realizations;  trial start=1 and end=50
foreach trial (`seq 1 30`) 
    foreach window ( `seq 10 10 50` )
	qsub -v run=$trial,win=$window -l nodes=1:x5672:ppn=1 -l walltime=20:00:00 saddle_node_persistence.sh;
    end
end








