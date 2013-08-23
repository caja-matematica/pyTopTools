#!/bin/csh 

cd $HOME/github/local/caja-matematica/pyTopTools/catastrophes/hopf

# run 50 realizations;  trial start=1 and end=50
foreach trial (`seq 1 50`) 
    foreach window (`seq 50 50 300`)
	echo $trial
	echo $window
	qsub -v run=$trial,win=$window -l nodes=1:x5672:ppn=1 -l walltime=10:00:00 persistence_hopfsub;
    end
end








