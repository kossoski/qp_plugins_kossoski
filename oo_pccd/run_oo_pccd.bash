#!/bin/bash

# This is a simple script that calls the oo_pccd program of Quantum Package.
# Each call represents one iteration in the orbital optimization procedure.
# When convergence is reached the file "converged" is written, and if the program crashes the file "failed" is written. In either case this script stops running.
# In the end, the file e_pccd_last contains the energy of the last iteration.

# Please modify the path to your Quantum Package directory, the name of the EZFIO file, and the initial (i) and final (f) indexes of the iteration.
###########################################################
source ~/qp2/quantum_package.rc
ezfio_file="ch+"
i=1
f=200
###########################################################

converged="converged"
if [ -f "$converged" ]; then
	rm $converged
fi
failed="failed"
if [ -f "$failed" ]; then
	rm $failed
fi

for (( j=$i; j<=$f; j++ ))
do
	if [ -f "$converged" ]; then
		echo "Converged..."
		exit
	fi	
	if [ -f "$failed" ]; then
		echo "Convergence failed..."
		exit
	fi	
	echo "$j"

	qp_run oo_pccd $ezfio_file > pcc_$j.out

	if grep -q 'Program exited' pcc_$j.out; then
  	  touch $failed
	fi

        grep EpCCD pcc_$j.out | sed 's/EpCCD//g' > e_pccd_last
        cp   pcc_$j.out pcc_last.out

done
