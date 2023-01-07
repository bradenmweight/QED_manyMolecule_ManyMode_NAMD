#!/bin/bash

cd TRAJ/
for t in {148..499}; do
    echo "Submitting traj-${t}"
    cd traj-${t};
        cp ../../NAMD.in .
        cp ../../submit.NAMD .
        sbatch submit.NAMD
        sleep 0.5
        cd ../
done

