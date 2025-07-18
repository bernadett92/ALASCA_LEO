#!/bin/bash

set -euo pipefail

module load harris/idl/88
module load nvidia/cuda/11.3

savedir='/home2/bstadler/ALASCA_LEO/'
params_file='params_alasca_leo_sh.pro'
atmos=('Tenerife' 'TURBO50')
seeds=(1)
#simuls=('noAO' 'TT' 'SATAO' 'LGSAO')
simuls=('noAO')
simul_time=10.0
aberrations=(0.0 160.3751 320.7501 481.1252)

for atmo in "${atmos[@]}"; do
    echo "$atmo"
    for seed in "${seeds[@]}"; do
        echo "$seed"
        for simul in "${simuls[@]}"; do
            echo "$simul"
            for aberration in "${aberrations[@]}"; do
                echo "$aberration"
                export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}
                time PYTHONPATH=/home/bstadler/Simulator/Durham/aotools/ idl -e "run_alasca_simulations, '$atmo','$params_file','$savedir',$seed,'$simul',$simul_time, $aberration"
            done
        done
    done
done
