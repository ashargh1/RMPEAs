#!/bin/bash
#SBATCH --ntasks-per-node 24
#SBATCH --mem-per-cpu 6G
#SBATCH --time 24:00:00
#SBATCH --partition=parallel 
#SBATCH -A jelawad1-condo
#SBATCH --job-name jupyter-notebook
#SBATCH --output jupyter-notebook-%J.log

ml anaconda
## or use your own python/conda enviromnent
#XDG_RUNTIME_DIR=””
#port=$(shuf -i8000-9999 -n1)
#echo $port
#node=$(hostname -s)
#user=$(whoami)
#jupyter-notebook –no-browser –port=${port} –ip=${node}

jupyter nbconvert --ExecutePreprocessor.timeout=172800 --to notebook --execute Untitled.ipynb
