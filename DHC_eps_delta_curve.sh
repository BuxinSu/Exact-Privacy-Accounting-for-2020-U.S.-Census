#!/bin/bash
#SBATCH --account=bingxin
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=04:00:00
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --array=1-43
#SBATCH --job-name=default
#SBATCH --output=/scratch/bell/bingxin/EHR_model/test_code/logs_file/DHC_eps_delta_curve_binary_search_%a.out
#SBATCH --chdir=/depot/bingxin/data/env/data_test

module purge
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/bell/bingxin/EHR_model/envs/ehr311

python /scratch/bell/bingxin/EHR_model/test_code/DHC_eps_delta_curve_binary_search.py --k "${SLURM_ARRAY_TASK_ID}"
