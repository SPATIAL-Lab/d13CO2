#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --account=bowen-np
#SBATCH --partition=bowen-np
#SBATCH --job-name=d13CO2_definitive
#SBATCH --array=0-6
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dustin.t.harper@utah.edu

profiles=(main gmst_scotese plate_torsvik2017 plate_merdith2021 plate_cao2024 cenozoic coupled)

export WORK_DIR="/uufs/chpc.utah.edu/common/home/bowen-group3/dth/d13CO2"
export RUN_PROFILE="${profiles[$SLURM_ARRAY_TASK_ID]}"
export D13CO2_MODEL_RUN="${D13CO2_MODEL_RUN:-definitive}"
export R_SCRIPT_PATH="${WORK_DIR}/R/model/d13CO2_RunAll.R"
export R_LIBS_USER="/uufs/chpc.utah.edu/common/home/u6041364/R/x86_64-redhat-linux-gnu-library/4.5"

module load gcc/8.5.0
module load jags
module load R

cd "$WORK_DIR"
echo "Starting ${D13CO2_MODEL_RUN}/${RUN_PROFILE} at $(date)"
Rscript "$R_SCRIPT_PATH"
echo "Finished ${D13CO2_MODEL_RUN}/${RUN_PROFILE} at $(date)"
