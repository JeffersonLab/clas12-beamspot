#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=beamspot_rgc
#SBATCH --output=./log/%x-%j-%N.out
#SBATCH --error=./log/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --time=24:00:00
#SBATCH --array=0-50

source /group/clas12/packages/setup.sh
module load clas12/pro

#FILES=(/volatile/clas12/rg-c/production/pass0v0-rf-raster/calib/recon/016066/rec_clas_016066.evio.0*hipo)
FILES=(/volatile/clas12/rg-c/production/pass0v0-beamoffsets/calib/recon/016044/rec_clas_016044.evio.00*hipo)
srun steer.sh $SLURM_ARRAY_TASK_ID ${FILES[@]:$(($SLURM_ARRAY_TASK_ID*20)):20}
