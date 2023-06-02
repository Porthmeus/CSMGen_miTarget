#!/bin/bash
#SBATCH --job-name=testJob
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=48:00:00
#SBATCH --output=testJob.out
#SBATCH --partition=all


echo "Starting matlab" &>> logs/extractConsistent_Matjes.colormore22.log; module load matlab/2018b; module load glpk; matlab -nodisplay -nodesktop -nosplash -nojvm -r "try; addpath(genpath('workflow/scripts')); extractConsistentModel('results/data/CBT_models/colormore22.mat', 'resources/diets/Matjes.csv', 'results/data/consistentModels/Matjes.colormore22_Consistent.csv','$HOME/Programs/CobraToolbox/'); exit; catch ME; disp(ME); exit; end;" &>>logs/extractConsistent_Matjes.colormore22.log
