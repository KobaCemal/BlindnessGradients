#!/bin/bash -l 
## task_name 
#SBATCH -J test 
## num nodes 
#SBATCH -N 1 
## tasks per node 
#SBATCH --ntasks-per-node=1
## cpus per taskq 
#SBATCH --cpus-per-task=16
## memory allocated per cpu 
#SBATCH --mem-per-cpu=6GB 
## max time 
#SBATCH --time=05:00:00 
## grant name 
#SBATCH -A plgsano5-cpu 
## partitionq 
#SBATCH --partition plgrid
## gpus allocation 
##SBATCH --gpus-per-node=4
## output file path 
#SBATCH --output="/net/ascratch/people/plgkoba/output_bal.out" 
## error file path 
#SBATCH --error="/net/ascratch/people/plgkoba/error_balt.err" 
## loading modules that are not in the environment 
## module add cudnn/8.4.1.50-cuda-11.6.0 
## bash commands to execute 



# for scan in $(cat scan_quality_batch1.txt )
# do 
# sub=$(echo $scan | cut -d "_" -f 1 | cut -d "-" -f2)
# ses=$(echo $scan | cut -d "_" -f 2 | cut -d "-" -f2)
 
# singularity run -B $SCRATCH  $SCRATCH/micalab_micapipe_v0.2.3-2024-01-18-4ee63843ce22.simg \
#     -bids $SCRATCH/datasets/baltimore \
#     -out $SCRATCH/datasets/baltimore/derivatives \
#     -tmpDir $SCRATCH/tmp \
#     -ses $ses \
#     -threads 16\
#     -fs_licence $SCRATCH/license.txt \
#     -sub $sub  -proc_structural -proc_surf -post_structural -atlas glasser-360 -GD -proc_func -mainScanStr task-rest_bold -noFIX -NSR -QC_subj
# done
# done 
#  

for scan in "sub-016"
do 
sub=$(echo $scan | cut -d "_" -f 1 | cut -d "-" -f2)
ses=$(echo $scan | cut -d "_" -f 2 | cut -d "-" -f2)
 
singularity run -B $SCRATCH  $SCRATCH/micalab_micapipe_v0.2.3-2024-01-18-4ee63843ce22.simg \
    -bids $SCRATCH/datasets/baltimore \
    -out $SCRATCH/datasets/baltimore/derivatives \
    -tmpDir $SCRATCH/tmp \
    -ses $ses \
    -threads 16\
    -fs_licence $SCRATCH/license.txt \
    -sub $sub  -proc_structural -proc_surf -post_structural -atlas glasser-360 -GD -proc_func -mainScanStr task-rest_bold -noFIX -NSR -QC_subj
done
 

