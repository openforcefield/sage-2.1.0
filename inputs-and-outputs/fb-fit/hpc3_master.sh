#!/bin/bash
#SBATCH -J lotr
#SBATCH -p standard
#SBATCH -t 200:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000mb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mail-user=dir@uci.edu
#SBATCH --constraint=fastscratch

rm -rf /tmp/$SLURM_JOB_NAME
source $HOME/.bashrc
eval "$(/opt/apps/anaconda/2020.07/bin/conda shell.bash hook)"

export SLURM_TMPDIR=/tmp
export TMPDIR=$SLURM_TMPDIR/$SLURM_JOB_NAME
mkdir  -p  $SLURM_TMPDIR/$SLURM_JOB_NAME
rsync -avzIi /dir/dir/dir/conda-env/fb_193.tar.gz $SLURM_TMPDIR/$SLURM_JOB_NAME
mkdir  -p  $SLURM_TMPDIR/$SLURM_JOB_NAME/fb_193
cd $SLURM_TMPDIR/$SLURM_JOB_NAME
tar xf fb_193.tar.gz -C fb_193
source fb_193/bin/activate
echo $(python -V)

rsync  -avzIi  $SLURM_SUBMIT_DIR/optimize.in  $SLURM_TMPDIR/$SLURM_JOB_NAME
rsync  -avzIi  $SLURM_SUBMIT_DIR/targets.tar.gz  $SLURM_TMPDIR/$SLURM_JOB_NAME
rsync  -avzIi  $SLURM_SUBMIT_DIR/forcefield  $SLURM_TMPDIR/$SLURM_JOB_NAME

tar -xzf targets.tar.gz

datadir=$(pwd)
mkdir -p $SLURM_SUBMIT_DIR/worker_logs
echo $(hostname) > $SLURM_SUBMIT_DIR/host

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

if ForceBalance.py optimize.in ; then
   tar -czf optimize.tmp.tar.gz optimize.tmp
   tar -czf result.tar.gz result
   mkdir -p ~/fit9/$SLURM_JOB_ID
   cp -rf result ~/fit9/$SLURM_JOB_ID
   find ./ -type d -exec chgrp dmobley_lab_share {} \;
   find ./ -type f -exec chgrp dmobley_lab_share {} \;
   find ./ -type d -exec chmod g+s {} \;
   find ./ -type f -exec chmod g+s {} \;
   rsync  -avzIi --no-o --no-g --recursive --exclude="optimize.tmp" --exclude="optimize.bak" --exclude="fb_193*" --exclude="targets*" $SLURM_TMPDIR/$SLURM_JOB_NAME/* $SLURM_SUBMIT_DIR --rsync-path="sudo rsync" > copy.log
   #rm -rf $SLURM_TMPDIR/$SLURM_JOB_NAME
fi
sleep 3h 
echo "All done"
