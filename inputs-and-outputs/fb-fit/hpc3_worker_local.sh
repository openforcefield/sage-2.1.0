#!/bin/bash

host=$1
port=$2
shift
shift

cmd=$(mktemp)
cat << EOF > $cmd
#!/usr/bin/env bash
#SBATCH -J wq-$port
#SBATCH -p free
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
# SBATCH --array=1-100
#SBATCH --account lab
# SBATCH --export ALL
#SBATCH -o /dev/null

#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
if ! diff -q /tmp/dir/fb_193.tar.gz /dir/dir/dir/conda-env/fb_193.tar.gz > /dev/null 2>&1;
then 
   rm -rf /tmp/dir/fb_*
fi

mkdir /tmp/dir -p
for i in \$(seq  \$SLURM_NTASKS ); do
        echo $i
        ./wq_worker_local.sh --cores 1 -s /tmp/dir --disk-threshold=0.002 --disk=3000 --memory-threshold=1000 -t 3600  -b 20 --memory=1000 $host:$port &
done
wait
EOF

sbatch $@ $cmd 
rm $cmd
