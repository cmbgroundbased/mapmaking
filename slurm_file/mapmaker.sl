#!/bin/bash -l
#
#SBATCH --job-name=mapmaking
#SBATCH --account=strip
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=20000
#SBATCH --time=2-00:00:00
#SBATCH --exclusive
#SBATCH --error=job_log/job.%J.err 
#SBATCH --output=job_log/job.%J.out

cd $HOME/mapmaking

cargo clean
cargo build
cargo run -- -t /home/projects/Strip/ATMMC/jan/

echo "SEFF REPORT:"
seff -d $SLURM_JOBID

