#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem=7G

word=$1

jobflag='-j 20'
kgflag=''
flags=''
drmaaflag=''
usedrmaa='false'

if [ "x$word" == 'xdry' ]; then
    flags="--dryrun --quiet" # --reason"
elif [ "x$word" == 'xunlock' ]; then
    flags='--unlock'
elif [ "x$word" == 'xtest' ]; then
    jobflag='-j 1'
    kgflag=''
elif [ "x$word" == 'xlocal' ]; then
    jobflag='-j 20'
    kgflag='--keep-going'
else
    usedrmaa='true'
    jobflag='-j 1000'
    kgflag='--keep-going'
    flags="--max-status-checks-per-second 0.1" # --restart-times 2"
fi

#module load slurm-drmaa

#export TMPDIR=$(realpath try3/tmp)
#echo "TMPDIR=$TMPDIR"
# tmpdir is necessary for using workflow.source to properly find
# scripts from within module calls.
# XXX: doesn't seem to be the case anymore?



# I can't get $drmaaflags to substitute properly because of the internal 's
if [ $usedrmaa == "true" ]; then
    snakemake $flags \
        --dir . \
        --latency-wait 60 $kgflag \
        --rerun-incomplete \
        -s snakemake/Snakefile \
        --max-threads 12 $jobflag \
        --drmaa ' -p priopark -A park_contrib --mem={resources.mem} -c {threads} -t 24:00:00 -o cluster-logs/slurm-%A.log'
else
    snakemake $flags \
        --dir . \
        --latency-wait 60 $kgflag \
        --rerun-incomplete \
        -s snakemake/Snakefile \
        -s snakemake/Snakefile \
        --max-threads 12 $jobflag \
        #--max-inventory-time 0 \
fi
