#!/usr/bin/env python3
# written by Bao Tram Vi, modified by julie Orjuela (IRD)
from snakemake.logging import logger

#logger.info("wrapper to launch CulebrONT on slurm HPC")

import os
import sys
from snakemake.utils import read_job_properties
from snakemake import load_configfile

jobscript = sys.argv[-1]
config = sys.argv[1]
cluster_config = sys.argv[2]

#read_job_propeties def reads the job properties defined in a snakemake jobscript and return a dict containing information about the job
job_properties = read_job_properties(jobscript)
config_properties = load_configfile(config)
cluster_properties = load_configfile(cluster_config)

rule = job_properties['rule']
jobid = job_properties['jobid']
#cluster = job_properties['cluster']
log = rule

#logger.info("cluster properties : ")
#logger.info(cluster_properties)

#"cluster": {"cpus-per-task": 4, "ntasks": 1, "mem-per-cpu": "2", "partition": "normal", "output": "logs/stdout/run_flye/fastq=5percentB1-1", "error": "logs/error/run_flye/fastq=5percentB1-1"}}

# recovery wildcards in variables
try: 
    references = job_properties['wildcards']['references']
    log += '_{}'.format(references)
except (AttributeError, KeyError):
    pass

try:
    fastq_name = job_properties['wildcards']['fastq_name']
    log += '_{}'.format(fastq_name)
except (AttributeError, KeyError):
    pass

try:
    treatment = job_properties['wildcards']['treatment']
    log += '_{}'.format(treatment)
except (AttributeError, KeyError):
    pass

try:
    ext = job_properties['wildcards']['ext']
    log += '_{}'.format(ext)
except (AttributeError, KeyError):
    pass

# recovery cpu per task from dict properties
#cpus_per_task = job_properties['threads']
outdir = config_properties['DATA']['directories']['out_dir']
logdir = os.path.join(outdir, "slurm_log")
os.makedirs(logdir, exist_ok=True)

#logger.info("cluster properties partition : ")
#if rule in cluster_properties :
#    logger.info("cluster properties partition : ")
#    logger.info(cluster_properties[rule]['partition'])
#    logger.info("cluster properties mem-per-cpu : ")
#    logger.info(cluster_properties[rule]['mem-per-cpu'])

# getting ressources from cluster config given by user
def get_ressources(rule):
    """
    define ressources from cluster_config file or get params default for rule
    """
    if rule in cluster_properties and 'partition' in cluster_properties[rule]:
       queue = cluster_properties[rule]['partition']
       mempercpu = cluster_properties[rule]['mem-per-cpu']
       cpus = cluster_properties[rule]['cpus-per-task']
       return f"--partition {queue} --mem-per-cpu {mempercpu}G --cpus-per-task {cpus} "
    #elif '__default__' in cluster_properties and 'partition' in cluster_properties['__default__']:
    else:
       queue = cluster_properties['__default__']['partition']
       mempercpu = cluster_properties['__default__']['mem-per-cpu']
       cpus = cluster_properties['__default__']['cpus-per-task']
       return f'--partition {queue} --mem-per-cpu {mempercpu}G --cpus-per-task {cpus}'
    
partition = get_ressources(rule)

sbatch = f'sbatch --parsable --job-name {rule} {partition} --ntasks 1 --output {logdir}/{log}.log_%j --error {logdir}/{log}.log_%j'

# read jobscript and insert information such as node, used and memory
with open(jobscript, "r") as j:
    scripts = j.readlines()

scripts.insert(1, "echo -e \"# sbatch parameters: \"{}\"\"\n".format(sbatch))
scripts.insert(2, "echo -e \"# Job running on node: $SLURM_JOB_NODELIST\"\n")
scripts.insert(3, "echo -e \"# Number of used CPUS: $SLURM_CPUS_ON_NODE\"\n")
scripts.insert(4, "echo -e \"# Memory per CPU in megabyte: $SLURM_MEM_PER_CPU\"\n")
scripts.insert(5, "echo -e \"# Partition: $SLURM_JOB_PARTITION\"\n")

with open(jobscript, "w") as j:
    j.writelines(scripts)

cmdline = " ".join([sbatch, jobscript])
logger.info(f'INFO : {cmdline}')

os.system(cmdline)

