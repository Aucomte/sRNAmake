#!/bin/bash
#SBATCH --job-name sRNA_pipeline
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#module purge

# cding into the path of the relvent pipeline instance
#cd /data/projects/rymv_xantho_smallrna/seb/pipeline_srnas

# Snakefile path
snkf="sRNA_pipeline.smk"

# run config file
#cfgf="config/config_BAI3.yaml"
cfgf="config/config.yaml"

# starting environments for the pipeline
module load system/Miniconda3/1.0
module load system/singularity/3.3.0
env=/home/$(whoami)/.conda/envs/snakemake
[ ! -d $env ] && echo -e "## [$(date) - sRNA_pipeline]\t Creating conda environment for snakemake" && conda env create -f envs/environment.yaml -n snakemake
source activate snakemake


# produit le graph du pipeline
#snakemake -s "$snkf" --configfile "$cfgf" --dag | dot -Tpdf > dag.pdf
#snakemake -s "$snkf" --configfile "$cfgf" --rulegraph | dot -Tpdf > rules.pdf

#snakemake --unlock --cores 1

# Run pipeline
echo -e "## [$(date) - sRNA_pipeline]\t Launching snakemake pipeline"
snakemake --nolock  --use-conda --use-singularity \
  --singularity-args "--bind /scratch:/tmp" \
  -s "$snkf" \
  --jobs 24 \
  --configfile "$cfgf" \
  -p --verbose \
  --latency-wait 60 --keep-going --restart-times 1 --rerun-incomplete  \
  --cluster "python3 script/slurm_wrapper.py "$cfgf" config/cluster_config_slurm.yaml" \
  --cluster-config "config/cluster_config_slurm.yaml" \
  --cluster-status "python3 script/slurm_status.py"

#  --local-cores $SLURM_CPUS_PER_TASK \

echo -e "## [$(date) - sRNA_pipeline]\t Creating snakemake report"
python3 script/snakemake_report.py -s "$snkf" --configfile "$cfgf"
