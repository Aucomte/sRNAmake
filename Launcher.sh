#!/bin/bash
#rm snakejob.*
module purge
module load system/python/3.7.3 #system/singularity/3.4.2

cluster_config="./cluster_config.yaml"
datas_config=."/config.yaml"

# produit le graph du pipeline
#snakemake --dryrun --configfile ${datas_config} -s sRNA_pipeline.smk
#exit
#
#snakemake --dag --configfile ${datas_config} -s sRNA_pipeline.smk  | dot -Tpdf > schema_pipeline_dag.pdf
#exit
#
#snakemake -s sRNA_pipeline.smk --configfile ${datas_config} --rulegraph | dot -Tpdf > schema_pipeline_global.pdf
#


snakemake --latency-wait 5184000 -s sRNA_pipeline.snake --jobs 100 --cluster "qsub {cluster.queue} {cluster.export_env} {cluster.cwd} \
{cluster.mem} {cluster.n_cpu} {cluster.logerror} {cluster.log} " --cluster-config ${cluster_config} --configfile ${datas_config}

#snakemake -s sRNA_pipeline.smk --configfile ${datas_config} --filegraph | dot -Tpdf > schema_pipeline_files.pdf
#snakemake -s sRNA_pipeline.smk --configfile ${datas_config} --dag | dot -Tpdf > schema_pipeline_samples.pdf
#snakemake -s sRNA_pipeline.smk --report REPORT.html

# add to clean files
#snakemake -s sRNA_pipeline.smk clean
