## Snakemake - sRNA_pipeline

import re, sys, glob, os, csv, pprint
import pandas as pd
from pathlib import Path
from snakemake import WorkflowError
#from snakemake import load_configfile
from os import listdir
import getpass

###############################################################################
# NOTE pas de caractere speciaux entre 2 wildcards

## ??????????? I DO NOT KNOW HOW TO DO THAT BUT WE WOULD NEED TO USE THE PATHS
## ??????????? PROVIDED IN THE SNAKEMAKE COMMAND LINE

# --- Importing Configuration Files --- #

#pprint.pprint(workflow.__dict__)
if len(workflow.overwrite_configfiles) == 0:
    logger.info("You need to use --configfile option to snakemake command line")
    raise ValueError("You have to use --configfile option to snakemake command line")
else:
    path_config = workflow.overwrite_configfiles[0]

#configfile: 'config/config.yaml'
cluster_config: "config/cluster_config_slurm.yaml"

###############################################################################

# helper function
def by_cond(cond, yes, no, cond_ext = '', no_ext = ''): # it's working but needs to be improved ...
    if not cond_ext:
        cond_ext = not cond
    if cond:
        return yes
    elif cond_ext:
        return no
    else:
        return no_ext

## cluster variables

#user = getpass.getuser()

#nasID = config['NASID']
#HOST_PREFIX = by_cond(nasID, user + '@' + nasID + ':', '')

##############################

# check configfile:
# existence of dir and validity of suffix

def check_config_file():
    if not os.path.exists(config["DATA"]["directories"]["fastq_dir"]):
        logger.info("CONFIG FILE CHECKING FAIL : in the DATA section, fastq_dir directory does not exists")
        raise ValueError("CONFIG FILE CHECKING FAIL : in the DATA section, fastq_dir directory does not exists")
    if config['DATA']['fasta_suffix'] not in ["fasta", "fna", "fa", "fsa"]:
        raise WorkflowError(f"\n\tERROR : {config['DATA']['fasta_suffix']} is not a good extention for fasta file (fasta, fna, fa, fsa).\nExiting...")

check_config_file()

# parse config file :
fastq_dir = Path(config["DATA"]["directories"]["fastq_dir"]).resolve().as_posix()
references_dir =  Path(config["DATA"]["directories"]["references_dir"]).resolve().as_posix()
out_dir = Path(config["DATA"]["directories"]["out_dir"]).resolve().as_posix()
log_dir = Path(config["DATA"]["directories"]["out_dir"]+"/LOGS/").resolve().as_posix()
miRNA_path = Path(config["DATA"]["files"]["miRNA_gff"]).resolve().as_posix()
samplefile = Path(config["DATA"]["files"]["sample_info"]).resolve().as_posix()
suffix_file = config["DATA"]["fasta_suffix"]

# to lunch separator
sep="#"

#############################################
# use threads define in cluster_config rule or rule default or default in snakefile
#############################################

def get_threads(rule, default):
    """
    use threads define in cluster_config rule or rule default or default in snakefile
    give threads or 'cpus-per-task from cluster_config rule : threads to SGE and cpus-per-task to SLURM
    """
    #cluster_config = load_configfile(cluster_config)
    if rule in cluster_config and 'threads' in cluster_config[rule]:
        return int(cluster_config[rule]['threads'])
    elif rule in cluster_config and 'cpus-per-task' in cluster_config[rule]:
        return int(cluster_config[rule]['cpus-per-task'])
    elif '__default__' in cluster_config and 'cpus-per-task' in cluster_config['__default__']:
        return int(cluster_config['__default__']['cpus-per-task'])
    elif '__default__' in cluster_config and 'threads' in cluster_config['__default__']:
        return int(cluster_config['__default__']['threads'])
    return default

#*###############################################################################
def final_return(wildcards):
    dico_final = {
                     "multiqc_fastqc" : expand(f"{out_dir}/6_MULTIQC/multiqc.html"),
                     "baminfo" : f"{out_dir}/2_mapping_sRNA/bamfile_info.txt",
                     "bam_merged_by_treatments" : expand(f"{out_dir}/3_merge_bam_sRNA/{{treatment}}_merge.bam",treatment=TREATMENT),
                     "ShortStack_gff" : expand(f"{out_dir}/4_ShortStack/ShortStack_All.gff3"),
                     "ShortStack_gff2" : expand(f"{out_dir}/4_ShortStack/ShortStack_All_miRNA.gff3"),
                     "ShortStack_html" : expand(f"{out_dir}/4_ShortStack/shortStack_analysis.html"),
                     #"sRNA_diff_exp_html" : f"{out_dir}/5_sRNA_loci_DE_analysis/sRNA_DE_analysis.html"
                     }
    return dico_final

#*###############################################################################

## ??????????  I AM NOT SURE WE NEED THAT? BECAUSE OF THIS THE ANNOTATION FILES (WHICH MUST BE GTF, WHY?) NEED
## ??????????  ALSO TO HAVE THE SAME NAME PREFIX, IS IT REALLY DESIRABLE?
REFERENCES, = glob_wildcards(f"{references_dir}/{{references}}.{suffix_file}", followlinks=True)

def unique(list1):
    # intilize a null list
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)

def checkIfDuplicates_1(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

samples = {}
header = ""
with open(samplefile, 'r') as f:
    spamreader = csv.reader(f, delimiter = "\t")
    for line in spamreader:
        header = line
        break
    for i in range(len(header)-1):
        with open(samplefile, 'r') as f:
            spamreader = csv.reader(f, delimiter = "\t")
            x = 0
            for line in spamreader:
                if line != header :
                    samples[header[i].lower(),x] = line[i]
                    x = x + 1

# on récupère la liste des traitements et des fastq_name en wildcard
treatments = []
sample_name = []
for key in samples:
    if key[0] == "treatment":
        treatments.append(samples[key])
    if key[0] == "samplename":
        sample_name.append(samples[key])
TREATMENT = unique(list(treatments))

def check_config_file():
    if checkIfDuplicates_1(sample_name):
        logger.info("CONFIG FILE CHECKING FAIL : in the DATA section, fastq_dir directory does not exists")
        raise ValueError("CONFIG FILE CHECKING FAIL : in the DATA section, fastq_dir directory does not exists")

## ????????? DO WE WANT TO RENAME THIS variable/wildcard TO SAMPLE_NAME?
## ????????? IDEALLY WOULD NEED TO CHECK THAT NONE OF THE VALUES OF samplename
## ????????? IN sampleInfo FILE ARE DUPLICATED BECAUSE THIS IS AN IMPLICIT
## ????????? ASSUMPTION AND IF IT HAPPENED IT WOULD BREAK THE PIPELINE
SAMPLE_NAME = unique(list(sample_name))

################################

# Récuperation des fastq_path
def get_fastq(wildcards):
    for key in samples:
        if key[0] == "samplename" and samples[key] == wildcards.fastq:
            row = key[1]
            for key1 in samples:
                if key1[0] == "filename":
                    if row == key1[1]:
                        return f"{samples[key1]}"

##############################
# --- Main Build Rules --- #
rule final:
    """
    construct a table of all resume files
    """
    input:
        unpack(final_return)


# -------------  0 QC:

rule run_Fastqc:
    """
        QC of fastq files on raw fastq files
    """
    threads: get_threads("run_Fastqc", 1)
    input:
        fastq = get_fastq
    output:
        html_fastqc = f"{out_dir}/0_QC/fastqc/{{fastq}}_raw_fastqc.html"
    log:
        error = f"{log_dir}/run_Fastqc/{{fastq}}.e",
        output = f"{log_dir}/run_Fastqc/{{fastq}}.o"
    singularity:
            config["SINGULARITY"]["MAIN"]
    conda:
            "envs/quality.yaml"
    shell:
         """
         fastqc -o {out_dir}/0_QC/fastqc -t {threads} {input.fastq}
         infilename=$(basename {input.fastq})
         fastqcOutFilePath="{out_dir}/0_QC/fastqc/${{infilename%%.*}}"_fastqc.html
         finalOutFilePath="{output.html_fastqc}"
         ## RUN ONLY IF SOURCE AND DEST ARE DIFFERENT PATHS 
         [[ "$fastqcOutFilePath" == "$finalOutFilePath" ]] || mv "$fastqcOutFilePath" "$finalOutFilePath"
         """


## FASTP DOES NOT ALLOW TO FINE TUNE ADAPTER TRIMMING = THIS IS A BIG LIMITATION
## RELATIVE TO OTHER TOOLS eg: https://adapterremoval.readthedocs.io/en/latest/manpage.html


rule run_fastp:
    """
        fastp on reads
    """
    threads: get_threads("run_fastp", 1)
    input:
        fastq = get_fastq
    params:
        options = config["PARAMS"]["FASTP"]["options"]
    output:
        output_html = f"{out_dir}/0_QC/{{fastq}}_fastp.html",
        output_json = f"{out_dir}/0_QC/{{fastq}}_fastp.json",
        fastq_trimmed =  f"{out_dir}/0_QC/{{fastq}}_trimmed.fastq"
    log:
        error = f"{log_dir}/run_fastp/{{fastq}}_fastp.e",
        output = f"{log_dir}/run_fastp/{{fastq}}_fastp.o"
    singularity:
        config["SINGULARITY"]["MAIN"]
    conda:
        "envs/quality.yaml"
    shell:
        """
            fastp --thread {threads} -i {input.fastq} -o {output.fastq_trimmed} -h {output.output_html} -j {output.output_json} {params.options}  1>>{log.output} 2>>{log.error}
        """

rule run_Fastqc_trimmed:
    """
        QC of fastq files on trimmed fastq files
    """
    threads: get_threads("run_Fastqc", 1)
    input:
        fastq = f"{out_dir}/0_QC/{{fastq}}_trimmed.fastq"
    output:
        html_fastqc = f"{out_dir}/0_QC/fastqc/{{fastq}}_trimmed_fastqc.html"
    log:
        error = f"{log_dir}/run_Fastqc/{{fastq}}.e",
        output = f"{log_dir}/run_Fastqc/{{fastq}}.o"
    singularity:
            config["SINGULARITY"]["MAIN"]
    conda:
            "envs/quality.yaml"
    shell:
         """
         fastqc -o {out_dir}/0_QC/fastqc -t {threads} {input.fastq}
         infilename=$(basename {input.fastq})
         fastqcOutFilePath="{out_dir}/0_QC/fastqc/${{infilename%%.*}}"_fastqc.html
         finalOutFilePath="{output.html_fastqc}"
         ## RUN ONLY IF SOURCE AND DEST ARE DIFFERENT PATHS 
         [[ "$fastqcOutFilePath" == "$finalOutFilePath" ]] || mv "$fastqcOutFilePath" "$finalOutFilePath"
         """

rule generate_fastqtrimmed_info:
    """
        generate_fastqtrimmed_info
    """
    threads: get_threads('generate_fastqtrimmed_info', 1)
    input:
        fastp_files = expand(f"{out_dir}/0_QC/{{fastq}}_trimmed.fastq", fastq = SAMPLE_NAME),
        samplefile = samplefile
    params:
        outdir = f"{out_dir}/0_QC/"
    output:
        out_file = f"{out_dir}/0_QC/fastq_trimmed_info.txt"
    script:
        "script/write_fastqtrimmed_info.py"


# -------------  1 preparing data:

rule cat_fasta :
    """
        cat all fasta in references dir
    """
    threads: get_threads('cat_fasta', 1)
    input:
        reference = expand(f"{references_dir}/{{references}}.{suffix_file}", references = REFERENCES)
    output:
        cat_ref = f"{out_dir}/1_cat_ref/all_ref.{suffix_file}"
    log:
         error = f"{log_dir}/cat_ref/all_ref.e",
         output = f"{log_dir}/cat_ref/all_ref.o"
    shell:
         """
            cat {input.reference} > {output.cat_ref}
         """

#rule cat_gtf :
#    """
#        cat all gtf in annotation dir
#    """
#    threads: get_threads('cat_gtf', 1)
#    input:
#        gtf = expand(f"{config['DATA']['directories']['annotation']}{{references}}.gtf", references = REFERENCES)
#    output:
#        cat_gtf = f"{out_dir}/1_cat_gtf/all_gtf.gtf"
#    log:
#        error = f"{log_dir}/cat_gtf/all_gtf.e",
#        output = f"{log_dir}/cat_gtf/all_gtf.o"
#    shell:
#         """
#            cat {input.gtf} > {output.cat_gtf}
#         """

# --------------------- 2 MAPPING

rule bwa_index:
    """
        make index with bwa for all reference files
    """
    threads: get_threads('bwa_index', 1)
    input:
            reference = f"{out_dir}/1_cat_ref/all_ref.{suffix_file}",
    output:
            sa_file = f"{out_dir}/1_cat_ref/all_ref.{suffix_file}.sa"
    log:
            error = f"{log_dir}/bwa_Index/all_ref.e",
            output = f"{log_dir}/bwa_Index/all_ref.o"
    message:
            f"""
            {sep*108}
            Execute {{rule}} for 
                Input:
                    - Fasta : {{input.reference}}
                Output:
                    - sa file : {{output.sa_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            {sep*108}"""
    singularity:
            config["SINGULARITY"]["MAIN"]
    conda:
            "envs/bam_treatment.yaml"
    shell:
            f"""
                bwa index {{input.reference}} 1>{{log.output}} 2>{{log.error}}
            """

###########

## THIS USED TO BE HOOKED TO THE INITIAL FASTQs NOT THE FASTP ONES
## DELIBERATE?
## SHOULD BE ABLE TO SPECIFY OPT ARGS TO bwa aln
## {output.bam_file} SHOULD BE TEMP

rule run_mapping:
    """
        make bwa aln for all samples on all merged references + sorting + indexing
    """
    threads: get_threads('run_mapping', 1)
    input:
            reference = rules.bwa_index.input.reference,
            index = rules.bwa_index.output.sa_file,
            fastq = rules.run_fastp.output.fastq_trimmed,
            fastq_trimmed_info = rules.generate_fastqtrimmed_info.output.out_file
    output:
            sai_file = f"{out_dir}/2_mapping_sRNA/{{fastq}}.sai",
            sam_file = f"{out_dir}/2_mapping_sRNA/{{fastq}}.sam",
            bam_file = f"{out_dir}/2_mapping_sRNA/{{fastq}}.bam",
            sorted_bam_file = f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam",
            sorted_bam_index = f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam.bai",
    params:
            options= config['PARAMS']['BWA_ALN']['options']
    log:
            error = f"{log_dir}/bwa_aln/{{fastq}}.e",
            output = f"{log_dir}/bwa_aln/{{fastq}}.o"
    message:
            f"""
            {sep*108}
            Execute {{rule}} for 
                Input:
                    - Fasta : {{input.reference}}
                    - Fastq : {{input.fastq}}
                Output:
                    - sai file : {{output.sai_file}}
                    - sorted bam file : {{output.sorted_bam_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            {sep*108}"""
    singularity:
            config["SINGULARITY"]["MAIN"]
    conda:
            "envs/bam_treatment.yaml"
    shell:
        """
            # Align reads to reference genome:
            bwa aln -t {threads} {params.options} {input.reference} {input.fastq} > {output.sai_file} 2>>{log.error}
            bwa samse -n 5 {input.reference} {output.sai_file} {input.fastq} > {output.sam_file} 2>>{log.error}
            samtools view -bSh -o {output.bam_file} {output.sam_file} 1>>{log.output} 2>>{log.error}
            samtools sort --threads {threads} -o {output.sorted_bam_file} {output.bam_file} 1>>{log.output} 2>>{log.error}
            samtools index {output.sorted_bam_file} 1>>{log.output} 2>>{log.error}
        """

rule generate_bamfile_info:
    """
        generate_bamfile_info + bamlist for each treatment (to be use for samtools merge by treatment)
    """
    threads: get_threads('generate_bamfile_info', 1)
    input:
        bam_files = expand(rules.run_mapping.output.sorted_bam_file, fastq = SAMPLE_NAME),
        samplefile = samplefile
    params:
        outdir = f"{out_dir}/2_mapping_sRNA/"
    output:
        out_file = f"{out_dir}/2_mapping_sRNA/bamfile_info.txt"
    script:
        "script/write_bamfile_info.py"


rule samtools_stats:
    """
        make stats on mappings
    """
    threads: get_threads('samtools_stats', 1)
    input:
            sorted_bam_file = rules.run_mapping.output.sorted_bam_file,
            sorted_bam_index = rules.run_mapping.output.sorted_bam_index
    output:
            bamstats = f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam.bamStats.txt",
            idxstats = f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam.idxstats.log",
            flagstat = f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam.flagstat.log"
    log:
            error = f"{log_dir}/samtools_stats/{{fastq}}.e",
            output = f"{log_dir}/samtools_stats/{{fastq}}.o"
    message:
            f"""
            {sep*108}
            Execute {{rule}} for 
                Input:
                    - sorted_bam_file : {{input.sorted_bam_file}}
                Output:
                    - bam stats : {{output.bamstats}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            {sep*108}"""
    singularity:
            config["SINGULARITY"]["MAIN"]
    conda:
            "envs/bam_treatment.yaml"
    shell:
            """
                samtools stats {input.sorted_bam_file} > {input.sorted_bam_file}.bamStats.txt 2>>{log.error}
                samtools idxstats --threads {threads} {input.sorted_bam_file} > {input.sorted_bam_file}.idxstats.log 2>>{log.error}
                samtools flagstat --threads {threads} {input.sorted_bam_file} > {input.sorted_bam_file}.flagstat.log 2>>{log.error}
            """

# --------------------- 3 Merge

rule samtools_merge_treatment:
    """
        make samtools merge + index for each treatment
    """
    threads: get_threads('samtools_merge_treatment', 1)
    input:
            lien = expand(f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam",fastq = SAMPLE_NAME),
            bamfile = rules.generate_bamfile_info.output.out_file,
    params:
            bam_list = f"{out_dir}/2_mapping_sRNA/List_bamfiles_{{treatment}}",
    output:
            bam_file = f"{out_dir}/3_merge_bam_sRNA/{{treatment}}_merge.bam",
            bam_index = f"{out_dir}/3_merge_bam_sRNA/{{treatment}}_merge.bam.bai",
    log:
            error = f"{log_dir}/samtools_merge_treatment/{{treatment}}_merge.e",
            output = f"{log_dir}/samtools_merge_treatment/{{treatment}}_merge.o"
    message:
            f"""
            {sep*108}
            Execute {{rule}} for 
                Input:
                    - bam list : {{input.lien}}
                Output:
                    - bam file merged : {{output.bam_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            {sep*108}"""
    singularity:
            config["SINGULARITY"]["MAIN"]
    conda:
            "envs/bam_treatment.yaml"
    shell:
           f"""
                samtools merge -@ {{threads}} -b {{params.bam_list}} {{output.bam_file}} 1>>{{log.output}} 2>>{{log.error}}
                samtools index {{output.bam_file}} 1>>{{log.output}} 2>>{{log.error}}
           """

rule samtools_merge:
    """
        make samtools merge + index for all bam
    """
    threads: get_threads('samtools_merge', 1)
    input:
            all_bam = expand(rules.run_mapping.output.sorted_bam_file, fastq = SAMPLE_NAME)
    params:
            bam_dir_ls  = f"ls {out_dir}/2_mapping_sRNA/*.sorted.bam > {out_dir}/2_mapping_sRNA/bamlist", ## WHY NOT DO THAT DIRECTLY IN THE SHELL CODE OF THE RULE?
            bam_list = f"{out_dir}/2_mapping_sRNA/bamlist"
    output:
            bam_file = f"{out_dir}/3_merge_bam_sRNA/all_merge.bam",
            bam_index = f"{out_dir}/3_merge_bam_sRNA/all_merge.bam.bai",
    log:
            error = f"{log_dir}/samtools_merge/all_merge.e",
            output = f"{log_dir}/samtools_merge/all_merge.o"
    message:
            f"""
            {sep*108}
            Execute {{rule}} for 
                Input:
                    - bam list : {{input.all_bam}}
                Output:
                    - bam file merged : {{output.bam_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            {sep*108}"""
    singularity:
            config["SINGULARITY"]["MAIN"]
    conda:
            "envs/bam_treatment.yaml"
    shell:
           f"""
                {{params.bam_dir_ls}}
                samtools merge -@ {{threads}} -b {{params.bam_list}} {{output.bam_file}} 1>>{{log.output}} 2>>{{log.error}}
                samtools index {{output.bam_file}} 1>>{{log.output}} 2>>{{log.error}}
            """

# --------------------- 4 SHORTSTACK

rule ShortStack:
    """
        Launch ShortStack
    """
    threads: get_threads('ShortStack', 1)
    input:
        reference = rules.bwa_index.input.reference,
        bamfile = rules.samtools_merge.output.bam_file
    output:
        all_gff3 = f"{out_dir}/4_ShortStack/ShortStack_All.gff3",
        resultsFile = f"{out_dir}/4_ShortStack/Results.txt",
        sstackLog =  f"{out_dir}/4_ShortStack/Log.txt",
    params:
        outdir = directory(f"{out_dir}/4_ShortStack/"),
        dicermin = config["PARAMS"]["SHORTSTACK"]["dicermin"],
        dicermax = config["PARAMS"]["SHORTSTACK"]["dicermax"],
        mincov = config["PARAMS"]["SHORTSTACK"]["mincov"],
        pad = config["PARAMS"]["SHORTSTACK"]["pad"],
        foldsize = config["PARAMS"]["SHORTSTACK"]["foldsize"],
        strand_cutoff = config["PARAMS"]["SHORTSTACK"]["strand_cutoff"],
    log:
        error = f"{log_dir}/ShortStack/ShortStack.e",
        output = f"{log_dir}/ShortStack/ShortStack.o"
    singularity:
        config["SINGULARITY"]["MAIN"]
    conda:
        "envs/shortstack.yaml"
    shell:
        """
            if [ -d "{params.outdir}" ]; then rm -Rf {params.outdir}; fi
            ShortStack --genomefile {input.reference} --bamfile {input.bamfile} --outdir {params.outdir}  --dicermin {params.dicermin} --dicermax {params.dicermax} --mincov {params.mincov} --pad {params.pad} --foldsize {params.foldsize} --strand_cutoff {params.strand_cutoff} 1>>{log.output} 2>>{log.error}
        """

rule shortStack_populateGFF:
    """
        Annotate gff from ShortStack with overlapping reference loci
    """
    threads: get_threads('shortStack_populateGFF', 1)
    input:
        all_gff3 = rules.ShortStack.output.all_gff3,
        miRNA_gff = miRNA_path
    params:
        outDirSStack = rules.ShortStack.params.outdir
    output:
        new_gff3 = f"{out_dir}/4_ShortStack/ShortStack_All_miRNA.gff3"
    singularity:
        config["SINGULARITY"]["MAIN"]
    shell:
       """
       Rscript script/shortStack_populateGFF.R --sstackDir={params.outDirSStack} --miRnaGffFile={input.miRNA_gff} --outputGffFile={output.new_gff3}
       """

rule shortStack_analysis:
    """
        ShortStack Analysis:
        Generate a html summary of Shortstack results
    """
    threads: get_threads('shortStack_analysis', 1)
    input:
        new_gff3 = rules.shortStack_populateGFF.output.new_gff3,
        resultsFile = rules.ShortStack.output.resultsFile
    params:
        outDirSStack = rules.ShortStack.params.outdir,
        outDir = directory(f"{out_dir}/4_ShortStack/")
    output:
        html_output = f"{out_dir}/4_ShortStack/shortStack_analysis.html"
    singularity:
        config["SINGULARITY"]["MAIN"]
    script:
        "script/shortStack_analysis.Rmd"

rule diff_exp_analysis:
    """
       Experimental sRNA Clusters Differential Expression Analysis
    """
    threads: get_threads('diff_exp_analysis', 4)
    input:
        bam_files_info_file = rules.generate_bamfile_info.output.out_file,
        genome_sequence_file = rules.cat_fasta.output.cat_ref, ## THIS MAY NOT BE NECESSARY ANY MORE!!!!!
        genome_annotation_file = config["DATA"]["files"]["annotation"],
        sRNA_loci_annot_file = rules.shortStack_populateGFF.output.new_gff3
    params:
        out_dir = lambda w, output: os.path.dirname(output.html_output),
        de_comparisons_file = config["DATA"]["files"]["de_comparisons_file"],
        filter_gff = config["DATA"]["files"]["filter_gff"],
        minRowSumTreshold = config["PARAMS"]["DE_ANALYSIS"]["minRowSumTreshold"],
        variableOfInterest = config["PARAMS"]["DE_ANALYSIS"]["variableOfInterest"],
        batch = config["PARAMS"]["DE_ANALYSIS"]["batch"],
        locfunc = config["PARAMS"]["DE_ANALYSIS"]["locfunc"],
        fitType = config["PARAMS"]["DE_ANALYSIS"]["fitType"],
        sizeFactEstimMethod = config["PARAMS"]["DE_ANALYSIS"]["sizeFactEstimMethod"],
        pAdjustMethod = config["PARAMS"]["DE_ANALYSIS"]["pAdjustMethod"],
        cooksCutoff = config["PARAMS"]["DE_ANALYSIS"]["cooksCutoff"],
        independentFiltering = config["PARAMS"]["DE_ANALYSIS"]["independentFiltering"],
        lfcThreshold = config["PARAMS"]["DE_ANALYSIS"]["lfcThreshold"],
        altHypothesis = config["PARAMS"]["DE_ANALYSIS"]["altHypothesis"],
        alpha = config["PARAMS"]["DE_ANALYSIS"]["alpha"],
        typeTrans = config["PARAMS"]["DE_ANALYSIS"]["typeTrans"]
    output:
        html_output = f"{out_dir}/5_sRNA_loci_DE_analysis/sRNA_DE_analysis.html"
    log:
        error = f"{log_dir}/diff_exp_analysis.e",
        output = f"{log_dir}/diff_exp_analysis.o"
    singularity:
        config["SINGULARITY"]["MAIN"]
    script:
        "script/sRNA_DE_analysis.Rmd"


rule multiqc:
    """
        multiqc on outdir directory
    """
    threads: get_threads("multiqc", 1)
    input:
        expand(f"{out_dir}/0_QC/{{fastq}}_fastp.json", fastq = SAMPLE_NAME),
        expand(f"{out_dir}/0_QC/fastqc/{{fastq}}_raw_fastqc.html", fastq = SAMPLE_NAME),
        expand(f"{out_dir}/0_QC/fastqc/{{fastq}}_trimmed_fastqc.html", fastq = SAMPLE_NAME),
        expand(f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam.bamStats.txt", fastq = SAMPLE_NAME),
        expand(f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam.idxstats.log", fastq = SAMPLE_NAME),
        expand(f"{out_dir}/2_mapping_sRNA/{{fastq}}.sorted.bam.flagstat.log", fastq = SAMPLE_NAME),
    output:
        f"{out_dir}/6_MULTIQC/multiqc.html"
    log:
        f"{log_dir}/multiqc/multiqc.log"
    singularity:
        config["SINGULARITY"]["MAIN"]
    wrapper:
        "0.65.0/bio/multiqc"


# TODO :
        # scratch
        # enlever cat_gtf pour l'instant --> un gtf faire a la main
        # fastqc --> pas la wildcard fastq mais le nom du chemin en sortie

        # trim parameter in config for the user to decide if trimming should be executed?
        # Logs: I believe it is more desirable to have both stdout and stderr in the same file, could use the trick employed in baseDmux to simplify coding

