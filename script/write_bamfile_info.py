import time
import csv
import re
import pandas as pd

bam = snakemake.input.bam_files
samplefile = snakemake.input.samplefile
bamfile = snakemake.output.out_file
outdir = snakemake.params.outdir


################################################################
def unique(list1):
    # intilize a null list
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)

################################################################

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

#on récupère la liste des traitements
treatments = []
for key in samples:
    if key[0] == "treatment":
        treatments.append(samples[key])
TREATMENT = unique(list(treatments))

for x in samples:
    if x[0] == "filename":
        for i in bam:
            if re.search(samples["samplename", x[1]],i):
                samples[x] = i

with open(bamfile, 'w') as f:
    writer = csv.writer(f, delimiter = "\t")
    writer.writerow(header)
    for i in range(len(header)):
        for key in samples.keys():
            if key[1] == i:
                f.write(f"{samples[key]}\t")
        f.write(f"\n")

####################################
## ecrire les listes de bam par traitement

for treatment in TREATMENT:
    samples = {}
    header = ""
    with open(bamfile, 'r') as f:
        spamreader = csv.reader(f, delimiter = "\t")
        for line in spamreader:
            header = line
            break
        for i in range(len(header)-1):
            with open(bamfile, 'r') as f:
                spamreader = csv.reader(f, delimiter = "\t")
                x = 0
                for line in spamreader:
                    if line != header :
                        samples[header[i].lower(),x] = [line[i]]
                        x = x + 1
        l = list()
        for key in samples:
            if key[0].lower() == "treatment":
                if samples[key][0] == treatment:
                    l.append(samples["filename", key[1]][0])
        treatmentfile = f"{outdir}List_bamfiles_{treatment}"
        with open(treatmentfile, 'w') as f:
            for key in l:
                f.write(f"{key}\n")

time.sleep(10)
