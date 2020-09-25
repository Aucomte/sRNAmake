#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @package countFromBam.py
# @author Sebastien Ravel
"""
	The countFromBam script
	=========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/11/2016
	:version: 0.1

	Script description
	------------------

	This Programme map fasta with bwa mem on reference file

	Example
	-------

	>>> countFromBam.py -f fasta/ -r referance.fasta -o output/

	Help Programm
	-------------

	optional arguments:
		- \-h, --help            show this help message and exit
		- \-v, --version         display ./countFromBam.py version number and exit
		- \-d, --debug          enter verbose/debug mode

	Input mandatory infos for running:
		- \-f <path/to/directory/fasta>, --fasta <path/to/directory/fasta>
						path to fasta files
		- \-r <path/to/dreference>, --ref <path/to/dreference>
						path to reference fasta files
		- \-o <path/to/directory>, --out <path/to/directory>
						Name of output file directory

	Input infos for running with default values:
		- \-j <int>, --nbjob <int>
						Number of job array lunch (default = 100)
		- \-th <int>, --thread <int>
						number of threads for mapping (default = 4)
"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
#from MODULES_SEB import *

## Python modules
import os
import sqlite3
import logging
from time import localtime

import argparse, pysam, traceback, sys, matplotlib
from pathlib import Path
from datetime import datetime
from collections import defaultdict, OrderedDict


import pandas as pd
# environment settings:
matplotlib.use('Agg')
pd.set_option('display.max_column',8)
pd.set_option('display.max_rows',50)
pd.set_option('display.max_seq_items',None)
pd.set_option('display.max_colwidth', 500)
pd.set_option('expand_frame_repr', True)
# pd.options.display.width = None
# remove pandas header styles
# this avoids the restriction that xlsxwriter cannot
# format cells where formatting was already applied
import pandas.io.formats.excel
pandas.io.formats.excel.header_style = None
from Bio.Seq import Seq
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
##################################################
## Variables Globales
version = "Release 0.1, 22nd of October, 2019"
shortVersion = 0.1
epilogTools = """******"""
descriptionTools = f"""
{"#"*80}
#
# More information:
#\tPipeline version: {version}
#
{"#"*80}
"""

##################################################
## Functions
def existant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open by default.

	:param x: a file path
	:type x: str()
	:rtype: PosixPath
	:return: PosixPath

	"""
	if not Path(x).exists():
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError(f"file '{x}' doesn't exist")
	elif not Path(x).is_file():
		raise argparse.ArgumentTypeError(f"{x} is not a valide file")

	return Path(x).resolve()

def getLibarySizeToDict(csv_file, str_size, sample_name):
	dico = pd.read_csv(csv_file, index_col=0, squeeze=True).to_dict(orient='index')
	return(dico[sample_name][str_size])

def normalizeData(value, library_size):
	return round((value/library_size)*1000000,2)

def toSTR(x):
	from ast import literal_eval
	if type(x) is float:
		print(f"FLOAT:{x}")
		return x
	if type(x) is str:
		print(f"str:{x}")
		return "/".join([f"{float(a):.0f}" for a in literal_eval(x)])
	if type(x) is list:
		print(f"List:{x}")
		return "/".join([f"{float(a):.0f}" for a in x])
	else:
		print(f"ELSE:{type(x)}")

def round_to_100_percent(number_set, digit_after_decimal=0):
	"""
		This function take a list of number and return a list of percentage, which represents the portion of each number in sum of all numbers
		Moreover, those percentages are adding up to 100%!!!
		Notice: the algorithm we are using here is 'Largest Remainder'
		The down-side is that the results won't be accurate, but they are never accurate anyway:)
	"""
	unround_numbers = [x / float(sum(number_set)) * 100 * 10 ** digit_after_decimal if float(sum(number_set)) != 0 else x/1 for x in number_set ]
	decimal_part_with_index = sorted([(index, unround_numbers[index] % 1) for index in range(len(unround_numbers))], key=lambda y: y[1], reverse=True)
	remainder = 100 * 10 ** digit_after_decimal - sum([int(x) for x in unround_numbers])
	index = 0
	while remainder > 0:
		unround_numbers[decimal_part_with_index[index][0]] += 1
		remainder -= 1
		index = (index + 1) % len(number_set)
	return [int(x) / float(10 ** digit_after_decimal) for x in unround_numbers]

###################################################
## Main code
##################################################
if __name__ == "__main__":
	# Parameters recovery
	parserOther = argparse.ArgumentParser(add_help=False)

	inOutOptional = parserOther.add_argument_group('Input infos not mandatory')
	inOutOptional.add_argument('-s', '--size', metavar="int", type = int, default=[20, 25],nargs = 2, required=False, dest = 'size_keep', help = 'size keep reads (default: 20 to 25)')
	# inOutOptional.add_argument('-p', '--nopdf',action='store_false', dest = 'build_PDF', help='do not build pdf of count if add arg')
	inOutOptional.add_argument('-t', '--notable',action='store_false', dest = 'build_table', help='do not build table if add arg')
	inOutOptional.add_argument('-v', '--version', action='version', version=version, help=f'Use if you want to know which version of {__file__} you are using')
	inOutOptional.add_argument('-h', '--help', action='help', help=f'show this help message and exit')
	inOutOptional.add_argument('-d', '--debug',action='store_true', help='enter verbose/debug mode')

	parserMandatory = argparse.ArgumentParser(
									 parents=[parserOther],
									 add_help=False,
									 prog=__file__,
									 formatter_class=argparse.RawDescriptionHelpFormatter,
									 description=descriptionTools,
									 epilog=epilogTools
									)

	inOutMandatory = parserMandatory.add_argument_group('Input mandatory infos for running')
	inOutMandatory.add_argument('-bf', '--bam', metavar="<path/to/file/bam>", type=existant_file, required=True, dest = 'bam_file_filter', help = 'Path to the bam filter file')
	inOutMandatory.add_argument('-ba', '--bamall', metavar="<path/to/file/bam>", type=existant_file, required=True, dest = 'bam_file', help = 'Path to the bam file not filter')
	inOutMandatory.add_argument('-l', '--lib', metavar="<path/to/file/csv>", type = existant_file, required=True, dest = 'csv_file', help = 'path to file with size library, use to normalize data')
	# inOutMandatory.add_argument('-o', '--outputDir', metavar="DIR", type = Path, required=True, dest = 'outputDir', help = 'Output folder name (it will be created)')

	# Check parameters
	args = parserMandatory.parse_args()

	#Welcome message
	print(f"""{"#"*80}\n#{Path(__file__).stem+" "+version:^78}#\n{"#"*80}\nStart time: {datetime.now():%d-%m-%Y at %H:%M:%S}\nCommande line run: {" ".join(sys.argv)}\n""")
	# resume to user
	print(" - Intput Info:")
	for k,v in vars(args).items():
		print(f"\t - {k}: {v}")


	# Récupère les arguments et initialisation des variables
	output_dir = Path(args.bam_file_filter).resolve().parent
	sample_name, reference_name ,*_ = Path(args.bam_file_filter).stem.split("_")

	size_keep_min, size_keep_max = args.size_keep
	list_size_keep = range(size_keep_min, size_keep_max+1)
	str_size = f"{size_keep_min}-{size_keep_max}"

	library_size = getLibarySizeToDict(args.csv_file, str_size, sample_name)

	print(library_size)
	print(sample_name)
	print(reference_name)

	listReadName = {}
	dicoResume = defaultdict(OrderedDict)
	dicoCountPos = defaultdict(OrderedDict)
	dicoCountFirstNT = defaultdict(OrderedDict)

	# index les bam si besoin
	if not os.path.exists(args.bam_file_filter.as_posix()+".bai"): pysam.index(args.bam_file_filter.as_posix())
	if not os.path.exists(args.bam_file.as_posix()+".bai"): pysam.index(args.bam_file.as_posix())


	## utilisation de pysam view pour compter les reads forward and reverse total pour 1-44
	nbReadsForTotal = int(pysam.view("-c","-F","16",args.bam_file.as_posix()))
	nbReadsRefTotal = int(pysam.view("-c","-f","16",args.bam_file.as_posix()))
	nbReadsTotal = nbReadsForTotal+nbReadsRefTotal
	dicoResume["size 1-44 forward/reverse/total"].update({f"{sample_name} {reference_name} Total":[nbReadsForTotal,nbReadsRefTotal,nbReadsTotal]})

	## utilisation de pysam view pour compter les reads forward and reverse total en perfect match pour 1-44
	samNbReadsForTotal_pm = pysam.view("-F","16",args.bam_file.as_posix())
	samNbReadsRefTotal_pm = pysam.view("-f","16",args.bam_file.as_posix())
	nbReadsForTotal_pm = len([read for read in samNbReadsForTotal_pm.split("\n") if "NM:i:0" in read])
	nbReadsRefTotal_pm = len([read for read in samNbReadsRefTotal_pm.split("\n") if "NM:i:0" in read])
	nbReadsTotal_pm = nbReadsForTotal_pm+nbReadsRefTotal_pm
	dicoResume["size 1-44 perfect match forward/reverse/total"].update({f"{sample_name} {reference_name} Total":[nbReadsForTotal_pm,nbReadsRefTotal_pm,nbReadsTotal_pm]})

	## utilisation de pysam view pour compter les reads forward and reverse sur 20-25
	nbReadsForSize = int(pysam.view("-c","-F","16",args.bam_file_filter.as_posix()))
	nbReadsRefSize = int(pysam.view("-c","-f","16",args.bam_file_filter.as_posix()))
	nbReadsSize = nbReadsForSize+nbReadsRefSize
	dicoResume[f"size {str_size} forward/reverse/total"].update({f"{sample_name} {reference_name} Total":[nbReadsForSize,nbReadsRefSize,nbReadsSize]})
	dicoResume[f"size {str_size} forward/reverse/total normalized"].update({f"{sample_name} {reference_name} Total":[normalizeData(nbReadsForSize,library_size),normalizeData(nbReadsRefSize,library_size),normalizeData(nbReadsSize,library_size)]})

	## utilisation de pysam view pour compter les reads forward and reverse sur 20-25 en perfect match
	samForSize = pysam.view("-F","16",args.bam_file_filter.as_posix())
	samRefSize = pysam.view("-f","16",args.bam_file_filter.as_posix())
	nbReadsForSize_pm = len([read for read in samForSize.split("\n") if "NM:i:0" in read])
	nbReadsRefSize_pm = len([read for read in samRefSize.split("\n") if "NM:i:0" in read])
	nbReadsSize_pm = nbReadsForSize_pm+nbReadsRefSize_pm

	dicoResume[f"size {str_size} perfect match forward/reverse/total"].update({f"{sample_name} {reference_name} Total": [nbReadsForSize_pm,nbReadsRefSize_pm,nbReadsSize_pm]})

	### Initialisation du dictionnaire pour résumer les infos
	for i in list_size_keep:
		dicoResume[f"size {str_size} perfect match forward/reverse/total"].update({f"{sample_name} {reference_name} {i}nt": [0,0,0]})
		dicoCountFirstNT[(f"Forward {i}nt")].update({"A":0,"T":0,"G":0,"C":0})
		dicoCountFirstNT[(f"Reverse {i}nt")].update({"A":0,"T":0,"G":0,"C":0})
		dicoCountFirstNT[(f"Forward {i}nt")].update({"Total":0})
		dicoCountFirstNT[(f"Reverse {i}nt")].update({"Total":0})

	# dataframe = pd.DataFrame.from_dict(dicoResume)
	# print(dataframe)

	#### Parcours du BAM avec pysam
	bamfile = pysam.AlignmentFile(args.bam_file_filter.as_posix(), "rb")

	dicoRefSize = dict([ (dico["SN"], dico["LN"]) for dico in bamfile.header["SQ"] ])					# lit le header et retourne une liste des séquences avec leurs tailles
	print(dicoRefSize)

	# exit()
	for ref in bamfile.references: 																		# parcours des séquences (chromosome, BAC, ...)
		print(ref)
		refSize = dicoRefSize[ref]																		# récupère la taille des séquences de référence
		try:																							# rentre en mode test du code
			for pos in range(1,refSize+1):																# parcours de la référence de la position 1 à fin

				iter = bamfile.fetch(ref,pos-1,pos, until_eof=True)										# génère un itérateur pour parcourir la liste des reads à la position donnée (pos-1, pos)

				# pour les sequence avec 50pbajoutees
				if "50pbajoutees" in ref and pos <=50 :
					oldPos= pos
					pos = refSize-(100-pos)
					# print(f"Inf50\tOLD POS:{oldPos}\t NEW POS:{pos}")
				elif "50pbajoutees" in ref and pos > 50:
					oldPos= pos
					pos = pos - 50
					# print(f"Sup50\tOLD POS:{oldPos}\t NEW POS:{pos}")
					# exit()

				# initialise le dico de count / pos à zero
				if args.build_table:
					for i in list_size_keep:
						dicoCountPos[f"Forward {i}nt"][pos] = 0
						dicoCountPos[f"Reverse {i}nt"][pos] = 0
						dicoCountPos["Total Forward"][pos] = 0
						dicoCountPos["Total Reverse"][pos] = 0
						dicoCountPos["Total Reverse-Forward"][pos] = 0

				# # affichage de l'avancement du processus#
				if ((pos % 10000 == 0) and (pos != 0)) or (float(pos) == refSize):
					percent = (float(pos)/float(refSize))*100
					print(f"Processed up to {percent:0.2f} at position {pos}", end='\n',flush=True)


				for read in iter:																		# Parcours les reads le long de la ref pour une position donnée
					readGroupFile = int(read.get_tags()[-1][-1].split("-")[-1].replace(f"size",""))		# regarde d'ou provient le read (file size dans le readGroup) donc taille du read
					readSequence = read.query_sequence												# récupère la séquence du read
					# print(f"readGroupFile: {readGroupFile}")
					# exit()

					dicoTags = dict(read.get_tags())
					cigarValue = str(read.cigarstring.replace("M",""))

					# test si en Perfect Match
					if cigarValue.isdigit() and dicoTags["MD"].isdigit() and int(cigarValue) == int(dicoTags["MD"]) :	# si en perfect match

						#regarde si le read à déjà était lu
						if read.query_name not in listReadName:
							listReadName[read.query_name] = ""									# ajoute a la liste des reads trouvés
							if read.is_reverse:														# test si en reverse
								firstNT = Seq(readSequence).reverse_complement()[0]					# créer un object Seq pour faire le reverse_complement et capturer le 1er nucléotide (cf Bio::SeqIO de biopython)
								if args.build_table:
									dicoCountPos[f"Reverse {readGroupFile}nt"][pos]+=1
									dicoCountPos["Total Reverse"][pos]+=1
									dicoCountPos["Total Reverse-Forward"][pos]+=1
								dicoResume[f"size {str_size} perfect match forward/reverse/total"][f"{sample_name} {reference_name} {readGroupFile}nt"][1]+=1
								dicoResume[f"size {str_size} perfect match forward/reverse/total"][f"{sample_name} {reference_name} {readGroupFile}nt"][2]+=1
								if firstNT != "N":
									dicoCountFirstNT[(f"Reverse {readGroupFile}nt")][firstNT] = (dicoCountFirstNT[(f"Reverse {readGroupFile}nt")].get(firstNT, 1))+1
							else:																	# si non en forward
								firstNT = readSequence[0]
								if args.build_table:
									dicoCountPos[f"Forward {readGroupFile}nt"][pos]+=1
									dicoCountPos["Total Forward"][pos]+=1
									dicoCountPos["Total Reverse-Forward"][pos]+=1
								dicoResume[f"size {str_size} perfect match forward/reverse/total"][f"{sample_name} {reference_name} {readGroupFile}nt"][0]+=1
								dicoResume[f"size {str_size} perfect match forward/reverse/total"][f"{sample_name} {reference_name} {readGroupFile}nt"][2]+=1
								if firstNT != "N":
									dicoCountFirstNT[(f"Forward {readGroupFile}nt")][firstNT] = (dicoCountFirstNT[(f"Forward {readGroupFile}nt")].get(firstNT, 1))+1

		except Exception as e:
			print("\n############ ERROR ############\n")
			print(traceback.format_exc())
			print("\n\n")
			print(read)
			print("query_name\t{}".format(read.query_name))
			print("get_reference_positions\t{}".format(read.get_reference_positions()))
			print("cigar\t{}".format(read.cigar))
			print("flag\t{}".format(read.flag))
			print("reference_id\t{}".format(read.reference_id))
			print("reference_start\t{}".format(read.reference_start))
			print("mapping_quality\t{}".format(read.mapping_quality))
			print("query_qualities\t{}".format(read.query_qualities))
			print("tags\t{}".format(read.tags))
			print("If is Reverse\t{}".format(read.is_reverse))
			print("cigarstring\t{}".format(read.cigarstring))
			print("Cigar TAG\t{}".format(read.get_tags()))
			print("readGroupFile\t{}".format(readGroupFile))
			exit()

	print("\n\nEnd bam read\n\n")
	######## Utilisation de Pandas pour générer les tableaux CSV à partir des dico python
	baseNameOutputCSV = f"{output_dir}/{sample_name}_{reference_name}"  # base de tous le nom de tableau


	# tableau des count des premiers nucléotides
	tableCountFirstNT = pd.DataFrame(dicoCountFirstNT)

	# génére la somme dans le tableau
	for i in list_size_keep:
		tableCountFirstNT[(f"Forward {i}nt")]["Total"] = tableCountFirstNT[f"Forward {i}nt"].sum(axis=0)
		tableCountFirstNT[(f"Reverse {i}nt")]["Total"] = tableCountFirstNT[f"Reverse {i}nt"].sum(axis=0)

		#initialise le dicoResume pour le passage en pourcentage
		dicoResume["5'nt %A/%T/%G/%C forward"][f"{sample_name} {reference_name} {i}nt"] = []
		dicoResume["5'nt %A/%T/%G/%C reverse"][f"{sample_name} {reference_name} {i}nt"] = []
		dicoResume["5'nt %A/%T/%G/%C total"][f"{sample_name} {reference_name} {i}nt"] = []

	# calcul pour chaque NT / taille a garder
	for letter in ["A","T","G","C"]:
		for i in list_size_keep:
			totalSizeFor = tableCountFirstNT[(f"Forward {i}nt")]["Total"]
			totalSizeRev = tableCountFirstNT[(f"Reverse {i}nt")]["Total"]
			totalSizeForRev = totalSizeFor+totalSizeRev
			if totalSizeForRev == 0:
				poucentageAFor = 0
				poucentageARev = 0
				poucentageAForRev = 0
			else:
				poucentageAFor = np.rint((tableCountFirstNT[(f"Forward {i}nt")][letter]/totalSizeForRev)*100)
				poucentageARev = np.rint((tableCountFirstNT[(f"Reverse {i}nt")][letter]/totalSizeForRev)*100)
				poucentageAForRev = np.rint(poucentageAFor+poucentageARev)

			if np.isnan(poucentageAFor):
				poucentageAFor = 0
			if np.isnan(poucentageARev):
				poucentageAFor = 0
			if np.isnan(poucentageAForRev):
				poucentageAFor = 0

			dicoResume["5'nt %A/%T/%G/%C forward"][f"{sample_name} {reference_name} {i}nt"].append(poucentageAFor)
			dicoResume["5'nt %A/%T/%G/%C reverse"][f"{sample_name} {reference_name} {i}nt"].append(poucentageARev)
			dicoResume["5'nt %A/%T/%G/%C total"][f"{sample_name} {reference_name} {i}nt"].append(poucentageAForRev)

	# appel la fonction pour arrondir a 100%
	for i in list_size_keep:
		dicoResume["5'nt %A/%T/%G/%C total"][f"{sample_name} {reference_name} {i}nt"] = round_to_100_percent(dicoResume["5'nt %A/%T/%G/%C total"][f"{sample_name} {reference_name} {i}nt"])

	dataframeResume = pd.DataFrame.from_dict(dicoResume)
	# print(dataframeResume[f"size {str_size} perfect match forward/reverse/total"])
	# on normalize les données
	dataframeResume[f"size {str_size} perfect match forward/reverse/total normalized"] = dataframeResume[f"size {str_size} perfect match forward/reverse/total"].map(lambda a: str( [ normalizeData(value,library_size) for value in a ]) )
	# On convertit le tout en string avec separateur /
	dataframeResume = dataframeResume.applymap(lambda x: toSTR(x))

	# sauvegarde dans un fichier CSV tabulé
	with open(f"{baseNameOutputCSV}_Resume.csv", "w") as libsizeFile:
		print(f"File Resume:\n{dataframeResume}\n")
		dataframeResume.to_csv(libsizeFile, index=True, sep="\t", float_format=".", header=True)


	#### Si tableau de count par position
	if args.build_table:
		dfCountPos = pd.DataFrame(dicoCountPos).sort_index()

		# on passe les valeurs de reverse en negative
		for i in list_size_keep:
			dfCountPos[(f"Reverse {i}nt")] = dfCountPos[(f"Reverse {i}nt")]*-1
		dfCountPos["Total Reverse"] = dfCountPos["Total Reverse"]*-1

		print(f"File table_pm :\n{dfCountPos}\n")
		with open(f"{baseNameOutputCSV}_table_pm_{str_size}.csv", "w") as countPosFile:
			dfCountPos.to_csv(countPosFile, index=True, index_label="position", sep="\t", float_format=".", header=True, columns=sorted(dfCountPos.columns))

		with open(f"{baseNameOutputCSV}_table_pm_{str_size}_normalised.csv", "w") as countPosFile:
			dfCountPosNormalized = dfCountPos.copy(deep=False).applymap(lambda value: str(round((value/library_size)*1000000,2)))
			dfCountPosNormalized.to_csv(countPosFile, index=True, index_label="position", sep="\t", float_format=".", header=True, columns=sorted(dfCountPos.columns))

	# if args.build_PDF:
		# with PdfPages(f"{baseNameOutputCSV}_countReads_{str_size}.pdf") as pdf:
			# for i in list_size_keep:
				# selectFor = dfCountPosNormalized[f"Forward {i}nt"].to_dict()
				# selectRev = dfCountPosNormalized[f"Reverse {i}nt"].to_dict()
				# maxFor = max(selectFor.values())
				# minRev = min(selectRev.values())
				# fig = plt.figure(figsize=(50, 8),dpi=800,facecolor='w', edgecolor='w') #
				# ax = plt.subplot(111)

				# ind = np.arange(len(selectRev))    # the x locations for the groups
				# width = 0.35       # the width of the bars: can also be len(x) sequence

				# p1 = plt.bar(ind, selectFor.values(), width, color='r')
				# p2 = plt.bar(ind, selectRev.values(), width, color="g")

				# plt.ylabel('Count Reads')
				# plt.title(f"{i} nt")
				# plt.xticks(np.arange(1, len(selectRev), 100))
				# #plt.yticks(np.arange(minRev-1,maxFor+1, 1))
				# plt.legend((p1[0], p2[0]), ('Forward', 'Reverse'))
				# pdf.savefig(fig)
				# plt.close()


	#print("\nStop time: ", datetime.now().strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")

	#exit()
