#!/home/mu42cuq/programs/Python-3.6.2/bin/python3.6
# -*- coding: utf-8 -*- 

import os
import sys
import argparse
import pysam
import multiprocessing
from mypymo import ls, sep3 #TODO: replace this dependency
import numpy as np #import array, zeros
from scipy.stats import ttest_ind_from_stats, gmean, norm, ttest_ind, gamma
from scipy.integrate import quad
from statsmodels.sandbox.stats.multicomp import multipletests
import pickle
import seaborn as sns
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import ColumnDataSource, CDSView, IndexFilter, HoverTool, BooleanFilter
from bokeh.layouts import gridplot




#####General Argument Parser
if(len(sys.argv) == 1 or (sys.argv[1] != 'count' and sys.argv[1] != 'analyze')):
	parser = argparse.ArgumentParser(description='MeRDE: MicroRNA Differential Expression Analysis.', epilog='Version 0.9 (Mai 2018)')
	parser.add_argument('Mode', help='Start MeRDE in COUNT mode to generate read count files of a given annotation and .bam mapping file(s) OR start MeRDE in ANALYZE mode to calculate differentially expressed miRNAs from given read count files.', choices=['count', 'analyze'])
	args = parser.parse_args()
	
#############


#####COUNT MODE
elif(sys.argv[1] == 'count'):
	parser = argparse.ArgumentParser(description='MeRDE: MicroRNA Differential Expression Analysis -- COUNT Mode.', epilog='Version 1.2 (September 2017)')
	parser.add_argument('Mode', help='Start MeRDE in COUNT mode to generate read count files of a given annotation and .bam mapping file(s) OR start MeRDE in ANALYZE mode to calculate differentially expressed miRNAs from given read count files.', choices=['count', 'analyze'])
	parser.add_argument('-a', '-annotation', help='Input annotation file in gff/gtf format.', metavar='<File>')
	parser.add_argument('-i', '-input_files', help='Input mapping files in bam format.', metavar='<File>', nargs='*')
	parser.add_argument('-s','-stranded', help='Input mapping data is strand specific.', action='store_true')
	parser.add_argument('-O','-readOverlap', help='Minimal amount of overlapping bases of a read and a feature to consider it for counting. [Default: 1]\n If a float number between 0.0 and 1.0 is given, it will be interpreted as the required overlapping fraction of read bases.', metavar='<Integer/Float>', default=1)
	parser.add_argument('-f','-feature_type', help='Type of feature in the 3rd column of the annotation file to consider as miRNA gene entries. [Default: miRNA]', default='miRNA')
	parser.add_argument('-p','-pre_miRNA', help='The given annotations refers to pre-mature miRNA genes, therefore the expression of their 5\' and 3\' ends will be counted individually.', action='store_true')
	parser.add_argument('-u','-unique_counting', help='Count only uniquely mapped reads.', action='store_true')
	parser.add_argument('-cpu' , help='Amount of CPUs to use for counting. [Default: 1]', metavar='<Integer>', default=1)
	parser.add_argument('-o', '-output_dir', help='Directory where output will be saved. [Default: current working directory]', metavar='<Path>', default=os.getcwd() +'/')
	parser.add_argument('-q', '-quiet', help='Do not write status information into console.', action='store_true')
	args = parser.parse_args()


	####FUNCTIONS

	def find_gene_name(lineSplit):
		if(lineSplit[8].find('gene_id=') != -1):
			return(lineSplit[8].split('gene_id=')[1].split(';')[0].strip())
		elif(lineSplit[8].find('gene_id "') != -1):
			return(lineSplit[8].split('gene_id "')[1].split('"')[0].strip())
		elif(lineSplit[8].find('gene_id ') != -1):
			return(lineSplit[8].split('gene_id ')[1].split(';')[0].strip())
		elif(lineSplit[8].find('geneid=') != -1):
			return(lineSplit[8].split('geneid=')[1].split(';')[0].strip())
		elif(lineSplit[8].find('geneid "') != -1):
			return(lineSplit[8].split('geneid "')[1].split('"')[0].strip())
		elif(lineSplit[8].find('geneid ') != -1):
			return(lineSplit[8].split('geneid ')[1].split(';')[0].strip())
		elif(lineSplit[8].find('ID=') != -1):
			return(lineSplit[8].split('ID=')[1].split(';')[0].strip())
		elif(lineSplit[8].find('ID "') != -1):
			return(lineSplit[8].split('ID "')[1].split('"')[0].strip())
		elif(lineSplit[8].find('ID ') != -1):
			return(lineSplit[8].split('ID ')[1].split(';')[0].strip())
		elif(lineSplit[8].find('gene_name=') != -1):
			return(lineSplit[8].split('gene_name=')[1].split(';')[0].strip())
		elif(lineSplit[8].find('gene_name "') != -1):
			return(lineSplit[8].split('gene_name "')[1].split('"')[0].strip())
		elif(lineSplit[8].find('gene_name ') != -1):
			return(lineSplit[8].split('gene_name ')[1].split(';')[0].strip())
		elif(lineSplit[8].find(';gene=') != -1):
			return(lineSplit[8].split(';gene=')[1].split(';')[0].strip())
		elif(lineSplit[8].find(';gene "') != -1):
			return(lineSplit[8].split(';gene "')[1].split('"')[0].strip())
		elif(lineSplit[8].find(';gene ') != -1):
			return(lineSplit[8].split(';gene ')[1].split(';')[0].strip())
		else:
			return('noNameFound_' +lineSplit[0] +'_'.join([lineSplit[0],lineSplit[3],lineSplit[4],lineSplit[6],]))

	def count_Reads(currentRead, currentExon):
		if(args.u):
			if(currentRead.is_secondary):               #check if secondary alignment, if unique mapping option is enabled
				return(False)
		if(currentRead.reference_length > 50):		#check if it is a split read
			return(False)
		if(args.s):						#TODO: allow reverse strand-specific counting
			if(currentRead.is_reverse):                  #check if strands match, if strand option is enabled
				if((currentRead.positions[-1]-int(currentExon[0])+1 < minOverlap) or (int(currentExon[1])-currentRead.positions[0]+1 < minOverlap)): #check read - feature overlap size)
					return(False)
				else:
					return(-1)
		
		if(type(minOverlap) == float):		#fractional read lenght overlap
			if((currentRead.positions[-1]-int(currentExon[0])+1 < (len(currentRead.positions)*minOverlap)) or (int(currentExon[1])-currentRead.positions[0]+1 < (len(currentRead.positions)*minOverlap))): #check read - feature overlap size)
				return(False)
			else:
				return(1)
		else:								#absolute bases overlap
			if((currentRead.positions[-1]-int(currentExon[0])+1 < minOverlap) or (int(currentExon[1])-currentRead.positions[0]+1 < minOverlap)): #check read - feature overlap size)
				return(False)
			else:
				return(1)
	
	#outdated function
	def perform_counting(bfile, outfile, annoLines):
		bamFile = pysam.AlignmentFile(bfile, 'rb')
		with open(outfile +'.count.tmp', 'w') as currentOutfile:
			for line in annoLines:
				lineSplit = line.split('\t')
				currentName = find_gene_name(lineSplit)
				if(args.p):					#given Annotations are pre-mature miRNAs
					currentCountPlusStrand5Prime = 0
					currentCountMinusStrand5Prime = 0
					currentCountPlusStrand3Prime = 0
					currentCountMinusStrand3Prime = 0
										
					try:					#split gene annotation and corresponding mapped reads to count for 5' and 3' part individually
						if(lineSplit[6] == '+'):
							fivePrimeReads = bamFile.fetch(lineSplit[0], int(lineSplit[3]), int(lineSplit[3]) + ((int(lineSplit[4]) - int(lineSplit[3]))//2))
							threePrimeReads = bamFile.fetch(lineSplit[0], (int(lineSplit[3]) + (int(lineSplit[4]) - int(lineSplit[3]))//2) + 1, int(lineSplit[4])) 
						else:
							threePrimeReads = bamFile.fetch(lineSplit[0], int(lineSplit[3]), int(lineSplit[3]) + ((int(lineSplit[4]) - int(lineSplit[3]))//2))
							fivePrimeReads = bamFile.fetch(lineSplit[0], (int(lineSplit[3]) + (int(lineSplit[4]) - int(lineSplit[3]))//2) + 1, int(lineSplit[4]))
					

						for read in fivePrimeReads:														#actual counting
							currentRes = count_Reads(read, (lineSplit[3], lineSplit[4], lineSplit[6], lineSplit[0]))
							if(currentRes):
								if(currentRes == 1):
									currentCountPlusStrand5Prime += 1
								else:
									currentCountMinusStrand5Prime += 1

						for read in threePrimeReads:														#actual counting
							currentRes = count_Reads(read, (lineSplit[3], lineSplit[4], lineSplit[6], lineSplit[0]))
							if(currentRes):
								if(currentRes == 1):
									currentCountPlusStrand3Prime += 1
								else:
									currentCountMinusStrand3Prime += 1

						if(args.s):																	#write count number
							currentOutfile.write(currentName +'_5Prime_Plus\t' +str(currentCountPlusStrand5Prime) +'\n')
							currentOutfile.write(currentName +'_5Prime_Minus\t' +str(currentCountMinusStrand5Prime) +'\n')
							currentOutfile.write(currentName +'_3Prime_Plus\t' +str(currentCountPlusStrand3Prime) +'\n')
							currentOutfile.write(currentName +'_3Prime_Minus\t' +str(currentCountMinusStrand3Prime) +'\n')
						else:
							currentOutfile.write(currentName +'_5Prime\t' +str(currentCountPlusStrand5Prime) +'\n')		#if counting is not strand specific, every read will be counted for plus strand (see count_Reads function)
							currentOutfile.write(currentName +'_3Prime\t' +str(currentCountPlusStrand3Prime) +'\n')
					except ValueError:
						if not(args.q):
							sys.stderr.write('Fetching reads went wrong for the following annotation line (so it was skipped): ' +line +'\n\n')

				else:
					currentCountPlusStrand = 0
					currentCountMinusStrand = 0
					try:
						currentReads = bamFile.fetch(lineSplit[0], int(lineSplit[3]), int(lineSplit[4]))                					#fetch mapped reads               

						for read in currentReads:															#actual counting
							currentRes = count_Reads(read, (lineSplit[3], lineSplit[4], lineSplit[6], lineSplit[0]))
							if(currentRes):
								if(currentRes == 1):
									currentCountPlusStrand += 1
								else:
									currentCountMinusStrand += 1
						if(args.s):																	#write count number
							currentOutfile.write(currentName +'_Plus\t' +str(currentCountPlusStrand) +'\n')
							currentOutfile.write(currentName +'_Minus\t' +str(currentCountMinusStrand) +'\n')
						else:
							currentOutfile.write(currentName +'\t' +str(currentCountPlusStrand) +'\n')				#if counting is not strand specific, every read will be counted for plus strand (see count_Reads function)
					except ValueError:
						if not(args.q):
							sys.stderr.write('Fetching reads went wrong for the following annotation line (so it was skipped): ' +line +'\n\n')

	def perform_counting2(bfile, currentCountDict, annoLines):
		bamFile = pysam.AlignmentFile(bfile, 'rb')
		for line in annoLines:
			lineSplit = line.split('\t')
			currentName = find_gene_name(lineSplit)
			if(args.p):					#given Annotations are pre-mature miRNAs
				currentCountPlusStrand5Prime, currentCountMinusStrand5Prime = 0, 0
				currentCountPlusStrand3Prime, currentCountMinusStrand3Prime = 0, 0
													
				try:					#split gene annotation and corresponding mapped reads to count for 5' and 3' part individually
					if(lineSplit[6] == '+'):
						fivePrimeReads = bamFile.fetch(lineSplit[0], int(lineSplit[3]), int(lineSplit[3]) + ((int(lineSplit[4]) - int(lineSplit[3]))//2))
						threePrimeReads = bamFile.fetch(lineSplit[0], (int(lineSplit[3]) + (int(lineSplit[4]) - int(lineSplit[3]))//2) + 1, int(lineSplit[4])) 
					else:
						threePrimeReads = bamFile.fetch(lineSplit[0], int(lineSplit[3]), int(lineSplit[3]) + ((int(lineSplit[4]) - int(lineSplit[3]))//2))
						fivePrimeReads = bamFile.fetch(lineSplit[0], (int(lineSplit[3]) + (int(lineSplit[4]) - int(lineSplit[3]))//2) + 1, int(lineSplit[4]))
				

					for read in fivePrimeReads:														#actual counting
						currentRes = count_Reads(read, (lineSplit[3], lineSplit[4], lineSplit[6], lineSplit[0]))
						if(currentRes == 1):
							currentCountPlusStrand5Prime += 1
						elif(currentRes == -1):
							currentCountMinusStrand5Prime += 1

					for read in threePrimeReads:														#actual counting
						currentRes = count_Reads(read, (lineSplit[3], lineSplit[4], lineSplit[6], lineSplit[0]))
						if(currentRes == 1):
							currentCountPlusStrand3Prime += 1
						elif(currentRes == -1):
							currentCountMinusStrand3Prime += 1

					if(args.s):
						with lock:															#write count number
							currentCountDict[currentName +'_5P+'] = currentCountPlusStrand5Prime
							currentCountDict[currentName +'_5P-'] = currentCountMinusStrand5Prime
							currentCountDict[currentName +'_3P+'] = currentCountPlusStrand3Prime
							currentCountDict[currentName +'_3P-'] = currentCountMinusStrand3Prime
					else:
						with lock:
							currentCountDict[currentName +'_5P'] = currentCountPlusStrand5Prime  	#if counting is not strand specific, every read will be counted for plus strand (see count_Reads function)
							currentCountDict[currentName +'_3P'] = currentCountPlusStrand3Prime
				except ValueError:
					if not(args.q):
						sys.stderr.write('Fetching reads went wrong for the following annotation line (so it was skipped): ' +line +'\n\n')

			else:
				currentCountPlusStrand, currentCountMinusStrand = 0, 0
				try:
					currentReads = bamFile.fetch(lineSplit[0], int(lineSplit[3]), int(lineSplit[4]))                					#fetch mapped reads               

					for read in currentReads:															#actual counting
						currentRes = count_Reads(read, (lineSplit[3], lineSplit[4], lineSplit[6], lineSplit[0]))
						if(currentRes == 1):
							currentCountPlusStrand += 1
						elif(currentRes == -1):
							currentCountMinusStrand += 1
					if(args.s):
						with lock:																	#write count number
							currentCountDict[currentName +'_+'] = currentCountPlusStrand
							currentCountDict[currentName +'_-'] = currentCountMinusStrand
					else:
						with lock:
							currentCountDict[currentName] = currentCountPlusStrand			#if counting is not strand specific, every read will be counted for plus strand (see count_Reads function)
				except ValueError:
					if not(args.q):
						sys.stderr.write('Fetching reads went wrong for the following annotation line (so it was skipped): ' +line +'\n\n')
			#print(currentCountDict)

	def split_Array(array, n):
		cache = [[] for i in range(n)]
		c = 0
		for x in array:
			cache[c].append(x)
			c+=1
			if(c > n-1):
				c=0
		return(cache)



#####PARAMETER ERROR HANDLING
	if(len(args.i) > 0):
		for file in args.i:
			if not(pysam.AlignmentFile(file, 'rb').is_bam):
				exit('The following file is either not a valid bam file or it is corrupted: ' +file)
	else:
		exit('You have to provide some .bam mapping files with the -i parameter for counting.')

	if not(args.a):
		exit('You have to provide an annotation file in .gff/.gtf format with the -a parameter for counting.')
	else:
		if not(os.path.isfile(args.a)):
			exit('Could not find the given annotation file: ' +args.a)	

	if not(os.path.isdir(args.o)):
		exit('The given ouput path does not exist:' +args.o)

	if not(os.access(args.o, os.W_OK)):
		exit('No permission to write to the given output path: ' +args.o)

	try:
		if(int(args.cpu) < 1):
			args.cpu = 1
			sys.stderr.write('Amount of used CPUs is set to 1.')
		else:
			args.cpu = int(args.cpu)
	except ValueError:
		args.cpu = 1
		sys.stderr.write('Amount of used CPUs is set to 1.')

	try:
		minOverlap = int(args.O)
		if(minOverlap <= 0):
			minOverlap = 1
			sys.stderr.write('Minimal amount of overlapping bases is set to 1.')
	except ValueError:
		try:
			if(float(args.O) > 0 and float(args.O) <= 1.0):
				minOverlap = float(args.O)
			else:
				exit('If you give the -O/-readOverlap parameter a float X, it has to be 0.0 < X <= 1.0 to be correctly interpreted as the fraction of the given read\'s length.')
		except ValueError:
			exit('You have to give the -O/-readOverlap parameter either a positive Integer or a Float between 0 and 1.')




#################
#################
#################



#####MAIN COUNTING
	if not(args.q):
		sys.stderr.write('\nStart counting for ' +str(len(args.i)) +' given bam files...' +'\n\n')
	annotation = []
	with open(args.a) as givenAnnotation:							#read in the annotation file
		for line in givenAnnotation:
			if(line[0] not in ['#','\n']):
				if(line.split('\t')[2] == args.f):
					annotation.append(line)					

	if(args.cpu > 1):
		splittedAnnotation = split_Array(annotation, args.cpu)
	else:
		splittedAnnotation = [annotation]

	for bamFile in args.i:
		currentCountDict = multiprocessing.Manager().dict()
		processes = []
		lock = multiprocessing.Lock()

		for x in splittedAnnotation:
			p = multiprocessing.Process(target=perform_counting2, args=(bamFile, currentCountDict, x,))		#call counting function
			processes.append(p)
			p.start()
		for p in processes:
			p.join()

		with open(args.o +bamFile.split('/')[-1] +'.count', 'w') as currentFinalOut:
			for name in sorted(currentCountDict.keys()):
				currentFinalOut.write(name +'\t' +str(currentCountDict[name]) +'\n')

		if not(args.q):
			sys.stderr.write('Counting done for file: ' +bamFile.split('/')[-1] +'\n')

	if not(args.q):
		sys.stderr.write('\nFinished counting.\n\n')
#############
#############
#############



#####ANALYZE MODE
elif(sys.argv[1] == 'analyze'):
	parser = argparse.ArgumentParser(description='MeRDE: MicroRNA Differential Expression Analysis -- ANALYZE Mode.', epilog='Version 1.2 (September 2017)')
	parser.add_argument('Mode', help='Start MeRDE in COUNT mode to generate read count files of a given annotation and .bam mapping file(s) OR start MeRDE in ANALYZE mode to calculate differentially expressed miRNAs from given read count files.', choices=['count', 'analyze'])
	parser.add_argument('-i', '-count_files', help='Input read count files. For count file format specifications see Documentation.', metavar='<File>', nargs='*')
	parser.add_argument('-c','-conditions', help='Define to which conditions the given input files belong.', metavar='<String>', nargs='*')
	parser.add_argument('-l', '-low_expression_filter', help='Ignore genes where the mean count of both conditions is lower than a certain threshold. [Default: 10]', metavar='<Integer>', default=10, type=int)
	parser.add_argument('-sigma', help='Set the sigma factor threshold for calculating the gene expression clusters. Should only be changed with caution. [Default: 2]', metavar='<Float>', default=2, type=float)
	parser.add_argument('-s', '-minimal_cluster_size', help='Set the minimal cluster size threshold for estimating the gene distribution parameters. Should only be changed with caution. [Default: 20]', metavar='<Integer>', default=20, type=int)
	parser.add_argument('-n', '-normalize', help='Normalize the count read files. Should only be applied to raw read counts. [Default: True]', default=True, choices=['True', 'False'])
	parser.add_argument('-e', '-outlier', help='Choose how to treat outliers, either "ignore", or "remove" or "replace". [Default="remove"]', default='remove', choices=['remove','ignore','correct', 'test'])
	parser.add_argument('-m', '-minimum_size', help='Minimum amount of Replicates per Condition.". [Default=4]', metavar='<Integer>', default=4, type=int)
	parser.add_argument('-no_html', help='Do not create a HTML results page. Will be faster since no pictures are created.', action='store_true')
	parser.add_argument('-cpu' , help='Amount of CPUs to use. [Default: 1]', metavar='<Integer>', default=1)
	parser.add_argument('-o', '-output_dir', help='Directory where output will be saved. [Default: current working directory]', metavar='<Path>', default=os.getcwd() +'/')
	parser.add_argument('-q', '-quiet', help='Do not write status information into console.', action='store_true')
	args = parser.parse_args()


	####FUNCTIONS
	def check_single_count_file(path):
		currentCountFile = open(path).readlines()
		geneNames = []
		geneCounts = []
		for line in currentCountFile:
			if(line[0] != '#'):
				lineSplit = line.split('\t')
				if(lineSplit[0] not in geneNames):
					geneNames.append(lineSplit[0])
					try:
						geneCounts.append(int(lineSplit[1].strip()))
					except ValueError:
							exit('The following line in file ' +path +'contains a non-integer value in a count column:\n' +line)
				else:
					exit('The gene name ' +lineSplit[0] +' in the file ' +path +' is duplicated. Make sure that every gene name is unique, or I will get confused during the analysis later on.')
		return((geneCounts, geneNames))

	
	def check_count_matrix_file(path):
		currentCountFile = open(path).readlines()
		geneNames = []
		geneCountsMatrix = [[] for i in range(len(args.c))]
		for line in currentCountFile:
			if(line[0] != '#'):
				lineSplit = line.split('\t')
				if(len(lineSplit)-1 != len(args.c)):
					exit('The amount of given conditions (' +str(len(args.c)) +') does not match the number of count columns of the following line:\n' +line)
				if(lineSplit[0] in geneNames):
					exit('The gene name ' +lineSplit[0] +' is duplicated. Make sure that every gene name is unique, or I will get confused during the analysis later on.')
				geneNames.append(lineSplit[0])
				for i in range(1, len(lineSplit)):
					try:
						geneCountsMatrix[i-1].append(int(lineSplit[i].strip()))
					except ValueError:
						exit('The following line in file ' +path +'contains a non-integer value in a count column:\n' +line)
		return(geneCountsMatrix, geneNames)

	
	def normalize_count_matrix(cMatrix):
		libSizeFactors = []
		for sample in range(len(cMatrix)):
			currSizeFactors = []
			for geneRow in range(len(cMatrix[sample])):
				if(int(cMatrix[sample][geneRow]) != 0):
					currRow = []
					for j in range(len(cMatrix)):
						if(int(cMatrix[j][geneRow]) != 0):
							currRow.append(int(cMatrix[j][geneRow]))
					currSizeFactors.append(int(cMatrix[sample][geneRow]) / gmean(currRow))
			libSizeFactors.append(np.median(currSizeFactors))
		
		normMatrix = [[] for i in range(len(cMatrix))]
		for i in range(len(cMatrix)):
			for j in range(len(cMatrix[i])):
				normMatrix[i].append(int(cMatrix[i][j]) / libSizeFactors[i])
		return(normMatrix, libSizeFactors)

	
	def filter_count_matrix(normCMatrix):
		filteredLowNormCMatrix = {conditionA:dict(), conditionB:dict()}
		lowExpressedNames = []

		for gene in normCMatrix[conditionA]:
			if(np.mean(normCMatrix[conditionA][gene]) >= args.l or np.mean(normCMatrix[conditionB][gene]) >= args.l):
				filteredLowNormCMatrix[conditionA][gene] = normCMatrix[conditionA][gene]
				filteredLowNormCMatrix[conditionB][gene] = normCMatrix[conditionB][gene]
			else:
				lowExpressedNames.append(gene)
		return(filteredLowNormCMatrix, lowExpressedNames)


	def calculate_gene_clusters(someCountMatrix):
		geneClusters = {conditionA:dict(), conditionB:dict()}

		for condition in someCountMatrix:
			for gene in someCountMatrix[condition]:
				if(someCountMatrix[conditionA][gene] != False and someCountMatrix[conditionB][gene] != False):		#If this gene was filtered in any of the conditions, it does not need a cluster
					geneClusters[condition][gene] = []
					mean = np.mean(someCountMatrix[condition][gene])
					std = np.std(someCountMatrix[condition][gene])

					clusterRange = (mean-(std*args.sigma/2), mean+(std*args.sigma))

					for gene2 in someCountMatrix[condition]:
						if(someCountMatrix[condition][gene2] != False):			# if a gene was completely filtered out, it should not be part of any cluster
							if(clusterRange[1] < args.l and np.mean(someCountMatrix[condition][gene2])+(np.std(someCountMatrix[condition][gene2])*args.sigma) < args.l):
								geneClusters[condition][gene].append(gene2)
							elif(np.mean(someCountMatrix[condition][gene2]) >= clusterRange[0] and np.mean(someCountMatrix[condition][gene2]) <= clusterRange[1]):
								geneClusters[condition][gene].append(gene2)
		return(geneClusters)

	
	def get_cluster_stats(geneCluster, condition, gene, someCountMatrix):
		clusterCounts = get_cluster_counts(geneCluster, condition, gene, someCountMatrix)
		return(np.mean(clusterCounts), np.std(clusterCounts))

	
	def get_cluster_counts(geneCluster, condition, gene, someCountMatrix):
		clusterCounts = []
		for geneID in geneCluster[condition][gene]:
			for count in someCountMatrix[condition][geneID]:
				clusterCounts.append(count)
		return(clusterCounts)

	
	def remove_outliers(filteredLowNormCMatrix, geneCluster):
		outMatrix = {x:{y:{} for y in filteredLowNormCMatrix[conditionA]} for x in [conditionA, conditionB]}
		'''
		for con in filteredLowNormCMatrix:
			for gen in filteredLowNormCMatrix[con]:
				outMatrix[con][gen] = own_cluster_outlier(gen, con, geneCluster, filteredLowNormCMatrix)
				outMatrix[con][gen] = own_samples_outlier(outMatrix[con][gen])
		return(outMatrix)
		'''

		'''
		print('Filter Sample outliers...')
		for con in filteredLowNormCMatrix:
			for gen in filteredLowNormCMatrix[con]:
				outMatrix[con][gen] = own_samples_outlier2(gen, con, filteredLowNormCMatrix)

		print('Calculate intermediate clusters...')
		intermediateClusters = calculate_gene_clusters(outMatrix)

		print('Filter cluster outliers...')
		outMatrix2 = {x:{y:{} for y in outMatrix[conditionA]} for x in [conditionA, conditionB]}
		for con in outMatrix:
			for gen in outMatrix[con]:
				outMatrix2[con][gen] = own_cluster_outlier(gen, con, intermediateClusters, outMatrix)
		return(outMatrix2)
		'''
		for con in filteredLowNormCMatrix:
			for gen in filteredLowNormCMatrix[con]:
				outMatrix[con][gen] = double_mad_outlier(gen, con, geneCluster, filteredLowNormCMatrix)
		return(outMatrix)


	def mean_std_ratio_test(filteredLowNormCMatrix):
		cache = {x:{y:{z:0 for z in ['mean', 'meanLog', 'msRatio', 'residual','zResidual']} for y in filteredLowNormCMatrix[conditionA]} for x in [conditionA, conditionB]}
		outMatrix = {x:{y:{} for y in filteredLowNormCMatrix[conditionA]} for x in [conditionA, conditionB]}
		valuesDict = {x:[] for x in ['mean', 'meanLog', 'msRatio', 'residual', 'name', 'std']}

		meanLog = []
		msRatio = []
		for con in filteredLowNormCMatrix:
			for gen in filteredLowNormCMatrix[con]:
				cache[con][gen]['mean'] = np.mean(filteredLowNormCMatrix[con][gen])

				if(np.sum(filteredLowNormCMatrix[con][gen]) > 0):
					valuesDict['name'].append(f'{gen}_{con}')
					valuesDict['mean'].append(np.mean(filteredLowNormCMatrix[con][gen]))
					valuesDict['meanLog'].append(np.log10(np.mean(filteredLowNormCMatrix[con][gen])))
					valuesDict['std'].append(np.std(filteredLowNormCMatrix[con][gen]))
					valuesDict['msRatio'].append(np.mean(filteredLowNormCMatrix[con][gen]) / np.std(filteredLowNormCMatrix[con][gen]))

				if(cache[con][gen]['mean'] >= args.l and np.sum(filteredLowNormCMatrix[con][gen]) > 0):
					cache[con][gen]['msRatio'] = np.mean(filteredLowNormCMatrix[con][gen]) / np.std(filteredLowNormCMatrix[con][gen])
					cache[con][gen]['meanLog'] = np.log10(np.mean(filteredLowNormCMatrix[con][gen]))
					meanLog.append(np.log10(cache[con][gen]['mean']))
					msRatio.append(cache[con][gen]['msRatio'])
				else:
					cache[con][gen]['msRatio'] = np.mean(filteredLowNormCMatrix[con][gen]) / (np.std(filteredLowNormCMatrix[con][gen])+0.1)


		linReg = np.polyfit(meanLog, msRatio, 1)
		linRegFunc = np.poly1d(linReg)

		negResiduen = []

		for con in cache:
			for gen in cache[con]:
				if(np.sum(filteredLowNormCMatrix[con][gen]) > 0):
					valuesDict['residual'].append(cache[con][gen]['msRatio'] - linRegFunc(np.log10(cache[con][gen]['mean'])))

				if(cache[con][gen]['mean'] >= args.l and np.sum(filteredLowNormCMatrix[con][gen]) > 0):
					cache[con][gen]['residual'] = cache[con][gen]['msRatio'] - linRegFunc(np.log10(cache[con][gen]['mean']))
					if(cache[con][gen]['residual'] < 0):
						negResiduen.append(cache[con][gen]['residual'])

		negResiduen = np.array(negResiduen)

		for con in cache:
			for gen in cache[con]:
				if(cache[con][gen]['residual'] < 0):
					currZnormedResidual = (cache[con][gen]['residual'] - np.mean(negResiduen)) / np.std(negResiduen) #Z-normalization of negative residuals
					if(currZnormedResidual < -2): 
						outMatrix[con][gen] = own_samples_outlier(filteredLowNormCMatrix[con][gen], 1)
					elif(currZnormedResidual < -1): 
						outMatrix[con][gen] = own_samples_outlier(filteredLowNormCMatrix[con][gen], 1.1)
					else:
						outMatrix[con][gen] = filteredLowNormCMatrix[con][gen]
				else:
					outMatrix[con][gen] = filteredLowNormCMatrix[con][gen]

		### Test output
		pickle.dump((cache,meanLog,msRatio), open(f'{args.o}/cache.pydict', 'wb'))

		draw_msRatio_plot(valuesDict, linRegFunc)

		return(cache, outMatrix)


	def own_samples_outlier(counts, fac):
		nonOutliers = []
		outliers = []
		cValues = np.array(counts)**(1/3)																	#Rough transform to normal distribution values
		samplesMean, samplesSTD = np.mean(cValues), np.std(cValues)

		for i,countValue in enumerate(cValues):
			if(countValue < (samplesMean-(samplesSTD*fac)) or countValue > (samplesMean+(samplesSTD*fac))):				#within 1 standard deviations of the mean
				outliers.append(counts[i])
			else:
				nonOutliers.append(counts[i])
		return(nonOutliers)


	def double_mad_outlier(gene, condition, geneCluster, filteredLowNormCMatrix):
		cValuesGene = filteredLowNormCMatrix[condition][gene]
		if(np.mean(cValuesGene) <= args.l):													 
			return(cValuesGene)
		cValues = get_cluster_counts(geneCluster, condition, gene, filteredLowNormCMatrix)
		med = np.median(cValues)
		absMedDev = [abs(x - med) for x in cValues]
		leftMad = np.median([x for i,x in enumerate(absMedDev) if cValues[i] <= med])
		rightMad = np.median([x for i,x in enumerate(absMedDev) if cValues[i] >= med])

		madDistances = []
		for x in cValues:
			if(x <= med):
				madDistances.append(abs((x - med)) / leftMad)
			else:
				madDistances.append(abs((x - med)) / rightMad)

		out = []
		for i,x in enumerate(madDistances):
			if(x > 3):
				out.append(cValues[i])
		
		newCValuesGene = [x for x in cValuesGene if x not in out]
		#if(len(newCValuesGene) < round(len(cValuesGene)*0.7)):					#If too much count values were removed (at the moment per default 4 replicates have to remain)
		if(len(newCValuesGene) < args.m):
			return(False)														#don't consider gene for diff analysis

		return(newCValuesGene)


	def estimate_distribution_parameters(outlierFilteredNormCMatrix, geneClustersRecalculated):
		outParameters = {conditionA:dict(), conditionB:dict()}
		for condition in outlierFilteredNormCMatrix:
			for gene in outlierFilteredNormCMatrix[condition]:
				if(outlierFilteredNormCMatrix[conditionA][gene] != False and outlierFilteredNormCMatrix[conditionB][gene] != False):	#If gene was filtered out in any condition, don't bother with estimating parameters
					cache = get_cluster_counts(geneClustersRecalculated, condition, gene, outlierFilteredNormCMatrix)
					currentClusterCounts = []
					for x in cache:
						if(x > 0):
							currentClusterCounts.append(x)
					
					if(len(currentClusterCounts) >= args.s):
						N = len(currentClusterCounts)
						sApprox = np.log((1/N)*sum(currentClusterCounts)) - ((1/N)*(sum(np.log(currentClusterCounts))))
						pApprox = (3 - sApprox + np.sqrt(((sApprox - 3)**2) + (24*sApprox))) / (12 * sApprox)
						bApprox = (1 / (pApprox * N)) * sum(currentClusterCounts)
						outParameters[condition][gene] = [sApprox, pApprox, bApprox, pApprox * bApprox, np.sqrt(pApprox * (bApprox**2)), N] #0: Approximation of s; 1: Approximation of shape parameter p; 2: Approximation of scale parameter b; 3: Approximation of mean; 4: Approximation of standard deviation; 5: Amount of counts in cluster
					else:
						outParameters[condition][gene] = False
				else:
					outParameters[condition][gene] = False
		
		return(outParameters)


	def boxplot(geneID, outlierFilteredNormCMatrix):
		fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
		bpCubic = axes.boxplot([outlierFilteredNormCMatrix[conditionA][geneID], outlierFilteredNormCMatrix[conditionB][geneID]], vert=True, patch_artist=True)

		colors = ['darkgreen', 'darkred']
		for patch, color in zip(bpCubic['boxes'], colors):
			patch.set_facecolor(color)

		axes.yaxis.grid(True)
		axes.set_xticks([y+1 for y in range(2)], )
		axes.set_ylabel('Normalized counts')

		plt.setp(axes, xticks=[y+1 for y in range(2)], xticklabels=[conditionA, conditionB])
		plt.savefig(f'{args.o}/boxplots/{geneID}.png')
		plt.close(fig)


	def boxplot2(geneID, outlierFilteredNormCMatrix, filteredLowNormCMatrix):
		fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7.5, 7.5))
		bpCubic = axes[0].boxplot([outlierFilteredNormCMatrix[conditionA][geneID], outlierFilteredNormCMatrix[conditionB][geneID]], vert=True, patch_artist=True)

		colors = ['darkgreen', 'darkred']
		for patch, color in zip(bpCubic['boxes'], colors):
			patch.set_facecolor(color)

		axes[0].yaxis.grid(True)
		axes[0].set_xticks([y+1 for y in range(2)], )
		axes[0].set_ylabel('Normalized counts')

		plt.setp(axes[0], xticks=[y+1 for y in range(2)], xticklabels=[conditionA, conditionB])

		bpCubic = axes[1].boxplot([filteredLowNormCMatrix[conditionA][geneID], filteredLowNormCMatrix[conditionB][geneID]], vert=True, patch_artist=True)

		colors = ['darkgreen', 'darkred']
		for patch, color in zip(bpCubic['boxes'], colors):
			patch.set_facecolor(color)

		axes[1].yaxis.grid(True)
		axes[1].set_xticks([y+1 for y in range(2)], )
		axes[1].set_ylabel('Normalized counts')

		plt.setp(axes[1], xticks=[y+1 for y in range(2)], xticklabels=[conditionA, conditionB])
		plt.tight_layout()
		plt.savefig(f'{args.o}/boxplots/{geneID}.png')
		plt.close(fig)


	def boxplot3(geneID, someCountMatrix, pV):
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7.5, 7.5) ,sharey=True)

		sns.set_style("whitegrid")
		plt.rcParams["font.family"] = "TeX Gyre Schola"

		vp = ax[0].violinplot(someCountMatrix[conditionA][geneID], vert=True, showmeans=True, showextrema=False, bw_method=.3)
		vp['bodies'][0].set_facecolor('#006d2c')
		vp['bodies'][0].set_alpha(.2)
		vp['cmeans'].set_edgecolor('#006d2c')
		ax[0].set_ylabel('Normalized Counts', fontname='TeX Gyre Schola', fontsize=15)
		ax[0].set_xlabel('annual', fontname='TeX Gyre Schola', fontsize=15)
		ax[0].set_xticks([0])
		ax[0].set_ylim([0, max(someCountMatrix[conditionA][geneID]+someCountMatrix[conditionB][geneID])*1.05])



		vp = ax[1].violinplot(someCountMatrix[conditionB][geneID], vert=True, showmeans=True, showextrema=False, bw_method=.3)
		vp['bodies'][0].set_facecolor('#a50f15')
		vp['bodies'][0].set_alpha(.2)
		vp['cmeans'].set_edgecolor('#a50f15')
		ax[1].set_xlabel('non-annual', fontname='TeX Gyre Schola', fontsize=15)
		ax[1].set_xticks([0])

		ax[0].plot([1]*len(someCountMatrix[conditionA][geneID]), someCountMatrix[conditionA][geneID], 'o', color='#006d2c')
		ax[1].plot([1]*len(someCountMatrix[conditionB][geneID]), someCountMatrix[conditionB][geneID], 'o', color='#a50f15')
		fig.text(0.5, 0.9, geneID, ha='center', va='center', fontname='TeX Gyre Schola', fontsize=15)
		fig.text(0.5, 0.05, f'adj. p-Value = {pV}', ha='center', va='center', fontname='TeX Gyre Schola', fontsize=10)

		plt.savefig(f'{args.o}/boxplots/{geneID}.png')
		plt.close(fig)


	def draw_msRatio_plot(dataDict, linRegFunc):
		output_file("/mnt/prostlocal/emanuel/DA/diapause/DEG/merde/test_scenario_1/meanLog_meanstd.html")
		source = ColumnDataSource(data=dataDict)
		view = CDSView(source=source, filters=[BooleanFilter([True if dataDict['meanLog'][i] >= 1 and dataDict['msRatio'][i] <= 1 else False for i in range(len(dataDict['name']))])])
		view2 = CDSView(source=source, filters=[BooleanFilter([True if dataDict['residual'][i] <= -0.8  else False for i in range(len(dataDict['name']))])])
		view3 = CDSView(source=source, filters=[BooleanFilter([True if dataDict['residual'][i] <= 0  else False for i in range(len(dataDict['name']))])])

		tools = ["box_select", "reset", 'box_zoom', 'pan']
		p = figure(title='Gene expression Mean versus Mean/STD ratio', plot_height=1000, plot_width=1000, tools=tools)
		p.circle(x='meanLog', y="msRatio", size=10, color='green', hover_color="blue", source=source)
		p.circle(x='meanLog', y="msRatio", size=10, color='yellow', hover_color="blue", source=source, view=view3)
		p.circle(x='meanLog', y="msRatio", size=10, color='orange', hover_color="blue", source=source, view=view)
		p.circle(x='meanLog', y="msRatio", size=10, color='red', hover_color="blue", source=source, view=view2)

		regLine = [linRegFunc(i) for i in np.linspace(0, 6, 30)] 
		p.line([i for i in np.linspace(0, 6, 30)], regLine, line_width=5, line_color='black', line_alpha=0.5)

		p.xaxis.axis_label = "Log10 Mean"
		p.yaxis.axis_label = "Mean / STD ratio"

		p.add_tools(HoverTool(tooltips = [
		    ("Name", "@name"),
		    ("mean", "@mean{0,0.00}"),
		    ("mean Log", "@meanLog{0,0.00}"),
		    ("STD", "@std{0,0.00}"),
		    ("Mean/STD ratio", "@msRatio{0,0.00}"),
		    ("Residual", "@residual{0,0.00}")]))

		show(p)


	def get_p_Value_gamma_test(geneID, outlierFilteredNormCMatrix, estimatedGeneDistributionParameters):
		#0: Approximation of s; 1: Approximation of shape parameter p; 2: Approximation of scale parameter b; 3: Approximation of mean; 4: Approximation of standard deviation; 5: Amount of counts in cluster
		shapeConditionA, scaleConditionA, meanConditionA, stdConditionA = estimatedGeneDistributionParameters[conditionA][geneID][1:5]
		shapeConditionB, scaleConditionB, meanConditionB, stdConditionB = estimatedGeneDistributionParameters[conditionB][geneID][1:5]
		meanConditionA = np.mean(outlierFilteredNormCMatrix[conditionA][geneID])
		meanConditionB = np.mean(outlierFilteredNormCMatrix[conditionB][geneID])

		zScoreConditionA = abs(((meanConditionA - meanConditionB) / stdConditionA) * np.sqrt(len(outlierFilteredNormCMatrix[conditionA][geneID])))
		zScoreConditionB = abs(((meanConditionB - meanConditionA) / stdConditionB) * np.sqrt(len(outlierFilteredNormCMatrix[conditionB][geneID]))) 

		pValueConditionA = abs(1 - quad(gamma.pdf, -zScoreConditionA+100, zScoreConditionA+100, args=(shapeConditionB, 100, scaleConditionB/stdConditionB))[0])
		pValueConditionB = abs(1 - quad(gamma.pdf, -zScoreConditionB+100, zScoreConditionB+100, args=(shapeConditionA, 100, scaleConditionA/stdConditionA))[0])

		if(pValueConditionA == np.inf and pValueConditionB == np.inf):
			return(get_p_value_cubic_Welsch_test(geneID, outlierFilteredNormCMatrix))

		return(min(pValueConditionA,pValueConditionB))


	def get_p_value_cubic_Welsch_test(geneID, outlierFilteredNormCMatrix):
		return( ttest_ind(np.array(outlierFilteredNormCMatrix[conditionA][gene])**(1/3), np.array(outlierFilteredNormCMatrix[conditionB][gene])**(1/3), equal_var=False)[1] )


##### PARAMETER ERROR HANDLING
	if not(args.c):
		exit('Conditions of the given input files not given, use the -c option to define them.')
	else:
		if(len(args.i) != len(args.c) and len(args.i) != 1):
			exit('Number of given input files (' +str(len(args.i)) +') does not match the number of given conditions (' +str(len(args.c)) +').') 
		conditionSet = set()
		for x in args.c:
			conditionSet.add(x)
		if(len(conditionSet) != 2):
			exit('You have to provide two different conditions, not more or less. Sorry!')
		conditionA = args.c[0]
		conditionB = list(conditionSet - {conditionA})[0]
		
	if not(os.path.isdir(args.o)):
		exit('The given ouput path does not exist:' +args.o)

	if not(os.access(args.o, os.W_OK)):
		exit('No permission to write to the given output path: ' +args.o)

	try:
		if(int(args.cpu) < 1):
			args.cpu = 1
			sys.stderr.write('Amount of used CPUs is set to 1.')
		else:
			args.cpu = int(args.cpu)
	except ValueError:
		args.cpu = 1
		sys.stderr.write('Amount of used CPUs is set to 1.')

	try:
		if(int(args.l) < 1):
			args.l = 0
			sys.stderr.write('Low expression filter set to 0 (Not recommended).')
		else:
			args.l = int(args.l)
	except ValueError:
		args.l = 10
		sys.stderr.write('Low expression filter set to 10.')

	if(args.n == 'True' or args.n == True):
		args.n = True
	else:
		args.n = False


#####################################################################################
#####################################################################################
#####################################################################################

#################### MAIN ANALYZE
	
############## Check Count Input format
	if not(args.q):
		sys.stderr.write('\nChecking count files... \n')
	if(len(args.i) == 1):
		countMatrix, geneNames = check_count_matrix_file(args.i[0])			#Input was a count matrix file. not single count files
	else:
		countMatrix = []													#Input were several single count files, each belonging to one sample
		geneNames = []
		for file in args.i:
			cache = check_single_count_file(file)							#check if the count file's format is correct
			countMatrix.append(cache[0])
			if(geneNames != []):
				if(geneNames != cache[1]):									#if the order of genes is not the same in all count files, abort
					exit('The order or naming of the genes of the following file does not fit the other files:\n' +file +'\nPlease make sure that the naming and order of genes is the same in every given coutn file.')
			geneNames = cache[1]
	if not(args.q):
		sys.stderr.write('Count files seem to be alright... \n')
##############
##############



############## Normalize Count Matrix
	if(args.n):
		if not(args.q):
			sys.stderr.write('Normalizing counts... \n')
		normalizedCountMatrix, libSizeFactors = normalize_count_matrix(countMatrix)				#unfiltered (may contain outliers and low expressed genes) normalized counts 
	else:
		normalizedCountMatrix, libSizeFactors = countMatrix, [0 for i in range(len(countMatrix))]		#unfiltered (may contain outliers and low expressed genes) non-normalized counts

	normCMatrix = {conditionA:dict(), conditionB:dict()}								#unfiltered (may contain outliers and low expressed genes) dictionary holding every gene and their counts separated by the two conditions
	for j in range(len(normalizedCountMatrix)):
		for i in range(len(normalizedCountMatrix[j])):
			if(geneNames[i] not in normCMatrix[args.c[j]]):
				normCMatrix[args.c[j]][geneNames[i]] = []
			normCMatrix[args.c[j]][geneNames[i]].append(normalizedCountMatrix[j][i])
##############
##############



############## Filter Low expressed Genes
	if not(args.q):
		sys.stderr.write('Filtering genes with a mean count of less than ' +str(args.l) +' in both conditions... \n')
	filteredLowNormCMatrix, lowExpressedNames = filter_count_matrix(normCMatrix)							#filtered (may contain outliers) dictionary holding every gene and their counts separeted by the two conditions

	### Test output
	pickle.dump(filteredLowNormCMatrix, open(f'{args.o}/filteredLowNormCMatrix.pydict', 'wb'))


##############
##############



############## Check for Genes with a severly disturbed mean/std ratio, they have to be handeled before building clusters
	if(args.e != 'ignore'):
		if not(args.q):
			sys.stderr.write('Checking for weird genes... \n')
		msRatioMatrix, filteredLowNormCMatrix = mean_std_ratio_test(filteredLowNormCMatrix)
		#pickle.dump(filteredLowNormCMatrix, open(f'{args.o}/filteredLowNormCMatrix.pydict', 'wb'))
		msRatioMatrix, filteredLowNormCMatrix = mean_std_ratio_test(filteredLowNormCMatrix)

##############
##############



############## Calculate Gene Clusters based on non-outlier-filtered Count Matrix
	if not(args.q):
		sys.stderr.write(str(len(normCMatrix[conditionA]) - len(filteredLowNormCMatrix[conditionA])) +' of ' +str(len(normCMatrix[conditionA])) +' genes were filtered... \n')
		sys.stderr.write('Calculating gene clusters... \n')
	geneClusters = calculate_gene_clusters(filteredLowNormCMatrix)

	### Test output
	pickle.dump(geneClusters, open(f'{args.o}/geneClusters.pydict', 'wb'))

##############
##############



############## Outlier filtering
	if(args.e == 'ignore'):
		if not(args.q):
			sys.stderr.write('Outlier detection is set to "ignore"... \n')

	if(args.e != 'ignore'):
		if not(args.q):
			sys.stderr.write('Removing outliers based on cluster properties... \n')
		
		if(args.e == 'remove'):
			outlierFilteredNormCMatrix = remove_outliers(filteredLowNormCMatrix, geneClusters)
			
			### Test output
			pickle.dump(outlierFilteredNormCMatrix, open(f'{args.o}/outlierFilteredNormCMatrix.pydict', 'wb'))

			if not(args.q):
				sys.stderr.write('Recalculating gene clusters... \n')
			geneClustersRecalculated = calculate_gene_clusters(outlierFilteredNormCMatrix) 

		if(args.e == 'correct'):
			outlierFilteredNormCMatrix = correct_outliers(filteredLowNormCMatrix, geneClusters)
			if not(args.q):
				sys.stderr.write('Recalculating gene clusters... \n')
			geneClustersRecalculated = calculate_gene_clusters(outlierFilteredNormCMatrix) 
	else:
		geneClustersRecalculated = geneClusters
##############
##############



############## Estimating Distribution Parameters
	if not(args.q):
		sys.stderr.write('Estimating distribution parameters... \n')

	estimatedGeneDistributionParameters = estimate_distribution_parameters(outlierFilteredNormCMatrix, geneClustersRecalculated)
##############
##############



############## Check for differential gene expresssion
	if not(args.q):
		sys.stderr.write('Checking for differentially expressed genes... \n')

	pValues = []
	geneNames = []
	outFilteredNames = []
	for gene in outlierFilteredNormCMatrix[conditionA]:
		if(outlierFilteredNormCMatrix[conditionA][gene] != False and outlierFilteredNormCMatrix[conditionB][gene] != False):	#If gene was filtered out in any condition, don't bother with estimating parameters
			geneNames.append(gene)
			if(estimatedGeneDistributionParameters[conditionA][gene] != False and estimatedGeneDistributionParameters[conditionB][gene] != False):
				pValues.append(get_p_Value_gamma_test(gene, outlierFilteredNormCMatrix, estimatedGeneDistributionParameters))
			else:
				pValues.append(get_p_value_cubic_Welsch_test(gene, outlierFilteredNormCMatrix))
		else:
			outFilteredNames.append(gene)

	adjustedPvalues = multipletests(pValues, alpha=0.1, method='fdr_bh')
##############
##############



############## Write Results 
	if not(args.q):
		sys.stderr.write(f'Writing output file... \n')

	result = []
	for i in range(len(geneNames)):
		result.append((geneNames[i], pValues[i], adjustedPvalues[1][i]))
	result.sort(key = lambda x: (x[2],x[1]))

	with open(f'{args.o}/merde3_results.csv', 'w') as f:
		f.write('Name\tlog2FC\tmean1\tmean2\tinv. Counts 1\tinv. Counts 2\tp-Value\tadj. p-Value\n')
		for x in result:
			meanA = np.mean(outlierFilteredNormCMatrix[conditionA][x[0]]) + 1
			meanB = np.mean(outlierFilteredNormCMatrix[conditionB][x[0]]) + 1
			logFC = np.log2(meanA / meanB)
			f.write(f'{x[0]}\t{logFC}\t{meanA}\t{meanB}\t{len(outlierFilteredNormCMatrix[conditionA][x[0]])}({len(normCMatrix[conditionA][x[0]])})\t{len(outlierFilteredNormCMatrix[conditionB][x[0]])}({len(normCMatrix[conditionB][x[0]])})\t{x[1]}\t{x[2]}\n')

		for x in outFilteredNames:
			meanA = np.mean(normCMatrix[conditionA][x]) + 1
			meanB = np.mean(normCMatrix[conditionB][x]) + 1
			logFC = np.log2(meanA / meanB)
			f.write(f'{x}\t{logFC}\t{meanA}\t{meanB}\t{len(filteredLowNormCMatrix[conditionA][x])}({len(normCMatrix[conditionA][x])})\t{len(filteredLowNormCMatrix[conditionB][x])}({len(normCMatrix[conditionB][x])})\tNA\tNA\n')

		for x in lowExpressedNames:
			meanA = np.mean(normCMatrix[conditionA][x]) + 1
			meanB = np.mean(normCMatrix[conditionB][x]) + 1
			logFC = np.log2(meanA / meanB)

			f.write(f'{x}\t{logFC}\t{meanA}\t{meanB}\t{len(normCMatrix[conditionA][x])}({len(normCMatrix[conditionA][x])})\t{len(normCMatrix[conditionB][x])}({len(normCMatrix[conditionB][x])})\tNA\tNA\n')
##############
##############

	### Test output
	pickle.dump(outlierFilteredNormCMatrix, open(f'{args.o}/outlierFilteredNormCMatrix.pydict', 'wb'))
	pickle.dump(geneClustersRecalculated, open(f'{args.o}/geneClustersRecalculated.pydict', 'wb'))
	pickle.dump(estimatedGeneDistributionParameters, open(f'{args.o}/estimatedGeneDistributionParameters.pydict', 'wb'))


############## Write HTML Output
	def get_results_data():
		e = open(f'{args.o}/merde3_results.csv').readlines()
		resDataDict = {x:[] for x in ['name', 'mean', 'meanA', 'meanB', 'fc', 'logFC', 'p', 'ap', 'noSamplesA', 'noSamplesB']}

		for lin in e[1:]:
			nil = lin.split('\t')
			resDataDict['name'].append(nil[0])
			resDataDict['mean'].append(np.log10(np.mean((float(nil[2]),float(nil[3])))+0.1))
			resDataDict['meanA'].append(float(nil[2]))
			resDataDict['meanB'].append(float(nil[3]))
			resDataDict['fc'].append(2**float(nil[1]) if nil[1] != 'NA' else 0)
			resDataDict['logFC'].append(float(nil[1]) if nil[1] != 'NA' else 0)
			resDataDict['p'].append(nil[-2])
			resDataDict['ap'].append(nil[-1].strip())
			resDataDict['noSamplesA'].append(nil[4])
			resDataDict['noSamplesB'].append(nil[5])

		return(resDataDict)

	def draw_mfc_plot(resData):
		source = ColumnDataSource(data=resData)
		sigGenes = CDSView(source=source, filters=[BooleanFilter([True if resData['ap'][i] != 'NA' and float(resData['ap'][i]) <= 0.05 else False for i in range(len(resData['ap']))])])

		tools = ["box_select", "reset", 'box_zoom', 'pan']
		p = figure(title='MA-Plot', plot_height=1000, plot_width=1000, tools=tools)
		p.circle(x='mean', y="logFC", size=10, color='black', hover_color="blue", source=source)
		p.circle(x='mean', y="logFC", size=10, color='green', hover_color="blue", source=source, view=sigGenes)

		p.xaxis.axis_label = "Log10 mean of normalized counts"
		p.yaxis.axis_label = "Log2 fold change"

		p.add_tools(HoverTool(tooltips = [
		    ("Name", "@name"),
		    ("Log10 Mean (total)", "@mean{0,0.00}"),
		    (f"Mean ({conditionA})", "@meanA{0,0.00}"),
		    (f"Mean ({conditionB})", "@meanB{0,0.00}"),
		    ("Log2 Fold change", "@logFC{0,0.00}"),
		    ("Fold change", "@fc{0,0.00}"),
		    ("adj. p-Value", "@ap{0,0.00}")]))
		return(p)


	if(args.no_html == False):
		if not(args.q):
			sys.stderr.write('Writing HTML output (this may take a while)... \n')
		
		from base64 import b64encode
		import matplotlib.pyplot as plt
		
		if not(os.path.exists(f'{args.o}/boxplots/')):
			os.mkdir(f'{args.o}/boxplots/')

		if not(os.path.exists(f'{args.o}/html/')):
			os.mkdir(f'{args.o}/html/')


		output_file(f'{args.o}/html/ma_plot.html')

		resData = get_results_data()
		save(draw_mfc_plot(resData))

		for gen in result:
			boxplot3(gen[0], outlierFilteredNormCMatrix, gen[2])
		#still to do
	##############
	##############
	##############


	if not(args.q):
		sys.stderr.write('Analysis finished!\n')