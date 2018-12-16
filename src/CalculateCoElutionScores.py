from __future__ import division

import mmap
import numpy as np
import scipy.stats
import sys
import math
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from scipy.spatial import distance
import operator
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support
from sklearn.model_selection import  cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import RFECV

from sklearn import svm, linear_model
from sklearn import metrics
import multiprocessing as mp

# The following imported libraries are for intergrating GeneMANIA data as Functional Evidence
import glob
import sys
from collections import defaultdict
from bs4 import BeautifulSoup
import urllib
import urllib2
import os 
import re
import gzip
import requests

#from keras.models import Sequential, Model
#from keras.layers import Dense, Dropout, Input
#from keras import regularizers

np.random.seed(1)


# GLobal static objects used across the files
# Load R packages globally in order to avoid loading each package again in each thread

# Load script for calculating Bayes correlation globally
r=robjects.r
r.source(os.path.realpath(__file__).rsplit(os.sep,1)[0] + os.sep + "Bayes_Corr.R")

cor1 = robjects.r["Bayes_Corr_Prior1"]
cor2 = robjects.r["Bayes_Corr_Prior2"]
cor3 = robjects.r["Bayes_Corr_Prior3"]
current_bayes_cor = ""

# Packages required for calculating WCC
rpackages.importr('wccsom')
r_wcc = robjects.r['wcc']

# array for storing elution matrices with poission noise for PCC + Noise co-elution freature
Poisson_cor_Mat = []
Bayes_cor_Mat	= []

# Default number of Threads is number of available cores
num_cores = mp.cpu_count()

#Global variables for calculating scores across mutliple threads
prot2profile = {}
ppiList = []

def lineCount(filename):
	lines = 0
	fh = open(filename)
	for l in fh:
		lines += 1
	return lines


# Helper functions for calculating co-elution scores on mutliple processors at the same time
def getScore(inputQ, outputList):
	for args in iter(inputQ.get, 'STOP'):
		global prot2profile
		ppi, score_index, fname, scoreFun = args
		protA, protB = ppi.split("\t")
		if protA in prot2profile[fname] and protB in prot2profile[fname]:
			profileA = prot2profile[fname][protA]
			profileB = prot2profile[fname][protB]
			score = scoreFun(profileA, profileB)
			outputList.put((score_index, score))
		else:
			outputList.put((score_index, 0.0))
		inputQ.task_done()





#Contacts:
#	Florian Goebels - florian.goebels@gmail.com
#
# @author: Florian Goebels
# This class stores and maintains co elution data
# It stores both the raw and normalizied co-elution data
class ElutionData():
	# @author: Florian Goebels
	# Class constructor
	# @Param:
	#		elutionProfileF flat file containing co elution data where each row represents a protein and collums contain spectral counts
	# @Class objects: 
	#		self.elutionMat numpy matric cnsisting of only the counts 
	#		self.normedElutionMat numpy matric containing normalized spectral counts 
	#		self.prot2Index dictonary mapping protein to it's respective row index in self.elutionMat and self.normedElutionMat
	def __init__(self, elutionProfileF, frac_count=2, max_frac_count=1):
		self.name = os.path.split(elutionProfileF)[-1]
		self.elutionMat, self.prot2Index  = self.loadElutionData(elutionProfileF, frac_count=frac_count, max_count_cutoff=max_frac_count)
		self.normedElutionMat = normalize_fracs(self.elutionMat)
		self.elutionMat = np.array(self.elutionMat)

	# @author: Florian Goebels
	# this methods load elution data stored as a flat file and returns the read in matrix and an index pointing each protein to it's rwo index
	# @Param:
	#		elutionProfileF elution data as flat file, tab separated
	def loadElutionData(self, elutionProfileF, frac_count = 2, max_count_cutoff=1):
		elutionProfileFH = open(elutionProfileF)
		elutionProfileFH.readline()
		i = 0
		elutionMat = []
		prot2Index = {}
		removed = 0
		pro_list = []
		for line in elutionProfileFH:
			line = line.rstrip()
			line = line.split("\t")
			protID = line[0]
			counts = map(float, line[1:])
			counts_array = np.array(counts)
			counts_array[np.isnan(counts_array)] = 0.0
			if len(list(set(np.where(np.array(counts_array) > 0.0)[0]))) >= frac_count and max(counts_array) >= max_count_cutoff:
				elutionMat.append(counts)
				prot2Index[protID] = i
				i += 1
				pro_list.append(protID)
			else:
				removed += 1
		print "finished processing %s\n removed %i (%.2f, total: %i, after filtering: %i) proteins found in less than %i fraction" % (elutionProfileF, removed, removed/(removed + len(prot2Index)), removed + len(prot2Index), len(prot2Index), frac_count)
		elutionProfileFH.close()
		elutionMat = np.nan_to_num(np.matrix(elutionMat))
		return elutionMat, prot2Index

	def orthmap(self, mapping):
		newMat = []
		newprot2Index = {}
		i = 0
		for prot, index in self.prot2Index.items():
			orth_prot = mapping.mapProtein(prot)
			if orth_prot == None: continue
			newMat.append(self.getElution(prot,normed=False))
			newprot2Index[orth_prot] = i
			i += 1
		self.elutionMat = np.array(newMat)
		self.prot2Index = newprot2Index
		self.normedElutionMat = normalize_fracs(self.elutionMat)
		

	# @author: Florian Goebels
	# return nurmalized co elution matrix, Normalization is based on previously published methods Pierre 2012 et al.
	def normalizeCoEulitionMat(self):
		self.elutionMat = normalize_fracs(self.elutionMat)
		

	# @author: Florian Goebels
	# returns elution profile of a given protein. If the protein has no profile the method returns NaN
	# @Param:
	#		normed if True return normalized counts, else raw counts
	def getElution(self, prot, normed=False):
		if prot not in self.prot2Index:
			return None
		if normed: 
			return self.normedElutionMat[self.prot2Index[prot]]
		else:
			return self.elutionMat[self.prot2Index[prot]]
	

	# @author: Florian Goebels
	# returns true if there is elution profile for a protein
	# @Param:
	#		prot  the protein in question
	def hasProt(self, prot):
		return prot in self.prot2Index

	# @author: Florian Goebels
	# reutrns a protein's row index in the co-elution table
	# @Param:
	#		prot  the protein in question
	def getProtIndex(self, prot):
		if prot in self.prot2Index:
			return self.prot2Index[prot]
		else:
			return None

	# @author: Florian Goebels
	# helper function for printing co-elution mat. Returns string
	# @Param:
	#		co elution mat as np matrix, which to be printed
	def printMat(self, mat):
		out = "ProtiID\tFraction_" + "\tFracion_".join(map(str, range(1,240)))
		for prot in self.prot2Index:
			index = self.prot2Index[prot]
			out += "\n%s\t%s" %  (prot, "\t".join(map(str,self.normedElutionMat[index])))
		return out

# @author: Florian Goebels
# helper function for normalization
def arr_norm(arr, axis=0):
	"""
	axis=0: normalize each column; 1: normalize each row
	"""
	mat = np.asmatrix(arr)
	return np.asarray(np.nan_to_num(mat / np.sum(mat, axis)))

def normalize_fracs(arr, norm_rows=True, norm_cols=True):
	if norm_cols:
		# Normalize columns first--seems correct for overall elution profile if
		# you're calculating correlation-type scores
		arr = arr_norm(arr, 0)
	if norm_rows:
		arr = arr_norm(arr, 1)
	return arr


# The following classes describe features for calculating Co elution scores
# each class needs to have two central functions: getScores and calculateScore
# whereas the first one allows the method to prepare the elution data how ever it needs and return int as two objects
# and the later one return the score for two given score objects
#Fueatures that can be calcuated in parrallel have always a static helper function, since pythons multiprocessing package does not support pickeling of objects

# @ author Florian Goebels
# returns apex similarity score which is based on previously published methods Pierre 2011 et al.
# returns 0 or 1 depending if the highest count peak of two proteins are in the same fraction or not
def calculateScore_apex(a,b):
	if a == b:
		return 1
	else:
		return 0

class Apex(object):
	def __init__(self):
		self.name="apex"
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = np.argmax(elutionData.getElution(prot))

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

	def get_name(self):
		return self.name

	calculateScore = staticmethod(calculateScore_apex)

# @ author Florian Goebels
# returns bayes correlation for two proteins
# reference on how the score is calculated is here: http://www.perkinslab.ca/sites/perkinslab.ca/files/Bayes_Corr.R
def calculateScoreBayes(a,b):
	global Bayes_cor_Mat
	index, a = a
	_, b = b
	return Bayes_cor_Mat[index][a][b]

class Bayes:
	def __init__(self, bayesType):
		self.name = "Bayes%i" % bayesType
		self.parallel = True
		self.bayesType = bayesType


	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		global current_bayes_cor, cor1, cor2, cor3, Bayes_cor_Mat
		if self.bayesType == 1:
			current_bayes_cor = cor1
		elif self.bayesType == 2:
			current_bayes_cor = cor2
		elif self.bayesType == 3:
			current_bayes_cor = cor3
		else:
			print("Invalid bayes selection")
			sys.exit()
			#TODO throw error instead of print
		elutionMat = elutionData.elutionMat
		dims = elutionMat.shape
		r_mat = robjects.r.matrix(robjects.FloatVector(elutionMat.flatten().tolist()), nrow=dims[0], ncol=dims[1])
		bayes_mat = np.array(current_bayes_cor(r_mat))
		Bayes_cor_Mat.append(bayes_mat)
		index = len(Bayes_cor_Mat) - 1
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = (index, elutionData.getProtIndex(prot))
			#prot2profile[fname][prot] = map( math.log, elutionData.getElution(prot)+1)


	def clear(self):
		global prot2profile, Bayes_cor_Mat
		Bayes_cor_Mat = []
		prot2profile = {}
		return True

	def get_name(self):
		return self.name

	# have this call always after __init__ since init initialize the right bayes function
	calculateScore = staticmethod(calculateScoreBayes)

# @ author Florian Goebels
# returns weighted cross correlation similarity score which is based on previously published methods Pierre 2011 et al.
def calculateScore_wcc(a,b):
	global r_wcc
	return r_wcc(robjects.FloatVector(a), robjects.FloatVector(b), 1)[0]

class Wcc:
	calculateScore = staticmethod(calculateScore_wcc)

	def __init__(self):
		self.name="wcc"
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = elutionData.getElution(prot)

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

	def get_name(self):
		return self.name


# @ author Florian Goebels
# returns travor correlation which is pearson correlation plus poisson noise to remove the influence of low counts
def calculateScore_PCCPN(a,b):
	global Poisson_cor_Mat
	index, a = a
	_, b = b
	return Poisson_cor_Mat[index][a][b]

def traver_corr(mat, repeat=1000, norm='columns', verbose=True):
	# As described in supplementary information in paper.
	# Randomly draw from poisson(C=A+1/M) for each cell
	# where A = the observed count and M is the total fractions
	# normalize each column to sum to 1
	# then correlate, and average together for repeat tries.
	def poisson_corr(mat, iteration_display, norm):
		if verbose: print(iteration_display)
		M = mat.shape[1]
		C = mat + 1/M
		poisson_mat = np.matrix(np.zeros(C.shape))
		for i in range(C.shape[0]):
			for j in range(M):
				poisson_mat[i,j] = np.random.poisson(C[i,j])
		if norm=='columns':
			poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
		elif norm=='rows': # seems to make no performance difference 1/25
			poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 1))
		corr = np.nan_to_num(np.corrcoef(poisson_mat))
		return corr
	avg_result = (reduce(operator.add, (poisson_corr(mat, i, norm=norm) for i in
										range(repeat))) / repeat)
	return avg_result

class Poisson:
	calculateScore = staticmethod(calculateScore_PCCPN)

	def __init__(self, repeat=100):
		self.name="poisson-%i" % (repeat)
		self.repeat=repeat
		self.parallel = True

	def init(self, elutionData):
		global Poisson_cor_Mat
		global prot2profile
		fname = "%s.%s" % (elutionData.name, self.name)
		noise_cor_mat = traver_corr(elutionData.elutionMat, self.repeat, 'columns', True)
		Poisson_cor_Mat.append(noise_cor_mat)
		index = len(Poisson_cor_Mat)-1
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = (index, elutionData.getProtIndex(prot))

	def setPoisson_cor_Mat(self, elutionData):
		global Poisson_cor_Mat
		Poisson_cor_Mat = traver_corr(elutionData.elutionMat, self.repeat, 'columns', True)

	def clear(self):
		global Poisson_cor_Mat
		Poisson_cor_Mat = ""
		global prot2profile
		prot2profile = {}

	def get_name(self):
		return self.name

# @ author Florian Goebels
# returns Mutual Information of two proteins
# mutual information is based on entropy calculation MI(x,y) = H(x) + H(y) - H(x,y) 
# wehre H means entropy
def calculateScore_MI(a, b):
	entropy_a, a_upper, a_lower, _ = a
	entropy_b, b_upper, b_lower, numFracs = b
	joint_probs = np.array(map(lambda x: len(x)/numFracs, [a_upper&b_upper, a_upper&b_lower, a_lower&b_upper, a_lower&b_lower]))
	joint_entropy_a_b = entropy(joint_probs, 2)
	mutual_information =  entropy_a  + entropy_b - joint_entropy_a_b
	return mutual_information

def bin_entropy(p):
	return entropy(np.array([p, 1 - p]))


def entropy(probs, base=0):
	if base == 0: base = len(probs)
	tmp_probs = probs
	tmp_probs[tmp_probs == 0] = 1
	return -sum(probs * map(lambda x: math.log(x, base), tmp_probs))


def getFracs(a, cutoff):
	upper = set([i for i, v in enumerate(a) if v > cutoff])
	lower = set([i for i, v in enumerate(a) if v <= cutoff])
	return (upper, lower)

class MutualInformation():
	calculateScore = staticmethod(calculateScore_MI)

	def __init__(self, minCounts = 2):
		self.name="MI"
		self.minCounts = minCounts
		self.parallel = True
	
	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			profile = elutionData.getElution(prot)
			(profile_upper, profile_lower) = getFracs(profile, self.minCounts)
			numFracs = len(profile)
			prot_entropy = bin_entropy(len(profile_upper) / numFracs)
			prot2profile[fname][prot] = (prot_entropy, profile_upper, profile_lower, numFracs)

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

	def get_name(self):
		return self.name



# @ author Florian Goebels
# calculates Jaccard overlap score which is as follows
# Jaccard(x,y) = sum(x!=1 and y!=1)/(sum(x!=1) + sum(y!=1))
def calculateScore_Jaccard(a_non_zero_fracs,b_non_zero_fracs):
	a_and_b_non_zero_fracs = len(a_non_zero_fracs & b_non_zero_fracs)
	a_or_b_non_zero_fracs = len(a_non_zero_fracs | b_non_zero_fracs)
	if a_or_b_non_zero_fracs == 0:
		return 0
	else:
		return a_and_b_non_zero_fracs/a_or_b_non_zero_fracs

class Jaccard():
	calculateScore = staticmethod(calculateScore_Jaccard)

	def __init__(self):
		self.name="Jaccard"
		self.non_zero_fracs_for_prot = {}
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = set(np.nonzero(elutionData.getElution(prot))[0])

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

	def get_name(self):
		return self.name



# @ author Florian Goebels
# return Pearson correlation of two proteins
def calculateScore_Pearson(a,b):
	score = scipy.stats.pearsonr(a, b)[0]
	if math.isnan(score): return 0.0
	return scipy.stats.pearsonr(a, b)[0]

class Pearson:
	calculateScore = staticmethod(calculateScore_Pearson)

	def __init__(self):
		self.name = "Pearson"
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = elutionData.getElution(prot)

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

	def get_name(self):
		return self.name
	
# @ author Lucas Ming Hu
# This is a helper class for getting GeneMANIA functional evidence for a given ElutionData object
class Genemania:
	
	def __init__(self, taxID):
		self.taxoID = taxID
		#create a protein_pair_mapping dictionary, GeneManiaName - UniProtName
		self.nameMappingDict = {}
		#create a geneName mapping dictionary based on Uniprot website database
		self.map_proteinNames()
		# Get all Genemania files
		self.files = []
		self.catchFile()
		# all functional evidence codes in GeneMANIA, excluding "Physical" and "complexes" and "Predicted" to eliminate circularity
		self.functionalEvidences = ['Co-expression', 'Genetic', 'Other', 'Shared']
		# ScoreCalc object contains edges and it's associated GeneMANIA scores
		self.scoreCalc = CalculateCoElutionScores("", "", "", num_cores=1, cutoff=0)
		# loads all of Worm Gene
		self.load_genemania()

	# @auothor Lucas Ming Hu
	# the catchFile function can help to download files from GeneMANIA website automatically.	
	def catchFile(self): 
		taxoIDspeciesDic = {'3702':'Arabidopsis_thaliana', '6239':'Caenorhabditis_elegans', '7955':'Danio_rerio', 
		                    '7227':'Drosophila_melanogaster','562':'Escherichia_coli','9606':'Homo_sapiens',
		                    '10090':'Mus_musculus','10116':'Rattus_norvegicus','4932':'Saccharomyces_cerevisiae'} 
		
		if self.taxoID not in taxoIDspeciesDic:
			return None #TODO throw illegal argument exception
	
		#urlbase = 'http://genemania.org/data/current'
		urlbase = 'http://genemania.org/data/archive/2014-10-15/'
		speciesURL = os.path.join(urlbase, taxoIDspeciesDic[self.taxoID])
		r = urllib.urlopen(speciesURL).read()
#		soup = BeautifulSoup(r)

		soup = BeautifulSoup(r, "html.parser")
    
		table = soup.find('table')
    
		allcell = []
		for row in table.find_all('tr'):
			for col in row.find_all('td'):
				allcell.append(col.getText())

		#filtering 
		#self.files = []
		for c in allcell:
			if '.txt' in c:
				self.files.append(os.path.join(speciesURL,c))

	# @author: Lucas Ming Hu        
	# a helper function to get the average of the GeneMANIA scores 
	# for each line of evidence
	def average(self, secondaryEvidenceDic):
		resultDict = defaultdict(float)		
		for key in secondaryEvidenceDic:
			resultDict[key] = sum(secondaryEvidenceDic[key]) * 1.0 / len(secondaryEvidenceDic[key])
		return resultDict

	# @author: Lucas Ming Hu
	# returns Functional anotation scores as a CalculateCoElutionScores Object
	def load_genemania(self):

#		data holders
		scores = {}
		self.header = []
		# read online database - by species taxoID

		for i, f_evidence in enumerate(self.functionalEvidences):
			this_evidence_scores = {}
			self.scoreCalc.header.append("GeneMania_%s" % f_evidence)
			for fp in self.files:                        #for de-bugging, I only used the first three files
				filename = str(fp.split('/')[-1])
				if filename.startswith(f_evidence):
					print "Processing: %s" % (filename)
					fh = urllib2.urlopen(fp)
					fh.readline()
					for line in fh:
						proteinA, proteinB, score = line.split()

						# transfer the GeneMANIA gene name to its corresponding UniPort name
						if ((proteinA in self.nameMappingDict) and (proteinB in self.nameMappingDict)):
							proteinAUniprot = self.nameMappingDict[proteinA]
							proteinBUniprot = self.nameMappingDict[proteinB]

							edge = "\t".join(sorted([proteinAUniprot, proteinBUniprot]))
							score = float(score)

							if edge not in this_evidence_scores:
								this_evidence_scores[edge] = [0, 0]
								(this_evidence_scores[edge])[0] = score
								(this_evidence_scores[edge])[1] = 1
							else:
								(this_evidence_scores[edge])[0] = (this_evidence_scores[edge])[0] + score
								(this_evidence_scores[edge])[1] = (this_evidence_scores[edge])[1] + 1

					fh.close()

			for edge in this_evidence_scores:
				score, counts = this_evidence_scores[edge]
				avg_score = score/counts
				if edge not in scores: scores[edge] = [0]*len(self.functionalEvidences)
				scores[edge][i] = avg_score

			self.scoreCalc.scores = np.zeros((len(scores.keys()), len(self.functionalEvidences)))
			i = 0
			for edge in scores:
				self.scoreCalc.ppiToIndex[edge] = i
				self.scoreCalc.IndexToPpi[i] = edge
				self.scoreCalc.scores[i,:] = scores[edge]
				i += 1

	def getScoreCalc(self):
		return self.scoreCalc


	# @author: Lucas Ming Hu
	# read name mapping database and put Uniprot and corresponding GeneMANIA name into a dictionary
	# key: GeneMANIA_name value: Uniprot_name
	def map_proteinNames(self):
		
		taxoIDurl = {'6239':'http://www.uniprot.org/uniprot/?query=taxonomy%3A6239&sort=score&columns=id,genes(ORF)&format=tab',
					 '3702':'http://www.uniprot.org/uniprot/?query=taxonomy%3A3702&sort=score&columns=id,genes(OLN)&format=tab',
					 '7955':'http://www.uniprot.org/uniprot/?query=taxonomy%3A7955&sort=score&columns=id,database(Ensembl)&format=tab',
					 '7227':'http://www.uniprot.org/uniprot/?query=taxonomy%3A7227&sort=score&columns=id,database(FlyBase)&format=tab',
					 '4932':'http://www.uniprot.org/uniprot/?query=taxonomy%3A4932&sort=score&columns=id,genes&format=tab'}

		response = urllib2.urlopen(taxoIDurl[self.taxoID])
		
		html = response.readlines() #return everything in a list, each item in the list is a line in original file
		
		unipront_geneNames_dic = {}
		
		for item in html[1:]:
			
			items = item.split("\t")
			uniprot = items[0]

			geneNames = items[1].strip(';') #some names has ; at the end
			geneNames = items[1].strip()
			
			genes_list = re.split('[;\s\|\/]', geneNames)
			
			new_list = list(genes_list) #make a new list to clone original list, otherwise, the code wont work.
			
			for gene_names in genes_list:
			    
				if "CELE" in gene_names:
					new_list.remove(gene_names)
			
			new_list = filter(None, new_list) #remove the empty item from the list
				
			if len(new_list) > 0 :
				unipront_geneNames_dic[uniprot] = new_list
			
		for key, value in unipront_geneNames_dic.iteritems():
			for i in range(0, len(value)):
				if value[i] not in self.nameMappingDict:
					self.nameMappingDict[value[i]] = key




# @ author Florian Goebels
# returns Euclidean distance of two proteins
def calculateScore_euclidean(a,b):
	return 1-distance.euclidean(a,b)

class Euclidiean:
	calculateScore = staticmethod(calculateScore_euclidean)

	def __init__(self):
		self.name = "Euclidiean"
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = elutionData.getElution(prot, normed=True)

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

	def get_name(self):
		return self.name

# @ author Florian Goebels
# This is a helper class for calculating co elution scores for a given ElutionData object
class CalculateCoElutionScores():
	# @author: Florian Goebels
	# this method inits the object
	# @Param:
	#		elutionData (optional) specify which data needs to be processed
	def __init__(self, scores, elutionData, outF, num_cores, cutoff = 0.5, verbose=True):
		self.num_cores = num_cores
		self.outF = outF
		self.verbose = verbose
		self.elutionData = elutionData
		self.scores = np.array([])
		self.header = ["ProtA","ProtB"]
		self.ppiToIndex = {}
		self.IndexToPpi = {}
		self.cutoff = cutoff
		self.features = scores
		self.to_predict = 0
		self.scoreF = ""
		self.fun_anno=""

		for eD in self.elutionData:
			for score in self.features:
				self.header.append("%s.%s" % (eD.name, score.name))
	###
	def get_val_proteins(self):
		val_proteins = set()

		for items in self.ppiToIndex.keys():
			proteins = items.split("\t")
			val_proteins.add(proteins[0])
			val_proteins.add(proteins[1])

		return val_proteins
	###

	def getShape(self):
		return self.scores.shape

	def get_all_scores(self):
		return self.scores

	def add_fun_anno(self, fun_anno):
		self.fun_anno = fun_anno
		self.merge(fun_anno, "l")

	def open(self):
		self.i = 0
		if self.scoreF != "":
			self.scoreFH = open(self.scoreF)
			self.scoreFH.readline() # remove header line
		else:
			self.to_predict = self.scores.shape[0]

	def has_edge(self, edge):
		return edge in self.ppiToIndex

	def get_score(self, edge):
		if not self.has_edge(edge): return None
		return self.scores[self.ppiToIndex[edge],:]

	def get_next(self):
		edge = ""
		scores = []
		# check if edge to predicted is in memory or saved on HD
		if self.scoreF != "":
			line = self.scoreFH.readline()
			if line == "": return None, None
			line = line.rstrip()
			line = line.split("\t")
			edge = "\t".join(sorted(line[0:2]))
			scores = map(float, line[2:])
			# add functional annotation if given
			if self.fun_anno!="":
				to_add = [0]*(len(self.fun_anno.header)-2)
				if self.fun_anno.has_edge(edge):
					to_add = self.fun_anno.get_score(edge)
				scores = np.append(scores, to_add)
		# Retrive edge and edge scores from local memory
		else:
			edge = self.IndexToPpi[self.i]
			scores = self.scores[self.i,:]
			self.i+=1

		return edge, np.array(scores)

	def close(self):
		if self.scoreF != "":
			self.scoreFH.close()

	# @author: Florian Goebels
	# this method combines who CalculateCoElutionScores objects onto one by comping the toMerge object into the self object
	# @Param:
	#		CalculateCoElutionScores toMerge a second CalculateCoElutionScores which should be combined with self object
	#		mode donates how to merge the sets, left (l), right (r), union (u), or  interaction(i)
	def merge(self, toMerge, mode):
		allPPIs = ""
		if mode == "u":
			allPPIs = set(self.ppiToIndex.keys()) | set(toMerge.ppiToIndex.keys())

		if mode == "l":
			allPPIs = set(self.ppiToIndex.keys())
		if mode == "r":
			allPPIs = set(toMerge.ppiToIndex.keys())
		if mode == "i":
			allPPIs = set(self.ppiToIndex.keys()) & set(toMerge.ppiToIndex.keys())
		#print len(set(toMerge.ppiToIndex.keys()))
		#print len(allPPIs)

		numFeature_in_merge = len(toMerge.header)-2

		numFeature_in_self = len(self.header)-2
		new_scores = np.zeros((len(allPPIs), numFeature_in_merge + numFeature_in_self))

		new_ppiToIndex = {}
		new_IndexToPpi = {}
		k = 0
		for ppi in allPPIs:
			scoresA = [0]*numFeature_in_self
			scoresB = [0]*numFeature_in_merge
			if ppi in self.ppiToIndex:
				scoresA = self.scores[self.ppiToIndex[ppi],:]
			if ppi in toMerge.ppiToIndex:
				scoresB = toMerge.scores[toMerge.ppiToIndex[ppi],:]
			new_score = np.array(list(scoresA) + list(scoresB))
			new_scores[k,:] = new_score
			new_IndexToPpi[k] = ppi
			new_ppiToIndex[ppi] = k
			k+=1
		self.scores = new_scores
		self.ppiToIndex = new_ppiToIndex
		self.IndexToPpi = new_IndexToPpi
		self.header.extend(toMerge.header[2:])

	# @author: Florian Goebels
	# create co elution table for a given list of scores and ppis
	# @Param:	
	#		scoreTypes	a list of co elutoin score objects
	#		PPIs		a list of ppis for which all scores in scoreTypes schoukd be calculated
	def calculateScores(self, toPred, tokeep=set([])):
		task_queue = mp.JoinableQueue()
		out_queue = mp.Queue()
		self.scores = []
		self.ppiToIndex = {}
		num_rows = len(tokeep)
		num_features = len(self.features)*len(self.elutionData)
		self.scores = np.zeros((num_rows, num_features))
		for i in range(self.num_cores):  # remove one core since listener is a sapareted process
			mp.Process(target=getScore, args=(task_queue, out_queue)).start()
		outFH = open(self.outF, "w")
		print >> outFH, "\t".join(self.header)
		k = 0
		ppi_index = 0
		write_buffer = ""
		for ppi in toPred:
			k += 1
			if k % 100000 == 0:
				if self.verbose: print(k)
				outFH.write(write_buffer)
				outFH.flush()
				write_buffer = ""
			i = 0
			for _ in self.elutionData:
				for score in self.features:
					fname = self.header[i+2]
					task = (ppi, i, fname, score.calculateScore)
					task_queue.put(task)
					i += 1
			task_queue.join()
			ppi_scores = [0]*num_features
			for _ in range(num_features):
				score_index, score = out_queue.get()
				ppi_scores[score_index] = score
			ppi_scores = np.nan_to_num(np.array(ppi_scores))

			if len(list(set(np.where(ppi_scores > self.cutoff)[0]))) > 0:
				self.to_predict += 1
				write_buffer +=  "%s\t%s\n" % (ppi, "\t".join(map(str, ppi_scores)))
				if ppi in tokeep:
					self.ppiToIndex[ppi] = ppi_index
					self.IndexToPpi[ppi_index] = ppi
					self.scores[ppi_index,:] = ppi_scores
					ppi_index += 1

		print >> outFH, write_buffer
		outFH.close()
		print("done calcualting co-elution scores")
		self.scores = self.scores[0:ppi_index,:]
		for i in range(self.num_cores):
			task_queue.put('STOP')
		self.scoreF = self.outF


	# @ author Florian Goebels
	# A filter for removing all possible protein pairs when predicting the network form elution data
	# here we decided to remove all candidate ppis with an Jaccard score of 0, since those interaction do not co-elute in any observed fraction
	def filter_interactions(self, eData, ppis):
		print("filtering")
		out = set([])
		del_Jaccard_scores = False
		global prot2profile
		if eData.name + ".Jaccard" not in prot2profile:
			this_Jaccard = Jaccard()
			this_Jaccard.init(eData)
			del_Jaccard_scores = True
		for ppi in ppis:
			protA, protB = ppi.split("\t")
			profileA = prot2profile[eData.name + ".Jaccard"][protA]
			profileB = prot2profile[eData.name + ".Jaccard"][protB]
			if not profileA.isdisjoint(profileB): out.add(ppi)

		if del_Jaccard_scores:
			del prot2profile[eData.name + ".Jaccard"]
		return out


	def getAllPairs(self):
		allfilteredPPIs = set([])
		for elution_data in self.elutionData:
			print("Filtering: %s" % (elution_data.name))
			allprots = elution_data.prot2Index.keys()
			candidatePPis = set([])
			for i in range(len(allprots)):
				for j in range(i + 1, len(allprots)):
					protA = allprots[i]
					protB = allprots[j]
					if protA == protB: continue
					protA, protB = sorted([protA, protB])
					candidatePPis.add("%s\t%s" % (protA, protB))
			print("Before filtering %i PPIs" % (len(candidatePPis)))
			filteredPPIs = self.filter_interactions(elution_data, candidatePPis)
			del candidatePPis
			print("After filtering %i PPIs" % (len(filteredPPIs)))
			allfilteredPPIs |= filteredPPIs
		print("Num of PPIs across all data sets after filtering %i" % (len(allfilteredPPIs)))
		return allfilteredPPIs

	def calculate_coelutionDatas(self, gs=""):
		toPred = self.getAllPairs()
		if gs !="":
			self.calculateScores(toPred, gs.get_edges())
		else:
			self.calculateScores(toPred)

	# @author: Florian Goebels
	# prints table
	# @Param:
	#		labels print class lable or not
	def toTable(self):
		valid_ppis = set(self.ppiToIndex.keys())
		out = ""
		out = "\t".join(self.header)
		for i in range(self.scores.shape[0]):
			ppi = self.IndexToPpi[i]
			if ppi not in valid_ppis: continue
			scores = "\t".join(map(str, self.scores[i, :]))
			line = "%s\t%s" % (ppi, scores)
			out += "\n" + line
		return out

	def readTable(self, scoreF, gs=""):
		self.scoreF = scoreF
		scoreFH = open(scoreF)
		tokeep = set([])
		self.header = scoreFH.readline().rstrip().split("\t")
		row_num = 0
		if gs !="":
			tokeep = gs.positive | gs.negative
			row_num = len(tokeep)
		else:
			row_num = lineCount(scoreF)
		self.scores = np.zeros((row_num, len(self.header)-2))
		i = 0
		self.ppiToIndex = {}
		for line in scoreFH:
			self.to_predict += 1
			line = line.rstrip()
			if line == "": continue
			linesplit = line.split("\t")
			edge = "\t".join(sorted(linesplit[0:2]))
			if gs != "" and edge not in tokeep: continue
			edge_scores = np.nan_to_num(np.array(map(float, linesplit[2:])))
			self.scores[i,:] = edge_scores
			self.IndexToPpi[i] = edge
			self.ppiToIndex[edge] = i
			i += 1
		scoreFH.close()
		self.scores = self.scores[0:i, :]

	# @author: Florian Goebels
	# return stored co elution scores for a given data set in form that is required for sklearn to learn and predict interaction
	def toSklearnData(self, gs):
		ids = []
		targets = []
		used_indeces = []
		for i in range(self.scores.shape[0]):
			ppi = self.IndexToPpi[i]
			label = "?"
			if ppi in gs.positive: label = 1
			if ppi in gs.negative: label = 0
			if label != "?":
				targets.append(label)
				ids.append(ppi)
				used_indeces.append(i)
		return ids, self.scores[used_indeces, :], np.array(targets)

# @ author Florian Goebels
# wrapper for machine learning
class CLF_Wrapper:
	# @author: Florian Goebels
	# class initializer, supports both random forest and svm
	# @Param:
	#		data matrix with features where each row is a data point
	#		targets list with class lables
	# 		forest if true ml is random forst, if false ml is svm
	def __init__(self, num_cores, forest=False, useFeatureSelection= False):
		self.num_cores = num_cores
		if forest:
			print("using Random forest")
			thisCLF = RandomForestClassifier(n_estimators=1000, n_jobs=self.num_cores, random_state=0)
		else:
			print("Using SVM")
			thisCLF =  svm.SVC(kernel="linear", probability=True)
		if useFeatureSelection:
			print("Using Feature selection RFECV")
			self.clf = Pipeline([
				('feature_selection', RFECV(estimator=thisCLF, step=1, scoring="accuracy")),
				('classification', thisCLF)
			])
		self.clf = thisCLF

	def fit(self, data, targets):
		self.clf.fit(data, targets)

	def eval(self, data, targets):
		probs = self.predict_proba(data)
		preds = self.predict(data)
		return self.get_metrics(probs, preds, targets)

	def get_metrics(self, probs, preds, targets):
		precision = metrics.precision_score(targets, preds, average=None)[1]
		recall = metrics.recall_score(targets, preds, average=None)[1]
		fmeasure = metrics.f1_score(targets, preds, average=None)[1]
		auc_pr = average_precision_score(targets, preds)
		auc_roc = roc_auc_score(targets, preds)
		curve_pr = precision_recall_curve(targets, probs)
		curve_roc = roc_curve(targets, probs)
		return [precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc]

	def cv_eval(self, data, targets, folds= 5):
		skf = StratifiedKFold(folds)
		probs = []
		preds = []
		this_targets = []
		i = 1
		for train, test in skf.split(data, targets):
			#print "Processing fold %i" % i
			print "Processing fold data..."
			i += 1
			self.fit(data[train], targets[train])
			probs.extend(self.predict_proba(data[test]))
			preds.extend(self.predict(data[test]))
			this_targets.extend(targets[test])
		return self.get_metrics(probs, preds, this_targets)


	# @author: Florian Goebels
	# @Param:
	#		toPred matric where each row is a data point and predicts interaction propability for a given set
	#		note trainer needs to be already trained to be able to predict
	def predict_proba(self, toPred):
		probas = self.clf.predict_proba(toPred)
		return probas[:,1]

	def predict(self, toPred):
		preds = self.clf.predict(toPred)
		return preds
"""
class MLP_wrapper(object):


	def __init__(self):
		print "Using MLP with Keras/tensorflow"
		self.model = Sequential()

	def fit(self, data, labels):
		num_features = data.shape[1]
		self.model.add(Dense(512, input_dim=num_features, activation='relu'))
		self.model.add(Dropout(0.5))
		self.model.add(Dense(512, activation='relu'))
		self.model.add(Dropout(0.5))
		self.model.add(Dense(1, activation='sigmoid'))
		self.model.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])
		self.model.fit(data, labels, epochs=500, batch_size=128)

	def eval(self, data, labels):
		return self.model.evaluate(data, labels)

	def predict_proba(self, toPred):
		return [x[0] for x in self.model.predict(toPred)]

	def predict(self, toPred):
		return [round(x[0]) for x in self.model.predict(toPred)]

class SAE_wrapper(MLP_wrapper):

	def __init__(self):
		print "Using stacked autoencoder"


	def fit(self, data, labels):
		print data.shape
		print len(labels)
		num_features = data.shape[1]
		input = Input(shape=(num_features,))

		hidden_layer1 = 1250
		hidden_layer2 = 600
		hidden_layer3 = 100
		hidden_layer4 = 1

		encoded = Dense(hidden_layer1, activation='relu')(input)
		encoded = Dense(hidden_layer2, activation='relu')(encoded)
		encoded = Dense(hidden_layer3, activation='relu', activity_regularizer=regularizers.l1(10e-5))(encoded)
		encoded = Dense(hidden_layer4, activation='sigmoid')(encoded)

		self.model = Model(input=input, output=encoded)

		self.model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=['accuracy'])

		self.model.fit(data, labels, epochs=20, batch_size=500, shuffle=True)
"""




# @ author: Lucas Ming Hu
# This is a helper function can help users to read their personal evidence file as extra functional evidence and integrate into the pipeline.
class ExternalEvidence:

	# author: Lucas Ming Hu
	# This methods initiates the object
	# FilenameWithDic is the directory with filename for the etxernal functional evidence data
	def __init__(self, FilenameWithDic):
		self.FilenameWithDic = FilenameWithDic
		self.scoreCalc = CalculateCoElutionScores("", "", "", 1)
		self.scoreCalc.readTable(FilenameWithDic) # added by flo using  CalculateCoElutionScores funtion for reading in functional scores
		self.scores = {}

	"""
	# @ author: Lucas Ming Hu
	# read external functional evidence, user can supply the external functional evidence as they want to integrate into the experimental data
	def readFile(self):

		#first, read the external functional evidence data into a dictionary
		#with open(self.FilenameWithDic) as fp:
		#	header = fp.readline()
		#	header = header.rstrip()
		#	print header
		#	for evidence in header.split("\t")[2:]:
		#		self.scoreCalc.header.append(evidence)
		#	for line in fp:
		#		names = line.rstrip().split("\t")
		#		if names[0] == "protein1":
		#			for evidence in names[2:]:
		#				self.scoreCalc.header.append(evidence)
		#		else:
		#			edge = "\t".join(sorted([names[0], names[1]]))
		#			self.scores[edge] = names[2:]
		###the code above always give errors, I changed back...

		with open(self.FilenameWithDic) as fp:
			for line in fp:
				names = line.rstrip().split("\t")
				if names[0] == "protein1" or names[0] == "ProtA":
					for evidence in names[2:]:
						self.scoreCalc.header.append(evidence)
				else:
					edge = "\t".join(sorted([names[0], names[1]]))
					self.scores[edge] = names[2:]

		i = 0
		self.scoreCalc.scores = np.zeros((len(self.scores.keys()), len(names[2:])))

		#second, from the dictionary, read the file into the
		for edge in self.scores:
			self.scoreCalc.ppiToIndex[edge] = i
			self.scoreCalc.IndexToPpi[i] = edge
			self.scoreCalc.scores[i, :] = self.scores[edge]
			i += 1
	"""

	def getScoreCalc(self):
		return self.scoreCalc

# @ author Lucas Ming Hu
# This is a helper class for getting STRING functional evidence for a given ElutionData object
# Make sure to leave 300MB in your disk for storing files from STRING, the script needs to download it first, since
# API in STRING doesn't support bulk data manipualation.
class STRING:
	def __init__(self, taxID, datadir = ""):
		self.datadir = datadir
		# the input is the Taxo ID for the given species
		self.TaxID = taxID
		# create a protein_pair_mapping dictionary, STRINGID - UniProtName
		self.nameMappingDict = {}
		# create a geneName mapping dictionary based on Uniprot website database
		self.nameMapping()
		# ScoreCalc object contains edges and it's associated STRING scores
		self.scoreCalc = CalculateCoElutionScores("", "", "", 1)
		# loads all of Worm Gene
		self.load_string()


	# @auothor Lucas Ming Hu
	# the nameMapping function can help to download name mapping file from STRING website automatically.
	def nameMapping(self):

		#this is the url for STRING_id and Uniprot_id mapping file
		url = "http://string-db.org/download/protein.aliases.v10.5/" + str(self.TaxID) + ".protein.aliases.v10.5.txt.gz"
		filename_protein = url.split("/")[-1]
		if self.datadir != "": filename_protein = self.datadir + os.sep + filename_protein

		with open(filename_protein, "wb") as f:
			r = requests.get(url)
			f.write(r.content)
		f.close()

		self.nameMappingDict = {}

		# read the local .gz file and store STRING_id and its uniprot_id to a dictionary
		with gzip.open(filename_protein, 'r') as fin:
			for line in fin:
				line.rstrip()
				items = line.split()
				# only read the uniprot ID into the dictionary...
				if len(items) >= 3 and "BLAST_UniProt_AC Ensembl_UniProt_AC" in line:
					self.nameMappingDict[items[0]] = items[1]  # key is STRING protein ID and value is Uniprot ID
		fin.close()

	# @auothor Lucas Ming Hu
	# the load_string function can help to download functional evidence data from STRING website automatically.
	# we exclude evidences from "exp", "database" and "combined"
	# the evidences left are: "neighbour", "fusion", "co-occurence", "co-expression" and "textmining"
	def load_string(self):

		# download the interaction data file from internet as compressed .gz file...
		url2 = "http://string-db.org/download/protein.links.detailed.v10.5/" + str(
			self.TaxID) + ".protein.links.detailed.v10.5.txt.gz"
		filename_interaction = url2.split("/")[-1]
		if self.datadir != "": filename_interaction = self.datadir + os.sep + filename_interaction

		with open(filename_interaction, "wb") as f:
			r = requests.get(url2)
			f.write(r.content)
		f.close()

		temp_score_dict = {}  # select the socres from the evidence we want and put them into a dictonary.
		self.header = [] # self.scoreCalc needs header row
		print ("finish reading the protein mapping file from STRING")

		# then read the .gz file
		with gzip.open(filename_interaction, 'r') as fin:
			for line in fin:
				line.rstrip()
				items = line.split()
				if "protein1 protein2" in line:

					# if changing the evidence, here is also changing
					for evidences in items[2:6]:
						self.scoreCalc.header.append(evidences)

					#self.scoreCalc.header.append(items[8]) #if including textmining, include this line.

				else:

					if items[0] in self.nameMappingDict and items[1] in self.nameMappingDict:
						proteinA = self.nameMappingDict[items[0]]
						proteinB = self.nameMappingDict[items[1]]
						edge = "\t".join(sorted([proteinA, proteinB]))

						score = items[2:]

						#select the slected evidence only, we exclude "exp", "database" and "combined_score"
						#usefulScores = [None] * 5
						#usefulScores[0:4] = score[0:4]
						#usefulScores[4] = score[6]

						#select the slected evidence only, we exclude "exp", "database" , "textmining" and "combined_score"
						usefulScores = [None] * 4
						usefulScores[0:4] = score[0:4]

						temp_score_dict[edge] = usefulScores
		print ("finish reading protein-protein functional evidence data from STRING")
		fin.close()

		self.scoreCalc.scores = np.zeros((len(temp_score_dict.keys()), len(usefulScores)))

		i = 0
		for edge in temp_score_dict:
			self.scoreCalc.ppiToIndex[edge] = i
			self.scoreCalc.IndexToPpi[i] = edge
			self.scoreCalc.scores[i,:] = temp_score_dict[edge]
			i += 1

	# @author: Lucas Ming Hu
	# the getScoreCalc function will return the self.scoreCalc .
	def getScoreCalc(self):
		return self.scoreCalc

