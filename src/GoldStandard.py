#!/usr/bin/python
from __future__ import division
import numpy as np
import urllib2
from xml.dom import minidom
import requests, sys
import os
import copy
import re
import math
import sys
import random as rnd
from sklearn.cluster import KMeans
rnd.seed(1)
np.random.seed(1)

# CORUM database got changed, so need to use new way to read the original database...
from StringIO import StringIO
from zipfile import ZipFile
from urllib import urlopen

class Goldstandard_from_Complexes():

	def __init__(self, name="unnamed", ratio = 5):
		self.complexes = Clusters(False)
		self.name = name
		self.ratio = ratio
		self.positive, self.negative = set([]), set([])

	# get all proteins from Goldstandard_from_Complexes object
	# a trial version added by Lucas HU
	def get_proteins(self):
		return self.complexes.get_all_prots()

	def set_lb(self, lb):
		self.complexes.lb = lb

	def set_ub(self, ub):
		self.complexes.ub = ub

	def make_reference_data(self, db_clusters, orthmap="", found_prots=""):
		total_complexes = 0
		i = 0
		for db_clust in db_clusters:
			tmp_clust = copy.deepcopy(db_clust.get_complexes())
			total_complexes += len(tmp_clust.complexes.keys())
			if tmp_clust.need_to_be_mapped == True and orthmap !="":
				orthmap.mapComplexes(tmp_clust)
			for compl in tmp_clust.complexes:
				self.complexes.addComplex("%i;%s;%s" % (i, db_clust.name, compl), tmp_clust.complexes[compl])
				i += 1

		print "Total number of complexes %i in %s" % (total_complexes, self.name)
		print "Number of complexes after ortholog mapping %i complexes in %s" % (len(self.complexes.complexes), self.name)

		if found_prots != "":
			self.complexes.remove_proteins(found_prots)
			print "After removing not indetified proteins %i number of complexes in % s" % (len(self.complexes.complexes), self.name)

		self.complexes.filter_complexes()
		print "After size filtering %i number of complexes in % s" % (len(self.complexes.complexes), self.name)

		self.complexes.merge_complexes()
		self.complexes.filter_complexes()
		print "After mergning %i number of complexes in % s" % (len(self.complexes.complexes), self.name)

		self.make_pos_neg_ppis()

	def make_pos_neg_ppis(self, val_ppis=""):
		self.positive, self.negative = self.complexes.getPositiveAndNegativeInteractions() #this is not right, PPI in complexes in hold_out might be considered as negative in trianing...
		if val_ppis!="":
			self.positive &= val_ppis
			self.negative &= val_ppis

	# a function added by Lucas HU to add PPIs into the self.positive and self.negative, inseatd of using the functions above
	# positive is also set object
	def add_positive(self, positives):
		self.positive = self.positive | positives

	def add_negative(self, negatives):
		self.negative = self.negative | negatives

	def get_complexes(self):
		return self.complexes

	def get_complex(self, complex_id):
		if complex_id in self.complexes:
			return self.complexes[complex_id]
		else:
			return ""

	def get_goldstandard(self):
		return self.positive, self.negative

	def get_positive(self):
		return self.positive

	def get_negative(self):
		return self.negative

	def get_edges(self):
		return self.positive | self.negative

	def get_overlap_scores(self, setA, setB):
		overlap_set = setA & setB
		overlap_score = (len(overlap_set) ** 2) / (len(setA) * len(setB))
		return overlap_score

	def get_overlap_proteins(self, setA, setB):
		overlap_set = setA & setB
		return overlap_set

	# a Funtion written by Lucas Hu,
	# to minimize the overlapped of the splitted complexes set.
	# only support num_folds = 2
	def n_fols_split(self, num_folds, overlap="False"):

		ref_cluster_ids = self.complexes.complexes.keys()
		out_folds = []

		# create a n by n zero matrix to store the binary overlap scores.
		overlap_score_matrix = np.zeros((len(self.complexes.complexes), len(self.complexes.complexes)))
		for i in range(len(self.complexes.complexes)):
			for j in range(i, len(self.complexes.complexes)):
				if i == j:
					overlap_score_matrix[i, j] = 1
				else:
					overlapScore = self.get_overlap_scores(self.complexes.complexes[i], self.complexes.complexes[j])
					overlap_score_matrix[i, j] = overlapScore
					overlap_score_matrix[j, i] = overlapScore

		estimator = KMeans(n_clusters=2) # A k-mean classifier to divide the whole set into two sets

		estimator.fit(overlap_score_matrix)
		labels = estimator.labels_

		# get the index of the item in the np.array...for two groups (group 0 and group 1)...
		itemindex_zero = (np.where(labels == 0)[0]).tolist()
		itemindex_one = (np.where(labels == 1)[0]).tolist()

		# do the while loop until the two sets are balanced ...
		while len(itemindex_one) > len(itemindex_zero):

			for individual_complex_one_index in itemindex_one:
				temp_protein_set_one = self.complexes.complexes[individual_complex_one_index]
				master_overlapped_set = set()
				max_overlap = 0
				max_index = ''

				for individual_complex_zero_index in itemindex_zero:
					temp_protein_set_zero = self.complexes.complexes[individual_complex_zero_index]
					overlapped_set = self.get_overlap_proteins(temp_protein_set_one, temp_protein_set_zero)
					master_overlapped_set = master_overlapped_set | overlapped_set

				if max_index == '':
					max_index = individual_complex_one_index
					max_overlap = len(master_overlapped_set)
				else:
					if max_overlap < len(master_overlapped_set):
						max_index = individual_complex_one_index
						max_overlap = len(master_overlapped_set)

			itemindex_one.remove(max_index)
			itemindex_zero.append(max_index)

		# further remove overlapped complexes...
		# identify the most overlapped complex from set1, remove...
		# idenifty the most overlapped complexes from set2, remove...
		# repeat for n rounds...
		round = 0
		while round < 10:

			for individual_complex_one_index in itemindex_one:
				temp_protein_set_one = self.complexes.complexes[individual_complex_one_index]
				master_overlapped_set = set()
				max_overlap = 0
				max_index = ''

				for individual_complex_zero_index in itemindex_zero:
					temp_protein_set_zero = self.complexes.complexes[individual_complex_zero_index]
					overlapped_set = self.get_overlap_proteins(temp_protein_set_one, temp_protein_set_zero)
					master_overlapped_set = master_overlapped_set | overlapped_set

				if max_index == '':
					max_index = individual_complex_one_index
					max_overlap = len(master_overlapped_set)
				else:
					if max_overlap < len(master_overlapped_set):
						max_index = individual_complex_one_index
						max_overlap = len(master_overlapped_set)

			itemindex_one.remove(max_index)


			for individual_complex_one_index in itemindex_zero:
				temp_protein_set_one = self.complexes.complexes[individual_complex_one_index]
				master_overlapped_set = set()
				max_overlap = 0
				max_index = ''

				for individual_complex_zero_index in itemindex_zero:
					temp_protein_set_zero = self.complexes.complexes[individual_complex_zero_index]
					overlapped_set = self.get_overlap_proteins(temp_protein_set_one, temp_protein_set_zero)
					master_overlapped_set = master_overlapped_set | overlapped_set

				if max_index == '':
					max_index = individual_complex_one_index
					max_overlap = len(master_overlapped_set)
				else:
					if max_overlap < len(master_overlapped_set):
						max_index = individual_complex_one_index
						max_overlap = len(master_overlapped_set)

			itemindex_zero.remove(max_index)

			round += 1

		print "length of complex set one: " + str(len(itemindex_one))
		print "length of complex set two: " + str(len(itemindex_zero))

		# randomize clusters
		#rnd.shuffle(ref_cluster_ids)

		fold_size = int(len(ref_cluster_ids) / num_folds)
		for i in range(num_folds): # only support two_folds cross validation here...

			# create training and evaluating complexes objects.
			evaluation = Goldstandard_from_Complexes("Evaluation")
			training = Goldstandard_from_Complexes("Training")

			if i == 0:
				eval_clust_ids = itemindex_one
				train_clust_ids = itemindex_zero
			else:
				eval_clust_ids = itemindex_zero
				train_clust_ids = itemindex_one

			for train_id in train_clust_ids: training.complexes.addComplex(train_id, self.complexes.complexes[train_id])
			for eval_id in eval_clust_ids: evaluation.complexes.addComplex(eval_id, self.complexes.complexes[eval_id])

			all_eval_prots = evaluation.complexes.get_all_prots()


			if not overlap:
				for comp in training.complexes.complexes:
					training.complexes.complexes[comp] = training.complexes.complexes[comp] - all_eval_prots


			training.make_pos_neg_ppis()
			evaluation.make_pos_neg_ppis()

			train = training.get_goldstandard()
			evaluate = evaluation.get_goldstandard()

			len_train_positive = len(train[0])
			len_eva_positive = len(train[0])
			len_train_negative = len(train[1])
			len_eva_negative = len(train[1])

			len_over_positive = len(train[0] & evaluate[0])
			len_over_negative = len(train[1] & evaluate[1])


			print "number of train and evaluation PPIs:"
			print len_train_positive + len_train_negative
			print "number of overlapped PPIs:"
			print len_over_positive + len_over_negative
			print "I am here"

			out_folds.append((training, evaluation))

		return out_folds

	# def n_fols_split(self, num_folds, overlap="False"):
    #
	# 	ref_cluster_ids = self.complexes.complexes.keys()
	# 	out_folds = []
    #
    #
	# 	# randomize clusters
	# 	rnd.shuffle(ref_cluster_ids)
    #
	# 	fold_size = int(len(ref_cluster_ids) / num_folds)
	# 	for i in range(num_folds):
	# 		# create training and evaluating complexes objects.
	# 		evaluation = Goldstandard_from_Complexes("Evaluation")
	# 		training = Goldstandard_from_Complexes("Training")
	# 		eval_clust_ids = ref_cluster_ids[fold_size * i: min(len(ref_cluster_ids), fold_size * (i + 1))]
	# 		train_clust_ids = list(set(ref_cluster_ids) - set(eval_clust_ids))
    #
	# 		for train_id in train_clust_ids: training.complexes.addComplex(train_id, self.complexes.complexes[train_id])
	# 		for eval_id in eval_clust_ids: evaluation.complexes.addComplex(eval_id, self.complexes.complexes[eval_id])
    #
	# 		all_eval_prots = evaluation.complexes.get_all_prots()
    #
    #
	# 		if not overlap:
	# 			for comp in training.complexes.complexes:
	# 				training.complexes.complexes[comp] = training.complexes.complexes[comp] - all_eval_prots
    #
    #
	# 		training.make_pos_neg_ppis()
	# 		evaluation.make_pos_neg_ppis()
    #
	# 		print ("debugging here")
	# 		train = training.get_goldstandard()
	# 		evaluate = evaluation.get_goldstandard()
    #
	# 		len_train_positive = len(train[0])
	# 		len_eva_positive = len(train[0])
	# 		len_train_negative = len(train[1])
	# 		len_eva_negative = len(train[1])
    #
	# 		len_over_positive = len(train[0] & evaluate[0])
	# 		len_over_negative = len(train[1] & evaluate[1])
    #
    #
	# 		print len_train_positive
	# 		print len_eva_positive
	# 		print len_over_positive
	# 		print len_train_negative
	# 		print len_eva_negative
	# 		print len_over_negative
	# 		sys.exit()
    #
	# 		out_folds.append((training, evaluation))
    #
	# 	return out_folds

	#new function added by Lucas HU to split into n-fold for n-fold cross-validation.
	#just a trial version to test if it works.
	def split_into_n_fold(self, n_fold, val_ppis, no_overlapp=False):  # what is vak_ppis

		tmp_clusters = self.complexes.complexes.keys()
		allPossiblePositive, allPossibleNegative = self.complexes.getPositiveAndNegativeInteractions()

		#get the neagtive PPIs which detected in co-fractionation experiments (the overlap of two sets)
		negativePPIs = val_ppis & allPossibleNegative

		trainingNegatives = set(list(negativePPIs)[:int(len(negativePPIs) / 2)])
		evaluationNegatives = set(list(negativePPIs)[int(len(negativePPIs) / 2):])

		rnd.shuffle(tmp_clusters)

		foldNumberComplex = int(len(tmp_clusters)/n_fold)

		training_evaluation_dictionary = {'turpleKey': []}

		# n_fold cross_validation, split the whole protein complexes set into n_fold, one is for training, and the rest is for validation.
		for i in range(n_fold):

			#create training and evaluating complexes objects.
			evaluation = Goldstandard_from_Complexes("Evaluation")
			training = Goldstandard_from_Complexes("Training")

			evaluatingComplexSet = tmp_clusters[foldNumberComplex * i : foldNumberComplex * (i + 1)]
			trainingComplexSet = list(set(tmp_clusters) - set(evaluatingComplexSet))

			#generate the positive PPIs in training complex set.
			for index in trainingComplexSet:

				tmp_complexes = Clusters(False)
				tmp_complexes.addComplex(index, self.complexes.complexes[index])
				tmp_p, _ = tmp_complexes.getPositiveAndNegativeInteractions()
				tmp_p = tmp_p & val_ppis
				training.complexes.addComplex(index, self.complexes.complexes[index])
				training.add_positive(tmp_p)
				training.add_negative(trainingNegatives)


			#generate the positive PPIs in evaluating complex set.
			for index in evaluatingComplexSet:

				tmp_complexes = Clusters(False)
				tmp_complexes.addComplex(index, self.complexes.complexes[index])
				tmp_p, _ = tmp_complexes.getPositiveAndNegativeInteractions()
				tmp_p = tmp_p & val_ppis
				evaluation.complexes.addComplex(index, self.complexes.complexes[index])
				evaluation.add_positive(tmp_p)
				evaluation.add_negative(evaluationNegatives)

			training.rebalance()
			evaluation.rebalance()

			training_evaluation_dictionary["turpleKey"].append((training, evaluation))

		return training_evaluation_dictionary

	#new function added by Lucas HU to split into n-fold for n-fold cross-validation.
	#In this case, negative PPIs are only generated with each fold of data, so negative PPIs might be positive PPIs in other fold of data
	# a trial version can be used for comparsion (suggested by Florian)
	def split_into_n_fold2(self, n_fold, val_ppis, no_overlapp=False):  # what is vak_ppis

		tmp_clusters = self.complexes.complexes.keys()
		allPossiblePositive, allPossibleNegative = self.complexes.getPositiveAndNegativeInteractions()

		#get the neagtive PPIs which detected in co-fractionation experiments (the overlap of two sets)
		negativePPIs = val_ppis & allPossibleNegative

		rnd.shuffle(tmp_clusters)

		foldNumberComplex = int(len(tmp_clusters)/n_fold)

		training_evaluation_dictionary = {'turpleKey': []}

		# n_fold cross_validation, split the whole protein complexes set into n_fold, one is for training, and the rest is for validation.
		for i in range(n_fold):

			#create training and evaluating complexes objects.
			evaluation = Goldstandard_from_Complexes("Evaluation")
			training = Goldstandard_from_Complexes("Training")

			evaluatingComplexSet = tmp_clusters[foldNumberComplex * i : foldNumberComplex * (i + 1)]
			trainingComplexSet = list(set(tmp_clusters) - set(evaluatingComplexSet))

			#generate the positive PPIs in training complex set.
			tmp_training_complexes = Clusters(False)
			for index in trainingComplexSet:

				tmp_training_complexes.addComplex(index, self.complexes.complexes[index])
				training.complexes.addComplex(index, self.complexes.complexes[index])

			tmp_p, tmp_n = tmp_training_complexes.getPositiveAndNegativeInteractions()
			tmp_p = tmp_p & val_ppis
			tmp_n = tmp_n & val_ppis
			training.add_positive(tmp_p)
			training.add_negative(tmp_n)

			#generate the positive PPIs in evaluating complex set.
			tmp_evaluation_complexes = Clusters(False)
			for index in evaluatingComplexSet:

				tmp_evaluation_complexes.addComplex(index, self.complexes.complexes[index])
				evaluation.complexes.addComplex(index, self.complexes.complexes[index])


			tmp_p2, tmp_n2 = tmp_evaluation_complexes.getPositiveAndNegativeInteractions()
			tmp_p2 = tmp_p2 & val_ppis
			tmp_n2 = tmp_n2 & val_ppis
			evaluation.add_positive(tmp_p2)
			evaluation.add_negative(tmp_n2)

			training.rebalance()
			evaluation.rebalance()

			training_evaluation_dictionary["turpleKey"].append((training, evaluation))

			print "the number of training negatives and positives for corss validation "
			print len(training.get_negative())
			print len(training.get_positive())

		return training_evaluation_dictionary




	def split_into_holdout_training(self, val_ppis, no_overlapp=False):

		holdout = Goldstandard_from_Complexes("Holdout")
		training = Goldstandard_from_Complexes("Training")

		tmp_clusters = self.complexes.complexes.keys()

		rnd.shuffle(tmp_clusters)


		val_negatives = list(self.negative & val_ppis)
		rnd.shuffle(val_negatives)
		t_n = set(val_negatives[:int(len(val_negatives)/2)])
		h_n = set(val_negatives[int(len(val_negatives)/2):])

		t_p, h_p = set([]), set([])
		# Balance data set on positive, since we have way more negative than positive
		i = 0;
		skipped_comp = []
		for complex in tmp_clusters:
			tmp_cluster = Clusters(False)
			tmp_cluster.addComplex(complex, self.complexes.complexes[complex])
			tmp_p, _ = tmp_cluster.getPositiveAndNegativeInteractions()
			tmp_p &= val_ppis #tmp_p = tmp_p & val_ppis
			if len(tmp_p) == 0: #should keep all the complexes
				skipped_comp.append(complex)
				#if i %2 == 0:
				#	training.complexes.addComplex(complex, self.complexes.complexes[complex])
				#else:
				#	holdout.complexes.addComplex(complex, self.complexes.complexes[complex])
				#i += 1
				continue


#
			if len(t_p)<=len(h_p):
				t_p |= tmp_p
				training.complexes.addComplex(complex, self.complexes.complexes[complex])
			else:
				h_p |= tmp_p
				holdout.complexes.addComplex(complex, self.complexes.complexes[complex])

#		for complex in skipped_comp:
#			if len(holdout.complexes.complexes) > len(training.complexes.complexes):
#				training.complexes.addComplex(complex, self.complexes.complexes[complex])
#			else:
#				holdout.complexes.addComplex(complex, self.complexes.complexes[complex])

		training.make_pos_neg_ppis(val_ppis)
		holdout.make_pos_neg_ppis(val_ppis)

		training.negative = t_n
		holdout.negative = h_n

		if no_overlapp:
			training.positive -= holdout.get_edges()

		training.rebalance()
		holdout.rebalance()
		return training, holdout

	# A function added by Lucas HU
	# used for 10_fold_cross_validation_on PPI level
	# a trial version
	def return_gold_standard_complexes(self, val_ppis):

		#holdout = Goldstandard_from_Complexes("Holdout")
		training = Goldstandard_from_Complexes("Training")

		tmp_clusters = self.complexes.complexes.keys()

		rnd.shuffle(tmp_clusters)

		val_negatives = list(self.negative & val_ppis)

		rnd.shuffle(val_negatives)
		t_n = set(val_negatives)

		t_p = set([])

		# Balance data set on positive, since we have way more negative than positive
		skipped_comp = []
		for complex in tmp_clusters:
			tmp_cluster = Clusters(False)
			tmp_cluster.addComplex(complex, self.complexes.complexes[complex])
			tmp_p, _ = tmp_cluster.getPositiveAndNegativeInteractions()
			tmp_p &= val_ppis  # tmp_p = tmp_p & val_ppis
			if len(tmp_p) == 0:  # should keep all the complexes
				skipped_comp.append(complex)

				continue

			t_p |= tmp_p
			training.complexes.addComplex(complex, self.complexes.complexes[complex])

		training.make_pos_neg_ppis(val_ppis)

		training.negative = t_n

		training.rebalance()

		return training

	# @author: Florian Goebels
	# this method combines who CalculateCoElutionScores objects onto one by comping the toMerge object into the self object
	# @Param:
	#		CalculateCoElutionScores toMerge a second CalculateCoElutionScores which should be combined with self object
	#		mode donates how to merge the sets, left (l), right (r), union (u), or  (i)
	def rebalance(self, ratio = 5):
		if len(self.positive) * self.ratio > len(self.negative):
			self.positive = set(rnd.sample(self.positive, int(len(self.negative) / self.ratio)))
			print("Warning: not enough negative data points in reference to create desired ratio pos:%s, neg:%s" % (len(self.positive), len(self.negative)))
		else:
			self.negative = set(rnd.sample(self.negative, len(self.positive)*self.ratio))


class Intact_clusters():

	def __init__(self, need_to_be_mapped, species = "homo_sapiens", name = "IntAct"):
		self.complexes = Clusters(need_to_be_mapped=need_to_be_mapped)
		self.need_to_be_mapped = True
		self.load_data(species)
		self.name = name


	def get_complexes(self):
		return self.complexes

	def load_data(self, species):
		intact_url = "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/2017-04-08/complextab/%s.tsv" % species
		intact_url_FH = urllib2.urlopen(intact_url)
		intact_url_FH.readline()
		i = 0
		for line in intact_url_FH:
			line = line.rstrip()
			linesplit = line.split("\t")
			evidence = linesplit[5]
			if not evidence.startswith("ECO:0000353"): continue
			members = linesplit[4]
			members = re.sub("\(\d+\)", "", members)
			members = re.sub("-\d+", "", members)
			members = re.sub("-PRO_\d+", "", members)
			members = set(members.split("|"))
			self.complexes.addComplex("%i;%s" % (i, linesplit[1]), members)
			i += 1

# @author Florian Goebels
# Wrapper class for downloading and handling CORUM protein complex information taken from here: http://mips.helmholtz-muenchen.de/genre/proj/corum/
class CORUM():
	# @author Florian Goebels
	# object constructor
	# @param
	#		lb lower bound complex should have at least lb members
	#		ub upper bound complex should have at most  ub members
	#		overlap_cutoff merge complexes that have an overlap_score > overlap_cutoff
	#		source_species select for which species the complexes should be maintained
	def __init__(self, need_to_be_mapped, source_species_regex = "(Human|Mammalia)", name = "CORUM"):
		self.complexes = Clusters(need_to_be_mapped=need_to_be_mapped)
		# static regex for identifying valid bochemical evidences codes
		self.biochemical_evidences_regex ="MI:(2193|2192|2191|2197|2195|2194|2199|2198|0807|0401|0400|0406|0405|0404|0089|0084|0081|0007|0006|0004|0513|1029|0979|0009|0008|0841|1312|2188|2189|0411|0412|0413|0928|0415|0417|0098|0729|0920|0921|0603|0602|0605|0604|0402|0095|0096|0606|0091|0092|1142|1145|1147|0019|1309|0696|0697|0695|0858|0698|0699|0425|0424|0420|0423|0991|0990|0993|0992|0995|0994|0997|0996|0999|0998|1028|1011|1010|1314|0027|1313|0029|0028|0227|0226|0225|0900|0901|0430|0434|0435|1008|1009|0989|1004|1005|0984|1007|1000|0983|1002|1229|1087|1325|0034|0030|0031|0972|0879|0870|1036|0678|1031|1035|1034|0676|0440|1138|1236|0049|0048|1232|0047|1137|0419|0963|1026|1003|1022|0808|0515|0514|1187|0516|0511|1183|0512|0887|0880|0889|0115|1006|1249|0982|0953|1001|0508|0509|0657|0814|1190|1191|0813|0066|0892|0899|1211|0108|1218|1352|1354|0949|0946|0947|0073|0071|1019|2168|0700|2167|1252|1017|0276|1189|1184)"
		self.source_species_regex = source_species_regex
		self.getCORUM()
		self.readCORUM()
		self.name = name


	def get_complexes(self):
		return self.complexes

	# @author Florian Goebels
	# downloads current version of corum and safe it to wd/data folder as corum.txt
	def getCORUM(self):
		self.corum_raw = {}

		#corum_url = "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt"
		#corum_url_FH = urllib2.urlopen(corum_url)
		#for line in corum_url_FH:
		#	line = line.rstrip()
		#	linesplit = np.array(line.split("\t"))
		#	corum_id = linesplit[0]
		#	self.corum_raw[corum_id] = linesplit[(2,5,7),]

		#Note, I changed the above part, coz the CORUM database deleted the above link, we need new ways to open it --- Lucas
		corum_url = "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip"

		url = urlopen(corum_url)
		zipfile = ZipFile(StringIO(url.read()))
		for line in zipfile.open("allComplexes.txt").readlines():
			line = line.rstrip()
			linesplit = np.array(line.split("\t"))
			corum_id = linesplit[0] + ";" + linesplit[1]
			self.corum_raw[corum_id] = linesplit[(2,5,7),]

	# @author Florian Goebels
	# reads in CORUM from flat file
	def readCORUM(self):
		for comp in self.corum_raw:
			(species, prots, evidence) = self.corum_raw[comp]
			# bool(...) returns true if evidence code is found => not bool(...) is true if not valid evidence is found, and thus skip this complex
			if not bool(re.search(self.biochemical_evidences_regex, evidence)):
				continue
			if not bool(re.search(self.source_species_regex, species)): continue
			prots = set(prots.split(";"))
			self.complexes.addComplex(comp, prots)


class FileClusters():
	def __init__(self, filepath, valprots, name="File"):
		self.complexes = Clusters(False)
		self.complexes.read_file(filepath)
		self.complexes.remove_proteins(valprots)
		self.name = name

	def get_complexes(self):
		return self.complexes

class Clusters():

	def __init__(self, need_to_be_mapped, overlap_cutoff = 0.8, lb = 3, ub = 50):
		self.complexes = {}
		self.overlap_cutoff = overlap_cutoff
		self.ub = ub
		self.lb = lb
		self.need_to_be_mapped = need_to_be_mapped

	def get_complexes(self):
		return self.complexes

	# added by Lucas Hu, to compare predicted complexes with gold standard complexes
	# to calculate sensitivity...
	def get_overlapped_complexes_set(self, another_complex_set): # complex_set is a Clusters project.
		Cluster_dict_another = another_complex_set.get_complexes()
		overlapped_complex_dic = {}
		for complex_another in Cluster_dict_another:
			for complex in self.complexes:
				if self.overlap(self.complexes[complex], Cluster_dict_another[complex_another]) >= 0.2:
					overlapped_complex_dic[complex_another] = Cluster_dict_another[complex_another]

		return overlapped_complex_dic

	def get_all_prots(self):
		out = set()
		for comp in self.complexes:
			out |= self.complexes[comp]
		return out

	def addComplex(self, complex, members):
		if complex not in self.complexes:
			self.complexes[complex] = set([])
		self.complexes[complex] = self.complexes[complex] | members

	def read_file(self, clusterF):
		clusterFH = open(clusterF)
		i = 0
		all_proteins_count = 0
		for line in clusterFH:
			line = line.rstrip()
			prots = set(line.split("\t"))
			self.addComplex(i, prots)
			i+=1
			all_proteins_count += len(prots)
		clusterFH.close()

		print "Average size of predicted complexes is: " + str((all_proteins_count)/i)

	def write_cuslter_file(self, outF):
		outFH = open(outF, "w")
		print >> outFH, self.to_string()
		outFH.close()

	def to_string(self):
		out = []
		for clust in self.complexes:
			prots = self.complexes[clust]
			out.append("\t".join(prots))
		return "\n".join(out)

	# @author Florian Goebels
	# creats all possible positive and negative protein interactions based on  co-complex membership
	def getPositiveAndNegativeInteractions(self):
		positive = set([])
		negative = set([])
		prot2cluster = self.getProtToComplexMap()
		for protA in prot2cluster:
			for protB in prot2cluster:
				if protA == protB: continue
				edge = "\t".join(sorted([protA, protB]))
				if len(prot2cluster[protA] & prot2cluster[protB]) > 0:
					positive.add(edge)

				else:
					negative.add(edge)
		return positive, negative

	def overlap(self, a, b):
		tmpa = set(a)
		tmpb = set(b)
		overlap = math.pow((len(tmpa & tmpb)), 2) / (len(tmpa) * len(tmpb))
		return overlap

	# @author Florian Goebels
	# merges complexes which have an overlapp score > overlap_cutoff, and continues to merge until there is nothing left to merge
	def merge_complexes(self):
		merged = set()
		allComplexes = self.complexes.keys()
		newComplexes = {}

		for i in range(len(allComplexes)):
			compI = allComplexes[i]
			if compI in merged: continue
			candidates = [compI]
			toMerge = set([compI])
			while(len(candidates)>0):
				this_complex = candidates.pop()
				for j in range(i + 1, len(allComplexes)):
					compJ = allComplexes[j]
					if compJ in merged: continue
					if self.overlap(self.complexes[compI], self.complexes[compJ]) > self.overlap_cutoff:
						toMerge.add(compJ)
						candidates.append(compJ)
				merged.add(this_complex)


			if len(toMerge) > 0:
				(newName, newProts) = set(), set()
				for name in toMerge:
					newName.add(name)
					newProts.update(self.complexes[name])
				newComplexes[(",".join(map(str,newName)),)] = newProts
				merged.add(compI)
			else:
				newComplexes[compI] = self.complexes[compI]
				merged.add(compI)

		self.complexes = newComplexes

	def filter_complexes(self):
		todel = set()
		for comp in self.complexes:
			if len(self.complexes[comp]) < self.lb or len(self.complexes[comp]) > self.ub:
				todel.add(comp)
		for delComp in todel:
			del self.complexes[delComp]
	"""
	#alt version for remove_proteins where complete complex gets removed when not all members have elution data
	def remove_proteins(self, to_keep):
		todel = set([])
		for comp in self.complexes:
			if len(self.complexes[comp] & to_keep)/len(self.complexes[comp])<1:
				todel.add(comp)
		for comp in todel:
			del self.complexes[comp]
		#	self.complexes[comp] = self.complexes[comp] & to_keep
	"""

	def remove_proteins(self, to_keep):
		todel = set([])
		for comp in self.complexes:
			new_prots = self.complexes[comp] & to_keep
			if len(new_prots)==0:
				todel.add(comp)
			else:
				self.complexes[comp] = new_prots

		for comp in todel:
			del self.complexes[comp]

	def getProtToComplexMap(self):
		out = {}
		for cluster in self.complexes:
			for prot in self.complexes[cluster]:
				if prot not in out: out[prot] = set([])
				out[prot].add(cluster)
		return out

	def getOverlapp(self, complexesB, cutoff = 0.5):
		out = 0
		for comp_ID_A in self.complexes.keys():
			protsA = self.complexes[comp_ID_A]
			matched = False
			for comp_ID_B in complexesB.complexes.keys():
				protsB = complexesB.complexes[comp_ID_B]
				if self.overlap(protsA, protsB) > cutoff:
					matched = True
					break

			if matched: out += 1
		return out

	def mmr(self, reference):
		matchingscores = {}
		for ref_complex in reference.complexes:
			if ref_complex not in matchingscores:
				matchingscores[ref_complex] = 0
			for complex in self.complexes:
				matchingscores[ref_complex] = max(matchingscores[ref_complex], self.overlap(self.complexes[complex],  reference.complexes[ref_complex]))
		mmr = matchingscores.values()
		mmr_score = sum(mmr)/len(mmr)
		return mmr_score

	def get_matching_complexes(self, reference):
		out_simcoe = set([])
		out_overlap = set([])
		out_comb = set([])

		def simco(a, b):
			tmpa = set(a)
			tmpb = set(b)
			return len(tmpa & tmpb) / (min(len(tmpa), len(tmpb)))

		for this_complex in self.complexes:
			for ref_complex in reference.complexes:
				overlap_score = self.overlap(self.complexes[this_complex], reference.complexes[ref_complex])
				simco_score = simco(self.complexes[this_complex], reference.complexes[ref_complex])
				if overlap_score >= 0.25:
					out_overlap.add(this_complex)
				if simco_score > 0.5:
					out_simcoe.add(this_complex)
				if (simco_score + overlap_score) / 2 >= 0.375:
					out_comb.add(this_complex)

		return out_simcoe, out_overlap, out_comb

	def frac_match_comp(self, reference):
		if len(self.complexes) == 0:
			return "0\t0\t0"
		else:
			out_simcoe, out_overlap, out_comb = self.get_matching_complexes(reference)
			out_overlap = len(out_overlap)/len(self.complexes)
			out_simcoe = len(out_simcoe)/len(self.complexes)
			out_comb = len(out_comb) / len(self.complexes)
			return "%f\t%f\t%f" % (out_overlap, out_simcoe, out_comb)

	def sensitivity(self, reference):
		sum_all_ref = 0.0
		sum_all_ref_overlaps = 0.0
		for reference_cluster in reference.complexes:
			sum_all_ref += len(reference.complexes[reference_cluster])
			max_overlap_for_this_ref = 0
			for predicted_cluster in self.complexes:
				overlap = len(self.complexes[predicted_cluster] & reference.complexes[reference_cluster])
				if overlap > max_overlap_for_this_ref:
					max_overlap_for_this_ref = overlap
			sum_all_ref_overlaps += max_overlap_for_this_ref
		return sum_all_ref_overlaps/sum_all_ref

	def ppv(self, reference):
		n = len(self.complexes)
		m = len(reference.complexes)
		overlap_mat = np.zeros((n,m))
		for n, predicted_cluster in enumerate(self.complexes):
			for m, reference_cluster in enumerate(reference.complexes):
				overlap = len(self.complexes[predicted_cluster] & reference.complexes[reference_cluster])
				overlap_mat[n,m] = overlap
		if np.sum(overlap_mat) == 0:
			return 0
		else:
			return np.sum(overlap_mat.max(axis=1))/np.sum(overlap_mat)


	def acc(self, reference, sn = None, ppv = None):
		if sn == None: sn = self.sensitivity(reference)
		if ppv == None: ppv = self.ppv(reference)
		return math.sqrt(sn*ppv)

	def clus_sep(self, reference):
		if len(self.complexes) == 0:
			return 0
		n = len(self.complexes)
		m = len(reference.complexes)
		row_F = np.zeros((n, m))
		col_F = np.zeros((n, m))

		for i, compA in enumerate(self.complexes):
			for j, compB in enumerate(reference.complexes):
				overlap = len(self.complexes[compA] & reference.complexes[compB])
				row_F[i,j] = overlap
				col_F[i,j] = overlap

		row_F = np.nan_to_num(row_F/np.sum(row_F, axis=1, keepdims=True))
		col_F = np.nan_to_num(col_F/np.sum(col_F, axis=0, keepdims=True))
		sep = np.nan_to_num(np.sum(row_F*col_F))
		sep_co = sep/m
		sep_cl = sep/n
		return math.sqrt(sep_co*sep_cl)

	def clus_eval(self, ref):
		mmr = self.mmr(ref)
		ppv = self.ppv(ref)
		sn = self.sensitivity(ref)
		acc = self.acc(ref, sn, ppv)
		sep = self.clus_sep(ref)
		prc = self.frac_match_comp(ref)
		return mmr, prc, sn, ppv, acc, sep


# @author Florian Goebels
# Wrapper class for retrieving go annotation for a given taxid from the QuickGo webservice : https://www.ebi.ac.uk/QuickGO/
class QuickGO():
	# @author Florian Goebels
	# object constructor makes internet connection and downloads go annotations into the wd/go_files folder as taxid.go file
	# @param
	#		taxid of species that go annotation should be downloaded
	def __init__(self, taxid, need_to_be_mapped, name ="QuickGO"):
		self.taxid = taxid
		self.complexes = Clusters(need_to_be_mapped=need_to_be_mapped)
		self.get_GO_complexes()
		self.name = name

	# @author LUCAS HU
	# Since QuickGO has been updated, format has been changed, we need to re-write this function to extract protein complexes information from GO
	# reads in go flat file if gaf 20 format as protein to go annotation mapping (as dictonary)
	# @param
	#		taxid species for which go annotation should be read into memory
	def get_GO_complexes(self):
		go_to_prot_map = {}
		prot_to_go_map = {}
		####
		requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?evidenceCode=ECO%3A0000353%2CECO%3A0000314%2CECO%3A0000269&goId=GO%3A0043234&taxonId="+str(self.taxid)+"&evidenceCodeUsage=descendants&evidenceCodeUsageRelationships=is_a"

		r = requests.get(requestURL, headers={"Accept": "text/gpad"})

		if not r.ok:
			r.raise_for_status()
			sys.exit()

		responseBody = r.text

		lines = responseBody.split("\n")
		for line in lines:
			linesplit = line.split("\t")
			if linesplit[0] == "UniProtKB":
				prot = linesplit[1]
				go_complex = linesplit[3]
				date = int(linesplit[8])
				if date > 20170512: continue
				# Adding prot to go map
				if prot not in prot_to_go_map: prot_to_go_map[prot] = set([])
				prot_to_go_map[prot].add(go_complex)
				# Adding go to prot map
				if go_complex not in go_to_prot_map: go_to_prot_map[go_complex] = set([])
				go_to_prot_map[go_complex].add(prot)

		i = 0
		for go_complex in go_to_prot_map:
			self.complexes.addComplex(go_complex, go_to_prot_map[go_complex])
			i+=1

	def get_complexes(self):
		return self.complexes


	# # @author Florian Goebels
	# # reads in go flat file if gaf 20 format as protein to go annotation mapping (as dictonary)
	# # @param
	# #		taxid species for which go annotation should be read into memory
	# def get_GO_complexes(self):
	# 	go_to_prot_map = {}
	# 	prot_to_go_map = {}
	# 	quickgoURL = "http://www.ebi.ac.uk/QuickGO-Old/GAnnotation?goid=GO:0043234&tax=%s&format=tsv&limit=1000000000&evidence=IDA,IPI,EXP," % (self.taxid)
	# 	print quickgoURL
	# 	print self.taxid
	# 	print "the url is: ..."
	# 	quickgoURL_FH = urllib2.urlopen(quickgoURL)
	# 	quickgoURL_FH.readline()
	# 	for line in quickgoURL_FH:
	# 		line = line.rstrip()
	# 		linesplit = line.split("\t")
	# 		prot = linesplit[1]
	# 		go_complex = linesplit[6] + ";" + linesplit[7]
	# 		date = int(linesplit[12])
	# 		if date > 20170512: continue
	# 		# Adding prot to go map
	# 		if prot not in prot_to_go_map: prot_to_go_map[prot] = set([])
	# 		prot_to_go_map[prot].add(go_complex)
	# 		# Adding go to prot map
	# 		if go_complex not in go_to_prot_map: go_to_prot_map[go_complex] = set([])
	# 		go_to_prot_map[go_complex].add(prot)
    #
	# 	quickgoURL_FH.close()
	# 	i = 0
	# 	for go_complex in go_to_prot_map:
	# 		self.complexes.addComplex(go_complex, go_to_prot_map[go_complex])
	# 		i+=1
    #
	# def get_complexes(self):
	# 	return self.complexes



# @author Florian Goebels
# Class for handling and getting data from INPARANOID http://inparanoid.sbc.su.se/cgi-bin/index.cgi
class Inparanoid():
	# @author Florian Goebels
	# init function fetches inparanoid data and read it in for creating CORUM reference data
	# thus it always parses and download human to target species mapping
	# @param
	#		taxid taxid of given target species from which should be mapped to human
	#		inparanoid_cutoff cut off used for accepting ortholog mappings TODO include bootstrapping score as well
	#		foundProts list of proteins that were found in the MS experiments. If given, in cases of uncertain ortholog mapping (e.g. more than one mapping if inparanoid score == 1), ortholog mappings to found prots are prefered
	def __init__(self, taxid, inparanoid_cutoff=1, foundProts = set([])):
		self.taxid2Name = self.getTaxId2Namemapping()
		self.inparanoid_cutoff = inparanoid_cutoff
		self.foundProts = foundProts
		if taxid != "":
			if taxid in self.taxid2Name:
				self.species = self.taxid2Name[taxid]
				xmldoc = self.getXML()
				self.orthmap, self.orthgroups = self.parseXML(xmldoc)
			else:
				print "Taxid:%s not supported" % taxid

	def mapProtein(self, prot):
		if prot not in self.orthmap: return None
		return self.orthmap[prot]


	# @author Florian Goebels
	# mappes protein interactions to their respectiv ortholog counterpart
	# here it maps human reference protein interactions from CORUM to given target species
	# @param
	def mapEdges(self, edges):
		mapped_edges = set([])
		for edge in edges:
			protA, protB = edge.split("\t")
			if protA not in self.orthmap or protB not in self.orthmap: continue
			edge = "\t".join(sorted([self.orthmap[protA], self.orthmap[protB]]))
			mapped_edges.add(tuple(edge))
		return mapped_edges

	def mapEdge(self, edge):
		protA, protB = edge.split("\t")
		if protA not in self.orthmap or protB not in self.orthmap: return ""
		return "\t".join(sorted([self.orthmap[protA], self.orthmap[protB]]))

	def mapComplexes(self, clusters):
		todel = set([])
		for clust in clusters.complexes:
			mapped_members = set([])
			for prot in clusters.complexes[clust]:
				if prot in self.orthmap:
					mapped_members.add(self.orthmap[prot])
#				else:
#					print "No map for %s" % prot

			if len(mapped_members)==0:
				todel.add(clust)
			else:
				clusters.complexes[clust] = mapped_members
		for clust in todel:
			del clusters.complexes[clust]

	# @author Florian Goebels
	# get taxid to inparanoid name mapping from the inparanoid website
	def getTaxId2Namemapping(self):
		#TODO do not hard code
		url_str = "http://inparanoid.sbc.su.se/download/current/sequences/species.mapping.inparanoid8"
		url_FH = urllib2.urlopen(url_str)
		taxid2Name = {}
		for line in url_FH:
			line = line.rstrip()
			(taxid, name) = line.split(".fasta\t")
			taxid2Name[taxid] = name
		return taxid2Name

	# @author Florian Goebels
	# reads in ortholog mapping from the given inparanoid XML object
	# @param:
	#		xmldoc parsed Inparanoid xml doc
	def parseXML(self, xmldoc):
		protID2prot = {}
		orthgroups = []
		targetGenes = set([])
		for species in xmldoc.getElementsByTagName('species'): 
			for protlist in species.getElementsByTagName('genes'):
				for prot in protlist.getElementsByTagName('gene'):
					if species.getAttribute('NCBITaxId') != "9606": targetGenes.add(str(prot.getAttribute('protId')))
					protID2prot[prot.getAttribute('id')] = prot.getAttribute('protId')
		for orthgroup in xmldoc.getElementsByTagName('orthologGroup'):
			protsInGroup = set([])
			for protRef in orthgroup.getElementsByTagName('geneRef'):
				prot = str(protID2prot[protRef.getAttribute('id')])
				inparanoidScore = float(protRef.getElementsByTagName('score')[0].getAttribute('value'))
	#			TODO integrate bootstraping score
	#			bootstrapScore = protRef.getElementsByTagName('score')[1].getAttribute('value')
				if inparanoidScore >= self.inparanoid_cutoff:
					protsInGroup.add(prot)
			if len(protsInGroup)>1:
				orthgroups.append(protsInGroup)

		if len(self.foundProts) != 0:
			toDel = targetGenes - self.foundProts
			for i in range(len(orthgroups)):
				orthgroups[i] = orthgroups[i] - toDel
		outmap = {}
		outgroups = []
		for orthgroup in orthgroups:
			if len(orthgroup) == 2:
				protA, protB = orthgroup
				if protA in targetGenes:
					outmap[protB] = protA
				else:
					outmap[protA] = protB
				outgroups.append(orthgroup)
		return outmap, outgroups

	# @author Florian Goebels
	# fetches human to target Imparanoid xml doc from the inparanoid web server
	def getXML(self):
		url_str = 'http://inparanoid.sbc.su.se/download/current/Orthologs_OrthoXML/%s/%s-%s.orthoXML'
		first = self.species
		second = "H.sapiens"	
		if self.species > "H.sapiens":
			first = "H.sapiens"
			second = self.species
		url_str = url_str % (first, first, second)
		xml_str = urllib2.urlopen(url_str).read()
		xmldoc = minidom.parseString(xml_str)
		return xmldoc

	def readTable(self, tableF, direction = 0):
		tableFH = open(tableF)
		targetGenes = set([])
		orthgroups = []
		def getids(ids_raw):
			out = []
			for i in range(int(len(ids_raw)/2)):
				this_id = ids_raw[i*2]
				this_score = float(ids_raw[(i*2)+1])

				if this_score < self.inparanoid_cutoff:continue
				out.append(this_id)
			return out
		tableFH.readline()
		for line in tableFH:
			if direction == 0:
				orthID, score, target_orth, source_orth = line.rstrip().split("\t")
			else:
				orthID, score, source_orth, target_orth = line.rstrip().split("\t")
			targetids = getids(target_orth.split())
			sourceids = getids(source_orth.split())
			group = []
			group.extend(targetids)
			group.extend(sourceids)
			orthgroups.append(set(group))
			targetGenes |= set(targetids)
		tableFH.close()

		if len(self.foundProts) != 0:
			toDel = targetGenes - self.foundProts
			for i in range(len(orthgroups)):
				orthgroups[i] = orthgroups[i] - toDel
		outmap = {}
		outgroups = []
		for orthgroup in orthgroups:
			if len(orthgroup) == 2:
				protA, protB = orthgroup
				if protA in targetGenes:
					outmap[protB] = protA
				else:
					outmap[protA] = protB
				outgroups.append(orthgroup)
		self.orthmap, self.orthgroups =  outmap, outgroups