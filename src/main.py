from __future__ import division

import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
import benchmark as bench
import sys
import os
import numpy as np



def Goldstandard_from_cluster_File(gsF, foundprots=""):
	clusters = GS.Clusters(need_to_be_mapped=False)
	clusters.read_file(gsF)
	if foundprots != "": clusters.remove_proteins(foundprots)
	gs = GS.Goldstandard_from_Complexes("All")
	gs.complexes = clusters
	gs.make_pos_neg_ppis()
	return gs

# def main():
# 	feature_combination, input_dir, use_rf, num_cores, mode, anno_source, anno_F, target_taxid, refF, output_dir = sys.argv[1:]
#
# 	#Create feature combination
# 	if feature_combination == "00000000": sys.exit()
# 	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
# 	this_scores = []
# 	for i, feature_selection in enumerate(feature_combination):
# 		if feature_selection == "1": this_scores.append(scores[i])
#
# 	print "\t".join([fs.name for fs in this_scores])
#
# 	# Initialize CLF
# 	use_rf = use_rf == "True"
# 	num_cores = int(num_cores)
# 	clf = CS.CLF_Wrapper(num_cores, use_rf)
#
# 	# Load elution data
# 	foundprots, elution_datas = utils.load_data(input_dir, this_scores)
#
# 	# Generate reference data set
# 	if refF == "none":
# 		all_gs = utils.create_goldstandard(target_taxid, foundprots)
# 	else:
# 		all_gs = Goldstandard_from_cluster_File(refF, foundprots)
#
# 	scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=num_cores, cutoff= 0.5)
# 	#scoreCalc.calculate_coelutionDatas(all_gs)
# 	scoreCalc.readTable(output_dir + ".scores.txt", all_gs)
#
# 	print "training ppis: %i" % len(set(scoreCalc.ppiToIndex.keys()))
#
# 	#n_fold cross validation to test the stability of preicted PPIs
# 	#utils.stability_evaluation(10, all_gs, scoreCalc, clf, output_dir, mode, anno_source, anno_F)
# 	#sys.exit()
#
# 	#n_fold cross validation to select the best features.
# 	#n_fold_cross_validation(10, all_gs, scoreCalc, clf, output_dir, mode, anno_source, anno_F)
#
# 	#sys.exit()
#
# 	###### actually predict the network using all data
# 	train, eval = all_gs.split_into_holdout_training(set(scoreCalc.ppiToIndex.keys()))
#
# 	print "All comp:%i" % len(all_gs.complexes.complexes)
# 	print "Train comp:%i" % len(train.complexes.complexes)
# 	print "Eval comp:%i" % len(eval.complexes.complexes)
#
# 	print "Num valid ppis in training pos: %i" % len(train.positive)
# 	print "Num valid ppis in training neg: %i" % len(train.negative)
# 	print "Num valid ppis in eval pos: %i" % len(eval.positive)
# 	print "Num valid ppis in eval neg: %i" % len(eval.negative)
#
# 	# Evaluate classifier
# 	#utils.bench_clf(scoreCalc, train, eval, clf, output_dir, verbose=True)
#
# 	functionalData = ""
# 	if mode != "exp":
# 		functionalData = utils.get_FA_data(anno_source, anno_F)
#
# 	print functionalData.scores.shape
#
# 	all_gs.rebalance()
# 	# Predict protein interaction
# 	network = utils.make_predictions(scoreCalc, mode, clf, all_gs, functionalData)
# 	outFH = open("%s.%s.pred.txt" % (output_dir, mode + anno_source), "w")
# 	print >> outFH, "\n".join(network)
# 	outFH.close()
#
# 	# Predicting clusters
# 	utils.predict_clusters("%s.%s.pred.txt" % (output_dir, mode + anno_source), "%s.%s.clust.txt" % (output_dir, mode + anno_source))
#
# 	# Evaluating predicted clusters
# 	pred_clusters = GS.Clusters(False)
# 	pred_clusters.read_file("%s.%s.clust.txt" % (output_dir, mode + anno_source))
# #	utils.clustering_evaluation(train.complexes, pred_clusters, "Train", True)
# 	clusterEvaluationScores = utils.clustering_evaluation(eval.complexes, pred_clusters, "", True)
# 	outFH = open("%s.%s.evaluation.txt" % (output_dir, mode + anno_source), "w")
#
# 	head = clusterEvaluationScores[1]
# 	cluster_scores = clusterEvaluationScores[0]
#
# 	tmp_head = head.split("\t")
# 	tmp_scores = cluster_scores.split("\t")
# 	for i in range(len(tmp_head)):
# 		outFH.write("%s\t%s" % (tmp_head[i], tmp_scores[i]))
# 		outFH.write("\n")

def main():

	input_dir, use_rf, num_cores, mode, anno_source, anno_F, target_taxid, refF, output_dir = sys.argv[1:]

	# read the whole MSB score files, and store it as CalculateCoElutionScores object
	# set all 0 to 1, so it reads all MSB scores
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5),
			  CS.Pearson(), CS.Apex()]
	this_scores = []
	for i, feature_selection in enumerate("11111111"):
		if feature_selection == "1": this_scores.append(scores[i])

	# read elution data
	foundprots, elution_datas = utils.load_data(input_dir, this_scores)

	# Generate reference data set
	if refF == "none":
		all_gs = utils.create_goldstandard(target_taxid, foundprots)
	else:
		all_gs = Goldstandard_from_cluster_File(refF, foundprots)



	all_scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=num_cores,
											cutoff=0.5)

	all_scoreCalc.readTable(output_dir + ".scores.txt", all_gs)

	feature_list = utils.generate_all_corr_combination(8)
	#feature_list = ["10000010", "00000001"] # just test this one

	# get all proteins in the reference data set.
	valprots = all_gs.get_proteins()

	# write the F1, auc_PR, auc_ROC to a file.
	outFH = open("%s.%s.evaluation.txt" % (output_dir, "PPI_feature_selection"), "w")
	outFH.write("%s\t%s\t%s\t%s" % ("features","F_measure", "Area_under_PR", "Area_under_ROC"))
	outFH.write("\n")

	for item_corr_combination in feature_list:

		# create feature combination
		if item_corr_combination == "00000000": sys.exit()
		scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5),
				  CS.Pearson(), CS.Apex()]
		this_scores = []
		for i, feature_selection in enumerate(item_corr_combination):
			if feature_selection == "1": this_scores.append(scores[i])

		print "\t".join([fs.name for fs in this_scores])

		# Initialize CLF
		if use_rf == "True":
			USE_rf = True
		else:
			USE_rf = False

		num_cores = int(num_cores)
		clf = CS.CLF_Wrapper(num_cores, USE_rf)

		# get the feature names and put into an empty set
		feature_names = set()
		for item_score in this_scores:
			feature_names.add(item_score.get_name())

		# cut the features to what we want to test
		test_scoreCalc = bench.feature_selector(feature_names, all_scoreCalc, valprots)

		print ("The size of chopped matrix for selected features")
		print np.shape(test_scoreCalc.get_scoreCalc().get_all_scores())

		print "training ppis: %i" % len(set(test_scoreCalc.ppiToIndex.keys()))

		train_gold_complexes = all_gs.return_gold_standard_complexes(set(test_scoreCalc.ppiToIndex.keys()))

		print "Train_gold comp:%i" % len(train_gold_complexes.complexes.complexes)

		print "Num valid ppis in pos: %i" % len(train_gold_complexes.positive)
		print "Num valid ppis in neg: %i" % len(train_gold_complexes.negative)

		# Evaluate classifier
		evaluation_results = utils.bench_by_PPI_clf(10, test_scoreCalc, train_gold_complexes, output_dir, clf, verbose=True)

		print evaluation_results

		outFH.write("%s\t%s\t%s\t%s" % (item_corr_combination,evaluation_results[0], evaluation_results[1],evaluation_results[2]))
		outFH.write("\n")


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)
