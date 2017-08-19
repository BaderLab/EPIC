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

def main():
	feature_combination, input_dir, use_rf, num_cores, mode, anno_source, anno_F, target_taxid, refF, output_dir = sys.argv[1:]

	#Create feature combination
 	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	this_scores = []
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])

	print "\t".join([fs.name for fs in this_scores])

	# Initialize CLF
 	use_rf = use_rf == "True"
	num_cores = int(num_cores)
	clf = CS.CLF_Wrapper(num_cores, use_rf)

	# Load elution data
 	foundprots, elution_datas = utils.load_data(input_dir, [])

	# Generate reference data set
 	if refF == "none":
		all_gs = utils.create_goldstandard(target_taxid, foundprots)
	else:
		all_gs = Goldstandard_from_cluster_File(refF, foundprots)

	scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=num_cores, cutoff= 0.5)
	#scoreCalc.calculate_coelutionDatas(all_gs)
 	scoreCalc.readTable(output_dir + ".scores.txt", all_gs)
	functionalData = ""
	if mode != "exp":
		functionalData = utils.get_FA_data(anno_source, anno_F)
		print functionalData.scores.shape

	all_gs.rebalance()
	# Predict protein interaction
#	network = utils.make_predictions(scoreCalc, mode, clf, all_gs, functionalData)
#	outFH = open("%s.%s.pred.txt" % (output_dir, mode + anno_source), "w")
#	print >> outFH, "\n".join(network)
#	outFH.close()

	# Predicting clusters
#	utils.predict_clusters("%s.%s.pred.txt" % (output_dir, mode + anno_source), "%s.%s.clust.txt" % (output_dir, mode + anno_source))

	# Evaluating predicted clusters
	pred_clusters = GS.Clusters(False)
	pred_clusters.read_file("%s.%s.clust.txt" % (output_dir, mode + anno_source))
	clusterEvaluationScores = utils.clustering_evaluation(all_gs.complexes, pred_clusters, "", True)
	outFH = open("%s.%s.evaluation.txt" % (output_dir, mode + anno_source), "w")

	head = clusterEvaluationScores[1]
	cluster_scores = clusterEvaluationScores[0]

	tmp_head = head.split("\t")
	tmp_scores = cluster_scores.split("\t")
	for i in range(len(tmp_head)):
		outFH.write("%s\t%s" % (tmp_head[i], tmp_scores[i]))
		outFH.write("\n")

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)
