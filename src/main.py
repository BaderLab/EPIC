from __future__ import division

import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
import benchmark as bench
import sys
import copy
import os
import numpy as np
import argparse

import warnings
warnings.filterwarnings('ignore')

def Goldstandard_from_cluster_File(gsF, foundprots=""):
	clusters = GS.Clusters(need_to_be_mapped=False)
	clusters.read_file(gsF)
	if foundprots != "": clusters.remove_proteins(foundprots)
	gs = GS.Goldstandard_from_Complexes("All")
	gs.complexes = clusters
	gs.make_pos_neg_ppis()
	return gs


def Goldstandard_from_PPI_File(gsF, foundprots=""):
	out = GS.Goldstandard_from_Complexes("gs")
	gsFH = open(gsF)
	for line in gsFH:
		line = line.rstrip()
		ida, idb, class_label = line.split("\t")[0:3]
		if foundprots !="" and (ida not in foundprots or idb not in foundprots): continue
		edge = "\t".join(sorted([ida, idb]))
		if class_label == "positive":
			out.positive.add(edge)
		else:
			out.negative.add(edge)
	gsFH.close()
	return out

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--feature_selection", type = str, help="Select which features to use. This is an 8 position long array of 0 and 1, where each position determines which co-elution feature to use. Features sorted by position are: MI, Bayes, Euclidean, WCC, Jaccard, PCCN, PCC, and Apex.  Each default=11101001", default="11101001")
	parser.add_argument("input_dir",  type = str, help="Directory containing the elution files for each experiment")
	parser.add_argument("-S", "--source", type = str, help="Flag for telling EPIC what input reference type to expect. TAXID: automatically download reference from GO,CORUM,INtACT. CLUST expect protein cluster flat file. PPI expects PPI flat file. Values: TAXID, CLUST, PPI, default: TAXID",
						default="TAXID")
	parser.add_argument("reference", type = str,help="Taxid, or file location of the used clust/ppi reference. When not using taxid, it is required to change the --source argument to either PPI or CLUST. default a taxonomic id")
	parser.add_argument("output_dir", type = str,help="Directory containing the output files")
	parser.add_argument("-o", "--output_prefix", type = str,help="Prefix name for all output Files", default="Out")

	parser.add_argument("-m", "--mode", type = str,help="Run EPIC with experimental, functional, or both evidences. Values: EXP, FA, COMB, default: EXP  ",
						default="EXP")
	parser.add_argument("-M", "--classifier", type = str,help="Select which classifier to use. Values: RF SVM, default RF",
						default="RF")
	parser.add_argument("-n", "--num_cores", type = int,help="Number of cores to be used, default 1",
						default=1)
	parser.add_argument("-f", "--fun_anno_source", type = str,help="Where to get functional annotaiton from. Values: STRING or GM or File, default= GM",
						default="GM")
	parser.add_argument("-F", "--fun_anno_file", type=str,
						help="Path to File containing functional annotation. This flag needs to be set when using FILE as fun_anno_source.",
						)
	parser.add_argument("-c", "--co_elution_cutoff", type = float,help="Co-elution score cutoff. default 0.5",
						default=0.5)
	parser.add_argument("-C", "--classifier_cutoff", type = float,help="Classifier confidence valye cutoff. default = 0.5",
						default=0.5)
	parser.add_argument("-e", "--elution_max_count", type = int,help="Removies protein that have a maximal peptide count less than the given value. default = 1",
						default=1)
	parser.add_argument("-p", "--precalcualted_score_file", type = str,help="Path to precalulated scorefile to read scores from for faster rerunning of EPIC. default = None",
						default="NONE")

	args = parser.parse_args()

	#Create feature combination
 	if args.feature_selection == "00000000":
		print "Select at least one feature"
		sys.exit()

	this_scores = utils.get_fs_comb(args.feature_selection)
	print "\t".join([fs.name for fs in this_scores])

	# Initialize CLF
 	use_rf = args.classifier == "RF"
	clf = CS.CLF_Wrapper(args.num_cores, use_rf)

	# Load elution data
 	foundprots, elution_datas = utils.load_data(args.input_dir, this_scores)

	gs = ""
	# Generate reference data set
 	if args.source == "TAXID":
		print "Geting PPIs from CORUM,GO,INTACT %s" % args.reference
		gs = utils.create_goldstandard(args.reference, foundprots)
	elif args.source == "CLUST":
		print "Reading cluster file from %s" % args.reference
		gs = Goldstandard_from_cluster_File(args.reference, foundprots)
	elif args.source == "PPI":
		gs = Goldstandard_from_PPI_File(args.reference, foundprots)
	else:
		print "Invalid reference source please select TAXID, CLUST, or PPI"
		sys.exit()

	output_dir = args.output_dir + os.sep + args.output_prefix

	scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=args.num_cores, cutoff= args.co_elution_cutoff)
	if args.precalcualted_score_file == "NONE":
		scoreCalc.calculate_coelutionDatas(gs)
	else:
 		scoreCalc.readTable(args.precalcualted_score_file, gs)

	print scoreCalc.scores.shape

	functionalData = ""
	gs.positive = set(gs.positive & set(scoreCalc.ppiToIndex.keys()))
	gs.negative = set(gs.negative & set(scoreCalc.ppiToIndex.keys()))
	gs.rebalance()

	print len(gs.positive)
	print len(gs.negative)

	if args.mode != "EXP":
		functionalData = utils.get_FA_data(args.fun_anno_source, args.fun_anno_file)
		print "Dimension of fun anno " + str(functionalData.scores.shape)
		tmp_sc = copy.deepcopy(scoreCalc)
		tmp_sc.add_fun_anno(functionalData)
		print "Start benchmarking"
		utils.cv_bench_clf(tmp_sc, clf, gs, output_dir, format="png")
	else:
		utils.cv_bench_clf(scoreCalc, clf, gs, output_dir, format="png")

	network = utils.make_predictions(scoreCalc, args.mode, clf, gs, fun_anno=functionalData)

	# Predict protein interaction
	outFH = open("%s.pred.txt" % (output_dir), "w")

	final_network = []
	for PPI in network:
		items = PPI.split("\t")
		if float(items[2]) >= args.classifier_cutoff:
			final_network.append(PPI)

	print >> outFH, "\n".join(final_network)
	outFH.close()

	# Predicting clusters
	utils.predict_clusters("%s.pred.txt" % (output_dir), "%s.clust.txt" % (output_dir))


	# Evaluating predicted clusters
	pred_clusters = GS.Clusters(False)
	pred_clusters.read_file("%s.clust.txt" % (output_dir))
	clust_scores, header = utils.clustering_evaluation(gs.complexes, pred_clusters, "", False)
	outFH = open("%s.eval.txt" % (output_dir), "w")
	header = header.split("\t")
	clust_scores = clust_scores.split("\t")
	for i, head in enumerate(header):
		print "%s\t%s" % (head, clust_scores[i])
		print >> outFH, "%s\t%s" % (head, clust_scores[i])
	outFH.close()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)
