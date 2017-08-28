from __future__ import division
import numpy as np
import copy, os, sys, re
import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.stats import zscore
import glob, math
import random as rnd

# a function added by lucas, to use n_fold cross_validation to help select features.
# a trial version though.
def n_fold_cross_validation(n_fold, all_gs, scoreCalc, clf, output_dir, overlap, local):
	out_scores = []
	out_head = []

	train_eval_container = all_gs.n_fols_split(n_fold, overlap)

	# create a matrix to store the computed complexes vealuation metrics
	complex_eval_score_vector = np.zeros((n_fold,10))
	val_ppis = set(scoreCalc.ppiToIndex.keys())

	print "Number of ppis with e-score>0.5: %i" % len(val_ppis)
	#the global cluster will contain all clusters predcited from n-fold-corss validation
	for index in range(n_fold):
		print "processinng fold " + str(index + 1)
		train, eval = train_eval_container[index]

		train.positive = train.positive & val_ppis
		train.negative = train.negative & val_ppis
		train.rebalance()

		print "All comp:%i" % len(all_gs.complexes.complexes)
		print "Train comp:%i" % len(train.complexes.complexes)
		print "Eval comp:%i" % len(eval.complexes.complexes)
		print "Num valid ppis in training pos: %i" % len(train.positive & val_ppis)
		print "Num valid ppis in training neg: %i" % len(train.negative & val_ppis)
		print "Num valid ppis in eval pos: %i" % len(eval.positive )
		print "Num valid ppis in eval neg: %i" % len(eval.negative 	)



		print "Overlap positive %i" % ( len(train.positive & eval.positive))
		print "Overlap negative %i" % ( len(train.negative & eval.negative))


		netF = "%s.fold_%s.pred.txt" % (output_dir, index)
		clustF = "%s.fold_%s.clust.txt" % (output_dir, index)

		if local :
		# Evaluate classifier
		# utils.bench_clf(scoreCalc, train, eval, clf, output_dir, verbose=True)
		# Predict protein interaction based on n_fold cross validation
			network = utils.make_predictions_cross_validation(scoreCalc, train, eval, clf)

		else:
			network = utils.predictInteractions(scoreCalc, clf, train, verbose=True)


		if len(network) == 0:
			print "No edges were predicted"
			tmp_scores = [0]*10
			prefix = "Fold %i " % (index+1)
			tmp_head = "\t".join(["%s%s" % (prefix, h) for h in
							  ["Num_pred_PPIS", "NUM_pred_CLUST", "mmr", "overlapp", "simcoe", "mean_simcoe_overlap", "sensetivity", "ppv", "accuracy",
							   "sep"]])
			out_head.append(tmp_head)
			out_scores.append("\t".join(map(str,tmp_scores)))
			complex_eval_score_vector[index, :] = tmp_scores
			continue

		tmp = []
		for ppi in network:
			prota, protb, score =  ppi.split("\t")
			if float(score)>0.5: # this is random forest confidence cut off
				tmp.append(ppi)
		network = tmp

		outFH = open(netF, "w")
		print >> outFH, "\n".join(network)
		outFH.close()

		# Predicting clusters
		utils.predict_clusters(netF, clustF)


		# Evaluating predicted clusters
		pred_clusters = GS.Clusters(False)
		pred_clusters.read_file(clustF)

		print "number of complexes"
		print len(pred_clusters.get_complexes())

		print "number of ppis"
		print len(network)

		fold_scores, fold_head = utils.clustering_evaluation(all_gs.complexes, pred_clusters, "Fold %i " % (index+1), True)
		out_scores.append("%i\t%i\t%s" % (len(network), len(pred_clusters.get_complexes()), fold_scores))
		out_head.append("Fold %i Num_pred_PPIS\tFold %i NUM_pred_CLUST\t%s" % ((index+1), (index+1), fold_head))

		tmp_scores = [len(network), len(pred_clusters.get_complexes())]
		tmp_scores.extend(map(float, fold_scores.split("\t")))
		tmp_scores = np.array(tmp_scores)
		complex_eval_score_vector[index, :] = tmp_scores

	averaged_complex_eval_metrics_vector = np.mean(complex_eval_score_vector, axis = 0)

	out_scores.append("\t".join(map(str, averaged_complex_eval_metrics_vector)))
	mean_head = re.sub("Fold \d", "Mean", "Fold %i Num_pred_PPIS\tFold %i NUM_pred_CLUST\t%s" % ((index+1), (index+1), fold_head))

	out_head.append(mean_head)
	return "\t".join(out_scores), "\t".join(out_head)
	#return averaged_complex_eval_metrics_vector, fold_head.split("\t")

def cut(args):
	fc, scoreF, outF = args
	if fc == "00000000": sys.exit()
	this_scores = get_fs_comb(fc)
	scoreCalc = CS.CalculateCoElutionScores("", "", "","", cutoff=0.5)
	empty_gs = GS.Goldstandard_from_Complexes()
	empty_gs.positive = set([])
	empty_gs.negative = set([])
	scoreCalc.readTable(scoreF, empty_gs)
	print scoreCalc.to_predict
	feature_comb = feature_selector([fs.name for fs in this_scores], scoreCalc)
	feature_comb.open()
	outFH = open(outF, "w")
	print >> outFH, "\t".join(feature_comb.scoreCalc.header)
	for i in range(feature_comb.to_predict):
		edge, edge_scores = feature_comb.get_next()
		if edge == "" or edge_scores == []: continue
		print >> outFH, "%s\t%s" % (edge, "\t".join(map(str, edge_scores)))
	outFH.close()
	feature_comb.close()


def merge_MS(args):
	def read_scores(scoreF, cutoff):
		num_prots = CS.lineCount(scoreF)
		scoreFH = open(scoreF)
		header = scoreFH.readline().rstrip()
		header = header.split("\t")
		out = CS.CalculateCoElutionScores("", "", "", 4)
		out.scores = np.zeros((num_prots , len(header[2:])))
		out.header = header
		i = 0
		for line in scoreFH:
			line = line.rstrip()
			if line == "":continue
			line = line.split("\t")
			edge = "\t".join(line[:2])
			this_score = np.array(map(float, line[2:]))
			if len(list(set(np.where(this_score >= cutoff)[0]))) > 0:
				out.ppiToIndex[edge] = i
				out.IndexToPpi[i] = edge
				out.scores[i, :] = this_score
				i += 1
		out.scores = out.scores[0:i, :]
		print i
		return out

	ms1_in, ms2_in, mode, outF = args
	ms1cutoff = 0
	ms2cutoff = 0

	if mode == "i":
		ms1cutoff = 0.5
		ms2cutoff = 0.5
	if mode == "u":
		ms1cutoff = 0
		ms2cutoff = 0
	if mode == "l":
		ms1cutoff = 0
		ms2cutoff = 0.5
	if mode == "r":
		ms1cutoff = 0.5
		ms2cutoff = 0

	ms1 = read_scores(ms1_in, ms1cutoff)
	print "Done reading in MS1"
	print ms1.scores.shape

	ms2 = read_scores(ms2_in, ms2cutoff)
	print "Done reading in MS2"
	print ms2.scores.shape


	ms2.merge(ms1, mode)
	print "Done merging MS1 and MS2"
	print ms2.scores.shape


	outFH = open(outF, "w")
	print >> outFH, "\t".join(ms2.header)
	for i in range(ms2.scores.shape[0]):
		if len(list(set(np.where(ms2.scores[i, :] > 0.5)[0]))) > 0:
			print >> outFH, "%s\t%s" % (ms2.IndexToPpi[i], "\t".join(map(str, ms2.scores[i, :])))
	outFH.close()

def exp_comb(args):
	FS, i, j, num_iter, input_dir, num_cores, ref_complexes, scoreF, mode, fun_anno_F, output_dir = args
	i,j, num_iter, num_cores = map(int, [i, j, num_iter, num_cores])


	search_engine = input_dir.split(os.path.sep)[-2]
	def get_eData_comb(data_dir, num_iex, num_beads):
		all_exp =  map(str, glob.glob(data_dir + "*.txt"))
		iex_exp = [f for f in all_exp if (f.split(os.sep)[-1].startswith("all"))]
		beads_exp = [ f for f in all_exp if ( not f.split(os.sep)[-1].startswith("all"))]
		if(i>len(iex_exp)):
			print "i is to large"
			sys.exit()
		if (j > len(beads_exp)):
			print "j is to large"
			sys.exit()

		sel_iex = rnd.sample(iex_exp, num_iex)
		sel_beads = rnd.sample(beads_exp, num_beads)
		return sel_iex + sel_beads


	# EPIC paramters
	if FS == "00000000": sys.exit()
	this_scores = get_fs_comb(FS)
	clf = CS.CLF_Wrapper(num_cores, True)

	ref_gs = Goldstandard_from_cluster_File(ref_complexes)

	scoreCalc = CS.CalculateCoElutionScores(this_scores, "", scoreF, num_cores=num_cores, cutoff=0.5)
	scoreCalc.readTable(scoreF, ref_gs)

	# the supplied functional evidence data needs to have the correct header row...
	functionalData = ""
	if mode == "comb":
		functionalData = utils.get_FA_data("FILE", fun_anno_F)

	if i == 0 and j == 0: sys.exit()

	out_head = ""
	all_scores = []

	for iter in range(num_iter):
		rnd.seed()
		this_eprofiles = get_eData_comb(input_dir, i, j)
		this_eprofiles_fnames = [f.rsplit(os.sep,1)[1] for f in this_eprofiles]
		rnd.seed(1)



		print this_eprofiles_fnames

		this_foundprots, _ = utils.load_data(this_eprofiles, [])
		print len(this_foundprots)

		feature_comb = feature_selector([fs.name for fs in this_scores], scoreCalc, valprots=this_foundprots, elution_file_names=this_eprofiles_fnames)
		if mode == "comb":
			feature_comb.add_fun_anno(functionalData)
		scores, head =  n_fold_cross_validation(5, ref_gs, feature_comb, clf, "%s_%i_%i" % (output_dir, i, j ), overlap = False, local = False)

	#	head, scores = run_epic_with_feature_combinations(this_scores, ref_gs, scoreCalc, clf, output_dir, valprots=this_foundprots)
		out_head = head
		all_scores.append("%s\t%i\t%i\t%s\t%i\t%s" % (FS, i,j,search_engine, len(this_foundprots), scores))
		print head
		print scores


	outFH = open(output_dir + ".%i_%i.all.eval.txt" % (i, j), "w")
	print >> outFH, "FS\tNum_iex\tNum_beads\tSearch_engine\tNum_Prots\t%s" % out_head
	for score in all_scores:
		print >> outFH, "%s" % (score)
	outFH.close()


def EPIC_cor(args):
	print "fuu"
	fs_eval_dir = args[0]
	vals = []
	fs = []
	for fs_eval_F in os.listdir(fs_eval_dir):
		if not fs_eval_F.endswith(".eval.txt"): continue
		fs_eval_F = fs_eval_dir + os.sep + fs_eval_F
		print fs_eval_F
		fs_eval_FH = open(fs_eval_F)
		print fs_eval_FH.readline().split("\t")[1:19]
		for line in fs_eval_FH:
			line = line.rstrip().split("\t")
			if (int(line[2])) < 100: continue
			vals.append(map(float, line[1:19]))
			fs.append(line[0])
		fs_eval_FH.close()
	vals = np.array(vals)
	for row in np.corrcoef(np.transpose(vals)):
		print "\t".join(map("{:.2f}".format, row))

def EPIC_eval_fs(args):
	in_dir, e_dir, scoreF, refF, outF = args
	ref_clusters = GS.Clusters(False)
	ref_clusters.read_file(refF)
	outFH = open(outF, "w")
	i = 0
	allFiles = paths = [os.path.join(in_dir,fn) for fn in next(os.walk(in_dir))[2]]
	for file in allFiles:
		if not file.endswith("clust.txt"): continue
		pred_clusters = GS.Clusters(False)
		pred_clusters.read_file(file)
		_, overlap, _ = pred_clusters.get_matching_complexes(ref_clusters)
		filesplit = file.split(".")[0:4]
		fs_comp = filesplit[0].split(os.sep)[-1]
		scores, head =  utils.clustering_evaluation(ref_clusters, pred_clusters, "Eval", False)
		if i == 0:
			print "FS_code\tCLF\tSE\tFS\tNum_complexes" + "\t".join(np.array(head.split("\t"))[[0,1,6]])
			print >> outFH, "FS_code\tCLF\tSE\tFS\tNum_complexes" + "\t".join(np.array(head.split("\t"))[[0,1,6]])
		print fs_comp + "\t" + "\t".join(filesplit[1:4]) + "\t"+ str(len(pred_clusters.complexes))+"\t" + "\t".join(np.array(scores.split("\t"))[[0, 1, 6]])
		print >> outFH, fs_comp + "\t" + "\t".join(filesplit[1:4]) + "\t"+ str(len(pred_clusters.complexes))+"\t" + "\t".join(np.array(scores.split("\t"))[[0, 1, 6]])
		i += 1
	outFH.close()

def EPIC_eval_fs_DIST(args):

	def getScore(scores):
		out = []
		for cat in [[2], [3], [8]]:
			out.append(sum(scores[cat])/len(cat))
		return out

	def dist(a,b):
		return distance.euclidean(a,b)

	def epic_read_eval(fs_eval_dir):

		def norm_score(scores, columns):
			for i in columns:
				min_score = min(scores[:, i])
				max_score = max(scores[:, i])
				for j in range(len(scores[:, i])):
					scores[j, i] = (scores[j, i] - min_score) / (max_score - min_score)

		vals = []
		fs = []
		header = ""
		for fs_eval_F in  os.listdir(fs_eval_dir):
			if not fs_eval_F.endswith(".eval.txt"): continue
			fs_eval_F = fs_eval_dir + os.sep + fs_eval_F
			fs_eval_FH = open(fs_eval_F)
			param = "-".join(fs_eval_F.split(os.sep)[-1].split(".")[0:2])
			header = np.array(fs_eval_FH.readline().strip().split("\t"))
			header =  np.append(header[1:3], header[11:])

			print fs_eval_F
			for line in fs_eval_FH:
				line = line.rstrip().split("\t")
				scores  = np.append(line[1:3], line[11:])
				vals.append(map(float, scores))
				fs.append("%s-%s" % (param, line[0]))

			fs_eval_FH.close()

		vals = np.array(vals)
		zvals = zscore(vals)

		return header, fs, vals, zvals

	all_scores = {}

	fs_eval_Files, outDir = args

	header, fs, vals, zvals = epic_read_eval(fs_eval_Files)
	fs = np.array(fs)

	def make_hists(x, header, outDir):
		fig = plt.figure()
		plt.rcParams.update({'font.size': 8})
		for i in range(len(header)):
			scores = x[:,i]

			ax = fig.add_subplot(3,4, i+1)
			ax.hist(scores)
			ax.set_title(header[i].replace(" ", "_"), fontsize=6)
			for tick in ax.get_xticklabels():
				tick.set_rotation(45)

		fig.tight_layout()
		fig.subplots_adjust(top=0.88)
		plt.savefig(outDir + ".hist.pdf")
		plt.close()

#	filtering feature selection with low number of predicted clusters

#	sel_vals = set(range(len(vals[:,1])))
#	for k in range(len(header)):
#		this_lb = np.percentile(vals[:, k], 5)
#		this_ub = np.percentile(vals[:, k], 95)
#		sel_vals &= set(np.where(vals[:, k] > this_lb)[0]) & set(np.where(vals[:, k] < this_ub)[0])


	print header[1]
	print header[2]
	print header[3]
	print header[8]

	sel_vals = np.where(vals[:, 1] > 100)[0]
	sel_vals = list(sel_vals)

	vals = vals[sel_vals,]
	zvals = zvals[sel_vals,]
	fs = np.array(fs)[sel_vals,]

	make_hists(vals, header, outDir + ".raw")
	make_hists(zvals, header, outDir + ".zscore")

	this_max_val_fc = []
	this_max_vals = []
	max_zvals = []
	for i in range(len(header)):
		max_index = np.argmax(np.array(zvals[:, i]))
		this_max_vals.append( vals[max_index, i])
		this_max_val_fc.append(fs[max_index])
		max_zvals.append(zvals[max_index, i])
	max_vals = np.array(getScore(np.array(max_zvals)))

#	print max_zvals
#	print max_vals


	composit_scores = {}
	scores = {}
	for i in range(len(fs)):
		this_f = fs[i]
		this_vals = getScore(zvals[i,:])
		this_dist = dist(this_vals,max_vals)
		summed_zscores = sum(this_vals)
		composit_scores[this_f] = summed_zscores
		scores[this_f] = this_dist

	scores_sorted = sorted(scores, key=scores.get)

	outFH = open(outDir + ".results.txt", "w")
	print >> outFH, "Artifical optimal vector"
	print >> outFH, "\t" + "\t".join(header)
	print >> outFH, "Scores:\t\t\t" + "\t".join(map(lambda x : "%.2f" % x, this_max_vals))
	print >> outFH, "Optimal FS per category:\t" + "\t".join(this_max_val_fc)

	print >> outFH, "Composit_scorte\tDistance\tFS\t" + "\t".join(header)
	for  f in scores_sorted:
		score = scores[f]
		c_score = composit_scores[f]
		f_scores = "\t".join(map(lambda x : "%.2f" % x, vals[np.where(fs == f)[0], :][0]))
		print >> outFH, "%.2f\t%.2f\t%s\t%s" % (c_score, score, f, f_scores)
	outFH.close()

def Goldstandard_from_cluster_File(gsF, foundprots = ""):
		clusters = GS.Clusters(need_to_be_mapped=False)
		clusters.read_file(gsF)
		if foundprots != "": clusters.remove_proteins(foundprots)
		gs = GS.Goldstandard_from_Complexes("All")
		gs.complexes = clusters
		gs.make_pos_neg_ppis()
		return gs

class feature_selector:
	def __init__(self, feature_names, scoreCalc, valprots=[],  elution_file_names=[]):
		self.valprots = valprots
		self.get_cols(scoreCalc.header, feature_names, elution_file_names)
		self.cutoff = scoreCalc.cutoff
		self.to_predict = scoreCalc.to_predict
		self.scoreCalc = self.filter_scoreCalc(scoreCalc)
		self.ppiToIndex = self.scoreCalc.ppiToIndex

	def getShape(self):
		return self.scoreCalc.getShape()

	def get_scoreCalc(self):
		return self.scoreCalc

	def set_cutoff(self, cutoff):
		self.cutoff = cutoff

	def get_cols(self, header, feature_names, elution_file_names= []):
		self.to_keep_header = [0, 1]
		self.to_keep_score = []
		for i in range(2, len(header)):
			colname = header[i]
			file_name, scorename = colname.rsplit(".",1)
			if scorename in feature_names and (file_name in elution_file_names or elution_file_names ==[]):
				self.to_keep_header.append(i)
				self.to_keep_score.append(i - 2)

	def valid_score(self, scores):
		return len(list(set(np.where(scores >= self.cutoff)[0])))>0

	def filter_score(self, scores):
		if self.valid_score(scores[self.to_keep_score]):
			return scores[self.to_keep_score]
		else:
			return []

	def filter_scoreCalc(self, scoreCalc):
		filtered_sc = CS.CalculateCoElutionScores("", "", "", 1)
		filtered_sc.scoreF = scoreCalc.scoreF
		filtered_sc.header = list(np.array(scoreCalc.header)[self.to_keep_header])
		filtered_sc.scores = np.zeros((len(scoreCalc.ppiToIndex.keys()),len(self.to_keep_score)))
		ppi_index = 0
		for i in range(scoreCalc.scores.shape[0]):
			ppi = scoreCalc.IndexToPpi[i]
			protA, protB = ppi.split("\t")
			if (protA not in self.valprots or protB not in self.valprots) and self.valprots != []: continue
			ppi_scores = self.filter_score(scoreCalc.scores[i, :])
			if ppi_scores ==[]: continue
			filtered_sc.ppiToIndex[ppi] = ppi_index
			filtered_sc.IndexToPpi[ppi_index] = ppi
			filtered_sc.scores[ppi_index, :] = ppi_scores
			ppi_index += 1
		filtered_sc.scores = filtered_sc.scores[0:ppi_index, :]
		return filtered_sc

	def get_next(self):
		edge, scores = self.scoreCalc.get_next()
		if edge =="":
			#print "recieved empty edge"
			return "", []
		protA, protB = edge.split("\t")

		if (protA not in self.valprots or protB not in self.valprots) and self.valprots != []:
			#print "no elution profile for edge %s\t%s" % (protA, protB)
			return "", []
		out_scores =  self.filter_score(scores)
		if self.scoreCalc.fun_anno != "" and out_scores !=[]:
			to_add = [0] * (len(self.scoreCalc.fun_anno.header) - 2)
			if self.scoreCalc.fun_anno.has_edge(edge):
				to_add = self.scoreCalc.fun_anno.get_score(edge)
			out_scores = np.append(out_scores, to_add)
		return edge, out_scores

	def toSklearnData(self, gs):
		return self.scoreCalc.toSklearnData(gs)

	def open(self):
		self.scoreCalc.open()

	def close(self):
		self.scoreCalc.close()

	def add_fun_anno(self, fun_anno):
		self.scoreCalc.add_fun_anno(fun_anno)

def write_reference(args):
	input_dir, taxid, output_dir = args
	foundprots, elution_datas = utils.load_data(input_dir, [])
	gs = utils.create_goldstandard(taxid, foundprots)
	out = gs.complexes.to_string()
	outFH = open(output_dir, "w")
	print >> outFH, out
	outFH.close()

def calc_feature_combination(args):
	feature_combination, se, input_dir, use_rf, overlap, local, cutoff, num_cores, scoreF, mode, anno ,faF, ref_complexes, output_dir = args
	#Create feature combination
	cutoff = float(cutoff)/100

	if feature_combination == "00000000": sys.exit()
	this_scores = get_fs_comb(feature_combination)
	num_cores = int(num_cores)
	use_rf = use_rf == "True"
	overlap = overlap == "True"
	local = local == "True"

	clf_name = "SVM"
	if use_rf: clf_name = "RF"

	clf = CS.CLF_Wrapper(num_cores, use_rf)

	ref_gs = Goldstandard_from_cluster_File(ref_complexes)


	scoreCalc = CS.CalculateCoElutionScores(this_scores, "", scoreF, num_cores=num_cores, cutoff=cutoff)
	scoreCalc.readTable(scoreF, ref_gs)
	feature_comb = feature_selector([fs.name for fs in this_scores], scoreCalc)
	print feature_comb.scoreCalc.scores.shape
	print scoreCalc.scores.shape
	if mode == "comb":
		if anno == "FILE":
			fa = utils.get_FA_data("FILE", faF)
		if anno == "GM":
			fa = utils.get_FA_data("GM", faF)
		if anno == "STRING":
			fa = utils.get_FA_data("STRING", faF)
			
		feature_comb.add_fun_anno(fa)

	print feature_comb.scoreCalc.scores.shape

	scores, head = n_fold_cross_validation(5, ref_gs, feature_comb, clf, output_dir, overlap, local)


	outFH = open(output_dir + ".eval.txt" , "w")
	print "FS\tSE\tCLF\t" + head
	print "%s\t%s\t%s\t" % (feature_combination, se, clf_name) + scores

	print >> outFH, "FS\tSE\tCLF\t" + head
	print >> outFH, "%s\t%s\t%s\t" % (feature_combination, se, clf_name) + scores
	outFH.close()

def get_fs_comb(comb_string):
	#Create feature combination
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	this_scores = []
	for i, feature_selection in enumerate(comb_string):
		if feature_selection == "1": this_scores.append(scores[i])
	return this_scores

def make_eval(args):
#	pred_clust_F, ref_clust_F, ppiF, cutoff, outF = args

	pred_clust_F, ref_clust_F = args


	#num_ppis = CS.lineCount(ppiF)
	pred_clusters = GS.Clusters(False)
	pred_clusters.read_file(pred_clust_F)

	ref_clusters = GS.Clusters(False)
	ref_clusters.read_file(ref_clust_F)

	#	utils.clustering_evaluation(train.complexes, pred_clusters, "Train", True)
	scores, head = utils.clustering_evaluation(ref_clusters, pred_clusters, "", True)

	#outFH = open(outF, "w")
	#outFH.write("%s\t%i\t%i\t%s\n" % (cutoff, num_ppis, len(pred_clusters.complexes), scores))
	#outFH.close()

def rf_cutoff(args):
	pred_clust_F, ref_clust_F, ppiF, cutoff, outF = args

	num_ppis = CS.lineCount(ppiF)
	pred_clusters = GS.Clusters(False)
	pred_clusters.read_file(pred_clust_F)

	ref_clusters = GS.Clusters(False)
	ref_clusters.read_file(ref_clust_F)

	#	utils.clustering_evaluation(train.complexes, pred_clusters, "Train", True)
	scores, head = utils.clustering_evaluation(ref_clusters, pred_clusters, "", True)

 	outFH = open(outF, "w")
 	outFH.write("%s\t%i\t%i\t%s\n" % (cutoff, num_ppis, len(pred_clusters.complexes), scores))
 	outFH.close()

def orth_map(args):
	clusterF, taxid, outF = args

	clust = GS.Clusters(False)
	clust.read_file(clusterF)

	orthmap = GS.Inparanoid(taxid=taxid)
	orthmap.mapComplexes(clust)

	clust.merge_complexes()
	clust.filter_complexes()

	outFH = open(outF,"w")
	outFH.write(clust.to_string())
	outFH.close()




def calc_scores(args):
	topred = []
	if args[0] == "-ref":
		_, refF, fs, numcores, cutoff, e_dir, outF = args
		gs = Goldstandard_from_cluster_File(refF)
		topred = list(gs.positive | gs.negative)
		print len(topred)
	else:
		fs, numcores, cutoff, e_dir, outF = args



	numcores = int(numcores)
	cutoff = float(cutoff)

	this_fs = get_fs_comb(fs)
	prots, edatas = utils.load_data(e_dir, this_fs)
	scoreCalc = CS.CalculateCoElutionScores(this_fs, edatas, outF, num_cores=numcores, cutoff=cutoff)
	if topred == []: topred = scoreCalc.getAllPairs()
	scoreCalc.calculateScores(topred)

def ppi_fs(args):
	fsc, scoreF, use_rf, se, num_cores, refComplexesF, output_dir = args
	num_cores = int(num_cores)
	use_rf = use_rf == "True"

	clf_name = "SVM"
	if use_rf: clf_name = "RF"
	clf = CS.CLF_Wrapper(num_cores, use_rf)

	this_fs = get_fs_comb(fsc)
	all_gs = Goldstandard_from_cluster_File(refComplexesF)
	valprots = all_gs.get_proteins()

	scoreCalc = CS.CalculateCoElutionScores(this_fs, "", scoreF, num_cores=num_cores, cutoff=-1)
	scoreCalc.readTable(scoreF, all_gs)
	print scoreCalc.scores.shape

	test_scoreCalc = feature_selector([fs.name for fs in this_fs], scoreCalc)

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

	outFH = open("%s.ppi_eva.txt" % (output_dir), "w")
	print >> outFH, "FS\tSE\tCLF\tFM\tauPR\tauROC\n%s\t%s\t%s\t%s" % (fsc, se, clf_name, "\t".join(map(str, evaluation_results)))
	outFH.close()


def main():
	mode = sys.argv[1]

	if mode == "-fs":
		calc_feature_combination(sys.argv[2:])

	elif mode == "-calc_s":
		calc_scores(sys.argv[2:])

	elif mode == "-make_ref":
		write_reference(sys.argv[2:])

	elif mode == "-cor_eval":
		EPIC_cor(sys.argv[2:])

	elif mode == "-best_fs":
		EPIC_eval_fs_DIST(sys.argv[2:])

	elif mode == "-exp_comb":
		exp_comb(sys.argv[2:])

	elif mode == "-merge_ms":
		merge_MS(sys.argv[2:])

	elif mode == "-cut":
		cut(sys.argv[2:])

	elif mode == "-best_fs2":
		EPIC_eval_fs(sys.argv[2:])

	elif mode == "-get_eval":
		make_eval(sys.argv[2:])

	elif mode == "-orthmap":
		orth_map(sys.argv[2:])

	elif mode == "-rfc":
		rf_cutoff(sys.argv[2:])

	elif mode == "-fs_ppi":
		ppi_fs(sys.argv[2:])


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)