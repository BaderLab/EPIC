# EPIC
## Elution Profile-Based Inference of Protein Complexes

### Installation
#### To install EPIC, first make sure you have Python 2.7 and scikit-learn package installed. Also one correlation score ("wcc") utilizes R to perform computation, thus R and rpy2 should be installed in your computer too.

#### First, open terminal and go to your desired directory. Then git clone at your desired directory:
git clone https://github.com/BaderLab/EPIC

#### Make a file folder and put all the elution profile files into this file folder. There are a couple of examples in the the folder: 
EPIC/test_data/elution_profiles/
#### Uniprot protein identifiers should be used in these elution profile files.




command = "python /Users/lucasminghu/Desktop/newest_EPIC/EPIC/src/main.py -s "+ feature + " /Users/lucasminghu/Desktop/Input_for_EPIC/Data_for_different_engines/MSB/ -c /Users/lucasminghu/Desktop/WormMap/All_FS_scores/For_selecting_best_features/All.complexes.txt /Users/lucasminghu/Desktop/WormMap/All_FS_scores/For_selecting_best_features/ -o " + feature + mode + anno + " -m " + mode + " -M RF -n 7 -f " + anno  + " -F /Users/lucasminghu/Desktop/WormMap/WormNet/WormNetV3_noZeros_no_physical_interactions.txt -P /Users/lucasminghu/Desktop/WormMap/All_FS_scores/For_selecting_best_features/MSB.11110001.scores.txt

	parser.add_argument("-s", "--feature_selection", type = str, help="Select which features to use. This is an 8 position long array of 0 and 1, where each position determines which co-elution feature to use. Features sorted by position are: MI, Bayes, Euclidean, WCC, Jaccard, PCCN, PCC, and Apex.  Each default=11101001", default="11101001")
	parser.add_argument("input_dir",  type = str, help="Directory containing the elution files for each experiment")

	parser.add_argument("-t", "--taxid", type = str, help="TAXID to automatically download reference from GO,CORUM,INtACT",
						default="")
	parser.add_argument("-c", "--cluster", type = str, help="Path to file containing protein clsuter reference",
						default="")
	parser.add_argument("-p", "--ppi", type = str, help="path to ppi File",
						default="")

	parser.add_argument("output_dir", type = str,help="Directory containing the output files")
	parser.add_argument("-o", "--output_prefix", type = str,help="Prefix name for all output Files", default="Out")

	parser.add_argument("-M", "--classifier", type = str,help="Select which classifier to use. Values: RF SVM, default RF",
						default="RF")
	parser.add_argument("-n", "--num_cores", type = int,help="Number of cores to be used, default 1",
						default=1)

	parser.add_argument("-m", "--mode", type=str,
						help="Run EPIC with experimental, functional, or both evidences. Values: EXP, FA, COMB, default: EXP  ",
						default="EXP")
	parser.add_argument("-f", "--fun_anno_source", type = str,help="Where to get functional annotaiton from. Values: STRING or GM or File, default= GM",
						default="GM")
	parser.add_argument("-F", "--fun_anno_file", type=str,
						help="Path to File containing functional annotation. This flag needs to be set when using FILE as fun_anno_source.",
						)
	parser.add_argument("-r", "--co_elution_cutoff", type = float,help="Co-elution score cutoff. default 0.5",
						default=0.5)
	parser.add_argument("-R", "--classifier_cutoff", type = float,help="Classifier confidence valye cutoff. default = 0.5",
						default=0.5)
	parser.add_argument("-e", "--elution_max_count", type = int,help="Removies protein that have a maximal peptide count less than the given value. default = 1",
						default=1)
	parser.add_argument("-E", "--frac_count", type = int,help="Number of fracrions a protein needs to be measured in. default = 2",
						default=2)

	parser.add_argument("-P", "--precalcualted_score_file", type = str,help="Path to precalulated scorefile to read scores from for faster rerunning of EPIC. default = None",
						default="NONE")
