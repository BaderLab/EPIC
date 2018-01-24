# EPIC
## Elution Profile-Based Inference of Protein Complexes

### Installation
##### To install EPIC, first make sure you have Python 2.7 and scikit-learn package installed. Also one correlation score ("wcc") utilizes R to perform computation, thus R and rpy2 should be installed in your computer too.

##### First, open terminal and go to your desired directory. Then git clone at your desired directory:
git clone https://github.com/BaderLab/EPIC

### Prepare Input Data
##### Make a file folder and put all the elution profile files into this file folder. There are a couple of examples in the the folder: 
EPIC/test_data/elution_profiles/
##### Uniprot protein identifiers should be used in these elution profile files. Assume the input elution profiles are stored at
../InputFolder/

### Run EPIC
##### Specify correlation scores to be used in EPIC. Eight different correlation socres were implemented in EPIC, in order they are: Mutual Information, Bayes Correlation, Euclidean Distance, Weighted Cross-Correlation, Jaccard Score, PCCN, Pearson Correlation Coefficient, and Apex Score. 
##### "0" indicates that we don't use this correlation score and "1" indicates that we use this correlation score. For example, 11101001 means we will use Mutual Information, Bayes Correlation, Euclidean Distance, Jaccard Score and Apex Score. To specify the correlation scores to use:
-s 11101001

##### There are two ways of generating gold standard protein complexes.
##### The first way is to curate gold standard protein complexes yourself. An exmple file (Worm_reference_complexes.txt) of gold standard protein complexes is given in EPIC/test_data/. The format of gold standard protein complexes file you curated should be same as this file. To specify using manually curated gold standard protein complexes, you can do (assume the protein complexes file is complexes.txt and stored at ../):
-c ../complexes.txt
##### The second way is to automatically download protein complexes from public databases (CORUM, IntAct and GO). In this case, you only need to specify the Taxonomy ID for the species you are studying on using "-t" (assume that we are studying C. elegans):
-t 6239
##### (You should only use one of it, that means, don't use -c and -t at the same time.)

 
##### Specify the directory of storing output files.
##### We need to specify a file folder that can store all the output files. This command needs to be given after the input file folder. Let's assume the output file folder is:
../OutputFolder/

##### Specify the prefix name of output files.
##### You can specify a prefix name for all the output files to distinguish between differnt runs. The defauly is "Out"
-o PrefixName

##### EPIC currently supports two machine learning classifiers: support vector machine and random forest. You can pick one here. "SVM" stands for support vector machine and "RF" stands for random forest. You can specify the machine learning classifier as:
-M RF
##### or
-M SVM

##### You need to specify the number of cores used to run EPIC, the default number is 1. Assume you want to use six cores to run EPIC, you can give the dollowing command:
-n 6

##### EPIC can generate the protein-protein interaction dataset using different resources. As disscussed in the paper, we can use experimental data only, functional evidence data only or the combination of experimental data and functional evidence data. You might want to only use experimental data or in most cases, you might want to use both experimental data and functional evidence data. You can specify which source of data you want to use to generate protein-protein interactions dataset. "EXP" stands for experimental data only, "FA" stands for functional evidence data only and "COMB" stands for using both experimental data and functional evidence data. for exmaple, if you want to use both experimental and functional evidence data, you can use the following command:
-m COMB

##### As mentioned above, EPIC can integrate functional evidence data with experimental data to generate final protein-protein interaction dataset. EPIC can automatically download functional evidence data from GeneMANIA and STRING databases. You can also specify your own curated functional evidence dataset. "GM" stands for GeneMANIA database, "STRING" stands for STRING database and "FILE" stands for the user curated database. For example, you want to use STRING as functional evidence, you can use the following command:
-f STRING
##### If you you use your own curated functional evidence data, you should also specify the directory of functional evidence data file. Assume the functional evidence data file is "fun_anno_file.txt" and it is stored in the file folder "../", you can command as:
-f FILE -F ../fun_anno_file.txt
##### There is an example file gives the format of "fun_anno_file.txt" here, you can find the file at:




### Contact
##### If you have problems, you could contact Lucas (lucasming.huATmail.utoronto.ca) or Florian (florian.goebelsATgooglemail.com).


command = "python /EPIC/src/main.py -s "+ feature + " InputFolder/ -c /Users/lucasminghu/Desktop/WormMap/All_FS_scores/For_selecting_best_features/All.complexes.txt /Users/lucasminghu/Desktop/WormMap/All_FS_scores/For_selecting_best_features/ -o " + feature + mode + anno + " -m " + mode + " -M RF -n 7 -f " + anno  + " -F /Users/lucasminghu/Desktop/WormMap/WormNet/WormNetV3_noZeros_no_physical_interactions.txt -P /Users/lucasminghu/Desktop/WormMap/All_FS_scores/For_selecting_best_features/MSB.11110001.scores.txt

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
