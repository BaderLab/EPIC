# WormMap
The predicted complexes Cytoscape file and elution profile files can be found at https://github.com/BaderLab/EPIC/tree/master/WormMap
(Make sure to update your Cytoscape software version to 3.6.0; elution profile files can be opened in JavaTreeView software.)
ELution files could be found at https://github.com/LucasHu/Worm-co-frac-data

# EPIC
## Elution Profile-Based Inference of Protein Complexes

### Installation
##### To install EPIC, first make sure you have Python 2.7 and scikit-learn package installed. Also one correlation score ("wcc") utilizes R to perform computation, thus R and rpy2 should be installed in your computer too. 
##### We recommend using the Anaconda environemnt to run EPIC and install associated libraries. Anaconda can be downloaded from "https://conda.io/docs/user-guide/install/download.html#anaconda-or-miniconda"
##### Create an Anaconda environment: type "conda create -n EPIC python=2.7 anaconda"
##### If the system says "conda: command not found", you probably need to type "export PATH=~/anaconda2/bin:$PATH" first
##### Activate your Anaconda EPIC environment: type "source activate EPIC"
##### Install "rpy2": type "conda install rpy2"
##### Install "requests": type "conda install requests"
##### Install "scikit-learn": type "conda install scikit-learn"
##### Install "beautifulsoup4": type "conda install beautifulsoup4"
##### Install "mock": type "conda install mock"
##### Install "R": type "conda install -c r r"
##### Install "kohonen": type "conda install -c r r-kohonen"
##### Install "numpy": type "pip install numpy"
##### If this step requires you to install "msgpack" or "argparse", just type "pip install msgpack" and "pip install argparse"
##### The package "wccsom" can be downloaded from https://cran.r-project.org/src/contrib/Archive/wccsom/wccsom_1.2.11.tar.gz
##### If the directory containing this package is "/wccsom_directory", go to r (type: R in command line), type: install.packages('/wccsom_directory/wccsom_1.2.11.tar.gz', type = 'source'). This will install the "wccsom" package.
##### (sometimes, you need to type: install.packages('class') and then type: install.packages('kohonen') in R before doing the step above.)
##### Install "matplotlib": type "python -mpip install -U matplotlib"
##### In case you missed any packages, you usaually can easily type "conda install missing_package_name" to install it...

##### First, open terminal and go to your desired directory. Then git clone at your desired directory:
git clone https://github.com/BaderLab/EPIC

### Prepare Input Data
##### Make a file folder and put all the elution profile files into this file folder. There are a couple of examples in the folder: 
EPIC/test_data/elution_profiles/
##### UniProt protein identifiers should be used in these elution profile files. Assume the input elution profiles (make sure no other folders are in this directory, only store elution profiles files there) are stored at
../InputFolder/

### Run EPIC
##### Specify correlation scores to be used in EPIC. Eight different correlation socres are implemented in EPIC, in order: Mutual Information, Bayes Correlation, Euclidean Distance, Weighted Cross-Correlation, Jaccard Score, PCCN, Pearson Correlation Coefficient, and Apex Score. 
##### "0" indicates that we don't use this correlation score and "1" indicates that we use this correlation score. For example, 11101001 means we will use Mutual Information, Bayes Correlation, Euclidean Distance, Jaccard Score and Apex Score. To specify the correlation scores to use:
-s 11101001

##### There are two ways of generating gold standard protein complexes.
##### The first way is to curate gold standard protein complexes yourself. An example file (Worm_reference_complexes.txt) of gold standard protein complexes is given in EPIC/test_data/. The format of the gold standard protein complex file you curated should be same as this file. To specify using manually curated gold standard protein complexes, you can run the following command (assuming the protein complex file is complexes.txt and stored at ../):
-c ../complexes.txt
##### The second way is to automatically download protein complexes from public databases (CORUM, IntAct and GO). In this case, you only need to specify the Taxonomy ID for the species you are studying on using "-t". For instnace, for C. elegans, we use:
-t 6239
##### (When you use both of the above two options, the reference will come from your input reference complexes file rather than curating from the Internet.)

 
##### Specify the directory of storing output files.
##### We need to specify a file folder that can store all the output files. This command needs to be given after the input file folder. Let's assume the output file folder is:
../OutputFolder/

##### Specify the prefix name of output files.
##### You can specify a prefix name for all the output files to distinguish between different runs. The default is "Out"
-o PrefixName

##### EPIC currently supports two machine learning classifiers: support vector machine and random forest. You can pick one here. "SVM" stands for support vector machine and "RF" stands for random forest. You can specify the machine learning classifier as:
-M RF
##### or
-M SVM

##### You need to specify the number of cores used to run EPIC, the default number is 1. Assume you want to use six cores to run EPIC, you can give the following command:
-n 6

##### EPIC can generate the protein-protein interaction dataset using different resources. As disscussed in the EPIC paper, we can use experimental data only, functional evidence data only or the combination of experimental data and functional evidence data. You might want to only use experimental data or in most cases, you might want to use both experimental data and functional evidence data. You can specify which source of data you want to use to generate protein-protein interactions dataset. "EXP" stands for experimental data only, "FA" stands for functional evidence data only and "COMB" stands for using both experimental data and functional evidence data. for exmaple, if you want to use both experimental and functional evidence data, you can use the following command:
-m COMB

##### As mentioned above, EPIC can integrate functional evidence data with experimental data to generate the final protein-protein interaction dataset. EPIC can automatically download functional evidence data from GeneMANIA and STRING databases. You can also specify your own curated functional evidence dataset. "GM" stands for GeneMANIA database, "STRING" stands for STRING database and "FILE" stands for the user curated database. For example, you want to use STRING as functional evidence, you can use the following command:
-f STRING
##### (Notice that STRING and GEMANIA only have limited number of species, an error will be reported if the species is not in these databases.)
##### If you you use your own curated functional evidence data, you should also specify the directory of the functional evidence data file. Assume the functional evidence data file is "fun_anno_file.txt" and it is stored in the file folder "../", you can use the following parameter:
-f FILE -F ../fun_anno_file.txt
##### There is an example file that shows the format of "fun_anno_file.txt" at:
EPIC/test_data/WormNetV3_noZeros_no_physical_interactions.txt.zip

##### You can also specify the cutoff of the correlation score of co-elution profile, which means if all the correlation socres of a particular protein-protein pair are smaller than the cutoff, this pair will be removed from the dataset. The default threshold is set as 0.5. Assume you want a cutoff at "xx" (this value should be from 0 to 1), you can specify it as:
-r xx

##### You can also set the cutoff for the machine leanring classifier, which means if the predicted protein-protein interaction has a score (assigned by the machine learning classifier) below the cutoff, it will be treated as a false protein-protein interaction. The default is 0.5, which is commonly used for many machine learning classifiers. Assume you want to set a cutoff at "##", you can do:
-R ##

##### For bottom-up proteomics, protein quantification can be determined by MS2 spectral counts. We can set a cutoff to determine the minimal MS2 spectral counts required for each protein. If a particular protein has MS2 spectral counts less than this cutoff across all fractions in a particular co-fractionation experiment, it will be removed from the input dataset. The default number here is 1, and if you want to set the cutoff at "##", you can do:
-e ##

##### If the number of fractions that a particular protein shows up in is less than 1, then the number is likely meaningless for calculating most correlation coefficients. You can set the cutoff for this value. If a proteins identified in n fractions in one co-fractionation experiment and n < cutoff, then this protein will be removed from the input dataset. if you want to set your own cutoff at "##", you can do:
-E ##

##### Calculating different correlation coefficents can be a very time consuming process (taking hours). After computation, a "scores.txt" will be generated and stored in your local computer, assume it is stored at the directory of "../". You can require EPIC to use the pre-calculated "scores.txt" as the input dataset (in which case EPIC will run likely in less than an hour for RF) and do:
-P ../scores.txt

##### Overall, you can run a command like (we take default values for many parameters):
python /EPIC/src/main.py -s 11101001 ../InputFolder/ -c ../complexes.txt ../OutputFolder/ -o PrefixName -M RF -n 6 -m COMB -f STRING 

### Output
##### Once EPIC finishes (which may take many hours), you will find all the output files in the file folder:
../OutputFolder/

### Contact
##### If you have problems, you could contact Lucas (lucasming.huATmail.utoronto.ca) or Florian (florian.goebelsATgooglemail.com).
