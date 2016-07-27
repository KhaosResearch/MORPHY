 -------------------------------
|                               |
|  MORPHY - README              |
|                               |
 -------------------------------

=======================================================================================
TABLE OF CONTENTS
=======================================================================================
1. Requirements 
2. Installing 
3. Executing 
4. Parameters
5. Structure
6. Results
=======================================================================================

=======================================================================================
1. Requirements 
=======================================================================================

MORPHY has been developed in Unix machines (Ubuntu and MacOS X) as well as in
Windows (MinGW) using the G++ (4.4.7) compiler. The make utility has been used to compile the
software package.

MORPHY is based on jMetalCpp and requires install two frameworks:

1) Bio++: a set of C++ libraries for Bioinformatics, including sequence analysis, phylogenetics, molecular evolution and population genetics. Home Page:  http://biopp.univ-montp2.fr/

2) PLL - Phylogenetic Likelihood Library: a highly optimized, parallized software library to ease the development of new software tools dealing with phylogenetic inference. 
Home Page: Phylogenetic Likelihood Library http://www.libpll.org/

Bio++ also requeries the CMake utility and Pll requeries Autoreconf utility.

=======================================================================================

=======================================================================================
2. Installing MORPHY
=======================================================================================

Copy the compressed file to the location where you want to install MORPHY and
unzip it.

Then, run the install.sh with the following command:
  %sh install.sh

It installs all the needed libraries: Bio++ (Bpp-Core, Bpp-Seq and Bpp-Phyl) and Pll.

=======================================================================================


=======================================================================================	
3. Executing MORPHY
=======================================================================================

The main binary is in the subfolder 'bin'. Enter this folder to execute MORPHY
	
  % cd bin
  % ./MORPHY param= parametersfile

The parameterfile is the filename wich has all the parameters requeried to execute the program. For example, to execute the program using one of our test problem 500 sequences dataset, type the following command:
	
  %./MORPHY param=parameters/params500.txt

The parameters to solve the state-of-the-art problems are defined in some files saved into the paramaters folder, and also are included the datasets (nucleotide sequences), partition model file and the initial user trees perfomanced by a bootstrap techniques using PhyML and DNAPars software. 

If you need run a simple example of the algorithm, only define the filename of the dataset of Amino-Acids sequences in Phylip format and the partition file, using the parameter sequencefile and partitionmodelfile

  %./MORPHY sequencefile=sequences.phy partitionmodelfile=sequences.model

All the parameters will be estimated, optimized and defined with the best default values obtanied in our experiments.

=======================================================================================


=======================================================================================
4. Structure 
=======================================================================================

* MORPHY: 
** lib
	*** Bpp: Source code of the Bio++ libraries: Bpp-Core, Bpp-Seq and Bpp-Phyl
	*** Pll: Source code of Phylogenetic Likelihood Library
	libjmetal.a: jMetalCpp library

** src: Source Code of MORPHY based in the jMetalCpp framework. Some changes were added to adjust to the Phylogenetic problem.
	
** bin: Main Binary of the software

** data
	*** sequences: Sequences file
	*** model: Partition Model file used by Pll 
	*** inputusertrees: Input Phylogenetic Trees used as repository to generate initial population

** parameters: Parameters file wich contains all the parameters needed to customize the software.

	
=======================================================================================
5. Parameters 
=======================================================================================

Parameters of the Metaheuristic:

experimentid	 = Experiment ID, the results file are renamed using this ID
DATAPATH   	 = data folder name 

populationsize 	 = Population Sizet
maxevaluations 	 = Number of the evaluations 

intervalupdateparameters = Interval of evaluations to Optimize the Branch-length and parameters of 				   the substitution model  

printtrace 		= Prints trace of the elapsed time and score value performed during 				  optimizing strategies   
printbestscores 	= Prints best scores found during the algorithm execution 

Genetic Operators

selection.method	= binarytournament and randomselection (for SMSEMOA) are Available.

mutation.method 		= NNI, SPR and TRB Topological mutations are available 
mutation.probability 		= Probability to execute

# ----------------------------------------------------------------------------------------
#                                     Phylogenetic Parameters
# ----------------------------------------------------------------------------------------

alphabet	    		= DNA or Protein are available.
sequencefile	 		= The sequence file to use (sequences must be aligned!).
		      	  	For example $(DATAPATH)/sequences/55.txt

input.sequence.format		= The alignment format, for exampĺe
			  	Phylip(order=sequential, type=extended, split=spaces)

partitionmodelfilepll 		= Model specifications file. For example: $(DATAPATH)/model/55.model

init.population 		= Method to generate or create the Initial phylogenetic trees: 
				user, random or stepwise
init.population.tree.file 	= FileName of Input User Trees, for example:
				$(DATAPATH)/inputusertrees/bootstrap1000_55

init.population.tree.format 	= Phylogenetic Trees Format: Newick

# ********** Parameters of the Evolutionary Model  ***********
model 				= Evolutionary Model Name, for example GTR+G

Frequences 

model.frequences 		= user or empirical
model.piA 			=  PiA
model.piC 			=  PiC 
model.piG 			=  PiG 
model.piT 			=  PiT 

GTR relative rate parameters 
model.AC 			= AC
model.AG 			= AG
model.AT 			= AT
model.CG 			= CG
model.CT 			= CT
model.GT 			= GT, always set 1
  
# ********** Distribution Gamma  ***********
rate_distribution 		= Gamma technique only
rate_distribution.ncat		= Number of categories, 4 by default
rate_distribution.alpha 	= Alpha value


# ********** Phylogenetic Optimization  ***********

High-level Strategies for searching the tree space: h2 y h1

H1: This strategy is based on the theoretical perspective of the strong relation between
parsimony and likelihood (minimizing the parsimony score is equivalent to
maximizing the likelihood under some assumptions), and searchs bi-objective phylogenetic trees applying topological moves using the Parametric Progressive Neighbourhood (PPN) technique (Go¨effon et al.,2008) to minimizing the parsimony and with this maximizing the likelihood
simultaneously; furthermore, with the aim of improve the likelihood, optimize all branch lengths affected for the topological moves. 

H2: This strategy is a parametric combination of two techniques that optimize both
criteria separately, for parsimony uses PPN technique applying topological
moves and adjusting the branch lengths affected, and for the likelihood uses a topological rearrange method provided by the Phylogenetic Likelihood Library (PLL).

Parameters:

optimization.method 			= h2 or h2
optimization.method.perc 		= Percentage to apply pllRearrangeSearch or PPN technique

optimization.ppn.numiterations 		= Number of iterations of PPN technique
optimization.ppn.maxsprbestmoves 	= Limit of best moves to be applied in PPN Technique


Parameters of the pllRearrangeSearch  function:

optimization.pll.percnodes 		= % of nodes to be used has root node in pllRearrangeSearch optimization.pll.mintranv		= Define the annulus in which the topological moves should 						take place.
optimization.pll.maxtranv		= Define the annulus in which the topological moves should 						take place.
optimization.pll.newton3sprbranch 	= Optimization of branch-length affected by topological moves


# ********** Branch Length Optimization ***********

Optimization of the Branch-Length of the initial population of phylogenetic trees
bl_optimization.starting 			= true or false
bl_optimization.starting.maxitevaluation	= 1000

Optimization of the Branch-Length during one of the techniques of tree-space exploration or during each Update Interval.
bl_optimization 				= true or false
bl_optimization.maxitevaluation			= 1000 

Optimization of the Branch-Length of the final population of phylogenetic trees perfomanced by algorithm
bl_optimization.final 				= true or false
bl_optimization.final.maxitevaluation		= 1000


# ********** Evolutionary Model Optimization ***********
Optimization of the Parameters of the Evolutionary Model in 
model_optimization 				= true or false
model_optimization.tolerance 			= 0.001

=======================================================================================
6. Results
=======================================================================================

The results of MORPHY are based in the jMetalCpp output format, 
creates two files called FUN+ExperimentID and VAR+ExperimentID. 
The FUN file contains the data of the pareto front approximation and 
the VAR file contains the optimized phylogenetic trees in newick format.

=======================================================================================


7.- News
=======================================================================================


=======================================================================================




