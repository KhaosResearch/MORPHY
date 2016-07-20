/** 
 * MO-PhyTree (version 1.0.0) a software tool for multi-objective 
 * phylogenetic inference. This software integrates features 
 * of the jMetalCpp, Bio++ and PLL frameworks.
 * 
 * Copyright (C) 2016 Cristian Zambrano-Vega, Antonio J. Nebro.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For any other enquiries send an Email to Cristian Zambrano
 * czambrano@uteq.edu.ec
 *
 * When publishing work that uses this software please cite us.
 * 
 * @file Phylogeny.h
 */


#ifndef PHYLOGENY_H_
#define PHYLOGENY_H_


#include <Problem.h>
#include <PhyloTreeSolutionType.h>
#include <Solution.h>
#include <SolutionSet.h>
#include <Ranking.h>
#include <CrowdingComparator.h>

#include <time.h>

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/Io/Fasta.h>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include <Bpp/Phyl/Model/Nucleotide/HKY85.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>

#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Node.h>

#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>


#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Io/FileTools.h>

#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/App/BppApplication.h>


#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/MetaOptimizer.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <pll/pll.h>

#include <Bpp/Text/KeyvalTools.h>

using namespace std;
using namespace bpp;
/**
  * @class Phylogeny
  * @brief Class representing problem Phylogeny
  **/

class Phylogeny : public Problem {
private:
    void PrintScores(double p,double l);
    
public:
  //Phylogeny(string solutionType);
  /*Phylogeny(string SequenceFile, double  kappa_, double alpha_, double beta_, int NumCat_ ,
                    double piA_, double piC_, double piG_, double piT_,
                    int bootstrapSize_,string BootStrapFilename_);*/
  
  Phylogeny(BppApplication *objApp);
  ~Phylogeny();

  void InitializePLLInstance();
  
  vector<string> getLeavesName();
  SiteContainer* getSites();
  SiteContainer* getSites2();
  string getSequenceFilename();
  double getMP();

  void evaluate(Solution *solution);
  void evaluate(Solution *solution,float p, float l);
  
  //void evaluate(Solution *solution, DRTreeParsimonyScore* ParsimonyEvaluator); 
  //void evaluate(Solution *solution, DRTreeParsimonyScore* ParsimonyEvaluator, NNIHomogeneousTreeLikelihood * LikelihoodEvaluator );
  
  double getParsimony(Solution *solution);
      

  void printParameters();
  void readParameters(BppApplication *objApp);
  void GenerateInitialTrees();
  void GenerateRandomTrees();
  void LoadUserTrees();
  boolean PLLisTreeValidate(TreeTemplate<Node> * tree);
  void PLLgenerateParsimonyTrees();
  void ReadSequences(BppApplication *objApp);
  void PrintSequences( pllAlignmentData * alignmentData);

  Alphabet* alphabet;
  IAlignment * seqReader;
  SiteContainer* sites_;
  string modelName;
  string model_frequences;
  DRTreeParsimonyScore* ParsimonyEvaluator ;
  Newick * newick;

 /********** PLL Data Objects ****/
  pllAlignmentData * alignmentData;
  pllQueue * partitionInfo;
  pllInstanceAttr attr;
  pllInstance * tr;
  partitionList * partitions;
  string rateheterogeneity_model;
  /*******************************/
  
  string SequenceFilename;
  map<int, double> neuclotidesfreqs;
  int numarbol;
  vector<Tree*> treesin;
  vector<Tree*> trees;
  int bootstrapSize;
  string bootstrapFilename;
  string initialTrees;
  string FilenameIntree;
  double mp;
  double ml;
  string PartitionModelFilePLL;

  double kappa, alpha, piA, piC, piG, piT,AC,AG,AT,CG,CT,GT;
  
  int NumPreliminarResults;
    
  size_t NumCat;
  
  int NumEvaluaciones;
  
  
  void OpenScores();
  void CloseScores();
  ofstream ComportamientoML,ComportamientoTime, ComportamientoMP , ComportamientoTimePar, ScoresML, ScoresMP;
  
  bool printtrace, printbestscores;
  
  bool printsequence;
  string  sequenceFormat;
  string  alphabetType;
 
 
  
  
  //Optimization 
  
  pair<double,boolean> PLLSearch(double iter);
  int PLLgetDeep(pllInstance * tr, nodeptr p);
  int PLLgetMaxDeep(pllInstance * tr, nodeptr p);
  double PLLOptimizarRamas2(Solution *  solution, int numevals);
  pair<double, string> PLLOptimizarRamas(Solution *  solution, int numevals);
  
  double PLLOptimizarParamModeloSust(double epsilon);
  
  double OptimizarParamModeloSust(Solution * solution);
  
  pair< pair<double *, int *>, Node *> SPR(Node * Nodo1, Node * Nodo2, int &NextIDNode);
  void SPR(Node * Nodo1, Node * Nodo2, double * b);
  //void ChangeBranchLenghth(Node * Nodo, double length, NNIHomogeneousTreeLikelihood * &NNILik  );
  int SPRvalide (Node* N1, Node* N2);
  void SPRreverse(Node * Nodo1, Node * Nodo2, double * b, int &NextIDNode);
  int getNivel(Node* nodo);
  bool isMov(vector< pair<int , int > > vIDS, int IDN1, int IDN2);
  //void PPNOptimiz(Solution * solution);
  pair<double,boolean> PPNSearch(Solution * solution,int NumMaxMovs, int PPN_NumIteraciones);
  
  void Optimization(Solution * solution);
  void OptimizationFull(Solution * solution,int iter);
  void OptimizacionHibrida(Solution * solution, int iter);
  void OptimizacionHibrida2(Solution * solution, int iter);
  
  void BranchLengthOptimization(Solution *  solution, string method, int NumEvals , double Tolerance);
  
  double BppGradientBLOptimization(Solution *  solution, int NumEvals , double Tolerancia);
      
  
  SolutionSet * generateInitialPopulation(int populationSize);
  void setPLLToTree(Solution *solution);
  void setTreeToPLL(Solution *solution, boolean hasToEvaluate) ;
  
  void PLLNewtonBLOptimization(int numevals);
  
  void updateInterval(SolutionSet * population);
  void printPreliminarResults(int NumExp, SolutionSet* population);
  
  void setPopulationNotNew(SolutionSet * population);
  int getNewSolutionsCount(SolutionSet * population);
  
  string OptimizationMethod;
  double porcentajecriterios;
  
  int PPN_NumIteraciones,PPN_MaxTopMoves;
  double pll_percnodes ,pll_mintranv,pll_maxtranv;
  bool pll_newton3sprbranch;
  
  
  
  bool StartingOptRamas;
  int StartingNumIterOptRamas;
  bool FinalOptRamas;
  int FinalNumIterOptRamas;
  
  bool OptimizationOptRamas;
  int OptimizationNumIterOptRamas;
  
  
  bool OptimizacionSubstModel;
  double TolerenciaOptSubstModel;
          

};

#endif /* PHYLOGENY_H_ */
