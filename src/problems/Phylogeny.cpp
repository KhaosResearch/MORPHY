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
 * @file Phylogeny.cpp
 */

#include <Phylogeny.h>
#include "vector"



/**
 * Class constructor
 */

Phylogeny::Phylogeny(BppApplication *objApp) {

        numberOfVariables_   = 1;
	numberOfObjectives_  = 2;
	numberOfConstraints_ = 0;
	problemName_="Phylogenetic Problem";
        solutionType_ = new PhyloTreeSolutionType(this);

	ParsimonyEvaluator=NULL;
        numarbol=0; 
        mp=0; ml=0;
        NumEvaluaciones =0;

        readParameters(objApp); //From Console or FileParameters
        
         /* Set the PLL instance attributes */
        if(rateheterogeneity_model=="GAMMA")
            attr.rateHetModel     = PLL_GAMMA;
        else
            attr.rateHetModel     = PLL_CAT;
        
        attr.fastScaling      = PLL_FALSE;
        attr.saveMemory       = PLL_TRUE;
        attr.useRecom         = PLL_FALSE;
        attr.randomNumberSeed = 0xDEADBEEF;
        attr.numberOfThreads  = 8; 


       ReadSequences(objApp);
       
       if (printsequence){         PrintSequences(alignmentData);       }
       
        partitionInfo = pllPartitionParse (PartitionModelFilePLL.c_str());
        if (!pllPartitionsValidate (partitionInfo, alignmentData)){
                cout << "Error: Partitions do not cover all sites in " << SequenceFilename << endl;	exit(-1);
        }
       
       InitializePLLInstance();
         
       newick = new Newick;
       GenerateInitialTrees();
        
       if(printtrace) OpenScores();
       
       
} // Phylogeny


void Phylogeny::InitializePLLInstance(){
    
/**Initialize PLL Instance and Partitons */ 
  tr = pllCreateInstance (&attr);
  partitions = pllPartitionsCommit (partitionInfo, alignmentData);
  
  pllTreeInitTopologyRandom (tr, alignmentData->sequenceCount, alignmentData->sequenceLabels);
  
  if (!pllLoadAlignment (tr, alignmentData, partitions)) {  cout << "InitializePLLInstance, Incompatible tree/alignment combination" << endl;	exit(-1); }
 
  pllInitParsimonyStructures (tr, partitions, PLL_TRUE);
  
  pllInitModel(tr, partitions);

  pllSetFixedAlpha(alpha, 0, partitions, tr); 
  
  if(model_frequences=="user"){
        double q[6] = {AC,AG,AT,CG,CT,GT}; 
        double f[4] = {piA,piC,piG,piT}; 
        pllSetFixedSubstitutionMatrix(q, 6, 0, partitions, tr);
        pllSetFixedBaseFrequencies(f, 4, 0, partitions, tr);
  }

  pllEvaluateParsimony(tr, partitions, tr->start, PLL_TRUE, PLL_TRUE);
  pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  
  
    
}

void Phylogeny::PrintSequences(pllAlignmentData * alignmentData){
    cout << "SEQUENCES" << endl;  
    cout << "Number of Sequences  " <<  alignmentData->sequenceCount << " Sequence Length " << alignmentData->sequenceLength << endl;
    
    //printf ("getNumberOfSequences %d\n", sites_->getNumberOfSequences());
       for (int i = 1; i <= alignmentData->sequenceCount; ++ i){
            cout << "Sequence N° " << i << " " << alignmentData->sequenceLabels[i] << endl << alignmentData->sequenceData[i] << endl << endl;
            //cout   << sites_->getSequence(i-1).getName()  << endl << sites_->getSequence(i-1).toString() << endl << endl;
        }
        
}

void Phylogeny::ReadSequences(BppApplication * objApp){
    
        alphabetType = ApplicationTools::getStringParameter("alphabet", objApp->getParams(), "DNA", "", false,false);
            
        if(alphabetType=="DNA")
            alphabet = new DNA();
        else
            alphabet = new ProteicAlphabet();
        
        
        sequenceFormat = ApplicationTools::getStringParameter("input.sequence.format", objApp->getParams(), "Phylip(order=interleaved, type=extended, split=spaces)", "", false,false);
         string format = "";
        std::map<std::string, std::string> Arguments;
        KeyvalTools::parseProcedure(sequenceFormat, format, Arguments);
        
        int SequenceTypePLL;
        
        if(format=="Phylip"){
            bool sequential = true, extended = true;
            string split = "  ";
            
            string order = ApplicationTools::getStringParameter("order", Arguments, "sequential", "", false, false);
            sequential = (order == "sequential")?true:false; //sequential or interleaved

            string type = ApplicationTools::getStringParameter("type", Arguments, "extended", "", false, false);
            if (type == "extended"){  // extended OR classic
                extended = true;
                split = ApplicationTools::getStringParameter("split", Arguments, "spaces", "", false, false);
                if (split == "spaces")
                  split = " ";
               else if (split == "tab")
                split = "\t";
           }else  extended = false;
   
            SequenceTypePLL = PLL_FORMAT_PHYLIP;
                    
            seqReader = new Phylip(extended,sequential,100,true,split);
            
        }else{ //Fasta
            
            int charsByLine = 100;
            bool checkSequenceNames = true, extended = false, strictSequenceNames = false;
            
            charsByLine = ApplicationTools::getIntParameter("charsbyline", Arguments, 100, "", false, false);

            checkSequenceNames = ApplicationTools::getBooleanParameter("checksequencenames", Arguments, true, "", false, false);
            extended = ApplicationTools::getBooleanParameter("extended", Arguments, false, "", false, false);
            strictSequenceNames = ApplicationTools::getBooleanParameter("strictsequencenames", Arguments, false, "", false, false);
            
            SequenceTypePLL = PLL_FORMAT_FASTA;
            
             seqReader = new Fasta(charsByLine, checkSequenceNames, extended, strictSequenceNames);
        }
       
       sites_ = seqReader->readAlignment(SequenceFilename, alphabet);
       
      //PLL Sequences Format
       alignmentData = pllParseAlignmentFile (SequenceTypePLL,SequenceFilename.c_str());
       if (!alignmentData){
           cout << "PLL Error while parsing " << SequenceFilename << endl;	exit(-1);
       }
        
       
}

void Phylogeny::GenerateInitialTrees(){
    
    if(initialTrees =="user")
        
        LoadUserTrees();
    
    else if(initialTrees =="stepwise")
        
        PLLgenerateParsimonyTrees();
        
   else  //Default Random

        GenerateRandomTrees();
}

void Phylogeny::LoadUserTrees(){
     if(!FileTools::fileExists(bootstrapFilename)){
            cout << "Can't open file of Initial Trees " << bootstrapFilename << endl;		exit(-1);
        } 
        
        IMultiTree* w = new Newick(true);
        w->read(bootstrapFilename, treesin);
        delete w;

        cout << "BootStrap Trees: " << treesin.size() << endl;
        
        if(treesin.size()<bootstrapSize){
            cout << "BootStrap size " << treesin.size() << " is lower than Population Size " << bootstrapSize << endl;		
            exit(-1);
        }
        
        int pos,i;
        TreeTemplate<Node> * tree;
        for(i=1;i<=bootstrapSize;i++){
            do{
                  pos=PseudoRandom::randInt(0, treesin.size() - 1);
                  tree = new TreeTemplate<Node>(*treesin[pos]);
            }while(!PLLisTreeValidate(tree));
            
            trees.push_back(tree);
            delete treesin[pos];
            treesin.erase(treesin.begin()+pos);
            
       }
           
        for(i=0;i<treesin.size();i++)  delete treesin[i]; 
        treesin.clear();
                
        cout << "Initial BootStrap Trees done!!! Size= " << trees.size() << endl;
}

void Phylogeny::GenerateRandomTrees(){
    
     TreeTemplate<Node> * tree;
  
     for(int k=1;k<=bootstrapSize;k++){
            vector<string> Leaves=getLeavesName();
            tree = TreeTemplateTools::getRandomTree(Leaves, false);
	    tree->setBranchLengths(0.05);
            trees.push_back(tree);
     }
     
     cout << "Random Trees done!!! Size= " << trees.size() << endl;
}


void Phylogeny::PLLgenerateParsimonyTrees(){
  
  pllInstance * tr;
  partitionList * partitions;
  
  tr = pllCreateInstance (&attr);
   
  partitions = pllPartitionsCommit (partitionInfo, alignmentData);
  
  pllTreeInitTopologyRandom (tr, alignmentData->sequenceCount, alignmentData->sequenceLabels); 
  
  if (!pllLoadAlignment (tr, alignmentData, partitions)) {  cout << "Evaluating, Incompatible tree/alignment combination" << endl;	exit(-1); }

  TreeTemplate<Node> * tree;
  
  for(int k=1;k<=bootstrapSize;k++){
      
    pllInitParsimonyStructures (tr, partitions, PLL_TRUE);  
    tr->randomNumberSeed = k*1234567;

    pllMakeParsimonyTreeFast(tr, partitions);
    cout << "Generating Parsimony Tree " << k << endl;
    //printf ("Score Parsimony %d : %d\n",k ,pllEvaluateParsimony(tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE) );

    pllTreeToNewick(tr->tree_string, tr, partitions, tr->start->back, false, true, false, false, false, 1, false, false);
    
    string arbol(tr->tree_string);   //fprintf(stderr, "%s\n", tr->tree_string );
    tree = TreeTemplateTools::parenthesisToTree(arbol);
    tree->setBranchLengths(0.05);
    
    trees.push_back(tree);
  }

 cout << "Initial Parsimony Trees ggenerated by StepWise Addtion done!!! Size= " << trees.size() << endl;
  
 //pllPartitionsDestroy (tr,&partitions);
 pllDestroyInstance (tr);
 
}

boolean Phylogeny::PLLisTreeValidate(TreeTemplate<Node> * tree){

  pllNewickTree * newick;
  boolean res=true;
  string treenewick = TreeTemplateTools::treeToParenthesis(*tree) ;
  newick = pllNewickParseString	(treenewick.c_str());
  
  if (!newick){
      res = false;
  }else {
      if (!pllValidateNewick (newick))  res=false;
  }
  pllNewickParseDestroy (&newick);
  return res;
}

void Phylogeny::readParameters(BppApplication *objApp){

    
    rateheterogeneity_model = ApplicationTools::getStringParameter("rateheterogeneity_model", objApp->getParams(), "GAMMA", "", false, false);
    
    //Sequences input File
    SequenceFilename = ApplicationTools::getStringParameter("sequencefile", objApp->getParams(), "55.txt", "", false, false);
    if(!FileTools::fileExists(SequenceFilename)){
            cout << "Can't open Sequences file " << SequenceFilename << endl;		exit(-1);
    }

    printtrace = ApplicationTools::getBooleanParameter("printtrace", objApp->getParams(), false, "", false, false);
    printbestscores = ApplicationTools::getBooleanParameter("printbestscores", objApp->getParams(), false, "", false, false);
    
    NumPreliminarResults = ApplicationTools::getIntParameter("NumPreliminarResults", objApp->getParams(), 0, "", false, false);
    
        
    printsequence = ApplicationTools::getBooleanParameter("input.sequence.print", objApp->getParams(), false, "", false, false);
     
   // Substitution Model Parameters
    modelName = ApplicationTools::getStringParameter("model", objApp->getParams(), "GTR", "", false, false); 

    model_frequences = ApplicationTools::getStringParameter("model.frequences", objApp->getParams(), "none", "", false, false); 
    if(model_frequences=="user"){
        piA =   ApplicationTools::getDoubleParameter("model.piA", objApp->getParams(), 0.25, "", false, false); 
        piC =   ApplicationTools::getDoubleParameter("model.piC", objApp->getParams(), 0.25, "", false, false); 
        piG =   ApplicationTools::getDoubleParameter("model.piG", objApp->getParams(), 0.25, "", false, false); 
        piT =   ApplicationTools::getDoubleParameter("model.piT", objApp->getParams(), 0.25, "", false, false); 
        AC =   ApplicationTools::getDoubleParameter("model.AC", objApp->getParams(), 1, "", false, false); 
        AG =   ApplicationTools::getDoubleParameter("model.AG", objApp->getParams(), 1, "", false, false); 
        AT =   ApplicationTools::getDoubleParameter("model.AT", objApp->getParams(), 1, "", false, false); 
        CG =   ApplicationTools::getDoubleParameter("model.CG", objApp->getParams(), 1, "", false, false); 
        CT =   ApplicationTools::getDoubleParameter("model.CT", objApp->getParams(), 1, "", false, false); 
        GT =   ApplicationTools::getDoubleParameter("model.GT", objApp->getParams(), 1, "", false, false); 
    }

    // Gamma Distribution Parameters
    NumCat = ApplicationTools::getDoubleParameter("rate_distribution.ncat", objApp->getParams(), 4, "", false, false); 
    alpha =  ApplicationTools::getDoubleParameter("rate_distribution.alpha", objApp->getParams(), 0.5, "", false, false); 

    //BootStrap Size mous to be equal to Population size
    bootstrapSize = ApplicationTools::getIntParameter("bootstrapSize", objApp->getParams(), 100, "", false, false);
    initialTrees = ApplicationTools::getStringParameter("init.population", objApp->getParams(), "stepwise", "", false, false);
    
    if(initialTrees=="user")
        bootstrapFilename = ApplicationTools::getStringParameter("init.population.tree.file", objApp->getParams(), "userinittree.txt", "", false, false);
    

    PartitionModelFilePLL = ApplicationTools::getStringParameter("partitionmodelfile", objApp->getParams(), "55.model", "", false, false);
  
    OptimizationMethod = ApplicationTools::getStringParameter("optimization.method", objApp->getParams(), "h2", "", false, false);
    if(OptimizationMethod=="h2")
        porcentajecriterios = ApplicationTools::getDoubleParameter("optimization.method.perc", objApp->getParams(), 0.5, "", false, false);

    //PPN Parameters
    PPN_NumIteraciones   = ApplicationTools::getDoubleParameter("optimization.ppn.numiterations", objApp->getParams(), 100, "", false, false);
    PPN_MaxTopMoves      = ApplicationTools::getDoubleParameter("optimization.ppn.maxsprbestmoves", objApp->getParams(), 20, "", false, false);
    
    //PLL Parameters
    pll_percnodes        = ApplicationTools::getDoubleParameter("optimization.pll.percnodes", objApp->getParams(), 0.2, "", false, false);
    pll_mintranv         = ApplicationTools::getDoubleParameter("optimization.pll.mintranv", objApp->getParams(), 1, "", false, false);
    pll_maxtranv         = ApplicationTools::getDoubleParameter("optimization.pll.maxtranv", objApp->getParams(), 20, "", false, false);
    pll_newton3sprbranch = ApplicationTools::getBooleanParameter("optimization.pll.newton3sprbranch", objApp->getParams(), true, "", false, false);
    
    //Branch Length Optimization
    StartingOptRamas           = ApplicationTools::getBooleanParameter("bl_optimization.starting", objApp->getParams(), true, "", false, false);
    StartingNumIterOptRamas    = ApplicationTools::getDoubleParameter("bl_optimization.starting.maxitevaluation", objApp->getParams(), 1000, "", false, false);
    
    FinalOptRamas           = ApplicationTools::getBooleanParameter("bl_optimization.final", objApp->getParams(), true, "", false, false);
    FinalNumIterOptRamas    = ApplicationTools::getDoubleParameter("bl_optimization.final.maxitevaluation", objApp->getParams(), 1000, "", false, false);

    OptimizationOptRamas           = ApplicationTools::getBooleanParameter("bl_optimization", objApp->getParams(), false, "", false, false);
    OptimizationNumIterOptRamas    = ApplicationTools::getDoubleParameter("bl_optimization.maxitevaluation", objApp->getParams(), 1000, "", false, false);

    //Substitution Model Optimization
    OptimizacionSubstModel           = ApplicationTools::getBooleanParameter("model_optimization", objApp->getParams(), true, "", false, false);
    TolerenciaOptSubstModel          = ApplicationTools::getDoubleParameter("model_optimization.tolerance", objApp->getParams(), 0.001, "", false, false);
   
}



void Phylogeny::printParameters(){
    
        cout <<"************** PARAMETERS OF THE PHYLOGENETIC PROBLEM  ***********************" << endl;
        cout << "Problem: " << getName() << endl;
        cout << "Sequence Filename: " << SequenceFilename << endl;
        cout << "Input Sequence Format " << sequenceFormat << endl;
	//cout << "Number of sequences Sites: " << sites_->getNumberOfSequences() << " Sites2: " << sites2_->getNumberOfSequences() << endl;
	//cout << "Number of sites: " << sites_->getNumberOfSites() << " Sites 2: " << sites2_->getNumberOfSites()  << endl;

        cout << "Number of sequences: " << alignmentData->sequenceCount << endl;
	cout << "Number of sites (sequenceLength): " << alignmentData->sequenceLength   << endl;
        
        cout << "Print input Sequences " << printsequence << endl;
        
        
        cout << "InitialTrees Mode: " << initialTrees << endl;
        if(initialTrees=="user")
                cout << "BootstrapFilename: " << bootstrapFilename << endl;
        cout << "BootstrapSize: " << bootstrapSize << endl;
        
        
        cout << "Rate Heterogeneity Model: " << rateheterogeneity_model << endl;  
        cout << "ModelName: " << modelName << endl;  
        //cout << "Transition / transversion ratio (Kappa): " << kappa << endl; 
        cout << "Num Categories: " << NumCat << endl; 
        cout << "Gamma Shape (alpha): " << alpha << endl; 
        cout << "Model Frequences: " << model_frequences << endl; 
        
        if(model_frequences=="user"){
                cout << "Frecuencies  piA: " << piA << " piC: " << piC << " piG: " << piG << " piT: " << piT << endl;  
                cout << "Relative rates AC: " << AC << " AG: " << AG << " AT: " << AT << " CG: " << CG << " CT: " << CT << " GT: " << GT << endl ;  
        }

        cout << "PartitionModelFilePLL: " << PartitionModelFilePLL << endl << endl;  
        
        
        cout << "************** Optimization Parameters ************ " << endl;
        
        cout << "Method: " <<  OptimizationMethod << endl;
        if(OptimizationMethod=="h2"){
         cout << "Criteria Percentaje: " <<  porcentajecriterios << endl;
        }
         
        cout << "PPN Number of Iterations " << PPN_NumIteraciones << endl;
        cout << "PPN_MaxTopMoves " << PPN_MaxTopMoves << endl;
        
        cout << "pll_percnodes " << pll_percnodes << endl;
        cout << "pll_mintranv " << pll_mintranv << endl;
        cout << "pll_maxtranv " << pll_maxtranv << endl;
        cout << "pll_newton3sprbranch " << pll_newton3sprbranch << endl;
        
        
        cout << "Starting Branch Length Optimization:  " <<  StartingOptRamas << endl;
        if(StartingOptRamas){
           cout << "    Numeric Method: Newton " <<  endl;
           cout << "    Max Iterations: " <<  StartingNumIterOptRamas << endl;
         }
        
        cout << "Final Branch Length Optimization:  " <<  FinalOptRamas << endl;
        if(FinalOptRamas){
           cout << "    Numeric Method: Newton" << endl;
           cout << "    Max Iterations: " <<  FinalNumIterOptRamas << endl;
         }
        
        cout << "Branch Length Optimization:  " <<  OptimizationOptRamas << endl;
        if(OptimizationOptRamas){
           cout << "    Numeric Method: Newton" <<  endl;
           cout << "    Max Iterations: " <<  OptimizationNumIterOptRamas << endl;
         }
        
        cout << "Optimization of Parameters of the Substitution Model:  " <<  OptimizacionSubstModel << endl;
        if(OptimizacionSubstModel){
                cout << "    Method: Brent" << endl;
                cout << "    Tolerance : " <<  TolerenciaOptSubstModel << endl;
        }
   
        cout << "Print Trace : " << printtrace << endl;
        cout << "Print Best Scores : " << printbestscores << endl;
        
        cout << "NumPreliminarResults: " << NumPreliminarResults << endl;
        
        
        
        cout << endl;
        
}

void Phylogeny::OpenScores(){
        ComportamientoML.open("TraceML.txt"); ComportamientoTime.open("TraceTime.txt");
        ComportamientoMP.open("TraceMP.txt"); ComportamientoTimePar.open("TraceTimePar.txt");
        ScoresMP.open("ScoresMP.txt"); ScoresML.open("ScoresML.txt");
}



void Phylogeny::CloseScores(){
      ComportamientoML.close();	ComportamientoTime.close();
      ComportamientoMP.close();	ComportamientoTimePar.close();
      ScoresMP.close();	ScoresML.close();
}

/**
 * Destructor
 */
Phylogeny::~Phylogeny() {

  delete seqReader;
  delete alphabet;
  delete sites_ ;
  delete newick;

  pllAlignmentDataDestroy (alignmentData);
  pllQueuePartitionsDestroy (&partitionInfo);
  pllPartitionsDestroy (tr,&partitions);
  pllDestroyInstance (tr);
  
  
  delete solutionType_;

   if(printtrace)  CloseScores();
   

} // ~Phylogeny


/**
  * Evaluates a solution
  */

void Phylogeny::evaluate(Solution *solution) {

//  double p,l;  p=0;l=0;
//  
//  NumEvaluaciones++;
//  
//  PhyloTree * Pt = (PhyloTree*) solution->getDecisionVariables()[0];
//  
//  l = getLikelihood(solution);  Pt->setLnLikelihoodValue(l);
//
//  p= getParsimony(solution);  Pt->setParsimonyScore(p);
// 
//  PrintScores(p,l); //Compare negative values
//
//  solution->setObjective(0,p);
//  solution->setObjective(1,l*-1); //The objective values are saved as Positive (jMetal format), but they are negative

} // evaluate

void Phylogeny::evaluate(Solution *solution, float p, float l) {
   
  NumEvaluaciones++;
   
  PhyloTree * Pt  = ((PhyloTree*)solution->getDecisionVariables()[0]);
  
  Pt->setParsimonyScore(p);  Pt->setLnLikelihoodValue(l);

  PrintScores(p,l);  //Compare negative values

  solution->setObjective(0,p);
  solution->setObjective(1,l*-1); //The objective values are saved as Positive (jMetal format), but they are negative

} // evaluate


double Phylogeny::getParsimony(Solution *solution){
  
  PhyloTree * Pt  = ((PhyloTree*)solution->getDecisionVariables()[0]);
  Tree * tree = Pt->getTree();
  
  ParsimonyEvaluator =	new DRTreeParsimonyScore(*tree, *sites_,false,true);
  double p = (double)ParsimonyEvaluator->getScore();
  delete ParsimonyEvaluator;
  return p;
}

void Phylogeny::PrintScores(double p,double l){
   
        if(ml==0) {   ml=l;  
                      if(printbestscores) cout << "First Likelihood " << ml  << endl; 
                      if(printtrace)  ScoresML << NumEvaluaciones << " " << l << endl; 
        }
        if(l>ml)  { 
                     if(printbestscores) cout << "Current Likelihood: "  << ml << " :: Now Best Likelihood: " << l << endl;
                     ml=l;  
                    if(printtrace)  ScoresML << NumEvaluaciones << " " << l << endl;
        }
        if(mp==0) {   
                     mp=p;  
                     if(printbestscores)  cout << "First Parsimony "  << mp << endl; 
                     if(printtrace)  ScoresMP << NumEvaluaciones << " " << p << endl;
        }
        if(p<mp)  {  
                    if(printbestscores)  cout << "Current Parsimony " << mp << " :: Now Best Parsimonia= " << p << endl;
                    mp=p;  
                    if(printtrace)  ScoresMP << NumEvaluaciones << " " << p << endl;
        }
}

vector<string> Phylogeny::getLeavesName(){
	return sites_->getSequencesNames();
}

SiteContainer* Phylogeny::getSites(){
	return sites_;
}

string Phylogeny::getSequenceFilename(){
	return SequenceFilename;
}


double Phylogeny::getMP(){ return mp;}

void Phylogeny::Optimization(Solution *solution) {
    
    if(OptimizationMethod=="h1")        
       OptimizacionHibrida(solution,1);
    else if(OptimizationMethod=="h2")   
       OptimizacionHibrida2(solution,1);
    
}

void Phylogeny::OptimizationFull(Solution *solution, int num) {
    
    if(OptimizationMethod=="h1")        
       OptimizacionHibrida(solution,num);
    else if(OptimizationMethod=="h2")   
       OptimizacionHibrida2(solution,num);
    
}

void Phylogeny::OptimizacionHibrida(Solution * solution, int iter){
    //cout <<"============================================================" << endl;
    //cout << "Parsimonia I " <<  problem_->getParsimony(solution)  << "  Likelihood I " << problem_->getLikelihood(solution) << endl;
    
    pair<double , boolean> Par; 
    
    Par = PPNSearch(solution,PPN_MaxTopMoves*iter, PPN_NumIteraciones * iter);
    
    setTreeToPLL(solution,true);
    
    if(Par.second){  
       if(OptimizationOptRamas) {
            pllOptimizeBranchLengths(tr, partitions, OptimizationNumIterOptRamas*iter);
            setPLLToTree(solution);
       }
    }
    
   evaluate(solution,(double)Par.first,(double)tr->likelihood); //Only to update scores in log files 
  
 //cout << "Parsimonia F " <<  Par.first << "  Likelihood F " << Lk  << endl;
 //cout <<"============================================================" << endl; 

}

void Phylogeny::OptimizacionHibrida2(Solution * solution, int iter){
  
    double Par;
    pair<double , boolean> Res; 
    
    //cout <<"============================================================" << endl;
   //cout << "Parsimonia I " <<  problem_->getParsimony(solution)  << "  Likelihood I " << problem_->getLikelihood(solution) << endl;

    if (PseudoRandom::randDouble() <= porcentajecriterios) {
        
        OptimizacionHibrida(solution, iter);
                
    }else {
        
         setTreeToPLL(solution,true);
         Res = PLLSearch((double)iter);
         
         if(Res.second){ //There is a topology change
             if(OptimizationOptRamas){ 
                  pllOptimizeBranchLengths(tr, partitions, OptimizationNumIterOptRamas*iter);
             }
            setPLLToTree(solution);
            Par =  getParsimony(solution);
            
            evaluate(solution,(double)Par,(double)tr->likelihood);   //Only to update scores in log files 
         }
         
         
    }
    
 
  
 //cout << "Parsimonia F " <<  Par << "  Likelihood F " << Lk  << endl;
 //cout <<"============================================================" << endl; 

}

boolean esHoja(int number, int maxTips)
{
  if(number <= maxTips)
    return PLL_TRUE;
  else
    return PLL_FALSE;
}

int Phylogeny::PLLgetMaxDeep(pllInstance * tr, nodeptr p){
  int d1,d2;
  
  d1= PLLgetDeep(tr,p);
  d2= PLLgetDeep(tr,p->back);
  if(d1>d2)  return d1; else return d2;
}

int Phylogeny::PLLgetDeep(pllInstance * tr, nodeptr p){
  nodeptr  p1, p2;
  int d1,d2;
  if (!esHoja (p->number, tr->mxtips)) {
     d1= PLLgetDeep(tr,p->next->back);
     d2= PLLgetDeep(tr,p->next->next->back);
     if(d1>d2)
         return d1+1;
     else 
         return d2+1;
 }else
     return 1;
     
}

pair<double,boolean>  Phylogeny::PLLSearch(double iter){
  
  pair<double,boolean> R;
  pllRearrangeList * rearrangeList;

  clock_t t_ini = clock();

 double lk1, lk2;
 lk1 = lk2 = R.first= tr->likelihood;
 R.second = false;
 
 tr->thoroughInsertion = pll_newton3sprbranch;
  
  nodeptr NodoInicial;
  int MaxNodePos=2*(tr->mxtips -1);
  int MaxD;
  
  for(int k=tr->mxtips + 1;k<=MaxNodePos;k++){
       if ( PseudoRandom::randDouble() <= (pll_percnodes * iter)) {
         
         rearrangeList = pllCreateRearrangeList (1);
           
         NodoInicial = tr->nodep[k];
          
         pllRearrangeSearch (tr, partitions, PLL_REARRANGE_SPR, NodoInicial,  pll_mintranv,  pll_maxtranv* iter,  rearrangeList);  
         pllEvaluateLikelihood (tr, partitions, tr->start, PLL_FALSE, PLL_FALSE); //Se aplican cambios en el scor
         
         if(rearrangeList->entries > 0 ) {
                pllRearrangeCommit (tr, partitions, &(rearrangeList->rearr[0]), PLL_TRUE);
                pllEvaluateLikelihood (tr, partitions, tr->start, PLL_FALSE, PLL_FALSE); //Se aplican cambios en el score
         
                 lk2 = tr->likelihood;
                 R.second = true;
                 if((lk1-lk2 ) > 0){
                          printf ("ERROR: POSITIVES DIFFERENCES Log-likelihood %f - %f \n", lk1 , tr->likelihood); 
                         exit(-1);
                }
         }
         
         pllDestroyRearrangeList (&rearrangeList);
         //cout << "Nodo Seleccionado " << k << " de " << MaxNodePos << " tr.likelid " << tr->likelihood << endl;
       }
  }
       
  //cout << " Tiempo " << ((double) (clock() - t_ini))/ CLOCKS_PER_SEC  << endl;
  
 if(printtrace){
        ComportamientoTime <<  ((double) (clock() - t_ini))/ CLOCKS_PER_SEC << endl; 
        ComportamientoML << lk1-lk2  << endl; 
 }
  
 R.first = lk2; 
 return R;
}


pair<double,boolean>  Phylogeny::PPNSearch(Solution * solution, int NumMaxMovs, int _PPN_NumIteraciones){
  
  pair<double,boolean> R;
  
  int NumApliedMovSPR=0; 
  vector< pair<int , int > > MovSPRPar; //Almacena todos los MovsSPR ya hechos entres ID1 y ID2
  pair<int , int> IDs;        //Mejor Mov SPR ID1 ID2 
    
  clock_t t_ini = clock();
  
    PhyloTree * Pt = (PhyloTree*) solution->getDecisionVariables()[0];
    Phylogeny * problem_ = (Phylogeny *) solution->getProblem();
    TreeTemplate<Node> * tree = Pt->getTree();
    int NextNodeID = tree->getNextId();
    
    //DRTreeParsimony usado para denotar los posibles buenos moviemitos basads en Parsimonia
    DRTreeParsimonyScore * Par =	new DRTreeParsimonyScore(*tree, *problem_->getSites() ,false,true);
    TreeTemplate<Node> * treeP = Par->getTreeP_();
    vector<Node*> nodes = treeP->getNodes();
    int NextNodeIDPar = treeP->getNextId();
    
    int ParsimonyScore;
    int MenorParsimonia, Par1;
    MenorParsimonia = Par1 = R.first =  Par->getScore();
    R.second = false;
    
    Node * Root;
    Node* Nodo1;
    Node* Nodo2;

  int MaxNivel, dmax,d,up,down,u, Iteracion;
  bool b;
  
  pair< pair<double *, int *>, Node *> pPar;
  pair<double *, int *> pnPar;
  
  Root=treeP->getRootNode();
  MaxNivel=0;
  for(int k=0;k<Root->getNumberOfSons();k++){      MaxNivel+= getNivel(Root->getSon(k));   }

  dmax =MaxNivel; 
  for(Iteracion=0; Iteracion<_PPN_NumIteraciones;Iteracion++){
      
        if (Iteracion%(_PPN_NumIteraciones/MaxNivel+1)==0) dmax--;
      
        if(dmax < 3) dmax=3;
        else if(dmax > MaxNivel) dmax=MaxNivel;

        do{
              b=true;
              do {
                  Nodo1 =  nodes[rand()%nodes.size()];
                  if(Nodo1->hasFather()){
                      if(Nodo1->getFather()->hasFather()) b=false;
                  }
              }while(b);

              d=rand()%dmax+1; //Distancia entre 1 y DMax
              u=up=rand()%d+2; //Subir desde 2 a D

              Nodo2=Nodo1;
              Node* prev;
              while(Nodo2->hasFather() && (u>0)) {
                 prev=Nodo2;
                 Nodo2=Nodo2->getFather();
                 u--;
              }

              down=d+2-up-u;

            while(Nodo2->getNumberOfSons()>1 && (down>0)) {
                  Nodo2 = Nodo2->getSon(rand()% Nodo2->getNumberOfSons());
                  if(Nodo1==Nodo2) {
                      Nodo2=Nodo2->getFather();
                  }else  down --;
            }
         }while(!SPRvalide (Nodo1,Nodo2));
         
         
         
         if(!isMov(MovSPRPar,Nodo1->getId() , Nodo2->getId())){
      
             IDs.first=Nodo1->getId();         IDs.second=Nodo2->getId();
             MovSPRPar.push_back(IDs);
             
            pPar=SPR(Nodo1, Nodo2, NextNodeIDPar); //Sobrr TreeParsimonia
            pnPar = pPar.first;

            Par->topologyChangeTested(*new TopologyChangeEvent());

            ParsimonyScore = Par->getScore();
            if( ParsimonyScore < MenorParsimonia){

                SPR(tree->getNode(Nodo1->getId()),tree->getNode(Nodo2->getId()), NextNodeID);

                MenorParsimonia = ParsimonyScore;
                R.second = true;
                
                if(++NumApliedMovSPR>NumMaxMovs) Iteracion=_PPN_NumIteraciones+1; //exit For Iteraciones

            }else {            
                 double *b = pnPar.first;
                if(b[0]==1)
                   SPR(Nodo1,pPar.second,b);
                else
                   SPRreverse(Nodo1,pPar.second,b,NextNodeIDPar);
            }
            delete[] pnPar.first;       delete[] pnPar.second;
       }
  }
  
  delete Par;
  
 
 if(printtrace){
        ComportamientoTimePar <<  ((double) (clock() - t_ini))/ CLOCKS_PER_SEC << endl; 
        ComportamientoMP << Par1 - MenorParsimonia   << endl; 
 }
  
 R.first = (double) MenorParsimonia;
 return R; 
}


bool Phylogeny::isMov(vector< pair<int , int > > vIDS, int IDN1, int IDN2){
    pair<int , int > eIDS;
     
    for (int i=0;i<vIDS.size();i++){
         eIDS = vIDS[i];
         if(eIDS.first == IDN1 and eIDS.second == IDN2)            return true;
    }
    return false;
}

pair< pair<double *, int *>, Node *> Phylogeny::SPR(Node * Nodo1, Node * Nodo2, int &NextIDNode){
    
    pair<pair<double *, int *>,Node*> p;
    int PosNodo;
    double distancetofather=0;
    double * b = new double[5]; //0 tipo, 1 brLenHermano, 2 brLenNodo2, 3 brLenNod1, 4 brLenPadre
    int * IDs = new int[3];  //0 IDNodoPadre, 1 IDNodo1, 2 IDNodo2
    Node * Padre;
    Node * Padre2;
    Node * GP;
    Node * Hermano;
  
    Padre=Nodo1->getFather();
    if (Nodo1->hasDistanceToFather()) b[3] = Nodo1->getDistanceToFather(); else b[3] =0 ;

    if(Padre->getNumberOfSons()==2){ //Si tiene 2 hijos Collapse Brother por Father
       PosNodo= Padre->getSonPosition(Nodo1);
       Hermano = Padre->getSon(PosNodo==0?1:0);

       if (Hermano->hasDistanceToFather()) {
           distancetofather = Hermano->getDistanceToFather();
           b[1] = distancetofather;
       }else      
           b[1] = 0;
       

       //Quito al Padre sin el hermano, y ubico al hermano en vez del Padre
       Padre->removeSon(Hermano);
       GP = Padre->getFather();
       GP->setSon(GP->getSonPosition(Padre),Hermano);
       
       if(Padre->hasDistanceToFather()) {
               distancetofather+=Padre->getDistanceToFather();
               b[4] = Padre->getDistanceToFather();
       }else b[4] =0;
           
       Hermano->setDistanceToFather(distancetofather);
       
       b[0] = 1; //Normal Nodos con dos Hijos
       p.second=Hermano;
      

     }else{ //Si tiene mas de un hermano, no se hace Collapse

       PosNodo= Padre->getSonPosition(Nodo1);
       Hermano = Padre->getSon(PosNodo==0?1:0);

       Padre->removeSon(Nodo1); //NO Elimina el NODO solo lo eliminar dle Vector de Sons

       Padre = new Node(NextIDNode++);
       Padre->addSon(Nodo1);
       
       b[0] = 2; //Nodo con mas de 2 Hijos SPRReverse
       p.second=Hermano;
     }

     distancetofather=0;
     
     Padre2 = Nodo2->getFather();
     
     //cout << "Hermano 1" << endl;
     
     if(Nodo2->hasDistanceToFather()) {
         distancetofather = Nodo2->getDistanceToFather();
         b[2]=distancetofather;
     }else   b[2]=0;
     
     Padre2->setSon(Padre2->getSonPosition(Nodo2),Padre);
     
     Padre->setDistanceToFather(distancetofather/2);
     
     //Agrego al Nodo2 como hijo del Padre
     Padre->addSon(Nodo2);
     
     Nodo2->setDistanceToFather(distancetofather/2);
     //cout << "PAdre 3 " << endl;
     
     IDs[0] = Padre->getId();      IDs[1] = Nodo1->getId();     IDs[2] = Nodo2->getId();
     
     /*if(Padre->hasDistanceToFather() and Nodo1->hasDistanceToFather() and Nodo2->hasDistanceToFather())
        if(Padre->getDistanceToFather()>0 and Nodo1->getDistanceToFather()>0 and Nodo2->getDistanceToFather()>0)
        cout << "Ramas a Optimizar dentro SPR " << Padre->getDistanceToFather() << " - " << Nodo1->getDistanceToFather() << " - " << Nodo2->getDistanceToFather() << endl;
     */
     pair<double *, int *> pn;
     pn.first=b; pn.second=IDs;
     
     p.first = pn;
     return p;
}

//SPR de Reserva a estado Anterior
void Phylogeny::SPR(Node * Nodo1, Node * Nodo2, double *b){
    
   //b -- 0 tipo, 1 brLenHermano, 2 brLenNodo2, 3 brLenNod1, 4 brLenPadre
    int PosNodo;
    Node * Padre;
    Node * Padre2;
    Node * GP;
    Node * Hermano;
  
    Padre=Nodo1->getFather();
    if (b[3]!= 0) Nodo1->setDistanceToFather(b[3]);

    PosNodo= Padre->getSonPosition(Nodo1);
   Hermano = Padre->getSon(PosNodo==0?1:0);
   if (b[2]!= 0) Hermano->setDistanceToFather(b[2]);

   //cout << "Hermano " << endl;

   //Quito al Padre sin el hermano, y ubico al hermano en vez del Padre
   Padre->removeSon(Hermano);
   GP = Padre->getFather();
   GP->setSon(GP->getSonPosition(Padre),Hermano);

   //cout << "Set Son Hermano" << endl;

  Padre2 = Nodo2->getFather();
     
  //cout << "PAdre " << endl;

   Padre2->setSon(Padre2->getSonPosition(Nodo2),Padre);
   if(b[4] != 0) Padre->setDistanceToFather(b[4]);

   //cout << "PAdre2 " << endl;

   //Agrego al Nodo2 como hijo del Padre
   Padre->addSon(Nodo2);
   if(b[1] != 0) Nodo2->setDistanceToFather(b[1]);
     //cout << "PAdre 3 " << endl;
     
     //cout << "Herman ID " << Hermano->getId() << endl;
      
}

void Phylogeny::SPRreverse(Node * Nodo1, Node * Nodo2 , double *b,int &NextIDNode){
    
    int PosNodo;
    Node * Padre;
    Node * Padre2;
    Node * GP;
    Node * Hermano;
  
    Padre=Nodo1->getFather();
    

    PosNodo= Padre->getSonPosition(Nodo1);
    Hermano = Padre->getSon(PosNodo==0?1:0);
    if (b[2]!= 0) Hermano->setDistanceToFather(b[2]);
   
     //Quito al Padre sin el hermano, y ubico al hermano en vez del Padre
     Padre->removeSon(Hermano);
     GP = Padre->getFather();
     GP->setSon(GP->getSonPosition(Padre),Hermano);


     delete Padre; NextIDNode--;
     Padre2 = Nodo2->getFather();
     //Agrego al Nodo2 como hijo del Padre
     Padre2->addSon(Nodo1);
     if (b[3]!= 0) Nodo1->setDistanceToFather(b[3]);
     //Nodo1->setDistanceToFather(distancetofather/2);
}


   
int Phylogeny::getNivel(Node* nodo){
    int Max=-1;
      
    if(nodo->isLeaf()){
        return 0;
    }else{
        
        vector<int> V;
        for(int i=0;i<nodo->getNumberOfSons();i++){
            V.push_back(getNivel(nodo->getSon(i)));
        }
        for(int i=0;i<nodo->getNumberOfSons();i++){
            if(V[i]>Max) Max=V[i];
        }
        return Max+1;
    }
}

int  Phylogeny::SPRvalide (Node* N1, Node* N2) {
    if (!N2->hasFather()) return 0;
    if (N1->getFather()==N2->getFather()) return 0;
    if (N1->getFather()==N2) return 0;
    
    return 1;
}

void Phylogeny::BranchLengthOptimization(Solution *  solution, string Method , int NumEvals, double Tolerance){
    if (Method==OptimizationTools::OPTIMIZATION_NEWTON)
         PLLNewtonBLOptimization(NumEvals);
    else
        BppGradientBLOptimization(solution,NumEvals,Tolerance);
        
}

double Phylogeny::BppGradientBLOptimization(Solution *  solution, int NumEvals , double Tolerancia){

//   PhyloTree * Pt = (PhyloTree*) solution->getDecisionVariables()[0];
//   Phylogeny * p = (Phylogeny *) solution->getProblem();
//   TreeTemplate<Node> * tree = Pt->getTree();
//  
//   // clock_t t_ini = clock();     
//
//   SubstitutionModel * Model = new GTR(p->alphabet,
//           Pt->CT/Pt->AG,Pt->AT/Pt->AG,Pt->GT/Pt->AG,Pt->AC/Pt->AG,Pt->CG/Pt->AG,
//           Pt->piA, Pt->piC, Pt->piG, Pt->piT); 
//  DiscreteDistribution * RateDist =  new GammaDiscreteDistribution(p->NumCat,Pt->alpha,Pt->beta);
//  RHomogeneousTreeLikelihood * tl1= new RHomogeneousTreeLikelihood(*tree, *p->sites2_, Model, RateDist, true,false,true);
//  tl1->initialize();
//   
//  //double lk1= tl1->getValue();
//  //cout << "Likelihood I " << lk1 << endl;
//    
//  tl1->enableFirstOrderDerivatives(true);
//  tl1->enableSecondOrderDerivatives(true);
//   
//  ParameterList pl = tl1->getBranchLengthsParameters();
//       
//  Optimizer* optimizer;
//   optimizer =new ConjugateGradientMultiDimensions(tl1);
//   //optimizer =new BfgsMultiDimensions(tl1);
//   //optimizer = new PseudoNewtonOptimizer(tl1);
//   optimizer->setMaximumNumberOfEvaluations(NumEvals);
//   optimizer->getStopCondition()->setTolerance(Tolerancia);
//   optimizer->setProfiler(NULL);
//   optimizer->setMessageHandler(NULL);
//   optimizer->setVerbose(0);
//   optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
//
//   optimizer->init(pl);
//   double lk= optimizer->optimize();
//   
//   TreeTemplate<Node>*  tLik = tl1->getTreeP_();
//   vector<Node *> nodos= tree->getNodes();
//   vector<Node *> nodos2= tLik->getNodes();
//   Node * nodo; 
//   Node * nodo2;
//  
//   if(nodos.size()!=nodos2.size()) { 
//       cout << "Warning: Number of Nodes are different: " << nodos.size() << " TLik " << nodos2.size() << endl;
//       
//       delete tree;
//       tree = new TreeTemplate<Node>();
//       tree ->setRootNode(TreeTemplateTools::cloneSubtree<Node>(*tLik->getRootNode()));
//            
//       //Newick * newickprinttree = new Newick;
//       //newickprinttree->write(*tree,"tree");
//       //newickprinttree->write(*tLik,"tlik");
//      //delete newickprinttree;
//       
//   }else{
//       
//        for(int i=0;i<nodos.size();i++){
//           nodo=nodos[i]; nodo2=nodos2[i];
//           if(nodo2->hasDistanceToFather())
//                nodo->setDistanceToFather(nodo2->getDistanceToFather());
//           else
//             nodo->deleteDistanceToFather();
//        }
//   }
//   
//   //ComportamientoTime <<  ((double) (clock() - t_ini))/ CLOCKS_PER_SEC << endl; 
//   //ComportamientoML << lk1 - lk   << endl; 
//   //cout <<"Tiempo " << ((double) (clock() - t_ini))/ CLOCKS_PER_SEC  << " DifereML " << lk1-lk << endl;
//   
//   delete optimizer;
//   delete tl1;
//   delete Model;
//   delete RateDist;
//
//   return lk;
}


void Phylogeny::PLLNewtonBLOptimization(int numevals){
    
    pllOptimizeBranchLengths (tr, partitions, numevals);
    
}

double  Phylogeny::PLLOptimizarParamModeloSust(double epsilon){
    
 pllOptimizeModelParameters(tr, partitions, epsilon);
  
 alpha = partitions->partitionData[0]->alpha;
  
  if(model_frequences=="user"){
    piA = partitions->partitionData[0]->frequencies[0];
    piC = partitions->partitionData[0]->frequencies[1];
    piG = partitions->partitionData[0]->frequencies[2];
    piT = partitions->partitionData[0]->frequencies[3];

    AC= partitions->partitionData[0]->substRates[0];
    AG= partitions->partitionData[0]->substRates[1];
    AT= partitions->partitionData[0]->substRates[2];
    CG= partitions->partitionData[0]->substRates[3];
    CT= partitions->partitionData[0]->substRates[4];
    GT= partitions->partitionData[0]->substRates[5];
  }
  
  return tr->likelihood;
}


double Phylogeny::OptimizarParamModeloSust(Solution * solution){

//    cout << "Optimizing Parameters of the Evolutionary Model" <<endl ;
//    Phylogeny * Prob = (Phylogeny *) solution->getProblem();
//    PhyloTree * Pt  = (PhyloTree*) solution->getDecisionVariables()[0];
//    TreeTemplate<Node> * t;  
//    t = Pt->getTree();
//    
//   SubstitutionModel * Model = new GTR(Prob->alphabet,
//           Pt->CT/Pt->AG,Pt->AT/Pt->AG,Pt->GT/Pt->AG,Pt->AC/Pt->AG,Pt->CG/Pt->AG,
//           Pt->piA, Pt->piC, Pt->piG, Pt->piT); 
//
//   DiscreteDistribution * RateDist =  new GammaDiscreteDistribution(Prob->NumCat,Pt->alpha,Pt->beta);
//      
//   DiscreteRatesAcrossSitesTreeLikelihood * tl1= new RHomogeneousTreeLikelihood(*t, *Prob->sites2_,Model, RateDist, true,false,true);
//   tl1->initialize();
//
//   double lk= tl1->getValue();
//
//    DerivableSecondOrder* f = tl1;
//    MetaOptimizerInfos* desc = new MetaOptimizerInfos();
//    MetaOptimizer* poptimizer = 0;
//    
//     ParameterList plsm,plrd;
//     
//    if (MetodoOptimizacionSubstModel == OptimizationTools::OPTIMIZATION_BRENT){
//               
//           plsm = tl1->getSubstitutionModelParameters();
//           desc->addOptimizer("Substitution model parameter", new SimpleMultiDimensions(f), plsm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
//            //for(int m=0;m<plsm.size();m++){ cout << "Par " << m << " " << plsm[m].getName() << " " << plsm[m].getValue() << endl;}
//          
//           if(OptimizacionRateDist){
//                plrd = tl1->getRateDistributionParameters();
//                desc->addOptimizer("Rate distribution parameter", new SimpleMultiDimensions(f), plrd.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
//           }
//          //for(int m=0;m<plrd.size();m++){cout << "Par " << m << " " << plrd[m].getName() << " " << plrd[m].getValue() << endl;}
//
//           poptimizer = new MetaOptimizer(f, desc, 0);
//
//     }   else { //OPTIMIZATION_BFGS 
//     
//               vector<string> vNameDer;
//               AbstractNumericalDerivative* fnum = new ThreePointsNumericalDerivative(f);
//               
//               plsm = tl1->getSubstitutionModelParameters();
//               vNameDer = plsm.getParameterNames();
//          
//               if(OptimizacionRateDist){
//                    plrd = tl1->getRateDistributionParameters();
//                    vector<string> vNameDer2 = plrd.getParameterNames();
//                    vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());
//               }
//               
//               fnum->setParametersToDerivate(vNameDer);
// 
//              desc->addOptimizer("Rate & model distribution parameters", new BfgsMultiDimensions(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);
//              poptimizer = new MetaOptimizer(fnum, desc, 100);
//        }
//    
//        poptimizer->setVerbose(0);
//        poptimizer->setMaximumNumberOfEvaluations(NumIterOptSubstModel);
//        poptimizer->getStopCondition()->setTolerance(TolerenciaOptSubstModel);
//        
//        poptimizer->setProfiler(NULL);
//        poptimizer->setMessageHandler(NULL);
//  
//        poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
//        
//        ParameterList pl2 = plsm;
//        pl2.addParameters(plrd);
//        
//        poptimizer->init(pl2); //tl1->getParameters()
//        
//        double lk2=poptimizer->optimize();
//        //cout << "Optimización SubstModel Lik " << lk << " OptimiLik " <<  lk2 << " Mejoramiento " << lk-lk2<< endl;
//        
//        ParameterList pl = poptimizer->getParameters();
//       
//        for(int m=0;m<pl.size();m++){ cout << pl[m].getName() << "=" << pl[m].getValue() << " - ";      }
//        cout << "NumberOfEvaluations " << poptimizer->getNumberOfEvaluations() << endl;
//   
//       double theta_  = pl.hasParameter("GTR.theta")?pl.getParameterValue("GTR.theta"):0;
//       double theta1_ = pl.hasParameter("GTR.theta1")?pl.getParameterValue("GTR.theta1"):0;
//       double theta2_ = pl.hasParameter("GTR.theta2")?pl.getParameterValue("GTR.theta2"):0;
//
//        if(theta_!=0 and theta1_!=0 and theta2_!=0){
//            Pt->piA = theta1_ * (1. - theta_);
//            Pt->piC = (1. - theta2_) * theta_;
//            Pt->piG = theta2_ * theta_;
//            Pt->piT =(1. - theta1_) * (1. - theta_);
//        }
//    
//        Pt->CT= pl.getParameterValue("GTR.a")*Pt->AG;
//        Pt->AT= pl.getParameterValue("GTR.b")*Pt->AG;       
//        Pt->GT= pl.getParameterValue("GTR.c")*Pt->AG;
//        Pt->AC= pl.getParameterValue("GTR.d")*Pt->AG;
//        Pt->CG= pl.getParameterValue("GTR.e")*Pt->AG;
//       
//        if(pl.hasParameter("Gamma.alpha"))  Pt->alpha = pl.getParameterValue("Gamma.alpha");
//        if(pl.hasParameter("Gamma.beta"))  Pt->beta = pl.getParameterValue("Gamma.beta");
//    
//         
//         delete tl1;
//         delete Model;
//         delete RateDist;
//         delete poptimizer;
//        
//         return lk2;
       
}



SolutionSet * Phylogeny::generateInitialPopulation(int populationSize){

  SolutionSet * population = new SolutionSet(populationSize);

  ApplicationTools::displayTask("Initial Population", true);

  Solution * newSolution;
  
  for (int i = 0; i < populationSize; i++) {
    
    newSolution = new Solution(this);
    
    double Parsimony=getParsimony(newSolution);
     
    setTreeToPLL(newSolution, true);
    if(StartingOptRamas)
        PLLNewtonBLOptimization(StartingNumIterOptRamas);
    
    if(OptimizacionSubstModel)
        PLLOptimizarParamModeloSust(TolerenciaOptSubstModel); 
    
    if(OptimizacionSubstModel or StartingOptRamas){
        setPLLToTree(newSolution);
    }
    
    evaluate(newSolution,Parsimony, (double)tr->likelihood );
    
    cout << "Initial tree " << i << " Scores P: " << newSolution->getObjective(0) << " L: " << newSolution->getObjective(1)*-1 << endl;
    
    population->add(newSolution);
    
  } 
  
  //Delete Temporal InitialTrees
  for (int i = 0; i < trees.size(); i++){
	  delete trees[i];
  }
  trees.clear();
  
  ApplicationTools::displayTaskDone();
  
  return population;
}


void Phylogeny::setTreeToPLL(Solution *solution, boolean hasToEvaluate) {
   
    PhyloTree * Pt = (PhyloTree*) solution->getDecisionVariables()[0];
    string strtree=TreeTemplateTools::treeToParenthesis(*Pt->getTree()); 
 
    pllNewickTree *newick = pllNewickParseString(strtree.c_str());
   pllTreeInitTopologyNewick(tr, newick, PLL_FALSE);
   pllNewickParseDestroy(&newick);

   if(hasToEvaluate){
        pllEvaluateLikelihood (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
   }
}

void Phylogeny::setPLLToTree(Solution *solution) {
    
    PhyloTree * Pt = (PhyloTree*) solution->getDecisionVariables()[0];
    
    pllTreeToNewick(tr->tree_string, tr, partitions, tr->start->back, true, true, false, false, false, 1, false, false);
    string strTree(tr->tree_string);
    
    delete Pt->getTree();
    Pt->setTree(TreeTemplateTools::parenthesisToTree(strTree));

}

void Phylogeny::setPopulationNotNew(SolutionSet * population){
      Solution * sol;
      PhyloTree * Pt;
      for(int i=0; i<population->size(); i++){
             sol=population->get(i);
             Pt = (PhyloTree*) sol->getDecisionVariables()[0];
             Pt->setNewSolution(false);
       }
    
}

int Phylogeny::getNewSolutionsCount(SolutionSet * population){
      int num=0 ; 
      Solution * sol;
      PhyloTree * Pt;
      for(int i=0; i<population->size(); i++){
             sol=population->get(i);
             Pt = (PhyloTree*) sol->getDecisionVariables()[0];
             if(Pt->isNewSolution()) num++;
       }
      return num;
    
}

void Phylogeny::updateInterval(SolutionSet * population){
      cout << "Updating and Optimizing Parameters.." << endl;
      Solution * sol;
      for(int i=0; i<population->size(); i++){
             sol=population->get(i);
             setTreeToPLL(sol, true);
             PLLNewtonBLOptimization(OptimizationNumIterOptRamas);
             setPLLToTree(sol);
             evaluate(sol,(double)sol->getObjective(0), (double)tr->likelihood);
                 
       }
      cout << "Done.." << endl;
    
}

void Phylogeny::printPreliminarResults(int NumExp, SolutionSet* population){
    ostringstream num;
    num<< NumExp;
        
   cout << "Printing front Approximation (only non dominated solutions). " << num.str() ;
   Ranking *ranking = new Ranking(population);
   SolutionSet * FrenteOP = new SolutionSet(ranking->getSubfront(0)->size());
   for (int i=0;i<ranking->getSubfront(0)->size();i++) {
            FrenteOP->add(new Solution(ranking->getSubfront(0)->get(i)));
   }
   delete ranking;
         
   cout << " Pareto front Size: " <<  FrenteOP->size() << endl;
   FrenteOP->sort( new DominanceComparator());
   FrenteOP->printVariablesToFile("VAR" + num.str());
   FrenteOP->printObjectivesToFile2("FUN" + num.str());
    
   delete FrenteOP;
}
