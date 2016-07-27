/** 
 * MORPHY (version 1.0.0) a software tool for multi-objective 
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
 * @file main.cpp
 */

#include <cstdlib>

#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

#include <Problem.h>
#include <Phylogeny.h>
#include <Solution.h>
#include <PhylogeneticMutation.h>
#include <BinaryTournament.h>
#include <RandomSelection.h>


#include <MORPHY.h>


using namespace std;
using namespace bpp;


void Message()
{
  (*ApplicationTools::message <<"* MORPHY: Multi-Objective softwaRe to Phylogenetic Inference              *").endLine();
  (*ApplicationTools::message  <<"*                                                                        *").endLine();
  (*ApplicationTools::message <<"* Authors:                                            Last Modif. 13/07/16*").endLine();
  (*ApplicationTools::message <<"* Cristian Zambrano-Vega                                                  *").endLine();
  (*ApplicationTools::message <<"* Antonio J. Nebro                                                        *").endLine();
  (*ApplicationTools::message <<"***************************************************************************").endLine();
  cout << endl;
}



int main(int argc, char** argv) {

  clock_t t_ini, t_fin;
  cout.precision(5);
  cout << fixed;

  Problem   * problem   ; // The problem to solve
  Algorithm * algorithm ; // The algorithm to use
  Operator  * mutation  ; // Mutation operator
  Operator  * selection ; // Selection operator
  
  srand (time(NULL));
   
  Message();
  
  BppApplication * objApp = new BppApplication(argc, argv, "Params");

  string NumExp =  ApplicationTools::getStringParameter("experimentid", objApp->getParams(), "1", "", false, false);

  //Problem
  problem =  new Phylogeny(objApp);
  
 //MOEA  
  int populationSize,maxEvaluations,IntervalOptSubsModel;
  
   populationSize= ApplicationTools::getIntParameter("populationsize", objApp->getParams(), 100, "", false, false);
   maxEvaluations = ApplicationTools::getIntParameter("maxevaluations", objApp->getParams(), 2500, "", false, false);
   IntervalOptSubsModel= ApplicationTools::getIntParameter("intervalupdateparameters", objApp->getParams(), 500, "", false, false);

   double umbral= ApplicationTools::getDoubleParameter("umbral", objApp->getParams(), 0.2, "", false, false);
  
    
  algorithm = new MORPHY(problem);
  algorithm->setInputParameter("umbral",&umbral);

  algorithm->setInputParameter("populationSize",&populationSize);
  algorithm->setInputParameter("maxEvaluations",&maxEvaluations);
  algorithm->setInputParameter("intervalupdateparameters",&IntervalOptSubsModel);


  //Operator Mutation
  map<string, void *> parameters;
  double mutationProbability = ApplicationTools::getDoubleParameter("mutation.probability", objApp->getParams(), 0.2, "", false, false);
  string OperadorMutacion= ApplicationTools::getStringParameter("mutation.method", objApp->getParams(), "NNI", "", false, false);
  
  parameters["probability"] = &mutationProbability;
  parameters["metodo"] = &OperadorMutacion;
  mutation = new PhylogeneticMutation(parameters);
  
  //Selection Operator
  string OperadorSeleccion = ApplicationTools::getStringParameter("selection.method", objApp->getParams(), "binarytournament", "", false, false);
  
  parameters.clear();
  if(OperadorSeleccion=="binarytournament"){
      selection = new BinaryTournament(parameters);
  }else if(OperadorSeleccion=="randomselection"){
        selection = new RandomSelection(parameters);     
  }else
      selection = new BinaryTournament(parameters);
      
  algorithm->addOperator("mutation",mutation);
  algorithm->addOperator("selection",selection);
    
  cout <<"********************PARAMETERS**********************************" << endl;
  cout << "NumExp: " << NumExp << endl << endl;
  
  cout <<"**************************MOEA**********************************" << endl ;
  cout << "populationSize: " << populationSize << endl;
  cout << "Umbral: " << umbral << endl;      
  cout << "maxEvaluations: " << maxEvaluations << endl << endl;
  cout << "IntervalUpdateParameters: " << IntervalOptSubsModel << endl;
  
  
  ((Phylogeny*)problem)->printParameters();
  ((PhylogeneticMutation*)mutation)->printParameters();
  
  cout << "********************* Selection Operator ********************* " << endl;
  cout << "Method: " << OperadorSeleccion << endl;
  cout << endl ;
  
  
   cout << "****************** Start of Algorithm ***************" << endl;
    
   time_t timer = time(NULL);    printf("Start:  %s\n", ctime(&timer));
   t_ini = clock();

   SolutionSet * population = algorithm->execute();

   
   Phylogeny* p = (Phylogeny*) problem;
   
   if(p->FinalOptRamas){
        Solution *sol;
        double Lk ;
        cout << "Optimizing Final Solutions - Population Size: " << population->size() << endl;
        
        for(int i=0; i<population->size(); i++){
                sol = population->get(i);
                Lk=sol->getObjective(1);

                cout << "Optimizing Solution " << i + 1 << " Likelihood = -" << Lk << " ";

                p->setTreeToPLL(sol, true);
                p->PLLNewtonBLOptimization(p->FinalNumIterOptRamas);
                Lk=p->tr->likelihood;
                p->setPLLToTree(sol);
                sol->setObjective(1,Lk*-1);

                cout << " final: " << Lk << endl;
                
        }
  
     timer = time(NULL);    printf("Final Optimization %s\n", ctime(&timer));
     
   }
    
    cout << "Optimal Pareto front Approximation (only non dominated solutions)." << endl;
    Ranking *ranking = new Ranking(population);
    SolutionSet * FrenteOP = new SolutionSet(ranking->getSubfront(0)->size());
    for (int i=0;i<ranking->getSubfront(0)->size();i++) {
            FrenteOP->add(new Solution(ranking->getSubfront(0)->get(i)));
    }
    delete ranking;
    cout << "Pareto front Size: " <<  FrenteOP->size() << endl;
    
    cout << "Ordering Solutions" << endl;
    FrenteOP->sort( new DominanceComparator());
    cout << "Variables values have been written to file VAR" << endl;
    FrenteOP->printVariablesToFile("VAR" + NumExp);
    cout << "Objectives values have been written to file FUN" << endl;
    FrenteOP->printObjectivesToFile2("FUN" + NumExp);

    if(p->consensus_tree){
        p->makeConsensus(FrenteOP, "CONSENSUS_TREE" + NumExp);
        cout << endl;
    }
    
    string filename="Parameters"+ NumExp;
    std::ofstream out(filename.c_str());
    out << "PARAMETERS" << endl  
        << "Gamma Shape (alpha): " << p->alpha << endl;
            
    if(p->model_frequences=="user"){
             out << "Frecuencias  piA: " << p->piA << " piC: " << p->piC << " piG: " << p->piG << " piT: " << p->piT << endl 
                 << "GTR relative rates AC: " << p->AC << " AG: " << p->AG << " AT: " << p->AT << " CG: " << p->CG << " CT: " << p->CT << " GT: " << p->GT << endl; 
    }
    out.close();
    
    t_fin = clock(); 
    double secs = ((double) (t_fin - t_ini))/ CLOCKS_PER_SEC;
    cout << "Total execution time: " << secs << "s" << endl;
    
    timer = time(NULL);    printf("End of Algorithm %s\n", ctime(&timer));
        
  
  delete objApp;
  delete mutation;
  delete selection;
  //delete problem;
  delete algorithm;
  
    
}



