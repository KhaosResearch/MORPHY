//  MORPHY.cpp
//
//  Author:
//       Cristian Zambrano-Vega <czambrano@uteq.edu.ec>
//
//  Copyright (C) 2017 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <MORPHY.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Phylogeny.h>

using namespace bpp;


MORPHY::MORPHY(Problem *problem) : Algorithm(problem) {
} // MORPHY


double MORPHY::RFDistance(Solution * sol1, Solution * sol2, bool normalized, int NumberOfTaxas) {

  PhyloTree *Pt1, *Pt2; 
  Pt1 = (PhyloTree *)sol1->getDecisionVariables()[0];
  Pt2 = (PhyloTree *)sol2->getDecisionVariables()[0];
 
  TreeTemplate<Node> * tree1 = Pt1->getTree();
  TreeTemplate<Node> * tree2 = Pt2->getTree();
  
  double distance = TreeTools::robinsonFouldsDistance(*tree1,*tree2, true);
  
  if(normalized)
      distance = distance / (2 * (NumberOfTaxas-3));
  
  return distance;
  
}

/*
 * Runs the MORPHY algorithm.
 * @return a <code>SolutionSet</code> that is a set of non dominated solutions
 * as a result of the algorithm execution
 */
SolutionSet * MORPHY::execute() {

  int populationSize;
  int maxEvaluations;
  int evaluations;
  int IntervalOptSubsModel;
  int NumberOfTaxas;
  double umbral;

//  QualityIndicator * indicators; // QualityIndicator object
  int requiredEvaluations; // Use in the example of use of the
                           // indicators object (see below)

  SolutionSet * population;
  SolutionSet * offspringPopulation;
  SolutionSet * unionSolution;

  Operator * mutationOperator;
  Operator * selectionOperator;

  Distance * distance = new Distance();

  //Read the parameters
  populationSize = *(int *) getInputParameter("populationSize");
  maxEvaluations = *(int *) getInputParameter("maxEvaluations");
  IntervalOptSubsModel = *(int *) getInputParameter("intervalupdateparameters");
  umbral = *(double *) getInputParameter("umbral");

  //Initialize the variables
  evaluations = 0;
  requiredEvaluations = 0;

  //Read the operators
  mutationOperator = operators_["mutation"];
  selectionOperator = operators_["selection"];
  
  Phylogeny * p = (Phylogeny *) problem_;
  population= p->generateInitialPopulation(populationSize);
  evaluations=populationSize;
      
  NumberOfTaxas = p->alignmentData->sequenceCount;
  
  time_t timer;
  // Generations
  //ApplicationTools::displayTask("Generations", true);
  while (evaluations < maxEvaluations) {
      
    // Create the offSpring solutionSet
    offspringPopulation = new SolutionSet(populationSize);

    timer = time(NULL);
    cout << "Evaluations  " <<  evaluations << " " << ctime(&timer) << endl;
    
    Solution * sol; PhyloTree *Pt;
    for (int i = 0; i < populationSize; i++) {
             sol = new Solution(population->get(i));
             Pt = (PhyloTree*) sol->getDecisionVariables()[0];
             Pt->setNewSolution(true);
            
             mutationOperator->execute(sol);
             p->Optimization(sol);
             offspringPopulation->add(sol);
             evaluations ++;
     }
        
    // Create the solutionSet union of solutionSet and offSpring
    unionSolution = population->join(offspringPopulation);
    delete offspringPopulation;

    // Ranking the union
    Ranking * ranking = new Ranking(unionSolution);
    
    int remain = populationSize;
    int index = 0;
    SolutionSet * front = NULL;
    for (int i=0;i<population->size();i++) {
      delete population->get(i);
    }
    population->clear();

    // Obtain the next front
    front = ranking->getSubfront(index);

    while ((remain > 0) && (remain >= front->size())) {
      //Assign crowding distance to individuals
      distance->crowdingDistanceAssignment(front, problem_->getNumberOfObjectives());

      //Add the individuals of this front
      for (int k = 0; k < front->size(); k++) {
        population->add(new Solution(front->get(k)));
      } // for

      //Decrement remain
      remain = remain - front->size();

      //Obtain the next front
      index++;
      if (remain > 0) {
        front = ranking->getSubfront(index);
      } // if
      
    } // while

    // Remain is less than front(index).size, insert only the best one
    if (remain > 0) {  // front contains individuals to insert
      distance->crowdingDistanceAssignment(front, problem_->getNumberOfObjectives());
      Comparator * c = new CrowdingComparator();
      front->sort(c);
      delete c;
      for (int k = 0; k < remain; k++) {
        population->add(new Solution(front->get(k)));
      } // for

      remain = 0;
    } // if

    delete ranking;
    delete unionSolution;

    int num=p->getNewSolutionsCount(population);
    //cout << "NewSolutionsCount Population " << num << " Evals: " << evaluations << endl;
    int value=(int)(umbral*populationSize);
    if(num< value){ //Reset Population
    
        Ranking *ranking = new Ranking(population);
        Solution *sol;
        for(int i=1; i< ranking->getNumberOfSubfronts();i++){
                for (int j=0;j<ranking->getSubfront(i)->size();j++) {
                    sol=ranking->getSubfront(i)->get(j);
                    p->OptimizationFull(sol,2);
                }
        }
        delete ranking;
        
    }
    
    
    p->setPopulationNotNew(population);
            
    //Update Interval
    if(IntervalOptSubsModel > 0){
         if(evaluations%IntervalOptSubsModel==0){ 
            p->updateInterval(population);
        }
    }
    
    
    if(p->NumPreliminarResults>0){
        // if(evaluations%p->NumPreliminarResults==0){ 
                p->printPreliminarResults((int)(evaluations/p->NumPreliminarResults), population);
        //}
    }
    
    
    
  } // while
  

  delete distance;

  // Return as output parameter the required evaluations
//  setOutputParameter("evaluations", &requiredEvaluations);

  // Return the first non-dominated front
  Ranking * ranking = new Ranking(population);
  SolutionSet * result = new SolutionSet(ranking->getSubfront(0)->size());
  for (int i=0;i<ranking->getSubfront(0)->size();i++) {
    result->add(new Solution(ranking->getSubfront(0)->get(i)));
  }
  delete ranking;
  delete population;

  return result;

} // execute
