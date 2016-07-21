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
 * @file PhylogeneticMutation.h
 */

#ifndef PhylogeneticMutation_H_
#define PhylogeneticMutation_H_

#include <Mutation.h>
#include <Solution.h>
#include <math.h>

#include <PseudoRandom.h>
#include <PhyloTree.h>
#include <Phylogeny.h>

#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>


           
using namespace bpp;

/**
  * @class Mutation
  * @brief This class implements a Phylogenetic Mutation operator.
**/
class PhylogeneticMutation : public Mutation {

public:

   PhylogeneticMutation(map<string, void *> parameters);
  ~PhylogeneticMutation();
  void * execute(void *);

  
  void printParameters();
  
 
  
private:
    
  string Metodo;
  double mutationProbability_;
  double distributionIndex_;

  //double randomgsl_alpha;  // shape paramter of gamma distribution
  //gsl_rng *random_gsl; // GSL (GNU Scientific Library) random generator for gamma numbers)

  
  
  void * doMutation(double mutationProbability_, Solution * solution);
  
  void NNI(Solution * solution);
  bool NNIValidate(Node * Nodo);
  
  void SPR(Solution * solution);
  int SPRvalide (Node* N1, Node* N2);
  
  void TBR(Solution * solution);
  
  //void ModificarRamasDistGamma(Solution * solution);
  
  // TODO: VALID_TYPES;

}; // PhylogeneticMutation

#endif /* PhylogeneticMutation */
