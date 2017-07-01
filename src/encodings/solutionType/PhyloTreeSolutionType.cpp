/** 
 * MORPHY (version 1.0.0) a software tool for multi-objective 
 * phylogenetic inference. This software integrates features 
 * of the jMetalCpp, Bio++ and PLL frameworks.
 * 
 * Copyright (C) 2017 Cristian Zambrano-Vega, Antonio J. Nebro.
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
 * @file PhyloTreeSolutionType.cpp
 */

#include <PhyloTreeSolutionType.h>
#include <cstddef>

/**
 * Constructor
 * @param problem
 */
PhyloTreeSolutionType::PhyloTreeSolutionType(Problem *problem)
: SolutionType(problem) { }

/**
 * Destructor
 */
PhyloTreeSolutionType::~PhyloTreeSolutionType(){

}

/**
 * Creates the variables of the solution
 * @param decisionVariables
 */
Variable **PhyloTreeSolutionType::createVariables(){
	 int i;
           
	  Variable **variables = new Variable*[problem_->getNumberOfVariables()]; //malloc(sizeof(Real) * problem->getNumberOfVariables());
	  if (variables ==  NULL) {
	    cout << "Error: Impossible to reserve memory for variable type" << endl;
	    exit(-1);
	  }

	  Phylogeny * ph = (Phylogeny*)problem_;

	  //if (ph->initialTrees =="user"){
               if (ph->numarbol >= ph->bootstrapSize) {
			cout << "Error: Impossible to get tree " << ph->numarbol<< " from BootStrap " << ph->bootstrapSize << endl;
                        exit(-1);
		}
                
		//Se envia UN NUEVO Tree a partir de Tree leido del InputTrees File
		variables[0] = new PhyloTree(new TreeTemplate<Node>(*ph->trees[ph->numarbol++]));
		 

//	  }else{
//
//		  variables[0] = new PhyloTree(ph->getLeavesName(),
//                                               ph->kappa, ph->piA, ph->piC, ph->piG, ph->piT, ph->alpha, ph->beta,
//                                               ph->AC,ph->AG,ph->AT,ph->CG, ph->CT,ph->GT);
//		  //cout << "LLamada a Creates Variables RANDOM TREE " << endl;
//	  }


	  return variables;

}//createVariables
