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
 * @file PhyloTree.h
 */


#ifndef PHYLOTREE_H_
#define PHYLOTREE_H_

#include <Variable.h>
#include <PseudoRandom.h>
#include <sstream>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>


using namespace bpp;

class PhyloTree : public Variable {

public:

  PhyloTree();
  //PhyloTree(string TreeDesc);
  PhyloTree(TreeTemplate<Node> * userTree);
  //PhyloTree(Tree * userTree);
  PhyloTree(vector<string> leavesnames);
  PhyloTree(string newickFile);
  PhyloTree (PhyloTree * phylotree);

  //PhyloTree(Variable * variable);
  ~PhyloTree();

  Variable * deepCopy();
  string toString() ;
  

  double getValue();
  void setValue(double value);
  double getLowerBound();
  double getUpperBound();

  Node * selectrandomnode();
  Node * getRandomLeaf();
  Node * selectrandomnodeToCross();
  void fixbugdropleaf(string leafname);
  TreeTemplate<Node> * cloneSubtree(int NodeID);
  int getNumberOfLeaves();
  void setRootNode(Node *nuevoroot);
  void resetNodesId();
  void writeTree(string TreeFilename);

  void setTree(TreeTemplate<Node> *tree);

  TreeTemplate<Node>* getTree();
  TreeTemplate<Node>* getTreeCopy();
  double getParsimonyScore(){return ParsimonyScore;}
  double getLnLikelihoodValue(){return LnLikelihoodValue;}

  void setParsimonyScore(double ParsimonyScore_){ ParsimonyScore=ParsimonyScore_;}
  void setLnLikelihoodValue(double LnLikelihoodValue_){LnLikelihoodValue = LnLikelihoodValue_;}

  bool isModificada() { return Modificada ; }
  void setModificada(bool Modificada_){ Modificada = Modificada_;}
  
  bool isNewSolution() { return newsolution ; }
  void setNewSolution(bool NewSolution_){ newsolution = NewSolution_;}

private:
  TreeTemplate<Node> * tree_;
  double ParsimonyScore ;
  double LnLikelihoodValue ;
  bool Modificada; //Indica si esta solución ha sido modifica entre Cruce y Mutaciópn
  bool newsolution;

};


#endif /* PHYLOTREE_H_ */
