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
 * @file PhyloTree.cpp
 */


#include <PhyloTree.h>
/**
  * Empty constructor.
  * It will only initialize all the variables.
 **/

PhyloTree::PhyloTree() {
	 tree_ = new TreeTemplate<Node>();

	 ParsimonyScore =0; LnLikelihoodValue=0;
         Modificada = false;
         newsolution=true;

} // PhyloTree

/*PhyloTree::PhyloTree(string TreeDesc) {
	tree_ = TreeTemplateTools::parenthesisToTree(TreeDesc);
} // PhyloTree
*/


PhyloTree::PhyloTree(vector<string> leavesnames){

	tree_ = TreeTemplateTools::getRandomTree(leavesnames, false);
	tree_->setBranchLengths(0.05);

	 ParsimonyScore =0; LnLikelihoodValue=0;
         Modificada = false; newsolution=true;
}

//Copia Referencia del userTree
//Debe enviarse UN Nuevo Tree
PhyloTree::PhyloTree(TreeTemplate<Node> * userTree){
    tree_  = userTree;
    ParsimonyScore =0; LnLikelihoodValue=0;
    Modificada = false; newsolution=true;
}


PhyloTree::PhyloTree(string newickFile) {
	Newick newick;
	tree_ = newick.read(newickFile);
        ParsimonyScore =0; LnLikelihoodValue=0;
        Modificada = false; newsolution=true;
} // PhyloTree

//Deep copya copian directamnete
PhyloTree::PhyloTree(PhyloTree * phylotree) {

	tree_ = new TreeTemplate<Node>();
	tree_ ->setRootNode(TreeTemplateTools::cloneSubtree<Node>(*phylotree->getTree()->getRootNode()));

        ParsimonyScore = phylotree->getParsimonyScore();
	LnLikelihoodValue = phylotree->getLnLikelihoodValue();
        Modificada = phylotree->isModificada();
        newsolution = phylotree->isNewSolution();
        

} // PhyloTree

/**
 * Destructor
 */
PhyloTree::~PhyloTree() {
    delete tree_;
}

/**
 * Returns a exact copy of the <code>PhyloTree</code> variable
 * @return the copy
 */
Variable *PhyloTree::deepCopy() {
    return new PhyloTree(this);
} // deepCopy


/**
 * Returns a string representing the object
 * @return The string
 */
string PhyloTree::toString(){
  stringstream ss;
  ss << TreeTemplateTools::treeToParenthesis(*tree_) << endl  ;
  return ss.str();
} // toString



TreeTemplate<Node> * PhyloTree::getTree() {
  return tree_;
} // getTree

TreeTemplate<Node>* PhyloTree::getTreeCopy() {
	TreeTemplate<Node>* Tree = new TreeTemplate<Node>();
	Tree = new TreeTemplate<Node>();
	Tree->setRootNode(TreeTemplateTools::cloneSubtree<Node>(*tree_->getRootNode()));
	return Tree;
} // getTree

void PhyloTree::setTree(TreeTemplate<Node> *tree) {
    tree_ = tree;  //Se hace una copia Integra de tree ->>> tree__
} // setTree



TreeTemplate<Node> * PhyloTree::cloneSubtree(int NodeID){
	return tree_->cloneSubtree(NodeID);
}

Node * PhyloTree::selectrandomnode()
{
	Node * nodo;
	vector<int> nodosIDs = tree_->getNodesId();

	if(nodosIDs.empty()){ cout << "Tree without Nodes " << endl;  exit(-1);	}

	do{
		nodo = tree_->getNode(RandomTools::pickOne(nodosIDs,true));
	}while(!nodo->hasFather() || !nodo->getFather()->hasFather() || nodo->isLeaf() );

	return nodo;
}

Node * PhyloTree::selectrandomnodeToCross()
{
	Node * nodo;
	vector<int> nodosIDs = tree_->getNodesId();
	do{
		nodo =  tree_->getNode(RandomTools::pickOne(nodosIDs,true));
	}while(!nodo->hasFather() || nodo->isLeaf() );

	return nodo;
}

//Retorna Una taxon Aleatorio, siempre y cuando el padre no sea ROOT del árbol
Node * PhyloTree::getRandomLeaf()
{
	Node * nodo;
	vector<int> nodosIDs = tree_->getLeavesId();
	do{
		do{
			nodo =  tree_->getNode(RandomTools::pickOne(nodosIDs,true));
		}while(!nodo->hasFather());
	}while(!nodo->getFather()->hasFather()); //Que el SubTree NO sea el Root

	return nodo;
}

void PhyloTree::fixbugdropleaf(string leafname)
{
	bool fix=false;
	Node* leaf = tree_->getNode(leafname);
	if (leaf->hasFather()){
		if(!leaf->getFather()->hasFather()){
			fix=true;
		}
	}


	TreeTemplateTools::dropLeaf(*tree_,leafname);

	if(fix) tree_ = new TreeTemplate<Node> (tree_->getRootNode());

}




int PhyloTree::getNumberOfLeaves(){
	return tree_->getNumberOfLeaves();
}

void PhyloTree::setRootNode(Node *nuevoroot){
	tree_->setRootNode(nuevoroot);
}

void PhyloTree::resetNodesId(){
	tree_->resetNodesId();
}


void PhyloTree::writeTree(string TreeFilename){
	Newick * newickprinttree = new Newick;
	newickprinttree->write(*tree_,TreeFilename);
	delete newickprinttree;
}



double PhyloTree::getValue(){ return 0;}
void PhyloTree::setValue(double value){   }
double PhyloTree::getLowerBound(){return 0;}
double PhyloTree::getUpperBound(){return 0; }

