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
 * @file PhylogeneticMutation.cpp
 */

#include <PhylogeneticMutation.h>


/**
 * Constructor
 * Creates a new instance of the PhylogeneticMutation  operator
 */
PhylogeneticMutation::PhylogeneticMutation(map<string, void *> parameters)
: Mutation(parameters) {

   mutationProbability_=0;  Metodo="NNI";
    
   if (parameters["metodo"] != NULL)              Metodo = *(string *) parameters["metodo"];
   if (parameters["probability"] != NULL)         mutationProbability_ = *(double *) parameters["probability"];
   
   //random_gsl = gsl_rng_alloc(gsl_rng_default);
   //randomgsl_alpha = 500;

   
} // PhylogeneticMutation



void PhylogeneticMutation::printParameters(){
   
   cout << "************** Parameters of the Mutation Operator  ************ " << endl;
   cout << "Mutation Method: " <<  Metodo << endl;
   cout << "mutationProbability: " << mutationProbability_ << endl;
  
   cout << endl;
}

/**
 * Destructor
 */
PhylogeneticMutation::~PhylogeneticMutation() { 

} // ~PhylogeneticMutation


/**
 * Perform the mutation operation
 * @param probability Mutation probability
 * @param solution The solution to mutate
 */
void * PhylogeneticMutation::doMutation(double mutationProbability_, Solution *solution) {
    
    if ( PseudoRandom::randDouble() <= mutationProbability_) {
        
                if(Metodo=="SPR") 
                    SPR(solution); 
                else if(Metodo=="TBR") 
                    TBR(solution);
                else
                     NNI(solution); 
                
    }

} // doMutation

void PhylogeneticMutation::TBR(Solution * solution){
    vector<int> nodosIDs;
    Variable **variables = solution->getDecisionVariables();
    PhyloTree * Pt = (PhyloTree*) variables[0];
    TreeTemplate<Node> * tree = Pt->getTree();

    Node * nodo;
    Node * nodoi;
    Node * nodoj;
    Node * nodoPadre;
    Node * nodoSubTree;

    nodosIDs = tree->getNodesId();
    do{
              do{
                      nodo = tree->getNode(RandomTools::pickOne(nodosIDs,true));
              }while(nodo->isLeaf() || nodo->getNumberOfSons() < 2 || !nodo->hasFather());

              if (RandomTools::flipCoin()) {
                      nodoi= nodo->getSon(0); nodoj= nodo->getSon(1);
              }else{
                      nodoi= nodo->getSon(1); nodoj= nodo->getSon(0);
              }

     }while(nodoi->isLeaf());

     TreeTemplate<Node> * subtree = new TreeTemplate<Node>(nodoi);
     subtree->resetNodesId();
     nodosIDs = subtree->getNodesId();
     do{
            nodoSubTree = subtree->getNode(RandomTools::pickOne(nodosIDs,true));
     }while(nodoSubTree->isLeaf());

     subtree->rootAt(nodoSubTree->getId());

     nodoPadre = nodo->getFather();
     nodoPadre->setSon(nodoPadre->getSonPosition(nodo),nodoj);
     tree->resetNodesId();
     nodosIDs = tree->getNodesId();
     do{
              nodo = tree->getNode(RandomTools::pickOne(nodosIDs,true));
     }while(nodo->isLeaf());


     int posSon;
     if (RandomTools::flipCoin()) posSon=0; else posSon=1;

     Node * nuevonodo = new Node();
     nuevonodo->addSon(nodo->getSon(posSon));
     nuevonodo->addSon(subtree->getRootNode());

     nodo->setSon(posSon,nuevonodo);
     tree->resetNodesId();

}


bool PhylogeneticMutation::NNIValidate(Node * Nodo){

    if (Nodo->getNumberOfSons()>1){
          if (!Nodo->getSon(0)->isLeaf() and !Nodo->getSon(1)->isLeaf()) return true;    
    }
    return false;

}

void PhylogeneticMutation::NNI(Solution * solution){
    
    PhyloTree * Pt = (PhyloTree*) solution->getDecisionVariables()[0];
    TreeTemplate<Node> * tree = Pt->getTree();

    Node * NodoSel;
    vector<Node *> nodes = tree->getNodes();

    do{
          NodoSel =  nodes[PseudoRandom::randInt(0, nodes.size() - 1)];
          
    }while(!NNIValidate(NodoSel));
    
    Node * Nodo1;
    Node * Nodo2;
    int Pos1, Pos2;

    Pos1 = PseudoRandom::randInt(0, NodoSel->getSon(0)->getNumberOfSons()-1);
    Pos2 = PseudoRandom::randInt(0, NodoSel->getSon(1)->getNumberOfSons()-1);

    Nodo1=  NodoSel->getSon(0)->getSon(Pos1);
    Nodo2=  NodoSel->getSon(1)->getSon(Pos2);

    NodoSel->getSon(0)->setSon(Pos1,Nodo2);
    NodoSel->getSon(1)->setSon(Pos2,Nodo1);

}



//void PhylogeneticMutation::ModificarRamasDistGamma(Solution * solution){
//    
//    /*Also modifies branch lengths in 	order to improve the tree likelihood value.
//    * Some branches, chosen at random, have their
//    * lengths multiplied by a factor obtained from a Γ-distribution (Lewis, 1998)
//    */
//    
//    TreeTemplate<Node> * tree = ((PhyloTree*) solution->getDecisionVariables()[0])->getTree();
//    
//    //if( PseudoRandom::randDouble() > 0.34 ){
//	 double gamma; double distance;
//	 vector<Node *> nodes = tree->getNodes();
//	 for(unsigned int i = 0; i < nodes.size(); i++){
//             if (PseudoRandom::randDouble() <= 0.5) {
//                if (nodes[i]->hasDistanceToFather()){
//			 distance = nodes[i]->getDistanceToFather();
//        		 gamma =  gsl_ran_gamma (random_gsl, randomgsl_alpha, (1.0/randomgsl_alpha) );
//			 if (distance * gamma != 0)	 nodes[i]->setDistanceToFather(distance * gamma);
//		 }
//             }
//	}
//   // }
//    
//}
		
void PhylogeneticMutation::SPR(Solution * solution){
    
    PhyloTree * Pt = (PhyloTree*) solution->getDecisionVariables()[0];
    TreeTemplate<Node> * tree = Pt->getTree();
    int NextIDNode = tree->getNextId();
   
    bool b;
    Node* Nodo1;
    Node* Nodo2;
    
     vector<Node*> nodes = tree->getNodes();
     do{
        b=true;
        do {
               Nodo1 =  nodes[rand()%nodes.size()];
               if(Nodo1->hasFather()){
                   if(Nodo1->getFather()->hasFather()) b=false;
                }
        }while(b);
     
        Nodo2 =  nodes[rand()%nodes.size()];
        
    }while(!SPRvalide (Nodo1,Nodo2));     
    
    
    int PosNodo;
    double distancetofather=0;
    Node * Padre;
    Node * Padre2;
    Node * GP;
    Node * Hermano;
  
    Padre=Nodo1->getFather();
    if(Padre->getNumberOfSons()==2){ //Si tiene 2 hijos Collapse Brother por Father
       PosNodo= Padre->getSonPosition(Nodo1);
       Hermano = Padre->getSon(PosNodo==0?1:0);

       if (Hermano->hasDistanceToFather()) {
           distancetofather = Hermano->getDistanceToFather();
       }
       
       //Quito al Padre sin el hermano, y ubico al hermano en vez del Padre
       Padre->removeSon(Hermano);
       GP = Padre->getFather();
       GP->setSon(GP->getSonPosition(Padre),Hermano);
       
       if(Padre->hasDistanceToFather()) {
               distancetofather+=Padre->getDistanceToFather();
       }
           
       Hermano->setDistanceToFather(distancetofather);

     }else{ //Si tiene mas de un hermano, no se hace Collapse

       PosNodo= Padre->getSonPosition(Nodo1);
       Hermano = Padre->getSon(PosNodo==0?1:0);

       Padre->removeSon(Nodo1); //NO Elimina el NODO solo lo eliminar dle Vector de Sons

       Padre = new Node(NextIDNode++);
       Padre->addSon(Nodo1);
       
     }

     distancetofather=0;
     
     Padre2 = Nodo2->getFather();
     
     if(Nodo2->hasDistanceToFather()) {
         distancetofather = Nodo2->getDistanceToFather();
     }

     Padre2->setSon(Padre2->getSonPosition(Nodo2),Padre);
     Padre->setDistanceToFather(distancetofather/2);
     
     //Agrego al Nodo2 como hijo del Padre
     Padre->addSon(Nodo2);
     
     Nodo2->setDistanceToFather(distancetofather/2);

     /*if(Padre->hasDistanceToFather() and Nodo1->hasDistanceToFather() and Nodo2->hasDistanceToFather())
        if(Padre->getDistanceToFather()>0 and Nodo1->getDistanceToFather()>0 and Nodo2->getDistanceToFather()>0)
        cout << "Ramas a Optimizar dentro SPR " << Padre->getDistanceToFather() << " - " << Nodo1->getDistanceToFather() << " - " << Nodo2->getDistanceToFather() << endl;
     */
}


int  PhylogeneticMutation::SPRvalide (Node* N1, Node* N2) {
    if (!N2->hasFather()) return 0;
    if (N1->getFather()==N2->getFather()) return 0;
    if (N1->getFather()==N2) return 0;
    if (N1 == N2) return 0;
    
    return 1;
}


void * PhylogeneticMutation::execute(void *object) {
  Solution *solution = (Solution *)object;
  // TODO: VALID_TYPES?
  //double probability = *(double *)getParameter("probability");
  doMutation(mutationProbability_, solution);
  return solution;
} // execute
 