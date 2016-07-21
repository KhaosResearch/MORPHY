#include "../../mem_alloc.h"
#include "../../axml.h"

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>



//#include "axml.h"


/** @file searchAlgo.c
    
    @brief Collection of routines for performing likelihood computation and branch optimization.

    Detailed description to appear soon.
*/




extern double accumulatedTime;

extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern char run_id[128];
extern double masterTime;
extern partitionLengths pLengths[MAX_MODEL];
extern char binaryCheckpointName[1024];
extern char binaryCheckpointInputName[1024];

boolean initrav (pllInstance *tr, partitionList *pr, nodeptr p)
{ 
  nodeptr  q;

  if (!isTip(p->number, tr->mxtips)) 
  {      
    q = p->next;

    do 
    {	   
      if (! initrav(tr, pr, q->back))  return PLL_FALSE;
      q = q->next;	
    } 
    while (q != p);  

    newviewGeneric(tr, pr, p, PLL_FALSE);
  }

  return PLL_TRUE;
} 


/** @brief Optimize the length of a specific branch

    Optimize the length of the branch connecting \a p and \a p->back
    for each partition (\a tr->numBranches) in library instance \a tr.
 
    @param tr
      The library instance

    @param pr
      Partition list
 
    @param p
      Endpoints of branch to be optimized 
*/
void update(pllInstance *tr, partitionList *pr, nodeptr p)
{       
  nodeptr  q; 
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

  #ifdef _DEBUG_UPDATE
    double 
      startLH;
  
    evaluateGeneric(tr, p);
  
    startLH = tr->likelihood;
  #endif

  q = p->back;   

  for(i = 0; i < numBranches; i++)
    z0[i] = q->z[i];    

  if(numBranches > 1)
    makenewzGeneric(tr, pr, p, q, z0, PLL_NEWZPERCYCLE, z, PLL_TRUE);
  else
    makenewzGeneric(tr, pr, p, q, z0, PLL_NEWZPERCYCLE, z, PLL_FALSE);

  for(i = 0; i < numBranches; i++)
  {         
    if(!tr->partitionConverged[i])
    {	  
      if(ABS(z[i] - z0[i]) > PLL_DELTAZ)  
      {	      
        tr->partitionSmoothed[i] = PLL_FALSE;
      }	 

      p->z[i] = q->z[i] = z[i];	 
    }
  }
 
  #ifdef _DEBUG_UPDATE
    evaluateGeneric(tr, p);
  
    if(tr->likelihood <= startLH)
      {
        if(fabs(tr->likelihood - startLH) > 0.01)
  	{
  	  printf("%f %f\n", startLH, tr->likelihood);
  	  assert(0);      
  	}
      }
  #endif
}

/** @brief Branch length optimization of specific branches

    Optimize the length of branches that have \a p as an endpoint 

    @param tr
      The library instance

    @param pr
      Partition list

    @param p
      Endpoint of branches to be optimized
*/
void smooth (pllInstance *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  update(tr, pr, p);    /*  Adjust branch */

  if (! isTip(p->number, tr->mxtips)) 
  {                                  /*  Adjust descendants */
    q = p->next;
    while (q != p) 
    {
      smooth(tr, pr, q->back);
      q = q->next;
    }	

    if(numBranches > 1 && !tr->useRecom)
      newviewGeneric(tr, pr,p, PLL_TRUE);
    else
      newviewGeneric(tr, pr,p, PLL_FALSE);
  }
} 

/**  @brief Check whether the branches in all partitions have been optimized
 
     Check if all branches in all partitions have reached the threshold for
     optimization. If at least one branch can be optimized further return \b PLL_FALSE.

     @param tr
       The library instance 

     @return
       If at least one branch can be further optimized return \b PLL_FALSE,
       otherwise \b PLL_TRUE.
             
*/
static boolean allSmoothed(pllInstance *tr, int numBranches)
{
  int i;
  boolean result = PLL_TRUE;

  for(i = 0; i < numBranches; i++)
  {
    if(tr->partitionSmoothed[i] == PLL_FALSE)
      result = PLL_FALSE;
    else
      tr->partitionConverged[i] = PLL_TRUE;
  }

  return result;
}


/** @brief Wrapper function for branch length optimization of the tree
  
    Perform \a maxtimes rounds of branch length optimization by running smooth()
    on all neighbour nodes of node \a tr->start.

    @param tr
      The library instance

    @param maxtimes
      Number of optimization rounds to perform
*/
/* do maxtimes rounds of branch length optimization */
void smoothTree (pllInstance *tr, partitionList *pr, int maxtimes)
{
	nodeptr  p, q;
	int i, count = 0;
    int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

	p = tr->start;
	for(i = 0; i < numBranches; i++)
		tr->partitionConverged[i] = PLL_FALSE;

	while (--maxtimes >= 0)
	{
		for(i = 0; i < numBranches; i++)
			tr->partitionSmoothed[i] = PLL_TRUE;

		smooth(tr, pr, p->back);
		if (!isTip(p->number, tr->mxtips))
		{
			q = p->next;
			while (q != p)
			{
				smooth(tr, pr, q->back);
				q = q->next;
			}
		}
		count++;

		if (allSmoothed(tr, numBranches)) break;
	}

	for(i = 0; i < numBranches; i++)
		tr->partitionConverged[i] = PLL_FALSE;
} 


/** @brief Optimize the branch length of edges around a specific node
    
    Optimize \a maxtimes the branch length of all (3) edges around a given node 
    \a p of the tree of library instance \a tr.

    @param tr
      The library instance

    @param p
      The node around which to optimize the edges

    @param maxtimes
      Number of optimization rounds to perform
*/
void localSmooth (pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes)
{ 
  nodeptr  q;
  int i;
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
  if (isTip(p->number, tr->mxtips)) return;

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionConverged[i] = PLL_FALSE;	

  while (--maxtimes >= 0) 
  {     
    for(i = 0; i < NUM_BRANCHES; i++)
      tr->partitionSmoothed[i] = PLL_TRUE;

    q = p;
    do 
    {
      update(tr, pr, q);
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr, numBranches))
      break;
  }

  for(i = 0; i < NUM_BRANCHES; i++)
  {
    tr->partitionSmoothed[i] = PLL_FALSE; 
    tr->partitionConverged[i] = PLL_FALSE;
  }
}




/** @brief Reset an \a infoList

    Resets an \a infoList by setting elements \a node and \a likelihood
    of each element of the \a bestInfo list structure to \b NULL and
    \a PLL_UNLIKELY, respectively.

    @param iList
      The given \a infoList.
*/
static void resetInfoList(infoList *iList)
{
  int 
    i;

  iList->valid = 0;

  for(i = 0; i < iList->n; i++)    
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = PLL_UNLIKELY;
  }    
}

/** @brief Initialize an \a infoList

    Initialize an \a infoList by creating a \a bestInfo list structure
    of \a n elements and setting the attributes \a node and \a likelihood
    of each element of the \a bestInfo list structure to \b NULL and
    \a PLL_UNLIKELY, respectively.

    @param iList
      The given \a infoList.

    @param n
      Number of elements to be created in the \a bestInfo list.
*/
static void initInfoList(infoList *iList, int n)
{
  int 
    i;

  iList->n = n;
  iList->valid = 0;
  iList->list = (bestInfo *)rax_malloc(sizeof(bestInfo) * (size_t)n);

  for(i = 0; i < n; i++)
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = PLL_UNLIKELY;
  }
}

/** @brief Deallocate the contents of an \a infoList
    
    Deallocate the contents of a given \a infoList by freeing
    the memory used by its \a bestInfo list structure.

    @param iList
      The \a infoList to be used.
*/
static void freeInfoList(infoList *iList)
{ 
  rax_free(iList->list);   
}


/** @brief Insert a record in an \a infoList

    Insert the pair \a likelihood and \node into the \a bestList \a list
    of \a iList \b only if there already exists a pair in \a list 
    which has the \a likelihood attribute smaller than the given \a 
    likelihoodby. The insertion is done by replacing the smallest
    likelihood pair with the new pair.

    @param node
      The given node

    @param likelihood
      The given likelihood

    @param iList
      The given \a infoList where the record will possibly be appended.
*/
static void insertInfoList(nodeptr node, double likelihood, infoList *iList)
{
  int 
    i,
    min = 0;

  double 
    min_l =  iList->list[0].likelihood;

  for(i = 1; i < iList->n; i++)
  {
    if(iList->list[i].likelihood < min_l)
    {
      min = i;
      min_l = iList->list[i].likelihood;
    }
  }

  if(likelihood > min_l)
  {
    iList->list[min].likelihood = likelihood;
    iList->list[min].node = node;
    iList->valid += 1;
  }

  if(iList->valid > iList->n)
    iList->valid = iList->n;
}


/** @brief  Optimize branch lengths of region

    Optimize the branch lenghts of only a specific region. The branch optimization starts
    at a node \a p and is carried out in all nodes with distance upto \a region from \a p.

    @param tr
      The library instance.
    
    @param p
      Node to start branch optimization from.

    @param region
      The allowed node distance from \p for which to still perform branch optimization.
*/
void smoothRegion (pllInstance *tr, partitionList *pr, nodeptr p, int region)
{ 
  nodeptr  q;

  update(tr, pr, p);   /* Adjust branch */

  if (region > 0)
  {
    if (!isTip(p->number, tr->mxtips)) 
    {                                 
      q = p->next;
      while (q != p) 
      {
        smoothRegion(tr, pr, q->back, --region);
        q = q->next;
      }	

      newviewGeneric(tr, pr,p, PLL_FALSE);
    }
  }
}

/** @brief Wrapper function for optimizing the branch length of a region \a maxtimes times

    Optimize the branch lengths of a specific region \a maxtimes times. The branch optimization
    starts at a given node \a p and is carried out in all nodes with distance upto \a region
    from \a p.

    @param tr
      The library instance.

    @param p
      Node to start branch optimization from.

    @param maxtimes
      Number of times to perform branch optimization.

    @pram region
      The allwed node distance from \p for which to still perform branch optimization.

    @todo
      In the previous version (before the model-sep merge) the loops were controlled by tr->numBranches,
      and now they are controlled by a constant NUM_BRANCHES. What is right?
*/
void regionalSmooth (pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes, int region)
{
  nodeptr  q;
  int i;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if (isTip(p->number, tr->mxtips)) return;            /* Should be an error */

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionConverged[i] = PLL_FALSE;

  while (--maxtimes >= 0) 
  {	
    for(i = 0; i < NUM_BRANCHES; i++)
      tr->partitionSmoothed[i] = PLL_TRUE;

    q = p;
    do 
    {
      smoothRegion(tr, pr, q, region);
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr, numBranches))
      break;
  }

  for(i = 0; i < NUM_BRANCHES; i++) {
    tr->partitionSmoothed[i] = PLL_FALSE;
    tr->partitionConverged[i] = PLL_FALSE;
  }
} 




/* @brief Split the tree into two components and optimize new branch length

   Split the tree into two components. The disconnection point is node \a p.
   First, a branch length is computed for the newly created branch between nodes
   \a p->next->back and \a p->next->next->back and then the two nodes are
   connected (hookup). Disconnection is done by setting \a p->next->next->back
   and \a p->next->back to \b NULL.

   @param tr
     The library instance

   @param p
     The node at which the tree should be decomposed into two components.

   @param numBranches
     Number of branches per partition

   @return q
     Node from the disconnected component

   @todo
     Why do we return this node?
*/
nodeptr  removeNodeBIG (pllInstance *tr, partitionList *pr, nodeptr p, int numBranches)
{  
  double   zqr[numBranches], result[numBranches];
  nodeptr  q, r;
  int i;

  q = p->next->back;
  r = p->next->next->back;

  for(i = 0; i < numBranches; i++)
    zqr[i] = q->z[i] * r->z[i];        

  makenewzGeneric(tr, pr, q, r, zqr, PLL_ITERATIONS, result, PLL_FALSE);

  for(i = 0; i < numBranches; i++)        
    tr->zqr[i] = result[i];

  hookup(q, r, result, numBranches); 

  p->next->next->back = p->next->back = (node *) NULL;

  return  q; 
}

/** @brief Split the tree into two components and recompute likelihood

    Split the tree into two component. The disconnection point is node \a p.
    Set the branch length of the new node between \a p->next->back and
    \a p->next->next->back to \a tr->currentZQR and then decompose the tree
    into two components by setting \a p->next->back and \a p->next->next->back
    to \b NULL.

    @param tr
      The library instance

    @param p
      The node at which the tree should be decomposed into two components.

    @return q
      the node after \a p

    @todo
      Why do we return this node? Why do we set to tr->currentZQR and not compute
      new optimized length? What is tr->currentZQR? 
*/
nodeptr  removeNodeRestoreBIG (pllInstance *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q, r;

  q = p->next->back;
  r = p->next->next->back;  

  newviewGeneric(tr, pr,q, PLL_FALSE);
  newviewGeneric(tr, pr,r, PLL_FALSE);

  hookup(q, r, tr->currentZQR, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  p->next->next->back = p->next->back = (node *) NULL;

  return  q;
}

/* @brief

   @todo
     What is tr->lzi ? What is thorough insertion?
*/
boolean insertBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;
  int i;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  r = q->back;
  s = p->back;

  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];

  if(tr->thoroughInsertion)
  { 
    double  zqr[numBranches], zqs[numBranches], zrs[numBranches], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;
    double defaultArray[numBranches];
    double e1[numBranches], e2[numBranches], e3[numBranches];
    double *qz;

    qz = q->z;

    for(i = 0; i < numBranches; i++)
      defaultArray[i] = PLL_DEFAULTZ;

    makenewzGeneric(tr, pr, q, r, qz, PLL_ITERATIONS, zqr, PLL_FALSE);
    /* the branch lengths values will be estimated using q, r and s
     * q-s are not connected, but both q and s have a valid LH vector , so we can call makenewzGeneric  to get a value for
     * lzsum, which is then use to generate reasonable starting values e1, e2, e3 for the new branches we create after the       insertion
     */

    makenewzGeneric(tr, pr, q, s, defaultArray, PLL_ITERATIONS, zqs, PLL_FALSE);
    makenewzGeneric(tr, pr, r, s, defaultArray, PLL_ITERATIONS, zrs, PLL_FALSE);


    for(i = 0; i < numBranches; i++)
    {
      lzqr = (zqr[i] > PLL_ZMIN) ? log(zqr[i]) : log(PLL_ZMIN); 
      lzqs = (zqs[i] > PLL_ZMIN) ? log(zqs[i]) : log(PLL_ZMIN);
      lzrs = (zrs[i] > PLL_ZMIN) ? log(zrs[i]) : log(PLL_ZMIN);
      lzsum = 0.5 * (lzqr + lzqs + lzrs);

      lzq = lzsum - lzrs;
      lzr = lzsum - lzqs;
      lzs = lzsum - lzqr;
      lzmax = log(PLL_ZMAX);

      if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
      else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
      else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          

      e1[i] = exp(lzq);
      e2[i] = exp(lzr);
      e3[i] = exp(lzs);
    }
    hookup(p->next,       q, e1, numBranches);
    hookup(p->next->next, r, e2, numBranches);
    hookup(p,             s, e3, numBranches);      		  
  }
  else
  {       
    double  z[numBranches];

    for(i = 0; i < numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      

      if(z[i] < PLL_ZMIN) 
        z[i] = PLL_ZMIN;
      if(z[i] > PLL_ZMAX)
        z[i] = PLL_ZMAX;
    }

    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }

  newviewGeneric(tr, pr,p, PLL_FALSE);

  if(tr->thoroughInsertion)
  {     
    localSmooth(tr, pr, p, PLL_MAX_LOCAL_SMOOTHING_ITERATIONS);
    for(i = 0; i < numBranches; i++)
    {
      tr->lzq[i] = p->next->z[i];
      tr->lzr[i] = p->next->next->z[i];
      tr->lzs[i] = p->z[i];            
    }
  }           

  return  PLL_TRUE;
}

boolean insertRestoreBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;

  r = q->back;
  s = p->back;

  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, numBranches);
    hookup(p->next->next, r, tr->currentLZR, numBranches);
    hookup(p,             s, tr->currentLZS, numBranches);
  }
  else
  {       
    double  z[NUM_BRANCHES];
    int i;

    for(i = 0; i < numBranches; i++)
    {
      double zz;
      zz = sqrt(q->z[i]);     
      if(zz < PLL_ZMIN) 
        zz = PLL_ZMIN;
      if(zz > PLL_ZMAX)
        zz = PLL_ZMAX;
      z[i] = zz;
    }

    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }   

  newviewGeneric(tr, pr,p, PLL_FALSE);

  return  PLL_TRUE;
}


static void restoreTopologyOnly(pllInstance *tr, bestlist *bt, int numBranches)
{ 
  nodeptr p = tr->removeNode;
  nodeptr q = tr->insertNode;
  double qz[NUM_BRANCHES], pz[NUM_BRANCHES], p1z[NUM_BRANCHES], p2z[NUM_BRANCHES];
  nodeptr p1, p2, r, s;
  double currentLH = tr->likelihood;
  int i;

  p1 = p->next->back;
  p2 = p->next->next->back;

  //memcpy(p1z, p1->z, numBranches*sizeof(double));
  //memcpy(p2z, p2->z, numBranches*sizeof(double));
  //memcpy(qz, q->z, numBranches*sizeof(double));
  //memcpy(pz, p->z, numBranches*sizeof(double));
  for(i = 0; i < numBranches; i++)
  {
    p1z[i] = p1->z[i];
    p2z[i] = p2->z[i];
  }

  hookup(p1, p2, tr->currentZQR, numBranches);

  p->next->next->back = p->next->back = (node *) NULL;             
  for(i = 0; i < numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }

  r = q->back;
  s = p->back;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, numBranches);
    hookup(p->next->next, r, tr->currentLZR, numBranches);
    hookup(p,             s, tr->currentLZS, numBranches);
  }
  else
  { 	
    double  z[NUM_BRANCHES];	
    for(i = 0; i < numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      
      if(z[i] < PLL_ZMIN)
        z[i] = PLL_ZMIN;
      if(z[i] > PLL_ZMAX)
        z[i] = PLL_ZMAX;
    }
    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }     

  tr->likelihood = tr->bestOfNode;

  saveBestTree(bt, tr, numBranches);

  tr->likelihood = currentLH;

  hookup(q, r, qz, numBranches);

  p->next->next->back = p->next->back = (nodeptr) NULL;

  if(tr->thoroughInsertion)    
    hookup(p, s, pz, numBranches);

  hookup(p->next,       p1, p1z, numBranches);
  hookup(p->next->next, p2, p2z, numBranches);
}


boolean testInsertBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{

  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;
  boolean doIt = PLL_TRUE;
  double startLH = tr->endLH;
  int i;

  r = q->back; 
  for(i = 0; i < numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }



  if(doIt)
  {     
    if (! insertBIG(tr, pr, p, q))       return PLL_FALSE;

    evaluateGeneric(tr, pr, p->next->next, PLL_FALSE, PLL_FALSE);

    if(tr->likelihood > tr->bestOfNode)
    {
      tr->bestOfNode = tr->likelihood;
      tr->insertNode = q;
      tr->removeNode = p;   
      for(i = 0; i < numBranches; i++)
      {
        tr->currentZQR[i] = tr->zqr[i];           
        tr->currentLZR[i] = tr->lzr[i];
        tr->currentLZQ[i] = tr->lzq[i];
        tr->currentLZS[i] = tr->lzs[i];      
      }
    }

    if(tr->likelihood > tr->endLH)
    {			  
      tr->insertNode = q;
      tr->removeNode = p;   
      for(i = 0; i < numBranches; i++)
        tr->currentZQR[i] = tr->zqr[i];      
      tr->endLH = tr->likelihood;                      
    }        

    hookup(q, r, qz, numBranches);

    p->next->next->back = p->next->back = (nodeptr) NULL;

    if(tr->thoroughInsertion)
    {
      nodeptr s = p->back;
      hookup(p, s, pz, numBranches);
    } 

    if((tr->doCutoff) && (tr->likelihood < startLH))
    {
      tr->lhAVG += (startLH - tr->likelihood);
      tr->lhDEC++;
      if((startLH - tr->likelihood) >= tr->lhCutoff)
        return PLL_FALSE;	    
      else
        return PLL_TRUE;
    }
    else
      return PLL_TRUE;
  }
  else
    return PLL_TRUE;  
}




void addTraverseBIG(pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
  {              
    if (! testInsertBIG(tr, pr, p, q))  return;

  }

  if ((!isTip(q->number, tr->mxtips)) && (--maxtrav > 0)) 
  {    
    addTraverseBIG(tr, pr, p, q->next->back, mintrav, maxtrav);
    addTraverseBIG(tr, pr, p, q->next->next->back, mintrav, maxtrav);
  }
} 





int rearrangeBIG(pllInstance *tr, partitionList *pr, nodeptr p, int mintrav, int maxtrav)
{  
  double   p1z[NUM_BRANCHES], p2z[NUM_BRANCHES], q1z[NUM_BRANCHES], q2z[NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;  
  boolean doP = PLL_TRUE, doQ = PLL_TRUE;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
  q = p->back;




  if (!isTip(p->number, tr->mxtips) && doP) 
  {     
    p1 = p->next->back;
    p2 = p->next->next->back;


    if(!isTip(p1->number, tr->mxtips) || !isTip(p2->number, tr->mxtips))
    {
      for(i = 0; i < numBranches; i++)
      {
        p1z[i] = p1->z[i];
        p2z[i] = p2->z[i];	   	   
      }

      if (! removeNodeBIG(tr, pr, p,  numBranches)) return PLL_BADREAR;

      if (!isTip(p1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, p, p1->next->back,
            mintrav, maxtrav);         

        addTraverseBIG(tr, pr, p, p1->next->next->back,
            mintrav, maxtrav);          
      }

      if (!isTip(p2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, p, p2->next->back,
            mintrav, maxtrav);
        addTraverseBIG(tr, pr, p, p2->next->next->back,
            mintrav, maxtrav);          
      }

      hookup(p->next,       p1, p1z, numBranches);
      hookup(p->next->next, p2, p2z, numBranches);
      newviewGeneric(tr, pr,p, PLL_FALSE);
    }
  }  

  if (!isTip(q->number, tr->mxtips) && maxtrav > 0 && doQ) 
  {
    q1 = q->next->back;
    q2 = q->next->next->back;

    /*if (((!q1->tip) && (!q1->next->back->tip || !q1->next->next->back->tip)) ||
      ((!q2->tip) && (!q2->next->back->tip || !q2->next->next->back->tip))) */
    if (
        (
         ! isTip(q1->number, tr->mxtips) && 
         (! isTip(q1->next->back->number, tr->mxtips) || ! isTip(q1->next->next->back->number, tr->mxtips))
        )
        ||
        (
         ! isTip(q2->number, tr->mxtips) && 
         (! isTip(q2->next->back->number, tr->mxtips) || ! isTip(q2->next->next->back->number, tr->mxtips))
        )
       )
    {

      for(i = 0; i < numBranches; i++)
      {
        q1z[i] = q1->z[i];
        q2z[i] = q2->z[i];
      }

      if (! removeNodeBIG(tr, pr, q, numBranches)) return PLL_BADREAR;

      mintrav2 = mintrav > 2 ? mintrav : 2;

      if (/*! q1->tip*/ !isTip(q1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, q, q1->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, pr, q, q1->next->next->back,
            mintrav2 , maxtrav);         
      }

      if (/*! q2->tip*/ ! isTip(q2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, q, q2->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, pr, q, q2->next->next->back,
            mintrav2 , maxtrav);          
      }	   

      hookup(q->next,       q1, q1z, numBranches);
      hookup(q->next->next, q2, q2z, numBranches);

      newviewGeneric(tr, pr,q, PLL_FALSE);
    }
  } 

  return  1;
} 





static double treeOptimizeRapid(pllInstance *tr, partitionList *pr, int mintrav, int maxtrav, analdef *adef, bestlist *bt, infoList *iList)
{
  int i, index;   

  nodeRectifier(tr);



  if (maxtrav > tr->mxtips - 3)  
    maxtrav = tr->mxtips - 3;  



  resetInfoList(iList);

  resetBestTree(bt);

  tr->startLH = tr->endLH = tr->likelihood;

  if(tr->doCutoff)
  {
    if(tr->bigCutoff)
    {	  
      if(tr->itCount == 0)    
        tr->lhCutoff = 0.5 * (tr->likelihood / -1000.0);    
      else    		 
        tr->lhCutoff = 0.5 * ((tr->lhAVG) / ((double)(tr->lhDEC))); 	  
    }
    else
    {
      if(tr->itCount == 0)    
        tr->lhCutoff = tr->likelihood / -1000.0;    
      else    		 
        tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));   
    }    

    tr->itCount = tr->itCount + 1;
    tr->lhAVG = 0;
    tr->lhDEC = 0;
  }

  /*
     printf("DoCutoff: %d\n", tr->doCutoff);
     printf("%d %f %f %f\n", tr->itCount, tr->lhAVG, tr->lhDEC, tr->lhCutoff);

     printf("%d %d\n", mintrav, maxtrav);
     */

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
  {           
    tr->bestOfNode = PLL_UNLIKELY;          

    //if(adef->permuteTreeoptimize)
    //index = perm[i];
      //else
      index = i;     

    if(rearrangeBIG(tr, pr, tr->nodep[index], mintrav, maxtrav))
    {    
      if(tr->thoroughInsertion)
      {
        if(tr->endLH > tr->startLH)                 	
        {			   	     
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }
        else
        { 		  
          if(tr->bestOfNode != PLL_UNLIKELY)
            restoreTopologyOnly(tr, bt, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }	   
      }
      else
      {
        insertInfoList(tr->nodep[index], tr->bestOfNode, iList);	    
        if(tr->endLH > tr->startLH)                 	
        {		      
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  
        }	    	  
      }
    }     
  }     

  if(!tr->thoroughInsertion)
  {           
    tr->thoroughInsertion = PLL_TRUE;  

    for(i = 0; i < iList->valid; i++)
    { 	  
      tr->bestOfNode = PLL_UNLIKELY;

      if(rearrangeBIG(tr, pr, iList->list[i].node, mintrav, maxtrav))
      {	  
        if(tr->endLH > tr->startLH)                 	
        {	 	     
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }
        else
        { 

          if(tr->bestOfNode != PLL_UNLIKELY)
          {	     
            restoreTopologyOnly(tr, bt, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
          }	
        }      
      }
    }       

    tr->thoroughInsertion = PLL_FALSE;
  }

  //if(adef->permuteTreeoptimize)
  //  rax_free(perm);

  return tr->startLH;     
}




boolean testInsertRestoreBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{    
  if(tr->thoroughInsertion)
  {
    if (! insertBIG(tr, pr, p, q))       return PLL_FALSE;

    evaluateGeneric(tr, pr, p->next->next, PLL_FALSE, PLL_FALSE);
  }
  else
  {
    if (! insertRestoreBIG(tr, pr, p, q))       return PLL_FALSE;

    {
      nodeptr x, y;
      x = p->next->next;
      y = p->back;

      if(! isTip(x->number, tr->mxtips) && isTip(y->number, tr->mxtips))
      {
        while ((! x->x)) 
        {
          if (! (x->x))
            newviewGeneric(tr, pr,x, PLL_FALSE);
        }
      }

      if(isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! y->x)) 
        {		  
          if (! (y->x))
            newviewGeneric(tr, pr,y, PLL_FALSE);
        }
      }

      if(!isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! x->x) || (! y->x)) 
        {
          if (! (x->x))
            newviewGeneric(tr, pr,x, PLL_FALSE);
          if (! (y->x))
            newviewGeneric(tr, pr,y, PLL_FALSE);
        }
      }				      	

    }

    tr->likelihood = tr->endLH;
  }

  return PLL_TRUE;
} 

void restoreTreeFast(pllInstance *tr, partitionList *pr)
{
  removeNodeRestoreBIG(tr, pr, tr->removeNode);
  testInsertRestoreBIG(tr, pr, tr->removeNode, tr->insertNode);
}

static void myfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, stream);

  assert(bytes_written == nmemb);
}

static void myfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t
    bytes_read;

  bytes_read = fread(ptr, size, nmemb, stream);

  assert(bytes_read == nmemb);
}



/** @brief Write tree to file

    Serialize tree to a file. 

    @todo Document this
*/
static void writeTree(pllInstance *tr, FILE *f)
{
  int 
    x = tr->mxtips + 3 * (tr->mxtips - 1);

  nodeptr
    base = tr->nodeBaseAddress;

  myfwrite(&(tr->start->number), sizeof(int), 1, f);
  myfwrite(&base, sizeof(nodeptr), 1, f);
  myfwrite(tr->nodeBaseAddress, sizeof(node), x, f);

}

int ckpCount = 0;

/** @brief Write a checkpoint

    Is checkpoint enabled?

    @todo fill this up
*/
static void writeCheckpoint(pllInstance *tr, partitionList *pr, int state)
{
  int   
    model; 

  char 
    extendedName[2048],
    buf[64];

  FILE 
    *f;

  strcpy(extendedName,  binaryCheckpointName);
  strcat(extendedName, "_");
  sprintf(buf, "%d", ckpCount);
  strcat(extendedName, buf);  

  ckpCount++;

  f = myfopen(extendedName, "w"); 

  /* cdta */   


  tr->ckp.accumulatedTime = accumulatedTime + (gettime() - masterTime);

  tr->ckp.state = state;

  tr->ckp.tr_optimizeRateCategoryInvocations = tr->optimizeRateCategoryInvocations;
  tr->ckp.tr_thoroughInsertion = tr->thoroughInsertion;
  tr->ckp.tr_startLH  = tr->startLH;
  tr->ckp.tr_endLH    = tr->endLH;
  tr->ckp.tr_likelihood = tr->likelihood;
  tr->ckp.tr_bestOfNode = tr->bestOfNode;

  tr->ckp.tr_lhCutoff = tr->lhCutoff;
  tr->ckp.tr_lhAVG    = tr->lhAVG;
  tr->ckp.tr_lhDEC    = tr->lhDEC;     
  tr->ckp.tr_itCount  = tr->itCount;
  tr->ckp.tr_doCutoff = tr->doCutoff;
  /* printf("Acc time: %f\n", tr->ckp.accumulatedTime); */

  /* user stupidity */


  tr->ckp.searchConvergenceCriterion = tr->searchConvergenceCriterion;
  tr->ckp.rateHetModel =  tr->rateHetModel;
  tr->ckp.maxCategories =  tr->maxCategories;
  tr->ckp.NumberOfModels = pr->numberOfPartitions;
  tr->ckp.numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;
  tr->ckp.originalCrunchedLength = tr->originalCrunchedLength;
  tr->ckp.mxtips = tr->mxtips;
  strcpy(tr->ckp.seq_file, seq_file);

  /* handle user stupidity */


  myfwrite(&(tr->ckp), sizeof(checkPointState), 1, f);

  myfwrite(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfwrite(tr->tree1, sizeof(char), tr->treeStringLength, f);

  myfwrite(tr->rateCategory, sizeof(int), tr->originalCrunchedLength, f);
  myfwrite(tr->patrat, sizeof(double), tr->originalCrunchedLength, f);
  myfwrite(tr->patratStored, sizeof(double), tr->originalCrunchedLength, f);

  /* need to store this as well in checkpoints, otherwise the branch lengths 
     in the output tree files will be wrong, not the internal branch lengths though */

  //TODO: We have to change the way to store the fracchanges
  //myfwrite(tr->fracchanges,  sizeof(double), pr->numberOfPartitions, f);
  myfwrite(&(tr->fracchange),   sizeof(double), 1, f);


  /* pInfo */

  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    int 
      dataType = pr->partitionData[model]->dataType;

    myfwrite(&(pr->partitionData[model]->numberOfCategories), sizeof(int), 1, f);
    myfwrite(pr->partitionData[model]->perSiteRates, sizeof(double), tr->maxCategories, f);
    myfwrite(pr->partitionData[model]->EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfwrite(pr->partitionData[model]->EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfwrite(pr->partitionData[model]->EI, sizeof(double),  pLengths[dataType].eiLength, f);

    myfwrite(pr->partitionData[model]->frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfwrite(pr->partitionData[model]->tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);
    myfwrite(pr->partitionData[model]->substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);
    myfwrite(&(pr->partitionData[model]->alpha), sizeof(double), 1, f);
    
    if(pr->partitionData[model]->protModels == LG4)
	{
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {
	      myfwrite(pr->partitionData[model]->EIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	      myfwrite(pr->partitionData[model]->EV_LG4[k], sizeof(double),  pLengths[dataType].evLength, f);
	      myfwrite(pr->partitionData[model]->EI_LG4[k], sizeof(double),  pLengths[dataType].eiLength, f);    
	      myfwrite(pr->partitionData[model]->frequencies_LG4[k], sizeof(double),  pLengths[dataType].frequenciesLength, f);
	      myfwrite(pr->partitionData[model]->tipVector_LG4[k], sizeof(double),  pLengths[dataType].tipVectorLength, f);  
	      myfwrite(pr->partitionData[model]->substRates_LG4[k], sizeof(double),  pLengths[dataType].substRatesLength, f);    
	    }
	}
	
  }



  writeTree(tr, f);

  fclose(f); 

  printBothOpen("\nCheckpoint written to: %s likelihood: %f\n", extendedName, tr->likelihood);
}

static void readTree(pllInstance *tr, partitionList *pr, FILE *f)
{
  int 
    nodeNumber,   
    x = tr->mxtips + 3 * (tr->mxtips - 1);





  nodeptr
    startAddress;

  myfread(&nodeNumber, sizeof(int), 1, f);

  tr->start = tr->nodep[nodeNumber];

  /*printf("Start: %d %d\n", tr->start->number, nodeNumber);*/

  myfread(&startAddress, sizeof(nodeptr), 1, f);

  /*printf("%u %u\n", (size_t)startAddress, (size_t)tr->nodeBaseAddress);*/



  myfread(tr->nodeBaseAddress, sizeof(node), x, f);

  {
    int i;    

    size_t         
      offset;

    boolean 
      addIt;

    if(startAddress > tr->nodeBaseAddress)
    {
      addIt = PLL_FALSE;
      offset = (size_t)startAddress - (size_t)tr->nodeBaseAddress;
    }
    else
    {
      addIt = PLL_TRUE;
      offset = (size_t)tr->nodeBaseAddress - (size_t)startAddress;
    }       

    for(i = 0; i < x; i++)
    {      	
      if(addIt)
      {	    
        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next + offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back + offset);
      }
      else
      {

        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next - offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back - offset);	   
      } 
    }

  }

  evaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

  printBothOpen("RAxML Restart with likelihood: %1.50f\n", tr->likelihood);
}


static void readCheckpoint(pllInstance *tr, partitionList *pr)
{
  int  
    restartErrors = 0,
                  model; 

  FILE 
    *f = myfopen(binaryCheckpointInputName, "r");

  /* cdta */   

  myfread(&(tr->ckp), sizeof(checkPointState), 1, f);



  if(tr->ckp.searchConvergenceCriterion != tr->searchConvergenceCriterion)
  {
    printf("restart error, you are trying to re-start a run where the ML search criterion was turned %s\n", (tr->ckp.searchConvergenceCriterion)?"ON":"OFF");
    restartErrors++;
  }  

  if(tr->ckp.rateHetModel !=  tr->rateHetModel)
  {
    printf("restart error, you are trying to re-start a run with a different model of rate heterogeneity, the checkpoint was obtained under: %s\n", (tr->ckp.rateHetModel == GAMMA)?"GAMMA":"PSR");
    restartErrors++;
  }  

  if(tr->ckp.maxCategories !=  tr->maxCategories)
  {
    printf("restart error, you are trying to re-start a run with %d per-site rate categories, the checkpoint was obtained with: %d\n", tr->maxCategories, tr->ckp.maxCategories);
    restartErrors++;
  }

  if(tr->ckp.NumberOfModels != pr->numberOfPartitions)
  {
    printf("restart error, you are trying to re-start a run with %d partitions, the checkpoint was obtained with: %d partitions\n", (int)pr->numberOfPartitions, tr->ckp.NumberOfModels);
    restartErrors++;      
  }

  if(tr->ckp.numBranches != pr->perGeneBranchLengths?pr->numberOfPartitions:1)
  {
    printf("restart error, you are trying to re-start a run where independent per-site branch length estimates were turned %s\n", (tr->ckp.numBranches > 1)?"ON":"OFF");
    restartErrors++;
  }

  if(tr->ckp.originalCrunchedLength != tr->originalCrunchedLength)
  {
    printf("restart error, you are trying to re-start a run with %d site patterns, the checkpoint was obtained with: %d site patterns\n", tr->ckp.originalCrunchedLength, tr->originalCrunchedLength);
    restartErrors++; 
  }

  if(tr->ckp.mxtips != tr->mxtips)
  {
    printf("restart error, you are trying to re-start a run with %d taxa, the checkpoint was obtained with: %d taxa\n", tr->mxtips, tr->ckp.mxtips);
    restartErrors++; 
  }

  if(strcmp(tr->ckp.seq_file, seq_file) != 0)
  {
    printf("restart error, you are trying to re-start from alignemnt file %s, the checkpoint was obtained with file: %s\n", tr->ckp.seq_file, seq_file);
    restartErrors++; 
  }

  printf("REstart errors: %d\n", restartErrors);

  if(restartErrors > 0)
  {
    printf("User induced errors with the restart from checkpoint, exiting ...\n");

    if(restartErrors > 4)
      printf(" ... maybe you should do field work instead of trying to use a computer ...\n");
    if(restartErrors > 6)
      printf(" ... kala eisai telios ilithios;\n");

    exit(-1);
  }

  tr->ntips = tr->mxtips;

  tr->startLH    = tr->ckp.tr_startLH;
  tr->endLH      = tr->ckp.tr_endLH;
  tr->likelihood = tr->ckp.tr_likelihood;
  tr->bestOfNode = tr->ckp.tr_bestOfNode;

  tr->lhCutoff   = tr->ckp.tr_lhCutoff;
  tr->lhAVG      = tr->ckp.tr_lhAVG;
  tr->lhDEC      = tr->ckp.tr_lhDEC;
  tr->itCount    = tr->ckp.tr_itCount;
  tr->thoroughInsertion       = tr->ckp.tr_thoroughInsertion;



  accumulatedTime = tr->ckp.accumulatedTime;

  /* printf("Accumulated time so far: %f\n", accumulatedTime); */

  tr->optimizeRateCategoryInvocations = tr->ckp.tr_optimizeRateCategoryInvocations;


  myfread(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfread(tr->tree1, sizeof(char), tr->treeStringLength, f);

  if(tr->searchConvergenceCriterion)
  {
    int bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 0) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 0))
    { 

#ifdef _DEBUG_CHECKPOINTING    
      printf("parsing Tree 0\n");
#endif

      treeReadTopologyString(tr->tree0, tr);   

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 0, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, PLL_FALSE, PLL_FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }

    bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 1) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 1))
    {

#ifdef _DEBUG_CHECKPOINTING
      printf("parsing Tree 1\n");
#endif

      treeReadTopologyString(tr->tree1, tr); 

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 1, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, PLL_FALSE, PLL_FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }
  }

  myfread(tr->rateCategory, sizeof(int), tr->originalCrunchedLength, f);
  myfread(tr->patrat, sizeof(double), tr->originalCrunchedLength, f);
  myfread(tr->patratStored, sizeof(double), tr->originalCrunchedLength, f);


  /* need to read this as well in checkpoints, otherwise the branch lengths 
     in the output tree files will be wrong, not the internal branch lengths though */

  //TODO: Same problem as writing the checkpoint
  //myfread(tr->fracchanges,  sizeof(double), pr->numberOfPartitions, f);
  myfread(&(tr->fracchange),   sizeof(double), 1, f);

  /* pInfo */

  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    int 
      dataType = pr->partitionData[model]->dataType;

    myfread(&(pr->partitionData[model]->numberOfCategories), sizeof(int), 1, f);
    myfread(pr->partitionData[model]->perSiteRates, sizeof(double), tr->maxCategories, f);
    myfread(pr->partitionData[model]->EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfread(pr->partitionData[model]->EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfread(pr->partitionData[model]->EI, sizeof(double),  pLengths[dataType].eiLength, f);

    myfread(pr->partitionData[model]->frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfread(pr->partitionData[model]->tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);
    myfread(pr->partitionData[model]->substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);
    myfread(&(pr->partitionData[model]->alpha), sizeof(double), 1, f);
    
    if(pr->partitionData[model]->protModels == LG4)
	{
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {
	      myfread(pr->partitionData[model]->EIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	      myfread(pr->partitionData[model]->EV_LG4[k], sizeof(double),  pLengths[dataType].evLength, f);
	      myfread(pr->partitionData[model]->EI_LG4[k], sizeof(double),  pLengths[dataType].eiLength, f);    
	      myfread(pr->partitionData[model]->frequencies_LG4[k], sizeof(double),  pLengths[dataType].frequenciesLength, f);
	      myfread(pr->partitionData[model]->tipVector_LG4[k], sizeof(double),  pLengths[dataType].tipVectorLength, f);  
	      myfread(pr->partitionData[model]->substRates_LG4[k], sizeof(double),  pLengths[dataType].substRatesLength, f);    
	    }
	}

    pllMakeGammaCats(pr->partitionData[model]->alpha, pr->partitionData[model]->gammaRates, 4, tr->useMedian);
  }

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  masterBarrier(THREAD_COPY_INIT_MODEL, tr, pr);
#endif

  updatePerSiteRates(tr, pr, PLL_FALSE);

  readTree(tr, pr, f);

  fclose(f); 

}

static void restoreTreeDataValuesFromCheckpoint(pllInstance *tr)
{
  tr->optimizeRateCategoryInvocations = tr->ckp.tr_optimizeRateCategoryInvocations;  
  tr->thoroughInsertion = tr->ckp.tr_thoroughInsertion;
  tr->likelihood = tr->ckp.tr_likelihood;              
  tr->lhCutoff = tr->ckp.tr_lhCutoff;
  tr->lhAVG    = tr->ckp.tr_lhAVG;
  tr->lhDEC    = tr->ckp.tr_lhDEC;   	 
  tr->itCount = tr->ckp.tr_itCount;
  tr->doCutoff = tr->ckp.tr_doCutoff;
}

void restart(pllInstance *tr, partitionList *pr)
{  
  readCheckpoint(tr, pr);

  switch(tr->ckp.state)
  {
    case REARR_SETTING:      
      break;
    case FAST_SPRS:
      break;
    case SLOW_SPRS:
      break;
    default:
      assert(0);
  }
}

int determineRearrangementSetting(pllInstance *tr, partitionList *pr, analdef *adef, bestlist *bestT, bestlist *bt)
{
  const 
    int MaxFast = 26;

  int 
    i,   
    maxtrav = 5, 
    bestTrav = 5;

  double 
    startLH = tr->likelihood; 

  boolean 
    impr   = PLL_TRUE,
           cutoff = tr->doCutoff;

  if(adef->useCheckpoint)
  {
    assert(tr->ckp.state == REARR_SETTING);

    maxtrav = tr->ckp.maxtrav;
    bestTrav = tr->ckp.bestTrav;
    startLH  = tr->ckp.startLH;
    impr     = tr->ckp.impr;      
    cutoff = tr->ckp.cutoff;

    adef->useCheckpoint = PLL_FALSE;
  }

  tr->doCutoff = PLL_FALSE;      

  resetBestTree(bt);    

#ifdef _DEBUG_CHECKPOINTING
  //  printBothOpen("MAXTRAV: %d\n", maxtrav);
#endif

  while(impr && maxtrav < MaxFast)
  {	
    recallBestTree(bestT, 1, tr, pr);
    nodeRectifier(tr);                      

    {
      tr->ckp.cutoff = cutoff;	
      tr->ckp.maxtrav = maxtrav;
      tr->ckp.bestTrav = bestTrav;
      tr->ckp.startLH  = startLH;
      tr->ckp.impr = impr;

      //writeCheckpoint(tr, pr, REARR_SETTING);
    }

    if (maxtrav > tr->mxtips - 3)  
      maxtrav = tr->mxtips - 3;    

    tr->startLH = tr->endLH = tr->likelihood;

    for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {                	         
      tr->bestOfNode = PLL_UNLIKELY;

      if(rearrangeBIG(tr, pr, tr->nodep[i], 1, maxtrav))
      {	     
        if(tr->endLH > tr->startLH)                 	
        {		 	 	      
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;		 
        }	         	       	
      }
    }

    treeEvaluate(tr, pr, 8 ); // 32 * 0.25
    saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

#ifdef _DEBUG_CHECKPOINTING
    //printBothOpen("TRAV: %d lh %f\n", maxtrav, tr->likelihood);
#endif

    if(tr->likelihood > startLH)
    {	 
      startLH = tr->likelihood; 	  	  	  
      //printLog(tr);	  
      bestTrav = maxtrav;	 
      impr = PLL_TRUE;
    }
    else	
      impr = PLL_FALSE;	



    if(tr->doCutoff)
    {
      tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));       

      tr->itCount =  tr->itCount + 1;
      tr->lhAVG = 0;
      tr->lhDEC = 0;
    }

    maxtrav += 5;


  }

  recallBestTree(bt, 1, tr, pr);

  tr->doCutoff = cutoff; 

#ifdef _DEBUG_CHECKPOINTING
  //printBothOpen("BestTrav %d\n", bestTrav);
#endif

  return bestTrav;     
}





void computeBIGRAPID_Test (pllInstance *tr, partitionList *pr, boolean estimateModel)
{   
  int
    i,
    impr, 
    bestTrav = 0, 
    rearrangementsMax = 0, 
    rearrangementsMin = 0,    
    thoroughIterations = 0,
    fastIterations = 0;

  analdef 
    *adef = (analdef *)rax_malloc(sizeof(analdef) * 1);

  double 
    lh = PLL_UNLIKELY, 
    previousLh = PLL_UNLIKELY, 
    difference, 
    epsilon;              

  bestlist 
    *bestT, 
    *bt;    

  infoList 
    *iList = (infoList*)rax_malloc(sizeof(infoList));

 

  adef->initialSet = PLL_FALSE;
  adef->useCheckpoint = PLL_FALSE;
  adef->stepwidth = 5;
  adef->max_rearrange = 21;
  

  bestT = (bestlist *) rax_malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);

  bt = (bestlist *) rax_malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  
  initInfoList(iList, 50); 

  difference = 10.0;
  epsilon = 0.01;    

  tr->thoroughInsertion = PLL_FALSE;     

  
  assert(estimateModel);

  modOpt(tr, pr, 10.0);

  printf("mo opt 1 %f\n", tr->likelihood);

  saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);  

  if(!adef->initialSet)   
    {
      if((!adef->useCheckpoint) || (adef->useCheckpoint && tr->ckp.state == REARR_SETTING))
	{
	  printf("K\n");
	  bestTrav = adef->bestTrav = determineRearrangementSetting(tr, pr, adef, bestT, bt);
	  printf("\nBest rearrangement radius: %d\n", bestTrav);
	}
    }
  else
    {
      bestTrav = adef->bestTrav = adef->initial;       
      printf("\nUser-defined rearrangement radius: %d\n", bestTrav);
    }
 
  modOpt(tr, pr, 5.0);
  
  printf("modopt %f\n", tr->likelihood);

  /* save the current tree again, while the topology has not changed, the branch lengths have changed in the meantime, hence
     we need to store them again */

  saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  /* set the loop variable to PLL_TRUE */

  impr = 1; 

  if(tr->doCutoff)
    tr->itCount = 0;

  while(impr)
    {    
      printf("Fast loop\n");

      recallBestTree(bestT, 1, tr, pr);        

      /* count how many fast iterations with so-called fast SPR moves we have executed */
      
      fastIterations++;	            
      
      treeEvaluate(tr, pr, 32);  // 32 * 1 = 32
     
      saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
    
      lh = previousLh = tr->likelihood;     

      treeOptimizeRapid(tr, pr, 1, bestTrav, adef, bt, iList);    

      impr = 0;    

      for(i = 1; i <= bt->nvalid; i++)
	{	    	
	  /* restore tree i from list generated by treeOptimizeRapid */
	  
	  recallBestTree(bt, i, tr, pr);
	  
	  /* optimize branch lengths of this tree */
	  
	  treeEvaluate(tr, pr, 8); // 0.25 * 32
	  
	  /* calc. the likelihood improvement */
	  
	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  
	  /* if the likelihood has improved save the current tree as best tree and continue */
	  /* note that we always compre this tree to the likelihood of the previous best tree */
	  
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;
	      printf("%f \n", tr->likelihood);
	      saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);	      
	    }	   	   
	}
    }

  
  tr->thoroughInsertion = PLL_TRUE;
  impr = 1;

  recallBestTree(bestT, 1, tr, pr);
  
  evaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
  
  modOpt(tr, pr, 1.0);


  while(1)
  {	 
    printf("slow loop\n");
    recallBestTree(bestT, 1, tr, pr);    

    if(impr)
      {
	
	rearrangementsMin = 1;      
	rearrangementsMax = adef->stepwidth;       
	thoroughIterations++;	  
      }			  			
    else
      {	  	  
	rearrangementsMax += adef->stepwidth;
	rearrangementsMin += adef->stepwidth; 	        	      
	
	/* if we have already tried them then abandon this loop, the search has converged */
	if(rearrangementsMax > adef->max_rearrange)	     	     	 
	  goto cleanup; 	   
      }
    
   

    treeEvaluate(tr, pr, 32 ); // 32 * 1

   
    previousLh = lh = tr->likelihood;	      
    saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
   

    treeOptimizeRapid(tr, pr, rearrangementsMin, rearrangementsMax, adef, bt, iList);

    impr = 0;			      		            

   

    for(i = 1; i <= bt->nvalid; i++)
    {		 
      recallBestTree(bt, i, tr, pr);

      treeEvaluate(tr, pr, 8); // 0.25	* 32

      difference = ((tr->likelihood > previousLh)? 
          tr->likelihood - previousLh: 
          previousLh - tr->likelihood); 	    
      if(tr->likelihood > lh && difference > epsilon)
      {
        impr = 1;	       
        lh = tr->likelihood;	  	     
        saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
      }	   	   
    }  


  }

cleanup: 

  /* do a final full tree traversal, not sure if this is required here */


  evaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);


  freeBestTree(bestT);
  rax_free(bestT);
  freeBestTree(bt);
  rax_free(bt);
  freeInfoList(iList);  
  rax_free(iList);

  rax_free(adef);

 

}



/* The number of maximum smoothing iterations is given explicitely */
boolean 
treeEvaluate (pllInstance *tr, partitionList *pr, int maxSmoothIterations)       /* Evaluate a user tree */
{
  smoothTree(tr, pr, maxSmoothIterations); /* former (32 * smoothFactor) */

  evaluateGeneric(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);

  return PLL_TRUE;
}

/** @brief Perform an NNI move

    Modify the tree topology of instance \a tr by performing an NNI (Neighbour Neighbor
    Interchange) move at node \a p. Perform one of the two possible NNI moves
    based on whether \a swap is set to 1 or 2.
*/
void NNI(pllInstance * tr, nodeptr p, int swap)
{
  nodeptr       q, tmp;

  q = p->back;
  assert(!isTip(q->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));


  if(swap == 1)
   {
     tmp = p->next->back;
     hookupFull(p->next, q->next->back, q->next->z);
     hookupFull(q->next, tmp,           p->next->z);
   }
  else
   {
      tmp = p->next->next->back;
      hookupFull(p->next->next, q->next->back, q->next->z);
      hookupFull(q->next,       tmp,           p->next->next->z);
   }
}
