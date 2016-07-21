/******************************************** genericParallelization.c ******************************************************/
#ifdef MEASURE_TIME_PARALLEL
static void reduceTimesWorkerRegions(pllInstance *tr, double *mins, double *maxs)
{
  int tid = tr->threadID; 
  int i,j ; 
  double reduction[NUM_PAR_JOBS * tr->numberOfThreads]; 

  ASSIGN_GATHER(reduction, timeBuffer, NUM_PAR_JOBS, DOUBLE, tr->threadID); 

#ifdef _USE_PTHREADS
  /* we'd need a proper barrier here... this evaluation is mostly interesting for MPI  */
  printf("\n\ncomment out MEASURE_TIME_PARALLEL\n\n");   
  assert(0); 
#else 
   MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* find min and max time */
  if(MASTER_P)
    {      
      for(j = 0; j < NUM_PAR_JOBS; ++j)
	{
	  boolean isFirst = PLL_TRUE; 
	  for(i = 0; i < tr->numberOfThreads; ++i)
	    {	      
	      double num = timeBuffer[i * NUM_PAR_JOBS + j]; 
	      if(isFirst || num < mins[j])
		mins[j] = num; 
	      if(isFirst || num > maxs[j])
		maxs[j] = num; 
	      isFirst = PLL_FALSE; 
	    }
	}	
    }  
}

static void printParallelTimePerRegion(double *mins, double *maxs)
{
  int i; 
  double allTime = 0; 
  double relTime[NUM_PAR_JOBS+1]; 
  for(i = 0; i < NUM_PAR_JOBS+1; ++i)
    allTime += timePerRegion[i]; 
  for(i = 0; i < NUM_PAR_JOBS+1; ++i)
    relTime[i] = (timePerRegion[i] / allTime) * 100 ; 

  printf("\n\nTime spent per region \nmasterTimeAbs\tmasterTimeRel\tloadBalance\tcommOverhead\n"); 
  for(i = 0; i < NUM_PAR_JOBS; ++i )
    if(timePerRegion[i] != 0)
      printf("%f\t%.2f\%\t%.2f\%\t%.2f\%\t%s\n", timePerRegion[i], relTime[i], (maxs[i] - mins[i]) * 100  / maxs[i] , (maxs[i] - timePerRegion[i]) * 100  / timePerRegion[i], getJobName(i) ); 
  printf("================\n%f\t%.2f\%\tSEQUENTIAL\n", timePerRegion[NUM_PAR_JOBS], relTime[i]); 
  printf("loadbalance: (minT - maxT) / maxT, \twhere minT is the time the fastest worker took for this region (maxT analogous) \n"); 
  printf("commOverhead: (maxWorker - masterTime) / masterTime, \t where maxWorker is the time the slowest worker spent in this region and masterTime is the time the master spent in this region (e.g., doing additional reduction stuff)\n"); 
}
#endif




/******************************************** utils.c ***********************************************************************/

void init_default(pllInstance *tr)
{

  /*********** tr inits **************/

  tr->numberOfThreads = 1; 
  tr->doCutoff = PLL_TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = PLL_FALSE;
  tr->rateHetModel = GAMMA;

  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->saveMemory = PLL_FALSE;

  tr->fastScaling = PLL_FALSE;

  tr->manyPartitions = PLL_FALSE;

  tr->startingTree = randomTree;

  tr->categories             = 25;

  tr->grouped = PLL_FALSE;
  tr->constrained = PLL_FALSE;

  tr->gapyness               = 0.0; 
  tr->useMedian = PLL_FALSE;
  /* recom */
  tr->useRecom = PLL_FALSE;
  tr->rvec = (recompVectors*)NULL;
  /* recom */

  /********* tr inits end*************/

}

void getDataTypeString(pllInstance *tr, pInfo *partitionInfo, char typeOfData[1024])
{
  switch(partitionInfo->dataType)
  {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
  }
}

void printResult(pllInstance *tr, partitionList *pr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "";

  strcpy(temporaryFileName, resultFileName);

  switch(adef->mode)
  {    
    case TREE_EVALUATION:
      pllTreeToNewick(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

      logFile = myfopen(temporaryFileName, "wb");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);

      if(adef->perGeneBranchLengths)
        printTreePerGene(tr, pr, adef, temporaryFileName, "wb");
      break;
    case BIG_RAPID_MODE:
      if(finalPrint)
      {
        switch(tr->rateHetModel)
        {
          case GAMMA:
          case GAMMA_I:

            pllTreeToNewick(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint,
                PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

            logFile = myfopen(temporaryFileName, "wb");
            fprintf(logFile, "%s", tr->tree_string);
            fclose(logFile);

            if(adef->perGeneBranchLengths)
              printTreePerGene(tr, pr, adef, temporaryFileName, "wb");
            break;
          case CAT:
            /*pllTreeToNewick(tr->tree_string, tr, pr, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint, adef,
              PLL_NO_BRANCHES, PLL_FALSE, PLL_FALSE);*/


            pllTreeToNewick(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE,
                PLL_TRUE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);




            logFile = myfopen(temporaryFileName, "wb");
            fprintf(logFile, "%s", tr->tree_string);
            fclose(logFile);

            break;
          default:
            assert(0);
        }
      }
      else
      {
        pllTreeToNewick(tr->tree_string, tr, pr, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint,
            PLL_NO_BRANCHES, PLL_FALSE, PLL_FALSE);
        logFile = myfopen(temporaryFileName, "wb");
        fprintf(logFile, "%s", tr->tree_string);
        fclose(logFile);
      }    
      break;
    default:
      printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
      exit(-1);
      break;
  }
}




/* removed the static keyword for using this function in the examples */
static boolean setupTree (pllInstance *tr, boolean doInit, partitionList *partitions)
{
  nodeptr  p0, p, q;
  int
    i,
    j;

  int
    tips,
    inter; 

  if(doInit)
    init_default(tr);

  tr->bigCutoff = PLL_FALSE;

  tr->maxCategories = MAX(4, tr->categories);

  tips  = (size_t)tr->mxtips;
  inter = (size_t)(tr->mxtips - 1);

  tr->treeStringLength = tr->mxtips * (PLL_NMLNGTH + 128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].executeModel = (boolean *)rax_malloc(sizeof(boolean) * (size_t)NUM_BRANCHES);
  tr->td[0].parameterValues = (double *)rax_malloc(sizeof(double) * (size_t)NUM_BRANCHES);

  tr->fracchange = -1.0;

  tr->constraintVector = (int *)rax_malloc((2 * (size_t)tr->mxtips) * sizeof(int));

  tr->nameList = (char **)rax_malloc(sizeof(char *) * (tips + 1));


  p0 = (nodeptr)rax_malloc((tips + 3 * inter) * sizeof(node));
  assert(p0);

  tr->nodeBaseAddress = p0;


  tr->nodep = (nodeptr *) rax_malloc((2* (size_t)tr->mxtips) * sizeof(nodeptr));
  assert(tr->nodep);    

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
  {
    p = p0++;

    p->hash   =  KISS32(); /* hast table stuff */
    p->x      =  0;
    p->xBips  = 0;
    p->number =  i;
    p->next   =  p;
    p->back   = (node *)NULL;
    p->bInf   = (branchInfo *)NULL;            
    tr->nodep[i] = p;
  }

  for (i = tips + 1; i <= tips + inter; i++)
  {
    q = (node *) NULL;
    for (j = 1; j <= 3; j++)
    {	 
      p = p0++;
      if(j == 1)
      {
        p->xBips = 1;
        p->x = 1;
      }
      else
      {
        p->xBips = 0;
        p->x =  0;
      }
      p->number = i;
      p->next   = q;
      p->bInf   = (branchInfo *)NULL;
      p->back   = (node *) NULL;
      p->hash   = 0;       
      q = p;
    }
    p->next->next->next = p;
    tr->nodep[i] = p;
  }

  tr->likelihood  = PLL_UNLIKELY;
  tr->start       = (node *) NULL;  

  tr->ntips       = 0;
  tr->nextnode    = 0;

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionSmoothed[i] = PLL_FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;

  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  for (i = 0; i < partitions->numberOfPartitions; i++) {
	partitions->partitionData[i] = (pInfo*)rax_malloc (sizeof(pInfo));
	partitions->partitionData[i]->partitionContribution = -1.0;
	partitions->partitionData[i]->partitionLH = 0.0;
	partitions->partitionData[i]->fracchange = 1.0;
  }

  return PLL_TRUE;
}


boolean modelExists(char *model, pllInstance *tr)
{
  /********** BINARY ********************/

  if(strcmp(model, "PSR") == 0)
  {
    tr->rateHetModel = CAT;
    return PLL_TRUE;
  }

  if(strcmp(model, "GAMMA") == 0)
  {
    tr->rateHetModel = GAMMA;
    return PLL_TRUE;
  }


  return PLL_FALSE;
}




/******************************************** treeIO.c ***********************************************************************/
void printTreePerGene(pllInstance *tr, partitionList *pr, analdef *adef, char *fileName, char *permission)
{  
  FILE *treeFile;
  char extendedTreeFileName[1024];
  char buf[16];
  int i;

  assert(adef->perGeneBranchLengths);

  int numberOfModels = pr->perGeneBranchLengths?pr->numberOfPartitions:1;
  for(i = 0; i < numberOfModels; i++)
    {
      strcpy(extendedTreeFileName, fileName);
      sprintf(buf,"%d", i);
      strcat(extendedTreeFileName, ".PARTITION.");
      strcat(extendedTreeFileName, buf);
      /*printf("Partitiuon %d file %s\n", i, extendedTreeFileName);*/
      pllTreeToNewick(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_TRUE, i, PLL_FALSE, PLL_FALSE);
      treeFile = myfopen(extendedTreeFileName, permission);
      fprintf(treeFile, "%s", tr->tree_string);
      fclose(treeFile);
    }  
    
}

void printTopology(pllInstance *tr, partitionList *pr, boolean printInner)
{
  if(!printInner)
  {
    boolean printBranchLengths = PLL_FALSE;
    pllTreeToNewick(tr->tree_string, tr, pr, tr->start->back, printBranchLengths, 0, 0, 0, 0, PLL_SUMMARIZE_LH, 0,0);
    fprintf(stderr, "%s", tr->tree_string);
  }
  else
  {
    TreeInner2StringREC(tr->tree_string, tr, pr, tr->start->back, PLL_FALSE, 0, 0, 0, 0, PLL_SUMMARIZE_LH, 0,0, PLL_TRUE);
    fprintf(stderr, "%s", tr->tree_string);
    fprintf(stderr, "Start was %d, pnb %d, pnnb %d, pback %d\n",
        tr->start->back->number,
        tr->start->back->next->back->number,
        tr->start->back->next->next->back->number,
        tr->start->number);
  }
}

void parsimonySPR(nodeptr p, partitionList *pr, pllInstance *tr)
{
  int i;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  double
    p1z[PLL_NUM_BRANCHES],
    p2z[PLL_NUM_BRANCHES];

  nodeptr
    p1 = p->next->back,
    p2 = p->next->next->back;

  //unsigned int score = evaluateParsimony(tr, pr, p, PLL_TRUE);

  //printf("parsimonyScore: %u\n", score);

  for(i = 0; i < numBranches; i++)
    {
      p1z[i] = p1->z[i];
      p2z[i] = p2->z[i];
    }

  tr->bestParsimony = INT_MAX;

  hookupDefault(p1, p2);

  p->next->next->back = p->next->back = (node *) NULL;

  if (p1->number > tr->mxtips)
    {
      addTraverseParsimony(tr, pr, p, p1->next->back, 0, 0, PLL_TRUE, PLL_TRUE);
      addTraverseParsimony(tr, pr, p, p1->next->next->back, 0, 0, PLL_TRUE, PLL_TRUE);
    }

  if(p2->number > tr->mxtips)
    {
      addTraverseParsimony(tr, pr, p, p2->next->back, 0, 0, PLL_TRUE, PLL_TRUE);
      addTraverseParsimony(tr, pr, p, p2->next->next->back, 0, 0, PLL_TRUE, PLL_TRUE);
    }

  //printf("best %u nodes %d %d\n",tr->bestParsimony, tr->insertNode->number, tr->insertNode->back->number);

  hookup(p1, p->next, p1z,       numBranches);
  hookup(p2, p->next->next, p2z, numBranches);
}

/******************************************** searchAlgo.c ***********************************************************************/

#ifdef _COMPUTEBIGRAPID
int determineRearrangementSetting(pllInstance *tr, partitionList *pr, bestlist *bestT, bestlist *bt)
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

  if(tr->useCheckpoint)
  {
    assert(tr->ckp.state == REARR_SETTING);

    maxtrav = tr->ckp.maxtrav;
    bestTrav = tr->ckp.bestTrav;
    startLH  = tr->ckp.startLH;
    impr     = tr->ckp.impr;      
    cutoff = tr->ckp.cutoff;

    tr->useCheckpoint = PLL_FALSE;
  }

  tr->doCutoff = PLL_FALSE;      

  resetBestTree(bt);    

#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("MAXTRAV: %d\n", maxtrav);
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

      writeCheckpoint(tr, pr, REARR_SETTING);
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
    printBothOpen("TRAV: %d lh %f\n", maxtrav, tr->likelihood);
#endif

    if(tr->likelihood > startLH)
    {	 
      startLH = tr->likelihood; 	  	  	  
      printLog(tr);	  
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
  printBothOpen("BestTrav %d\n", bestTrav);
#endif

  return bestTrav;     
}
#endif



#ifdef _COMPUTEBIGRAPID
/** @brief Hill climbing search algorithm for exploring the tree space

    Detailed description coming
    
    @param tr
      Tree instance

    @param pr
      List of partitions

    @todo
      Remove adef

    @param estimateModel
      Dont know yet

*/
void computeBIGRAPID (pllInstance *tr, partitionList *pr, analdef *adef, boolean estimateModel)
{   
  int
    i,
    impr, 
    bestTrav = 0, 
    rearrangementsMax = 0, 
    rearrangementsMin = 0,    
    thoroughIterations = 0,
    fastIterations = 0;

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

  /* now here is the RAxML hill climbing search algorithm */


  /* initialize two lists of size 1 and size 20 that will keep track of the best 
     and 20 best tree topologies respectively */

  bestT = (bestlist *) rax_malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);

  bt = (bestlist *) rax_malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  /* initialize an additional data structure used by the search algo, all of this is pretty 
     RAxML-specific and should probably not be in the library */

  initInfoList(iList, 50);

  /* some pretty atbitrary thresholds */

  difference = 10.0;
  epsilon = 0.01;    

  /* Thorough = 0 means that we will do fast SPR inbsertions without optimizing the 
     three branches adjacent to the subtree insertion position via Newton-Raphson 
     */

  tr->thoroughInsertion = PLL_FALSE;     

  /* if we are not using a checkpoint and estimateModel is set to PLL_TRUE we call the function 
     that optimizes model parameters, such as the CAT model assignment, the alpha paremeter
     or the rates in the GTR matrix. Otherwise we just optimize the branch lengths. Note that 
     the second parameter of treeEvaluate() controls how many times we will iterate over all branches 
     of the tree until we give up, provided that, the br-len opt. has not converged before.
     */

  if(!adef->useCheckpoint)
  {
    if(estimateModel)
      modOpt(tr, pr, 10.0);
    else
      treeEvaluate(tr, pr, 64); // 32 * 2
  }

  /* print some stuff to the RAxML_log file */

  printLog(tr); 

  /* save the current tree (which is the input tree parsed via -t in the bestT list */

  saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  /* if the rearrangmenet radius has been set by the user ie. adef->initailSet == PLL_TRUE 
     then just set the apppropriate parameter.
     Otherwise, call the function  determineRearrangementSetting() that seeks 
     for the best radius by executing SPR moves on the initial tree with different radii
     and returns the smallest radius that yields the best log likelihood score after 
     applying one cycle of SPR moves to the tree 
     */

  if(!adef->initialSet)   
  {
    if((!adef->useCheckpoint) || (adef->useCheckpoint && tr->ckp.state == REARR_SETTING))
    {
      bestTrav = adef->bestTrav = determineRearrangementSetting(tr, pr, adef, bestT, bt);
      printBothOpen("\nBest rearrangement radius: %d\n", bestTrav);
    }
  }
  else
  {
    bestTrav = adef->bestTrav = adef->initial;       
    printBothOpen("\nUser-defined rearrangement radius: %d\n", bestTrav);
  }


  /* some checkpointing noise */
  if(!(adef->useCheckpoint && (tr->ckp.state == FAST_SPRS || tr->ckp.state == SLOW_SPRS)))
  {      

    /* optimize model params more thoroughly or just optimize branch lengths */
    if(estimateModel)
      modOpt(tr, pr, 5.0);
    else
      treeEvaluate(tr, pr, 32);   // 32 * 1
  }

  /* save the current tree again, while the topology has not changed, the branch lengths have changed in the meantime, hence
     we need to store them again */

  saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  /* set the loop variable to PLL_TRUE */

  impr = 1;

  /* this is for the additional RAxML heuristics described imn this paper here:

     A. Stamatakis,  F. Blagojevic, C.D. Antonopoulos, D.S. Nikolopoulos: "Exploring new Search Algorithms and Hardware for Phylogenetics: RAxML meets the IBM Cell". 
     In Journal of VLSI Signal Processing Systems, 48(3):271-286, 2007.

     This is turned on by default 
     */


  if(tr->doCutoff)
    tr->itCount = 0;

  /* figure out where to continue computations if we restarted from a checkpoint */

  if(adef->useCheckpoint && tr->ckp.state == FAST_SPRS)
    goto START_FAST_SPRS;

  if(adef->useCheckpoint && tr->ckp.state == SLOW_SPRS)
    goto START_SLOW_SPRS;

  while(impr)
  {
START_FAST_SPRS:
    /* if re-starting from checkpoint set the required variable values to the 
       values that they had when the checkpoint was written */

    if(adef->useCheckpoint && tr->ckp.state == FAST_SPRS)
    {	    
      impr = tr->ckp.impr;	  
      bestTrav = tr->ckp.bestTrav;	  
      rearrangementsMax = tr->ckp.rearrangementsMax;
      rearrangementsMin = tr->ckp.rearrangementsMin;
      thoroughIterations = tr->ckp.thoroughIterations;
      fastIterations = tr->ckp.fastIterations;  
      lh = tr->ckp.lh;
      previousLh = tr->ckp.previousLh;
      difference = tr->ckp.difference;
      epsilon    = tr->ckp.epsilon;  

      restoreTreeDataValuesFromCheckpoint(tr);	  

      adef->useCheckpoint = PLL_FALSE;
    }
    else
      /* otherwise, restore the currently best tree */
      recallBestTree(bestT, 1, tr, pr);

    /* save states of algorithmic/heuristic variables for printing the next checkpoint */


    tr->ckp.impr = impr;	
    tr->ckp.bestTrav = bestTrav;      
    tr->ckp.rearrangementsMax = rearrangementsMax;
    tr->ckp.rearrangementsMin = rearrangementsMin;
    tr->ckp.thoroughIterations = thoroughIterations;
    tr->ckp.fastIterations = fastIterations;  
    tr->ckp.lh = lh;
    tr->ckp.previousLh = previousLh;
    tr->ckp.difference = difference;
    tr->ckp.epsilon    = epsilon;              
    tr->ckp.bestTrav = bestTrav;       
    tr->ckp.impr = impr;                 

    /* write a binary checkpoint */
    writeCheckpoint(tr, pr, FAST_SPRS);

    /* this is the aforementioned convergence criterion that requires computing the RF,
       let's not worry about this right now */

    if(tr->searchConvergenceCriterion)
    {
      int bCounter = 0;	  	      	 	  	  	

      if(fastIterations > 1)
        cleanupHashTable(tr->h, (fastIterations % 2));		
      
      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, fastIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, PLL_FALSE, PLL_FALSE, tr->threadID);	    

      {
        char 
          *buffer = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
#ifdef _DEBUG_CHECKPOINTING
        printf("Storing tree in slot %d\n", fastIterations % 2);
#endif

        pllTreeToNewick(buffer, tr, pr, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

        if(fastIterations % 2 == 0)	      
          memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
        else
          memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    

        rax_free(buffer);
      }


      assert(bCounter == tr->mxtips - 3);	    	   

      if(fastIterations > 0)
      {
        double rrf = convergenceCriterion(tr->h, tr->mxtips);

        if(rrf <= 0.01) /* 1% cutoff */
        {
          printBothOpen("ML fast search converged at fast SPR cycle %d with stopping criterion\n", fastIterations);
          printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
          cleanupHashTable(tr->h, 0);
          cleanupHashTable(tr->h, 1);
          goto cleanup_fast;
        }
        else		    
          printBothOpen("ML search convergence criterion fast cycle %d->%d Relative Robinson-Foulds %f\n", fastIterations - 1, fastIterations, rrf);
      }
    }


    /* count how many fast iterations with so-called fast SPR moves we have executed */

    fastIterations++;	

    /* optimize branch lengths */

    treeEvaluate(tr, pr, 32);  // 32 * 1 = 32

    /* save the tree with those branch lengths again */

    saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

    /* print the log likelihood */

    printLog(tr);    

    /* print this intermediate tree to file */

    printResult(tr, pr, adef, PLL_FALSE);

    /* update the current best likelihood */

    lh = previousLh = tr->likelihood;

    /* in here we actually do a cycle of SPR moves */

    treeOptimizeRapid(tr, pr, 1, bestTrav, bt, iList);

    /* set impr to 0 since in the immediately following for loop we check if the SPR moves above have generated 
       a better tree */

    impr = 0;

    /* loop over the 20 best trees generated by the fast SPR moves, and check if they improve the likelihood after all of their branch lengths
       have been optimized */

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
        saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

      }
    }
#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("FAST LH: %f\n", lh);
#endif


  }

  /* needed for this RF-based convergence criterion that I actually describe in here:

     A. Stamatakis: "Phylogenetic Search Algorithms for Maximum Likelihood". In M. Elloumi, A.Y. Zomaya, editors. 
     Algorithms in Computational Biology: techniques, Approaches and Applications, John Wiley and Sons

     a copy of this book is in my office */

  if(tr->searchConvergenceCriterion)
  {
    cleanupHashTable(tr->h, 0);
    cleanupHashTable(tr->h, 1);
  }

cleanup_fast:  
  /*
     now we have jumped out of the loop that executes 
     fast SPRs, and next we will execute a loop that executes thorough SPR cycles (with SPR moves 
     that optimize via newton-Raphson all adjacent branches to the insertion point) 
     until no thorough SPR move can be found that improves the likelihood further. A classic 
     hill climbing algo.
     */

  tr->thoroughInsertion = PLL_TRUE;
  impr = 1;

  /* restore the currently best tree. this si actually required, because we do not know which tree
     is actually stored in the tree data structure when the above loop exits */

  recallBestTree(bestT, 1, tr, pr);

  {
    /* RE-TRAVERSE THE ENTIRE TREE */

    evaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("After Fast SPRs Final %f\n", tr->likelihood);   
#endif
  }

  /* optimize model params (including branch lengths) or just 
     optimize branch lengths and leave the other model parameters (GTR rates, alhpa) 
     alone */

  if(estimateModel)
    modOpt(tr, pr, 1.0);
  else
    treeEvaluate(tr, pr, 32 ); //32 * 1

  /* start loop that executes thorough SPR cycles */

  while(1)
  {	 
    /* once again if we want to restart from a checkpoint that was written during this loop we need
       to restore the values of the variables appropriately */
START_SLOW_SPRS:
    if(adef->useCheckpoint && tr->ckp.state == SLOW_SPRS)
    {	        
      impr = tr->ckp.impr;	 
      bestTrav = tr->ckp.bestTrav;	 
      rearrangementsMax = tr->ckp.rearrangementsMax;
      rearrangementsMin = tr->ckp.rearrangementsMin;
      thoroughIterations = tr->ckp.thoroughIterations;
      fastIterations = tr->ckp.fastIterations;     
      lh = tr->ckp.lh;
      previousLh = tr->ckp.previousLh;
      difference = tr->ckp.difference;
      epsilon    = tr->ckp.epsilon;                    

      restoreTreeDataValuesFromCheckpoint(tr);	  

      adef->useCheckpoint = PLL_FALSE;
    }
    else
      /* otherwise we restore the currently best tree and load it from bestT into our tree data 
         structuire tr */
      recallBestTree(bestT, 1, tr, pr);

    /* now, we write a checkpoint */

    tr->ckp.impr = impr;      
    tr->ckp.bestTrav = bestTrav;      
    tr->ckp.rearrangementsMax = rearrangementsMax;
    tr->ckp.rearrangementsMin = rearrangementsMin;
    tr->ckp.thoroughIterations = thoroughIterations;
    tr->ckp.fastIterations = fastIterations;	  
    tr->ckp.lh = lh;
    tr->ckp.previousLh = previousLh;
    tr->ckp.difference = difference;
    tr->ckp.epsilon    = epsilon;              
    tr->ckp.bestTrav = bestTrav;       
    tr->ckp.impr = impr;                 

    /* write binary checkpoint to file */
    writeCheckpoint(tr, pr, SLOW_SPRS);

    if(impr)
    {
      /* if the logl has improved write out some stuff and adapt the rearrangement radii */
      printResult(tr, pr, adef, PLL_FALSE);
      /* minimum rearrangement radius */
      rearrangementsMin = 1;
      /* max radius, this is probably something I need to explain at the whiteboard */
      rearrangementsMax = adef->stepwidth;	

      /* once again the convergence criterion */

      if(tr->searchConvergenceCriterion)
      {
        int bCounter = 0;	      

        if(thoroughIterations > 1)
          cleanupHashTable(tr->h, (thoroughIterations % 2));		

        bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, thoroughIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
            &bCounter, 1, PLL_FALSE, PLL_FALSE, tr->threadID);	    


        {
          char 
            *buffer = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));

#ifdef _DEBUG_CHECKPOINTING		
          printf("Storing tree in slot %d\n", thoroughIterations % 2);
#endif

          pllTreeToNewick(buffer, tr, pr, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

          if(thoroughIterations % 2 == 0)	      
            memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
          else
            memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    

          rax_free(buffer);
        }

        assert(bCounter == tr->mxtips - 3);

        if(thoroughIterations > 0)
        {
          double rrf = convergenceCriterion(tr->h, tr->mxtips);

          if(rrf <= 0.01) /* 1% cutoff */
          {
            printBothOpen("ML search converged at thorough SPR cycle %d with stopping criterion\n", thoroughIterations);
            printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
            goto cleanup;
          }
          else		    
            printBothOpen("ML search convergence criterion thorough cycle %d->%d Relative Robinson-Foulds %f\n", thoroughIterations - 1, thoroughIterations, rrf);
        }
      }



      thoroughIterations++;	  
    }
    else
    {

      /* if the lnl has not imrpved by the current SPR cycle adapt the min and max rearrangemnt radii and try again */

      rearrangementsMax += adef->stepwidth;
      rearrangementsMin += adef->stepwidth; 	        	      

      /* if we have already tried them then abandon this loop, the search has converged */
      if(rearrangementsMax > adef->max_rearrange)	     	     	 
        goto cleanup; 	   
    }

    /* optimize branch lengths of best tree */

    treeEvaluate(tr, pr, 32 ); // 32 * 1

    /* do some bokkeeping and printouts again */
    previousLh = lh = tr->likelihood;	      
    saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
    printLog(tr);

    /* do a cycle of thorough SPR moves with the minimum and maximum rearrangement radii */

    treeOptimizeRapid(tr, pr, rearrangementsMin, rearrangementsMax, bt, iList);

    impr = 0;			      		            

    /* once again get the best 20 trees produced by the SPR cycle, load them from the bt tree list into tr
       optimize their branch lengths and figure out if the LnL of the tree has improved */

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

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("SLOW LH: %f\n", lh);              
#endif
  }

cleanup: 

  /* do a final full tree traversal, not sure if this is required here */

  {
    evaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("After SLOW SPRs Final %f\n", tr->likelihood);   
#endif
  }

  /* free data structures */

  if(tr->searchConvergenceCriterion)
  {
    freeBitVectors(tr->bitVectors, 2 * tr->mxtips);
    rax_free(tr->bitVectors);
    freeHashTable(tr->h);
    rax_free(tr->h);
  }

  freeBestTree(bestT);
  rax_free(bestT);
  freeBestTree(bt);
  rax_free(bt);
  freeInfoList(iList);  
  rax_free(iList);

  printLog(tr);

  printResult(tr, pr, adef, PLL_TRUE);

  /* and we are done, return to main() in axml.c  */

}
#endif

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

    if(pr->partitionData[model]->protModels == PLL_LG4)
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

