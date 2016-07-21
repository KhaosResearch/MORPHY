#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pll/pll.h>
#include <time.h>

static void printModelParameters(partitionList *pr)
{
  int
    i;

  for(i = 0; i < pr->numberOfPartitions; i++)
    {
      int
    j,
    states = pr->partitionData[i]->states,
    rates = (states * states - states) / 2;

      printf("Partition %d\n", i);

      printf("alpha: %f\n", pr->partitionData[i]->alpha);

      printf("Frequencies: \n");

      for(j = 0; j < states; j++)
    printf("%f ",  pr->partitionData[i]->frequencies[j]);
      printf("\n");

      printf("SubstRates: \n");

      for(j = 0; j < rates; j++)
    printf("%f ",  pr->partitionData[i]->substRates[j]);

     printf("\n");


     printf ("Gamma rates.....: %f %f %f %f\n\n\n",
          pr->partitionData[i]->gammaRates[0],
          pr->partitionData[i]->gammaRates[1],
          pr->partitionData[i]->gammaRates[2],
          pr->partitionData[i]->gammaRates[3]);
    }

}



int main (int argc, char * argv[])
{
  pllAlignmentData * alignmentData;
  pllInstance * tr;
  pllNewickTree * newick;
  partitionList * partitions;
  pllQueue * partitionInfo;
  int i;
  pllInstanceAttr attr;
  pllRearrangeList * rearrangeList;

#ifdef _FINE_GRAIN_MPI
  pllInitMPI (&argc, &argv);
#endif

  if (argc < 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file] [threads]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  clock_t  t_fin, t_ini= clock();

  /* Set the PLL instance attributes */
  attr.rateHetModel     = PLL_GAMMA;
  attr.fastScaling      = PLL_FALSE;
  attr.saveMemory       = PLL_FALSE;
  attr.useRecom         = PLL_FALSE;
  attr.randomNumberSeed = 0xDEADBEEF;
  attr.numberOfThreads  = (argc > 4) ? (atoi(argv[4]) > 0 ? atoi(argv[4]) : 8) : 8;            /* This only affects the pthreads version */

  /* Create a PLL tree */
  tr = pllCreateInstance (&attr);

  /* Parse a PHYLIP file */
  alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, argv[1]);


  if (!alignmentData)
   {
     fprintf (stderr, "Error while parsing %s\n", argv[1]);
     return (EXIT_FAILURE);
   }

  /* Parse a NEWICK file */
  newick = pllNewickParseFile (argv[2]);

  if (!newick)
   {
     fprintf (stderr, "Error while parsing newick file %s\n", argv[2]);
     return (EXIT_FAILURE);
   }
  if (!pllValidateNewick (newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
   {
     fprintf (stderr, "Invalid phylogenetic tree\n");
     printf ("%d\n", errno);
     return (EXIT_FAILURE);
   }

  /* Parse the partitions file into a partition queue structure */
  partitionInfo = pllPartitionParse (argv[3]);
  
  /* Validate the partitions */
  if (!pllPartitionsValidate (partitionInfo, alignmentData))
   {
     fprintf (stderr, "Error: Partitions do not cover all sites\n");
     return (EXIT_FAILURE);
   }

  /* Commit the partitions and build a partitions structure */
  partitions = pllPartitionsCommit (partitionInfo, alignmentData);

  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&partitionInfo);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllAlignmentRemoveDups (alignmentData, partitions);

  /* Set the topology of the PLL tree from a parsed newick tree */
  //PLL_FALSE: the branch lengths of the parsed tree are reset to some default values.
  pllTreeInitTopologyNewick (tr, newick, PLL_FALSE);

   /*Or instead of the previous function use the next commented line to create
     a random tree topology */
  //pllTreeInitTopologyRandom (tr, alignmentData->sequenceCount, alignmentData->sequenceLabels); 

  /* Connect the alignment and partition structure with the tree structure */
  if (!pllLoadAlignment (tr, alignmentData, partitions))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  
  
  /*********************** Parsimony Evaluation **************************/

  /* Before testing
   printf ("Max tips: %d bestParsimony %d\n", tr->mxtips , tr->bestParsimony);
   allocateParsimonyDataStructures(tr, partitions);
   int p= evaluateParsimony(tr, partitions, tr->start, PLL_TRUE);
   pllMakeParsimonyTreeFast(tr,partitions);
   printf(" Parsiomona %d \n" ,p);
  */

  /*int j;
   pllInitParsimonyStructures (tr, partitions, PLL_TRUE);   // the last parameter is that you want to allocate buffers for the per site computation 
   printf ("Score Parsimony : %d\n", pllEvaluateParsimony(tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE)); 
   for (i = 0; i < partitions->numberOfPartitions; ++i) 
   { 
    printf ("Per-site scores for partition %d\n", i); 
    for (j = 0; j < partitions->partitionData[i]->width; ++j) printf ("%d ", partitions->partitionData[i]->perSiteParsScores[j]); 
    printf ("\n"); 
  } */
  
  pllInitParsimonyStructures (tr, partitions, PLL_TRUE);
  tr->bestParsimony = pllEvaluateParsimony(tr, partitions, tr->start, PLL_TRUE, PLL_TRUE);

  printf ("Valor Parsimony: %d\n", tr->bestParsimony);

/* Initialize the model. Note that this function will also perform a full
     tree traversal and evaluate the likelihood of the tree. Therefore, you
     have the guarantee that tr->likelihood the valid likelihood */

    pllInitModel(tr, partitions);
    printf ("Log-likelihood INITIAL of topology: %f\n", tr->likelihood);

    //55
	double alpha=0.365;
	double q[6] = {1.7283,5.1263,0.6618,1.9542,8.3688,1}; 
	double f[4] = {0.3104,0.1588 ,0.1669 ,0.3639 }; 
   
  
	//186
	/*double alpha=0.05;
  	double q[6] = {1.5849,58.5552,1.1812,1.6618,39.6149,1};
        double f[4] = {0.3104,0.3166, 0.1296, 0.2434}; */

	//218
	/* double alpha=0.350863;
  	double q[6] = {1,1,1,1,1,1};
	double f[4] = {0.24484,0.23583,0.31868, 0.20065}; 
	*/

       //500
	/*double alpha=0.867;
	double q[6] = {1.24303,4.02694,0.34805,1.73304,3.25578,1}; //500
        double f[4] = {0.28759,0.17822,0.14682,0.38737}; 
  */
  
  //pllSetFixedAlpha(alpha , 0, partitions, tr);    
  //pllSetFixedSubstitutionMatrix(q, 6, 0, partitions, tr);
  //pllSetFixedBaseFrequencies(f, 4, 0, partitions, tr);

  pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("Valor de Log-likelihood AFTER MODELPARAMETSR AND FREQUENCIES: %f\n", tr->likelihood);

  //pllOptimizeBranchLengths (tr, partitions, 100000);
  //printf ("Valor de Log-likelihood After Branch Length  Optimization : %f\n", tr->likelihood);



   //printf ("Valor de bestParsimony: %d\n", tr->bestParsimony);

  //printModelParameters(partitions);
  //pllOptimizeModelParameters(tr, partitions, 0.000001);
  //printf ("Valor de Log-likelihood After pllOptimizeModelParameters: %f\n", tr->likelihood);
  //printModelParameters(partitions);


  /* Now let's create another list and compute 30 rearrangement moves */
  rearrangeList = pllCreateRearrangeList (30);

  /* The next flag specifies that the length of the new branch that is created
     by an SPR move need not be optimized */
  tr->thoroughInsertion = PLL_TRUE;


  printf ("Computing the best 30 SPR in radius (1,30)\n");
  pllRearrangeSearch (tr, partitions, 
                      PLL_REARRANGE_SPR, 
                      tr->nodep[tr->mxtips + 1], 
                      1, 
                      30, 
                      rearrangeList);

  pllEvaluateLikelihood(tr, partitions, tr->start, PLL_FALSE, PLL_FALSE);
  printf ("Current log-likelihood: %f\n\n", tr->likelihood);

  printf ("Number of computed rearrangements: %d\n", rearrangeList->entries);
  printf ("------------------------------------\n");
  for (i = 0; i < rearrangeList->entries; ++ i)
   {
     printf ("%2d  Type: SPR  Likelihood: %f \n", i, rearrangeList->rearr[i].likelihood);
   }

  printf ("Current log-likelihood BEFORE COMMIT SPRMOVE: %f\n\n", tr->likelihood);


  printf ("Committing rearrangeList->rearr[0]\n");
  pllRearrangeCommit (tr, partitions, &(rearrangeList->rearr[0]), PLL_TRUE);
  //pllOptimizeBranchLengths (tr, partitions, 64);
  pllEvaluateLikelihood (tr, partitions, tr->start, PLL_FALSE, PLL_FALSE);
  int p = pllEvaluateParsimony(tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE);
  printf ("New log-likelihood: %f Parsimony %d \n\n", tr->likelihood,p);

  /* Rolling back to the previous topology. Note that if we evaluate the
     likelihood with a partial traversal we might get an invalid log likelihood.
     This is due to the fact that the likelihood vectors no longer correspond
     to the old topology, hence we need to do full traversal. I left the
     partial traversal here as an example */
  
  printf ("Rolling back...\n");
  pllRearrangeRollback (tr, partitions);
  pllEvaluateLikelihood (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New log-likelihood AFTER ROLLBAKC TRUU FULL TRANSVERSAL: %f\n\n", tr->likelihood);

  /* We do one more rollback to get to the original topology, but this time we
     do a full traversal to fix the log-likelihood to the correct value plus we
     do branch-length optimization */
  printf ("Rolling back...\n");
  pllRearrangeRollback (tr, partitions);
  pllEvaluateLikelihood (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  //pllOptimizeBranchLengths (tr, partitions, 64);
  printf ("New log-likelihood: %f\n\n", tr->likelihood);

  /* DEallocate the rearrange list */
  pllDestroyRearrangeList (&rearrangeList);

  t_fin = clock();
  double secs = ((double) (t_fin - t_ini))/ CLOCKS_PER_SEC;
  printf("Total execution time: %f \n",secs);

  /* Do some cleanup */
  pllAlignmentDataDestroy (alignmentData);
  pllNewickParseDestroy (&newick);
  pllPartitionsDestroy (tr, &partitions);
  pllDestroyInstance (tr);
  
  return (EXIT_SUCCESS);
}
