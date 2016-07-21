#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pll.h>

int main (int argc, char * argv[])
{
  pllAlignmentData * alignmentData;
  pllInstance * tree;
  pllNewickTree * newick;
  partitionList * partitions;
  pllQueue * partitionInfo;
  pllInstanceAttr attr;
  pllRearrangeList * rearrangeList;

  int i,
      useDefaultBranchLengths = PLL_TRUE,
      fullTraversal = PLL_TRUE,
      perSiteLikelihoods = PLL_FALSE,
      saveRollbackInfo = PLL_TRUE;

#ifdef _FINE_GRAIN_MPI
  pllInitMPI (&argc, &argv);
#endif

  if (argc < 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file] [threads]\n", argv[0]);
     return (EXIT_FAILURE);
   }

#ifdef _USE_PTHREADS
   /* The number of threads must be greater than 0 */
   if ((argc > 4) && (atoi(argv[4]) <= 0)) {
     fprintf (stderr, "Error: The number of threads must be greater than 0\n");
     return (EXIT_FAILURE);
   }
#endif

  /* Set the PLL instance attributes */
  attr.rateHetModel     = PLL_GAMMA;
  attr.fastScaling      = PLL_FALSE;
  attr.saveMemory       = PLL_FALSE;
  attr.useRecom         = PLL_FALSE;
  attr.randomNumberSeed = 0xDEADBEEF;
  attr.numberOfThreads  = (argc > 4) ? atoi(argv[4]) : 2;  /* This only affects the pthreads version */

  /* Create a PLL tree */
  tree = pllCreateInstance (&attr);

  if (!tree) {
    fprintf (stderr, "Error creating PLL instance\n");
    return (EXIT_FAILURE);
  }


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

  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
  if (!pllValidateNewick (newick))  
   {
     fprintf (stderr, "Invalid phylogenetic tree\n");
     printf ("%d\n", errno);
     return (EXIT_FAILURE);
   }

  /* Parse the partitions file into a partition queue structure */
  partitionInfo = pllPartitionParse (argv[3]);
  
  if (!pllPartitionsValidate (partitionInfo, alignmentData))
   {
     fprintf (stderr, "Error: Partitions do not cover all sites\n");
     return (EXIT_FAILURE);
   }

  /* Commit the partitions and build a partitions structure */
  partitions = pllPartitionsCommit (partitionInfo, alignmentData);
  partitions->perGeneBranchLengths = PLL_TRUE;

  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&partitionInfo);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllAlignmentRemoveDups (alignmentData, partitions);

  /* Set the topology of the PLL tree from a parsed newick tree */
  pllTreeInitTopologyNewick (tree, newick, useDefaultBranchLengths);

  /* Or instead of the previous function use the next commented line to create
     a random tree topology 
  pllTreeInitTopologyRandom (tree, alignmentData->sequenceCount, alignmentData->sequenceLabels); */

  /* Connect the alignment and partition structure with the tree structure */
  if (!pllLoadAlignment (tree, alignmentData, partitions))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  
  /* Initialize the model and starts PTHREADS (if it applies). Note that this function 
     will also perform a full tree traversal and evaluate the likelihood of the tree. 
     Therefore, you have the guarantee that tr->likelihood contains the valid likelihood */
  pllInitModel(tree, partitions);

  pllOptimizeBranchLengths (tree, partitions, 64);

  printf ("Log-likelihood of topology: %f\n", tree->likelihood);

  /* Create a list that will hold information for at most 20 rearrangement moves */
  rearrangeList = pllCreateRearrangeList (20);


  /* The next flag specifies that PLL optimizes the length of the new branch
     that is created by an SPR move */
  tree->thoroughInsertion = PLL_FALSE;

  /* Note that the following commands will fill the list with at most 20 
     SPR and NNI rearrangement moves, i.e. the best 20 will appear in the
     list */
  printf ("Computing the best 20 SPR and NNI rearrangements in radius (1,20)\n");
  pllRearrangeSearch (tree, 
                      partitions, 
                      PLL_REARRANGE_SPR, 
                      tree->nodep[tree->mxtips + 1], 
                      1, 
                      20,
                      rearrangeList);

  pllRearrangeSearch (tree, 
                      partitions, 
                      PLL_REARRANGE_NNI, 
                      tree->nodep[tree->mxtips + 1], 
                      1, 
                      20, 
                      rearrangeList);

  printf ("Number of computed rearrangements: %d\n", rearrangeList->entries);
  printf ("------------------------------------\n");

  for (i = 0; i < rearrangeList->entries; ++ i)
   {
     printf ("%2d  Type: %s  Log-likelihood: %f\n", 
             i, 
             rearrangeList->rearr[i].rearrangeType == PLL_REARRANGE_SPR ? 
                "SPR" : "NNI", 
             rearrangeList->rearr[i].likelihood);
   }


  printf ("Committing move 0\n");
  pllRearrangeCommit(tree, partitions, &(rearrangeList->rearr[0]), saveRollbackInfo);

  fullTraversal = PLL_TRUE;
  pllEvaluateLikelihood (tree, partitions, tree->start, fullTraversal, perSiteLikelihoods);
  printf ("New log-likelihood: %f\n\n", tree->likelihood);

  /* We don't need the rearrange list anymore */
  pllDestroyRearrangeList (&rearrangeList);

  /* Now let's create another list and compute 30 rearrangement moves */
  rearrangeList = pllCreateRearrangeList (30);

  /* The next flag specifies that the length of the new branch that is created
     by an SPR move need not be optimized */
  tree->thoroughInsertion = PLL_TRUE;

  printf ("Computing the best 30 SPR in radius (1,30)\n");
  pllRearrangeSearch (tree, partitions, 
                      PLL_REARRANGE_SPR, 
                      tree->nodep[tree->mxtips + 1], 
                      1, 
                      30, 
                      rearrangeList);

  printf ("Number of computed rearrangements: %d\n", rearrangeList->entries);
  printf ("------------------------------------\n");
  for (i = 0; i < rearrangeList->entries; ++ i)
   {
     printf ("%2d  Type: SPR  Likelihood: %f\n", i, rearrangeList->rearr[i].likelihood);
   }

  printf ("Committing rearrangeList->rearr[0]\n");
  pllRearrangeCommit (tree, partitions, &(rearrangeList->rearr[0]), saveRollbackInfo);

  fullTraversal = PLL_FALSE;
  pllEvaluateLikelihood (tree, partitions, tree->start, fullTraversal, perSiteLikelihoods);
  printf ("New log-likelihood: %f\n\n", tree->likelihood);

  /* Rolling back to the previous topology. Note that if we evaluate the
     likelihood with a partial traversal we might get an invalid log likelihood.
     This is due to the fact that the likelihood vectors no longer correspond
     to the old topology, hence we need to do full traversal. I left the
     partial traversal here as an example */
  printf ("Rolling back...\n");
  pllRearrangeRollback (tree, partitions);
  pllEvaluateLikelihood (tree, partitions, tree->start, fullTraversal, perSiteLikelihoods);
  printf ("New log-likelihood: %f\n\n", tree->likelihood);

  /* We do one more rollback to get to the original topology, but this time we
     do a full traversal to fix the log-likelihood to the correct value plus we
     do branch-length optimization */
  printf ("Rolling back...\n");
  pllRearrangeRollback (tree, partitions);

  fullTraversal = PLL_TRUE;
  pllEvaluateLikelihood (tree, partitions, tree->start, fullTraversal, perSiteLikelihoods);
  pllOptimizeBranchLengths (tree, partitions, 64);
  printf ("New log-likelihood: %f\n\n", tree->likelihood);

  /* DEallocate the rearrange list */
  pllDestroyRearrangeList (&rearrangeList);

  /* Do some cleanup */
  pllAlignmentDataDestroy (alignmentData);
  pllNewickParseDestroy (&newick);
  /* This function also stops PTHREADS if it applies */
  pllPartitionsDestroy (tree, &partitions);
  pllDestroyInstance (tree);
  
  return (EXIT_SUCCESS);
}
