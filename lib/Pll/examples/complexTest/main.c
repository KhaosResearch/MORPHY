#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../../axml.h"
#include "../../utils.h"
#include "../../lexer.h"
#include "../../hash.h"

static nodeptr pickMyRandomSubtree(pllInstance *tr)
{
  nodeptr p;
  //do
  {
    /* select a random inner node */
    p = tr->nodep[(rand() % (tr->mxtips - 2)) + 1 + tr->mxtips];

    /* select a random orientation */
    int exitDirection = rand() % 3;
    switch(exitDirection)
    {
      case 0:
        break;
      case 1:
        p = p->next;
        break;
      case 2:
        p = p->next->next;
        break;
      default:
        assert(0);
    }
  }
  //while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));
  return p;
}


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
      printf("\n\n");
      
    }

}

static void testProteinStuff()
{
  pllAlignmentData * alignmentData;
  pllInstance * tr;
  pllNewickTree * newick;
  
  partitionList * partitions;
  
  struct pllQueue * parts;
  
  int i;
  
  for(i = 0; i < 5; i++)
    {
      //write a simple partition file with 3 partitions 
      //for dataset dna.phy.dat contained 
      //in this source directory 
      
      FILE *f = fopen("proteinPartitions", "w");
      
      switch(i)
	{
	case 0:
	  fprintf(f, "WAG, p1 = 1-200\n");
	  fprintf(f, "WAG, p2 = 201-600\n");
	  fprintf(f, "WAG, p3 = 601-1104\n");
	  break;
	case 1:
	  fprintf(f, "LG, p1 = 1-200\n");
	  fprintf(f, "LG, p2 = 201-600\n");
	  fprintf(f, "LG, p3 = 601-1104\n");
	  break;
	case 2:
	  fprintf(f, "JTT, p1 = 1-200\n");
	  fprintf(f, "JTT, p2 = 201-600\n");
	  fprintf(f, "JTT, p3 = 601-1104\n");
	  break;
	case 3:	  
	case 4:
	  fprintf(f, "GTR, p1 = 1-200\n");
	  fprintf(f, "GTR, p2 = 201-600\n");
	  fprintf(f, "GTR, p3 = 601-1104\n");
	  break;
	default:
	  assert(0);
	}
      
      fclose(f);
      
      tr = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
      
      alignmentData = pllParsePHYLIP ("prot.phy");

      /* or alternatively, parse a FASTA file */
      // alignmentData = pllParseFASTA ("prot.phy");
      
      newick = pllNewickParseFile("parsimonyTree");
      
      parts = pllPartitionParse ("proteinPartitions");
      
      /* Validate the partitions */
      if (!pllPartitionsValidate (parts, alignmentData))
	{
	  fprintf (stderr, "Error: Partitions do not cover all sites\n");
	  return (EXIT_FAILURE);
	}
      
      /* commit the partitions and build a partitions structure */
      partitions = pllPartitionsCommit (parts, alignmentData);
      
      /* destroy the  intermedia partition queue structure */
      pllQueuePartitionsDestroy (&parts);
      
      /* eliminate duplicate sites from the alignment and update weights vector */
      pllPhylipRemoveDuplicate (alignmentData, partitions);
      
      pllTreeInitTopologyNewick (tr, newick, PLL_TRUE);
      if (!pllLoadAlignment (tr, alignmentData, partitions, PLL_DEEP_COPY))
	{
	  fprintf (stderr, "Incompatible tree/alignment combination\n");
	  return (EXIT_FAILURE);
	}
      //pllInitModel(tr, PLL_TRUE, alignmentData, partitions);
      pllInitModel(tr, alignmentData, partitions);
      
      switch(i)
	{
	case 0:
	  //all params unlinked 
	  
	  pllLinkAlphaParameters("0,1,2", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);	
	  break;
	case 1:
	  //link params in another way 
	  
	  pllLinkAlphaParameters("0,0,0", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);    
	  break;
	case 2:
	  //link params in yet another way 
	  
	  pllLinkAlphaParameters("0,0,0", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,0", partitions);    	
	  break;

	case 3:
	  //also fiddle around with the Q matrices, make them to be non-GTR, but simpler
	  
	  pllLinkAlphaParameters("0,1,2", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);    
	  
	  //these are GTR models
	  //pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,5", partitions, 0);	  
	  //pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,5", partitions, 1);

	  //this is a simpler model with 5 parameters, parameter a and f have 
	  //the same value
	  //pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,0", partitions, 2);
	  break;

	case 4:
	  {
	    //test case to show how the model parameters can be set to fixed values

	    // set up arrays of user-defined base frequencies 
	    // and a user defined q matrix 
	    double 
	      f[4] = {0.25, 0.25, 0.25, 0.25},
	      q[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.5};
	    
	      //unlink alpha parameters base frequencies and Q matrices 
	      //across all partitions
	    pllLinkAlphaParameters("0,1,2", partitions);
	    pllLinkFrequencies("0,1,2", partitions);
	    pllLinkRates("0,0,0", partitions);
	    
	    //set alpha to a fixed value of 1.0 for partition 0 and 
	    //parition 1
	    //pllSetFixedAlpha(1.0, 0, partitions, tr);
	    //pllSetFixedAlpha(1.0, 1, partitions, tr);

	    //fix the base frequencies to 0.25 for 
	    //partitions 0 and 1
	    //pllSetFixedBaseFrequencies(f, 4, 0, partitions, tr);
	    //pllSetFixedBaseFrequencies(f, 4, 1, partitions, tr);
	    
	    //set the Q matrix to fixed values for partition 
	    //0
	    //pllSetFixedSubstitutionMatrix(q, 6, 0, partitions, tr);	    	    
	  }	  
	  break;
	default:
	  assert(0);
	}
      
      evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
      printf("%f \n", tr->likelihood);
      pllOptimizeModelParameters(tr, partitions, 1.0);

      //print the model parameters 

      //printModelParameters(partitions);

      printf("%f \n", tr->likelihood); 
      //cleanup
      pllAlignmentDataDestroy (alignmentData);
      pllNewickParseDestroy (&newick);
      
      pllPartitionsDestroy (tr, &partitions);
      pllTreeDestroy (tr);      
    }
}

int main (int argc, char * argv[])
{
  pllAlignmentData *alignmentData1, *alignmentData2;
  pllInstance * tr, *tr2;
  pllNewickTree * newick;
  partitionList * partitions, *partitions2;
  struct pllQueue * parts;
  int i;

  if (argc != 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* Create a PLL tree */
  tr = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
  tr2 = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);

  /* Parse a PHYLIP file */
  alignmentData1= pllParsePHYLIP (argv[1]);
  alignmentData2 = pllParsePHYLIP (argv[1]);

  if (!alignmentData1)
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
     return (EXIT_FAILURE);
   }

  /* Parse the partitions file into a partition queue structure */
  parts = pllPartitionParse (argv[3]);

  /* Validate the partitions */
  if (!pllPartitionsValidate (parts, alignmentData1))
   {
     fprintf (stderr, "Error: Partitions do not cover all sites\n");
     return (EXIT_FAILURE);
   }

  /* commit the partitions and build a partitions structure */
  partitions = pllPartitionsCommit (parts, alignmentData1);
  partitions2 =  pllPartitionsCommit (parts, alignmentData2);

  /* destroy the  intermedia partition queue structure */
  pllQueuePartitionsDestroy (&parts);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllPhylipRemoveDuplicate (alignmentData1, partitions);
  pllPhylipRemoveDuplicate (alignmentData2, partitions2);


  /* Set the topology of the PLL tree from a parsed newick tree */
  //pllTreeInitTopologyNewick (tr, newick, PLL_TRUE);
  /* Or instead of the previous function use the next commented line to create
     a random tree topology
  pllTreeInitTopologyRandom (tr, phylip->nTaxa, phylip->label); */

  pllTreeInitTopologyForAlignment(tr, alignmentData1);

  /* Connect the alignment with the tree structure */
  if (!pllLoadAlignment (tr, alignmentData1, partitions, PLL_DEEP_COPY))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }

  /* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
 pllInitModel(tr, alignmentData1, partitions);

  /* TODO transform into pll functions !*/

  /*
     allocateParsimonyDataStructures(tr, partitions);
     pllMakeParsimonyTreeFast(tr, partitions);
     pllFreeParsimonyDataStructures(tr, partitions);
  */

  pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
  pllTreeToNewick (tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
  printf ("Tree: %s %d\n", tr->tree_string, tr->start->number);
  evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

  double
    firstTree = tr->likelihood;

  printf("%f \n", tr->likelihood);
  //computeBIGRAPID_Test(tr, partitions, PLL_TRUE);
  printf("final like %f\n", tr->likelihood);
  //pllInitModel(tr, PLL_TRUE, phylip, partitions);

  pllTreeInitTopologyNewick (tr2, newick, PLL_TRUE);
  if (!pllLoadAlignment (tr2, alignmentData2, partitions2, PLL_DEEP_COPY))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  pllInitModel(tr2, alignmentData2, partitions2);

  pllTreeToNewick (tr2->tree_string, tr2, partitions2, tr2->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
  printf ("Tree: %s %d\n", tr2->tree_string, tr2->start->number);
  evaluateGeneric(tr2, partitions2, tr2->start, PLL_TRUE, PLL_FALSE);

  printf("%f \n", tr2->likelihood);

  double
    secondTree = tr2->likelihood;

  assert(firstTree == secondTree);

  pllOptimizeModelParameters(tr2, partitions2, 10.0);

  printf("%f \n", tr2->likelihood);

  pllAlignmentDataDestroy (alignmentData1);
  pllNewickParseDestroy (&newick);

  pllPartitionsDestroy (tr, &partitions);
  pllTreeDestroy (tr);

  pllAlignmentDataDestroy (alignmentData2); 
  pllPartitionsDestroy (&partitions2, tr2->mxtips);
  pllTreeDestroy (tr2);



  for(i = 0; i < 5; i++)
    {
      //write a simple partition file with 3 partitions 
      //for dataset dna.phy.dat contained 
      //in this source directory 
      
      FILE *f = fopen("dummy", "w");
      
      fprintf(f, "DNA, p1 = 1-200\n");
      fprintf(f, "DNA, p1 = 201-400\n");
      fprintf(f, "DNA, p1 = 401-705\n");
      
      fclose(f);
      
      tr = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
      
      alignmentData1= pllParsePHYLIP (argv[1]);
      
      newick = pllNewickParseFile (argv[2]);
      
      parts = pllPartitionParse ("dummy");
      
      /* Validate the partitions */
      if (!pllPartitionsValidate (parts, alignmentData1))
	{
	  fprintf (stderr, "Error: Partitions do not cover all sites\n");
	  return (EXIT_FAILURE);
	}
      
      /* commit the partitions and build a partitions structure */
      partitions = pllPartitionsCommit (parts, alignmentData1);
      
      /* destroy the  intermedia partition queue structure */
      pllQueuePartitionsDestroy (&parts);
      
      /* eliminate duplicate sites from the alignment and update weights vector */
      pllPhylipRemoveDuplicate (alignmentData1, partitions);
      
      pllTreeInitTopologyNewick (tr, newick, PLL_TRUE);
      if (!pllLoadAlignment (tr, alignmentData1, partitions, PLL_DEEP_COPY))
	{
	  fprintf (stderr, "Incompatible tree/alignment combination\n");
	  return (EXIT_FAILURE);
	}
      pllInitModel(tr, alignmentData1, partitions);
      
      switch(i)
	{
	case 0:
	  //link params in one way 
	  
	  pllLinkAlphaParameters("0,1,2", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);	
	  break;
	case 1:
	  //link params in another way 
	  
	  pllLinkAlphaParameters("0,0,0", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);    
	  break;
	case 2:
	  //link params in yet another way 
	  
	  pllLinkAlphaParameters("0,0,0", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,0", partitions);    	
	  break;

	case 3:
	  //also fiddle around with the Q matrices, make them to be non-GTR, but simpler
	  
	  pllLinkAlphaParameters("0,1,2", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);    
	  
	  //these are GTR models
	  pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,5", partitions, 0);	  
	  pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,5", partitions, 1);

	  //this is a simpler model with 5 parameters, parameter a and f have 
	  //the same value
	  pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,0", partitions, 2);
	  break;

	case 4:
	  {
	    //test case to show how the model parameters can be set to fixed values

	    // set up arrays of user-defined base frequencies 
	    // and a user defined q matrix 
	    double 
	      f[4] = {0.25, 0.25, 0.25, 0.25},
	      q[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.5};
	    
	      //unlink alpha parameters base frequencies and Q matrices 
	      //across all partitions
	    pllLinkAlphaParameters("0,1,2", partitions);
	    pllLinkFrequencies("0,0,1", partitions);
	    pllLinkRates("0,1,2", partitions);
	    
	    //set alpha to a fixed value of 1.0 for partition 0 and 
	    //parition 1
	    pllSetFixedAlpha(1.0, 0, partitions, tr);
	    pllSetFixedAlpha(1.0, 1, partitions, tr);

	    //fix the base frequencies to 0.25 for 
	    //partitions 0 and 1
	    pllSetFixedBaseFrequencies(f, 4, 0, partitions, tr);
	    pllSetFixedBaseFrequencies(f, 4, 1, partitions, tr);
	    
	    //set the Q matrix to fixed values for partition 
	    //0
	    pllSetFixedSubstitutionMatrix(q, 6, 0, partitions, tr);	    	    
	  }	
	  break;
	default:
	  assert(0);
	}
      
      evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
      printf("%f \n", tr->likelihood);
      pllOptimizeModelParameters(tr, partitions, 10.0);

      //print the model parameters 

      printModelParameters(partitions);

      printf("%f \n", tr->likelihood); 
      //cleanup
      pllAlignmentDataDestroy (alignmentData1);
      pllNewickParseDestroy (&newick);
      
      pllPartitionsDestroy (&partitions, tr->mxtips);
      pllTreeDestroy (tr);      
    }
  
  testProteinStuff();

  return (EXIT_SUCCESS);
}
