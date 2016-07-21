
/******************************************** axml.h ******************************************************/
#define LIKELIHOOD_EPSILON                      0.0000001

#define INVAR_MIN                               0.0001
#define INVAR_MAX                               0.9999

#define TT_MIN                                  0.0000001
#define TT_MAX                                  1000000.0

#define FREQ_MIN                                0.001
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))
#define LOGF(x) logf(x)
#define PLL_GAMMA_I                             2
#define PLL_MEM_APROX_OVERHEAD     1.3 /* TODOFER can we measure this empirically? */

#define  TREE_EVALUATION            0
#define  BIG_RAPID_MODE             1
#define  CALC_BIPARTITIONS          3
#define  SPLIT_MULTI_GENE           4
#define  CHECK_ALIGNMENT            5
#define  PER_SITE_LL                6
#define  PARSIMONY_ADDITION         7
#define  CLASSIFY_ML                9
#define  DISTANCE_MODE              11
#define  GENERATE_BS                12
#define  COMPUTE_ELW                13
#define  BOOTSTOP_ONLY              14
#define  COMPUTE_LHS                17
#define  COMPUTE_BIPARTITION_CORRELATION 18
#define  THOROUGH_PARSIMONY         19
#define  COMPUTE_RF_DISTANCE        20
#define  MORPH_CALIBRATOR           21
#define  CONSENSUS_ONLY             22
#define  MESH_TREE_SEARCH           23
#define  FAST_SEARCH                24
#define  MORPH_CALIBRATOR_PARSIMONY 25
#define  SH_LIKE_SUPPORTS           28

#define  GPU_BENCHMARK              29

#define M_GTRCAT         1
#define M_GTRGAMMA       2
#define M_BINCAT         3
#define M_BINGAMMA       4
#define M_PROTCAT        5
#define M_PROTGAMMA      6
#define M_32CAT          7
#define M_32GAMMA        8
#define M_64CAT          9
#define M_64GAMMA        10


#define BIPARTITIONS_ALL       0
#define GET_BIPARTITIONS_BEST  1
#define DRAW_BIPARTITIONS_BEST 2
#define BIPARTITIONS_BOOTSTOP  3

/* bootstopping stuff */

//#define BOOTSTOP_PERMUTATIONS 100
//#define START_BSTOP_TEST      10

//#define FC_THRESHOLD          99
//#define FC_SPACING            50
//#define FC_LOWER              0.99
//#define FC_INIT               20

//#define FREQUENCY_STOP 0
//#define MR_STOP        1
//#define MRE_STOP       2
//#define MRE_IGN_STOP   3

//#define MR_CONSENSUS 0
//#define MRE_CONSENSUS 1
//#define STRICT_CONSENSUS 2



/* bootstopping stuff end */

/** @brief ???Expected likelihood weight
 * @todo add explanation, is this ever used?  */
typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

/*
  typedef uint64_t parsimonyNumber;

  #define PCF 16


typedef unsigned char parsimonyNumber;

#define PCF 2
*/

/** @brief LH entries @warning UNUSED???
  */
typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;
/** @brief LH list. @warning UNUSED???
  */
typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;


#define randomTree    0
#define givenTree     1 
#define parsimonyTree 2

extern int determineRearrangementSetting ( pllInstance *tr, partitionList *pr, analdef *adef, bestlist *bestT, bestlist *bt );
extern void computeBIGRAPID ( pllInstance *tr, partitionList *pr, analdef *adef, boolean estimateModel);
extern void computeBIGRAPID_Test (pllInstance *tr, partitionList *pr, boolean estimateModel);
extern boolean treeEvaluatePartition ( pllInstance *tr, double smoothFactor, int model );
extern void treeEvaluateRandom (pllInstance *tr, double smoothFactor);
extern void treeEvaluateProgressive(pllInstance *tr);
extern boolean compatible(entry* e1, entry* e2, unsigned int bvlen);
extern boolean issubset(unsigned int* bipA, unsigned int* bipB, unsigned int vectorLen);
extern void testGapped(pllInstance *tr);
extern void readBinaryModel(pllInstance *tr);
extern void writeBinaryModel(pllInstance *tr);
extern int *permutationSH(pllInstance *tr, int nBootstrap, long _randomSeed);
void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile);
boolean modelExists(char *model, pllInstance *tr);
void reorder( double *x, int n, int span );
void reorder_back( double *x, int n, int span );
void init_default(pllInstance *tr);
nodeptr pickRandomSubtree (pllInstance *tr);
static boolean setupTree (pllInstance *tr, boolean doInit, partitionList *partitions);
void allocNodex(pllInstance *tr, int tid, int n);

void threadMakeVector(pllInstance *tr, int tid);
void threadComputeAverage(pllInstance *tr, int tid);
void threadComputePearson(pllInstance *tr, int tid);
extern FILE *getNumberOfTrees(pllInstance *tr, char *fileName, analdef *adef);
extern hashtable *copyHashTable(hashtable *src, unsigned int vectorLength);
extern void compareBips(pllInstance *tr, char *bootStrapFileName, analdef *adef);
extern void computeRF(pllInstance *tr, char *bootStrapFileName, analdef *adef);
extern void printPartitions(pllInstance *tr);
extern void parseSecondaryStructure(pllInstance *tr, analdef *adef, int sites);
extern void reductionCleanup(pllInstance *tr, int *originalRateCategories, int *originalInvariant);
extern void computeNextReplicate(pllInstance *tr, long *seed, int *originalRateCategories, int *originalInvariant, boolean isRapid, boolean fixRates);
/*extern void computeNextReplicate(pllInstance *tr, analdef *adef, int *originalRateCategories, int *originalInvariant);*/
extern void parseProteinModel(analdef *adef);
extern void catToGamma(pllInstance *tr, analdef *adef);
extern void gammaToCat(pllInstance *tr);
extern void calculateModelOffsets(pllInstance *tr);
extern void fixModelIndices(pllInstance *tr, int endsite, boolean fixRates);
extern void categorizeIterative(pllInstance *, int startIndex, int endIndex);
extern void evaluateGenericVectorIterative(pllInstance *, int startIndex, int endIndex);
extern double evaluateGenericInitravPartition(pllInstance *tr, nodeptr p, int model);

extern void determineFullTraversal(nodeptr p, pllInstance *tr);
/*extern void optRateCat(pllInstance *, int i, double lower_spacing, double upper_spacing, double *lhs);*/
extern boolean bootStop(pllInstance *tr, hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, unsigned int vectorLength);
extern void computeBootStopOnly(pllInstance *tr, char *bootStrapFileName, analdef *adef);
extern void computeConsensusOnly(pllInstance *tr, char* treeSetFileName, analdef *adef);

extern void printTreePerGene(pllInstance *tr, partitionList *pr, analdef *adef, char *fileName, char *permission);
extern boolean treeReadLenMULT ( FILE *fp, pllInstance *tr, analdef *adef );
extern void meshTreeSearch(pllInstance *tr, analdef *adef, int thorough);
extern void computeBOOTRAPID (pllInstance *tr, analdef *adef, long *radiusSeed);
extern void optimizeRAPID ( pllInstance *tr, analdef *adef );
extern void thoroughOptimization ( pllInstance *tr, analdef *adef, topolRELL_LIST *rl, int index );
extern int treeOptimizeThorough ( pllInstance *tr, int mintrav, int maxtrav);
extern int checker ( pllInstance *tr, nodeptr p );
extern boolean tipHomogeneityChecker ( pllInstance *tr, nodeptr p, int grouping );
extern void doAllInOne ( pllInstance *tr, analdef *adef );
extern void printBootstrapResult ( pllInstance *tr, analdef *adef, boolean finalPrint );
extern void printBipartitionResult ( pllInstance *tr, analdef *adef, boolean finalPrint );
extern unsigned int genericBitCount(unsigned int* bitVector, unsigned int bitVectorLength);
extern void computePlacementBias(pllInstance *tr, analdef *adef);
