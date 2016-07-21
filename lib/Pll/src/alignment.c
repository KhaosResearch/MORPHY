/** 
 * PLL (version 1.0.0) a software library for phylogenetic inference
 * Copyright (C) 2013 Tomas Flouri and Alexandros Stamatakis
 *
 * Derived from 
 * RAxML-HPC, a program for sequential and parallel estimation of phylogenetic
 * trees by Alexandros Stamatakis
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
 * For any other enquiries send an Email to Tomas Flouri
 * Tomas.Flouri@h-its.org
 *
 * When publishing work that uses PLL please cite PLL
 * 
 * @file alignment.c
 *
 * @brief Collection of routines for reading alignments
 *
 * Auxiliary functions for storing alignments read from predefined file formats
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pll.h"
#include "pllInternal.h"

/** @defgroup alignmentGroup Reading and parsing multiple sequence alignments
    
    This set of functions handles the reading and parsing of several file formats that describe multiple sequence alignments. They are also responsible for storing the alignment in an internal structure
*/
static pllAlignmentData * pllParsePHYLIP (const char * filename);
static pllAlignmentData * pllParseFASTA (const char * filename);
static int read_phylip_header (int * inp, int * sequenceCount, int * sequenceLength);
static inline int parsedOk (int * actLen, int sequenceCount, int sequenceLength);
static int parse_phylip (pllAlignmentData * alignmentData, int input);
//static int getFastaAlignmentInfo (int * inp, int * seqCount, int * seqLen);
//static int parseFastaAlignment (pllAlignmentData * alignmentData, int input);

static int query_getnext(char ** head, int * head_len, 
                         char ** seq, int * seq_len, int * qno);
static int query_open(const char * filename);
static void query_close(void);

#define PLL_MEMCHUNK    4096
#define PLL_LINEALLOC   1048576

static FILE * query_fp;
static char query_line[PLL_LINEALLOC];

static int query_no = -1;
static char * query_head = 0;
static char * query_seq = 0;

static long query_head_len = 0;
static long query_seq_len = 0;

static long query_head_alloc = 0;
static long query_seq_alloc = 0;

static long query_filesize = 0;
static int query_lineno;

static unsigned int chrstatus[256] =
  {
    /*

      How to handle input characters

      0=stripped, 1=legal, 2=fatal, 3=silently stripped

    @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
    P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

    0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,  2,  2,  3,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,
    1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
    0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,
    0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
  };

#ifdef __PLL_DEBUG_PARSER
static int
printTokens (int input)
{
  pllLexToken token;

  do
   {
     NEXT_TOKEN

     /* begin of parser */
     switch (token.tokenType)
      {
        case PLL_TOKEN_NUMBER:
          printf ("PLL_TOKEN_NUMBER (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case PLL_TOKEN_STRING:
          printf ("PLL_TOKEN_STRING (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case PLL_TOKEN_EOF:
          printf ("PLL_TOKEN_EOF\n");
          break;
        case PLL_TOKEN_WHITESPACE:
          printf ("PLL_TOKEN_WHITESPACE\n");
          break;
        case PLL_TOKEN_NEWLINE:
          printf ("PLL_TOKEN_NEWLINE\n");
          break;
        case PLL_TOKEN_UNKNOWN:
          printf ("PLL_TOKEN_UNKNOWN (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        default:
          break;
      }
     /* end of parser */


   }
  while (token.tokenType != PLL_TOKEN_EOF && token.tokenType != PLL_TOKEN_UNKNOWN);

  if (token.tokenType == PLL_TOKEN_UNKNOWN) return (0);

  return (1);
}
#endif

/** @ingroup alignmentGroup
    @brief Initialize alignment structure fields

    Allocates memory for the data structure that will hold the alignment and
    initializes it. It requires the number of sequences \a sequenceCount and
    the length of sequences \a sequenceLength. It returns a pointer to the
    initialized data structure.

    @param sequenceCount
      Number of sequences in the alignment
    
    @param sequenceLength
      Length of the sequences

    @param 
      Initialized alignment data structured
*/
pllAlignmentData *
pllInitAlignmentData (int sequenceCount, int sequenceLength)
 {
   int i;
   pllAlignmentData * alignmentData;
   void * mem;
   
   /** TODO */
   alignmentData               =  (pllAlignmentData *) malloc (sizeof (pllAlignmentData));
   alignmentData->sequenceData = (unsigned char **) malloc ((sequenceCount + 1) * sizeof (unsigned char *));
   mem = (void *) malloc (sizeof (unsigned char) * (sequenceLength + 1) * sequenceCount);
   for (i = 1; i <= sequenceCount; ++i)
    {
      alignmentData->sequenceData[i]                 = (unsigned char *) (mem + (i - 1) * (sequenceLength + 1) * sizeof (unsigned char));
      alignmentData->sequenceData[i][sequenceLength] = 0;
    }
   alignmentData->sequenceData[0] = NULL;
    
   alignmentData->sequenceLabels = (char **) calloc ((sequenceCount + 1), sizeof (char *));

   alignmentData->sequenceCount  = sequenceCount;
   alignmentData->sequenceLength = sequenceLength;
   alignmentData->originalSeqLength = sequenceLength;

   /** TODO: remove siteWeights from alignment */
   alignmentData->siteWeights    = NULL;

   return (alignmentData);
 }

/** @ingroup alignmentGroup
    @brief Deallocates the memory associated with the alignment data structure
    
    Deallocates the memory associated with the alignment data structure \a alignmentData.

    @param alignmentData
      The alignment data structure
*/
void
pllAlignmentDataDestroy (pllAlignmentData * alignmentData)
{
  int i;

  for (i = 1; i <= alignmentData->sequenceCount; ++ i)
   {
     rax_free (alignmentData->sequenceLabels[i]);
   }
  rax_free (alignmentData->sequenceData[1]);
  rax_free (alignmentData->sequenceLabels);
  rax_free (alignmentData->sequenceData);
  rax_free (alignmentData->siteWeights);
  rax_free (alignmentData);
}


/** @ingroup alignmentGroup
    @brief Prints the alignment to the console

    @param alignmentData
      The alignment data structure
*/
void pllAlignmentDataDumpConsole (pllAlignmentData * alignmentData)
 {
   int i;

   printf ("%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
   for (i = 1; i <= alignmentData->sequenceCount; ++ i)
    {
      printf ("%s %s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
    }
 }



static void dump_fasta_content(FILE * fp, pllAlignmentData * alignmentData)
{
  int i;

  for (i = 1; i <= alignmentData->sequenceCount; ++i)
     fprintf (fp, ">%s\n%s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
}

static void dump_phylip_content(FILE * fp, pllAlignmentData * alignmentData)
{
  int i;

  for (i = 1; i <= alignmentData->sequenceCount; ++i)
     fprintf (fp, "%s %s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
}

/** @ingroup alignmentGroup
    @brief Dump the alignment to a file of format \a fileFormat

    Dumps the alignment contained in \a alignmentData to file \a filename of type \a fileFormat.

    @note If \a filename exists, all contents will be erased

    @param alignmentData
      Alignment data structure

    @param fileFormat
      Format of output file. Can take the value \b PLL_FORMAT_PHYLIP or \b PLL_FORMAT_FASTA

    @param filename
      Output filename

    @return
      Returns \b PLL_TRUE on success, otherwise \b PLL_FALSE.
*/
int
pllAlignmentDataDumpFile (pllAlignmentData * alignmentData, int fileFormat, const char * filename)
{
  FILE * fp;
  void (*outfun)(FILE *, pllAlignmentData *);
  
  if (fileFormat != PLL_FORMAT_PHYLIP && fileFormat != PLL_FORMAT_FASTA) return (PLL_FALSE);

  outfun = (fileFormat == PLL_FORMAT_PHYLIP) ? dump_phylip_content : dump_fasta_content;

  fp = fopen (filename, "w");
  if (!fp) return (PLL_FALSE);
  
  /* if PHYLIP print the silly header at the beginning */
  if (fileFormat == PLL_FORMAT_PHYLIP)
   {
     fprintf (fp, "%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
   }
  
  outfun(fp, alignmentData);

  fclose (fp);
  return (PLL_TRUE);
}



/* ROUTINES FOR PHYLIP PARSING */
/** @ingroup alignmentGroup
    @brief Parse the PHYLIP file header
*/
static int
read_phylip_header (int * inp, int * sequenceCount, int * sequenceLength)
{
  pllLexToken token;
  int input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

  if (token.tokenType != PLL_TOKEN_NUMBER) return (0);

  *sequenceCount = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
  if (token.tokenType != PLL_TOKEN_NUMBER) return (0);

  *sequenceLength = atoi (token.lexeme);

  *inp = input;

  return (*sequenceCount && *sequenceLength);
}

static inline int
parsedOk (int * actLen, int sequenceCount, int sequenceLength)
{
  int i;

  for (i = 1; i <= sequenceCount; ++ i)
   {
     if (actLen[i] != sequenceLength) return (0);
   }
  
  return (1);
}


/** @ingroup alignmentGroup
    @brief Parse the PHYLIP file body
*/
static int
parse_phylip (pllAlignmentData * alignmentData, int input)
{
  int i,j;
  pllLexToken token;
  int * sequenceLength;
  int rc;

  sequenceLength = (int *) rax_calloc (alignmentData->sequenceCount + 1, sizeof (int));

  NEXT_TOKEN
  for (i = 0; ; ++i)
  {
    j = i % alignmentData->sequenceCount;
    if (i < alignmentData->sequenceCount) 
     {
       if (token.tokenType == PLL_TOKEN_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.tokenType == PLL_TOKEN_UNKNOWN)
        {
          rax_free (sequenceLength);
          return (0);
        }

       CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)


       if (token.tokenType != PLL_TOKEN_STRING && token.tokenType != PLL_TOKEN_NUMBER && token.tokenType != PLL_TOKEN_FLOAT)
        {
          rax_free (sequenceLength);
          return (0);
        }
       alignmentData->sequenceLabels[i + 1] = strndup (token.lexeme, token.len);
       NEXT_TOKEN
       CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
     }
    
    while (1)
     {
       if (token.tokenType == PLL_TOKEN_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.tokenType == PLL_TOKEN_UNKNOWN)
        {
         rax_free (sequenceLength);
         return (0);
        }
       
       if (token.tokenType == PLL_TOKEN_NEWLINE) break;

       if (token.tokenType != PLL_TOKEN_STRING)
        {
          rax_free (sequenceLength);
          return (0);
        }

       if (sequenceLength[j + 1] + token.len > alignmentData->sequenceLength) 
        {
          fprintf (stderr, "Sequence %d is larger than specified\n", j + 1);
          rax_free (sequenceLength);
          return (0);
        }
       memmove (alignmentData->sequenceData[j + 1] + sequenceLength[j + 1], token.lexeme, token.len);
       sequenceLength[j + 1] += token.len;

       NEXT_TOKEN
       CONSUME (PLL_TOKEN_WHITESPACE)
     }
    CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE);
  }
}

/* Phylip parsers. Use the following attributed grammar 
 * 
 *        S -> HEADER ENDL DATA
 *   HEADER -> PLL_TOKEN_NUMBER PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER ENDL |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER ENDL
 *     ENDL -> PLL_TOKEN_WHITESPACE PLL_TOKEN_NEWLINE | PLL_TOKEN_NEWLINE
 *     DATA -> PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING ENDL DATA |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING ENDL DATA | 
 *             PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_EOF |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_EOF
 */

/** @ingroup alignmentGroup
    @brief Parse a PHYLIP file

    Parses the PHYLIP file \a filename and returns a ::pllAlignmentData structure
    with the alignment.

    @param filename
      Name of file to be parsed

    @return
      Returns a structure of type ::pllAlignmentData that contains the alignment, or \b NULL
      in case of failure.
*/
static pllAlignmentData *
pllParsePHYLIP (const char * filename)
{
  int i, input, sequenceCount, sequenceLength;
  char * rawdata;
  long filesize;
  pllAlignmentData * alignmentData;

  rawdata = pllReadFile (filename, &filesize);
  if (!rawdata)
   {
     errno = PLL_ERROR_FILE_OPEN;
     return (NULL);
   }
  
  init_lexan (rawdata, filesize);
  input = get_next_symbol();

  /* parse the header to obtain the number of taxa and sequence length */
  if (!read_phylip_header (&input, &sequenceCount, &sequenceLength))
   {
     rax_free (rawdata);
     fprintf (stderr, "Error while parsing PHYLIP header (number of taxa and sequence length)\n");
     errno = PLL_ERROR_PHYLIP_HEADER_SYNTAX;
     return (NULL);
   }

  lex_table_amend_phylip();

  /* allocate alignment structure */
  alignmentData = pllInitAlignmentData (sequenceCount, sequenceLength);

  if (! parse_phylip (alignmentData, input))
   {
     errno = PLL_ERROR_PHYLIP_BODY_SYNTAX;
     pllAlignmentDataDestroy (alignmentData);
     lex_table_restore();
     rax_free (rawdata);
     return (NULL);
   }
  
  lex_table_restore();
  rax_free (rawdata);

  alignmentData->siteWeights  = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i) 
    alignmentData->siteWeights[i] = 1;

  return (alignmentData);
}


/* FASTA routines */

static void query_close(void)
{
  fclose(query_fp);
  
  if (query_seq)
    free(query_seq);
  if (query_head)
    free(query_head);

  query_head = 0;
  query_seq = 0;
}

static int query_open(const char * filename)
{
  query_head = NULL;
  query_seq = NULL;

  query_head_len = 0;
  query_seq_len = 0;

  query_head_alloc = PLL_MEMCHUNK;
  query_seq_alloc = PLL_MEMCHUNK;

  query_head = (char *) malloc((size_t)query_head_alloc);
  if (!query_head) return 0;
  query_seq = (char *) malloc((size_t)query_seq_alloc);
  if (!query_seq) 
  {
    free(query_head);
    return 0;
  }

  query_no = -1;

  /* open query file */
  query_fp = NULL;
  query_fp = fopen(filename, "r");
  if (!query_fp) 
  {
    free(query_head);
    free(query_seq);
    return 0;
  }

  if (fseek(query_fp, 0, SEEK_END)) 
  {
    free(query_head);
    free(query_seq);
    return 0;
  }

  query_filesize = ftell(query_fp);
  
  rewind(query_fp);

  query_line[0] = 0;
  fgets(query_line, PLL_LINEALLOC, query_fp);
  query_lineno = 1;

  return 1;
}

static int query_getnext(char ** head, int * head_len,
                  char ** seq, int * seq_len, int * qno)
{
  while (query_line[0])
    {
      /* read header */

      if (query_line[0] != '>')
        return -1;
      
      if (strlen(query_line) + 1 == PLL_LINEALLOC)
        return -1;

      /* terminate header at first space or end of line */

      char * z0 = query_line + 1;
      char * z = z0;
      while (*z)
        {
          if (*z == '\n')
            break;
          z++;
        }
      long headerlen = z - z0;
      query_head_len = headerlen;

      /* store the header */

      if (headerlen + 1 > query_head_alloc)
        {
          query_head_alloc = headerlen + 1;
          query_head = (char *) realloc(query_head, (size_t)query_head_alloc);
          if (!query_head) return -1;
        }

      memcpy(query_head, query_line + 1, (size_t)headerlen);
      query_head[headerlen] = 0;

      /* get next line */

      query_line[0] = 0;
      fgets(query_line, PLL_LINEALLOC, query_fp);
      query_lineno++;

      /* read sequence */

      query_seq_len = 0;

      while (query_line[0] && (query_line[0] != '>'))
        {
          char c;
          char m;
          char * p = query_line;

          while((c = *p++))
            {
              m = chrstatus[(int)c];
              switch(m)
                {
                case 0:
                  /* fatal character */
                  return -1;
                  break;

                case 1:
                  /* legal character */
                  if (query_seq_len + 1 > query_seq_alloc)
                    {
                      query_seq_alloc += PLL_MEMCHUNK;
                      query_seq = (char *) realloc(query_seq, (size_t)query_seq_alloc);
                      if (!query_seq) return -1;
                    }
                  *(query_seq + query_seq_len) = c;
                  query_seq_len++;

                  break;

                case 2:
                  /* silently stripped chars */
                  break;
                default:
                  return -1;

                }
            }

          query_line[0] = 0;
          fgets(query_line, PLL_LINEALLOC, query_fp);
          query_lineno++;
        }

      /* add zero after sequence */

      if (query_seq_len + 1 > query_seq_alloc)
        {
          query_seq_alloc += PLL_MEMCHUNK;
          query_seq = (char *) realloc(query_seq, (size_t)query_seq_alloc);
          if (!query_seq) return -1;
        }
      *(query_seq + query_seq_len) = 0;

      query_no++;
      *head = query_head;
      *seq = query_seq;
      *head_len = query_head_len;
      *seq_len = query_seq_len;
      *qno = query_no;

      return 1;
    }
  
  return 0;
}

static void _pll_free_temp_alignment(pllAlignmentData * alignmentData, int count)
{
  int i;

  for (i = 1; i <= count; ++i)
  {
    free(alignmentData->sequenceData[i]);
    free(alignmentData->sequenceLabels[i]);
  }
  free(alignmentData->sequenceData);
  free(alignmentData->sequenceLabels);
  free(alignmentData);
}

static pllAlignmentData *
pllParseFASTA (const char * filename)
{
  int i,j;
  int status;

  pllAlignmentData * alignmentData;

  if (!query_open(filename))
  {
    errno = PLL_ERROR_FILE_OPEN;
    return NULL;
  }

  alignmentData = (pllAlignmentData *) malloc(sizeof(pllAlignmentData));
  alignmentData->sequenceData = NULL;
  alignmentData->sequenceLabels = NULL;

  int prev_qseqlen = -1;
  while(1)
  {
    char * qhead;
    int query_head_len;
    char *qseq;
    int qseqlen;
//    int query_no = -1;

    if ((status = query_getnext(&qhead, &query_head_len,
                      &qseq, &qseqlen, &query_no)) > 0)
    {
      alignmentData->sequenceData = (unsigned char **)realloc(alignmentData->sequenceData, 
                                                              (query_no + 2)*sizeof(unsigned char *));
      alignmentData->sequenceLabels = (char **)realloc(alignmentData->sequenceLabels,
                                                                (query_no + 2)*sizeof(char *));

      /* remove trailing whitespace from sequence names */
      j = query_head_len-1;
      while(j>0 && (qhead[j] == ' ' || qhead[j] == '\t'))
      {
        qhead[j] = 0;
      }
      alignmentData->sequenceData[query_no+1] = (unsigned char *)strdup(qseq);
      alignmentData->sequenceLabels[query_no+1] = strdup(qhead);
      
      if (prev_qseqlen != -1)
      {
        /* fasta sequences not aligned, free everything except last read
           which will be freed by query_close() */
        if (qseqlen != prev_qseqlen)
        {
          errno = PLL_ERROR_FASTA_SYNTAX;
          _pll_free_temp_alignment(alignmentData, query_no+1);
          query_close();
          return NULL;
        }
      }
      else
      {
        alignmentData->sequenceLength = qseqlen;
        alignmentData->originalSeqLength = qseqlen;
        prev_qseqlen = qseqlen;
      }
    }
    else if (status == -1)
    {
      errno = PLL_ERROR_FASTA_SYNTAX;
      _pll_free_temp_alignment(alignmentData, query_no+1);
      query_close();
      return NULL;
    }
    else break;
  }
  alignmentData->sequenceCount  = query_no+1;
  query_close();

  alignmentData->siteWeights = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i)
    alignmentData->siteWeights[i] = 1;

  /* ugly hack to turn it to one contiguous block of memory. This should be redesigned */
  void * mem = malloc((alignmentData->sequenceCount)*(alignmentData->sequenceLength+1)*sizeof(unsigned char));
  for (i = 1; i <= alignmentData->sequenceCount; ++i)
  {
    void * tmp = alignmentData->sequenceData[i];
    alignmentData->sequenceData[i] = (unsigned char *) (mem + (i - 1) * (alignmentData->sequenceLength + 1) * sizeof (unsigned char));
    memcpy(alignmentData->sequenceData[i], tmp, alignmentData->sequenceLength);
    alignmentData->sequenceData[i][alignmentData->sequenceLength] = 0;
    free(tmp);
  }
  alignmentData->sequenceData[0] = NULL;

  return (alignmentData);
}


/** @ingroup alignmentGroup
    @brief Parse a file that contains a multiple sequence alignment

    Parses the file \a filename of type \a fileType which contains a multiple sequence alignment.
    The supported file types are the sequential and interleaved versions of PHYLIP format, and
    the FASTA format. The parsed alignment is returned as a pointer to a structure of type
    ::pllAlignmentData

    @param fileType
      Type of file to parse. Can be either \b PLL_FORMAT_PHYLIP or \b PLL_FORMAT_FASTA

    @param filename
      Name of file to parse

    @return
      Returns a structure of type ::pllAlignmentData that contains the multiple sequence alignment,
      otherwise returns \b NULL in case of failure.
*/
pllAlignmentData * pllParseAlignmentFile (int fileType, const char * filename)
{
  switch (fileType)
   {
     case PLL_FORMAT_PHYLIP:
       return (pllParsePHYLIP (filename));
     case PLL_FORMAT_FASTA:
       return (pllParseFASTA (filename));
     default:
       /* RTFM */
       errno = PLL_ERROR_INVALID_FILETYPE;
       return (NULL);
   }

}
