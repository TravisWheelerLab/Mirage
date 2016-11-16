/* TransSW.c - Alex Nord - 2016
 *
 * USAGE:  ./TransSW <protein.fa> <dna.fa> [-debug]
 *
 * ABOUT:  This program takes a file containing a protein sequence
 *         and a file containing a DNA sequence (both FASTA-format)
 *         and outputs data used to identify highly identical
 *         alignments of the two sequences.  The program is built
 *         around an implementation of the Smith-Waterman local 
 *         sequence alignment algorithm.
 *
 *  NOTE:  The generic output of this program is extremely unfriendly
 *         to human readers---use the '-debug' option to see the
 *         same content with a little annotation.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Diagonals.h"


// Maximum and minimum macros
#define MAX(a,b) a > b ? a : b
#define MIN(a,b) a < b ? a : b


// Gap cost during Smith-Waterman
#define GAPCOST (-1)



/* FUNCTION: SortHits
 *
 *    ABOUT: This function generates an index that sorts the
 *           rightmost column of the Smith-Waterman matrix,
 *           so that we can quit after we've identified the
 *           highest-scoring run(s) through this segment.
 *           The backbone of the function is a simple mergesort
 *           implementation.
 */
int
SortHits
(
 int ** Source,
 int  * Target,
 int    x,
 int    TargetLength
)
{
  int i,j,k;
  
  for (i=0; i<TargetLength; i++)
    Target[i] = i;
  
  int * tempTarget;
  if ((tempTarget = malloc(TargetLength*sizeof(int))) == NULL) return 1;

  int blocksize,blocknum;
  int A,B,A_stop,B_stop;
  int placer;
  blocksize = 1;
  while (blocksize < TargetLength) {
    blocknum = 0;
    while ((blocknum+1) * blocksize < TargetLength) {
      A = blocksize * blocknum;
      B = A + blocksize;
      A_stop = B;
      B_stop = MIN(B+blocksize,TargetLength);
      placer = A;
      while (A < A_stop && B < B_stop) {
	if (Source[Target[A]][x] > Source[Target[B]][x]) {
	  tempTarget[placer] = Target[A];
	  A++;
	} else {
	  tempTarget[placer] = Target[B];
	  B++;
	}
	placer++;
      }
      while (A < A_stop) {
	tempTarget[placer] = Target[A];
	A++;
	placer++;
      }
      while (B < B_stop) {
	tempTarget[placer] = Target[B];
	B++;
	placer++;
      }
      blocknum += 2;
    }
    memcpy(Target,tempTarget,TargetLength*sizeof(int));
    blocksize *= 2;
  }

  if (tempTarget) free(tempTarget);
  
  return 0;
}



/* FUNCTION: TraceBack
 *
 *    ABOUT: This function traces a path backwards through a
 *           filled-out Smith-Waterman table, based on the
 *           scoring scheme.  It returns 1 if a path has been
 *           identified, and a 0 otherwise.
 */
void
TraceBack
(
 int ** DPtable,
 char * Protein,
 int    ProtEnd,
 char * TransLine,
 int    NuclEnd,
 int  * goodhit,
 int  * ProtStart,
 int  * NuclStart,
 int    debug
)
{
  *goodhit = 0;
  
  int i = NuclEnd;
  int j = ProtEnd;
  int k;
  
  int score = DPtable[i][j];
  int alignScore;
  
  //if (debug) printf("  ");
  
  while (score && i > 2) {
    
    *ProtStart = j;
    *NuclStart = i;

    // Overwriting this position to avoid redundant hits
    // NOTE: In order to avoid missing overlapping hits,
    //       we'll need to check for starts between the
    //       start+1 and end+1 (that end at least at end+1)
    //       during the stitching phase.
    DPtable[i][j] = 0;
    
    alignScore = CalcScore(TransLine[i],Protein[j-1]);
    
    //if (debug) printf("(%d,%d):%d  ",i,j,score);
    
    // 1. Clean match?
    if (score == DPtable[i-3][j-1] + alignScore) {
      i -= 3;
      j--;
      
      // 2. Unclean match - one nucleotide inserted
    } else if (score == DPtable[i-4][j-1] + alignScore + GAPCOST) {
      i -= 4;
      j--;
      
      // 3. Unclean match - two nucleotides inserted
    } else if (score == DPtable[i-5][j-1] + alignScore + GAPCOST) {
      i -= 5;
      j--;
      
      // 4. Prot. seq. insertion
    } else if (score == DPtable[i][j-1] + GAPCOST) {
      j--;
      
      // 5. Codon insertion
    } else if (score == DPtable[i-3][j] + GAPCOST) {
      i -= 3;
      
      // 6. Codon + 1 nucleotide insertion
    } else if (score == DPtable[i-4][j] + GAPCOST) {
      i -= 4;
      
      // 7. Codon + 2 nucleotides insertion
    } else if (score == DPtable[i-5][j] + GAPCOST) {
      i -= 5;
      
    } else { 

      // We're all washed up! :(
      //if (debug) printf("\n");
      return;

    }
    
    score = DPtable[i][j];
    
  }
  
  //if (debug) printf("\n\n\n");
  
  *goodhit = 1;

  return;

}




/* FUNCTION: 
 *
 *    ABOUT: 
 *
 */
int
TriangleSW
(
 char * Protein,
 int    ProtLength,
 int    debug
)
{
  int i,j,k;

  // Note that this COULD use reduced space (although
  // memory isn't expected to be a problem...  
  int ** DPtable;
  if ((DPtable = (int **)malloc(ProtLength*sizeof(int *))) == NULL) return 1;
  for (i=0; i<ProtLength; i++) 
    if ((DPtable[i] = (int *)malloc((ProtLength-i)*sizeof(int))) == NULL) return 1;
  
  // We're just looking for the max. value (as an indicator
  // of repetitiveness in the sequence)
  int max = 0;
  int max_i,max_j;
  
  // Initializing the topmost row
  DPtable[0][0] = 0;
  for (i=1; i<ProtLength; i++) {
    if (Protein[0] == Protein[i]) { // 0 is lame
      DPtable[0][i] = 1;
      max = 1;
    } else {
      DPtable[0][i] = 0;
    }
  }
  if (debug) {
    printf("     ");
    for(j=0; j<ProtLength; j++)
      printf("%c ",Protein[j]);
    printf("\n\n");
    printf("%c    ",Protein[0]);
    for (j=0; j<ProtLength; j++)
      printf("%d ",DPtable[0][j]);
    printf("\n");
  }
  
  // Running the recurrence <- Do we want to double-down on gaps/mismatches?
  //                           Currently we do (I think sensibly), but...
  int match;
  for (i=1; i<ProtLength; i++) {
    DPtable[i][0] = 0;
    for (j=1; j<ProtLength-i; j++) {

      if (Protein[i] == Protein[j+i]) match = 1;
      else                            match = -1;

      DPtable[i][j] = MAX(DPtable[i-1][j+1]-1,DPtable[i][j-1]-1);
      DPtable[i][j] = MAX(DPtable[i][j],match+DPtable[i-1][j]);
      DPtable[i][j] = MAX(DPtable[i][j],0);
      if (DPtable[i][j] > max) {
	max = DPtable[i][j];
	max_i = i;
	max_j = j;
      }
    }
    if (debug) {
      printf("%c    ",Protein[i]);
      for (j=ProtLength-i; j<ProtLength; j++)
	printf("  ");
      for (j=0; j<ProtLength-i; j++)
	printf("%d ",DPtable[i][j]);
      printf("\n");
    }
  }

  // Swag!  Tell us what ya saw
  printf("%d\n",max);

  // Mister Smith-Watermanchov, tear down that table!
  for (i=0; i<ProtLength; i++)
    free(DPtable[i]);
  free(DPtable);

  return 0;
}



/* FUNCTION: 
 *
 *    ABOUT: 
 *
 */
int
TranslatedSW
(
 char * Protein,
 int    ProtLength,
 char * NuclString,
 int    NuclLength,
 int    JustRepeatChecking,
 int    debug
)
{
  int i,j,k;
  
  // Allocate the dynamic programming table.
  int ** DPtable;
  if ((DPtable = (int **)malloc((3+NuclLength) * sizeof(int *))) == NULL)   return 1;
  for (i=0; i < 3+NuclLength; i++)
    if ((DPtable[i] = (int *)malloc((1+ProtLength) * sizeof(int))) == NULL) return 1;
  
  // Initialize the top row and leftmost column with 0s
  for (i=0; i < 3+NuclLength; i++)
    DPtable[i][0] = 0;
  for (i=0; i <= ProtLength; i++) {
    DPtable[0][i] = 0;
    DPtable[1][i] = 0;
    DPtable[2][i] = 0;
  }

  // We also want to allocate an array to hold the translated string corresponding
  // to the DP table.
  // Note that the amino acid corresponding to a nucleotide triple is placed at the
  // position following the end of the triple.  This is so that we can start at
  // the 'NuclLength'-th position during traceback, to make things more intuitive.
  char * TransLine;
  if ((TransLine = (char *)malloc((3+NuclLength) * sizeof(char))) == NULL) return 1;
  for (i=3; i <= NuclLength; i++)
    TransLine[i] = DNAtoAA(NuclString[i-3],NuclString[i-2],NuclString[i-1]);  
  
  
  // Fill in the dynamic programming table, such that each
  // position corresponds to the translated sequence 
  char  nextchar;
  int   nextval,lastval;
  int   max, max_i, max_j;
  max = max_i = max_j = 0;
  for (i=3; i <= NuclLength; i++) {
    if (debug) printf("  ");
    for (j=1; j <= ProtLength; j++) {
      nextchar = TransLine[i];
      // Because it's possible for single-nucleotide insertions to occur,
      // we have to consider each possible 'mistake' offset
      lastval = 0;

      // (1.) Amino acid insertion in protein string
      nextval = DPtable[i][j-1]+GAPCOST; 

      // (2.) Match states, accepting the possibility of single-nucleotide insertions,
      //      where 'k' is the number of nucleotides inserted into this codon
      for (k=0; k<3; k++) {
	if ((i-3)-k < 0) break;
	if (k) {
	  nextval = MAX(nextval, DPtable[(i-3)-k][j-1]+CalcScore(nextchar,Protein[j-1])+GAPCOST);
	} else {
	  nextval = MAX(nextval, DPtable[(i-3)-k][j-1]+CalcScore(nextchar,Protein[j-1]));
	}
	
	// (3.) Nucleotide insertions in DNA string.
	nextval = MAX(nextval, DPtable[(i-3)-k][j]+GAPCOST);
      }

      // Because this is Smith-Waterman, we max. with 0
      lastval = MAX(nextval, lastval);
      
      DPtable[i][j] = lastval;
      if (lastval > max) {
	max   = nextval;
	max_i = i;
	max_j = j;
      }
      if (debug) printf("%d ",DPtable[i][j]);
    }
    if (debug) printf("\n");
  }

  if (debug) {
    printf("\n");
    printf("  MAX: %d at (%d,%d)\n",max,max_i,max_j);
    printf("\n\n");
  }

  // If all we care about is checking whether the sequences have some good match, then
  // we only have to report the maximum here and be done with it.
  int * HitGuide;
  if (JustRepeatChecking) {
    
    printf("%d\n",max); 
    HitGuide = NULL;

  } else {

    // Next, we generate an index that sorts the final column in decreasing order,
    // such that we can try to track down any full-protein-length runs.
    if ((HitGuide = malloc((NuclLength-2)*sizeof(int))) == NULL) return 1;
    
    int no_hits = 1;
    int goodhit,ProtStart,NuclStart,highscore;

    
    k = ProtLength;
    while (k >= 0) {
      
      if (SortHits(DPtable,HitGuide,k,NuclLength-2)) return 1;
      
      // By order of terminal score, try to fill out a run through the protein sequence
      // segment.
      
      i = 0;
      highscore = DPtable[HitGuide[i]][k];
      if (highscore == 0) break;
      
      TraceBack(DPtable,Protein,k,TransLine,HitGuide[i],&goodhit,&ProtStart,&NuclStart,debug);
      while (i < NuclLength && DPtable[HitGuide[i]][k] && !goodhit) {
	i++;
	highscore = DPtable[HitGuide[i]][k];
	TraceBack(DPtable,Protein,ProtLength,TransLine,HitGuide[i],&goodhit,&ProtStart,&NuclStart,debug);
      }
      
      if (goodhit) {
	printf("%d-%d,%d-%d,%d\n",ProtStart-1,k-1,NuclStart-2,HitGuide[i],highscore);
	no_hits = 0;
      }
      
      // If we get a hit, we also try to run out all other hits //of equal score.
      j = 1;
      //while (i+j < NuclLength) { //&& DPtable[HitGuide[i+j]][k]) { // ONE
      while (i+j < NuclLength && DPtable[HitGuide[i+j]][k] == highscore) { // TWO
	TraceBack(DPtable,Protein,k,TransLine,HitGuide[i+j],&goodhit,&ProtStart,&NuclStart,debug);
	if (goodhit) printf("%d-%d,%d-%d,%d\n",ProtStart-1,k-1,(NuclStart-2),HitGuide[i+j],highscore);
	j++;
      }
      
      k--;
      
    }
    
  }

  // Relinquish your grasp on that memory! 
  // It has been fouled by the absence of utility!
  for (i=0; i<3+NuclLength; i++) {
    if (DPtable[i]) free(DPtable[i]);
  }
  if (HitGuide)  free(HitGuide);
  if (DPtable)   free(DPtable);
  if (TransLine) free(TransLine);

  return 0;
}




/* FUNCTION: PrintUsageMsg
 *
 *    ABOUT: This function gives the user a nudge in the right direction.
 *
 */
int PrintUsageMsg ()
{
  printf("\n\n");
  printf("\tUSAGE:  ./TransSW  <protein file>  <DNA file>  [OPTS]\n\n");
  printf("\tABOUT:  This program performs the basic Smith-Waterman algorithm\n");
  printf("\t        on a protein sequence and a translated DNA sequence.\n");
  printf("\t        Note: If using for against-self search, the 'DNA file' argument\n");
  printf("\t        still needs to be supplied, but will be ignored ('-' is fine).\n\n");
  printf("\tOPT.S:  -debug : Provide human-friendly output\n");
  printf("\t        -match : We're just interested in ID-ing the highest score in the table\n");
  printf("\t        -range : We're just searching some range of the protein\n");
  printf("\t        -self  : Search protein against itself, implies '-repeats'\n");
  printf("\n\n");
  return 1;
}



/* FUNCTION: main
 *
 *    ABOUT: Parse the input files, run the functions, do good things forever.
 *
 */
int main (int argc, char ** argv) 
{

  // In case the user seems confused
  if (argc < 3) return PrintUsageMsg();

  int i,j,k;

  // Checking for additional options
  i = 3; 
  int debug   = 0;
  int repeats = 0;
  int bounded = 0;
  int against_self = 0;
  int ProtStart, ProtEnd; // If we're using bounded search
  while (i < argc) {
    if (!strcmp(argv[i],"-debug")) {
      debug = 1;
    } else if (!strcmp(argv[i],"-match")) {
      repeats = 1;
    } else if (!strcmp(argv[i],"-range")) {
      bounded   = 1;
      ProtStart = atoi(argv[i+1]);
      ProtEnd   = atoi(argv[i+2]);
      i += 2;
    } else if (!strcmp(argv[i],"-self")) {
      repeats = 1;
      against_self = 1;
    } else {
      printf("  Unrecognized option '%s' ignored.\n",argv[i]);
    }
    i++;
  }

  // Open the protein FASTA file
  FILE * protfile;
  protfile = fopen(argv[1],"r");
  if (protfile == NULL) {
    printf("  ERROR:  Could not open file '%s'\n",argv[1]);
    return 1;
  }
  
  // A string used to rip lines out of the FASTA file.
  char * next_line;
  if ((next_line = malloc(128*sizeof(char))) == NULL) {
    printf("  ERROR:  Could not allocate string 'next_read'...?\n");
    return 1;
  }
  
  // Because we don't want to risk weird formatting breaking everything,
  // we eat the header line char by char.
  while (fgetc(protfile) >= 32);
  
  // Prep the protein string  
  int ProtLength;
  char * Protein;
  int capacity = 512;

  // If the sequence is bounded, we know exactly how large our array will be
  if (bounded) { 
  
    ProtLength = (ProtEnd-ProtStart)+1;
    if ((Protein = malloc(ProtLength*sizeof(char))) == NULL) {
      printf("  ERROR:  Could not allocate string 'Protein'!\n");
      return 1;
    } else {
      for (i=0; i<ProtLength; i++)
	Protein[i] = 0;
    }
    capacity = ProtLength;
    
  } else { // More like the nucleotide sequence rip
    
    if ((Protein = malloc(capacity*sizeof(char))) == NULL) {
      printf("  ERROR: Could not allocate string 'Protein'!\n");
      return 1;
    } else {
      for (i=0; i<capacity; i++) 
	Protein[i] = 0;
    }

  }
  
  // Rip the protein file, keeping an eye on all the junk you need to.
  char * ptemp = NULL;
  int pos = 0;
  ProtLength = 0;
  while (!feof(protfile)) {    
    fscanf(protfile,"%s\n",next_line);
    i = 0;
    while (next_line[i] > 32) {
      if (!bounded || (pos >= ProtStart && pos <= ProtEnd)) {
	Protein[ProtLength] = next_line[i];
	ProtLength++;
      }
      if (!bounded && ProtLength == capacity) {
	if (ptemp != NULL) {
	  if ((ptemp = realloc(ptemp,capacity*sizeof(char))) == NULL) {
	    printf("  ERROR:  Failed to reallocate string 'ptemp'\n");
	    return 1;
	  }
	} else {
	  if ((ptemp = malloc(capacity*sizeof(char))) == NULL) {
	    printf("  ERROR:  Failed to allocate string 'ptemp'\n");
	    return 1;
	  }
	}
	memcpy(ptemp,Protein,capacity*sizeof(char));
	if ((Protein = realloc(Protein,capacity*2*sizeof(char))) == NULL) {
	  printf("  ERROR: Failed to reallocate string 'Protein'!\n");
	  return 1;
	}
	memcpy(Protein,ptemp,capacity*sizeof(char));
	for (j=capacity; j<capacity*2; j++) 
	  Protein[j] = 0;
	capacity *= 2;
      }
      i++;
      pos++;
    }
  }
  if (ptemp) free(ptemp);
  fclose(protfile);

  // If we're searching against ourself, we can just re-use the Protein string
  int NuclLength = 0;
  char * NuclString;
  if (against_self) {
  
    NuclString = NULL;
    NuclLength = ProtLength;

  } else {

    // Open up the DNA/RNA string
    FILE * nuclfile;
    nuclfile = fopen(argv[2],"r");
    if (nuclfile == NULL) {
      printf("  ERROR:  Could not open file '%s'\n",argv[2]);
      fclose(protfile);
      return 1;
    }
    
    // Because we don't want to risk weird formatting breaking everything,
    // we eat the header line char by char.
    while (fgetc(nuclfile) >= 32);
    
    // Prep the DNA string
    capacity   = 512;
    if ((NuclString = malloc(capacity*sizeof(char))) == NULL) {
      printf("  ERROR:  Could not allocate string 'NuclString'!\n");
      return 1;
    } else {
      for (i=0; i<capacity; i++)
	NuclString[i] = 0;
    }
    
    // Rip that DNA file!
    char * ntemp = NULL;
    while (!feof(nuclfile)) {    
      fscanf(nuclfile,"%s\n",next_line);
      i = 0;
      while (next_line[i] > 32) {
	NuclString[NuclLength] = next_line[i];
	NuclLength++;
	i++;
	if (NuclLength == capacity) {	
	  if (ntemp != NULL) {
	    if ((ntemp = realloc(ntemp,capacity*sizeof(char))) == NULL) {
	      printf("  ERROR:  Failed to reallocate string 'ntemp'\n");
	      return 1;
	    }
	  } else {
	    if ((ntemp = malloc(capacity*sizeof(char))) == NULL) {
	      printf("  ERROR:  Could not allocate string 'ntemp'\n");
	      return 1;
	    }
	  }	
	  memcpy(ntemp,NuclString,capacity*sizeof(char));
	  if ((NuclString = realloc(NuclString,capacity*2*sizeof(char))) == NULL) {
	    printf("  ERROR:  Failed to reallocate string 'NuclString'!\n");
	    return 1;
	  }
	  memcpy(NuclString,ntemp,capacity*sizeof(char));	
	  for (j=capacity; j<capacity*2; j++)
	    NuclString[j] = 0;
	  capacity *= 2;
	}
      }
    }    
    if (ntemp) free(ntemp);
    fclose(nuclfile);
  }

  // No longer need next_line
  free(next_line);
  
  // If we're debugging, we want to make sure that we read the files in correctly
  if (debug) {
    printf("\n");
    printf(" Input Protein: %s\n",Protein);
    if (!against_self) 
      printf("     Input DNA: %s\n",NuclString);
    printf("\n");
  }
  
  // Because single nucleotides can be inserted, we need to do a DP
  // table including all possible offsets -- NOT ONE DP TABLE PER!
  if (against_self) {
    if (TriangleSW(Protein,ProtLength,debug)) {
      printf("\n\tFAILURE at search 1\n\n");
      return 1;
    } 
  } else {
    if (TranslatedSW(Protein,ProtLength,NuclString,NuclLength,repeats,debug)) {
      printf("\n\tFAILURE at search 1\n\n");
      return 1;
    }
  }
  
  // Get outta here, you rascals!
  if (Protein)    free(Protein);
  if (NuclString) free(NuclString);
  
  // No problem / 2ez / gg
  return 0;

}
