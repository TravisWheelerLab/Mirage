/* MultiSeqNW.c - Alex Nord - 2016
 *
 * USAGE:  ./MultiSeqNW <MSA one> <MSA two> [-debug]
 *
 * ABOUT:  This program takes two multiple sequence alignments (MSAs)
 *         and performs a profile-profile global alignment based on
 *         Needleman-Wunsch.
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MultiSeqNW.h"


// Maximum and minimum macros -- UNSAFE
//#define MAX(a,b) a > b ? a : b
//#define MIN(a,b) a < b ? a : b


// The costs of generating gaps
/*
#define INTRON_GAP_BASE (-5.0)   // Wher we start our intron gapstart cost
#define INTRON_GAP_MAX  (-200.0) // The most we'll charge for an intron gap
#define INTRON_GAP_CONT (0.0)    // Cost of continuing a gap over an intron 
#define GAP_START       (-11.0)  // Affine gap open cost
#define GAP_CONTINUE    (-1.0)   // Affine gap continuation
// The costs of matching when working with special chars
#define INTRON_TO_INTRON (0.0)
#define NEGATIVE_INF     (0.0-INFINITY)
*/


// The BLOSUM 62 scoring matrix.  Because we're performing global alignment,
// our choice of scoring matrix is more or less arbitrary.  Scores are in
// half-bits.
const float BLOSUM62[400] = {
   4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,
  -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3,
  -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,
  -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,
   0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
  -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,
  -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,
   0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,
  -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,
  -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,
  -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,
  -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,
  -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,
  -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,
  -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,
   1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,
   0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,
  -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,
  -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,
   0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4,
};




float Min (float a, float b)
{
  if (a < b) return a;
  return b;
}


float Max (float a, float b)
{
  if (a > b) return a;
  return b;
}






/* FUNCTION: InitTupleSet
 *
 *    ABOUT: This function initializes a 'TUPLE_SET' (i.e., profile
 *           position) object, returning the object on successful
 *           initialization or else NULL.
 *
 */
TUPLE_SET *
InitTupleSet 
(
 int num_letters, 
 int is_non_letter
)
{
  TUPLE_SET * ts;

  if ((ts = malloc(sizeof(TUPLE_SET))) == NULL)
    return NULL;

  ts->MarksIntron    = 0;
  ts->NumTuples      = 0;
  ts->NearestIntron  = 0;

  if (is_non_letter) {

    ts->TupleChars  = NULL;
    ts->TupleRatios = NULL;
    if (is_non_letter == 1) 
      ts->MarksIntron = 1;

  } else {

    if ((ts->TupleChars = malloc(num_letters * sizeof(char))) == NULL) {
      free(ts);
      return NULL;
    }
    if ((ts->TupleRatios = malloc(num_letters * sizeof(float))) == NULL) {
      free(ts->TupleChars);
      free(ts);
      return NULL;
    }
    ts->NumTuples = num_letters;

  }

  return ts;

}



/* FUNCTION: ConvertToTupleSet
 *
 *    ABOUT: This function takes a set of letters and builds a
 *           TUPLE_SET object representing the expression of
 *           those characters.  These characters correspond to
 *           some column in one of the MSAs.
 *
 */
TUPLE_SET *
ConvertToTupleSet 
(
 char * letters, 
 int num_letters
)
{
  int i;

  TUPLE_SET * ts;

  // We do a start-up check to make sure that this isn't a
  // gap or intron position.
  for (i=0; i<num_letters; i++) {
    if (letters[i] == '*') {
      ts = InitTupleSet(0,1);
      return ts;
    }
  }

  // Allocate a full alphabet of space
  int * Alphabet;
  if ((Alphabet = malloc(26*sizeof(int))) == NULL)
    return NULL;
  for (i=0; i<26; i++)
    Alphabet[i] = 0;
  
  // How many actual (non-gap) characters occur?
  int actual_chars = 0;
  for (i=0; i<num_letters; i++) {
    if (letters[i] == '-')
      continue;
    if (letters[i] >= 97) 
      letters[i] -= 32;
    Alphabet[letters[i]-65]++;
    actual_chars++;
  }
  
  // Figure out some stats for our profile
  int distinct_letters = 0;
  for (i=0; i<26; i++) {
    if (Alphabet[i]) 
      distinct_letters++;
  }

  // And initialize this position of the profile
  ts = InitTupleSet(distinct_letters,0);
  int placer = 0;
  for (i=0; i<26; i++) {
    if (Alphabet[i]) {
      ts->TupleChars[placer]  = i+65;
      ts->TupleRatios[placer] = ((float)Alphabet[i])/((float)actual_chars);
      placer++;
    }
  }

  // Peace out
  free(Alphabet);
  return ts;

}






/* FUNCTION: DestroyTupleSet
 *
 *    ABOUT: See, you take this here TUPLE_SET object here and make
 *           it not be here anymore, see.
 *
 */
void DestroyTupleSet (TUPLE_SET * ts)
{
  if (ts->TupleChars)  free(ts->TupleChars);
  if (ts->TupleRatios) free(ts->TupleRatios);
  free(ts);
}





/* FUNCTION: ConvertToIndex
 *
 *    ABOUT: Take a character and figure out where its index is
 *           in BLOSUM62.
 *
 */
int
ConvertToIndex (char letter)
{
  // Capitalize (should be unnecessary, but better safe than sorry)
  if (letter >= 97) 
    letter -= 32;

  // Convert
  switch(letter) {
  case 'A':
    return 0;
  case 'R':
    return 1;
  case 'N':
    return 2;
  case 'D':
    return 3;
  case 'C':
    return 4;
  case 'Q':
    return 5;
  case 'E':
    return 6;
  case 'G':
    return 7;
  case 'H':
    return 8;
  case 'I':
    return 9;
  case 'L':
    return 10;
  case 'K':
    return 11;
  case 'M':
    return 12;
  case 'F':
    return 13;
  case 'P':
    return 14;
  case 'S':
    return 15;
  case 'T':
    return 16;
  case 'W':
    return 17;
  case 'Y':
    return 18;
  case 'V':
    return 19;
  default:
    return -1;
  }
}




/* FUNCTION: CalcGap
 *
 *    ABOUT: Calculate the cost of starting a gap in our DP matrix
 *
 */
float
CalcGap
(
 float       basescore,
 TUPLE_SET * Tuple1,
 TUPLE_SET * Tuple2,
 int         starting,
 int         is_insert, // i.e., Tuple1 is actually a position back
 float       intron_gap_base,
 float       intron_gap_max,
 float       intron_gap_cont,
 float       intron_gap_mult,
 float       gap_start,
 float       gap_continue,
 float       negative_inf
)
{

  // Are we kicking off a gap?
  if (starting) {

    // Are either of these tuples introns?
    if (Tuple1->MarksIntron || Tuple2->MarksIntron) {
      
      // Are BOTH introns?
      if (Tuple1->MarksIntron && Tuple2->MarksIntron) {
	return basescore+gap_start;
      }
      
      // Didn't think so.  Need to ID intron and calc. dist. for friendo.
      // Note that 'INTRON_GAP_MAX' is the largest cost in the sense of "damage done"
      // rather than lowest-value.
      if (Tuple1->MarksIntron && !is_insert) {

	if (Tuple2->NearestIntron > 10 && intron_gap_base > 0) 
	  return basescore+intron_gap_max; // Avoiding overflow
	return basescore+Max(intron_gap_max,intron_gap_base * pow(intron_gap_mult,Tuple2->NearestIntron-1));

      }

      // Either Tuple2 is the intron or we're going the wrong way
      if (Tuple2->MarksIntron && is_insert) {
	if (Tuple1->NearestIntron > 10 && intron_gap_base > 0) 
	  return basescore+intron_gap_max; // Avoiding overflow
	return basescore+Max(intron_gap_max,(intron_gap_base*pow(intron_gap_mult,Tuple1->NearestIntron-1))-1);

      }

    }
  
    // If neither of the special cases panned out, we do a standard start
    return basescore+gap_start;

  }

  // Nope, just continuing a gap.
  // How dreadfuly dull.
  return basescore+gap_continue;

}





/* FUNCTION: CalcMatch
 *
 *    ABOUT: Calculate the score of aligning two positions in
 *           a profile, with special considerations for gaps and
 *           introns.
 *
 */
float
CalcMatch 
(
 float       basescore,
 TUPLE_SET * Tuple1, 
 TUPLE_SET * Tuple2,
 float       intron_to_intron,
 float       negative_inf
)
{
  // Quick catch 1: Introns
  if (Tuple1->MarksIntron) {
    if (Tuple2->MarksIntron) {
      return basescore+intron_to_intron;
    } else {
      return negative_inf;
    }
  } else if (Tuple2->MarksIntron) {
    return negative_inf;
  }

  // Otherwise, use our scoring matrix
  float score = 0.0;
  int i,j,x,y;
  for (i=0; i<Tuple1->NumTuples; i++) {
    for (j=0; j<Tuple2->NumTuples; j++) {
      x = 20*ConvertToIndex(Tuple1->TupleChars[i]);
      y = ConvertToIndex(Tuple2->TupleChars[j]);
      if (x < 0 || y < 0) continue; // Something's weird, but it'll be weird across the board
      score += (BLOSUM62[x+y] * Tuple1->TupleRatios[i] * Tuple2->TupleRatios[j]);
    }
  }

  return basescore+score;

}




/* FUNCTION: NWTraceback
 *
 *    ABOUT: Once we've completed our Needleman-Wunsch dynamic
 *           programming matrix, we identify an optimal path
 *           back through the table.  We compute the coordinates
 *           of each step along that path and return the length
 *           of the path.
 *
 */
int
NWTraceback
(
 float     ** Match,
 float     ** Insert,
 float     ** Delete,
 TUPLE_SET ** TS1,
 int          TS1Length,
 TUPLE_SET ** TS2,
 int          TS2Length,
 int        * OptimalX,
 int        * OptimalY,
 float        intron_gap_base,
 float        intron_gap_max,
 float        intron_gap_cont,
 float        intron_gap_mult,
 float        gap_start,
 float        gap_continue,
 float        negative_inf
)
{
  int i,j,k;
  i = TS1Length;
  j = TS2Length;
  k = 0;

  int state = 1;

  while (i > 1 && j > 1) {

    if (state == 1) {
      
      if (Match[i-1][j-1] > Max(Insert[i-1][j-1],Delete[i-1][j-1])) {
	state = 1;
      } else if (Insert[i-1][j-1] > Delete[i-1][j-1]) {
	state = 2;
      } else {
	state = 3;
      }      
      i--;
      j--;

    } else if (state == 2) {
    
      if (Insert[i][j] == CalcGap(Insert[i][j-1],TS1[i-1],TS2[j-2],0,1,
				  intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
				  gap_start,gap_continue,negative_inf)) {
	state = 2;
      } else {
	state = 1;
      }
      j--;
      
    } else {

      if (Delete[i][j] == CalcGap(Delete[i-1][j],TS1[i-2],TS2[j-1],0,0,
				  intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
				  gap_start,gap_continue,negative_inf)) {
	state = 3;
      } else {
	state = 1;
      }
      i--;

    }
     
    OptimalX[k] = i;
    OptimalY[k] = j;
    k++;

  }
  
  while (i > 1) {
    i--;
    OptimalX[k] = i;
    OptimalY[k] = j;
    k++;
  }

  while (j > 1) {
    j--;
    OptimalX[k] = i;
    OptimalY[k] = j;
    k++;
  }

  // Because (1,1) has to be a match, we just build this in to our
  // MSA construction
  if (OptimalX[k] != 1 || OptimalY[k] != 1) {
    OptimalX[k] = 1;
    OptimalY[k] = 1;
    k++;
  }
  /*
  printf("\n");
  for (i = k-1; i >= 0; i--) {
    printf("%d, %d\n",OptimalX[i],OptimalY[i]);
  }
  printf("\n\n");
  */
  // By now, k is the length of the path we're tracing back.
  return k;

}





/* FUNCTION: IntronAdjustment
 *
 *    ABOUT: Try to make smart corrections to our final MSA
 *           so that we have the thing on the right, instead
 *           of the thing on the left:
 *
 *                          AB*-C     AB*C
 *                          AB*-C --> AB*C
 *                          A-*BC     AB*C
 *
 *           We do a similar final check in FinalMSA.pl, but
 *           it's worth making these adjustments here, too.
 *
 */
int
IntronAdjustment
(
 char ** MSA,
 int     numSeqs,
 int     seqLength
)
{
  int i,j,k;

  int  intron;
  char compare;
  int  yay;
  int  nay;
  int  Lgaps,Rgaps;

  int alreadySwapped[numSeqs];

  int placer = 1;
  int runner = 1;
  while (runner < seqLength-2) {

    // Check if we're in an intron 'danger zone'
    intron = 0;
    for (i=0; i<numSeqs; i++) {
      if (MSA[i][runner+1] == '*') {
	runner++;
	intron = 1;
	break;
      }
    }
    
    // Might we get to do some movin' around? Please?
    if (intron) {

      // For each position, check if it has a gap on the left
      // side.  If so and if it makes sense to swap with the gap
      // then we swap.
      for (i=0; i<numSeqs; i++) {

	alreadySwapped[i] = 0;

	if (MSA[i][runner-1] == '-') {

	  compare = MSA[i][runner+1];
	  yay     = 0;
	  nay     = 0;

	  for (j=0; j<numSeqs; j++) {

	    if (i == j)                      continue;
	    if (MSA[j][runner-1] == compare) yay++;
	    if (MSA[j][runner+1] == compare) nay++;

	  }

	  if (yay > nay || nay == 0) {
	    MSA[i][runner+1]  = '-';
	    MSA[i][runner-1]  = compare;
	    alreadySwapped[i] = 1;
	  }
	}
      }

      // For each position, check if it has a gap on the right
      // side.  If so and if it makes sense to swap with the gap
      // then we swap.  Note that we don't do 'no complaints'
      // swaps from the left to the right.
      for (i=0; i<numSeqs; i++) {

	if (!alreadySwapped[i] && MSA[i][runner+1] == '-') {

	  compare  = MSA[i][runner-1];
	  yay      = 0;
	  nay      = 0;

	  for (j=0; j<numSeqs; j++) {

	    if (i == j)                      continue;
	    if (MSA[j][runner+1] == compare) yay++;
	    if (MSA[j][runner-1] == compare) nay++;
	    
	  }

	  if (yay > nay) {
	    MSA[i][runner-1] = '-';
	    MSA[i][runner+1] = compare;
	  }
	}
      }
 
      // Now we check the columns surrounding the intron to make
      // sure there aren't any 'all-gap' comlumns.  If there are,
      // we cover them up.
	
      Lgaps = 1;
      for (i=0; i<numSeqs; i++) {
	if (MSA[i][runner-1] != '-') {
	  Lgaps = 0;
	  break;
	}
      }
      
      Rgaps = 1;
      for (i=0; i<numSeqs; i++) {
	if (MSA[i][runner+1] != '-') {
	  Rgaps = 0;
	  break;
	}
      }
      
      if (Lgaps) {
	for (i=0; i<numSeqs; i++) {
	  MSA[i][placer] = MSA[i][runner];
	}
	placer++;
      } else {
	for (i=0; i<numSeqs; i++) {
	  MSA[i][placer]   = MSA[i][runner-1];
	  MSA[i][placer+1] = MSA[i][runner];
	}
	placer += 2;
      }
      
      runner++;
      if (!Rgaps) {
	for (i=0; i<numSeqs; i++) {
	  MSA[i][placer] = MSA[i][runner];
	}
	placer++;
      }
      runner++;
      
      // We SHOULD be ready to rock and roll some more!
      
    } else {

      // Copy and increment
      for (i=0; i<numSeqs; i++)
	MSA[i][placer] = MSA[i][runner];
      placer++;
      runner++;
      
    }

  }

  // Note that it's possible for runner to extend past segLength-2,
  // such that we might have already filled in the final column
  if (runner == seqLength-2) {
    for (i=0; i<numSeqs; i++) {
      MSA[i][placer]   = MSA[i][runner];
      MSA[i][placer+1] = MSA[i][runner+1];
    }
    placer += 2;
  } else if (runner == seqLength-1) {
    for (i=0; i<numSeqs; i++)
      MSA[i][placer] = MSA[i][runner];
    placer++;
  }

  // Now "placer" is our MSA length.
  seqLength = placer;

  // We do a second scan looking for "jigsaw" situations (possibly
  // introduced during the last phase).  For example (via OBSCN):
  //
  //               S*-SV         S*SV
  //               S*-SV   -->   S*SV
  //               S*S-V         S*SV
  //
  // Because of the work done during the first scan, these should
  // be restricted to one side of the splice site or the other.
  // We also assume they'll be single gaps.

  placer = 1;
  runner = 1;
  while (runner < seqLength-1) {
    
    // This "problem" should be isolated to places where we see splice
    // sites, so we're intron hunting again.  Also, because we shift left
    // above, we just need to worry about what's right of the splice site
    intron = 0;
    for (i=0; i<numSeqs; i++) {
      if (MSA[i][runner] == '*') {
	intron = 1;
	break;
      }
    }

    // Record and move forward
    for (i=0; i<numSeqs; i++)
      MSA[i][placer] = MSA[i][runner];
    placer++;
    runner++;
    
    // Did we see a splice site?  Are we reasonable to think this could
    // be a "jigsaw" area?
    if (intron) {
     
      // It'd be surprising if this sort of issue occurred for more than
      // one pair of columns, but we might as well cover all our bases
      while (intron && runner < seqLength-2 
	     && (MSA[0][runner] == '-' || MSA[0][runner+1] == '-')
	     && MSA[0][runner] != MSA[0][runner+1]) {

	// To play it safe with post-hoc alignment changing, we'll only
	// perform this switch if there's unanimous consensus on the 
	// "puzzle piece"
	if (MSA[0][runner] == '-') { 
	  alreadySwapped[0] = 1;
	  compare = MSA[0][runner+1];
	} else { 
	  alreadySwapped[0] = 0;
	  compare = MSA[0][runner];   
	}
	
	// Now we just make sure this fits the profile (well, really two profiles)
	yay = 1; // Are we happy?
	for (i=1; i<numSeqs; i++) {
	  
	  // We didn't like the looks of this
	  if ((MSA[i][runner] != '-' && MSA[i][runner+1] != '-')
	      || MSA[i][runner] == MSA[i][runner+1]) {
	    yay = 0;
	    break;
	  }
	  
	  // Are we on the path to unanimous agreement?
	  if (MSA[i][runner] == '-') {
	    alreadySwapped[i] = 1;
	    if (MSA[i][runner+1] != compare) {
	      yay = 0;
	      break;
	    }
	  } else {
	    alreadySwapped[i] = 0;
	    if (MSA[i][runner] != compare) {
	      yay = 0;
	      break;
	    }
	  }
	  
	}

	// Ready to get weird with things?
	if (yay) {
	  
	  // Oh, things are getting weird alright
	  for (i=0; i<numSeqs; i++) 
	    MSA[i][placer] = MSA[i][runner+alreadySwapped[i]];
	  placer++;
	  runner += 2;
	  
	} else {

	  // One of these days you're going to get weird.
	  intron = 0;

	}
	
      }

    }

  }

  // Copy over the final position (should be a splice site) and bail
  for (i=0; i<numSeqs; i++) 
    MSA[i][placer] = MSA[i][seqLength-1];
  placer++;

  return placer;

}





/* FUNCTION: MSNeedlemanWunsch
 *
 *    ABOUT: Perform multiple-sequence Needleman-Wunsch dynamic
 *           programming algorithm on our profiles.
 *
 */
int
MSNeedlemanWunsch
(
 TUPLE_SET ** TS1,
 char      ** MSA1,
 char      ** MSA1Names,
 int          MSA1Size,
 int          MSA1Length,
 TUPLE_SET ** TS2,
 char      ** MSA2,
 char      ** MSA2Names,
 int          MSA2Size,
 int          MSA2Length,
 float        intron_gap_base,
 float        intron_gap_max,
 float        intron_gap_cont,
 float        intron_gap_mult,
 float        gap_start,
 float        gap_continue,
 float        intron_to_intron,
 float        negative_inf,
 int          debug
)
{
  int i,j,k;

  float ** Match;
  if ((Match = malloc((MSA1Length+1)*sizeof(float *))) == NULL) {
    printf("  ERROR:  Failed to initialize Match\n");
    return 1;
  }
  for (i=0; i<=MSA1Length; i++) {
    if ((Match[i] = malloc((MSA2Length+1)*sizeof(float))) == NULL) {
      printf("  ERROR:  Failed to complete initialization of Match\n");
      return 1;
    }
  }

  float ** Insert;
  if ((Insert = malloc((MSA1Length+1)*sizeof(float *))) == NULL) {
    printf("  ERROR:  Failed to initialize Insert\n");
    return 1;
  }
  for (i=0; i<=MSA1Length; i++) {
    if ((Insert[i] = malloc((MSA2Length+1)*sizeof(float))) == NULL) {
      printf("  ERROR:  Failed to complete initialization of Insert\n");
      return 1;
    }
  }

  float ** Delete;
  if ((Delete = malloc((MSA1Length+1)*sizeof(float *))) == NULL) {
    printf("  ERROR:  Failed to initialize Delete\n");
    return 1;
  }
  for (i=0; i<=MSA1Length; i++) {
    if ((Delete[i] = malloc((MSA2Length+1)*sizeof(float))) == NULL) {
      printf("  ERROR:  Failed to complete initialization of Delete\n");
      return 1;
    }
  }

  // Setting the corners
  Match[0][0]  = 0.0;
  Insert[0][0] = 0.0;
  Delete[0][0] = 0.0;
  Match[1][0]  = negative_inf;
  Insert[1][0] = negative_inf;
  Delete[1][0] = gap_start;
  Match[0][1]  = negative_inf;
  Insert[0][1] = gap_start;
  Delete[0][1] = negative_inf;

  // Setting the margins
  for (i=2; i<=MSA1Length; i++) {
    Match[i][0]  = negative_inf;
    Insert[i][0] = negative_inf;
    Delete[i][0] = Delete[i-1][0]+gap_continue;
  }
  for (j=2; j<=MSA2Length; j++) {
    Match[0][j]  = negative_inf;
    Insert[0][j] = Insert[0][j-1]+gap_continue;
    Delete[0][j] = negative_inf;
  }

  // Last very special cases: (1,1), (2,1), and (1,2)
  Match[1][1]  = intron_to_intron;
  Insert[1][1] = negative_inf;
  Insert[1][1] = negative_inf;

  for (i=2; i<=MSA1Length; i++) {

    if (TS1[i-1]->MarksIntron) Match[i][1] = intron_to_intron;
    else                       Match[i][1] = negative_inf;
    Insert[i][1] = Max(Match[i][0]+gap_start,Insert[i][0]+gap_continue);
    Delete[i][1] = Max(Match[i-1][1]+gap_start,Delete[i-1][1]+gap_continue);

  }

  for (j=2; j<=MSA2Length; j++) {

    if (TS2[j-1]->MarksIntron) Match[1][j] = intron_to_intron;
    else                       Match[1][j] = negative_inf;
    Delete[1][j] = Max(Match[0][j]+gap_start,Delete[0][j]+gap_continue);
    Insert[1][j] = Max(Match[1][j-1]+gap_start,Insert[1][j-1]+gap_continue);

  }

  // The std. recurrence.  Note that because we're starting at 1,1 in the
  // dynamic programming table, everything that references the profiles 
  // (TS1,TS2) requires subtracting an additional 1.
  for (i=2; i<=MSA1Length; i++) {
    for (j=2; j<=MSA2Length; j++) {

      // Match State
      Match[i][j] = Max(Insert[i-1][j-1],Delete[i-1][j-1]);
      Match[i][j] = Max(Match[i][j],Match[i-1][j-1]);
      Match[i][j] = CalcMatch(Match[i][j],TS1[i-1],TS2[j-1],intron_to_intron,negative_inf);
      
      // Note that we're worried about intron gaps in different ways when
      // we're aligning to intron-containing MSAs and when we're aligning
      // a raw sequence to an intron-containing MSA, so we need to check this.
      // >> In brief, when we're 
      // Insert State
      Insert[i][j] = Max(CalcGap(Insert[i][j-1],TS1[i-1],TS2[j-2],0,1,
				 intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
				 gap_start,gap_continue,negative_inf),
			 CalcGap(Match[i][j-1],TS1[i-1],TS2[j-2],1,1,
				 intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
				 gap_start,gap_continue,negative_inf));
      
      // Delete Stat
      Delete[i][j] = Max(CalcGap(Delete[i-1][j],TS1[i-2],TS2[j-1],0,0,
				 intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
				 gap_start,gap_continue,negative_inf),
			 CalcGap(Match[i-1][j],TS1[i-2],TS2[j-1],1,0,
				 intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
				 gap_start,gap_continue,negative_inf));
      
    }

  }

  /*
   * Debugging - Not recommended for large tests :P
   *
  printf("\nMATCH:\n");
  for (j=0; j<=MSA2Length; j++) {
    for (i=0; i<=MSA1Length; i++) {
      printf("%.0f\t",Match[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");
  printf("\nDELETE:\n");
  for (j=0; j<=MSA2Length; j++) {
    for (i=0; i<=MSA1Length; i++) {
      printf("%.0f\t",Delete[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");
  printf("\nINSERT:\n");
  for (j=0; j<=MSA2Length; j++) {
    for (i=0; i<=MSA1Length; i++) {
      printf("%.0f\t",Insert[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");
  */

  // Prep the arrays that will store the coordinates resulting from traceback
  int * OptimalX;
  int * OptimalY;
  if ((OptimalX = malloc((MSA1Length+MSA2Length)*sizeof(int))) == NULL) {
    printf("  ERROR:  Failed to allocate 'OptimalX'\n");
    return 1;
  }
  if ((OptimalY = malloc((MSA1Length+MSA2Length)*sizeof(int))) == NULL) {
    printf("  ERROR:  Failed to allocate 'OptimalY'\n");
    return 1;
  }

  
  int OutputSeqLength;
  OutputSeqLength = NWTraceback(Match,Insert,Delete,TS1,MSA1Length,
				TS2,MSA2Length,OptimalX,OptimalY,
				intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
				gap_start,gap_continue,negative_inf);


  char ** FinalMSA;
  if ((FinalMSA = malloc((MSA1Size+MSA2Size) * sizeof(char *))) == NULL) {
    printf("  ERROR:  Failed to initialize 'tempMSA'\n");
    return 1;
  }
  for (i=0; i<MSA1Size+MSA2Size; i++) {
    if ((FinalMSA[i] = malloc(OutputSeqLength * sizeof(char))) == NULL) {
      printf("  ERROR:  Failed to finish initializing 'tempMSA'\n");
      return 1;
    }
  }

  
  // Prepare to generate the final MSA
  for (k=0; k<MSA1Size; k++) {
    
    i = 0;
    j = MSA1Length-1;
    while (i < OutputSeqLength-1 && j > 0) {

      if (i && OptimalX[i] == OptimalX[i-1]) {
	FinalMSA[k][(OutputSeqLength-1)-i] = '-';
      } else {
	FinalMSA[k][(OutputSeqLength-1)-i] = MSA1[k][j];
	j--;
      }
      i++;
      
    }

    while (i < OutputSeqLength-1) {
      FinalMSA[k][(OutputSeqLength-1)-i] = '-';
      i++;
    }
    FinalMSA[k][0] = '*';

  }

  
  // Prepare for the final MSA (part 2)
  for (k=MSA1Size; k<MSA1Size+MSA2Size; k++) {
    
    i = 0;
    j = MSA2Length-1;
    while (i < OutputSeqLength-1 && j > 0) {
      
      if (i && OptimalY[i] == OptimalY[i-1]) {
	FinalMSA[k][(OutputSeqLength-1)-i] = '-';
      } else {
	FinalMSA[k][(OutputSeqLength-1)-i] = MSA2[k-MSA1Size][j];
	j--;
      }
      i++;
      
    }
    
    while (i < OutputSeqLength-1) {
      FinalMSA[k][(OutputSeqLength-1)-i] = '-';
      i++;
    }
    FinalMSA[k][0] = '*';

  }

  // Nothing more to see here, folks
  free(OptimalX);
  free(OptimalY);
  for (i=0; i<=MSA1Length; i++) {
    free(Match[i]);
    free(Insert[i]);
    free(Delete[i]);
  }
  free(Match);
  free(Insert);
  free(Delete);

  
  // We run a quick adjustment around the intronic regions
  OutputSeqLength = IntronAdjustment(FinalMSA,MSA1Size+MSA2Size,OutputSeqLength);

  // Print out the final MSA!
  for (k=0; k<MSA1Size; k++) {
    printf("%s",MSA1Names[k]);
    i = 0;
    j = 0;
    while (i < OutputSeqLength) {
      printf("%c",FinalMSA[k][i]);
      i++;
      j++;
      if (j == 50) {
	printf("\n");
	j = 0;
      }
    }
    if (j) printf("\n");
    printf("\n");
  }
  
  for (k=MSA1Size; k<MSA1Size+MSA2Size; k++) {
    printf("%s",MSA2Names[k-MSA1Size]);
    i = 0;
    j = 0;
    while (i < OutputSeqLength) {
      printf("%c",FinalMSA[k][i]);
      i++;
      j++;
      if (j == 50) {
	printf("\n");
	j = 0;
      }
    }
    if (j) printf("\n");
    printf("\n");
  }

  for (i=0; i<MSA1Size+MSA2Size; i++)
    free(FinalMSA[i]);
  free(FinalMSA);

  return 0;
}




/* FUNCTION: PrintUsageMsg
 *
 *    ABOUT: This function gives the user a nudge in the right direction.
 *
 */
int PrintUsageMsg ()
{
  printf("\n");
  printf("\tUSAGE:  ./MultiSeqNW  <msa1>  <msa1 size>  <msa2>  <msa2 size>  [-debug]\n");
  printf("\n");
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
  if (argc < 5) return PrintUsageMsg();
  
  int i,j,k;
  
  // gap penalties (command-line customizable) <- changing to a struct might be cleaner
  float intron_gap_base  = -5.0;
  float intron_gap_max   = -200.0;
  float intron_gap_cont  = 0.0;
  float intron_gap_mult  = 2.0;
  float gap_start        = -11.0; // shouldn't change
  float gap_continue     = -1.0;  // shouldn't change
  float intron_to_intron = 0.0;   // eh...
  float negative_inf     = 0.0-INFINITY; // SHOULDN'T CHANGE

  int debug = 0;
  i = 5;
  while (i < argc) {
    if (!strcmp(argv[i],"-debug")) {
      debug = 1;
    } else if (!strcmp(argv[i],"-igBase")) {
      i++;
      intron_gap_base = atof(argv[i]);
    } else if (!strcmp(argv[i],"-igMax")) {
      i++;
      intron_gap_max = atof(argv[i]);
    } else if (!strcmp(argv[i],"-igMult")) {
      i++;
      intron_gap_mult = atof(argv[i]);
      if (intron_gap_mult < 0.0) {
	printf("\n  Cannot have a negative multiplier\n\n");
	intron_gap_mult = 2.0;
      }
    } else if (!strcmp(argv[i],"-igCont")) {
      i++;
      intron_gap_cont = atof(argv[i]);
    } else if (!strcmp(argv[i],"-itoi")) {
      i++;
      intron_to_intron = atof(argv[i]);
    } else {
      printf("  Unrecognized option '%s' ignored\n",argv[i]);
    }
    i++;
  }
  
  // How many sequences are in each of the MSAs?
  int MSA1Size = atoi(argv[2]);
  int MSA2Size = atoi(argv[4]);

  // Allocate for names (assuming no file has a name > 256 chars)
  char ** MSA1Names;
  if ((MSA1Names = malloc(MSA1Size*sizeof(char *))) == NULL) {
    printf("  ERROR:  Failed to allocate for 'MSA1Names'!\n");
    return 1;
  }
  for (i=0; i<MSA1Size; i++) {
    if ((MSA1Names[i] = malloc(256*sizeof(char))) == NULL) {
      printf("  ERROR:  Failed to allocate for 'MSA1Names'\n");
      return 1;
    }
    for (j=0; j<256; j++)
      MSA1Names[i][j] = 0;
  }
  char ** MSA2Names;
  if ((MSA2Names = malloc(MSA2Size*sizeof(char *))) == NULL) {
    printf("  ERROR:  Failed to allocate for 'MSA2Names'!\n");
    return 1;
  }
  for (i=0; i<MSA2Size; i++) {
    if ((MSA2Names[i] = malloc(256*sizeof(char))) == NULL) {
      printf("  ERROR:  Failed to allocate for 'MSA2Names'\n");
      return 1;
    }
    for (j=0; j<256; j++)
      MSA2Names[i][j] = 0;
  }

  // Open the protein FASTA file
  FILE * MSA1File;
  MSA1File = fopen(argv[1],"r");
  if (MSA1File == NULL) {
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
  int name_pos = 0;
  while ((MSA1Names[0][name_pos] = fgetc(MSA1File)) >= 32) {
    name_pos++;
  }
  
  // Prep for the first MSA
  int    MSA1Length;
  char * FirstSeq;
  int    capacity = 512;

  if ((FirstSeq = malloc(capacity*sizeof(char))) == NULL) {
    printf("  ERROR: Could not allocate string 'Seq1'!\n");
    return 1;
  } else {
    for (i=0; i<capacity; i++) 
      FirstSeq[i] = 0;
  }

  // Rip the protein file, keeping an eye on all the junk you need to.
  char * tempSeq = NULL;
  MSA1Length = 0;
  while (!feof(MSA1File)) {
    fscanf(MSA1File,"%s",next_line);
    if (next_line[0] == '>') break;
    i = 0;
    while (next_line[i] > 32) {
      FirstSeq[MSA1Length] = next_line[i];
      next_line[i] = 0;
      MSA1Length++;
      if (MSA1Length == capacity) {
	if (tempSeq != NULL) {
	  if ((tempSeq = realloc(tempSeq,capacity*sizeof(char))) == NULL) {
	    printf("  ERROR:  Failed to reallocate string 'tempSeq'\n");
	    return 1;
	  }
	} else {
	  if ((tempSeq = malloc(capacity*sizeof(char))) == NULL) {
	    printf("  ERROR:  Failed to allocate string 'tempSeq'\n");
	    return 1;
	  }
	}
	memcpy(tempSeq,FirstSeq,capacity*sizeof(char));
	if ((FirstSeq = realloc(FirstSeq,capacity*2*sizeof(char))) == NULL) {
	  printf("  ERROR: Failed to reallocate string 'Seq1'!\n");
	  return 1;
	}
	memcpy(FirstSeq,tempSeq,capacity*sizeof(char));
	for (j=capacity; j<capacity*2; j++) 
	  FirstSeq[j] = 0;
	capacity *= 2;
      }
      i++;
    }
  }


  // It's possible that a sequence might be passed in without any splice
  // site markers, which would cause problems during traceback, so we just
  // do a quick check
  int needs_splices = 0;
  if (FirstSeq[0] != '*' && FirstSeq[0] != '-') {
    
    // First off, do we have the room to add the flanking splice site characters?
    if (MSA1Length+2 >= capacity) {
      if (tempSeq != NULL) {
	if ((tempSeq = realloc(tempSeq,capacity*sizeof(char))) == NULL) {
	  printf("  ERROR:  Failed to reallocate string 'tempSeq'\n");
	  return 1;
	}
      } else {
	if ((tempSeq = malloc(capacity*sizeof(char))) == NULL) {
	  printf("  ERROR:  Failed to allocate string 'tempSeq'\n");
	  return 1;
	}
      }
      memcpy(tempSeq,FirstSeq,capacity*sizeof(char));
      if ((FirstSeq = realloc(FirstSeq,capacity*2*sizeof(char))) == NULL) {
	printf("  ERROR: Failed to reallocate string 'Seq1'!\n");
	return 1;
      }
      memcpy(FirstSeq,tempSeq,capacity*sizeof(char));
      for (j=capacity; j<capacity*2; j++) 
	FirstSeq[j] = 0;
      capacity *= 2;
    }
    
    // Now we can actual get our sequence set up
    for (i=MSA1Length; i>0; i--)
      FirstSeq[i] = FirstSeq[i-1];
    FirstSeq[0]            = '*';
    FirstSeq[MSA1Length+1] = '*';
    MSA1Length += 2;

    // Record that we're going to be adding splice sites to this batch
    needs_splices = 1;

  }


  // Now we know the length of the first batch of seqs, so our job becomes a bit easier
  // for the following sequences in the MSA.
  char ** MSA1;
  if ((MSA1 = malloc(MSA1Size*sizeof(char *))) == NULL) {
    printf("  ERROR: Failed to allocate 'MSA1'!\n");
    return 1;
  }
  for (k=0; k<MSA1Size; k++) {
    if ((MSA1[k] = malloc(MSA1Length*sizeof(char))) == NULL) {
      printf("  ERROR: Failed to allocate 'MSA1'\n");
      return 1;
    }
  }

  // Copying over the first sequence in the MSA
  for (i=0; i<MSA1Length; i++)
    MSA1[0][i] = FirstSeq[i];


  // Hold onto the string containing the sequence name
  char * nameholder;
  if ((nameholder = malloc(128*sizeof(char))) == NULL) {
    printf("  ERROR:  Failed to allocate string 'nameholder'\n");
    return 1;
  }
  

  // Grabbing the remaining MSAs
  int pos;
  for (k=1; k<MSA1Size; k++) {
  
    // Read in the name of the next sequence    
    if (next_line[0] != '>') {
      fscanf(MSA1File,"%s",nameholder);
      while (nameholder[0] != '>') { 
	fscanf(MSA1File,"%s",nameholder); 
      }
    } else {
      i = 0;
      while (next_line[i] >= 32) {
	nameholder[i] = next_line[i];
	i++;
      }
      nameholder[i] = 0;
    }

    // Record the sequence name
    name_pos = 0;
    while (nameholder[name_pos] >= 32) {
      MSA1Names[k][name_pos] = nameholder[name_pos];
      name_pos++;
    }    
    
    // Make sure it ends in a linebreak
    if (MSA1Names[k][name_pos] != '\n') {
      MSA1Names[k][name_pos] = '\n';
      name_pos++;
    }
    MSA1Names[k][name_pos] = 0;
    
    // Rip the protein file, keeping an eye on all the junk you need to.
    pos = 0;
    
    // Do we need to add splice sites?
    if (needs_splices) {
      MSA1[k][pos] = '*';
      pos++;
    }
    
    while (!feof(MSA1File)) {
      fscanf(MSA1File,"%s",next_line);
      if (next_line[0] == '>') break;
      i = 0;
      while (next_line[i] > 32) {
	MSA1[k][pos] = next_line[i];
	next_line[i] = 0;
	pos++;
	i++;
      }
    }

    if (needs_splices) {
      MSA1[k][pos] = '*';
      pos++;
    }

  }
  fclose(MSA1File);
  
  // Prep for the second MSA
  int MSA2Length= 0;

  // Open up the DNA/RNA string
  FILE * MSA2File;
  MSA2File = fopen(argv[3],"r");
  if (MSA2File == NULL) {
    printf("  ERROR:  Could not open file '%s'\n",argv[2]);
    fclose(MSA2File);
    return 1;
  }
    
  // Because we don't want to risk weird formatting breaking everything,
  // we eat the header line char by char.
  name_pos = 0;
  while ((MSA2Names[0][name_pos] = fgetc(MSA2File)) >= 32) {
    name_pos++;
  }

  // Re-prep the FirstSeq string
  for (i=0; i<capacity; i++)
    FirstSeq[i] = 0;
  
  // Rip that DNA file!
  while (!feof(MSA2File)) {    
    fscanf(MSA2File,"%s",next_line);
    if (next_line[0] == '>') break;
    i = 0;
    while (next_line[i] > 32) {
      FirstSeq[MSA2Length] = next_line[i];
      next_line[i] = 0;
      MSA2Length++;
      i++;
      if (MSA2Length == capacity) {	
	if (tempSeq != NULL) {
	  if ((tempSeq = realloc(tempSeq,capacity*sizeof(char))) == NULL) {
	    printf("  ERROR:  Failed to reallocate string 'tempSeq'\n");
	    return 1;
	  }
	} else {
	  if ((tempSeq = malloc(capacity*sizeof(char))) == NULL) {
	    printf("  ERROR:  Could not allocate string 'tempSeq'\n");
	    return 1;
	  }
	}	
	memcpy(tempSeq,FirstSeq,capacity*sizeof(char));
	if ((FirstSeq = realloc(FirstSeq,capacity*2*sizeof(char))) == NULL) {
	  printf("  ERROR:  Failed to reallocate string 'NuclString'!\n");
	  return 1;
	}
	memcpy(FirstSeq,tempSeq,capacity*sizeof(char));	
	for (j=capacity; j<capacity*2; j++)
	  FirstSeq[j] = 0;
	capacity *= 2;
      }
    }
  }
  if (tempSeq) free(tempSeq);


  // It's possible that this sequence might be passed in without splice sites, too
  needs_splices = 0;
  if (FirstSeq[0] != '*' && FirstSeq[0] != '-') {
    
    // First off, do we have the room to add the flanking splice site characters?
    if (MSA2Length+2 >= capacity) {
      if (tempSeq != NULL) {
	if ((tempSeq = realloc(tempSeq,capacity*sizeof(char))) == NULL) {
	  printf("  ERROR:  Failed to reallocate string 'tempSeq'\n");
	  return 1;
	}
      } else {
	if ((tempSeq = malloc(capacity*sizeof(char))) == NULL) {
	  printf("  ERROR:  Failed to allocate string 'tempSeq'\n");
	  return 1;
	}
      }
      memcpy(tempSeq,FirstSeq,capacity*sizeof(char));
      if ((FirstSeq = realloc(FirstSeq,capacity*2*sizeof(char))) == NULL) {
	printf("  ERROR: Failed to reallocate string 'Seq1'!\n");
	return 1;
      }
      memcpy(FirstSeq,tempSeq,capacity*sizeof(char));
      for (j=capacity; j<capacity*2; j++) 
	FirstSeq[j] = 0;
      capacity *= 2;
    }
    
    // Now we can actual get our sequence set up
    for (i=MSA2Length; i>0; i--)
      FirstSeq[i] = FirstSeq[i-1];
    FirstSeq[0]            = '*';
    FirstSeq[MSA2Length+1] = '*';
    MSA2Length += 2;

    // Record that we're going to be adding splice sites to this batch
    needs_splices = 1;

  }

  // Now we know the length of the second batch of seqs, so etc.
  char ** MSA2;
  if ((MSA2 = malloc(MSA2Size*sizeof(char *))) == NULL) {
    printf("  ERROR: Failed to allocate 'MSA2'!\n");
    return 1;
  }
  for (k=0; k<MSA2Size; k++) {
    if ((MSA2[k] = malloc(MSA2Length*sizeof(char))) == NULL) {
      printf("  ERROR: Failed to allocate 'MSA2'\n");
      return 1;
    }
  }

  // Copying over the first sequence in the MSA and freeing FirstSeq
  for (i=0; i<MSA2Length; i++)
    MSA2[0][i] = FirstSeq[i];
  free(FirstSeq);

  // Grabbing the remaining MSAs
  for (k=1; k<MSA2Size; k++) {
  
    // Read in the name of the next sequence
    if (next_line[0] != '>') {
      fscanf(MSA2File,"%s",nameholder);
      while (nameholder[0] != '>') { 
	fscanf(MSA2File,"%s",nameholder); 
      }
    } else {
      i = 0;
      while (next_line[i] >= 32) {
	nameholder[i] = next_line[i];
	i++;
      }
      nameholder[i] = 0;
    }

    // Record the sequence name
    name_pos = 0;
    while (nameholder[name_pos] >= 32) {
      MSA2Names[k][name_pos] = nameholder[name_pos];
      name_pos++;
    }

    // Make sure it ends in a linebreak
    if (MSA2Names[k][name_pos] != '\n') {
      MSA2Names[k][name_pos] = '\n';
      name_pos++;
    }
    MSA2Names[k][name_pos] = 0;

    // Do we need to add splice sites?
    pos = 0;
    if (needs_splices) {
      MSA2[k][pos] = '*';
      pos++;
    }
    
    // Rip the protein file, keeping an eye on all the junk you need to.
    while (!feof(MSA2File)) {
      fscanf(MSA2File,"%s",next_line);
      if (next_line[0] == '>') break;
      i = 0;
      while (next_line[i] > 32) {
	MSA2[k][pos] = next_line[i];
	next_line[i] = 0;
	pos++;
	i++;
      }
    }

    if (needs_splices) {
      MSA2[k][pos] = '*';
      pos++;
    }    

  }
  fclose(MSA2File);

  // No longer need next_line
  free(next_line);


  /*
   *
   * Now we condense our sequences down into tuples
   *
   */

  
  // The starting point for our tuples
  char * tuple_base;
  if ((tuple_base = malloc(Max(MSA1Size,MSA2Size) * sizeof(char))) == NULL) {
    printf("  ERROR:  Failed to initialize 'tuple_base'\n");
    return 1;
  }


  // MSA1 -> TS1
  TUPLE_SET ** TS1;
  if ((TS1 = malloc(MSA1Length * sizeof(TUPLE_SET *))) == NULL) {
    printf("  ERROR:  Failed to initialize 'TS1'!\n");
    return 1;
  }
  for (j=0; j<MSA1Length; j++) {
    for (i=0; i<MSA1Size; i++)
      tuple_base[i] = MSA1[i][j];
    TS1[j] = ConvertToTupleSet(tuple_base,MSA1Size);
  }

  // Scan along, marking proximity to nearest intron
  for (j=1; j<MSA1Length; j++) {
    if (TS1[j]->MarksIntron) continue;
    TS1[j]->NearestIntron = TS1[j-1]->NearestIntron+1;
  }
  for (j=MSA1Length-2; j>0; j--) {
    if (TS1[j]->MarksIntron) continue;
    TS1[j]->NearestIntron = Min(TS1[j]->NearestIntron,TS1[j+1]->NearestIntron+1);
  }

  // MSA2 -> TS2
  TUPLE_SET ** TS2;
  if ((TS2 = malloc(MSA2Length * sizeof(TUPLE_SET *))) == NULL) {
    printf("  ERROR:  Failed to initialize 'TS1'!\n");
    return 1;
  }
  for (j=0; j<MSA2Length; j++) {
    for (i=0; i<MSA2Size; i++)
      tuple_base[i] = MSA2[i][j];
    TS2[j] = ConvertToTupleSet(tuple_base,MSA2Size);    
  }

  // Scan along, marking proximity to nearest intron
  for (j=1; j<MSA2Length; j++) {
    if (TS2[j]->MarksIntron) continue;
    TS2[j]->NearestIntron = TS2[j-1]->NearestIntron+1;
  }
  for (j=MSA2Length-2; j>0; j--) {
    if (TS2[j]->MarksIntron) continue;
    TS2[j]->NearestIntron = Min(TS2[j]->NearestIntron,TS2[j+1]->NearestIntron+1);
  }

  // You have outlived your utility, tuple_base.  Prepare to die.
  free(tuple_base);

  // Get to work!
  i = MSNeedlemanWunsch(TS1,MSA1,MSA1Names,MSA1Size,MSA1Length,
                        TS2,MSA2,MSA2Names,MSA2Size,MSA2Length,
			intron_gap_base,intron_gap_max,intron_gap_cont,intron_gap_mult,
			gap_start,gap_continue,intron_to_intron,negative_inf
			,debug);

  
  // YOU RASCALLY TUPLE_SETS! GET OFF MY LAWN!
  for (i=0; i<MSA1Size; i++)
    DestroyTupleSet(TS1[i]);
  free(TS1);
  for (i=0; i<MSA2Size; i++)
    DestroyTupleSet(TS2[i]);
  free(TS2);

  
  // Don't even get me started on you, MSA1 and MSA2
  for (i=0; i<MSA1Size; i++)
    free(MSA1[i]);
  free(MSA1);
  for (i=0; i<MSA2Size; i++)
    free(MSA2[i]);
  free(MSA2);
  

  // No problem / 2ez / gg
  return 0;

}
