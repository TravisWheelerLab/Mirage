/* Diagonals.c - Alex Nord - 2016
 *
 * USAGE:  N/A
 *
 * ABOUT:  This file contains a handful of simple functions used
 *         by "FastDiagonals."
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Diagonals.h"



/* FUNCTION: InitDiagonals
 *
 *    ABOUT: This function generates the DIAGONALS struct used during
 *           alignment.  It takes the length of the (protein) sequence,
 *           a minimum score threshold for starting a diagonal, and a
 *           list of diagonal scores (of length seq_length), and generates
 *           a DIAGONALS struct listing only those positions that reached
 *           the score threshold.
 */
DIAGONALS *
InitDiagonals
(
 int   seq_length,
 int   start_threshold,
 int * diag_scores,
 int   diag_length
)
{
  int i,j,k;
  int capacity;
  
  // 1. Allocate the struct itself
  DIAGONALS * diags;  
  if ((diags = malloc(sizeof(DIAGONALS))) == NULL) return NULL;
  
  // Allocate for and record the diagonal starting indices and scores.
  j = 0;

  capacity = 32; // This should be plenty of capacity, since we're starting 2 in
  if ((diags->diagonal_starts  = malloc(capacity * sizeof(int))) == NULL) return NULL;
  if ((diags->diagonal_scores  = malloc(capacity * sizeof(int))) == NULL) return NULL;
  if ((diags->diagonal_strikes = malloc(capacity * sizeof(int))) == NULL) return NULL;

  for (i=diag_length-1; i<seq_length; i++) {

    if (diag_scores[i] > start_threshold) {    
      
      diags->diagonal_starts[j]  = i-1;
      diags->diagonal_scores[j]  = diag_scores[i];
      diags->diagonal_strikes[j] = 0;
      j++;
      
      // Hopefully we won't need to resize, but just in case
      if (j == capacity) {
	
	capacity *= 2;
	
	if ((diags->diagonal_starts = realloc(diags->diagonal_starts,capacity*sizeof(int))) == NULL) 
	  return NULL;
	
	if ((diags->diagonal_scores = realloc(diags->diagonal_scores,capacity*sizeof(int))) == NULL) 
	  return NULL;

	if ((diags->diagonal_strikes = realloc(diags->diagonal_strikes,capacity*sizeof(int))) == NULL) 
	  return NULL;

      }
    }
  }

  diags->num_diagonals   = j;
  diags->diagonal_length = diag_length;
  
  // It was tough, but together we can overcome anything!
  return diags;

}


/* FUNCTION: DestroyDiagonals
 *
 *    ABOUT: This function frees up all of the 'interior' components
 *           of an instance of a DIAGONALS struct.  At an earlier time
 *           there was more to this struct, so making a whole function
 *           seemed slightly more prudent than it looks now.
 *
 */
void
DestroyDiagonals
( DIAGONALS * diags )
{
  free(diags->diagonal_starts);
  free(diags->diagonal_scores);
  free(diags->diagonal_strikes);
}



/* FUNCTION: DNAtoAA
 *
 *    ABOUT: MASSIVE conversion table.  Nothing more.
 *
 */
char 
DNAtoAA
(
 char first,
 char second,
 char third
)
{
  if (first  > 96) first  -= 32;
  if (second > 96) second -= 32;
  if (third  > 96) third  -= 32;

  if (first == 'A') {
    if (second == 'A') {
      if (third == 'A') {
	return 'K'; // AAA
      } else if (third == 'C') {
	return 'N'; // AAC
      } else if (third == 'G') {
	return 'K'; // AAG	
      } else if (third == 'T' || third == 'U') {
	return 'N'; // AAT	
      } else {
	return 0;   // ERROR
      }      
    } else if (second == 'C') {
      if (third == 'A') {
	return 'T'; // ACA	
      } else if (third == 'C') {
	return 'T'; // ACC	
      } else if (third == 'G') {
	return 'T'; // ACG	
      } else if (third == 'T' || third == 'U') {
	return 'T'; // ACT	
      } else {
	return 0;   // ERROR
      }
    } else if (second == 'G') {
      if (third == 'A') {
	return 'R'; // AGA
      } else if (third == 'C') {
	return 'S'; // AGC	
      } else if (third == 'G') {
	return 'R'; // AGG	
      } else if (third == 'T' || third == 'U') {
	return 'S'; // AGT	
      } else {
	return 0;   // ERROR	
      }
    } else if (second == 'T' || second == 'U') {
      if (third == 'A') {
	return 'I'; // ATA	
      } else if (third == 'C') {
	return 'I'; // ATC	
      } else if (third == 'G') {
	return 'M'; // ATG	
      } else if (third == 'T' || third == 'U') {
	return 'I'; // ATT
      } else {
	return 0;   // ERROR
      }
    } else {
      return 0;     // ERROR
    }    
  } else if (first == 'C') {
    if (second == 'A') {
      if (third == 'A') {
	return 'Q'; // CAA	
      } else if (third == 'C') {
	return 'H'; // CAC	
      } else if (third == 'G') {
	return 'Q'; // CAG	
      } else if (third == 'T' || third == 'U') {
	return 'H'; // CAT	
      } else {
	return 0;   // ERROR	
      }
    } else if (second == 'C') {
      if (third == 'A') {
	return 'P'; // CCA	
      } else if (third == 'C') {
	return 'P'; // CCC	
      } else if (third == 'G') {
	return 'P'; // CCG
      } else if (third == 'T' || third == 'U') {
	return 'P'; // CCT
      } else {
	return 0;   // ERROR	
      }      
    } else if (second == 'G') {
      if (third == 'A') {
	return 'R'; // CGA	
      } else if (third == 'C') {
	return 'R'; // CGC	
      } else if (third == 'G') {
	return 'R'; // CGG
      } else if (third == 'T' || third == 'U') {
	return 'R'; // CGT	
      } else {
	return 0;   // ERROR
      }      
    } else if (second == 'T' || second == 'U') {
      if (third == 'A') {
	return 'L'; // CTA	
      } else if (third == 'C') {
	return 'L'; // CTC
      } else if (third == 'G') {
	return 'L'; // CTG
      } else if (third == 'T' || third == 'U') {
	return 'L'; // CTT
      } else {
	return 0;   // ERROR
      }      
    } else {
      return 0;     // ERROR
    }    
  } else if (first == 'G') {
    if (second == 'A') {
      if (third == 'A') {
	return 'E'; // GAA
      } else if (third == 'C') {
	return 'D'; // GAC
      } else if (third == 'G') {
	return 'E'; // GAG
      } else if (third == 'T' || third == 'U') {
	return 'D'; // GAT
      } else {
	return 0;   // ERROR
      }      
    } else if (second == 'C') {
      if (third == 'A') {
	return 'A'; // GCA
      } else if (third == 'C') {
	return 'A'; // GCC (the codon used to compile C programs)
      } else if (third == 'G') {
	return 'A'; // GCG
      } else if ( third == 'T' || third == 'U') {
	return 'A'; // GCT
      } else {
	return 0;   // ERROR
      }      
    } else if ( second == 'G') {
      if (third == 'A') {
	return 'G'; // GGA
      } else if (third == 'C') {
	return 'G'; // GGC
      } else if (third == 'G') {
	return 'G'; // GGG
      } else if (third == 'T' || third == 'U') {
	return 'G'; // GGT
      } else {
	return 0;   // ERROR
      }      
    } else if (second == 'T' || second == 'U') {
      if (third == 'A') {
	return 'V'; // GTA (Rated M for Mature)
      } else if (third == 'C') {
	return 'V'; // GTC
      } else if (third == 'G') {
	return 'V'; // GTG
      } else if (third == 'T' || third == 'U') {
	return 'V'; // GTT
      } else {
	return 0;   // ERROR
      }      
    } else {
      return 0;     // ERROR
    }    
  } else if (first == 'T' || first == 'U') {
    if (second == 'A') {
      if (third == 'A') {
	return 'x'; // TAA // STOP
      } else if (third == 'C') {
	return 'Y'; // TAC
      } else if (third == 'G') {
	return 'x'; // TAG // STOP
      } else if (third == 'T' || third == 'U') {
	return 'Y'; // TAT
      } else {
	return 0;   // ERROR
      }      
    } else if (second == 'C') {
      if (third == 'A') {
	return 'S'; // TCA
      } else if (third == 'C') {
	return 'S'; // TCC
      } else if (third == 'G') {
	return 'S'; // TCG
      } else if (third == 'T' || third == 'U') {
	return 'S'; // TCT
      } else {
	return 0;   // ERROR
      }      
    } else if (second == 'G') {
      if (third == 'A') { 
	return 'x'; // TGA // STOP 
      } else if (third == 'C') {
	return 'C'; // TGC
      } else if (third == 'G') {
	return 'W'; // TGG
      } else if (third == 'T' || third == 'U') { 
	return 'C'; // TGT
      } else {
	return 0;   // ERROR
      }      
    } else if (second == 'T' || second == 'U') {
      if (third == 'A') {
	return 'L'; // TTA
      } else if (third == 'C') {
	return 'F'; // TTC
      } else if (third == 'G') {
	return 'L'; // TTG
      } else if (third == 'T' || third == 'U') {
	return 'F'; // TTT
      } else {
	return 0;   // ERROR
      }      
    } else {
      return 0;     // ERROR
    }    
  } else {
    return 0;       // ERROR
  }
}



/* FUNCTION: CalcScore
 *
 *    ABOUT: This function takes two chars and computes their
 *           alignment score.  It's pretty challenging, so I
 *           won't even try to explain what's going on here.
 *
 */
int 
CalcScore
(char one, char two)
{
  // Wow!
  // Much calculations!
  // So programmer!
  if (one == 'x' || two == 'x') 
    return -100;
  if (one == two) 
    return 1;
  return -3;
  // Numbers R GOOD!
}


