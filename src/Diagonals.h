/* Diagonals.h - Alex Nord - 2016
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef __DIAGONALS_H__
#define __DIAGONALS_H__


/* STRUCT: DIAGONALS
 *
 * Used to store information about high-scoring diagonals in our alignment
 * of the input protein sequence to one of the three DNA translations.
 *
 */
typedef struct _DIAGONALS {
  int   diagonal_length;  // How far have we extended the recorded diagonals?
  int   num_diagonals;    // How many diagonals have been recorded, so far?
  int * diagonal_starts;  // What is the starting point of each diagonal (in the protein)?
  int * diagonal_scores;  // What is the score of each diagonal?
  int * diagonal_strikes; // How many mismatches has a diagonal had?
} DIAGONALS;



DIAGONALS * InitDiagonals (int seq_length, int start_threshold, int * diag_scores, int diag_length);
void DestroyDiagonals (DIAGONALS * diags);
char DNAtoAA (char first, char second, char third);
int CalcScore (char one, char two);


#endif
