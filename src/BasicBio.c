#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "BasicBio.h"

/* List of Functions
 *
 *   MBB_MinInt
 *   MBB_MaxInt
 *   MBB_MinFloat
 *   MBB_MaxFloat
 *   MBB_SortInts
 *   MBB_SortFloats
 *   MBB_DNAtoNum
 *   MBB_TranslateCodon
 *   MBB_AminoAliScore
 *
 */

// Function: MBB_MinInt
int MBB_MinInt (int a, int b) {
  if (a < b) return a;
  return b;
}

// Function: MBB_MaxInt
int MBB_MaxInt (int a, int b) {
  if (a > b) return a;
  return b;
}

// Function: MBB_MinFloat
float MBB_MinFloat (float a, float b) {
  if (a < b) return a;
  return b;
}

// Function: MBB_MaxFloat
float MBB_MaxFloat (float a, float b) {
  if (a > b) return a;
  return b;
}



/*
 * Function: MBB_SortInts
 *
 */
void MBB_SortInts (int * Vals, int * Index, int num_vals) {

  int i,j,k,j_lim,k_lim;

  int * Read  = malloc(num_vals*sizeof(int));
  int * Write = malloc(num_vals*sizeof(int));
  int * Swap;
  
  // Initialize the index
  for (i=0; i<num_vals; i++)
    Read[i] = i;

  // Merge sort
  int blocksize  = 1;
  int num_blocks = num_vals / blocksize;
  int pos        = 0;
  while (num_blocks) {

    for (i=0; i+blocksize<num_vals; i+=2*blocksize) {

      j     = i;
      j_lim = j+blocksize;
      k     = j_lim;
      k_lim = k+blocksize;
      if (k_lim>num_vals)
	k_lim = num_vals;
      
      while (j<j_lim && k<k_lim) {
	if (Vals[Read[j]]<Vals[Read[k]])
	  Write[pos++] = Read[j++];
	else
	  Write[pos++] = Read[k++];
      }

      while (j<j_lim)
	Write[pos++] = Read[j++];

      while (k<k_lim)
	Write[pos++] = Read[k++];
      
    }

    // If we didn't have enough to fill up a full block, we'll still need to
    // copy over the final values into 'Write'
    while (i<num_vals)
      Write[pos++] = Read[i++];

    // Read becomes write, write becomes read, the first shall be last, etc!
    Swap  = Read;
    Read  = Write;
    Write = Swap;

    // Brace yourselves for another iteration!
    num_blocks /= 2;
    blocksize  *= 2;
    pos         = 0;
    
  }

  // Our final sorted index is stored in 'Read'
  for (i=0; i<num_vals; i++)
    Index[i] = Read[i];
  
  free(Read);
  free(Write);

}



/*
 * Function: MBB_SortFloats
 *
 */
void MBB_SortFloats (float * Vals, int * Index, int num_vals) {

  int i,j,k,j_lim,k_lim;

  int * Read  = malloc(num_vals*sizeof(int));
  int * Write = malloc(num_vals*sizeof(int));
  int * Swap;
  
  // Initialize the index
  for (i=0; i<num_vals; i++)
    Read[i] = i;

  // Merge sort
  int blocksize  = 1;
  int num_blocks = num_vals / blocksize;
  int pos        = 0;
  while (num_blocks) {

    for (i=0; i+blocksize<num_vals; i+=2*blocksize) {

      j     = i;
      j_lim = j+blocksize;
      k     = j_lim;
      k_lim = k+blocksize;
      if (k_lim>num_vals)
	k_lim = num_vals;
      
      while (j<j_lim && k<k_lim) {
	if (Vals[Read[j]]<Vals[Read[k]])
	  Write[pos++] = Read[j++];
	else
	  Write[pos++] = Read[k++];
      }

      while (j<j_lim)
	Write[pos++] = Read[j++];

      while (k<k_lim)
	Write[pos++] = Read[k++];
      
    }

    // If we didn't have enough to fill up a full block, we'll still need to
    // copy over the final values into 'Write'
    while (i<num_vals)
      Write[pos++] = Read[i++];

    // Read becomes write, write becomes read, the first shall be last, etc!
    Swap  = Read;
    Read  = Write;
    Write = Swap;

    // Brace yourselves for another iteration!
    num_blocks /= 2;
    blocksize  *= 2;
    pos         = 0;
    
  }

  // Our final sorted index is stored in 'Read'
  for (i=0; i<num_vals; i++)
    Index[i] = Read[i];
  
  free(Read);
  free(Write);

}



/*
 * Function: MBB_DNAtoNum
 *
 */
int MBB_DNAtoNum (char residue) {

  // Set to uppercase
  if (residue > 96)
    residue -= 32;

  // Find character
  if (residue < 69) {
    if (residue == 65) return 0;
    if (residue == 67) return 1;
  } else {
    if (residue == 71) return 2;
    if (residue == 84) return 3;
    if (residue == 85) return 3;
  }
  return -1;
  
}

/*
 * Function: MBB_AminoToIndex
 *
 */
int MBB_AminoToIndex (char residue) {
  if (residue >= 97) residue -= 32;
  if (residue < 65 || residue > 90)
    return -1;
  return MBB_AMINO_INDEX[residue];
}


/*
 * Function: MBB_TranslateCodon
 *
 */
char MBB_TranslateCodon (char * Codon) {
  int index = MBB_DNAtoNum(*Codon)<<4;
  index    += MBB_DNAtoNum(*(Codon+1))<<2;
  index    += MBB_DNAtoNum(*(Codon+2));
  return MBB_CODON_TABLE[index];
}


/*
 * Function: MBB_AminoAliScore
 *
 */
float MBB_AminoAliScore (char amino1, char amino2) {

  // Make sure we're uppercase
  if (amino1 > 96) amino1 -= 32;
  if (amino2 > 96) amino2 -= 32;

  int index1 = MBB_AMINO_INDEX[amino1-65];
  int index2 = MBB_AMINO_INDEX[amino2-65];

  // Invalid chars get the big -inf
  if (index1 < 0 || index2 < 0)
    return MBB_NINF;

  return MBB_BLOSUM62[21 * index1 + index2];
  
}
