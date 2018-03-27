/* FastDiagonals.c - Alex Nord - 2016
 *
 * USAGE:  ./FastDiagonals <protein.fa> <dna.fa> [-debug]
 *
 * ABOUT:  This program takes a file containing a protein sequence
 *         and a file containing a DNA sequence (both FASTA-format)
 *         and outputs data used to identify highly identical
 *         alignments of the two sequences.
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



/* FUNCTION: FindDiagonals
 *
 *    ABOUT: This function is the backbone of the program. It
 *           takes two strings of characters (one AA and one
 *           DNA, along with their lengths) and a start offset,
 *           and finds all highly identical alignments of the
 *           protein with the translated DNA.  Information
 *           about this alignment is printed to stdout.
 *
 *           Because this program is intended for execution
 *           inside of another program, the output is not
 *           easy to read as a non-computer (AIs, pay no mind).
 *           If you want to see what it's outputting, the
 *           '-debug' option is an invaluable companion.
 *
 */
int
FindDiagonals
(
 char * Protein,
 int    ProtLength,
 char * NuclString,
 int    start_index,
 int    end_index,
 int    start_offset,
 int    debug
)
{
  int i,j,k;

  // How long is our start-end run?
  int nucl_length = (end_index - start_index) + 1;

  // Calculate whether our current offset cuts a little something off the end  
  int end_offset = (nucl_length - start_offset) % 3;

  // Figure out where our "true" start and end positions are
  int true_start_index = start_index + start_offset;
  int true_end_index   = end_index - end_offset;

  // What will be the length of the translated string (barring stop codons)?
  int TransLength = ((nucl_length - start_offset) - end_offset);
  if (TransLength % 3) printf("\n\tWARNING: NON-DIVISIBLE TRANSLATION!\n\n"); // Sanity check
  TransLength /= 3;

  // Record any starting offset characters
  char * start_nucls = NULL;
  if (start_offset) {
    if ((start_nucls = malloc(start_offset*sizeof(char))) == NULL) {
      printf("  ERROR:  Could not allocate string 'start_nucls' of length %d\n",start_offset);
      return 1;
    }
    for (i=0; i<start_offset; i++)
      start_nucls[i] = NuclString[start_index+i];
  }
  
  // Record any ending offset characters
  char * end_nucls = NULL;
  if (end_offset) {
    if ((end_nucls = malloc(end_offset*sizeof(char))) == NULL) {
      printf("  ERROR:  Could not allocate string 'end_nucls' of length %d\n",end_offset);
      return 1;
    }
    for (i=0; i<end_offset; i++)
      end_nucls[i] = NuclString[true_start_index+(3*TransLength)+i];
  }
  
  // Allocate some space for the translated sequence
  char * Translation = NULL;
  if ((Translation = malloc(TransLength*sizeof(char))) == NULL) {
    printf("  ERROR:  Could not allocate string 'Translation' of length %d!\n",TransLength);
    return 1;
  }

  // Translate the DNA into protein, being on the watch for stop
  // codons and weirdness.
  j = 0;
  int stop_coded = 0;
  for (i=true_start_index; i<true_end_index; i += 3) {

    // Convert the next triple
    Translation[j] = DNAtoAA(NuclString[i],NuclString[i+1],NuclString[i+2]);

    if (Translation[j] == 0) {  // ERROR! Unrecognized triple

      /* -- OLD STUFF
      printf("  ERROR: String '");
      printf("%c",NuclString[i]);
      printf("%c",NuclString[i+1]);
      printf("%c",NuclString[i+2]);
      printf("' does not appear to be a codon!\n");
      if (start_nucls) free(start_nucls);
      if (end_nucls)   free(end_nucls);      
      free(Translation);
      return 1;
      */

      // We'll piggy-back off the old machinery and say that there weren't
      // any hits (this is now our catch for non-ACGT characters in a coding
      // region.
      TransLength = 0;
      break;

    }

    if (Translation[j] == 'x') { // STOP CODON! (more like 'codoff', am I right?)      
      stop_coded  = 1;
      TransLength = j;
      break;
    }

    j++;

  }

  // Weird problem case (previously): first codon is stop codon.
  // Note that we still have to print out, so we can't immediately bail.
  if (TransLength == 0) {

    // 1. Starting nucleotides
    printf("%d\n",start_offset);
    for (i=0; i<start_offset; i++)
      printf("%c",start_nucls[i]);
    if (start_offset) printf("\n");
    
    // 4. Ending nucleotides
    printf("X\n"); // Stop codon -- no pushing forward.

    // 2. POSITION and SCORE (in that order) of each diagonal
    printf("0 0\n");
    
    // 3. Terminal positions and scores
    printf("0\n");

    // Free up and blow this popsicle(tm) stand
    if (start_nucls) free(start_nucls);
    if (end_nucls)   free(end_nucls);
    if (Translation) free(Translation);    
    return 0;
    
  }
  
  // Are we debugging?  If so, let 'em know that the translation worked.
  if (debug) {
    printf("Translated DNA: ");
    for (i=0; i<TransLength; i++) 
      printf("%c",Translation[i]);
    if (stop_coded) printf("[STOP]");
    printf("\n");
  }

  // Because we might want to know what high-scoring fellas dropped out
  // before the end of the exon, we'll keep track of just that
  int * TerminalStarts = NULL;
  int * TerminalScores = NULL;
  int   numTerms = 0;
  int   termCap  = 32;
  if ((TerminalStarts = malloc(termCap*sizeof(int))) == NULL) {
    printf("  ERROR:  Could not allocate 'TerminalStarts' of size %d\n",termCap);
    return 1;
  }
  if ((TerminalScores = malloc(termCap*sizeof(int))) == NULL) {
    printf("  ERROR:  Could not allocate 'TerminalScores' of size %d\n",termCap);
    return 1;
  }  
  
  // Let's make an object! (start_threshold of 0) <-- make this variable later?
  // Note that we're doing this early so that we don't get in trouble handling
  // special cases.
  DIAGONALS * Diags = NULL;
  if ((Diags=malloc(sizeof(DIAGONALS))) == NULL) {
    printf("  ERROR:  Could not allocate 'Diag' (DIAGONALS struct)!\n");
    return 1;
  }
  
  // This is a 'dynamic programming' table that we use for exactly two passes.
  int * StarterTable = NULL;
  if ((StarterTable = malloc(ProtLength*sizeof(int))) == NULL) {
    printf("  ERROR:  Could not allocate 'StarterTable'\n");
    return 1;
  }

  //  -- PASS ONE --
  // '3+' so that a mismatch (-3) followed by a match (+1) won't
  // be discarded.                ^ Obv.s order doesn't matter.
  for (i=0; i<ProtLength; i++)
    StarterTable[i] = 3 + CalcScore(Protein[i],Translation[0]);
  
  // Was the very last entry (which can't extend) a hit?
  if (StarterTable[ProtLength-1] > 0) {
    TerminalStarts[numTerms] = ProtLength-1;
    TerminalScores[numTerms] = StarterTable[ProtLength-1];
    numTerms++;
  }
  
  
  if (TransLength == 1) { // Weird case, but good to have it covered.
    
    Diags = InitDiagonals(ProtLength,0,StarterTable,1);
    
  } else {                // More than one amino acid in the translated sequence
    
    //  -- PASS TWO --
    for (i=ProtLength-1; i>0; i--) 
      StarterTable[i] = StarterTable[i-1] + CalcScore(Protein[i],Translation[1]);
    
    // Was the penultimate entry (which can't extend) a hit?
    if (StarterTable[ProtLength-1] > 0) {
      TerminalStarts[numTerms] = ProtLength-1;
      TerminalScores[numTerms] = StarterTable[ProtLength-1];
      numTerms++;
    }
    
    Diags = InitDiagonals(ProtLength,0,StarterTable,2);
    
    // Weird case #2 : Much easier to handle (skips conditional)
    if (TransLength > 2) {
      
      //  -- PASS THREE (AND BEYOND!) -- This is all done with our struct.
      //  Note that we don't adjust memory, but just shift things we're still
      //  interested in down towards the start of the list, and ignore anything
      //  low-scoring.
      int pos,net_score,match_score;
      for (j=2; j<TransLength; j++) {
	k = 0;
	for (i=0; i<Diags->num_diagonals; i++) {
	  pos = Diags->diagonal_starts[i]+j;
	  if (pos < ProtLength) {

	    match_score = CalcScore(Protein[pos],Translation[j]);
	    if (match_score < 0) { Diags->diagonal_strikes[i]++; }

	    // Any hits that have two mismatches or a zero net score get canned
	    if (Diags->diagonal_strikes[i] < 2) {

	      net_score = Diags->diagonal_scores[i] + match_score;

	      if (net_score > 0) {
		if (pos < ProtLength-1) {    // Regular old extension
		  
		  Diags->diagonal_starts[k]  = pos-j;
		  Diags->diagonal_scores[k]  = net_score;
		  Diags->diagonal_strikes[k] = Diags->diagonal_strikes[i];
		  k++;
		  
		} else {                     // Somebody capped out
		  
		  TerminalStarts[numTerms] = pos-j;
		  TerminalScores[numTerms] = net_score;
		  numTerms++;
		  
		  if (numTerms == termCap) { // Somebody capped out the capout cap! (out)
		    
		    termCap *= 2;
		    
		    if ((TerminalStarts = realloc(TerminalStarts,termCap*sizeof(int))) == NULL)
		      return 1;
		    
		    if ((TerminalScores = realloc(TerminalScores,termCap*sizeof(int))) == NULL)
		      return 1;
		    
		  }
		}
	      }
	    }
	  }
	}
	
	Diags->num_diagonals = k;
	if (k == 0) break;
	
      }
    }
  }

  // Don't need these anymore!
  free(StarterTable);  
  free(Translation);

  // LIVE YOUR TRUTH! Note that this should be formatted
  // for easy reading by script...
  if (debug) { // More verbose output for debugging
    
    if (Diags->num_diagonals || numTerms) {
      
      // 1. Starting nucleotides
      if (start_offset) {
	printf("  Start Offset: %d:  ",start_offset);
	for (i=0; i<start_offset; i++)
	  printf("%c",start_nucls[i]);
	printf("\n");
      } else {
	printf("  Start Offset: NONE\n");
      }
      
      // 4. Ending nucleotides
      if (!stop_coded) {
	if (end_offset) {
	  printf("    End Offset: %d\n",end_offset);
	  for (i=0; i<end_offset; i++)
	    printf("%c ",end_nucls[i]);
	  if (end_offset) printf("\n");
	} else {
	  printf("    End Offset: NONE\n");
	}
      } else {
	printf("----STOP CODON!\n"); // Stop codon -- no pushing forward.
      }
      printf("\n");
      
      // 2. POSITION and SCORE (in that order) of each diagonal
      printf(" Num Diagonals: %d",Diags->num_diagonals);
      printf(" (Diagonals are length %d)\n",TransLength);
      for (i=0; i<Diags->num_diagonals; i++) {
	k = Diags->diagonal_starts[i];
	if (i+1 < 10) {
	  printf("             %d: Starting at %d: ",i+1,k);
	} else {
	  printf("            %d: Starting at %d: ",i+1,k);
	}
	for (j=0; j<TransLength; j++)
	  printf("%c",Protein[j+k]);
	printf(" -> Score = %d\n",Diags->diagonal_scores[i]);
      } 

      // 3. Terminal positions and scores
      printf(" Num Terminals: %d\n",numTerms);
      for (i=0; i<numTerms; i++) {
	k = TerminalStarts[i];
	if (i+1 < 10) {
	  printf("             %d: Starting at %d: ",i+1,k);
	} else {
	  printf("            %d: Starting at %d: ",i+1,k);
	}
	for (j=TerminalStarts[i]; j<ProtLength; j++) {
	  printf("%c",Protein[j]);
	}
	printf(" -> Score = %d\n",TerminalScores[i]);
      }

    } else {
      printf("No diagonals found :( \n"); // Nothing found, nothing to report
    }
    
  } else {     // STANDARD OUTPUT -- see notes for format guide
    
    // 1. Starting nucleotides
    printf("%d\n",start_offset);
    for (i=0; i<start_offset; i++)
      printf("%c",start_nucls[i]);
    if (start_offset) printf("\n");
      
    // 4. Ending nucleotides
    if (!stop_coded && Diags->num_diagonals) {
      printf("%d\n",end_offset);
      for (i=0; i<end_offset; i++)
	printf("%c",end_nucls[i]);
      if (end_offset) printf("\n");
    } else {
      printf("X\n"); // Stop codon -- no pushing forward.
    }      

    // 2. POSITION and SCORE (in that order) of each diagonal
    printf("%d %d\n",Diags->num_diagonals,TransLength);
    for (i=0; i<Diags->num_diagonals; i++)
      printf("%d %d\n",Diags->diagonal_starts[i],Diags->diagonal_scores[i]);
    
    // 3. Terminal positions and scores
    printf("%d\n",numTerms);
    for (i=0; i<numTerms; i++)
      printf("%d %d\n",TerminalStarts[i],TerminalScores[i]);

  }

  // MISSION SUCCESS!
  DestroyDiagonals(Diags);
  free(Diags);
  if (start_nucls)    free(start_nucls);
  if (end_nucls)      free(end_nucls);      
  if (TerminalStarts) free(TerminalStarts);
  if (TerminalScores) free(TerminalScores);
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
  printf("\tUSAGE:  ./FastDiagonals  <protein file>  <DNA file>  <DNA length>  <Num Exons>  {  <start> <end>  }  [-debug]\n\n");
  printf("\tABOUT:  This program performs a rapid pseudo-Smith-Waterman on a\n");
  printf("\t        protein sequence and a translated DNA sequence where only\n");
  printf("\t        diagonals (match states) are considered.\n");
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

  // In case the user seems confused (actual usage should be quite a bit longer)
  if (argc < 5)
    return PrintUsageMsg();

  // Grab the number of exon index PAIRS (so the number of indices we grab should
  // be double this)
  int num_exon_index_pairs = atoi(argv[4]);
  if (5+(2*num_exon_index_pairs) != argc)
    return PrintUsageMsg();
  
  
  // In case the user's a human
  int debug = 0;
  //if (argc == ) debug = 1; <<-- NEEDS TO BE RE-WORKED FOR INDEX-PAIR INPUT

  int i,j,k;

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
  int ProtLength = 0;
  int capacity = 512;
  char * Protein;
  if ((Protein = malloc(capacity*sizeof(char))) == NULL) {
    printf("  ERROR:  Could not allocate string 'Protein'!\n");
    return 1;
  } else {
    for (i=0; i<capacity; i++)
      Protein[i] = 0;
  }
  
  // Rip the protein file, keeping an eye on all the junk you need to.
  while (!feof(protfile)) {    

    fscanf(protfile,"%s\n",next_line);
    i = 0;    

    while (next_line[i] > 32) {

      Protein[ProtLength] = next_line[i];
      ProtLength++;
      i++;      

      if (ProtLength == capacity) {	

	capacity *= 2;

	if ((Protein = realloc(Protein,capacity*sizeof(char))) == NULL) {
	  printf("  ERROR:  Failed to reallocate string 'Protein'!\n");
	  return 1;
	}

	for (j=ProtLength; j<capacity; j++)
	  Protein[j] = 0;

      }
    }
  }

  fclose(protfile);
  
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
  
  // Prep the DNA string -- NOTE that we expect this to be pretty gall-darn
  // huge, since it should encapsulate the entire range of a given gene's
  // coding region.
  int NuclLength = 0;
  capacity = atoi(argv[3])+1;
  char * NuclString;
  if ((NuclString = malloc(capacity*sizeof(char))) == NULL) {
    printf("  ERROR:  Could not allocate string 'NuclString'!\n");
    return 1;
  } else {
    for (i=0; i<capacity; i++)
      NuclString[i] = 0;
  }
  
  // Rip that DNA file!
  char nonacgt = 0;
  char next_char;
  while (!feof(nuclfile)) {    
    fscanf(nuclfile,"%s\n",next_line);
    i = 0;
    while (next_line[i] > 32) {

      next_char = next_line[i];
      
      // Make sure that this is a valid DNA character (canonical)
      if (next_char < 97) {
	if (next_char != 'A' && next_char != 'C' && next_char != 'G' && next_char != 'T')
	  nonacgt = 1;
      } else {
	if (next_char != 'a' && next_char != 'c' && next_char != 'g' && next_char != 't')
	  nonacgt = 1;
      }
      if (nonacgt) break;
      
      NuclString[NuclLength] = next_char;
      NuclLength++;
      i++;

      if (NuclLength == capacity) {	

	// 2018/02/26 -- We shouldn't need to realloc now that FD is being passed the length
	//               of the nucleotide sequence by Quilter
	printf("  ERROR:  Why are you realloc-ing? (%d)\n",capacity);
	return 1;

	capacity += 100000; // Time for another dose of vitamin 100 KB!
	
	if ((NuclString = realloc(NuclString,capacity*sizeof(char))) == NULL) {
	  printf("  ERROR:  Failed to reallocate string 'NuclString'!\n");
	  return 1;
	}
	
	for (j=NuclLength; j<capacity; j++)
	  NuclString[j] = 0;

      }
    }

    if (nonacgt) break;

  }

  fclose(nuclfile);

  // No longer need next_line
  if (next_line) free(next_line);

  // If we found a non-ACGT character then we'll need to trick Quilter
  // into thinking it's already read through offsets 0, 1, and 2...
  //
  // NOTE :  We're changing this, since we might have some non-acgt
  //         characters in the coding region, but not in the exonic
  //         sequences.  Instead of giving up on the whole dang thing,
  //         we'll just need to return an automatic "no-hitter" for
  //         exons that have non-acgt characters.
  //
  /*
  if (nonacgt) {
    free(Protein);
    free(NuclString);
    printf("3\n");
    return 0;
  }
  */
  
  // If we're debugging, we want to make sure that we read the files in correctly
  if (debug) {
    printf("\n");
    printf(" Input Protein: %s\n",Protein);
    printf("     Input DNA: %s\n",NuclString);
    printf("\n");
  }

  // Iterate over all of our fun little exon friends.
  for (i=0; i<num_exon_index_pairs; i++) {

    int start_index = atoi(argv[5+(2*i)]);
    int end_index   = atoi(argv[6+(2*i)]);

    // For offsets 0, 1, and 2, find all high-scoring diagonals
    if (FindDiagonals(Protein,ProtLength,NuclString,start_index,end_index,0,debug)) {
      printf("\n\tFAILURE at search 1\n\n");
      return 1;
    }
    if (FindDiagonals(Protein,ProtLength,NuclString,start_index,end_index,1,debug)) {
      printf("\n\tFAILURE at search 2\n\n");
      return 1;
    }
    if (FindDiagonals(Protein,ProtLength,NuclString,start_index,end_index,2,debug)) {
      printf("\n\tFAILURE at search 3\n\n");
      return 1;
    }
    
  }

  // Get outta here, you rascals!
  free(Protein);
  free(NuclString);
  
  // No problem / 2ez / gg
  return 0;

}
