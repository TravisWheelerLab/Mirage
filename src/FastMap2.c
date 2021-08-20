/*  
 * FastMap2.c
 * 
 * About: This program is a potential replacement for the implementation of 'FastMap'
 *        from the original release of Mirage.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "BasicBio.h"

// FUNCTION LIST
//
// + PrintUsage
// + ParseNucleotideFile
// + ParseProteinFile
// + TranslateRF
// + PrintMap
// + BlockScan
// + OriginalFastMap

// How many mismatches are allowed (on each side) as we extend a diagonal?
const int MAX_MISMATCH = 2;

// What blocksize do we use during diagonal identification?
const int BLOCK_SIZE = 6;
  
// PrintUsage
int PrintUsage () {
  printf("\n");
  printf("  USAGE:  ./FastMap2.c [protein-seq-file] [num-protein-seqs] [chromosomal-seq-file] {exon-start-index exon-end-index}\n");
  printf("\n");
  printf("  NOTES:  The 'exon-...-index' values are wrt genome, not relative to\n");
  printf("          the range of chromosomal sequence present in the file.\n");
  printf("          The range is treated as inclusive.\n");
  printf("\n");
  printf("          For reverse strand mapping, start indices should be the\n");
  printf("          larger of the two numbers.\n");
  printf("\n");
  printf("          OUTPUT amino acid ranges are [0..protein_length) and reading\n");
  printf("          frame values are [0..2].  Nucleotide coordinates don't include\n");
  printf("          the flanking 'non-coding' nucleotides on either side of the\n");
  printf("          mapped region.\n");
  printf("\n");
  return 1;
}


/*
 *  Function:  ParseNucleotideFile
 *
 *  About:  This function parses a FASTA-formatted sequence file
 *          containing a single DNA sequence extracted from a chromosome
 *          using sfetch (or at least using the naming convention
 *          of '>[chr]/[start]-[end]')
 */
char * ParseNucleotideFile
(
 char * fname,
 int  * nucl_seq_len,
 int  * nucl_start_index,
 int  * nucl_end_index
 ){

  int i,j;

  FILE * inf;
  inf = fopen(fname,"r");
  if (inf == NULL) {
    fprintf(stderr,"\n  ERROR:  Failed to open nucleotide sequence file '%s'\n\n",fname);
    return NULL;
  }

  char * line;
  if ((line=malloc(512*sizeof(char)))==NULL) {
    fprintf(stderr,"\n  ERROR:  Malloc error [loc1]\n\n");
    return NULL;
  }

  // First line is going to be the sequence name, formatted using Easel standards.
  // We can figure out the extracted region since the name should be structured as
  // >[chr]/[start]-[end]
  fscanf(inf,"%s\n",line);
  i=1;
  while (line[i] != '/')
    i++;
  
  int start = line[++i] - '0';
  while (line[++i] != '-')
    start = start * 10 + (line[i] - '0');
  *nucl_start_index = start;

  int end = line[++i] - '0';
  while (line[++i] >= '0')
    end = end * 10 + (line[i] - '0');
  *nucl_end_index = end;

  // Now that we know the range, we can compute the sequence length and allocate
  *nucl_seq_len = abs(end-start)+1;
  char * NuclSeq;
  if ((NuclSeq=malloc((*nucl_seq_len+1) * sizeof(char)))==NULL) {
    fprintf(stderr,"\n  ERROR:  Malloc error [loc2]\n\n");
    return NULL;
  }

  // Read in the rest of the sequence
  i=0;
  while (!feof(inf)) {
    fscanf(inf,"%s\n",line);
    strcpy(NuclSeq+i,line);
    i += strlen(line);
  }

  // Wrap it on up!
  free(line);
  fclose(inf);

  return NuclSeq;
  
}



/*
 *  Function:  ParseProteinFile
 *  
 *  About:  Duh
 *
 */
char ** ParseProteinFile (char * fname, int num_seqs)
{
  FILE * inf;
  inf = fopen(fname,"r");
  if (inf == NULL) {
    fprintf(stderr,"\n  ERROR:  Failed to open protein sequence file '%s'\n\n",fname);
    return NULL;
  }

  char * line;
  if ((line=malloc(512*sizeof(char)))==NULL) {
    fprintf(stderr,"\n  ERROR:  Malloc error [loc3]\n\n");
    return NULL;
  }

  // We shouldn't ever need more than 40kb for a protein sequence,
  // unless there's something bigger than titin that's somehow
  // escaped discovery...
  char * seqbuffer;
  if ((seqbuffer=malloc(40000*sizeof(char)))==NULL) {
    fprintf(stderr,"\n  ERROR:  Malloc error [loc4]\n\n");
    return NULL;
  }

  // Eat our way up to the start of the first sequence
  while (!feof(inf)) {
    fscanf(inf,"%s\n",line);
    if (line[0] == '>')
      break; 
  }

  // If we hit the end of the file something's not right...
  if (feof(inf)) {
    fprintf(stderr,"\n  ERROR:  No sequences found in protein sequence file\n\n");
    return NULL;
  }

  // Where we'll put the sequences we read in.
  char ** ProtSeqs;
  if ((ProtSeqs=malloc(num_seqs*sizeof(char *)))==NULL) {
    fprintf(stderr,"\n  ERROR:  Malloc error [loc5]\n\n");
    return NULL;
  }

  // Let's rock 'n' roll!
  num_seqs   = 0;
  int seqlen = 0;
  while (!feof(inf)) {

    fscanf(inf,"%s\n",line);

    if (line[0] == '>') {

      // Allocate and store the last sequence
      if ((ProtSeqs[num_seqs]=malloc((seqlen+1)*sizeof(char)))==NULL) {
	fprintf(stderr,"\n  ERROR:  Malloc error [loc6,seq%d]\n\n",num_seqs);
	return NULL;
      }

      seqbuffer[seqlen] = 0;
      strcpy(ProtSeqs[num_seqs],seqbuffer);

      // Make sure the whole dang sequence is uppercase
      while (seqlen) {
	seqlen--;
	if (ProtSeqs[num_seqs][seqlen] > 96)
	  ProtSeqs[num_seqs][seqlen] -= 32;
      }

      num_seqs++;
      
    } else {
      strcpy(seqbuffer+seqlen,line);
      seqlen += strlen(line);
    }

  }

  // There should always be a final sequence, but we'll handle the possibility
  // of a name without a sequence just in case...
  //
  // This is the same as the above code, so no comments
  //
  if (seqlen) {
    if ((ProtSeqs[num_seqs]=malloc((seqlen+1)*sizeof(char)))==NULL) {
      fprintf(stderr,"\n  ERROR:  Malloc error [loc6,seq%d]\n\n",num_seqs);
      return NULL;
    }
    seqbuffer[seqlen] = 0;
    strcpy(ProtSeqs[num_seqs],seqbuffer);
    num_seqs++;
  } else {
    fprintf(stderr,"\n  ERROR:  Empty sequence name in protein file\n\n");
    return NULL;
  }
  
  // Clean it up
  free(seqbuffer);
  free(line);
  fclose(inf);

  // Done!
  return ProtSeqs;
  
}



/*
 *  Function: TranslateRF
 *
 */
void TranslateRF (char * NuclSeq, int start, int end, char * RF)
{
  int i=0;
  while (start<end) {
    RF[i++] = MBB_TranslateCodon(&NuclSeq[start]);
    start += 3;
  }
  RF[i] = 0;
}




/*
 *  Function: PrintMap
 *
 */
void PrintMap
( char * ProtSeq,
  int    prot_id,
  int    prot_start_index,
  int    exon_id,
  int    rf_index,
  char * ORF,
  int    map_start_index,
  int    map_len,
  char * NuclSeq,
  int    nucl_seq_len,
  int    rel_start_nucl,
  int    chr_start_nucl,
  float  score,
  int    revcomp
 ){

  int i;

  // ALL BOUNDS ARE INCLUSIVE!
  int prot_end_index = prot_start_index + (map_len-1);
  int map_end_index  = map_start_index  + (map_len-1);

  rel_start_nucl  += 3 * map_start_index;
  int rel_end_nucl = rel_start_nucl + 3*map_len - 1;

  int chr_end_nucl;
  if (revcomp) {
    chr_start_nucl -= 3 * map_start_index;
    chr_end_nucl    = chr_start_nucl - 3*map_len + 1;
  } else {
    chr_start_nucl += 3 * map_start_index;
    chr_end_nucl    = chr_start_nucl + 3*map_len - 1;
  }
  
  // NOTE: I'm choosing for there to be 17 buffer nucleotides so
  //       that, when stitching hits together, we can have
  //       (up to) an extension of 5 full aminos from each hit,
  //       and also observe whether we have AG/GT at full 5aa extension.
  int num_buffer_nucls = 17; // *17*

  //       ... and, as you might expect, this puts a little restriction
  //       on the hits we can report (mainly applies to nonstandard
  //       chromosomal sequences)
  if (rel_start_nucl - num_buffer_nucls < 0
      || rel_end_nucl + num_buffer_nucls > nucl_seq_len)
    return;

  // Metadata
  printf("\n");
  printf("Protein Num   : %d\n",prot_id);
  printf("Exon Num      : %d\n",exon_id);
  printf("Reading Frame : %d\n",rf_index);
  printf("Score         : %f\n",score);
  printf("Amino Range   : %d..%d\n",prot_start_index,prot_end_index);
  printf("Mapped Nucl.s : %d..%d\n",chr_start_nucl,chr_end_nucl);

  // Query Protein
  printf("Query Protein :  ");
  for (i=0; i<num_buffer_nucls; i++) printf(" ");
  while (prot_start_index <= prot_end_index)
    printf(" %c ",ProtSeq[prot_start_index++]);
  printf("\n");

  // Nucls
  printf("Nucleotide Seq:  ");
  for (i=num_buffer_nucls; i>0; i--) {
    if (NuclSeq[rel_start_nucl-i] >= 96) printf("%c",NuclSeq[rel_start_nucl-i]);
    else                                 printf("%c",NuclSeq[rel_start_nucl-i]+32);
  }
  while (rel_start_nucl <= rel_end_nucl) {
    if (NuclSeq[rel_start_nucl] >= 96) printf("%c",NuclSeq[rel_start_nucl]-32);
    else                               printf("%c",NuclSeq[rel_start_nucl]);
    rel_start_nucl++;
  }
  for (i=0; i<num_buffer_nucls; i++) {
    if (NuclSeq[rel_start_nucl+i] >= 96) printf("%c",NuclSeq[rel_start_nucl+i]);
    else                                 printf("%c",NuclSeq[rel_start_nucl+i]+32);
  }
  printf("\n");

  // Translation
  printf("Translated Seq:  ");
  for (i=0; i<num_buffer_nucls; i++) printf(" ");
  while (map_start_index <= map_end_index)
    printf(" %c ",ORF[map_start_index++]);
  printf("\n\n");

}



/*
 *  Function: BlockScan
 *
 *  About: Experimental approach to mapping to an exon.
 *
 *         We'll first break the exon into blocks of a given
 *         length, find which blocks map perfectly to parts
 *         of the protein, and then see if we can extend a
 *         lil' bit into the adjacent (imperfectly mapped) blocks.
 *
 *         The reason why we'd be interested in doing this is
 *         essentially to catch ARFs, since they aren't going to
 *         necessarily cover the ORFs from SRFs.
 *
 */
void BlockScan
( char * ProtSeq,
  int    prot_id,
  int    exon_id,
  int    rf_index,
  char * ORF,
  int    orf_len,
  int    micro_exon,
  char * NuclSeq,
  int    nucl_seq_len,
  int    rel_start_nucl,
  int    chr_start_nucl,
  int    revcomp
 ){

  int i,j;

  int prot_len = strlen(ProtSeq);

  // No idea what the right number would be for this -- too small
  // and we might as well fill in all the diagonals (necessarily
  // a bad idea? -- prot_len * orf_len...), but too big and we
  // might miss something we don't really want to miss...
  int blocksize = BLOCK_SIZE;
  if (micro_exon)
    blocksize = 2;
  else if (blocksize > orf_len)
    blocksize = orf_len;
  
  int num_blocks = orf_len / blocksize; // OK to risk being one short

  // Rather than having a chart of blocks to look at, we'll highlight
  // where there are positive diagonals (since this should be a *sparse*
  // matrix).
  //
  // The value we store is equal to (prot_len*block_num + i)
  //
  int bp_size = prot_len;
  int * BlockPositives = NULL;
  BlockPositives = malloc(bp_size * sizeof(int));
  int num_positives = 0;

  // Now we'll check, treating each protein position as a candidate
  // start position, whether we can get a mapping to this block
  int block_num = 0;
  while (block_num < num_blocks) {

    char Block[blocksize];
    for (i=0; i<blocksize; i++)
      Block[i] = ORF[block_num*blocksize+i];

    // Run through the length of the protein, asking the perennial question
    // "you good?"
    for (i=0; i<prot_len-blocksize; i++) {

      j=0;
      while (j<blocksize && Block[j]==ProtSeq[i+j])
	j++;

      if (j==blocksize) {

	// If we need to realloc, do that right quick
	if (num_positives == bp_size-1) {
	  bp_size *= 2;
	  BlockPositives = realloc(BlockPositives,bp_size*sizeof(int));
	}
	BlockPositives[num_positives++] = prot_len * block_num + i;

      }
      
      // TODO
      // Realloc BlockPositives, or jump ship on account of repetitiveness
      // now...

    }

    block_num++;
    
  }

  // Who would've guessed all those blocks would take so long to fill out?
  // Oh well, now that we've burned a thousand years let's do some more
  // work!
  //
  // And by "do some more work" I mean
  //
  // Note that we can precompute what the next diagonal would be if it also
  // hit, and that we're working with ordered values -- use binary search to
  // make this a bit more efficient.
  //
  int AlreadyUsed[num_positives];
  for (i=0; i<num_positives; i++)
    AlreadyUsed[i] = 0;

  // Let's zag!
  for (i=0; i<num_positives; i++) {
    
    // If we've already looked at this diagonal, we're over it
    if (AlreadyUsed[i]) continue;
    AlreadyUsed[i] = 1;
    
    // Reverse engineer the block_num and starting position
    block_num = BlockPositives[i] / prot_len;
    int prot_start_index = BlockPositives[i] - block_num * prot_len;
    
    // Let's see how far this diagonal goes!
    int next_block_num   = block_num + 1;
    int next_start_index = prot_start_index + blocksize;
    int next_positive    = prot_len * next_block_num + next_start_index;
    int search_high = num_positives-1;
    int search_low  = i+1;

    while (search_low <= search_high && next_block_num < num_blocks) {

      int search_bin = (search_high + search_low) / 2;
	
      if (BlockPositives[search_bin] == next_positive) {
	
	// THE DREAM LIVES ON!
	AlreadyUsed[search_bin] = 1;
	next_block_num++;
	next_start_index += blocksize;
	next_positive = prot_len * next_block_num + next_start_index;
	search_high = num_positives-1;
	search_low = search_bin+1;

      } else if (BlockPositives[search_bin] < next_positive) {
	
	search_low = search_bin+1;
	
      } else {
	
	search_high = search_bin-1;
	
      }
    }

    // We now have the full run of this diagonal!
    // That is, the full run of this diagonal according to the BLOCKS -- but
    // let's see if we can extend a bit further!
    int prot_end_index  = next_start_index - 1;
    int orf_start_index = block_num * blocksize;
    int orf_end_index   = next_block_num * blocksize - 1;

    // We'll track the score of the diagonal as we run extend it.
    //
    // While the current implementation of this could theoretically
    // be done after we have the full extension, since this might
    // give us an X-drop-y alternative way to decide when to stop
    // extending, I'll build the basic infrastructure of tracking
    // score in that way here.
    float score = 0.0;
    for (j=0; j<=prot_end_index-prot_start_index; j++)
      score += MBB_AminoAliScore(ProtSeq[prot_start_index+j],ORF[orf_start_index+j]);

    int num_mismatches = 0;
    while (prot_start_index && orf_start_index) {

      // Check for mismatches (possibly jumping ship)
      if (ProtSeq[prot_start_index-1] != ORF[orf_start_index-1]) {
	num_mismatches++;
	if (num_mismatches > MAX_MISMATCH)
	  break;
      }

      // I'm still standing!
      score += MBB_AminoAliScore(ProtSeq[--prot_start_index],ORF[--orf_start_index]);

    }
    
    num_mismatches = 0;
    while (prot_end_index+1 < prot_len && orf_end_index+1 < orf_len) {

      // Check for mismatches (possibly jumping ship)
      if (ProtSeq[prot_end_index+1] != ORF[orf_end_index+1]) {
	num_mismatches++;
	if (num_mismatches > MAX_MISMATCH)
	  break;
      }

      // Boost that stinky score!
      score += MBB_AminoAliScore(ProtSeq[++prot_end_index],ORF[++orf_end_index]);

    }
    
    // Time to shout 'n' scream about that mapping!
    int map_len = prot_end_index - prot_start_index + 1;
    PrintMap(ProtSeq,prot_id,prot_start_index,exon_id,rf_index,ORF,orf_start_index,
	     map_len,NuclSeq,nucl_seq_len,rel_start_nucl,chr_start_nucl,score,revcomp);
    
    
  }
  
  if (BlockPositives) free(BlockPositives);
  
}




/*
 *  Function: OriginalFastMap
 *
 *  About: Scan a protein for mapping regions to an ORF
 *
 // NOTE: Right now this is somewhat simplistic -- we're assuming that use
 //       of an ARF would start right after a stop codon.
 *
 *
 */
void OriginalFastMap
( char * ProtSeq,
  int    prot_id,
  int    exon_id,
  int    rf_index,
  char * ORF,
  int    orf_len,
  char * NuclSeq,
  int    nucl_seq_len,
  int    rel_start_nucl,
  int    chr_start_nucl,
  int    revcomp
  ){

  int prot_len = strlen(ProtSeq);

  // The maximimum number of acceptable misses
  int max_misses = 1;

  // If we hit the maximum number of misses, but have a diagonal of this length
  // then we'll override and output the mapping anyways (set to prot_len to disable)
  int output_threshold = 10;
  output_threshold--; // Because we start at 'map_len = 0'

  // We'll require 5 aminos to a valid exon, so my condolences to the last 4 positions
  int num_live_diags = prot_len - 4;
  int * ProtStartPos = malloc(num_live_diags * sizeof(int));
  int * NumMisses    = malloc(num_live_diags * sizeof(int));

  // Because we allow 'max_misses' = 0, this is just an initialization (not a priming)
  int i;
  for (i=0; i<num_live_diags; i++) {
    ProtStartPos[i] = i;
    NumMisses[i] = 0;
  }

  // We'll keep track of where we are in nucleotides
  // UPDATE: This can be figured out when we need it later.
  //
  //int rel_end_nucl = rel_start_nucl + 3;
  //int chr_end_nucl = chr_start_nucl;
  //if (revcomp) chr_end_nucl -= 3;
  //else         chr_end_nucl += 3;
  
  // Run that scan
  //
  // Because this is currently just around for posterity, I'm being lazy and setting
  // the score output for any hit to 0.  If I were ever to make this an option (why?)
  // then I'd want to change that, but (why?)
  //
  int map_len  = 0;
  while (map_len < orf_len && num_live_diags) {

    // Who lives? WHO DIES?!
    //
    // Because we may be front-filling as diagonals die off,
    // we'll use a 'while' loop rather than a 'for' loop
    //
    i=0;
    while (i<num_live_diags) {

      // First off, is this diagonal finished?
      int prot_index = ProtStartPos[i]+map_len;
      if (prot_index == prot_len) {

	// It's the end of the road!
	PrintMap(ProtSeq,prot_id,ProtStartPos[i],exon_id,rf_index,ORF,0,map_len,
		 NuclSeq,nucl_seq_len,rel_start_nucl,chr_start_nucl,0.0,revcomp);

	// Front-fill, decrement num_live_diags
	num_live_diags--;
	ProtStartPos[i] = ProtStartPos[num_live_diags];
	NumMisses[i] = NumMisses[num_live_diags];
	
      } else {

	// This diagonal could continue, but is it good enough?
	if (ProtSeq[prot_index] != ORF[map_len]) NumMisses[i]++;

	// If the number of misses exceeds our maximum number of misses,
	// we'll kill this diagonal.
	if (NumMisses[i] > max_misses) {

	  // If the diagonal was behaving well up until this point, cut it some slack!
	  if (map_len >= output_threshold) {
	    PrintMap(ProtSeq,prot_id,ProtStartPos[i],exon_id,rf_index,ORF,0,map_len,
		     NuclSeq,nucl_seq_len,rel_start_nucl,chr_start_nucl,0.0,revcomp);
	  }

	  // Front-fill, decrement num_live_diags
	  num_live_diags--;
	  ProtStartPos[i] = ProtStartPos[num_live_diags];
	  NumMisses[i] = NumMisses[num_live_diags];
	
	} else {

	  // You live to diagonal another day!
	  i++;
	  
	}
	
      }
      
    }

    // Every amino we step forward in the mapping process is 3 nucleotides, duh
    map_len++;
    //rel_end_nucl += 3;
    //if (revcomp) chr_end_nucl -= 3;
    //else         chr_end_nucl += 3;
    
  }

  // Finally, if we still have any live diagonals after exhausting our ORF, you'd
  // better believe we want to hear about it!
  for (i=0; i<num_live_diags; i++) {
    int prot_index = ProtStartPos[i]+map_len;
    PrintMap(ProtSeq,prot_id,ProtStartPos[i],exon_id,rf_index,ORF,0,map_len,
	     NuclSeq,nucl_seq_len,rel_start_nucl,chr_start_nucl,0.0,revcomp);
  }
  
}



/*
 *  M A I N
 *
 */
int main (int argc, char ** argv) {

  // Check usage
  if (argc < 6) return PrintUsage();

  // Parse the nucleotide file
  // Note that in the case of reverse strand the start_index will be larger
  int nucl_seq_len;
  int nucl_start_index;
  int nucl_end_index;
  char * NuclSeq = ParseNucleotideFile(argv[3],&nucl_seq_len,&nucl_start_index,&nucl_end_index);
  
  // Make sure things went well 
  if (NuclSeq == NULL) {
    fprintf(stderr,"\n  ERROR:  Failed to parse nucleotide file '%s'\n\n",argv[2]);
    return 1;
  }

  // Parse the protein file
  int num_seqs = atoi(argv[2]);
  char ** ProtSeqs = ParseProteinFile(argv[1],num_seqs);

  // Make sure things went not poorly
  if (ProtSeqs == NULL) {
    fprintf(stderr,"\n  ERROR:  Failed to parse protein file '%s'\n\n",argv[1]);
    return 1;
  }

  // Alright then!  If we made it this far, then it would seem the game is afoot!
  int i,j,k;

  // Our first task is just to process the exon locations so that they index into
  // the right places in our nucleotide sequence array.
  int num_exons    = (argc - 3) / 2;
  int * ExonStarts = malloc(num_exons*sizeof(int));  
  int * ExonEnds   = malloc(num_exons*sizeof(int));

  // We'll start off by just pulling in the coordinates
  j=4;
  int revcomp = 0;
  if (nucl_start_index > nucl_end_index) {
    revcomp = 1;
    for (i=0; i<num_exons; i++) {
      ExonStarts[i] = nucl_start_index - atoi(argv[j++]);
      ExonEnds[i]   = nucl_start_index - atoi(argv[j++]);
    }
  } else {
    for (i=0; i<num_exons; i++) {
      ExonStarts[i] = atoi(argv[j++]) - nucl_start_index;
      ExonEnds[i]   = atoi(argv[j++]) - nucl_start_index;
    }
  }
  // We should now be able to index directly in / out of exons

  // Now we'll expand out a couple codons, for the hallibut!
  // NOTE that we don't expand if the GTF is telling us that this is a smallish
  // exon -- in this case we do a much narrower search, which will be more expensive
  // but will help us pick up smaller bits
  int exon_ext_len = 6 * BLOCK_SIZE;
  int micro_thresh = 3 * (3 * BLOCK_SIZE / 2); // 150% of block size
  for (i=0; i<num_exons; i++) {
    if (ExonEnds[i] - ExonStarts[i] >= micro_thresh) {
      ExonStarts[i] -= exon_ext_len;
      ExonEnds[i]   += exon_ext_len;
    } else {
      ExonStarts[i] -= 2;
      ExonEnds[i]   += 2;
    }
  }
  
  // Moreover, since things are nice 'n' easy, let's figure out
  // how much space we'll need to store our ORFs.
  // Note that we're going to be one short of the actual length,
  // but that doesn't matter too much...
  int max_exon_len = 0;
  for (i=0; i<num_exons; i++) {
    if (ExonEnds[i] - ExonStarts[i] > max_exon_len)
      max_exon_len = ExonEnds[i] - ExonStarts[i];
  }
  char * RF = malloc((2+max_exon_len/3)*sizeof(char));
  
  // Now for the fun part!  Let's do this as an exon-by-exon thing,
  // where we translate out the 3 ORFs and do an all-vs-all against
  // the proteins.
  for (i=0; i<num_exons; i++) {

    // Where is this actually situated on the genome?
    int chr_start_nucl = nucl_start_index;
    if (revcomp) chr_start_nucl -= ExonStarts[i];
    else         chr_start_nucl += ExonStarts[i];
    
    // Is this a short exon?  This will impact the blocksize we use...
    int micro_exon = 0;
    if (abs(ExonStarts[i]-ExonEnds[i])-4 < micro_thresh)
      micro_exon = 1;

    // Iterate through the reading frames
    int rf_index;
    for (rf_index=0; rf_index<3; rf_index++) {

      // Because there could be stop codons sprinkled through this
      // frame, we'll distinguish what we're getting out of this
      // function as an RF rather than as an ORF.
      TranslateRF(NuclSeq,ExonStarts[i]+rf_index,ExonEnds[i]+rf_index,RF);
      int rf_len = strlen(RF);

      // We'll scan through the frame and pull out each open reading
      // frame (w/ at least 5 aminos)
      int min_orf_aminos  = 5;
      int orf_start_index = 0;
      int orf_end_index;
      while (orf_start_index < rf_len) {

	// In case we start with an 'X' we'll want to scan ahead
	while (orf_start_index < rf_len && RF[orf_start_index] == 'X')
	  orf_start_index++;
	if (orf_start_index == rf_len)
	  break;

	// Scan to the end of the reading frame or the next stop codon,
	// filling out our ORF.
	// Note that we'll end our looping at the value that doesn't work,
	// which makes calculating 'len' easier, but requires a quick
	// decrement.
	orf_end_index = orf_start_index+1;
	while (orf_end_index < rf_len && RF[orf_end_index] != 'X')
	  orf_end_index++;
	int orf_len = orf_end_index - orf_start_index;
	orf_end_index--;

	// If our ORF is long enough to be worth checking, run through our
	// protein sequences and see if we have any mappings worth remarking on
	if (orf_len >= min_orf_aminos || (micro_exon && orf_len >= 2)) {

	  // For reporting we'll need to know some more detailed intel about the
	  // ORF's location, both relative to the sequence we have and the whole
	  // dang chromosome.
	  //
	  // Keep in mind that these are inherently oriented w.r.t. codon starts
	  // rather than centers
	  //
	  int orf_offset   = rf_index + 3*orf_start_index;
	  int orf_rel_nucl = ExonStarts[i] + orf_offset;
	  int orf_chr_nucl = chr_start_nucl;
	  if (revcomp) orf_chr_nucl -= orf_offset;
	  else         orf_chr_nucl += orf_offset;;

	  // Two methods for doing our mapping
	  //
	  // [1]
	  //
	  // 'BlockScan' fills out a sparse matrix of k-mer matches, and tries
	  // to stitch those together and then extend them out.
	  //
	  // 'BlockScan' will be more computationally expensive and is (as
	  // it's currently implemented) intolerant of mismatches, but lets us
	  // find internal parts of exons that match to our proteins.
	  //
	  for (j=0; j<num_seqs; j++)
	    BlockScan(ProtSeqs[j],j,i,rf_index,&RF[orf_start_index],orf_len,micro_exon,
		      NuclSeq,nucl_seq_len,orf_rel_nucl,orf_chr_nucl,revcomp);

	  // [2]
	  //
	  // 'OriginalFastMap' runs through the open reading frame start to finish,
	  // looking for hits with a limited number of mismatches.
	  //
	  // The primary shortcoming with 'OriginalFastMap' is that it won't be
	  // able to capture ARFs very well, since they'll often be embedded in
	  // the middle of an exon, and thus won't be discovered when we
	  // start at the start.
	  //
	  //for (j=0; j<num_seqs; j++)
	  //OriginalFastMap(ProtSeqs[j],j,i,rf_index,&RF[orf_start_index],orf_len,
	  //NuclSeq,orf_rel_nucl,orf_chr_nucl,revcomp);
	  
	}
	
	// Prep for next round
	orf_start_index = orf_end_index+1;
	
      }

      // Iterate through the protein sequences
      //for (j=0; j<num_seqs; j++)
      //FindMaps(ProtSeqs[j],ORF,chr_start_nucl);

    }
    
  }

  // Release the sequences!
  free(NuclSeq);
  for (i=0; i<num_seqs; i++)
    free(ProtSeqs[i]);
  free(ProtSeqs);
  free(RF);

  // These things are also irrelevant now!
  free(ExonStarts);
  free(ExonEnds);

  return 0;
  
}
