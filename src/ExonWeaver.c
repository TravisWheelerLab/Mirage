#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "BasicBio.h"


// FUNCTION LIST
// -------------
//  + BarfNodeInfo
//  + IncreaseInEdgeCapacity
//  + IncreaseOutEdgeCapacity
//  + DestroyGraph
//  + SpliceProbToScore
//  + ReadNodesFromFile
//  + UpdateHitExtScore
//  + AttemptConnection
//  + SortFloats
//  + OrganizeEdgeList
//  + ConnectGraph
//  + RecursivePathEval
//  + GetHitExonChrStart
//  + GetHitExonChrEnd
//  + PrintHitMeta
//  + PrintExonMeta
//  + PrintSplicedAlignment
//  + ReportSplicing
//  + ReportMaximalPaths
//  + ParseCommandArgs
//  + PrintUsage



// NOTE: I'm currently assuming that the input nucleotide sequence
//       provides 17 nucleotides on either side of the 'coding' sequence,
//       which has some impacts on assumptions I make about coding in
//       how far we're willing to consider extending a hit, and stuff
//       like that.
//
//       Because this might not always be the case, I'm going to add
//       a comment reading '*17*' wherever this assumption is being
//       explicitly applied.


// There's a small automatic penalty for splicing, so that everything else
// equal we'll prefer alignments with fewer exons.
// There's also a maximum score for looking like a great 3' or 5' splice site!
const float SPLICE_COST      = -6.0;
const float MAX_SPLICE_SCORE =  5.0;


/*
 *  Struct: HW_NODE
 *
 */
typedef struct _HW_NODE_T_ {

  // The data that we pull directly from the file
  int start_amino;
  int end_amino;
  int start_nucl;
  int end_nucl;
  int nucl_str_len; // Includes buffer nucleotides
  float score;
  char  * Nucls; // All uppercase, with 17 buffer nucleotides on each side of coding region
  float * TPSS;  // 3' SS score (prob. of being last  nucl. in intron)
  float * FPSS;  // 5' SS score (prob. of being first nucl. in intron)

  struct _HW_NODE_T_ ** Incoming;
  int   * IncomingID;
  float * InEdgeScore;       // What does this add to the scores of the nodes?
  int   * FirstCodingNucl;   // What position in 'Nucls' follows the splice signal?
  int num_incoming;
  int max_incoming; // For reallocing

  struct _HW_NODE_T_ ** Outgoing;
  int   * OutgoingID;
  float * OutEdgeScore;     // What does this add to the scores of the nodes?
  int   * LastCodingNucl;   // What position in 'Nucls' precedes the splice signal?
  float * CumScoreMemo;     // Cumulative score memoization (avoid re-treading graph)
  int num_outgoing;
  int max_outgoing; // For reallocing

} HW_NODE;




/*
 *  Struct: HW_OPTS
 *
 */
typedef struct _HW_OPTS_T_ {
  int req_consistency; // Do we require hits follow a logical order on the chromosome?
  int output_format;   // Complete alignments (1), or exon-by-exon (2)?
  int report_singles;  // Do hits need to have more than one exon?
} HW_OPTS;




/*
 *  DEBUGGING: BarfNodeInfo
 *
 */
void BarfNodeInfo (HW_NODE * N, int node_id) {
  int i;
  printf("\n");
  printf("NODE %d\n",node_id);
  printf("  Node Score : %f\n",N->score);
  printf("  Amino Range: %d..%d\n",N->start_amino,N->end_amino);
  printf("  Nucl. Range: %d..%d\n",N->start_nucl,N->end_nucl);
  printf("  Nucleotides: %s\n",N->Nucls);
  printf("  In  Edges  : %d",N->num_incoming);
  if (N->num_incoming) {
    printf("  (%d:%f",N->IncomingID[0],N->InEdgeScore[0]);
    for (i=1; i<N->num_incoming; i++)
      printf("/%d:%f",N->IncomingID[i],N->InEdgeScore[i]);
    printf(")");
  }
  printf("\n");
  printf("  Out Edges  : %d",N->num_outgoing);
  if (N->num_outgoing) {
    printf("  (%d:%f",N->OutgoingID[0],N->OutEdgeScore[0]);
    for (i=1; i<N->num_outgoing; i++)
      printf("/%d:%f",N->OutgoingID[i],N->OutEdgeScore[i]);
    printf(")");
  }
  printf("\n");
  printf("\n");
}





/*
 *  Functions: Increase[In/Out]EdgeCapacity
 *
 */
void IncreaseInEdgeCapacity (HW_NODE * N) {

  // We'll start off by creating some temporary arrays, to guarantee everything's
  // copied correctly
  HW_NODE ** TempNode = malloc(N->num_incoming * sizeof(HW_NODE *));
  int TempIID[N->num_incoming];
  int TempFCN[N->num_incoming];
  float TempIES[N->num_incoming];

  int i;
  for (i=0; i<N->num_incoming; i++) {
    TempNode[i] = N->Incoming[i];
    TempIID[i]  = N->IncomingID[i];
    TempFCN[i]  = N->FirstCodingNucl[i];
    TempIES[i]  = N->InEdgeScore[i];
  }

  N->max_incoming *= 2;

  N->Incoming = realloc(N->Incoming, N->max_incoming * sizeof(HW_NODE *));
  N->IncomingID = realloc(N->IncomingID, N->max_incoming * sizeof(int));
  N->FirstCodingNucl = realloc(N->FirstCodingNucl, N->max_incoming * sizeof(int));
  N->InEdgeScore = realloc(N->InEdgeScore, N->max_incoming * sizeof(float));
  
  for (i=0; i<N->num_incoming; i++) {
    N->Incoming[i] = TempNode[i];
    N->IncomingID[i] = TempIID[i];
    N->FirstCodingNucl[i] = TempFCN[i];
    N->InEdgeScore[i] = TempIES[i];
  }

  free(TempNode);

}
//
//
void IncreaseOutEdgeCapacity (HW_NODE * N) {

  HW_NODE ** TempNode = malloc(N->num_outgoing * sizeof(HW_NODE *));
  int TempOID[N->num_outgoing];
  int TempLCN[N->num_outgoing];
  float TempOES[N->num_outgoing];

  int i;
  for (i=0; i<N->num_outgoing; i++) {
    TempNode[i] = N->Outgoing[i];
    TempOID[i]  = N->OutgoingID[i];
    TempLCN[i]  = N->LastCodingNucl[i];
    TempOES[i]  = N->OutEdgeScore[i];
  }

  N->max_outgoing *= 2;

  N->Outgoing = realloc(N->Outgoing, N->max_outgoing * sizeof(HW_NODE *));
  N->OutgoingID = realloc(N->OutgoingID, N->max_outgoing * sizeof(int));
  N->LastCodingNucl = realloc(N->LastCodingNucl, N->max_outgoing * sizeof(int));
  N->OutEdgeScore = realloc(N->OutEdgeScore, N->max_outgoing * sizeof(float));
  N->CumScoreMemo = realloc(N->CumScoreMemo, N->max_outgoing * sizeof(float));
  
  for (i=0; i<N->num_outgoing; i++) {
    N->Outgoing[i] = TempNode[i];
    N->OutgoingID[i] = TempOID[i];
    N->LastCodingNucl[i] = TempLCN[i];
    N->OutEdgeScore[i] = TempOES[i];
    N->CumScoreMemo[i] = MBB_NINF;
  }

  free(TempNode);
  
}




/*
 *  Function: DestroyGraph
 *
 */
void DestroyGraph (HW_NODE ** Graph, int num_exons) {
  HW_NODE * N;
  int i;
  for (i=0; i<num_exons; i++) {
    N = Graph[i];
    free(N->Nucls);
    free(N->TPSS);
    free(N->FPSS);
    free(N->Incoming);
    free(N->IncomingID);
    free(N->InEdgeScore);
    free(N->FirstCodingNucl);
    free(N->Outgoing);
    free(N->OutgoingID);
    free(N->OutEdgeScore);
    free(N->LastCodingNucl);
    free(N->CumScoreMemo);
    free(N);
  }
  free(Graph);
}



/*
 *  Function: SpliceProbToScore
 *
 */
float SpliceProbToScore (float prob) {
  // For now, we'll do a simple quadratic function, with a max score of 10
  return MAX_SPLICE_SCORE * (prob * prob);
}




/*
 *  Function: ReadNodesFromFile
 *
 */
void ReadNodesFromFile (HW_NODE ** Graph, int num_exons, FILE * inf) {

  // We already know how much reading we have assigned!
  int i,j;
  for (i=0; i<num_exons; i++) {
    
    // Start off by allocating space for this node
    Graph[i] = (HW_NODE *)malloc(sizeof(HW_NODE));
    HW_NODE * Node = Graph[i];

    fscanf(inf,"\n"); // There should be an empty line...
    fscanf(inf,"Hit Score   : %f\n",&(Node->score));
    fscanf(inf,"Amino Range : %d..%d\n",&(Node->start_amino),&(Node->end_amino));
    fscanf(inf,"Nucl Range  : %d..%d\n",&(Node->start_nucl) ,&(Node->end_nucl));

    // We need to do a bit of computation to make sure we've nailed down the
    // actual length of our nucleotide string: (diff + 1) + 17 left + 17 right
    //
    // *17*
    //
    Node->nucl_str_len = abs(Node->end_nucl - Node->start_nucl) + 35;
    Node->Nucls = malloc((Node->nucl_str_len+1)*sizeof(char));
    Node->TPSS  = malloc((Node->nucl_str_len+1)*sizeof(float));
    Node->FPSS  = malloc((Node->nucl_str_len+1)*sizeof(float));

    // Reading in the nucleotide stuff is gonna prove my mastery of fscanf...
    fscanf(inf,"Nucleotides : %s\n",Node->Nucls);
    fscanf(inf,"3' SS Str.  :");
    float splice_prob;
    for (j=0; j<Node->nucl_str_len; j++) {
      fscanf(inf," %f",&splice_prob);
      Node->TPSS[j] = SpliceProbToScore(splice_prob);
    }
    fscanf(inf,"\n5' SS Str.  :");
    for (j=0; j<Node->nucl_str_len; j++) {
      fscanf(inf," %f",&splice_prob);
      Node->FPSS[j] = SpliceProbToScore(splice_prob);
    }
    fscanf(inf,"\n");

    // Finish up this node by initializing edge-connection stuff
    int init_max_in    = 4;
    Node->num_incoming = 0;
    Node->max_incoming = init_max_in;
    Node->Incoming        = malloc(init_max_in * sizeof(HW_NODE *));
    Node->IncomingID      = malloc(init_max_in * sizeof(int));
    Node->InEdgeScore     = malloc(init_max_in * sizeof(float));
    Node->FirstCodingNucl = malloc(init_max_in * sizeof(int));

    int init_max_out   = 4;
    Node->num_outgoing = 0;
    Node->max_outgoing = init_max_out;
    Node->Outgoing       = malloc(init_max_out * sizeof(HW_NODE *));
    Node->OutgoingID     = malloc(init_max_out * sizeof(int));
    Node->OutEdgeScore   = malloc(init_max_out * sizeof(float));
    Node->LastCodingNucl = malloc(init_max_out * sizeof(int));
    Node->CumScoreMemo   = malloc(init_max_out * sizeof(float));

  }
  

}





/*
 *  Function: UpdateHitExtScore
 *
 */
void UpdateHitExtScore
(
 char  * TransSeq,
 char  * RefSeq,
 float * AliScores,
 int     array_index,
 int     array_len,
 char    next_amino,
 int   * num_stop_codons,
 float * sum_score
 ){

  // Whatever we do, the next amino is the new amino.
  char trans_char = TransSeq[array_index];
  TransSeq[array_index] = next_amino;

  // It'd be super cool if we were fixing an existing stop codon...
  if (trans_char == 'X') {

    if (next_amino != 'X') {

      // YES! At least one stop codon patched over!
      *num_stop_codons -= 1;

      // If we still have stop codons, the score is still -inf
      if (*num_stop_codons) return;

      // OH BOOOIIIIIII, we're finally free from the tyranny of stop codons!
      *sum_score = 0.0;
      int i;
      for (i=0; i<array_len; i++) {
	AliScores[i] = MBB_AminoAliScore(TransSeq[i],RefSeq[i]);
	*sum_score  += AliScores[i];
      }

    }
      
  } else if (next_amino == 'X') {

    *num_stop_codons += 1;
    *sum_score = MBB_NINF;
    
  } else if (*num_stop_codons == 0) {

    *sum_score -= AliScores[array_index];
    AliScores[array_index]
      = MBB_AminoAliScore(TransSeq[array_index],RefSeq[array_index]);
    *sum_score += AliScores[array_index];

  }

}





/*
 *  Function: AttemptConnection
 *
 */
void AttemptConnection
(
 HW_NODE ** Graph,
 int left_index,
 int right_index,
 char * Seq
 ){

  int i,j,k;
  
  HW_NODE *  LeftNode = Graph[left_index];
  HW_NODE * RightNode = Graph[right_index];

  int left_end_amino = LeftNode->end_amino;
  int right_start_amino = RightNode->start_amino;

  // To avoid the headache of indirection, we're going to actually
  // extract the potentially overlapping sequences into two arrays
  // that we can walk along with a single index.
  //
  // Keep in mind that there are 17 buffer nucleotides on either
  // side of each "exon," of which 15 (5 aminos) can be extended
  // into. 

  // How far do we need to unravel into each on account of overlap?
  // Note that this will work fine if we're extending into the buffer
  // zone, rather than unraveling the original exons.
  //
  int overlap_len = left_end_amino - right_start_amino + 1;

  // So how many coding nucleotides are overlapping?
  int nucl_overlap = overlap_len * 3;

  // If we have negative overlap (i.e., we need to extend into buffer)
  // then we'll want to have that flagged
  int overlapped = 1;
  if (overlap_len < 0) {
    overlapped    =  0;
    overlap_len  *= -1;
    nucl_overlap *= -1;
  }

  // Next we'll define some boundaries for our search space.
  //
  // NOTE that just to guarantee that we get a good splice signal check
  // in all cases (esp. overlap_len==0), we're going to extend one
  // additional codon into the buffer zone.
  //
  overlap_len  += 1;
  nucl_overlap += 3;
  //
  // *17* IN A MAJOR WAY
  //
  int lx_lb_amino, lx_rb_amino, rx_lb_amino, rx_rb_amino;
  int lx_lb_nucl , lx_rb_nucl , rx_lb_nucl , rx_rb_nucl;
  if (overlapped) {

    lx_rb_amino = LeftNode->end_amino + 1;
    lx_lb_amino = lx_rb_amino - overlap_len;

    rx_lb_amino = RightNode->start_amino - 1;
    rx_rb_amino = rx_lb_amino + overlap_len;

    lx_rb_nucl  = LeftNode->nucl_str_len - 15;   // Last  nucl of codon
    lx_lb_nucl  = lx_rb_nucl - 2 - nucl_overlap; // First nucl of codon

    rx_lb_nucl  = 14;                            // First nucl of codon
    rx_rb_nucl  = rx_lb_nucl + 2 + nucl_overlap; // Last  nucl of codon
    
  } else {

    lx_lb_amino = LeftNode->end_amino - 1;
    lx_rb_amino = lx_lb_amino + overlap_len;

    rx_rb_amino = RightNode->start_amino + 1;
    rx_lb_amino = rx_rb_amino - overlap_len;

    lx_lb_nucl  = LeftNode->nucl_str_len - 20;   // First nucl of codon
    lx_rb_nucl  = lx_lb_nucl + 2 + nucl_overlap; // Last  nucl of codon

    rx_rb_nucl  = 19;                            // Last  nucl of codon
    rx_lb_nucl  = rx_rb_nucl - 2 - nucl_overlap; // First nucl of codon

  }

  // Pull out our reference aminos
  char RefSeq[overlap_len];
  for (i=0; i<overlap_len; i++)
    RefSeq[i] = Seq[lx_lb_amino+i];

  // Build up the extracted nucleotide-y data
  //
  // NOTE: We're going to pull in the splicing data so that it's representing the
  //       probability of a nucleotide being the last in the exon, rather than first
  //       in the intron.
  //
  char   LeftNucls[nucl_overlap];
  char  RightNucls[nucl_overlap];
  float   LeftFPSS[nucl_overlap];
  float  RightTPSS[nucl_overlap];
  for (i=0; i<nucl_overlap; i++) {
    LeftNucls[i]  = LeftNode->Nucls[lx_lb_nucl + i];
    LeftFPSS[i]   = LeftNode->FPSS[lx_lb_nucl + i + 1];
    RightNucls[i] = RightNode->Nucls[rx_lb_nucl + i];
    RightTPSS[i]  = RightNode->TPSS[rx_lb_nucl + i - 1];
  }


  // Compute the cost of unraveling the aminos at each overlapping position
  float UnravelCost[overlap_len];
  if (overlapped) {

    float LeftUnravelCost[overlap_len];
    float RightUnravelCost[overlap_len];
    float pos_score;

    // Compute the score of each position independently for the two exons
    for (i=0; i<overlap_len; i++) {

      pos_score = MBB_AminoAliScore(RefSeq[i],MBB_TranslateCodon(&LeftNucls[i*3]));
      if (pos_score != MBB_NINF) LeftUnravelCost[i] = 0.0 - pos_score;
      else                       LeftUnravelCost[i] = 0.0;

      pos_score = MBB_AminoAliScore(RefSeq[i],MBB_TranslateCodon(&RightNucls[i*3]));
      if (pos_score != MBB_NINF) RightUnravelCost[i] = 0.0 - pos_score;
      else                       RightUnravelCost[i] = 0.0;

    }

    // Accumulate cost for the left-side node, from right-to-left
    for (i=overlap_len-2; i>=0; i--)
      LeftUnravelCost[i] += LeftUnravelCost[i+1];

    // Accumulate cost for the right-side node, from left-to-right
    for (i=1; i<overlap_len; i++)
      RightUnravelCost[i] += RightUnravelCost[i-1];
    
    // Transfer over into a combined unraveling cost
    for (i=0; i<overlap_len; i++)
      UnravelCost[i] = LeftUnravelCost[i] + RightUnravelCost[i];
    
  } else {
    for (i=0; i<overlap_len; i++)
      UnravelCost[i] = 0.0;
  }
  
  

  // Finally, initialize our translated sequence space.
  //
  // We'll start the case where a single nucleotide is contributed by
  // the left exon, so that the work of our loop is just scanning until
  // we only have a single nucleotide contributed by the right exon.
  //
  // To avoid having the addition of negative infinities break
  // everything forever, we'll keep count of how many stop codons
  // are in the translated sequence, and have a trigger to check
  // 'sum_score' only when the counter's at zero.
  //
  char  TransSeq[overlap_len];
  float AminoAliScores[overlap_len];
  float sum_score = 0.0;
  int   num_stop_codons = 0;
  for (i=1; i<overlap_len; i++) {
    TransSeq[i] = MBB_TranslateCodon(&RightNucls[i*3]);
    if (TransSeq[i] == 'X') {
      num_stop_codons++;
    } else {
      AminoAliScores[i] = MBB_AminoAliScore(RefSeq[i],TransSeq[i]);
      sum_score += AminoAliScores[i];
    }
  }

  // Load up the first splicing
  char Codon[3];
  Codon[0] = LeftNucls[0];
  Codon[1] = RightNucls[1];
  Codon[2] = RightNucls[2];
  TransSeq[0] = MBB_TranslateCodon(Codon);
  AminoAliScores[0] = MBB_AminoAliScore(RefSeq[0],TransSeq[0]);
  if (TransSeq[0] == 'X')
    num_stop_codons++;

  // Could still be -inf, but whatever...
  sum_score += AminoAliScores[0];

  // If we started out with any stop codons, our sum score is -inf
  if (num_stop_codons) sum_score = MBB_NINF;

  // Track the highest score and what positions defined the ends of
  // our coding regions.  Currently, these will always be one off from
  // one another, but if we ever add gappery these could split apart...
  int top_left_break  = 0;
  int top_right_break = 1;
  float top_score = sum_score + LeftFPSS[0] + RightTPSS[1] + UnravelCost[0];

  // i tracks the position that will be the first drawing from
  //   the left node's nucleotide sequence,
  // j tracks the position in TransSeq that we write to
  // k is equal to 'i % 3' (avoids running mod)
  j = 0;
  k = 1;
  for (i=1; i<nucl_overlap-1; i++) {

    // Exactly what we're switching over (from Right to Left)
    // will depend on how far into the current codon we are
    if (k==0) {
      Codon[0] = LeftNucls[i];
      Codon[1] = RightNucls[i+1];
      Codon[2] = RightNucls[i+2];
      j++;
      k=1;
    } else if (k==1) {
      Codon[1] = LeftNucls[i];
      k=2;
    } else { // if (k==2)
      Codon[2] = LeftNucls[i];
      k=0;
    }

    // Well, that wasn't toooooo nasty!  Let's run the translation and
    // compare scores.
    //
    // Just because this code ends up looking a little ugly, I'm offloading
    // it to a different function.
    //
    char next_amino = MBB_TranslateCodon(Codon);
    UpdateHitExtScore(TransSeq,RefSeq,AminoAliScores,
		      j,overlap_len,next_amino,&num_stop_codons,&sum_score);

    // If we've beaten out the previous champion, then we'll seize the throne!
    if (num_stop_codons == 0) {
      float score_check = sum_score + LeftFPSS[i] + RightTPSS[i+1] + UnravelCost[j];
      if (score_check > top_score) {
	top_left_break  = i;
	top_right_break = i+1;
	top_score       = score_check;
      }
    }
    
  }

  // Now we *should* have a top_score (if it's negative, bail!), and we
  // can reverse engineer the actual positions in our sequences that we'd
  // use for splicing
  if (top_score == MBB_NINF) return;

  // WE'RE SPLICING! (so you gotta pay that lil' tax...)
  top_score += SPLICE_COST;

  // Before we record any position info., check if we need to bump up the
  // space for recording edges 
  if (LeftNode->num_outgoing == LeftNode->max_outgoing)
    IncreaseOutEdgeCapacity(LeftNode);

  if (RightNode->num_incoming == RightNode->max_incoming)
    IncreaseInEdgeCapacity(RightNode);

  // Soooooo, where were we?
  int l_out_edges = LeftNode->num_outgoing;
  int r_in_edges  = RightNode->num_incoming;
  
  LeftNode->num_outgoing += 1;
  LeftNode->Outgoing[l_out_edges] = RightNode;
  LeftNode->OutgoingID[l_out_edges] = right_index;
  LeftNode->LastCodingNucl[l_out_edges] = lx_lb_nucl + top_left_break;
  LeftNode->OutEdgeScore[l_out_edges] = top_score;
  LeftNode->CumScoreMemo[l_out_edges] = MBB_NINF;

  RightNode->num_incoming += 1;
  RightNode->Incoming[r_in_edges] = LeftNode;
  RightNode->IncomingID[r_in_edges] = left_index;
  RightNode->FirstCodingNucl[r_in_edges] = rx_lb_nucl + top_right_break;
  RightNode->InEdgeScore[r_in_edges] = top_score;

  // That's an edge!

}





/*
 *  Function: SortFloats
 *
 */
void SortFloats (float * Vals, int * Index, int num_vals) {

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
  while (num_blocks > 1) {

    for (i=0; i+blocksize<num_vals; i+=2*blocksize) {

      j     = i;
      j_lim = j+blocksize;
      k     = j_lim;
      k_lim = k+blocksize;
      if (k_lim>num_vals)
	k_lim = num_vals;
      
      while (j<j_lim && k<k_lim) {
	if (Vals[Read[j]]>Vals[Read[k]])
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
 *  Function: OrganizeEdgeList
 *
 */
void OrganizeEdgeList (HW_NODE * N) {

  int i;

  // We'll make an array of edge indices, and an equally-sized
  // temporary array
  int index_size = MBB_MaxInt(N->num_incoming,N->num_outgoing);
  if (!index_size)
    return;
  
  int      * Index      = malloc(index_size*sizeof(int));
  HW_NODE ** TempNodes  = malloc(index_size*sizeof(HW_NODE *));
  int      * TempInts   = malloc(index_size*sizeof(int));
  float    * TempFloats = malloc(index_size*sizeof(float *));

  // 1. INCOMING EDGES
  SortFloats(N->InEdgeScore,Index,N->num_incoming);

  for (i=0; i<N->num_incoming; i++) TempNodes[i] = N->Incoming[i];
  for (i=0; i<N->num_incoming; i++) N->Incoming[i] = TempNodes[Index[i]];

  for (i=0; i<N->num_incoming; i++) TempInts[i] = N->IncomingID[i];
  for (i=0; i<N->num_incoming; i++) N->IncomingID[i] = TempInts[Index[i]];

  for (i=0; i<N->num_incoming; i++) TempInts[i] = N->FirstCodingNucl[i];
  for (i=0; i<N->num_incoming; i++) N->FirstCodingNucl[i] = TempInts[Index[i]];

  for (i=0; i<N->num_incoming; i++) TempFloats[i] = N->InEdgeScore[i];
  for (i=0; i<N->num_incoming; i++) N->InEdgeScore[i] = TempFloats[Index[i]];

  // 2. OUTGOING EDGES  
  SortFloats(N->OutEdgeScore,Index,N->num_outgoing);

  for (i=0; i<N->num_outgoing; i++) TempNodes[i] = N->Outgoing[i];
  for (i=0; i<N->num_outgoing; i++) N->Outgoing[i] = TempNodes[Index[i]];

  for (i=0; i<N->num_outgoing; i++) TempInts[i] = N->OutgoingID[i];
  for (i=0; i<N->num_outgoing; i++) N->OutgoingID[i] = TempInts[Index[i]];

  for (i=0; i<N->num_outgoing; i++) TempInts[i] = N->LastCodingNucl[i];
  for (i=0; i<N->num_outgoing; i++) N->LastCodingNucl[i] = TempInts[Index[i]];

  for (i=0; i<N->num_outgoing; i++) TempFloats[i] = N->OutEdgeScore[i];
  for (i=0; i<N->num_outgoing; i++) N->OutEdgeScore[i] = TempFloats[Index[i]];

  // SWEET RELEASE
  free(Index);
  free(TempNodes);
  free(TempInts);
  free(TempFloats);

}





/*
 *  Function: ConnectGraph
 *
 */
void ConnectGraph (HW_NODE ** Graph, int num_exons, char * Seq, HW_OPTS * Opts) {

  int i,j;

  // If there's a short-ish gap between the two exons, we'll extend
  // into the putative intron to see if we can get a decent splicing.
  // We're currently assuming 17 nucleotides on either side for our
  // buffer, although I don't know if we want to go much further than
  // 4 aminos...
  int max_extension = 4; // *17* (sort of -- could go up to 10...)

  // We're going to assume that all of our hits are mappings to
  // the same chromosome, so we'll do a quick strand check
  int revcomp = 0;
  if (Graph[0]->start_nucl > Graph[0]->end_nucl)
    revcomp = 1;

  // At least for now, what we'll do is scan through our graph,
  // and at each node we'll run ahead to find the range of other
  // nodes that it could splice into (benefiting from the fact
  // that the nodes are pre-sorted by 'start_amino' values).
  for (i=0; i<num_exons-1; i++) {

    int left_end_amino = Graph[i]->end_amino;

    // Wowee! Let's see how many friends we can hook up with ;)
    for (j=i+1; j<num_exons; j++) {

      int right_start_amino = Graph[j]->start_amino;
      int overlap_len = left_end_amino - right_start_amino + 1;

      // We don't want to overextend ourselves, now
      if (overlap_len < -1 * max_extension)
	break;

      // Is there a sensible splice to be found?
      if (Graph[i]->start_amino == right_start_amino
	  || left_end_amino >= Graph[j]->end_amino)
	continue;
      
      // Remember to make sure these hits are compatibly located on the chromosome!
      // (but only if that's what you're into!)
      if (Opts->req_consistency &&
	  ((!revcomp && Graph[i]->end_nucl - (3*overlap_len) < Graph[j]->start_nucl) ||
	   ( revcomp && Graph[i]->end_nucl + (3*overlap_len) > Graph[j]->start_nucl)))
	AttemptConnection(Graph,i,j,Seq);
      
    }
    
  }

  // Now that we've got all our connections drawn in, let's
  // sort each node's edge list by score (descending, obv.s)
  for (i=0; i<num_exons; i++) OrganizeEdgeList(Graph[i]);

}



/*
 *  Function: ExtendTransMap
 *
 */
float ExtendTransMap
(
 HW_NODE * N,
 int start_nucl,
 int end_nucl,
 char * Seq,
 int * amino_depth,
 char * MapNucls,
 int * map_nucl_len
){

  float mapscore = 0.0;
  
  int i;
  for (i=start_nucl; i<=end_nucl; i++) {
    MapNucls[*map_nucl_len] = N->Nucls[i];
    *map_nucl_len += 1;
    if (*map_nucl_len % 3 == 0) {
      mapscore += MBB_AminoAliScore(Seq[*amino_depth],
				    MBB_TranslateCodon(&MapNucls[*map_nucl_len-3]));
      *amino_depth += 1;
    }
  }

  return mapscore;
  
}






/*
 *  Function: RecursivePathEval
 *
 */
float RecursivePathEval
(
 HW_NODE ** Graph,
 char * Seq,
 int node_id,
 int source_id,
 float * TopScoreThroughNode,
 int * TopScoreSourceNode,
 int * TopScoreTargetNode,
 int target_amino,
 char * MapNucls,
 int start_map_nucl_len,
 int start_amino_depth,
 float map_cum_score
){

  int i,j;

  HW_NODE * N = Graph[node_id];

  // We'll want to know what starting nucleotide we'd be using
  // (according to the node we arrived from) so we can (a) ensure
  // that we contribute at least one codon and (b) score the
  // mapping
  int start_nucl;
  if (source_id >= 0) {
    i=0;
    while (N->IncomingID[i] != source_id)
      i++;
    start_nucl = N->FirstCodingNucl[i];
  } else {
    start_nucl = 17; // *17*
  }

  // What score do we get from splicing in at this point?
  float tp_score = N->TPSS[start_nucl-1];

  // Our second special catch is if this exon hits the target
  // amino we know ahead of time that this is as far as we can
  // go.
  if (N->end_amino == target_amino) {

    // We'll still need to get the score of the mapping under the
    // given splice position
    int amino_depth = start_amino_depth;
    int map_nucl_len = start_map_nucl_len;    
    float map_score = ExtendTransMap(N,start_nucl,strlen(N->Nucls)-17,Seq,
				     &amino_depth,MapNucls,&map_nucl_len);
    map_score += map_cum_score + tp_score;

    if (map_score > TopScoreThroughNode[node_id]) {
      TopScoreThroughNode[node_id] = map_score;
      TopScoreSourceNode[node_id]  = source_id;
    }
    return map_score;
    
  }

  // Rats!  Why can't life be easy all the time?
  // Oh well, let's get recursive with it...
  int top_target;
  float top_score = MBB_NINF;
  for (i=0; i<N->num_outgoing; i++) {

    // Confirm that this is a reasonable splice pairing (I'm a splice-mmelier)
    int end_nucl = N->LastCodingNucl[i];
    if (end_nucl - start_nucl < 2) continue;

    float fp_score = N->FPSS[end_nucl+1];

    int amino_depth = start_amino_depth;
    int map_nucl_len = start_map_nucl_len;
    float map_score = ExtendTransMap(N,start_nucl,end_nucl,Seq,&amino_depth,MapNucls,
				     &map_nucl_len);
    map_score += map_cum_score + tp_score + fp_score + SPLICE_COST;

    // Wait a minute there, fella -- I ain't seen you before, now, have I?
    if (N->CumScoreMemo[i] == MBB_NINF) {
      
      // RECURSE!
      float recur_score =
	RecursivePathEval(Graph,Seq,N->OutgoingID[i],node_id,TopScoreThroughNode,
			  TopScoreSourceNode,TopScoreTargetNode,target_amino,MapNucls,
			  map_nucl_len,amino_depth,map_score);

      // Phew, that was some hard work! Sure don't want to have to do it again!
      N->CumScoreMemo[i] = recur_score - map_score;
      map_score = recur_score;

    } else {

      // Aha! I thought you looked familiar!
      map_score += N->CumScoreMemo[i];
      
    }

    // Put it all together and what do you get? A score!
    if (map_score > top_score) {
      top_score  = map_score;
      top_target = N->OutgoingID[i];
    }

  }


  // If this is where we're returning, we've spliced!
  if (top_score > TopScoreThroughNode[node_id]) {
    TopScoreThroughNode[node_id] = top_score;
    TopScoreSourceNode[node_id]  = source_id;
    TopScoreTargetNode[node_id]  = top_target;
  }
  return top_score;

}




/*
 *  Function: GetHitExonChrStart
 *
 */
void GetHitExonChrStart
(
 HW_NODE * N,
 int   in_node_id,
 int * chr_start_nucl,
 int * in_edge_id
 ){

  // What's the start nucleotide of the full exon?
  int nucl = N->start_nucl;

  // Find where this edge lives
  int edge_id = 0;
  while (N->IncomingID[edge_id] != in_node_id)
    edge_id++;
  *in_edge_id = edge_id;

  // How far over do we shift to get to the first nucleotide used
  // in the spliced exon?
  int offset = N->FirstCodingNucl[edge_id] - 17; // *17*

  // Because we allow for the possibility of inconsistent graphs, we'll
  // check revcomp for each exon individually
  if (N->start_nucl > N->end_nucl) nucl -= offset; // Reverse strand
  else                             nucl += offset; // Forward strand
  *chr_start_nucl = nucl;
  
}



/*
 *  Function: GetHitExonChrEnd
 *
 */
void GetHitExonChrEnd
(
 HW_NODE * N,
 int   out_node_id,
 int * chr_end_nucl,
 int * out_edge_id
 ){

  // What's the first nucleotide of the full exon?
  int nucl = N->start_nucl;

  int edge_id = 0;
  while (N->OutgoingID[edge_id] != out_node_id)
    edge_id++;
  *out_edge_id = edge_id;

  // How far over do we shift to get to the last nucleotide used
  // in the spliced exon?
  int offset = N->LastCodingNucl[edge_id] - 17; // *17*

  // Because we allow for the possibility of inconsistent graphs, we'll
  // check revcomp for each exon individually
  if (N->start_nucl > N->end_nucl) nucl -= offset; // Reverse strand
  else                             nucl += offset; // Forward strand
  *chr_end_nucl = nucl;
  
}



/*
 *  Function: PrintHitMeta
 *
 */
void PrintHitMeta
(
 int * ExonIDs,
 int   start_amino,
 int   end_amino,
 int   start_nucl,
 int   end_nucl,
 int   num_exons,
 float score,
 HW_OPTS * Opts
 ){

  int i;

  printf("\n  ----- Exons Spliced -----\n\n");

  // We'll ID the hit by the set of exons it used, using biologist-friendly indexing
  printf("  Spliced Exon IDs : %d",ExonIDs[0]+1);
  for (i=1; i<num_exons; i++)
    printf(",%d",ExonIDs[i]+1);
  printf("\n");

  // Give the score of the hit
  printf("  Alignment Score  : %.2f\n",score);

  // Give the range of aminos, as well as the broad span of nucleotides
  // If consistency requirements have been disabled, we'll note that here
  printf("  Amino Acid Range : %d..%d\n",start_amino,end_amino);
  if (Opts->req_consistency)
    printf("  Nucleotide Range : %d..%d\n",start_nucl,end_nucl);
  else
    printf("  Nucleotide Range : Undetermined (inconsistencies permitted)\n");
  printf("\n");
    

}




/*
 *  Function: PrintExonMeta
 *
 *  NOTE that this seems dumb, but in case we want to reformat it'll be nice to
 *  have this output stuck in one place
 *
 */
void PrintExonMeta
(
 int exon_id,
 int start_amino,
 int end_amino,
 int start_nucl,
 int end_nucl
 ){
  printf("  Exon %d: Aminos %d..%d, Nucls %d..%d\n",
	 exon_id+1,start_amino,end_amino,start_nucl,end_nucl);
}





/*
 *  Function: PrintSplicedAlignment
 *
 */
void PrintSplicedAlignment
(
 int  * ExonIDs,
 int  * StartNucls,
 int  * EndNucls,
 int  * StartAminos,
 int  * EndAminos,
 char * TransSeq,
 char * NuclSeq,
 char * RefSeq,
 int    num_exons,
 float  score,
 int    out_str_len,
 HW_OPTS * Opts
 ){

  int i,j,k;
  int line_len = 60;

  // No matter which option we pick, we'll be printing out the hit metadata
  // at the start.
  PrintHitMeta(ExonIDs,StartAminos[0],EndAminos[num_exons-1],StartNucls[0],
	       EndNucls[num_exons-1],num_exons,score,Opts);
  
  // If we're doing output method 1, we'll dump all of the metadata first,
  // then print the alignment as one big thing
  if (Opts->output_format == 1) {
    for (i=0; i<num_exons; i++) {
      PrintExonMeta(ExonIDs[i],StartAminos[i],EndAminos[i],StartNucls[i],EndNucls[i]);
    }
    printf("\n");
  }

  char * Str;
  int exon_num = 0;
  int read_pos = 0;
  int exon_adv; // How many exons are we about to advance?
  while (exon_num < num_exons) {

    // If we're doing output method 2, we'll dump individual metadata lines as
    // we go
    if (Opts->output_format == 2) {
      if (exon_num) printf("\n");
      PrintExonMeta(ExonIDs[exon_num],StartAminos[exon_num],EndAminos[exon_num],
		    StartNucls[exon_num],EndNucls[exon_num]);
      printf("\n");
    }

    Str = TransSeq; 
    i = 0; // Which string are we writing out?
    while (i<3) {

      exon_adv = 0; 
      j = read_pos; // Our writing position      
      printf("  "); // Just some white space to head things off

      while (j<read_pos+line_len) {

	printf("%c",Str[j++]);

	// If we're at the end of our exon, we'll do some careful stuff
	if (Str[j]=='|' || j==out_str_len) {
	  exon_adv++;
	  if (Opts->output_format == 2 || exon_num+exon_adv == num_exons)
	    break;
	}

      }
      
      // Move it forward!
      printf("\n");
      if (i==0) Str = NuclSeq;
      else      Str = RefSeq;
      i++;
      
    }

    printf("\n");

    // Set up for the next set!
    exon_num += exon_adv;
    read_pos  = j;
    exon_adv  = 0;

  }
  
  printf("\n");
    
}




/*
 *  Function: ReportSplicing
 *
 */
void ReportSplicing
(
 HW_NODE ** Graph,
 char * Seq,
 int  * ExonIDs,
 int num_exons,
 int num_aminos,
 float score,
 HW_OPTS * Opts
 ){

  int i,j;

  // We'll load up a bunch of arrays with the data we want to put out,
  // and then we'll do the actual outputery at the end.
  
  // We'll start off by nabbing the chromosome start/end positions, and
  // building up a list of internal Outgoing indices.
  int StartNuclOut[num_exons];
  int EndNuclOut[num_exons];
  int InEdgeIndex[num_exons];
  int OutEdgeIndex[num_exons];

  // We'll treat the first and last a little bit specially, since they don't
  // get spliced.
  HW_NODE * N = Graph[ExonIDs[0]];
  StartNuclOut[0] = N->start_nucl;
  InEdgeIndex[0]  = -1;
  for (i=0; i<num_exons-1; i++) {
    GetHitExonChrEnd(N,ExonIDs[i+1],&EndNuclOut[i],&OutEdgeIndex[i]);
    N = Graph[ExonIDs[i+1]];
    GetHitExonChrStart(N,ExonIDs[i],&StartNuclOut[i+1],&InEdgeIndex[i+1]);
  }
  EndNuclOut[i]   = N->end_nucl;
  OutEdgeIndex[i] = -1;

  // Now we can actually build up the output (and figure out those pesky amino indices)
  int  outseqbuffer = 5*num_aminos; // We could be more precise, but it doesn't matter
  char TransSeqOut[outseqbuffer];
  char NuclSeqOut[outseqbuffer];
  char RefSeqOut[outseqbuffer];
  int  StartAminoOut[num_exons];
  int  EndAminoOut[num_exons];

  // We'll need to be able to translate -- DUH!
  char Codon[3];
  int codon_pos = 0;
  int write_pos = 0;
  int amino_index = Graph[ExonIDs[0]]->start_amino;
  int last_ref_write; // In case we split a codon and need to write to a distant pos.
  for (i=0; i<num_exons; i++) {

    N = Graph[ExonIDs[i]];

    // NOTE that we want to friendly to biologists, so we'll index starting at '1'
    // Moreover, we'll keep this as '+1' so we can do a straight amino_count as the
    // EndAminoOut at the end of the loop.
    StartAminoOut[i] = amino_index + 1;

    int start_index;
    if (i) start_index = N->FirstCodingNucl[InEdgeIndex[i]];
    else   start_index = 17; // *17*

    int end_index;
    if (i<num_exons-1) end_index = N->LastCodingNucl[OutEdgeIndex[i]];
    else               end_index = N->nucl_str_len - 18; // *17*
    
    // Grab the leading dinucls (if applicable)
    // This is also a great opportunity to pull in the score of this edge!
    if (i) {

      start_index -= 2;

      TransSeqOut[write_pos] = ' ';
      NuclSeqOut[write_pos]  = N->Nucls[start_index++];
      RefSeqOut[write_pos]   = ' ';
      if (NuclSeqOut[write_pos] < 96) NuclSeqOut[write_pos] += 32;
      write_pos++;

      TransSeqOut[write_pos] = ' ';
      NuclSeqOut[write_pos]  = N->Nucls[start_index++];
      RefSeqOut[write_pos]   = ' ';
      if (NuclSeqOut[write_pos] < 96) NuclSeqOut[write_pos] += 32;
      write_pos++;

    }

    // WRITE OUT THIS EXON!!!
    while (start_index <= end_index) {

      NuclSeqOut[write_pos] = N->Nucls[start_index];
      if (NuclSeqOut[write_pos] > 96) NuclSeqOut[write_pos] -= 32;
      Codon[codon_pos] = NuclSeqOut[write_pos];

      // What we do next depends on the fill level of our codon -- if we've
      // filled one up, then we translate; if we're at the middle nucl we record
      // the ref; if we're starting one off, we take a rest
      if (codon_pos == 2) {
	TransSeqOut[last_ref_write] = MBB_TranslateCodon(Codon);
	TransSeqOut[write_pos] = ' ';
	RefSeqOut[write_pos] = ' ';
	codon_pos = 0;
      } else if (codon_pos == 1) {
	TransSeqOut[write_pos] = ' ';
	RefSeqOut[write_pos] = Seq[amino_index++];
	last_ref_write = write_pos;
	codon_pos = 2;
      } else {
	TransSeqOut[write_pos] = ' ';
	RefSeqOut[write_pos] = ' ';
	codon_pos = 1;
      }

      start_index++;
      write_pos++;
      
    }

    // Toss in the tailing dinucls (if applicable)
    if (i<num_exons-1) {
      
      TransSeqOut[write_pos] = ' ';
      NuclSeqOut[write_pos]  = N->Nucls[start_index++];
      RefSeqOut[write_pos]   = ' ';
      if (NuclSeqOut[write_pos] < 96) NuclSeqOut[write_pos] += 32;
      write_pos++;

      TransSeqOut[write_pos] = ' ';
      NuclSeqOut[write_pos]  = N->Nucls[start_index++];
      RefSeqOut[write_pos]   = ' ';
      if (NuclSeqOut[write_pos] < 96) NuclSeqOut[write_pos] += 32;
      write_pos++;

      // Put in an indicator that we're splicing!
      TransSeqOut[write_pos] = '|';
      NuclSeqOut[write_pos]  = '|';
      RefSeqOut[write_pos]   = '|';
      write_pos++;
      
    }

    EndAminoOut[i] = amino_index;
    
  }

  // TIME TO ACTUALLY PRINT THIS STUFF OUT!
  //
  // Just to have a quick guide to the variables we're interested in
  //
  //
  // Coordinates:
  //    ExonIDs
  //    StartNuclOut
  //    EndNuclOut
  //    StartAminoOut
  //    EndAminoOut
  //
  // Sequences:
  //    TransSeqOut
  //    NuclSeqOut
  //    RefSeqOut
  //
  PrintSplicedAlignment(ExonIDs,StartNuclOut,EndNuclOut,StartAminoOut,EndAminoOut,
			TransSeqOut,NuclSeqOut,RefSeqOut,num_exons,score,write_pos,
			Opts);
  
}





/*
 *  Function: ReportMaximalPaths
 *
 */
void ReportMaximalPaths
(
 HW_NODE ** Graph,
 int num_exons,
 char * Seq,
 HW_OPTS * Opts
 ){

  int i,j;
  
  // We'll begin by doing a big ol' scan through the graph
  // to figure out how many connected components there are
  // and which nodes belong to which components
  //
  // Note that '0' is 'unassigned'
  //
  int ComMembership[num_exons];
  for (i=0; i<num_exons; i++)
    ComMembership[i] = 0;
  int num_components = 0;

  // We'll do our traversal with a stack, which we'll initialize...
  // NOW!
  HW_NODE * NodeStack[num_exons];
  HW_NODE * NextNode;
  int nodestacksize = 0;

  // Simple DFS-y approach
  i=0;
  while (i<num_exons) {

    // Ooooh, boy, here comes another component!
    if (!ComMembership[i]) {
      num_components++;
      ComMembership[i] = num_components;
      NodeStack[0]  = Graph[i];
      nodestacksize = 1;
    }

    while (nodestacksize) {

      // Pop the next node off the stack
      nodestacksize--;
      NextNode = NodeStack[nodestacksize];

      // See if any of this node's incoming friends are (as-of-yet)
      // unaccounted for
      for (j=0; j<NextNode->num_incoming; j++) {
	
	// Looks like somebody's joining the team!
	if (!ComMembership[NextNode->IncomingID[j]]) {
	  ComMembership[NextNode->IncomingID[j]] = num_components;
	  NodeStack[nodestacksize++] = NextNode->Incoming[j];
	}
	
      }

      // See if any of this node's outgoing pals are (at-this-time)
      // unassigned
      for (j=0; j<NextNode->num_outgoing; j++) {

	// You're hired!
	if (!ComMembership[NextNode->OutgoingID[j]]) {
	  ComMembership[NextNode->OutgoingID[j]] = num_components;
	  NodeStack[nodestacksize++] = NextNode->Outgoing[j];
	}
	
      }

    }

    i++;

  }

  // Well, that was surprisingly easy -- we can now figure out the
  // amino acid ranges covered by each connected component.
  //
  // NOTE: I'm fairly certain it's impossible for there to be a
  //       pathological graph where the min_amino and max_amino for
  //       a connected component can't be reached through forward
  //       movement alone. FAIRLY certain.
  //
  // +1 because '0' is unassigned, and easy indexing is better than
  // ultra-nitpicker-level memory efficiency.
  //
  int ComMinAmino[num_components+1];
  int ComMaxAmino[num_components+1];
  for (i=1; i<=num_components; i++) {
    ComMinAmino[i] = 1000000;
    ComMaxAmino[i] = -1;
  }

  for (i=0; i<num_exons; i++) {
    j = ComMembership[i];
    ComMinAmino[j] = MBB_MinInt(ComMinAmino[j],Graph[i]->start_amino);
    ComMaxAmino[j] = MBB_MaxInt(ComMaxAmino[j],Graph[i]->end_amino);
  }

  // So now we're going to do some nasty business.
  //
  // Basically, the plan is to find the highest-scoring path of 
  // longest coverage through each component, and then report it.
  //
  // I'm just now realizing that the work we do to sort scores in
  // 'ConnectGraph' might be unnecessary... I'm going to leave that
  // in for now, but it may not survive through to release.
  //
  // We'll take a DFS-based approach to determining this path, but
  // this time using recursive function calls, for easy retrieval
  // of scores.

  float TopScoreThroughNode[num_exons];
  int   TopScoreSourceNode[num_exons]; // the prev node on the highest-scoring path
  int   TopScoreTargetNode[num_exons]; // the next node on the highest-scoring path
  float ToppestScores[num_components+1];
  int   ToppestSources[num_components+1]; // the first node on the highest-scoring path

  // We'll need to build up the nucleotide sequence as we examine splicings
  // to determine the actual scores.
  char MapNucls[strlen(Seq)*3];

  // Initialize!
  for (i=1; i<=num_components; i++)
    ToppestScores[i] = MBB_NINF;

  // Initialize more!
  for (i=0; i<num_exons; i++) {
    TopScoreSourceNode[i] = -1;
    TopScoreTargetNode[i] = -1;
    TopScoreThroughNode[i] = MBB_NINF;
  }

  // Now run through our exons looking for possible starts, those
  // being exons whose start aminos equal ComMinAmino for their
  // component
  for (i=0; i<num_exons; i++) {

    // Figure out which component we're in, and see if we're a valid
    // starting position.
    j = ComMembership[i];
    if (Graph[i]->start_amino == ComMinAmino[j]) {

      // Are you my dad?
      float pathscore
	= RecursivePathEval(Graph,Seq,i,-1,TopScoreThroughNode,TopScoreSourceNode,
			    TopScoreTargetNode,ComMaxAmino[j],MapNucls,0,
			    Graph[i]->start_amino,0.0);

      if (pathscore > ToppestScores[j]) {
	ToppestScores[j]  = pathscore;
	ToppestSources[j] = i;
      }
      
    }
  }

  // We can now report the highest-scoring path through each connected
  // component among the set of paths that go from the min amino for the
  // component to the max amino!
  for (i=1; i<=num_components; i++) {

    // Like I've mentioned before, I'm fairly certain we're guaranteed at least
    // one path per component, but I'll add a check for positive scores just to
    // be certain.
    if (ToppestScores[i] <= 0.0)
      continue;

    // How many exons are in this path?
    j = TopScoreTargetNode[ToppestSources[i]];
    int num_com_exons = 1;
    while (j != -1) {
      num_com_exons++;
      j = TopScoreTargetNode[j];
    }

    // While we may not be especially interested in single-exon outputs, we are
    // provided the additional context that they're part of their own (unconnected)
    // components, which may be informative for figuring out how to string together
    // partial spliced mappings (esp. when building a bestiary of weirdness)
    if (num_com_exons == 1 && Opts->report_singles == 0)
      continue;
    
    // UGH, just list them for me
    int ComExons[num_com_exons];
    ComExons[0] = ToppestSources[i];
    for (j=1; j<num_com_exons; j++)
      ComExons[j] = TopScoreTargetNode[ComExons[j-1]];

    // Also, how many aminos (not to pry)?
    int num_com_aminos = ComMaxAmino[i] - ComMinAmino[i] + 1;

    // To keep things a little bit cleaner than they'd be otherwise, let's offload
    // the next bit of work to a dedicated output prep. function.
    ReportSplicing(Graph,Seq,ComExons,num_com_exons,num_com_aminos,ToppestScores[i],Opts);

  }

  // DONE!
  
}




/*
 *  Function: ParseCommandArgs
 *
 */
void ParseCommandArgs (HW_OPTS * Opts, int argc, char ** argv) {

  // Initialize the option variables

  // Do we require that hits follow one another on the chromosome?
  // The reason we'd be interested in turning this off is if we're investigating
  // strange splicing phenomena.
  Opts->req_consistency = 1;

  // Output options
  //  == 1 : We'll print out our metadata first, followed by the complete alignment
  //  == 2 : We'll print out each exon individually, as paired metadata / data
  Opts->output_format = 1;

  // Only print spliced alignments with multiple exons (turning on allows checking
  // for single-node unconnected components of graph)
  Opts->report_singles = 0;

  int i=1;
  while (i<argc-1) {
    if (!strcmp(argv[i],"--allow-inconsistency")) {
      Opts->req_consistency = 0;
    } else if (!strcmp(argv[i],"-outfmt")) {
      Opts->output_format = atoi(argv[++i]);
    } else if (!strcmp(argv[i],"--report-singles")) {
      Opts->report_singles = 1;
    } else {
      fprintf(stderr,"  ExonWeaver Warning: Unrecognized option '%s' ignored\n",argv[i]);
    }
    i++;
  }
  
}





/*
 *  Function: PrintUsage
 *
 */
int PrintUsage () {
  printf("\n  USAGE: ./ExonWeaver {OPTS} [weaver-input-file]\n\n");
  return 1;
}





/*
 *         M  A  I  N
 *
 */
int main (int argc, char ** argv) {

  if (argc < 2) return PrintUsage();

  // Parse commandline arguments, stored to an options struct
  HW_OPTS * Opts = malloc(sizeof(HW_OPTS));
  ParseCommandArgs(Opts,argc,argv);
  
  // We'll kick things off by reading that big stinky file
  FILE * inf;
  inf = fopen(argv[argc-1],"r");

  // First line tells us the number of exons that we'll be considering
  int num_exons;
  fscanf(inf,"Num Hits : %d\n",&num_exons);

  // Second line tells us the length of the protein sequence
  int seq_len;
  fscanf(inf,"Seq Len  : %d\n",&seq_len);

  // Third line is the protein sequence itself (all uppercase)
  char * Seq = malloc((seq_len+1)*sizeof(char));
  fscanf(inf,"Sequence : %s\n",Seq);
  
  // Rad! Now that we've covered the metadata, let's get down to the data!
  HW_NODE ** Graph = malloc(num_exons*sizeof(HW_NODE *));
  ReadNodesFromFile(Graph,num_exons,inf);

  // Reading time is over -- now it's time to get graphical!
  fclose(inf);

  // Draw up connections between nodes wherever possible.
  // These edges end up being sorted by score (highest first)
  ConnectGraph(Graph,num_exons,Seq,Opts);

  // Identify the connected components of the graph, and reduce
  // down to non-overlapping (in protein seq.) connected components
  // with maximum coverage.
  ReportMaximalPaths(Graph,num_exons,Seq,Opts);

  // Once we're all through, we are become Free, Destroyer of Arrays
  DestroyGraph(Graph,num_exons);
  free(Opts);
  free(Seq);

  // And that's called ExonWeaver!
  return 0;
  
}


