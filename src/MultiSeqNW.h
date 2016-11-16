#ifndef __MULTI_SEQ_NW_H__
#define __MULTI_SEQ_NW_H__

typedef struct _tuple_set_ {
  int     MarksIntron;
  int     NumTuples;
  int     NearestIntron;
  char  * TupleChars;
  float * TupleRatios;
} TUPLE_SET;

TUPLE_SET * InitTupleSet(int num_letters, int is_non_letter);
TUPLE_SET * ConvertToTupleSet(char * letters, int num_letters);
void DestroyTupleSet(TUPLE_SET * tuple_set);

#endif
