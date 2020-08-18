#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _BasicBio_H_
#define _BasicBio_H_

static float MBB_INF  =  1.0/0.0;
static float MBB_NINF = -1.0/0.0;

int MBB_MinInt(int a, int b);
int MBB_MaxInt(int a, int b);
float MBB_MinFloat(float a, float b);
float MBB_MaxFloat(float a, float b);

int MBB_DNAtoNum(char residue);
char MBB_TranslateCodon(char * Codon);
float MBB_AminoAliScore(char amino1, char amino2);

#endif
