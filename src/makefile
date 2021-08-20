final_build: build_FastMap2 build_ExonWeaver build_MultiSeqNW clean

build_FastMap2: FastMap2.c BasicBio.c
	gcc -O2 -o FastMap2 FastMap2.c BasicBio.c -lm

build_ExonWeaver: ExonWeaver.c BasicBio.c
	gcc -O2 -o ExonWeaver ExonWeaver.c BasicBio.c -lm

build_MultiSeqNW: MultiSeqNW.c
	gcc -O2 -o MultiSeqNW MultiSeqNW.c -lm

clean:
	rm -f *.o *~
