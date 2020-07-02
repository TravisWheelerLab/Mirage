final_build: build_FastMap2 build_HitWeaver build_MultiSeqNW clean

build_FastMap2: FastMap2.c BasicBio.c
	gcc -O2 -o FastMap2 FastMap2.c BasicBio.c -lm

build_HitWeaver: HitWeaver.c BasicBio.c
	gcc -O2 -o HitWeaver HitWeaver.c BasicBio.c -lm

build_MultiSeqNW: MultiSeqNW.c
	gcc -O2 -o MultiSeqNW MultiSeqNW.c -lm

clean:
	rm -f *.o *~
