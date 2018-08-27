final_build: build_FastMap build_MultiSeqNW build_TransSW clean

build_FastMap: FastMap.c Diagonals.c
	gcc -O2 -o FastMap FastMap.c Diagonals.c -lm

build_MultiSeqNW: MultiSeqNW.c
	gcc -O2 -o MultiSeqNW MultiSeqNW.c -lm

build_TransSW: TransSW.c Diagonals.c
	gcc -O2 -o TransSW TransSW.c Diagonals.c -lm

clean:
	rm -f *.o *~
