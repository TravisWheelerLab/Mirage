final_build: build_FindDiagonals build_MultiSeqNW build_TransSW clean

build_FindDiagonals: FindDiagonals.c Diagonals.c
	gcc -O2 -o FindDiagonals FindDiagonals.c Diagonals.c -lm

build_MultiSeqNW: MultiSeqNW.c
	gcc -O2 -o MultiSeqNW MultiSeqNW.c -lm

build_TransSW: TransSW.c Diagonals.c
	gcc -O2 -o TransSW TransSW.c Diagonals.c -lm

clean:
	rm -f *.o *~
