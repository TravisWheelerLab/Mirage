final_build: build_FastDiagonals build_MultiSeqNW build_TransSW clean

build_FastDiagonals: FastDiagonals.c Diagonals.c
	gcc -O2 -o FastDiagonals FastDiagonals.c Diagonals.c -lm

build_MultiSeqNW: MultiSeqNW.c
	gcc -O2 -o MultiSeqNW MultiSeqNW.c -lm

build_TransSW: TransSW.c Diagonals.c
	gcc -O2 -o TransSW TransSW.c Diagonals.c -lm

clean:
	rm -f *.o *~
