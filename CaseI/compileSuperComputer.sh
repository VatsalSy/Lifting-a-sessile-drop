qcc -source -D_MPI=1 $file.c
mpicc -Wall -std=c99 -O2 -D_MPI=1 _$file.c -o $file -lm
