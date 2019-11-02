#Compiles Basilisk scripts without bview stuff.
qcc -O2 -Wall $file.c -o $file -lm
