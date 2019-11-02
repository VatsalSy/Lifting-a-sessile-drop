#Compiles Basilisk scripts with bview stuff.
qcc -O2 -Wall $file.c -o $file -lm \
   -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm
