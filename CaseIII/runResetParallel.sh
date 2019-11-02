if [ -f "log" ];
then
rm log
fi
if [ -f "dump" ];
then
rm dump
fi
if [ -d "intermediate" ]
then
rm -r intermediate
fi
mkdir intermediate

qcc -fopenmp -Wall -O2 dropOnDropImpact.c -o dropOnDropImpact -lm

export OMP_NUM_THREADS=8
./dropOnDropImpact
