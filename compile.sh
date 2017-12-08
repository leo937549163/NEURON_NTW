#!/bin/sh

#just to clean all old object files !
rm -f *.o *.so

#In case of without MPI
gcc -Wall -fPIC -c connectivity.c cppInterface.cpp species.cpp reaction.cpp SV.cpp Main.cpp -pg
g++ -shared -o connectivity.so connectivity.o cppInterface.o species.o reaction.o SV.o Main.o

gcc -Wall -fPIC -c reaction_function.c -pg
g++ -shared -o reaction_function.so reaction_function.o species.o reaction.o

#In case of with MPI
# mpicc.openmpi -Wall -fPIC -c connectivity.c cppInterface.cpp species.cpp reaction.cpp SV.cpp Main.cpp
# mpic++.openmpi -shared -o connectivity.so connectivity.o cppInterface.o species.o reaction.o SV.o Main.o

#to run CaBuf model
# python CaBuf_stochastic_low.py
# python CaBuf_stochastic_high.py
# python CaBuf_deterministic_low.py
# python CaBuf_deterministic_high.py

# Profiler
python setup.py build_ext --inplace -lprofiler
# pprof --text connectivity.so output.prof