CC = g++
LIBS = -lblas -llapacke -lstdc++fs
FLAGS = -fopenmp
FEHD: driver_new.cpp FEHD.cpp GCopenmp.cpp kernelsCPU.cpp mkAR.cpp timeSeriesOPs.cpp utility.cpp
	$(CC) $^ -o $@ $(LIBS) $(FLAGS)
