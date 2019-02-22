all: km_pthreads km_openmp

km_pthreads: km_pthreads.o utils.o 
	g++ -o km_pthreads -pthread km_pthreads.o utils.o

km_openmp: km_openmp.o utils.o
	g++ -o km_openmp -fopenmp km_openmp.o utils.o

utils.o: utils.cpp utils.hh
	g++ -c utils.cpp

km_pthreads.o: km_pthreads.cpp pkmeans.hh
	g++ -c km_pthreads.cpp

km_openmp.o: km_openmp.cpp
	g++ -c km_openmp.cpp

clean:
	rm -f km_pthreads.o km_openmp.o utils.o
