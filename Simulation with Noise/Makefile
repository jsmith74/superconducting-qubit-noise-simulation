all: TEvolve

TEvolve: main.o importHFuncs.o timeEvolutionGenFuncs.o
	mpic++ -O3 main.o importHFuncs.o timeEvolutionGenFuncs.o -o TEvolve

main.o: main.cpp
	mpic++ -O3 -I EIGEN main.cpp -c -march=native

timeEvolutionGenFuncs.o: timeEvolutionGenFuncs.cpp
	mpic++ -O3 -I EIGEN timeEvolutionGenFuncs.cpp -c -march=native

importHFuncs.o: importHFuncs.cpp
	mpic++ -O3 -I EIGEN importHFuncs.cpp -c -march=native

clean:
	rm *o TEvolve
