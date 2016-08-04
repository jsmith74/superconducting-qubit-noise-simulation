all: RTCGen

EXTRALIBS = -lm -lpthread -lrt -lstdc++ -ldl -luuid
WSTPLIB = WSTP64i4

RTCGen: main.o ftm.o finiteDifferenceFuncs.o RTCGen.o
	g++ -O3  main.o ftm.o finiteDifferenceFuncs.o RTCGen.o -I WSTP -L WSTP -l ${WSTPLIB} ${EXTRALIBS} -o RTCGen
	
main.o: main.cpp
	g++ -O3 -I WSTP -I EIGEN main.cpp -c -march=native

RTCGen.o: reducedTaylorCoefficientGen.cpp
	g++ -O3 -I WSTP -I EIGEN reducedTaylorCoefficientGen.cpp -c -march=native -o RTCGen.o


finiteDifferenceFuncs.o: finiteDifferenceFuncs.cpp
	g++ -O3 -I WSTP -I EIGEN finiteDifferenceFuncs.cpp -c

.c.o :
	g++ -O3 -c -I WSTP $<

ftm.c: f.tm
	./WSTP/wsprep f.tm -o ftm.c

clean:
	rm *o RTCGen

