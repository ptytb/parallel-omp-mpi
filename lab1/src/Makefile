
integ: Taylor.o fobj.o  main.o
	g++ -fopenmp -lgomp Taylor.o fobj.o main.o -o integ

.cpp.o:
	g++ -fopenmp -g -O0 -c -o $@ $*.cpp

Taylor.o: Taylor.cpp
fobj.o: fobj.cpp
main.o: main.cpp

clean:
	rm -f Taylor.o fobj.o main.o integ

.PHONY: clean

