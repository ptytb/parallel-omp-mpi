
heatcond: heatcond.o
	mpic++ heatcond.o -o heatcond

.cpp.o:
	mpic++ -g -O0 -c -o $@ $*.cpp

heatcond.o: heatcond.cpp

clean:
	rm -f heatcond.o heatcond

.PHONY: clean

