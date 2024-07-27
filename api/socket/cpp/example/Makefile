#CXX=g++
CXX=g++ -std=c++11
#CXX=ccache g++
all: example 

example: CppAPIExample.o
	$(CXX) -o example CppAPIExample.o -L../lib -lgstoreconnector
	#$(CXX) -o example CppAPIExample.o -L../lib -lgstoreconnector -L/usr/local/lib -lboost_filesystem

CppAPIExample.o: CppAPIExample.cpp
	$(CXX) -c -I../src/ CppAPIExample.cpp -o CppAPIExample.o

#QueryGenerator: QueryGenerator.cpp
#	$(CXX) -o QueryGenerator QueryGenerator.cpp

test: test.o
	$(CXX) -o test test.o -L../lib -lgstoreconnector -pthread

test.o: test.cpp
	$(CXX) -c -I../src/ test.cpp -o test.o -pthread

clean:
	rm -rf *.o example QueryGenerator test
