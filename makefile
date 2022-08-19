CC = g++

CFLAGS= -O3 -march=native -std=c++17 -Wall
LDFLAGS = -LSIMDCompressionAndIntersection -l:libSIMDCompressionAndIntersection.a -LTurboPFor-Integer-Compression -l:libic.a -pthread -lpthread

EXEC = index
OBJ = index.o main.o utils.o fastDelta.o bloomFilter.o

all : $(EXEC) 

index : $(OBJ)
	$(CC) -o index $^ $(LDFLAGS)

index.o: sources/index.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: sources/main.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: sources/utils.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

fastDelta.o: include/fastDelta.c
	$(CC) -o $@ -c $< $(CFLAGS)

bloomFilter.o: sources/bloomfilter.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(EXEC) $(OBJ) *.dat

rebuild: clean $(EXEC)