CC = g++
CFLAGS = -Wall -O3 --std=c++11
main.exe : main.o ../../Library/PRNG/random.o ../../Library/blockAverage/blockAverage.hpp ../../Library/library.hpp 
	$(CC) ../../Library/PRNG/random.o main.o -o main.exe $(CFLAGS)
main.o : main.cpp ../../Library/PRNG/random.hpp ../../Library/blockAverage/blockAverage.hpp ../../Library/library.hpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
../../Library/PRNG/random.o : ../../Library/PRNG/random.cpp 
	$(CC) -c ../../Library/PRNG/random.cpp -o ../../Library/PRNG/random.o $(CFLAGS)
run : main.exe
	./main.exe
clean :
	rm *.o main.exe