CC = g++
CFLAGS = -Wall -O3 --std=c++11
main.exe : main.o ../../Library/PRNG/random.o
	$(CC) ../../Library/PRNG/random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
run : main.exe
	./main.exe
clean :
	rm *.o main.exe seed.out