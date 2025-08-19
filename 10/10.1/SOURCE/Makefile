# Usa il compilatore MPI wrapper invece di g++ direttamente
MPICC = mpic++
CC = g++
# Rimuovi i sanitizer problematici e usa flag compatibili
CFLAGS = -Wall -O2 --std=c++14
# Se vuoi il debug, usa questi flag invece:
# CFLAGS = -Wall -g -O0 --std=c++14

# Per linking con MPI, usa mpic++ che include automaticamente le librerie necessarie
main.exe : main.o ../../../Library/PRNG/random.o ../../../Library/blockAverage/blockAverage.hpp ../../../Library/library.hpp geneticLibrary.hpp cities.hpp
	$(MPICC) ../../../Library/PRNG/random.o main.o -o main.exe $(CFLAGS)

main.o : main.cpp ../../../Library/PRNG/random.hpp ../../../Library/blockAverage/blockAverage.hpp ../../../Library/library.hpp geneticLibrary.hpp cities.hpp
	$(MPICC) -c main.cpp -o main.o $(CFLAGS)

../../../Library/PRNG/random.o : ../../../Library/PRNG/random.cpp
	$(CC) -c ../../../Library/PRNG/random.cpp -o ../../../Library/PRNG/random.o $(CFLAGS)

run : main.exe
	mpirun -np 5 ./main.exe

clean :
	rm -f *.o main.exe