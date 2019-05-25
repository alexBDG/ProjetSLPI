FFLAGS = -g -O0 -w -std=c++11 -I Eigen/

main.exe: methodes.h methodes.o main.o
	g++ $(FFLAGS) methodes.o main.o -o main.exe

main.o: main.cpp methodes.h methodes.cpp
	g++ $(FFLAGS) -c methodes.cpp main.cpp

methodes.o: methodes.cpp methodes.h
	g++ $(FFLAGS) -c methodes.cpp

clean:
	rm -f *~ *.exe *.out
