FFLAGS = -O0 -std=c++11 -I Eigen/Eigen -w

TableRase: run
	@ rm *.o
	clear

run: Methodes.o main.o Lecture.o
	@ g++ $(FFLAGS) Methodes.o main.o Lecture.o -o run

main.o: Methodes.cpp main.cc
	@ g++ $(FFLAGS) -c Methodes.cpp Lecture.cpp main.cc

Methodes.o: Methodes.cpp
	@ g++ $(FFLAGS) -c Methodes.cpp

Lecture.o: Lecture.cpp
	@ g++ $(FFLAGS) -c Lecture.cpp

clean:
	@ rm -f *.o run *~ *.txt
