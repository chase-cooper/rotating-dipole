all : main.o
	g++ -o main.o main.cpp

clean :
	rm -r *.o