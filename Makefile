CC = g++
FLAGS = -pedantic -Wextra -Wall -Weffc++

all: compr uncompr

compr : compr.cc bin globaldefs.h stream.o
	$(CC) $(FLAGS) -lboost_program_options-mt compr.cc -o bin/compr

uncompr : uncompr.cc bin globaldefs.h
	$(CC) $(FLAGS) -lboost_program_options-mt uncompr.cc -o bin/uncompr

stream.o : stream.h stream.cc
	$(CC) $(FLAGS) stream.cc -c -o bin/stream.o

bin :
	mkdir bin

clean :
	rm bin/*
	rmdir bin
