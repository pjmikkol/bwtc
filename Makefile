CC = g++
FLAGS = -pedantic -Wextra -Wall -Weffc++

all: compr uncompr

compr : compr.cc bin globaldefs.h stream.o preprocessor.o
	$(CC) $(FLAGS) -lboost_program_options-mt compr.cc bin/stream.o \
	bin/preprocessor.o -o bin/compr

uncompr : uncompr.cc bin globaldefs.h
	$(CC) $(FLAGS) -lboost_program_options-mt uncompr.cc -o bin/uncompr

stream.o : stream.h stream.cc
	$(CC) $(FLAGS) stream.cc -c -o bin/stream.o

preprocessor.o : preprocessor.h preprocessor.cc stream.h
	$(CC) $(FLAGS) preprocessor.cc -c -o bin/preprocessor.o

bin :
	mkdir bin

clean :
	rm bin/*
	rmdir bin
