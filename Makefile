CC = g++
FLAGS = -pedantic -Wextra -Wall

all: compr uncompr

compr : compr.cc bin globaldefs.h
	$(CC) $(FLAGS) -lboost_program_options-mt compr.cc -o bin/compr

uncompr : uncompr.cc bin globaldefs.h
	$(CC) $(FLAGS) -lboost_program_options-mt uncompr.cc -o bin/uncompr

bin :
	mkdir bin

clean :
	rm bin/*
	rmdir bin