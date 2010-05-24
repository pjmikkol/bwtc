CC = g++
FLAGS = -pedantic -Wextra -Wall

main : main.cc bin
	$(CC) $(FLAGS) -lboost_program_options-mt main.cc -o bin/main

bin :
	mkdir bin

clean :
	rm bin/*
	rmdir bin