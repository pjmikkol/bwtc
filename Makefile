CC = g++
FLAGS = -pedantic -Wextra -Wall

main : main.o
	$(CC) $(FLAGS) -lboost_program_options-mt main.cc -o main