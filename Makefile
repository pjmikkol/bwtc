CC = g++
FLAGS = -pedantic -Wextra -Wall -Weffc++ -g

all: compr uncompr

compr : compr.cc bin globaldefs.h bin/stream.o bin/preprocessor.o bin/block.o \
	bin/block_manager.o bin/coders.o 
	$(CC) $(FLAGS) -lboost_program_options-mt compr.cc bin/stream.o \
	bin/preprocessor.o bin/block.o bin/block_manager.o bin/coders.o \
	bin/rl_compress.o -o bin/compr

uncompr : uncompr.cc bin globaldefs.h
	$(CC) $(FLAGS) -lboost_program_options-mt uncompr.cc -o bin/uncompr

bin/stream.o : stream.h stream.cc
	$(CC) $(FLAGS) stream.cc -c -o bin/stream.o

bin/preprocessor.o : preprocessor.h preprocessor.cc stream.h bin/block.o
	$(CC) $(FLAGS) preprocessor.cc -c -o bin/preprocessor.o

bin/block.o : block.h block.cc
	$(CC) $(FLAGS) block.cc -c -o bin/block.o

bin/coders.o : coders.cc coders.h rl_compress.h bin/rl_compress.o
	$(CC) $(FLAGS) coders.cc -c -o bin/coders.o

bin/block_manager.o : block_manager.cc block_manager.h block.h
	$(CC) $(FLAGS) block_manager.cc -c -o bin/block_manager.o

bin/rl_compress.o : rl_compress.cc rl_compress.h globaldefs.h
	$(CC) $(FLAGS) rl_compress.cc -c -o bin/rl_compress.o 
bin :
	mkdir bin

clean :
	rm bin/*
	rmdir bin
