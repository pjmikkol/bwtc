CC = g++
FLAGS = -pedantic -Wextra -Wall -Weffc++ -g

all: bin/compr bin/uncompr

bin/compr : compr.cc globaldefs.h bin/stream.o bin/preprocessor.o bin/block.o \
	bin/block_manager.o bin/coders.o 
	$(CC) $(FLAGS) -lboost_program_options-mt compr.cc bin/stream.o \
	bin/preprocessor.o bin/block.o bin/block_manager.o bin/coders.o \
	bin/rl_compress.o -o bin/compr

bin/uncompr : uncompr.cc globaldefs.h bin/coders.o bin/rl_compress.o \
	bin/stream.o
	$(CC) $(FLAGS) -lboost_program_options-mt bin/coders.o uncompr.cc \
	bin/rl_compress.o bin/stream.o -o bin/uncompr

bin/stream.o : stream.h stream.cc 
	$(CC) $(FLAGS) stream.cc -c -o bin/stream.o

bin/preprocessor.o : preprocessor.h preprocessor.cc stream.h bin/block.o 
	$(CC) $(FLAGS) preprocessor.cc -c -o bin/preprocessor.o

bin/block.o : block.h block.cc 
	$(CC) $(FLAGS) block.cc -c -o bin/block.o

bin/coders.o : coders.cc coders.h probmodels/base_prob_model.h rl_compress.h \
	bin/rl_compress.o 
	$(CC) $(FLAGS) coders.cc -c -o bin/coders.o

bin/block_manager.o : block_manager.cc block_manager.h block.h 
	$(CC) $(FLAGS) block_manager.cc -c -o bin/block_manager.o

bin/rl_compress.o : rl_compress.cc rl_compress.h globaldefs.h 
	$(CC) $(FLAGS) rl_compress.cc -c -o bin/rl_compress.o 

clean :
	rm -f bin/*
	rm -f test/coderstest
	rm -f test/streamtest
	rm -f test/preproctest
	rm -f test/testfile.txt

# Rest of the file is for tests:
tests : test/streamtest test/preproctest test/coderstest
	./test/streamtest
	./test/preproctest
	./test/coderstest

test/streamtest : test/stream_test.cc test/testdefs.h bin/stream.o
	$(CC) $(FLAGS) bin/stream.o test/stream_test.cc -lboost_filesystem-mt \
	-o test/streamtest

test/preproctest : test/preproc_test.cc test/testdefs.h bin/block.o \
	bin/preprocessor.o bin/stream.o bin/block_manager.o
	$(CC) $(FLAGS) bin/block.o bin/preprocessor.o bin/stream.o \
	bin/block_manager.o test/preproc_test.cc -o test/preproctest

test/coderstest : test/coders_test.cc test/testdefs.h bin/coders.o \
	bin/stream.o bin/rl_compress.o
	$(CC) $(FLAGS) bin/coders.o bin/rl_compress.o bin/stream.o \
	test/coders_test.cc -o test/coderstest
