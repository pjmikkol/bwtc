CC = g++
FLAGS = -pedantic -Wextra -Wall -Weffc++ -g

tests : streamtest preproctest
	./streamtest
	./preproctest

streamtest : stream_test.cc testdefs.h
	$(CC) $(FLAGS) ../bin/stream.o stream_test.cc -lboost_filesystem-mt -o streamtest

preproctest : preproc_test.cc testdefs.h
	$(CC) $(FLAGS) ../bin/block.o ../bin/preprocessor.o ../bin/stream.o preproc_test.cc -o preproctest

clean :
	rm streamtest
	rm preproctest
	rm testfile.txt

