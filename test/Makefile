CC = g++
FLAGS = -pedantic -Wextra -Wall -Weffc++ -g

tests : streamtest preproctest coderstest
	./streamtest
	./preproctest
	./coderstest

streamtest : stream_test.cc testdefs.h
	$(CC) $(FLAGS) ../bin/stream.o stream_test.cc -lboost_filesystem-mt \
	-o streamtest

preproctest : preproc_test.cc testdefs.h
	$(CC) $(FLAGS) ../bin/block.o ../bin/preprocessor.o ../bin/stream.o \
	../bin/block_manager.o preproc_test.cc -o preproctest

coderstest : coders_test.cc testdefs.h
	$(CC) $(FLAGS) ../bin/coders.o ../bin/rl_compress.o ../bin/stream.o \
	coders_test.cc -o coderstest

clean :
	rm coderstest
	rm streamtest
	rm preproctest
	rm testfile.txt

