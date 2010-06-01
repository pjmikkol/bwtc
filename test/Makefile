CC = g++
FLAGS = -pedantic -Wextra -Wall -Weffc++ -g

#tests : streamtest
#	./streamtest

streamtest : stream_test.cc
	$(CC) $(FLAGS) ../bin/stream.o stream_test.cc -lboost_filesystem-mt -o streamtest

clean :
	rm streamtest
	rm temptest.txt

