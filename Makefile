CC = g++
DFLAGS = -g
FLAGS = -pedantic -Wextra -Wall $(DFLAGS)
TFLAGS = -Wall $(DFLAGS) # because of template assertions

all: bin/compr bin/uncompr

bin/compr : compr.cc globaldefs.h bin/stream.o bin/preprocessor.o bin/block.o \
	bin/block_manager.o bin/coders.o bin/dcbwt.o bin/bw_transform.o \
	bin/difference_cover.o bin/prob_models.o
	$(CC) $(FLAGS) -lboost_program_options-mt compr.cc bin/stream.o \
	bin/preprocessor.o bin/block.o bin/block_manager.o bin/coders.o \
	bin/rl_compress.o bin/dcbwt.o bin/bw_transform.o bin/prob_models.o \
	bin/difference_cover.o -o bin/compr

bin/uncompr : uncompr.cc globaldefs.h bin/coders.o bin/rl_compress.o \
	bin/stream.o bin/inverse_bwt.o bin/prob_models.o
	$(CC) $(FLAGS) -lboost_program_options-mt bin/coders.o uncompr.cc \
	bin/rl_compress.o bin/stream.o bin/inverse_bwt.o bin/prob_models.o \
	-o bin/uncompr

# Streams
bin/stream.o : stream.h stream.cc 
	$(CC) $(FLAGS) stream.cc -c -o bin/stream.o

# Preprocessor
bin/preprocessor.o : preprocessors/preprocessor.h preprocessors/preprocessor.cc\
	 stream.h block.h block_manager.h
	$(CC) $(FLAGS) preprocessors/preprocessor.cc -c -o bin/preprocessor.o

bin/testpreprocessor.o : preprocessors/preprocessor.h block.h block_manager.h \
	stream.h preprocessors/test_preprocessor.h \
	preprocessors/test_preprocessor.cc 
	$(CC) $(FLAGS) preprocessors/test_preprocessor.cc -c -o \
	bin/testpreprocessor.o

# Blocks and related things
bin/block.o : block.h block.cc 
	$(CC) $(FLAGS) block.cc -c -o bin/block.o

bin/block_manager.o : block_manager.cc block_manager.h block.h 
	$(CC) $(FLAGS) block_manager.cc -c -o bin/block_manager.o

# Arithmetic range coding
bin/coders.o : coders.cc coders.h bin/rl_compress.o probmodels/base_prob_model.h
	$(CC) $(FLAGS) coders.cc -c -o bin/coders.o

bin/rl_compress.o : rl_compress.cc rl_compress.h globaldefs.h 
	$(CC) $(FLAGS) rl_compress.cc -c -o bin/rl_compress.o 

# Probability models for arithmetic coding
bin/prob_models.o : probmodels/base_prob_model.cc probmodels/base_prob_model.h
	$(CC) $(FLAGS) probmodels/base_prob_model.cc -c -o bin/prob_models.o

# Burrows-wheeler transforms
bin/bw_transform.o : bwtransforms/bw_transform.cc bwtransforms/bw_transform.h \
	bwtransforms/dcbwt.h
	$(CC) $(FLAGS) bwtransforms/bw_transform.cc -c -o bin/bw_transform.o

bin/dcbwt.o : bwtransforms/bw_transform.h bwtransforms/dcbwt.h \
	bwtransforms/difference_cover.h bwtransforms/difference_cover-inl.h \
	bwtransforms/dcbwt.cc 
	$(CC) $(TFLAGS) bwtransforms/dcbwt.cc -c -o bin/dcbwt.o

bin/difference_cover.o :  bwtransforms/difference_cover-inl.h \
	bwtransforms/difference_cover.h bwtransforms/difference_cover.cc
	$(CC) $(FLAGS) bwtransforms/difference_cover.cc -c -o \
	bin/difference_cover.o

bin/inverse_bwt.o : bwtransforms/inverse_bwt.h bwtransforms/inverse_bwt.cc
	$(CC) $(FLAGS) bwtransforms/inverse_bwt.cc -c -o bin/inverse_bwt.o

clean :
	rm -f bin/*
	rm -f test/coderstest
	rm -f test/streamtest
	rm -f test/preproctest
	rm -f test/testfile.txt

# Rest of the file is for tests:
tests : test/streamtest test/preproctest test/coderstest test/dcbwttest \
	test/preprocalgotest
	./test/streamtest
	./test/preproctest
	./test/coderstest
	./test/dcbwttest
	./test/preprocalgotest

test/streamtest : test/stream_test.cc test/testdefs.h bin/stream.o
	$(CC) $(FLAGS) bin/stream.o test/stream_test.cc -lboost_filesystem-mt \
	-o test/streamtest

test/preproctest : test/preproc_test.cc test/testdefs.h bin/block.o \
	bin/preprocessor.o bin/stream.o bin/block_manager.o
	$(CC) $(FLAGS) bin/block.o bin/preprocessor.o bin/stream.o \
	bin/block_manager.o test/preproc_test.cc -o test/preproctest

test/coderstest : test/coders_test.cc test/testdefs.h bin/coders.o \
	bin/stream.o bin/rl_compress.o bin/prob_models.o
	$(CC) $(FLAGS) bin/coders.o bin/rl_compress.o bin/stream.o \
	test/coders_test.cc bin/prob_models.o -o test/coderstest

test/dcbwttest : test/dcbwt_test.cc block.h bwtransforms/dcbwt.h \
	bin/bw_transform.o bin/dcbwt.o bin/block.o bin/difference_cover.o
	$(CC) $(FLAGS) bin/bw_transform.o bin/dcbwt.o test/dcbwt_test.cc \
	bin/block.o bin/difference_cover.o -o test/dcbwttest

test/preprocalgotest : test/preproc_algo_test.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/block.o \
	bwtransforms/dcbwt.h bin/bw_transform.o bin/dcbwt.o \
	bin/difference_cover.o
	$(CC) $(FLAGS) test/preproc_algo_test.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/block.o \
	bin/bw_transform.o bin/dcbwt.o bin/difference_cover.o -o \
	test/preprocalgotest
