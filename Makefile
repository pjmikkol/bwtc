CC = g++
DFLAGS = -g -O3 -DNDEBUG -fomit-frame-pointer #-fno-inline
FLAGS = -pedantic -Wextra -Wall $(DFLAGS)
TFLAGS = -Wall $(DFLAGS) # less flags because of template assertions

all: bin bin/compr bin/uncompr 

bin/compr : compr.cc globaldefs.h bin/stream.o bin/preprocessor.o bin/MainBlock.o \
	bin/block_manager.o bin/coders.o bin/dcbwt.o bin/bw_transform.o \
	bin/difference_cover.o bin/prob_models.o bin/utils.o bin/sa-is-bwt.o
	$(CC) $(FLAGS) -lboost_program_options compr.cc bin/stream.o \
	bin/preprocessor.o bin/MainBlock.o bin/block_manager.o bin/coders.o \
	bin/rl_compress.o bin/dcbwt.o bin/bw_transform.o bin/prob_models.o \
	bin/difference_cover.o bin/utils.o bin/sa-is-bwt.o -o bin/compr

bin/uncompr : uncompr.cc globaldefs.h bin/coders.o bin/rl_compress.o \
	bin/stream.o bin/inverse_bwt.o bin/prob_models.o bin/utils.o
	$(CC) $(FLAGS) -lboost_program_options bin/coders.o uncompr.cc \
	bin/rl_compress.o bin/stream.o bin/inverse_bwt.o bin/prob_models.o \
	bin/utils.o -o bin/uncompr

bin :
	mkdir bin

# Streams
bin/stream.o : stream.h stream.cc 
	$(CC) $(FLAGS) stream.cc -c -o bin/stream.o

# Pre- and postprocessors
bin/preprocessor.o : preprocessors/preprocessor.h preprocessors/preprocessor.cc\
	stream.h MainBlock.hpp block_manager.h utils.h
	$(CC) $(FLAGS) preprocessors/preprocessor.cc -c -o bin/preprocessor.o

bin/testpreprocessor.o : preprocessors/preprocessor.h MainBlock.hpp \
	block_manager.h stream.h preprocessors/test_preprocessor.h \
	preprocessors/test_preprocessor.cc 
	$(CC) $(FLAGS) preprocessors/test_preprocessor.cc -c -o \
	bin/testpreprocessor.o

bin/postprocessor.o : preprocessors/postprocessor.cc \
	preprocessors/postprocessor.h utils.h
	$(CC) $(FLAGS) preprocessors/postprocessor.cc -c -o bin/postprocessor.o

bin/sequence_detector.o : preprocessors/sequence_detector.h \
	preprocessors/sequence_detector.cc preprocessors/longsequences.h
	$(CC) $(FLAGS) preprocessors/sequence_detector.cc -c -o \
	bin/sequence_detector.o

bin/expander.o : preprocessors/expander.cc preprocessors/expander.h \
	preprocessors/longsequences.h
	$(CC) $(FLAGS) preprocessors/expander.cc -c -o bin/expander.o

#bin/longsequences.o : preprocessors/longsequences.cc \
#	preprocessors/longsequences.h preprocessors/sequence_detector.h 
#	$(CC) $(FLAGS) preprocessors/longsequences.cc -c -o bin/longsequences.o

bin/long_sequences.o : preprocessors/long_sequences.cc \
	preprocessors/longsequences.h preprocessors/long_sequences.h \
	preprocessors/sequence_detector-inl.h 
	$(CC) $(FLAGS) preprocessors/long_sequences.cc -c -o \
	bin/long_sequences.o

# Blocks and related things
bin/MainBlock.o : MainBlock.hpp MainBlock.cpp 
	$(CC) $(FLAGS) MainBlock.cpp -c -o bin/MainBlock.o

bin/block_manager.o : block_manager.cc block_manager.h MainBlock.hpp 
	$(CC) $(FLAGS) block_manager.cc -c -o bin/block_manager.o

BASICUTILITIES = bin/stream.o bin/preprocessor.o bin/block_manager.o \
	bin/MainBlock.o bin/utils.o

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

bin/sa-is-bwt.o : bwtransforms/sa-is-bwt.cc bwtransforms/sa-is-bwt.h \
	bwtransforms/bw_transform.h bwtransforms/sais.hxx
	$(CC) $(FLAGS) bwtransforms/sa-is-bwt.cc -c -o bin/sa-is-bwt.o

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

bin/utils.o : utils.cc utils.h
	$(CC) $(FLAGS) utils.cc -c -o bin/utils.o

clean :
	rm -f bin/*
	rm -f test/*test
	rm -f test/testfile.txt

# Rest of the file is for tests:
tests : test/preproctest test/coderstest test/dcbwttest  \
	test/preprocalgotest test/streamtest test/sa-is_test \
	test/prime_sequence_detector_test test/mask_sequence_detector_test \
	test/expander_test test/long_sequences_test test/utils_test \
	test/speedtest
	#test/longsequencetest #test/speedtest 
	./test/streamtest
	./test/preproctest
	./test/coderstest
	./test/dcbwttest
	./test/preprocalgotest
	#./test/longsequencetest
	./test/sa-is_test
	./test/prime_sequence_detector_test
	./test/mask_sequence_detector_test
	./test/expander_test
	./test/utils_test

test/utils_test : test/utils_test.cc globaldefs.h bin/utils.o
	$(CC) $(FLAGS) bin/utils.o test/utils_test.cc \
	-lboost_unit_test_framework -o test/utils_test

test/long_sequences_test : bin/long_sequences.o test/long_sequences_test.cc \
	bin/sequence_detector.o $(BASICUTILITIES)
	$(CC) $(FLAGS) bin/long_sequences.o test/long_sequences_test.cc \
	bin/sequence_detector.o $(BASICUTILITIES) -o test/long_sequences_test

test/expander_test : bin/expander.o  test/expander_test.cc
	$(CC) $(FLAGS) -lboost_unit_test_framework test/expander_test.cc \
	bin/expander.o -o test/expander_test

test/prime_sequence_detector_test : preprocessors/sequence_detector-inl.h \
	test/seq_det_test.cc bin/sequence_detector.o bin/utils.o
	$(CC) $(FLAGS) -lboost_unit_test_framework test/seq_det_test.cc \
	bin/sequence_detector.o bin/utils.o -o \
	test/prime_sequence_detector_test

test/mask_sequence_detector_test : preprocessors/sequence_detector-inl.h \
	test/mask_seq_det_test.cc bin/sequence_detector.o bin/utils.o
	$(CC) $(FLAGS) -lboost_unit_test_framework test/mask_seq_det_test.cc \
	bin/sequence_detector.o bin/utils.o -o \
	test/mask_sequence_detector_test

test/sa-is_test : test/sa-is_test.cc bwtransforms/sa-is-bwt.h \
	bwtransforms/sais.hxx
	$(CC) $(FLAGS) test/sa-is_test.cc -o test/sa-is_test

test/streamtest : test/stream_test.cc test/testdefs.h bin/stream.o
	$(CC) $(FLAGS) bin/stream.o test/stream_test.cc -lboost_filesystem \
	-o test/streamtest

test/preproctest : test/preproc_test.cc test/testdefs.h bin/MainBlock.o \
	bin/preprocessor.o bin/stream.o bin/block_manager.o bin/utils.o
	$(CC) $(FLAGS) bin/MainBlock.o bin/preprocessor.o bin/stream.o \
	bin/block_manager.o bin/utils.o test/preproc_test.cc -o \
	test/preproctest

test/coderstest : test/coders_test.cc test/testdefs.h bin/coders.o \
	bin/stream.o bin/rl_compress.o bin/prob_models.o bin/utils.o
	$(CC) $(FLAGS) bin/coders.o bin/rl_compress.o bin/stream.o \
	test/coders_test.cc bin/prob_models.o bin/utils.o -o test/coderstest

test/dcbwttest : test/dcbwt_test.cc MainBlock.hpp bwtransforms/dcbwt.h \
	bin/sa-is-bwt.o \
	bin/bw_transform.o bin/dcbwt.o bin/MainBlock.o bin/difference_cover.o
	$(CC) $(FLAGS) bin/bw_transform.o bin/dcbwt.o test/dcbwt_test.cc \
	bin/MainBlock.o bin/difference_cover.o bin/sa-is-bwt.o -o test/dcbwttest

# When long sequences ready plug into this
test/speedtest : test/bwt_and_preproctest.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/MainBlock.o \
	bin/utils.o bin/postprocessor.o bin/dcbwt.o bin/bw_transform.o \
	bin/difference_cover.o bin/sa-is-bwt.o
	$(CC) $(FLAGS) test/bwt_and_preproctest.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/MainBlock.o \
	bin/utils.o bin/postprocessor.o bin/dcbwt.o bin/bw_transform.o \
	bin/difference_cover.o bin/sa-is-bwt.o -o test/speedtest

test/preprocalgotest : test/preproc_algo_test.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/MainBlock.o \
	bwtransforms/dcbwt.h bin/bw_transform.o bin/dcbwt.o bin/sa-is-bwt.o \
	bin/difference_cover.o bin/postprocessor.o bin/utils.o
	$(CC) $(FLAGS) test/preproc_algo_test.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/MainBlock.o \
	bin/bw_transform.o bin/sa-is-bwt.o bin/dcbwt.o bin/difference_cover.o \
	bin/postprocessor.o bin/utils.o -o test/preprocalgotest

test/longsequencetest : test/longsequence_test.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/MainBlock.o \
	bin/longsequences.o bin/utils.o bin/postprocessor.o
	$(CC) $(FLAGS) test/longsequence_test.cc bin/testpreprocessor.o \
	bin/block_manager.o bin/preprocessor.o bin/stream.o bin/MainBlock.o \
	bin/longsequences.o bin/utils.o bin/postprocessor.o \
	-o test/longsequencetest
