all: test_subset_sum

test_subset_sum: test_subset_sum.cc subset_sum.o pointer.o
	g++ -o test_subset_sum test_subset_sum.cc subset_sum.o pointer.o

subset_sum.o: subset_sum.cc subset_sum.h
	g++ -c subset_sum.cc

pointer.o: pointer.cc pointer.h
	g++ -c pointer.cc
