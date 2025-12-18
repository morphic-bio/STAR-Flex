# Makefile for trimvalidate - CLI validator tool

CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra -I../../source

LIBTRIM := ../../source/libtrim/libtrim.a

all: trimvalidate

trimvalidate: trimvalidate.cpp $(LIBTRIM)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBTRIM)

$(LIBTRIM):
	$(MAKE) -C ../../source/libtrim

test: trimvalidate
	@./run_parity_test.sh

clean:
	rm -f trimvalidate

.PHONY: all clean test
