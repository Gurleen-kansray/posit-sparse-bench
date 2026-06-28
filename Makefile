CXX = g++
CXXFLAGS = -std=c++20 -O2
INCLUDES = -I../universal/include
LDFLAGS = -lm

all: generic_ladder cg_compare

generic_ladder: src/generic_ladder.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

cg_compare: src/cg_compare.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

clean:
	rm -f generic_ladder cg_compare
