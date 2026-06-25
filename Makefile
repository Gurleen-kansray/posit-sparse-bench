CXX = g++
CXXFLAGS = -std=c++20 -O2
INCLUDES = -I../universal/include
LDFLAGS = -lm

all: bcsstk03_ladder bcsstk14_ladder bcsstk38_ladder mhd4800b_ladder \
     bcsstk03_posit bcsstk14_posit bcsstk14_posit64

bcsstk03_ladder: src/bcsstk03_ladder.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk14_ladder: src/bcsstk14_ladder.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk38_ladder: src/bcsstk38_ladder.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

mhd4800b_ladder: src/mhd4800b_ladder.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk03_posit: src/bcsstk03_posit.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk14_posit: src/bcsstk14_posit.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk14_posit64: src/bcsstk14_posit64.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

clean:
	rm -f bcsstk03_ladder bcsstk14_ladder bcsstk38_ladder mhd4800b_ladder \
	      bcsstk03_posit bcsstk14_posit bcsstk14_posit64
