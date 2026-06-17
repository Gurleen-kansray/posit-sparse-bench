CXX = g++
CXXFLAGS = -std=c++20 -O2
INCLUDES = -I../universal/include
LDFLAGS = -lm

all: add32_posit bcsstk03_posit bcsstk14_posit bcsstk14_posit64

add32_posit: src/add32_posit.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk03_posit: src/bcsstk03_posit.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk14_posit: src/bcsstk14_posit.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

bcsstk14_posit64: src/bcsstk14_posit64.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

run_all: all
	./add32_posit
	./bcsstk03_posit
	./bcsstk14_posit
	./bcsstk14_posit64

clean:
	rm -f add32_posit bcsstk03_posit bcsstk14_posit bcsstk14_posit64
